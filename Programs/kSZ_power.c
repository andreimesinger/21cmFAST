/*************************************************************************************
Written mostly by Matt McQuinn.

This program calculates the kSZ from reionization boxes.  The script that executes this program should have the
format:

FORMAT: kSZ_power <list of density boxes> <list of xH boxes> <list of velocity boxes> <ouput filename>

old version: kSZ_power <filename_that_lists_reionization_data_files. <number_of_files_to_read_in $output_file_name>

$filename_with_list_of_reionization_data_files has the following format:
  
On each of the $number_of_files_to_read_in lines in this file, the format is

        $path/$base_file_name    $starting_redshift

with the files in order of increasing redshift.
and each $base_file_name comes with three files:

          base_file_name.dens = overdensity field
          base_file_name.vel = velocity field (in units of km/s)
          base_file_name.ion = ionization field (ones and zeros)

Each of these files is $data_grid_size^3 floats in row major order and the los direction is taken to be 
in k where float_number = i*GRIDHII_DIM*GRIDHII_DIM + j*GRIDHII_DIM +k.  Note that floats are wasteful, but
you cannot convert everything to float blindly because fftwf by default wants float format. I can help 
if you would rather do floats.

See kSZ.sh for an example script file, and filenames.txt is an example of $number_of_files_to_read_in.

Currently this code assumes hydrogen and 1 electron of helium were reionized simultaneously.

I do not put in exp(-tau) factor that attenuates the anisotropies.  We need to eventually put this in.  
Should not bee too hard to do this. I note where it needs to go in this code.

There are two outputs:
The temperature field field: $output_file_name.dat
The angular power spectrum of kSZ: $output_file_name.pk   (l l^2 C_l/2pi error)
Because this is the fourier transform of a 2D field, .pk is noisy when plotted in l spaced at the nyquist frequency.
We will eventually want to bin it.

Note that at low-l you should see a large bump in the power spectrum that looks like a bug.  It is not, and called
the doppler anisotropy. 

Compiler flags:
TAU_POWER: this makes it so that this program calculates
the tau field and its power spectrum rather than the kSZ.  
PARALLEL_APPROX: assumes that rays run parallel to line-of-sight vector,
which is a reasonable approximation at high z
NO_DELTA: Calculate ignoring density fluctuations
ASCII: .dat file is ASCII
*************************************************************************************/
#include "fftCMB.c"

#define PARALLEL_APPROX (int) 0 // use parallel light ray approximation

void projectArray(float *Tcmb, float *Tcmb_3d, float *dtau_3d);
int readInArrays(float*, float*, FILE*, FILE*, FILE*);
void printTField(float *Tcmb);
void printStatistics(float *K, float *Pk, int *Count);
unsigned long long position(int i, int j, int ti, int tj);
void subtract_avg(float *Tcmb);
double hubbleZ(float z);
double confDist(float z);

float START_REDSHIFT = -1; //lowest redshift

float REDSHIFT;
double currZ;
double mean_taue_curr_z = 0;
FILE *LOG;
char OUTFILE[400];
double *taue_arry;
void doFFT(float *R, fftwf_complex *Rf, long Nf, int dim);
void makePk(float *K, float *Pk, int *Count, fftwf_complex *Rf);
void initialize(float **K,float **Pk, int **count, int flag);
  long N, Nf;
  fftwf_complex *Tcmbf = NULL;
  float *Tcmb, *Tcmb_3D, *dtau_3D;
  float *Cl = NULL,*L = NULL;
  int *Count = NULL;


int main(int argc, char *argv[])
{
  unsigned short seed[3] = {-11, 444, 3};//random number 
  unsigned long long i, ct;
  char filename[400];
  void doInverseFFT(fftwf_complex *Rf, float *R, long Nf, int dim);
  FILE *IN, *delta_filelist, *v_filelist, *xH_filelist;


  if (argc != 5){
    fprintf(stderr, "kSZ_power: FORMAT: kSZ_power <list of density boxes> \
          <list of xH boxes> <list of velocity boxes> <ouput filename>\nAborting\n");
    return -1;
  }
  strcpy(OUTFILE, argv[4]);
  seed48(seed);//initialize random number seed
  //size of arrays for 2Dfft
  N = HII_DIM*HII_DIM;
  Nf = (HII_DIM/2 +1)*HII_DIM;

  taue_arry = (double *) malloc(sizeof(double)*HII_DIM*HII_DIM);

  /****  OPEN FILES  ****/
  if (!(LOG=fopen("../Log_files/kSZ_power_LOG", "w"))){
    fprintf(stderr, "kSZ_power: WARNING: unable to open LOG file at \
                       ../Log_files/kSZ_power_LOG.\n");
  }
  if (!(delta_filelist=fopen(argv[1], "r"))){
    fprintf(stderr, "kSZ_power: ERROR: unable to open delta filelist %s.\n", argv[1]);
    fprintf(LOG, "kSZ_power: ERROR: unable to open delta filelist %s.\n", argv[1]);
    fclose(LOG); return -1;
  }
  if (!(xH_filelist=fopen(argv[2], "r"))){
    fprintf(stderr, "kSZ_power: ERROR: unable to open xHI filelist %s.\n", argv[2]);
    fprintf(LOG, "kSZ_power: ERROR: unable to open xHI filelist %s.\n", argv[2]);
    fclose(LOG); fclose(delta_filelist); return -1;
  }
  if (!(v_filelist=fopen(argv[3], "r"))){
    fprintf(stderr, "kSZ_power: ERROR: unable to open xHI filelist %s.\n", argv[3]);
    fprintf(LOG, "kSZ_power: ERROR: unable to open xHI filelist %s.\n", argv[3]);
    fclose(LOG); fclose(delta_filelist); fclose(xH_filelist); return -1;
  }
  fprintf(stderr, "kSZ_power: Opened filelist files\n");
  fprintf(LOG, "kSZ_power: Opened filelist files\n");
  
  /****  ALLOCATE MEMORY  ****/
  //allocates 2 and 3 temperature field and initializes 2D 
  dtau_3D = (float*)  fftwf_malloc(sizeof(float) * HII_TOT_NUM_PIXELS);
  Tcmb = (float*) fftwf_malloc(sizeof(float) * N);
  Tcmb_3D = (float*)  fftwf_malloc(sizeof(float) * N*HII_DIM);
  if (!dtau_3D || !Tcmb || !Tcmb_3D){ // sloppy and memory bleeding.. i know. but how often do you get allocation errors?
    fprintf(stderr, "kSZ_power: ERROR allocating memory.\nAborting\n");
    fprintf(LOG, "kSZ_power: ERROR allocating memory.\nAborting\n");
    fclose(LOG); fclose(delta_filelist); fclose(xH_filelist); fclose(v_filelist); return -1;
  }
  for(i = 0; i< N; i++)
    Tcmb[i] = 0.;
  //fourier conjugate
  Tcmbf = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nf);
  initialize(&L, &Cl, &Count, 0);


  /************************  END INITIALIZATION    *****************************************/

  //    readInFileNames(filename, 0);//initializes function that provides filename to read in
  //    for(i = 0; i < NUM_FILES; i++)//each file is a different redshift interval
  while (!feof(delta_filelist)){      
    fprintf(stderr, "Working on the next set of boxes...\n");
    fprintf(LOG, "Working on the next set of boxes...\n");

    //read in boxes
    if (readInArrays(Tcmb_3D, dtau_3D, delta_filelist, xH_filelist, v_filelist) < 0 ){
      fclose(LOG); fclose(delta_filelist); fclose(xH_filelist); fclose(v_filelist);
      return -1;	
    }

    fprintf(stderr, "Initialized the arrays, now doing the projection\n");
    fprintf(LOG, "Initialized the arrays, now doing the projection\n");

    //calculates kSZ (or tau signal)
    projectArray(Tcmb, Tcmb_3D, dtau_3D); 
    /*
    //print power spectrum
    doFFT(Tcmb, Tcmbf, Nf, 2);
    makePk(L, Cl, Count, Tcmbf); //calculates power spectrum
    printStatistics(L, Cl, Count); //prints power spectrum
    initialize(&L, &Cl, &Count, 3);
    */
  }
  fprintf(stderr, "Done with all boxes\n Now doing FFT and printing the power spectrum\n");
  fprintf(LOG, "Done with all boxes\n Now doing FFT and printing the power spectrum\n");
  fclose(delta_filelist); fclose(xH_filelist); fclose(v_filelist);
  sprintf(filename, "cp %s_z%f.pk %s.pk", OUTFILE, REDSHIFT, OUTFILE);
  fprintf(stderr, "Copying the last power spectrum file: %s\n", filename);
  fprintf(LOG, "Copying the last power spectrum file: %s\n", filename);
  system(filename);
  fftwf_free(Tcmb_3D); //frees grid
  fftwf_free(dtau_3D);
  printTField(Tcmb); //prints the kSZ (tau) field


  //frees all the rest of the gllocated memory
  initialize(&L, &Cl, &Count, 1);
  fftwf_free(Tcmbf);
  fftwf_free(Tcmb);
  fclose(LOG);
  free(taue_arry);

  return 0;
}


/********************************************************************
Calculates 2D kSZ field
*******************************************************************/
void projectArray(float *Tcmb, float *Tcmb_3d, float *dtau_3d)
{
  int i, j, k, ii, jj;
  double A, dz, inc;
  int tx, ty;
  double taue, DA_zstart = confDist(START_REDSHIFT);
  double DA_zstartcurrbox = confDist(REDSHIFT);
  double dR = (BOX_LEN / (double) HII_DIM); // comoving length of single cell in cMpc
  double last_slice_taue = 0;
  unsigned long long ct;


  //this coeficinet is used to calculate kSZ
  A = N_b0*SIGMAT*BOX_LEN*CMperMPC/(double)HII_DIM; // N*sigma*dR (comoving units)

  //randomly chooses location with which to start the projection
  //different for each sub-box so structures do not align across interface
  tx = (int )(HII_DIM*drand48());   //Put in your favorite random # generator if you like such as ran1(&idum) in numerical recipes
  ty = (int )(HII_DIM*drand48());


  currZ = REDSHIFT;
  taue = mean_taue_curr_z; //optical depth to thompson scattering along this k sightline
  for(k = 0; k < HII_DIM; k++){
    inc = (k*dR + DA_zstartcurrbox) / DA_zstart; //ratio of the angular diameter distances	  
    dz = - dR * CMperMPC / drdz(currZ);
    currZ += dz;

    for(i = 0; i < HII_DIM; i++){
      //      fprintf(stderr, "%i ", i);
      for(j = 0; j < HII_DIM; j++){ 

	if (PARALLEL_APPROX){
	  ii = i; jj = j; //assumes rays are parallel
	}
	else{
	  ii = ((int) ((i+.5)*inc)) % HII_DIM; //does not
	  jj = ((int) ((j+.5)*inc)) % HII_DIM;
	}

	  //adding missing conversion factors for length units
	taue_arry[position(i,j,0,0)] += dtau_3d[HII_R_INDEX(ii,jj,k)] * (1+currZ) * (1+currZ);
	    
	  Tcmb[position(i,j, tx, ty)] 
	    += A*(1+currZ)*Tcmb_3d[HII_R_INDEX(ii,jj,k)] * exp(-taue_arry[position(i,j,0,0)]);
	  if ((i==0)&&(j==0)&&(k<110) && (100<k)){
	    fprintf(stderr, "taue along sightline is %e\n", taue_arry[position(i,j,0,0)]);
	    fprintf(LOG, "taue along sightline is %e\n", taue_arry[position(i,j,0,0)]);
	  }
	  /**** AM: I took out one factor of (1+z) since the velocity is in comoving units **/
	}

    }// done with all x,y

    if (k%50 ==0){
      //print power spectrum
      doFFT(Tcmb, Tcmbf, Nf, 2);
      makePk(L, Cl, Count, Tcmbf); //calculates power spectrum
      printStatistics(L, Cl, Count); //prints power spectrum
      initialize(&L, &Cl, &Count, 3);
    }
  } // end look through box

  mean_taue_curr_z =0;
  for (ct=0; ct<HII_DIM*HII_DIM; ct++)
    mean_taue_curr_z += taue_arry[ct];
  mean_taue_curr_z /=  pow(HII_DIM, 2); // updates the mean taue for next box

  //#ifdef TAU_POWER  //removes mean--kSZ naturally has zero mean, so only relevant for tau
   subtract_avg(Tcmb);  
//#endif


  //prints out starting and ending redshift
  fprintf(stderr, "zstart =%f zend = %f, <taue(zend)>=%f\n", REDSHIFT, currZ, mean_taue_curr_z);
  fprintf(LOG, "zstart =%f zend = %f, <taue(zend)>=%f\n", REDSHIFT, currZ, mean_taue_curr_z);
}

/************************************************************************
Returns position, where ti and tj are the center position
***********************************************************************/
unsigned long long position(int i, int j, int ti, int tj)
{
    return ((i+ti)%HII_DIM)*HII_DIM + (j +tj)%HII_DIM; 
}


/********************************************************************
Reads in necessary arrays to calculate kSZ/tau power spectrum.
**********************************************************************/
int readInArrays(float *Tcmb_3d, float *dtau_3d, 
		 FILE *delta_filelist, FILE *xH_filelist, FILE *v_filelist)
{
  unsigned long long i, ct;
    int x,y,z;
    float v, delta, xi;
    FILE *delta_IN, *xH_IN, *v_IN;
    char filename[200], *token;

    //density field in overdensity
    if (fscanf(delta_filelist, "%s\n", filename) <= 0){
      fprintf(stderr, "kSZ_power: ERROR: early termination in delta filelist\n.");
      fprintf(LOG, "kSZ_power: ERROR: early termination in delta filelist\n.");
      return -1;
    }
    if( !(delta_IN = fopen(filename, "rb"))){
	fprintf(stderr, "Could not read in %s.\n", filename);
	fprintf(LOG, "Could not read in %s.\n", filename);
	return -1;
    }

    //neutral fraction field
    if (fscanf(xH_filelist, "%s\n", filename) <= 0){
      fprintf(stderr, "kSZ_power: ERROR: early termination in delta filelist\n.");
      fprintf(LOG, "kSZ_power: ERROR: early termination in delta filelist\n.");
      return -1;
    }
    if( !(xH_IN = fopen(filename, "rb"))){
	fprintf(stderr, "Could not read in %s.\n", filename);
	fprintf(LOG, "Could not read in %s.\n", filename);
	fclose (delta_IN); return -1;
    }

    //velocity field
    if (fscanf(v_filelist, "%s\n", filename) <= 0){
      fprintf(stderr, "kSZ_power: ERROR: early termination in delta filelist\n.");
      fprintf(LOG, "kSZ_power: ERROR: early termination in delta filelist\n.");
      return -1;
    }
    if( !(v_IN = fopen(filename, "rb"))){
	fprintf(stderr, "Could not read in %s.\n", filename);
	fprintf(LOG, "Could not read in %s.\n", filename);
	fclose(delta_IN); fclose(xH_IN); return -1;
    }

    // get start redshift
    if (START_REDSHIFT < 0){
      token = strtok(filename, "_");
      token = strtok(NULL, "_");
      token = strtok(NULL, "t");
      token = strtok(NULL, "t");
      token = strtok(NULL, "_");
      START_REDSHIFT = atof(token);
      fprintf(stderr, "Starting redshift: %f\n", START_REDSHIFT);
      fprintf(LOG, "Starting redshift: %f\n", START_REDSHIFT);
      REDSHIFT = START_REDSHIFT;
      mean_taue_curr_z = tau_e(0, START_REDSHIFT, NULL, NULL, 0);
      for (ct=0; ct< HII_DIM*HII_DIM; ct++)
	taue_arry[ct] = mean_taue_curr_z;
    }
    else{ // get the starting one from this current box
      token = strtok(filename, "_");
      token = strtok(NULL, "_");
      token = strtok(NULL, "t");
      token = strtok(NULL, "t");
      token = strtok(NULL, "_");
      REDSHIFT= atof(token);
      fprintf(stderr, "Starting redshift of new box: %f\n", REDSHIFT);
      fprintf(LOG, "Starting redshift of new box: %f\n", REDSHIFT);
    }

    // now load the arrays
    for(i = 0; i < HII_TOT_NUM_PIXELS; i++){  
      if (fread(&delta, sizeof(float), 1, delta_IN) != 1){
	fprintf(stderr, "kSZ_power: ERROR: reading from delta box\n.");
	fprintf(LOG, "kSZ_power: ERROR: reading from delta box\n.");
	fclose(delta_IN); fclose(xH_IN); fclose(v_IN); return -1;
      }
      if (fread(&v, sizeof(float), 1, v_IN)	 != 1){
	fprintf(stderr, "kSZ_power: ERROR: reading from velocity box\n.");
	fprintf(LOG, "kSZ_power: ERROR: reading from velocity box\n.");
	fclose(delta_IN); fclose(xH_IN); fclose(v_IN); return -1;
      }

      if (fread(&xi, sizeof(float), 1, xH_IN) != 1) {
	fprintf(stderr, "kSZ_power: ERROR: reading from xH box\n.");
	fprintf(LOG, "kSZ_power: ERROR: reading from xH box\n.");
	fclose(delta_IN); fclose(xH_IN); fclose(v_IN); return -1;
      }
	xi = 1.0-xi; // input is neutral fraction not ionized
	v *= CMperMPC/C; //in units of C
	//*****AM:  my velocity fields are in comoving units
	//*****will convert to proper units in the projection
	// check for wierd rounding errors
	if (delta < -1){
	  fprintf(stderr, "kSZ_power: flagged delta=%e < -1 at index=%llu\nChanging to -1+FRACT_FLOAT_ERR\n", delta, i);
	  fprintf(LOG, "kSZ_power: flagged delta=%e < -1 at index=%llu\nChanging to -1+FRACT_FLOAT_ERR\n", delta, i);
	  delta = -1 + FRACT_FLOAT_ERR;
	}
	if (xi < 0){
	  fprintf(stderr, "kSZ_power: flagged xi=%e < 0 at index=%llu\nChanging to 0\n", xi, i);
	  fprintf(LOG, "kSZ_power: flagged xi=%e < 0 at index=%llu\nChanging to 0\n", xi, i);
	  xi = 0;
	}
	if (xi > 1){
	  fprintf(stderr, "kSZ_power: flagged xi=%e > 1 at index=%llu\n.Changing to 1\n", xi, i);
	  fprintf(LOG, "kSZ_power: flagged xi=%e > 1 at index=%llu\n.Changing to 1\n", xi, i);
	  xi = 1;
	}

	Tcmb_3d[i] = v*xi*(1.0+delta);


	// and the e^dtau array, i.e.  optical depth contribution of the current cell
	dtau_3d[i] =  N_b0*(1.0+delta)*xi*SIGMAT*CMperMPC*BOX_LEN/(double)HII_DIM;
	//AM: NOTE THAT it is missing a factor of (1+z)^2 from proper to comoving conversions
    }


    fclose(delta_IN); fclose(v_IN); fclose(xH_IN);
}


//Prints dT_kSZ/tau field
void printTField(float *Tcmb)
{
    char filename[400];
    FILE *outfile;
    sprintf(filename, "%s.dat", OUTFILE);
    outfile = fopen(filename, "wb");
    /*
    int i,j;
    for(i=0; i< HII_DIM;i++)
      for(j=0; j< HII_DIM;j++)
	{
	  fprintf(outfile, "%le\n", Tcmb[i*HII_DIM+j]);
	}   
    */
    fwrite(Tcmb, sizeof(float), HII_DIM*HII_DIM, outfile);
    fprintf(stderr, "Just called fwrite to write file to %s\n", filename);
    fprintf(LOG, "Just called fwrite to write file to %s\n", filename);
    fclose(outfile);
}


//prints angular power spectrum
void printStatistics(float *K, float *Pk, int *Count) 
{
    int i;
    FILE *outfile;
    float coef, x;
    char filename[400];
   

    sprintf(filename, "%s_z%f.pk", OUTFILE, currZ);

    if((outfile = fopen(filename, "w")) == NULL)
    {
	fprintf(stderr, "Could not open output file.\n");
	exit(1);

    }

    x = confDist(START_REDSHIFT);//I binned in terms of size of cell at START_REDSHIFT 
    fprintf(outfile, "#prints l, l^2 C_l/2pi, L^2 Cl/2pi sqrt(2/Nmodes)\n");
    for(i =  1; i< HII_DIM*sqrt(2)/2.; i++)
    {
	coef = K[i]*K[i]/(2.*M_PI);
	fprintf(outfile, "%f %e %e\n", x*K[i], coef*Pk[i],
		coef*Pk[i]/sqrt(Count[i]/2.)); //prints l, l^2 C_l/2pi, L^2 Cl/2pi sqrt(2/Nmodes)
    }
    fclose(outfile);
}



/****************************************************
Subtracts the mean of field (zero mode)
needed to calculate power of delta tau
****************************************************/
void subtract_avg(float *Tcmb)
{
    int i, j;
    double avg = 0.;
    for(i = 0; i < HII_DIM; i++)
	for(j = 0; j < HII_DIM; j++)
	    avg += Tcmb[i*HII_DIM + j];

    avg/= (HII_DIM*HII_DIM);
    fprintf(stderr, "mean is %e\n", avg);
    for(i = 0; i < HII_DIM; i++)
	for(j = 0; j < HII_DIM; j++)
	  Tcmb[i*HII_DIM + j] -= (float)avg;
}

/*********************************************************
Conformal distance :

Andrei, feel free to replace this with more accurate formula
Didn't want to include my cosmology calculator in this code.
This assues OmegaM = 0.27
********************************************************/
double confDist(float z)
{
  return comovingdistance(0, z)/CMperMPC;
  //return 6689./hlittle*pow((1.+z)/10., .2); //cMpc //AM: replaced
}
//Hubble(z)/H0
double hubbleZ(float z)
{
  return hubble(z)/(double)Ho;
  //   return sqrt(OMm*pow(1.+z, 3.) + (1.-OMm)); //AM: replaced
}
