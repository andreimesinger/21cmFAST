#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "filter.c"

FILE *LOG;

/*
  USAGE: find_halos <redshift>

  Program FIND_HALOS takes in a k_space box of the linear overdensity field
  and filters it on decreasing scales in order to find virialized halos.  
  Virialized halos are defined according to the linear critical overdensity.
  FIND_HALOS outputs a list of virialized halos, <halos_N_LMpc>, where
  N^3 is the spacial resolution and L is the comoving length of the box.
  The columns in <halos_N_LMpc> contains each halo's:
  mass (in M_sun)
  x location (in box size, [0 to 1) )
  y location (in box size, [0 to 1) )
  z location (in box size, [0 to 1) )

  NOTE: Relevant parameters are taken from INIT_PARAMS.H and ANAL_PARAMS.H

  Author: Andrei Mesinger
  Date: 10/29/06
*/


/*
  Funtion OVERLAP_HALO checks if the would be halo with radius R
  and centered on (x,y,z) overlaps with a preesisting halo
*/
int overlap_halo(char * in_halo, float R, int x, int y, int z){
  int x_curr, y_curr, z_curr, x_min, x_max, y_min, y_max, z_min, z_max, R_index;
  float Rsq_curr_index, xsq, xplussq, xminsq, ysq, yplussq, yminsq, zsq, zplussq, zminsq;
  int x_index, y_index, z_index;

    fprintf(LOG, "begin overlap halo (%i, %i, %i), clock=%.2f\n", x,y,z,(double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);

  // scale R to a effective overlap size, using R_OVERLAP_FACTOR
  R *= R_OVERLAP_FACTOR;

  // convert R to index units
  R_index = ceil(R/BOX_LEN*DIM);
  Rsq_curr_index = pow(R/BOX_LEN*DIM, 2); // convert to index

  // set parameter range
  x_min = x-R_index;
  x_max = x+R_index;
  y_min = y-R_index;
  y_max = y+R_index;
  z_min = z-R_index;
  z_max = z+R_index;

  //    printf("min %i, %i, %i\n", x_min, y_min, z_min);
  //printf("max %i, %i, %i\n", x_max, y_max, z_max);
  for (x_curr=x_min; x_curr<=x_max; x_curr++){
    for (y_curr=y_min; y_curr<=y_max; y_curr++){
      for (z_curr=z_min; z_curr<=z_max; z_curr++){
	x_index = x_curr;
	y_index = y_curr;
	z_index = z_curr;
	// adjust if we are outside of the box
	if (x_index<0) {x_index += DIM;}
	else if (x_index>=DIM) {x_index -= DIM;}
	if (y_index<0) {y_index += DIM;}
	else if (y_index>=DIM) {y_index -= DIM;}
	if (z_index<0) {z_index += DIM;}
	else if (z_index>=DIM) {z_index -= DIM;}

	// remember to check all reflections
	xsq = pow(x-x_index, 2);
	ysq = pow(y-y_index, 2);
	zsq = pow(z-z_index, 2);
	xplussq = pow(x-x_index+DIM, 2);
	yplussq = pow(y-y_index+DIM, 2);
	zplussq = pow(z-z_index+DIM, 2);
	xminsq = pow(x-x_index-DIM, 2);
	yminsq = pow(y-y_index-DIM, 2);
	zminsq = pow(z-z_index-DIM, 2);
	if ( in_halo[R_INDEX(x_index, y_index, z_index)] && 
	     ( (Rsq_curr_index > (xsq + ysq + zsq)) || // AND pixel is within this halo
	       (Rsq_curr_index > (xsq + ysq + zplussq)) || 
	       (Rsq_curr_index > (xsq + ysq + zminsq)) || 
			   
	       (Rsq_curr_index > (xsq + yplussq + zsq)) || 
	       (Rsq_curr_index > (xsq + yplussq + zplussq)) || 
	       (Rsq_curr_index > (xsq + yplussq + zminsq)) || 

	       (Rsq_curr_index > (xsq + yminsq + zsq)) || 
	       (Rsq_curr_index > (xsq + yminsq + zplussq)) || 
	       (Rsq_curr_index > (xsq + yminsq + zminsq)) || 


	       (Rsq_curr_index > (xplussq + ysq + zsq)) || 
	       (Rsq_curr_index > (xplussq + ysq + zplussq)) || 
	       (Rsq_curr_index > (xplussq + ysq + zminsq)) || 

	       (Rsq_curr_index > (xplussq + yplussq + zsq)) || 
	       (Rsq_curr_index > (xplussq + yplussq + zplussq)) || 
	       (Rsq_curr_index > (xplussq + yplussq + zminsq)) || 

	       (Rsq_curr_index > (xplussq + yminsq + zsq)) || 
	       (Rsq_curr_index > (xplussq + yminsq + zplussq)) || 
	       (Rsq_curr_index > (xplussq + yminsq + zminsq)) || 


	       (Rsq_curr_index > (xminsq + ysq + zsq)) || 
	       (Rsq_curr_index > (xminsq + ysq + zplussq)) || 
	       (Rsq_curr_index > (xminsq + ysq + zminsq)) || 

	       (Rsq_curr_index > (xminsq + yplussq + zsq)) || 
	       (Rsq_curr_index > (xminsq + yplussq + zplussq)) || 
	       (Rsq_curr_index > (xminsq + yplussq + zminsq)) || 

	       (Rsq_curr_index > (xminsq + yminsq + zsq)) || 
	       (Rsq_curr_index > (xminsq + yminsq + zplussq)) || 
	       (Rsq_curr_index > (xminsq + yminsq + zminsq))
	       ) ){
	    
	  // this pixel already belongs to a halo, and would want to become part of this halo as well
	  return 1;
	}
      }
    }
  }

  return 0;
}


/*
  Funtion UPDATE_IN_HALO takes in a box <in_halo> and flags all points
  which fall within radius R of (x,y,z).
*/
void update_in_halo(char * in_halo, float R, int x, int y, int z){
  int x_curr, y_curr, z_curr, x_min, x_max, y_min, y_max, z_min, z_max, R_index;
  float Rsq_curr_index, xsq, xplussq, xminsq, ysq, yplussq, yminsq, zsq, zplussq, zminsq;
  int x_index, y_index, z_index;

    fprintf(LOG, "begin update halo (%i, %i, %i), clock=%.2f\n", x,y,z,(double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);

  // convert R to index units
  R_index = ceil(R/BOX_LEN*DIM);
  Rsq_curr_index = pow(R/BOX_LEN*DIM, 2); // convert to index

  // set parameter range
  x_min = x-R_index;
  x_max = x+R_index;
  y_min = y-R_index;
  y_max = y+R_index;
  z_min = z-R_index;
  z_max = z+R_index;

  //printf("min %i, %i, %i\n", x_min, y_min, z_min);
  //printf("max %i, %i, %i\n", x_max, y_max, z_max);
  for (x_curr=x_min; x_curr<=x_max; x_curr++){
    for (y_curr=y_min; y_curr<=y_max; y_curr++){
      for (z_curr=z_min; z_curr<=z_max; z_curr++){
	x_index = x_curr;
	y_index = y_curr;
	z_index = z_curr;
	// adjust if we are outside of the box
	if (x_index<0) {x_index += DIM;}
	else if (x_index>=DIM) {x_index -= DIM;}
	if (y_index<0) {y_index += DIM;}
	else if (y_index>=DIM) {y_index -= DIM;}
	if (z_index<0) {z_index += DIM;}
	else if (z_index>=DIM) {z_index -= DIM;}

	// now check
	if (!in_halo[R_INDEX(x_index, y_index, z_index)]){ // untaken pixel (not part of other halo)
	  // remember to check all reflections
	  xsq = pow(x-x_index, 2);
	  ysq = pow(y-y_index, 2);
	  zsq = pow(z-z_index, 2);
	  xplussq = pow(x-x_index+DIM, 2);
	  yplussq = pow(y-y_index+DIM, 2);
	  zplussq = pow(z-z_index+DIM, 2);
	  xminsq = pow(x-x_index-DIM, 2);
	  yminsq = pow(y-y_index-DIM, 2);
	  zminsq = pow(z-z_index-DIM, 2);
	  if ( (Rsq_curr_index > (xsq + ysq + zsq)) || 
	       (Rsq_curr_index > (xsq + ysq + zplussq)) || 
	       (Rsq_curr_index > (xsq + ysq + zminsq)) || 

	       (Rsq_curr_index > (xsq + yplussq + zsq)) || 
	       (Rsq_curr_index > (xsq + yplussq + zplussq)) || 
	       (Rsq_curr_index > (xsq + yplussq + zminsq)) || 

	       (Rsq_curr_index > (xsq + yminsq + zsq)) || 
	       (Rsq_curr_index > (xsq + yminsq + zplussq)) || 
	       (Rsq_curr_index > (xsq + yminsq + zminsq)) || 


	       (Rsq_curr_index > (xplussq + ysq + zsq)) || 
	       (Rsq_curr_index > (xplussq + ysq + zplussq)) || 
	       (Rsq_curr_index > (xplussq + ysq + zminsq)) || 

	       (Rsq_curr_index > (xplussq + yplussq + zsq)) || 
	       (Rsq_curr_index > (xplussq + yplussq + zplussq)) || 
	       (Rsq_curr_index > (xplussq + yplussq + zminsq)) || 

	       (Rsq_curr_index > (xplussq + yminsq + zsq)) || 
	       (Rsq_curr_index > (xplussq + yminsq + zplussq)) || 
	       (Rsq_curr_index > (xplussq + yminsq + zminsq)) || 


	       (Rsq_curr_index > (xminsq + ysq + zsq)) || 
	       (Rsq_curr_index > (xminsq + ysq + zplussq)) || 
	       (Rsq_curr_index > (xminsq + ysq + zminsq)) || 

	       (Rsq_curr_index > (xminsq + yplussq + zsq)) || 
	       (Rsq_curr_index > (xminsq + yplussq + zplussq)) || 
	       (Rsq_curr_index > (xminsq + yplussq + zminsq)) || 

	       (Rsq_curr_index > (xminsq + yminsq + zsq)) || 
	       (Rsq_curr_index > (xminsq + yminsq + zplussq)) || 
	       (Rsq_curr_index > (xminsq + yminsq + zminsq))
	     ){
	    
	    // we are within the sphere defined by R, so change flag in in_halo array
	    in_halo[R_INDEX(x_index, y_index, z_index)] = 1;
	    //	    printf("%i, %i, %i\n", x_index, y_index, z_index);
	  }
	}
      }
    }
  }
}



int main(int argc, char ** argv){
  fftwf_complex *box;
  fftwf_plan plan;
  FILE *IN, *OUT, *F;
  float growth_factor, R, delta_m, dm, dlnm, M, Delta_R, delta_crit, REDSHIFT;
  double fgrtm, dfgrtm;
  unsigned long long ct;
  char filename[80], *in_halo, *forbidden;
  int x,y,z,dn, n;
  float R_temp, x_temp, y_temp, z_temp, dummy, M_MIN;

  /************  BEGIN INITIALIZATION ****************************/
  if (argc != 2){
    fprintf(stderr, "USAGE: find_halos <redshift>\nAborting...\n");
    return -1;
  }
  REDSHIFT = atof(argv[1]);

  system("mkdir ../Log_files");
  system("mkdir ../Output_files");
  system("mkdir ../Output_files/DNDLNM_files");
  system("mkdir ../Output_files/FgtrM_files");
  system("mkdir ../Output_files/Halo_lists");
  
  // initialize power spectrum 
  init_ps(0, 1e10);
  growth_factor = dicke(REDSHIFT); // normalized to 1 at z=0
  delta_crit = Deltac; // for now set to spherical; check if we want elipsoidal later

  //set the minimum source mass
  M_MIN = M_TURNOVER;

  // open log file
  system("mkdir ../Log_files");
  LOG = fopen("../Log_files/zscroll_log_file", "w");
  if (!LOG){
    fprintf(stderr, "find_halos.c: Unable to open log file\n Aborting...\n");
    return -1;
  }

  // allocate array for the k-space box
  box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);
  if (!box){
    fprintf(stderr, "find_halos.c: Error allocating memory for box.\nAborting...\n");
    return -1;
  }

  // allocate memory for the boolean in_halo box
  in_halo = (char *) malloc(sizeof(char)*TOT_NUM_PIXELS);
  if (!in_halo){
    fprintf(stderr, "find_halos.c: Error allocating memory for in_halo box\nAborting...\n");
    fftwf_free(box);
    return -1;
  }
  // initialize
  memset(in_halo, 0, sizeof(char)*TOT_NUM_PIXELS);

  if (OPTIMIZE){
    forbidden = (char *) malloc(sizeof(char)*TOT_NUM_PIXELS);
    if (!forbidden){
      fprintf(stderr, "find_halos.c: Error allocating memory for forbidden box\nAborting...\n");
      fftwf_free(box);
      free(in_halo);
      return -1;
    }
  }

  // open k-space box to read in
  sprintf(filename, "../Boxes/deltak_z0.00_%i_%.0fMpc", DIM, BOX_LEN);
  IN = fopen(filename, "rb");
  if (!IN){
    fprintf(stderr, "find_halos.c: Unable to open file %s for reading\nAborting...\n", filename);
    fftwf_free(box);
    free(in_halo);
    return -1;
  }

  // remove conflicting files
  sprintf(filename, "../Output_files/FgtrM_files/hist_halos_z%.2f_%i_%.0fMpc_b%.3f_c%.3f", REDSHIFT, DIM, BOX_LEN, SHETH_b, SHETH_c);
  system(filename);
  sprintf(filename, "rm ../Output_files/DNDLNM_files/hist_halos_z%.2f_%i_%.0fMpc_b%.3f_c%.3f", REDSHIFT, DIM, BOX_LEN,  SHETH_b, SHETH_c);
  system(filename);

  // open the output files
  sprintf(filename, "../Output_files/Halo_lists/halos_z%.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
  OUT = fopen(filename, "w");
  if (!OUT){
    fprintf(stderr, "Unable to open file %s for writting!\n", filename);
    fftwf_free(box);
    free(in_halo);
    fclose(IN);
    return -1;
  }
  /************  END INITIALIZATION ****************************/

  // lets filter it now
  // set initial R value
  R = MtoR(M_MIN*1.01); // one percent higher for rounding
  while (R < L_FACTOR*BOX_LEN)
    R*=DELTA_R_FACTOR;
  fgrtm=dfgrtm=0;
  n=0;
  Delta_R = L_FACTOR*2*BOX_LEN/(DIM+0.0);
  while ((R > 0.5*Delta_R) && (RtoM(R) >= M_MIN)){ // filter until we get to half the pixel size or M_MIN
    M = RtoM(R);
    if (DELTA_CRIT_MODE == 1){ // use sheth tormen correction
      delta_crit = growth_factor*sheth_delc(Deltac/growth_factor, sigma_z0(M));
    }
    fprintf(stderr, "Now doing R=%.3f Mpc, M=%.2e...\n", R, M);
    fprintf(LOG, "Now doing R=%.3f Mpc, M=%.2e, clock=%.2f\n", R, M, (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);

    // first let's check if virialized halos of this size are rare enough 
    // that we don't have to worry about them (let's define 7 sigma away, as in Mesinger et al 05)
    if ((sigma_z0(M)*growth_factor*7) < delta_crit){
      fprintf(stderr, "Too rare! Skipping...\n");
      R /= DELTA_R_FACTOR;
      continue;
    }

    fprintf(LOG, "begin read, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    // read in the box
    rewind(IN);
    if (mod_fread(box, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS, 1, IN)!=1){
      fprintf(stderr, "find_halos.c: Read error occured!\n");
      fftwf_free(box);
      fclose(IN);
      fclose(OUT);
      free(in_halo);
      return -1;
    }
    fprintf(LOG, "end read, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);

    // now filter the box on scale R
    // 0 = top hat in real space, 1 = top hat in k space
    fprintf(LOG, "begin filter, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    filter(box, HALO_FILTER, R);
    fprintf(LOG, "end filter, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);

    // do the FFT to get delta_m box
    fprintf(LOG, "begin fft, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    plan = fftwf_plan_dft_c2r_3d(DIM, DIM, DIM, (fftwf_complex *)box, (float *)box, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    fprintf(LOG, "end fft, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);

    /*****************  BEGIN OPTIMIZATION *****************************/
    // to optimize speed, if the filter size is large (switch to collapse fraction criteria later)
    // 
    if (OPTIMIZE && (M > OPTIMIZE_MIN_MASS)){
      fprintf(LOG, "begin initialization of forbidden, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
      memset(forbidden, 0, sizeof(char)*TOT_NUM_PIXELS);
      // now go through the list of existing halos and paint on the no-go region onto <forbidden>
      sprintf(filename, "../Output_files/Halo_lists/halos_z%.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
      F = fopen(filename, "r");
	fscanf(F, "%e %f %f %f", &dummy, &x_temp, &y_temp, &z_temp);
      while (!feof(F)){
	R_temp = MtoR(dummy);
	update_in_halo(forbidden, R_temp+R_OVERLAP_FACTOR*R, rint(x_temp*DIM), rint(y_temp*DIM), rint(z_temp*DIM));
	fscanf(F, "%e %f %f %f", &dummy, &x_temp, &y_temp, &z_temp);
      }
      fclose(F);
      fprintf(LOG, "end initialization of forbidden, clock=%.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
    }
    /****************  END OPTIMIZATION  *******************************/

    // now lets scroll through the box, flagging all pixels with delta_m > delta_crit
    dn=0;
    for (x=0; x<DIM; x++){
      for (y=0; y<DIM; y++){
	for (z=0; z<DIM; z++){
	  delta_m = *((float *)box + R_FFT_INDEX(x,y,z)) * growth_factor / VOLUME;       // don't forget the factor of 1/VOLUME!
	  // if not within a larger halo, and radii don't overlap print out stats, and update in_halo box
	  /*********  BEGIN OPTIMIZATION **********/
	  if (OPTIMIZE && (M > OPTIMIZE_MIN_MASS)){
	    if ( (delta_m > delta_crit) && !forbidden[R_INDEX(x,y,z)]){
	    fprintf(stderr, "Found halo #%i, delta_m = %.3f at (x,y,z) = (%i,%i,%i)\n", n+1, delta_m, x,y,z);
	    fprintf(OUT, "%e\t%f\t%f\t%f\n", M, x/(DIM+0.0), y/(DIM+0.0), z/(DIM+0.0));
	    fflush(NULL);
	    update_in_halo(in_halo, R, x,y,z); // flag the pixels contained within this halo
	    update_in_halo(forbidden, (1+R_OVERLAP_FACTOR)*R, x,y,z); // flag the pixels contained within this halo
	    dn++; // keep track of the number of halos
	    n++;
	    }
	  }
	  /*********  END OPTIMIZATION **********/

	  else if ((delta_m > delta_crit) && !in_halo[R_INDEX(x,y,z)] && !overlap_halo(in_halo, R, x,y,z)){ // we found us a "new" halo!
	    fprintf(stderr, "Found halo #%i, delta_m = %.3f at (x,y,z) = (%i,%i,%i)\n", n+1, delta_m, x,y,z);
	    fprintf(OUT, "%e\t%f\t%f\t%f\n", M, x/(DIM+0.0), y/(DIM+0.0), z/(DIM+0.0));
	    fflush(NULL);
	    update_in_halo(in_halo, R, x,y,z); // flag the pixels contained within this halo
	    dn++; // keep track of the number of halos
	    n++;
	  }
	}
      }
    }

    if (dn > 0){
      // now lets print out the mass functions (FgrtR)
      fgrtm += M/(RHOcrit*OMm)*dn/VOLUME;
      dfgrtm += pow(M/(RHOcrit*OMm)*sqrt(dn)/VOLUME, 2);
      sprintf(filename, "../Output_files/FgtrM_files/hist_halos_z%.2f_%i_%.0fMpc_b%.3f_c%.3f", REDSHIFT, DIM, BOX_LEN, SHETH_b, SHETH_c);
      F = fopen(filename, "a");
      fprintf(F, "%e\t%e\t%e\t%e\t%e\t%e\n", M, fgrtm, sqrt(dfgrtm), FgtrM(REDSHIFT, M), FgtrM_st(REDSHIFT, M), FgtrM_bias(REDSHIFT, M, 0, sigma_z0(RtoM(L_FACTOR*BOX_LEN))) );
      fclose(F);

      // and the dndlnm files
      sprintf(filename, "../Output_files/DNDLNM_files/hist_halos_z%.2f_%i_%.0fMpc_b%.3f_c%.3f", REDSHIFT, DIM, BOX_LEN,  SHETH_b, SHETH_c);
      F = fopen(filename, "a");
      //dm = RtoM(DELTA_R_FACTOR*R)-M;
      dlnm = log(RtoM(DELTA_R_FACTOR*R)) - log(M);
      fprintf(F, "%e\t%e\t%e\t%e\t%e\t%e\n", M, dn/VOLUME/dlnm, sqrt(dn)/VOLUME/dlnm, M*dNdM(REDSHIFT, M), M*dNdM_st(REDSHIFT, M), M*dnbiasdM(M, REDSHIFT, RtoM(L_FACTOR*BOX_LEN), 0) );
      //      fprintf(F, "%e\t%e\t%e\t%e\t%e\t%e\n", M, M*dn/VOLUME/dm, M/dm/VOLUME*sqrt(dn), M*dNdM(REDSHIFT, M), M*dNdM_st(REDSHIFT, M), M*dnbiasdM(M, REDSHIFT, RtoM(BOX_LEN), 0) );
      fclose(F);
    }

    R /= DELTA_R_FACTOR;
  }


  // deallocate 
  fclose(OUT);
  fclose(IN);
  fclose(LOG);
  fftwf_free(box);

  /*
  // print in_halo box
  sprintf(filename, "../Boxes/in_halo_z%.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
  OUT = fopen(filename, "wb");
  fprintf(stderr, "Now writting in_halo box at %s\n", filename);
  if (mod_fwrite(in_halo, sizeof(char)*TOT_NUM_PIXELS, 1, OUT)!=1){
    fprintf(stderr, "find_halos.c: Write error occured while writting in_halo box.\n");
  }
  fclose(OUT);
  */

  free(in_halo);

  if (OPTIMIZE)
    free(forbidden);

  return 0;
}
