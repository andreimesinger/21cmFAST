#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"

#define NUM_TRIALS (unsigned long long) 1e7 // num of times we will perform the test
#define dr (float)(1.0/(HII_DIM*5.0))
#define SIZE_BINSTEP (float) 2*BOX_LEN/(HII_DIM+0.0) // single step in histogram construction
#define MAX_MAX (float) 1000.0 // truncate distributions after this value
#define VOXEL_NF_CUTOFF (float) 0.5 // value demarkating "neutral" and "ionized" voxels
/*

  USAGE: gen_size_distr <REDSHIFT> <REGION> <IN_BUBBLE BOX filename>

  PROGRAM GEN_SIZE_DISTR generates size distributions by randomly choosing a pixel and randomly
  generating a directional vector, then distance along that directional vector to a phase transition

  <REGION> is an integer:
     0 = genarate size distributions of ionized bubbles
     1 = genarate size distributions of neutral regions

*/


int main (int argc, char ** argv){
  char *in_bubble, filename[100];
  FILE *F, *LOG;
  float REDSHIFT, r_mag, r_x, r_y, r_z, bin_floor, R, MAX, dpdR, P, *dist, nf, xH;
  int REGION_FLAG, x_0,y_0,z_0, x_curr, y_curr,z_curr, wrap,j,k;
  unsigned long long i, bin_ct;
  gsl_rng * r;

  /***************   BEGIN INITIALIZATION   **************************/

  // check usage
  if (argc != 4){
    fprintf(stderr, "USAGE: gen_size_distr <REDSHIFT> <REGION> <IN_BUBBLE BOX filename>\nAborting...\n");
    return -1;
  }
  REDSHIFT = atof(argv[1]);
  REGION_FLAG = atoi(argv[2]);
  if ((REGION_FLAG != 1) && (REGION_FLAG != 0)){
    fprintf(stderr, "USAGE: gen_size_distr <REDSHIFT> <REGION> <IN_BUBBLE BOX filename>\nAborting...\n");
    return -1;
  }
  // allocate
  dist = (float *) malloc(sizeof(float)*NUM_TRIALS);
  if (!dist){
    fprintf(stderr, "gen_size_distr: Unable to allocate enough memory\nAborting\n");
    return -1;
  }

  // seed the random number generator
  //  gaus_seed(SIZE_RANDOM_SEED);
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, SIZE_RANDOM_SEED);


  // read in the bubble box
  // allocate memory for the boolean in_bubble
  in_bubble = (char *) malloc(sizeof(char)*HII_TOT_NUM_PIXELS);
  if (!in_bubble){
    fprintf(stderr, "gen_size_distr: Error allocating memory for in_bubble box\nAborting...\n");
    gsl_rng_free (r); return -1;
  }
  F = fopen(argv[3], "rb");
  if (!F){
    fprintf(stderr, "gen_size_distr: Error opening file %s for reading\nAborting...\n", argv[3]);
    free(in_bubble);
    gsl_rng_free (r); return -1;
  }
  fprintf(stderr, "Reading in xH box\n");
  nf = 0;
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
        if (fread(&xH, sizeof(float), 1, F)!=1){
          fprintf(stderr, "delta_T: Read error occured while reading neutral_fraction box.\n");
	  fclose(F); free(in_bubble);
	  gsl_rng_free (r); return -1;
        }

	nf += xH;
	if (xH < VOXEL_NF_CUTOFF)
	  in_bubble[HII_R_INDEX(i,j,k)] = 1;
	else
	  in_bubble[HII_R_INDEX(i,j,k)] = 0;    
      }
    }
  }
  fclose(F);
  nf /= (double)HII_TOT_NUM_PIXELS;

  // check if the ionization field is fully neutral or ionized. if so calling this function is retarded
  if ((nf<FRACT_FLOAT_ERR) || (nf>(1-FRACT_FLOAT_ERR))){
    fprintf(stderr, "gen_size_distr: The ionization field is only a single phase.  Aborting gen_size_distr.\n");
    return 0;
  }


  // open log file
  LOG = fopen("../Log_files/gen_size_log_file", "w");
  if (!LOG){
    fprintf(stderr, "gen_size_distr: Error opening log file\nAborting...\n");
    free(in_bubble);
    gsl_rng_free (r); return -1;
  }

  // open output file
  if (REGION_FLAG){
    if (USE_HALO_FIELD)
      sprintf(filename, "../Output_files/Size_distributions/neutral_nf%f_z%06.2f_%i_%.0fMpc", nf, REDSHIFT, HII_DIM, BOX_LEN);
    else
      sprintf(filename, "../Output_files/Size_distributions/neutral_no_halos_nf%f_z%06.2f_%i_%.0fMpc", nf, REDSHIFT, HII_DIM, BOX_LEN);
  }
  else{
    if (USE_HALO_FIELD)
      sprintf(filename, "../Output_files/Size_distributions/ionized_nf%f_z%06.2f_%i_%.0fMpc", nf, REDSHIFT, HII_DIM, BOX_LEN);
    else
      sprintf(filename, "../Output_files/Size_distributions/ionized_no_halos_nf%f_z%06.2f_%i_%.0fMpc", nf, REDSHIFT, HII_DIM, BOX_LEN);
  }
  F=fopen(filename, "w");
  if (!F){
    fprintf(stderr, "gen_size_distr: Error opening output file\nAborting...\n");
    free(in_bubble);
    fclose(LOG);
    gsl_rng_free (r); return -1;
  }
  MAX = 0;
  /***************   END INITIALIZATION   **************************/


  // go through the trials
  for (i=0; i<NUM_TRIALS; i++){

    // find appropriate pixel
    x_0 = gsl_rng_uniform (r)*HII_DIM;
    y_0 = gsl_rng_uniform (r)*HII_DIM;
    z_0 = gsl_rng_uniform (r)*HII_DIM;
    if ( ((REGION_FLAG==0) && !in_bubble[HII_R_INDEX(x_0, y_0, z_0)] ) ||
	 ((REGION_FLAG==1) && in_bubble[HII_R_INDEX(x_0, y_0, z_0)] ) ){
      // we are not in the desired phase, try again
      i--;
      continue;
    }


    // construct unit directional vector
    r_x = gsl_rng_uniform (r)-0.5;
    r_y = gsl_rng_uniform (r)-0.5;
    r_z = gsl_rng_uniform (r)-0.5;
    r_mag = sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
    r_x /= r_mag; r_y /= r_mag; r_z /= r_mag;
    //    fprintf(stderr, "Found %llu at (%i, %i, %i), going in direction (%.3f, %.3f, %.3f)\n", i, x_0, y_0, z_0, r_x, r_y, r_z);
    //fprintf(LOG, "Found %llu at (%i, %i, %i), going in direction (%.3f, %.3f, %.3f)\n", i, x_0, y_0, z_0, r_x, r_y, r_z);
    wrap = 0;     // flag if we wrapped around
    r_mag = 0;
    // now follow the directional vector until we reach a phase transition
    while (1){
      r_mag += dr;
      x_curr = (r_mag * r_x) * HII_DIM  + x_0+0.5;
      y_curr = (r_mag * r_y) * HII_DIM  + y_0+0.5;
      z_curr = (r_mag * r_z) * HII_DIM  + z_0+0.5;

      // check if we stepped out of the box
      while (x_curr >= HII_DIM){ x_curr -= HII_DIM; wrap=1;}
      while (x_curr < 0){ x_curr += HII_DIM; wrap=1;}
      while (y_curr >= HII_DIM){ y_curr -= HII_DIM; wrap=1;}
      while (y_curr < 0){ y_curr += HII_DIM; wrap=1;}
      while (z_curr >= HII_DIM){ z_curr -= HII_DIM; wrap=1;}
      while (z_curr < 0){ z_curr += HII_DIM; wrap=1;}

      // check if we crossed a boundary
      if ( in_bubble[HII_R_INDEX(x_0, y_0, z_0)] != in_bubble[HII_R_INDEX(x_curr, y_curr, z_curr)] ){ // got it!
	//	fprintf(stderr, "Got it at d=%.3f Mpc, (%i, %i, %i)\n", r_mag*BOX_LEN, x_curr, y_curr, z_curr);
	//fprintf(LOG, "Got it at d=%.3f Mpc, (%i, %i, %i)\n", r_mag*BOX_LEN, x_curr, y_curr, z_curr);
       	dist[i] = r_mag*BOX_LEN;
	if (dist[i] > MAX)
	  MAX = dist[i];
	break;
      }
      else if (wrap && (x_0==x_curr) && (y_0==y_curr) && (z_0==z_curr)){ // went around and back to first pixel
	// try again
	fprintf(stderr, "We wrapped around to the start pixel\n");
	fprintf(LOG, "We wrapped around to the start pixel\n");
	i--;
	break;
      }
    } // end while
  } // end for
  fprintf(stderr, "\nDone with Monte-Carlo\n\nNow printing histogram\n");
  fprintf(LOG, "\nDone with Monte-Carlo\n\nNow printing histogram\n");


  // print out the histogram
  // do more efficientl if we need super long trials
  bin_floor = 0;
  P=0;
  if (MAX > MAX_MAX){ MAX = MAX_MAX;}
  while (bin_floor < (MAX+SIZE_BINSTEP)){ // go through each bin
    // cycle through all the entries counting the ones that fall in the bin
    bin_ct = 0;
    for (i=0; i<NUM_TRIALS; i++){
      if ( (dist[i] >= bin_floor) && (dist[i] < (bin_floor+SIZE_BINSTEP)) ){ bin_ct++;}
    }

    // now print our results: mean value of bin followed by num of files that fall in bin
    R =  bin_floor+0.5*SIZE_BINSTEP;
    dpdR = (float)bin_ct/(float)NUM_TRIALS;
    dpdR /= SIZE_BINSTEP;
    //    printf("%f\t%f\t%f\t%u\n", dpdR, (float)bin_ct/(float)NUM_TRIALS/SIZE_BINSTEP, SIZE_BINSTEP, bin_ct);
    fprintf(F, "%f\t%e\n", R, R*dpdR );
    P+=dpdR;
    bin_floor += SIZE_BINSTEP;
  }
  fprintf(stderr, "Done %f!\n", P*SIZE_BINSTEP);
  fprintf(LOG, "Done!\n");


  // deallocate
  free(in_bubble);
  free(dist);
  fclose(LOG);
  fclose(F);

  gsl_rng_free (r); return 0;
}
