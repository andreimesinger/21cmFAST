#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"

/*
  USAGE: delta_ps <deltax filename> <output filename>

  box is assumed to be of the HII dimension defined in ANAL_PARAM.H
*/

#define FORMAT (int) 0 /* 0= unpadded binary box; 1= FFT padded binary box (outdated) */
#define CONVERT_TO_DELTA (int) 0 /* 1= convert the field to a zero-mean delta; 
		     be careful not to do this with fields which are already zero mean */

int main(int argc, char ** argv){
  char filename[100];
  FILE *F;
  float REDSHIFT;
  int x,y,z, format;
  fftwf_complex *deltax;
  fftwf_plan plan;
  float k_x, k_y, k_z, k_mag, k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor;
  int i,j,k, n_x, n_y, n_z, NUM_BINS;
  double dvdx, ave, new_ave, *p_box, *k_ave;
  unsigned long long ct, *in_bin_ct;

  // check arguments
  if (argc != 3){
    fprintf(stderr, "USAGE: delta_ps <deltax filename> <output filename>\nAborting\n");
    return -1;
  }
  // initialize and allocate thread info
  if (fftwf_init_threads()==0){
    fprintf(stderr, "init: ERROR: problem initializing fftwf threads\nAborting\n.");
    return -1;
  }
  fftwf_plan_with_nthreads(NUMCORES); // use all processors for init

  ave=0;
  //allocate and read-in the density array
  deltax = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!deltax){
    fprintf(stderr, "delta_T: Error allocating memory for deltax box\nAborting...\n");
    fftwf_cleanup_threads(); return -1;
  }
  F = fopen(argv[1], "rb");
  switch (FORMAT){
    // FFT format
  case 1:
    fprintf(stderr, "Reading in FFT padded deltax box\n");
    if (mod_fread(deltax, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS, 1, F)!=1){
      fftwf_free(deltax);
      fprintf(stderr, "deltax_ps.c: unable to read-in file\nAborting\n");
      fftwf_cleanup_threads(); return -1;
    }
    break;

    // unpaded format
  case 0:
    fprintf(stderr, "Reading in unpadded box\n");
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if (fread((float *)deltax + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	    fprintf(stderr, "init.c: Read error occured!\n");
	    fftwf_free(deltax);
	    fftwf_cleanup_threads(); return -1;	    
	  }
      	  ave += *((float *)deltax + HII_R_FFT_INDEX(i,j,k));
	}
      }
    }
    ave /= (double)HII_TOT_NUM_PIXELS;
    fprintf(stderr, "Average is %e\n", ave);
    break;

  default:
    fprintf(stderr, "Wrong format code\naborting...\n");
    fftwf_free(deltax);
    fftwf_cleanup_threads(); return -1;	    
  }
  fclose(F);


  if (CONVERT_TO_DELTA){
    new_ave = 0;
    fprintf(stderr, "Now converting field to zero-mean delta\n");
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  *((float *)deltax + HII_R_FFT_INDEX(i,j,k)) /= ave;
	  *((float *)deltax + HII_R_FFT_INDEX(i,j,k)) -= 1;
	  new_ave += *((float *)deltax + HII_R_FFT_INDEX(i,j,k));
	}
      }
    }
    new_ave /= (double) HII_TOT_NUM_PIXELS;
    fprintf(stderr, "The mean value of the field is now %e\n", new_ave);
  }



  // do the FFTs
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax, (fftwf_complex *)deltax, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
     deltax[ct] *= VOLUME/(HII_TOT_NUM_PIXELS+0.0);
  }


  /******  PRINT OUT THE POWERSPECTRUM  *********/

  k_factor = 1.4;
  k_first_bin_ceil = DELTA_K;
  k_max = DELTA_K*HII_DIM;
  // initialize arrays
   NUM_BINS = 0;
  k_floor = 0;
  k_ceil = k_first_bin_ceil;
  while (k_ceil < k_max){
    NUM_BINS++;
    k_floor=k_ceil;
    k_ceil*=k_factor;
  }

  p_box =  (double *)malloc(sizeof(double)*NUM_BINS);
  k_ave =  (double *)malloc(sizeof(double)*NUM_BINS);
  in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*NUM_BINS);
  if (!p_box || !in_bin_ct || !k_ave){ // a bit sloppy, but whatever..
    fprintf(stderr, "delta_T.c: Error allocating memory.\nAborting...\n");
    fftwf_free(deltax);
    fftwf_cleanup_threads(); return -1;
  }
  for (ct=0; ct<NUM_BINS; ct++){
    p_box[ct] = k_ave[ct] = 0;
    in_bin_ct[ct] = 0;
  }


  // now construct the power spectrum file
  for (n_x=0; n_x<HII_DIM; n_x++){
    if (n_x>HII_MIDDLE)
      k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
    else
      k_x = n_x * DELTA_K;

    for (n_y=0; n_y<HII_DIM; n_y++){
      if (n_y>HII_MIDDLE)
	k_y =(n_y-HII_DIM) * DELTA_K;
      else
	k_y = n_y * DELTA_K;

      for (n_z=0; n_z<=HII_MIDDLE; n_z++){ 
	k_z = n_z * DELTA_K;
	
	k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

	// now go through the k bins and update
	ct = 0;
	k_floor = 0;
	k_ceil = k_first_bin_ceil;
	while (k_ceil < k_max){
	  // check if we fal in this bin
	  if ((k_mag>=k_floor) && (k_mag < k_ceil)){
	    in_bin_ct[ct]++;
	    p_box[ct] +=  pow(k_mag,3)*pow(cabs(deltax[HII_C_INDEX(n_x, n_y, n_z)]), 2) / (2.0*PI*PI*VOLUME);
	    // note the 1/VOLUME factor, which turns this into a power density in k-space

	    k_ave[ct] += k_mag;
	    break;
	  }

	  ct++;
	  k_floor=k_ceil;
	  k_ceil*=k_factor;
	}
      }
    }
  } // end looping through k box
  fftwf_free(deltax);

  // now lets print out the k bins
  F = fopen(argv[2], "w");
  if (!F){
    fprintf(stderr, "delta_T.c: Couldn't open file %s for writting!\n", filename);
    fftwf_cleanup_threads(); return -1;
  }
  for (ct=1; ct<NUM_BINS; ct++){
    fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
  }
  fclose(F);

  /****** END POWER SPECTRUM STUFF   ************/

  free(p_box); free(k_ave); free(in_bin_ct);

  fftwf_cleanup_threads(); return 0;
}
