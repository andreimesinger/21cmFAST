#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "bubble_helper_progs.c"

#define   bin_factor (float) 1.2
#define   bin_first_bin_ceil (float) 0.001
#define   bin_max (float) 1000

/*
  Function FILTER_DEN_HIST prints out a histogram of the overdensities in the

  USAGE: filter_den_hist <R> <1=old format with FFT padding, 0= new format without padding> <deltax filename> <output filename>

  if R <= 0, the filtering step is skipped

  NOTE:  the relevant box parameters are taken from INIT_PARAMS.H
*/

// returns gaussian integrand
float dgausdx(float x, float mean, float sigma){
  float ex;

  ex = -0.5*pow((x-mean)/sigma, 2);
  return 1.0/(sigma*sqrt(TWOPI)) * pow(E, ex);
}

int main(int argc, char ** argv){
  fftwf_plan plan;
  fftwf_complex *box;
  float M, *delta_m, floor, ciel, del, R;
  double MAX, MIN;
  unsigned long long index, ct, *in_bin_ct;
  FILE *F;
  int x,y,z, format, i,j,k, NUM_BINS;
  char filename[500];
  double ave, *p_box, *bin_ave, value;
  float  bin_floor, bin_ceil;


  /********************* INITIALIZATION **************************/
  if (argc != 5){
    fprintf(stderr, "USAGE: filter_den_hist <R> <1=my format, 0=hy format> <deltax filename> <output filename>\nAborring...\n");
    return 0;
  }
  R = atof(argv[1]);
  format = atoi(argv[2]);

  // initialize power spectrum crap
  init_ps();

  /*
  printf("%f\n", sigma_z0(1e8));
  printf("%f\n", sigma_z0(1e9));
  printf("%f\n", sigma_z0(1e10));
  printf("%f\n", sigma_z0(1e11));
  printf("%f\n", sigma_z0(1e12));
  printf("%f\n", sigma_z0(1e15));
  printf("%f\n", MtoR(1e15));
  printf("%f\n", MtoR(1e16));
  return 0;
  */

  //allocate and read-in the density array
  box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!box){
    fprintf(stderr, "delta_T: Error allocating memory for box box\nAborting...\n");
    return -1;
  }
  F = fopen(argv[3], "rb");
  fprintf(stderr, "Reading in box box of HII_DIM=%i\n", HII_DIM);
  switch (format){
    // my format
  case 1:
    if (mod_fread(box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS, 1, F)!=1){
      fftwf_free(box); fclose(F);
      fprintf(stderr, "box_ps.c: unable to read-in file\nAborting\n");
      return -1;
    }
    
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  *((float *)box + HII_R_FFT_INDEX(i,j,k)) += 1;
	  if (*((float *)box + HII_R_FFT_INDEX(i,j,k)) < 0){
	    fprintf(stderr, "Less than 0???\n");
	  }
	}
      }
    }
    
    //convert \delta -> \Delta
    //    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){ box[ct] += 1;}
    break;

    // Hy format
  case 0:
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if (fread((float *)box + HII_R_FFT_INDEX(k,j,i), sizeof(float), 1, F)!=1){
	    fprintf(stderr, "init.c: Read error occured!\n");
	    fftwf_free(box); fclose(F);
	    return -1;	    
	  }
	  	  *((float *)box + HII_R_FFT_INDEX(k,j,i)) += 1; // convert to Ddelta

	  //	  *((float *)box + HII_R_FFT_INDEX(k,j,i)) *= *((float *)box + HII_R_FFT_INDEX(k,j,i));

	  if (*((float *)box + HII_R_FFT_INDEX(k,j,i)) < 0){
	    fprintf(stderr, "Less than 0???, %e\n", *((float *)box + HII_R_FFT_INDEX(k,j,i)));
	  }

	}
      }
    }
    break;

  default:
    fprintf(stderr, "Wrong format code\naborting...\n");
    fftwf_free(box);
    return -1;	    
  }
  fclose(F);


  // output file
  F = fopen(argv[4], "w");
  if (!F){
    fprintf(stderr, "Unable to open file: %s for writing\nAborting...\n", argv[4]);
    fftwf_free(box);
    return -1;
  }

  // allocate memory for the filtered histogram
  delta_m = (float *)malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!delta_m){
    fprintf(stderr, "Filter_den_hist.c: Error allocating memory for histogram array\nAborting...\n");
    fftwf_free(box); fclose(F);
    return -1;
  }

  /************  END INITIALIZATION ******************/

  if (R>0){
    //convert to k-space to filter
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)box, (fftwf_complex *)box, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
    //  real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
      box[ct] /= (HII_TOT_NUM_PIXELS+0.0);
    }

    // filter the field
    fprintf(stderr, "Read done, now filtering on scale R=%e Mpc\n", R);
    HII_filter(box, 0, R);
  
    // do the FFT to get delta_m box
    fprintf(stderr, "begin fft\n");
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)box, (float *)box, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    fprintf(stderr, "end fft\n");
    /*
      for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
      if (*((float *)box + HII_R_FFT_INDEX(k,j,i)) < 0){
      fprintf(stderr, "Less than 0???\n");
      }
      
      }
      }
      }
    */
  } // end filtering

    // copy the delta values over to the histogram array
  MIN = 1000;
  MAX = -1000;
  for (x=0; x<HII_DIM; x++){
    for (y=0; y<HII_DIM; y++){
      for (z=0; z<HII_DIM; z++){
	// get the filtered mass on this scale
	//delta_m[HII_R_INDEX(x,y,z)] = pow(*((float *)box + HII_R_FFT_INDEX(x,y,z)), 2);
		delta_m[HII_R_INDEX(x,y,z)] = *((float *)box + HII_R_FFT_INDEX(x,y,z));
	if (delta_m[HII_R_INDEX(x,y,z)] < MIN){ MIN = delta_m[HII_R_INDEX(x,y,z)];}
	if (delta_m[HII_R_INDEX(x,y,z)] > MAX){ MAX = delta_m[HII_R_INDEX(x,y,z)];}
      }
    }
  }
  fprintf(stderr, "Min = %.2f, Max = %.2f\n", MIN, MAX);

  // ghetto counting (lookup how to do logs of arbitrary bases in c...)
  NUM_BINS = 0;
  bin_floor = 0;
  bin_ceil = bin_first_bin_ceil;
  while (bin_ceil < bin_max){
    NUM_BINS++;
    bin_floor=bin_ceil;
    bin_ceil*=bin_factor;
  }

  p_box =  (double *)malloc(sizeof(double)*NUM_BINS);
  bin_ave =  (double *)malloc(sizeof(double)*NUM_BINS);
  in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*NUM_BINS);
  if (!p_box || !in_bin_ct || !bin_ave){ // a bit sloppy, but whatever..
    fprintf(stderr, "delta_T.c: Error allocating memory.\nAborting...\n");
    fftwf_free(box);
    return -1;
  }
  for (ct=0; ct<NUM_BINS; ct++){
    p_box[ct] = bin_ave[ct] = 0;
    in_bin_ct[ct] = 0;
  }

  for (x=0; x<HII_DIM; x++){
    for (y=0; y<HII_DIM; y++){
      for (z=0; z<HII_DIM; z++){
	value = delta_m[HII_R_INDEX(x,y,z)];
	if (value < 0)
	  value = 0;
	// now go through the k bins and update
	ct = 0;
	bin_floor = 0;
	bin_ceil = bin_first_bin_ceil;
	while (bin_ceil < bin_max){
	  // check if we fal in this bin
	  if ((value>=bin_floor) && (value < bin_ceil)){
	    in_bin_ct[ct]++;
	    bin_ave[ct] += value;
	    break;
	  }

	  ct++;
	  bin_floor=bin_ceil;
	  bin_ceil*=bin_factor;
	}
      }
    }
  } // end loop through box


  // now print
  for (ct=0; ct<NUM_BINS; ct++){
    value = bin_ave[ct]/(double)in_bin_ct[ct];
    fprintf(F, "%e\t%e\n", value, in_bin_ct[ct]/(double)HII_TOT_NUM_PIXELS/(bin_factor-1.0));
  }
 
  // deallocate memory
  fclose(F);
  free(delta_m);
  fftwf_free(box);
  free(p_box); free(bin_ave); free(in_bin_ct);
}
