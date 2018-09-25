#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"

/*
  USAGE: boxcar_smooth_field <0=no fft padding; 1=fft padding> <input highres filename> <output lowres filename>

  Date 31.12.2009
  Author:  Andrei Mesinger
*/

int main(int argc, char ** argv){
  fftwf_complex *box, *smoothed_box;
  fftwf_plan plan;
  unsigned long long ct;
  int format, pixel_factor,i,j,k;
  float mass_factor;
  FILE *F;

  if (argc != 4){
    fprintf(stderr, "USAGE: boxcar_smooth_field <0=no fft padding; 1=fft padding> <input highres filename> <output lowres filename>\nAborting...\n");
    return 0;
  }
  format = atoi(argv[1]);
  // find factor of HII pixel size / deltax pixel size
  pixel_factor = DIM/HII_DIM +0.5;
  mass_factor = pow(pixel_factor, 3);

  // allocate memory
  box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);
  smoothed_box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!box || !smoothed_box){
    fprintf(stderr, "smooth_field.c: Error allocating memory for box.\nAborting...\n");
    return -1;
  }
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
    smoothed_box[ct]=0;
  }

  // open file and read-in
  F=fopen(argv[2], "rb");
  if (!F){
    fprintf(stderr, "smooth_field.c: Error open binary file %s for reading\nAborting...\n", argv[2]);
    fftwf_free(box);
    return -1;
  }
  if (format==0){ // box has no fft padding
    for (i=0; i<DIM; i++){
      for (j=0; j<DIM; j++){
	for (k=0; k<DIM; k++){
	  if (fread((float *)box + R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	    fprintf(stderr, "smooth_field.c: Error reading-in binary file %s\nAborting...\n", argv[2]);
	    fftwf_free(box), fclose(F);
	    return -1;
	  }
	}
      }
    }
  }
  else if (format==1){ // box has fft padding
    if (mod_fread(box, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "smooth_field.c: Error reading-in binary file %s\nAborting...\n", argv[2]);
      fftwf_free(box), fclose(F);
      return -1;
    }
  }
  else{
    fprintf(stderr, "smooth_field.c: Incorrect format specifier %i\nAborting...\n", format);
    fftwf_free(box), fclose(F);
    return -1;
  }
  fclose(F);


  // go through the high-res box, mapping the mass onto the low-res (updated) box
  for (i=0; i<DIM;i++){
    for (j=0; j<DIM;j++){
      for (k=0; k<DIM;k++){

	*( (float *) smoothed_box + HII_R_FFT_INDEX(i/pixel_factor, j/pixel_factor, k/pixel_factor)) +=
	  *( (float *) box + R_FFT_INDEX(i,j,k))/mass_factor;
      }
    }
  }
  

  // now sample and print to file
  F=fopen(argv[3], "wb");
  if (!F){
    fprintf(stderr, "smooth_field.c: Error open binary file %s for writting\nAborting...\n", argv[3]);
    fftwf_free(box);
    return -1;
  }
  if (format==0){ // no fft padding
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if( fwrite( (float *)smoothed_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	    fprintf(stderr, "smooth_field.c: Error writting binary file %s\nAborting...\n", argv[3]);
	    fftwf_free(box), fclose(F);
	    return -1;
	  }
	}
      }
    }
  }
  else{ //fft padding
   if (mod_fwrite(smoothed_box, sizeof(float)*HII_TOT_FFT_NUM_PIXELS, 1, F)!=1){
     fprintf(stderr, "smooth_field.c: Error writting binary file %s\n", argv[3]);
   }
  }

   fftwf_free(smoothed_box); fftwf_free(box), fclose(F);
  return 0;
}
