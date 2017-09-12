#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "filter.c"

#define FLIP_BOX (int) (0)

/*****************************************************************************************
   USAGE: redshift_interpolate_boxes <BOX TYPE> 
           <filename containing list of boxes to be interpolated in increasing redshift order>

   BOX TYPE is type int and can take the following values:
     0 = boxes without FFT padding such as xH and deltaT boxes
     1 = OUTDATED boxes with FFT padding, such as older version overdensity delta box
     2 = velocity boxes WITHOUT FFT PADDING

   for the velocity field, the filename list (argv[2]) should just list one velocity component
   (vx, vy, or vz); it doesn't matter which one; the code will open the appropriate box.

   Outputed boxes are always oriented along the z-axis! They will be placed (without FFT
   padding) in ../Boxes, begging with the same prefix as the input and ending with "lighttravel"

   Date: 7.4.2011
   Author: Andrei Mesinger
 *****************************************************************************************/

FILE *LOG;

int read_box(char *box_filename, fftwf_complex *box, int format, int LOS_direction){
  int i,j,k;
  char *token, input_filename_prefix[300], input_filename_sufix[300];
  FILE *F;
  //fprintf(stderr, "%s\n", box_filename);
  // and read-in the first box
  if (format == 2){ // velocity box; we need to decide which component to read
    token = strtok(box_filename, "v");
    strcpy(input_filename_prefix, token);
    token = strtok(NULL, "_");
    token = strtok(NULL, "v");
    strcpy(input_filename_sufix, token);
    if (LOS_direction == 0){
      sprintf(box_filename, "%svx_%s", input_filename_prefix, input_filename_sufix);
    }
    else if (LOS_direction == 1){
      sprintf(box_filename, "%svy_%s", input_filename_prefix, input_filename_sufix);
    }
    else{
      sprintf(box_filename, "%svz_%s", input_filename_prefix, input_filename_sufix);
    }
  }
  fprintf(stderr, "Reading-in box: %s\n", box_filename);
  fprintf(LOG, "Reading-in box: %s\n", box_filename);
  if (!(F = fopen(box_filename, "rb"))){
    fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to open %s.\nAborting.\n", box_filename);
    fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to open %s.\nAborting.\n", box_filename);
    fclose(LOG);
    return -1;
  }

  // read in array
  if (format == 1){ // box has fft padding
    fprintf(stderr, "redshift_interpolate_boxes: WARNING: you should not be using box format code 1 for >=v1.1, since boxes are outputed without FFT padding.\n");
    if (mod_fread(box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS, 1, F)!=1)
      return -1;
  }
  else { // box has no fft padding
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if (fread((float *)box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	    return -1;
	  }
	}
      }
    }
  }

    
  fclose(F);
  fprintf(stderr, "Finished read\n");
  fprintf(LOG, "Finished read\n");

  return 0;
}

void copy_slice(fftwf_complex *box_interpolate, fftwf_complex *box_z1, fftwf_complex *box_z2, 
		double z, double z1, double z2, int slice_ct, int LOS_direction){
  int i,j,k;
  float fz1, fz2, fz;

  if (LOS_direction == 0){ // LOS is x - axis
    for (j=0;j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
	fz1 = *((float *)box_z1 + HII_R_FFT_INDEX(slice_ct, j,k));
	fz2 = *((float *)box_z2 + HII_R_FFT_INDEX(slice_ct, j,k));
	fz = (fz2 - fz1) / (z2 - z1) * (z - z1) + fz1; // linearly interpolate in z (time actually)
	*((float *)box_interpolate + HII_R_FFT_INDEX(j, k, slice_ct)) = fz;
	//	if ((j==10) && (k==10))
	//fprintf(stderr, "slice %i: interpolating between f(z1=%f)=%f and f(z2=%f)=%f to get f(z=%f)=%f\n", slice_ct, z1, fz1, z2, fz2, z, fz);
      }
    }
  }
  else if (LOS_direction == 1){ // LOS is y - axis
    for (i=0;i<HII_DIM; i++){
      for (k=0; k<HII_DIM; k++){
	fz1 = *((float *)box_z1 + HII_R_FFT_INDEX(i, slice_ct, k));
	fz2 = *((float *)box_z2 + HII_R_FFT_INDEX(i, slice_ct, k));
	fz = (fz2 - fz1) / (z2 - z1) * (z - z1) + fz1; // linearly interpolate in z
	*((float *)box_interpolate + HII_R_FFT_INDEX(i, k, slice_ct)) = fz;
      }
    }
  }
  else{
    for (i=0;i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	fz1 = *((float *)box_z1 + HII_R_FFT_INDEX(i, j, slice_ct));
	fz2 = *((float *)box_z2 + HII_R_FFT_INDEX(i, j, slice_ct));
	fz = (fz2 - fz1) / (z2 - z1) * (z - z1) + fz1; // linearly interpolate in z
	*((float *)box_interpolate + HII_R_FFT_INDEX(i, j, slice_ct)) = fz;
      }
    }
  }
}


int main(int argc, char ** argv){
  char box_filename[300], box_filename_z1[300], box_filename_z2[300], output_filename_prefix[300], output_filename_sufix[300], output_filename[300];
  char input_filename_prefix[300], input_filename_sufix[300];
  FILE *BOX_LIST, *F;
  fftwf_complex *box_z1, *box_z2, *box_interpolate; 
  int format, LOS_direction, slice_ct;
  double z1, z2, z, dR;
  float start_z, end_z;
  int i,j,k;
  char *token;
  unsigned long long ct;

  /*
   char line[200];
  start_z = 5.986;
  sprintf(line, "%06.2f", start_z);
  printf("%s\n", line);
  printf("%f\n", atof(line));


   sprintf(line, "../Boxes/updated_vz_z5.60_450_500Mpc");//../Boxes/xH_nohalos_nf0.000000_eff31.5_HIIfilter1_Mmin1.0e+08_RHIImax30_z5.60_450_500Mpc");

  token = strtok(line, "v");
  strcpy(input_filename_prefix, token);
  token = strtok(NULL, "_");
  token = strtok(NULL, "v");
  strcpy(input_filename_sufix, token);
  printf("%s\n%s\n", input_filename_prefix, input_filename_sufix);
  return 0;
  while (!isdigit(token[0]))
    token = strtok(NULL, "z");
  z1 = atof(strtok(token, "_"));
  printf("%f\n", z1);
  return 0;
  */

  if (argc != 3){
    fprintf(stderr, "USAGE: redshift_interpolate_boxes <box type> <filename containing list of boxes to be interpolated in increasing redshift order>\nAborting\n");
    return -1;
  }
  format = atoi(argv[1]);
  LOS_direction = 2; // 0 = x-axis, 1= y-axis, 2= z-axis; starting value if flip boxes is initialized
  dR = (BOX_LEN / (double) HII_DIM) * CMperMPC; // size of cell (in comoving cm)

  // open the box list and log files
  if (!(BOX_LIST = fopen(argv[2], "r"))){
    fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to open box filelist %s.\nAborting.\n", argv[2]);
    return -1;
  }
  system("mkdir ../Log_files");
  if (!(LOG = fopen("../Log_files/redshift_interpolate_boxes_log", "w"))){
    fprintf(stderr, "WARNING: redshift_interpolate_boxes: Unable to open log file\n");
  }

  // allocate memory for the three boxes
  if (!(box_z1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to allocate memory\nAborting\n");
    fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to allocate memory\nAborting\n");
    fclose(LOG); fclose(BOX_LIST); fclose(F);
  }
  if (!(box_z2 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to allocate memory\nAborting\n");
    fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to allocate memory\nAborting\n");
    fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fclose(F);
  }
  if (!(box_interpolate = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to allocate memory\nAborting\n");
    fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to allocate memory\nAborting\n");
    fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); fclose(F);
  }

  // read in the first box
  fscanf(BOX_LIST, "%s\n", box_filename);
  strcpy(box_filename_z1, box_filename);
  // and read-in the first box
  if (read_box(box_filename, box_z1, format, LOS_direction) != 0){
    fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename);
    fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename);
    fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); return -1;
  }

  // and get the redshift from the filename
  strcpy(box_filename, box_filename_z1);
  token = strtok(box_filename, "z");
  strcpy(output_filename_prefix, token);
  while (!isdigit(token[0]))
    token = strtok(NULL, "z");
  z1 = atof(strtok(token, "_"));
  fprintf(stderr, "Output filename prefix is %s.\nRead-in box at redshift %f\n", output_filename_prefix, z1);
  fprintf(LOG, "Output filename prefix is %s.\nRead-in box at redshift %f\n", output_filename_prefix, z1);
  /***************************  END INITIALIZATIONS  *****************************************************/


  // scroll through all of the boxes
  z = start_z = z1;
  slice_ct = 0;
  while (!feof(BOX_LIST)){
    // read in the next box
    fscanf(BOX_LIST, "%s\n", box_filename);
    strcpy(box_filename_z2, box_filename);
    if (read_box(box_filename, box_z2, format, LOS_direction) != 0){
      fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename);
      fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename);
      fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); fftwf_free(box_interpolate); fclose(F);
      return -1;
    }
    // and get the redshift from the filename
    strcpy(box_filename, box_filename_z2);
    token = strtok(box_filename, "z");
    while (!isdigit(token[0]))
      token = strtok(NULL, "z");
    z2 = atof(strtok(token, "_"));
    fprintf(stderr, "Read-in box at redshift %f\nInterpolating in the direction of %i\n", z2, LOS_direction);
    fprintf(LOG, "Read-in box at redshift %f\nInterpolating in the direction of %i\n", z2, LOS_direction);

    // now do the interpolation
    while (z < z2){ // until we move to the next set of boxes
 
     // check if we filled-up our array and write-out
      if (slice_ct == HII_DIM){
	end_z = z;
	// smooth the field?
	// maybe?
	// nah...

	// write out the box
	sprintf(output_filename, "%s_zstart%09.5f_zend%09.5f_FLIPBOXES%i_%i_%.0fMpc_lighttravel", 
		output_filename_prefix, start_z, end_z, FLIP_BOX, HII_DIM, BOX_LEN);
	if (!(F=fopen(output_filename, "wb"))){
	  fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to open file %s.\nAborting\n", output_filename);
	  fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to open file %s.\nAborting\n", output_filename);
	  fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); fftwf_free(box_interpolate);
	  return -1;
	}
	// don't include FFT padding
	for (i=0; i<HII_DIM; i++){
	  for (j=0; j<HII_DIM; j++){
	    for (k=0; k<HII_DIM; k++){

	      if (fwrite((float *)box_interpolate + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F) != 1){
		fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to open file %s.\nAborting\n", output_filename);
		fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to open file %s.\nAborting\n", output_filename);
		fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); fftwf_free(box_interpolate); fclose(F);
		return -1;
	      }
	      
	    }
	  }
	}
	fclose(F);

	fprintf(stderr, "Written light travel box at %s.\n", output_filename);
	fprintf(LOG, "Written light travel box at %s.\n", output_filename);

	// update quantities
	slice_ct=0;
	start_z = end_z;
	memset(box_interpolate, 0, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); //not needed but easier to debug
	if (FLIP_BOX){
	  LOS_direction = ++LOS_direction % 3;
	  if (format == 2){  // we are rotating the cube, so we need to read-in the appropriate velocity components
	    if (read_box(box_filename_z1, box_z1, format, LOS_direction) != 0){
	      fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename_z1);
	      fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename_z1);
	      fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); fftwf_free(box_interpolate); fclose(F);
	      return -1;
	    }
	    if (read_box(box_filename_z2, box_z2, format, LOS_direction) != 0){
	      fprintf(stderr, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename_z2);
	      fprintf(LOG, "ERROR: redshift_interpolate_boxes: Unable to read-in box %s\nAborting\n", box_filename_z2);
	      fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); fftwf_free(box_interpolate); fclose(F);
	      return -1;
	    }
	  }
	}
      } // we are now continuing with a new interpolation box

      copy_slice(box_interpolate, box_z1, box_z2, gettime(z), gettime(z1), gettime(z2), slice_ct, LOS_direction);
      // note the interpolation in copy_slice is done linearly in time using the "gettime" function
      slice_ct++;
      z -= dR / drdz(z);
    } // done with this pair of boxes, moving on to the next redshift pair

    // copy the z2 box now to z1, since it is the lower bound, and read the next filename at z2
    memcpy(box_z1, box_z2, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    z1 = z2;
    strcpy(box_filename_z1, box_filename_z2);
    fprintf(stderr, "\n");
    fprintf(LOG, "\n");
  }
  

  fclose(LOG); fclose(BOX_LIST); fftwf_free(box_z1); fftwf_free(box_z2); fftwf_free(box_interpolate);
  return 0;
}
