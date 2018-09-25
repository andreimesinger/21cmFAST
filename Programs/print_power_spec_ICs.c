#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"


/*
  USAGE: print_power_spec <kbox_file_name> <output filename> [<k_max> <k_first_bin_ceil> <k_factor>]

  Function PRINT_POWER_SPEC reads in a delta_k box from <kbox_file_name>,
  and outputs:
     k-magnitude   <power_spectrum>   (1/V)<|delta_k|^2>
  to file with name <output filename>.
  
  The bins in <OUT> are generated using <k_max> and <k_factor>, and
  <|delta_k|^2> is the average value in a particular bin from <box>,
  whereas <power_spectrum> is the average input power spectrum
  (see COSMOLOGY.H)

  <k_max>, <k_first_bin_ceil>, and <k_factor> are optional and default to:
  DIM*TWOPI/BOX_LEN, TWOPI/BOX_LEN, 1.5, respectively

  NOTE:  the relevant box parameters are taken from INIT_PARAMS.H
*/
int main(int argc, char ** argv){
  fftwf_complex *box;
  unsigned long *in_bin_ct;
  int n_x, n_y, n_z, NUM_BINS, ct;
  float k_x, k_y, k_z, k_mag, k, p;
  float k_floor, k_ceil, k_max, growth_factor, k_first_bin_ceil, k_factor;
  double *p_true, *p_box, *k_ave;
  FILE *OUT, *IN;
  unsigned long long longct;

  /********************* INITIALIZATION **************************/
  // initialize power spectrum crap
  init_ps();
  growth_factor = dicke(INITIAL_REDSHIFT); // normalized to 1 at z=0
  growth_factor = 1;
  // read in arguemnts
  if (argc == 6){
    k_factor = atof(argv[5]);
    k_first_bin_ceil = atof(argv[4]);
    k_max = atof(argv[3]);
  }
  else if (argc == 3){
    k_factor = 1.5;
    k_first_bin_ceil = DELTA_K;
    k_max = DELTA_K*DIM;
  }
  else{ // bad call
    fprintf(stderr, "USAGE: print_power_spec <box_file_name> <output filename> [<k_max> <k_first_bin_ceil> <k_factor>]\nAborting...\n");
    return -1;
  }


  // read in the box
  // allocate array for the k-space and real-space boxes
  box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);
  if (!box){
    fprintf(stderr, "print_power_spec.c: Error allocating memory for box.\nAborting...\n");
    return -1;
  }
  fprintf(stderr, "now opening file for reading\n");
  if (!(IN = fopen(argv[1], "rb"))){
    fprintf(stderr, "Error opening IC file %s file for read\nAborting\n", argv[1]);
    fftwf_free(box); return -1;
  }
  else if (mod_fread(box, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS, 1, IN)!=1){
    fprintf(stderr, "print_power_spec.c: Read error occured!\n");
    return -1;
  }
  fclose(IN);
  /*
  fprintf(stderr, "VALUES AFTER FULL READ\n");
  fprintf(stderr, "%f+%f*I\n", creal(box[C_INDEX(0,0,0)]), cimag(box[C_INDEX(0,0,0)]));
  fprintf(stderr, "%f+%f*I\n", creal(box[C_INDEX(0,10,0)]), cimag(box[C_INDEX(0,10,0)]));
  fprintf(stderr, "%f+%f*I\n", creal(box[C_INDEX(0,0,100)]), cimag(box[C_INDEX(0,0,100)]));
  fprintf(stderr, "%f+%f*I\n", creal(box[C_INDEX(100,0,0)]), cimag(box[C_INDEX(100,0,0)]));
  fprintf(stderr, "%f+%f*I\n", creal(box[C_INDEX(DIM-1,DIM-1,MIDDLE)]), cimag(box[C_INDEX(DIM-1,DIM-1,MIDDLE)]));
  */

  // output file
  OUT = fopen(argv[2], "w");
  if (!OUT){
    fprintf(stderr, "Unable to open file: %s for writing\nAborting...\n", argv[2]);
    return -1;
  }

  // initialize arrays
  // ghetto counting (lookup how to do logs of arbitrary bases in c...)
  NUM_BINS = 0;
  k_floor = 0;
  k_ceil = k_first_bin_ceil;
  while (k_ceil < k_max){
    NUM_BINS++;
    k_floor=k_ceil;
    k_ceil*=k_factor;
  }

  p_true = (double *)malloc(sizeof(double)*NUM_BINS);
  p_box =  (double *)malloc(sizeof(double)*NUM_BINS);
  k_ave =  (double *)malloc(sizeof(double)*NUM_BINS);
  in_bin_ct = (unsigned long *)malloc(sizeof(unsigned long)*NUM_BINS);
  if (!p_true || !p_box || !in_bin_ct || !k_ave){ // a bit sloppy, but whatever..
    fprintf(stderr, "print_power_spec: Error allocating memory.\nAborting...\n");
    return -1;
  }

  for (ct=0; ct<NUM_BINS; ct++){
    p_true[ct] = p_box[ct] = k_ave[ct] = 0;
    in_bin_ct[ct] = 0;
  }
  /************  END INITIALIZATION ******************/


  // now construct the power spectrum file
  for (n_x=0; n_x<DIM; n_x++){
    // convert index to numerical value for this component of the k-mode: k = (2*pi/L) * n
    if (n_x>MIDDLE)
      k_x =(n_x-DIM) * DELTA_K;  // wrap around for FFT convention
    else
      k_x = n_x * DELTA_K;

    for (n_y=0; n_y<DIM; n_y++){
      if (n_y>MIDDLE)
	k_y =(n_y-DIM) * DELTA_K;
      else
	k_y = n_y * DELTA_K;

      // since physical space field is real, only half contains independent modes
      for (n_z=0; n_z<=floor(DIM/2); n_z++){ 
	k_z = n_z * DELTA_K;
	
	k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
	p = power_in_k(k_mag)*growth_factor*growth_factor;


	// now go through the k bins and update
	ct = 0;
	k_floor = 0;
	k_ceil = k_first_bin_ceil;
	while (k_ceil < k_max){
	  // check if we fal in this bin
	  if ((k_mag>=k_floor) && (k_mag < k_ceil)){
	    in_bin_ct[ct]++;
	    p_true[ct] += p;
	    p_box[ct] += pow(cabs(box[C_INDEX(n_x, n_y, n_z)])*growth_factor, 2)/VOLUME;
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

  // now lets print out the k bins
  for (ct=1; ct<NUM_BINS; ct++){
    fprintf(OUT, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_true[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
  }


  // deallocate memory
  fclose(OUT);
  free(p_true);
  free(p_box);
  free(k_ave);
  free(in_bin_ct);
  fftwf_free(box);
}
