#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"


/*
  USAGE: delta_T [-p <NUM THREADS>] <redshift> <xH filename> [<Ts filename>]

  generates the 21-cm temperature offset from the CMB field and power spectrum 
  the spin temperature filename is optional

  NOTE: the optional argument of thread number including the -p flag, MUST
  be the first two arguments.  If these are omitted, num_threads defaults
  to NUMCORES in INIT_PARAMS.H
*/

int main(int argc, char ** argv){
  fftwf_complex *deldel_T;
  fftwf_plan plan;
  char filename[1000], psoutputdir[1000], *token;
  float *deltax, REDSHIFT, growth_factor, dDdt, pixel_x_HI, pixel_deltax, *delta_T, *v, H, dummy;
  FILE *F, *LOG;
  int i,j,k, n_x, n_y, n_z, NUM_BINS, curr_Pop, arg_offset,num_th;
  double dvdx, ave, *p_box, *k_ave, max_v_deriv;
  unsigned long long ct, *in_bin_ct, nonlin_ct, temp_ct;
  float nf, max, maxi, maxj, maxk, maxdvdx, min, mini, minj, mink, mindvdx;
  float k_x, k_y, k_z, k_mag, k_sq, k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor;
  float *xH, const_factor, *Ts, T_rad, pixel_Ts_factor, curr_alphaX, curr_TvirX;
  double ave_Ts, min_Ts, max_Ts, temp, curr_zetaX;


  /************  BEGIN INITIALIZATION ****************************/
  max = -1e3;
  min = 1e3;
  ave = 0;
  nonlin_ct=0;

  // check arguments
  if (argc < 3){
    fprintf(stderr, "USAGE: delta_T [-p <NUM THREADS>] <redshift> <xH filename> [<Ts filename>]\nAborting!\n");
    return -1;
  }
  if ( (argv[1][0]=='-') && ((argv[1][1]=='p') || (argv[1][1]=='P')) ){
    // user specified num proc
    num_th = atoi(argv[2]);
    fprintf(stderr, "delta_T: threading with user-specified %i threads\n", num_th);
    arg_offset = 2;
  }
  else{
    num_th = NUMCORES;
    fprintf(stderr, "delta_T: threading with default %i threads\n", num_th);
    arg_offset = 0;
  }
  if (fftwf_init_threads()==0){
    fprintf(stderr, "delta_T: ERROR: problem initializing fftwf threads\nAborting\n.");
    return -1;
  }
  fftwf_plan_with_nthreads(num_th);    

  // open LOG file
  REDSHIFT = atof(argv[1+arg_offset]);
  system("mkdir ../Log_files");
  sprintf(filename, "../Log_files/delta_T_log_file_%d", getpid());
  LOG = fopen(filename, "w");
  if (!LOG){ fprintf(stderr, "delta_T.c: Error opening log file %s\n", filename);}
  T_rad = T_cmb*(1+REDSHIFT);
  H = hubble(REDSHIFT);
  const_factor = 27 * (OMb*hlittle*hlittle/0.023) * 
    sqrt( (0.15/OMm/hlittle/hlittle) * (1+REDSHIFT)/10.0 );
  system("mkdir ../Output_files/");
  system("mkdir ../Output_files/Deldel_T_power_spec");
  system("mkdir ../Log_files");

  // get the neutral fraction and HII filter from the filename
  strcpy(filename, argv[2+arg_offset]);
  strtok(filename, "f");
  printf("strtok(filename, 'f') %s\n",strtok(filename, "f"));
  nf = atof(strtok(NULL, "_"));
  if (argc == (4+arg_offset)){
    // get zetaX and HII filter from the filename
    strcpy(filename, argv[3+arg_offset]);
    token = strtok(filename, "X");
    curr_zetaX = atof(strtok(NULL, "_"));
    token = strtok(NULL, "X");
    curr_alphaX = atof(strtok(NULL, "_"));
    token = strtok(NULL, "X");
    curr_TvirX = atof(strtok(NULL, "_"));
    strcpy(filename, argv[3+arg_offset]);
    token = strtok(filename, "P");
    token = strtok(NULL, "p");
    curr_Pop = atoi(strtok(NULL, "_"));
    if (curr_Pop>4){curr_Pop=2;} // due to a previous typo
  }
  else{
    curr_zetaX = -1;
    curr_alphaX = -1;
    curr_TvirX = -1;
    curr_Pop = -1;
  }
  // initialize power spectrum 
  init_ps();
  growth_factor = dicke(REDSHIFT); // normalized to 1 at z=0


  // allocate memory for xH box and read it in
  xH = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!xH){
    fprintf(stderr, "delta_T: Error allocating memory for xH box\nAborting...\n");
    fprintf(LOG, "delta_T: Error allocating memory for xH box\nAborting...\n");
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  if (!(F = fopen(argv[2+arg_offset], "rb"))){
    fprintf(stderr, "delta_T: unable to open xH box at %s\nAborting...\n", argv[2+arg_offset]);
    fprintf(LOG, "delta_T: unable to open xH box at %s\nAborting...\n", argv[2+arg_offset]);
    free(xH);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  fprintf(stderr, "Reading in xH box\n");
  fprintf(LOG, "Reading in xH box\n");
  if (mod_fread(xH, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    fprintf(stderr, "delta_T: Read error occured while reading neutral_fraction box.\n");
    fprintf(LOG, "delta_T: Read error occured while reading neutral_fraction box.\n");
    fclose(F); free(xH);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  fclose(F);
 
  // allocate memory for deltax box and read it in
  deltax = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  if (!deltax){
    fprintf(stderr, "delta_T: Error allocating memory for deltax box\nAborting...\n");
    fprintf(LOG, "delta_T: Error allocating memory for deltax box\nAborting...\n");
    free(xH);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  if (!(F = fopen(filename, "rb"))){
    fprintf(stderr, "delta_T: Error openning deltax box for reading at %s\n", filename);
    fprintf(LOG, "delta_T: Error openning deltax box for reading at %s\n", filename);
    free(xH); free(deltax);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  fprintf(stderr, "Reading in deltax box\n");
  fprintf(LOG, "Reading in deltax box\n");
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
	if (fread((float *)deltax + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	  fprintf(stderr, "delta_T: Read error occured while reading deltax box.\n");
	  fprintf(LOG, "delta_T: Read error occured while reading deltax box.\n");
	  fclose(F); free(xH); free(deltax);
	  fclose(LOG); fftwf_cleanup_threads(); return -1;
	}
      }
    }
  }
  fclose(F);


  // allocate memory for our delta_T box
  delta_T = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!delta_T){
    fprintf(stderr, "delta_T: Error allocating memory for delta_T box\nAborting...\n");
    fprintf(LOG, "delta_T: Error allocating memory for delta_T box\nAborting...\n");
    free(xH); free(deltax);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }

  // allocate memory for the velocity box and read it in
  v = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  if (!v){
    fprintf(stderr, "delta_T: Error allocating memory for velocity box\nAborting...\n");
    fprintf(LOG, "delta_T: Error allocating memory for velocity box\nAborting...\n");
    free(xH); free(deltax); free(delta_T);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  switch(VELOCITY_COMPONENT){
  case 1:  sprintf(filename, "../Boxes/updated_vx_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
    break;
  case 3:  sprintf(filename, "../Boxes/updated_vz_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
    break;
  default: sprintf(filename, "../Boxes/updated_vy_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  }
  if (T_USE_VELOCITIES){
    if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "delta_T: Error opening velocity file at %s\n", filename);
      fprintf(LOG, "delta_T: Error opening velocity file at %s\n", filename);
      free(xH); free(deltax); free(delta_T); free(v);
      fclose(LOG); fftwf_cleanup_threads(); return -1;
    }
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if (fread((float *)v + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	    fprintf(stderr, "delta_T: Read error occured while reading velocity box.\n");
	    fprintf(LOG, "delta_T: Read error occured while reading velocity box.\n");
	    fclose(F); free(xH); free(deltax); free(delta_T); free(v);
	    fclose(LOG); fftwf_cleanup_threads(); return -1;
	  }
	}
      }
    }
    fclose(F);
  }

  if (USE_TS_IN_21CM){
    // and allocate memory for the spin temperature box and read it in
    if (!(Ts = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
      fprintf(stderr, "delta_T.c: Error in memory allocation for Ts box\nAborting...\n");
      fprintf(LOG, "delta_T.c: Error in memory allocation for Ts box\nAborting...\n");
      free(xH); free(deltax); free(delta_T); free(v);
      fclose(LOG); fftwf_cleanup_threads(); return -1;
    }
    if (!(F = fopen(argv[3+arg_offset], "rb") )){
      fprintf(stderr, "delta_T.c: Error openning Ts file %s to read from\nAborting...\n", argv[3+arg_offset]);
      fprintf(LOG, "delta_T.c: Error openning Ts file %s to read from\nAborting...\n", argv[3+arg_offset]);
      free(xH); free(deltax); free(delta_T); free(v); free(Ts);
      fclose(LOG); fftwf_cleanup_threads(); return -1;
    }
    if (mod_fread(Ts, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
      fprintf(LOG, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
      free(xH); free(deltax); free(delta_T); free(v); free(Ts); fclose(F);
      fclose(LOG); fftwf_cleanup_threads(); return -1;
    }
    fclose(F);
  }

  /************  END INITIALIZATION ****************************/

  // ok, lets fill the delta_T box; which will be the same size as the bubble box
  ave_Ts = max_Ts = 0;
  min_Ts = 1e5;
  temp=0;
  temp_ct=0;
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){

	pixel_deltax = deltax[HII_R_FFT_INDEX(i,j,k)];
	pixel_x_HI = xH[HII_R_INDEX(i,j,k)];

	if (pixel_x_HI > TINY){
	  temp = pixel_deltax;
	  temp_ct++;
	}

	delta_T[HII_R_INDEX(i,j,k)] = const_factor*pixel_x_HI*(1+pixel_deltax);

	if (USE_TS_IN_21CM){
	  pixel_Ts_factor = (1 - T_rad / Ts[HII_R_INDEX(i,j,k)]);
	  delta_T[HII_R_INDEX(i,j,k)] *= pixel_Ts_factor;
	  ave_Ts += Ts[HII_R_INDEX(i,j,k)];
	  if (min_Ts > Ts[HII_R_INDEX(i,j,k)]) { min_Ts = Ts[HII_R_INDEX(i,j,k)];}
	  if (max_Ts < Ts[HII_R_INDEX(i,j,k)]) { max_Ts = Ts[HII_R_INDEX(i,j,k)];}
	}

	if (max < delta_T[HII_R_INDEX(i,j,k)]){ max = delta_T[HII_R_INDEX(i,j,k)];}
	if (min > delta_T[HII_R_INDEX(i,j,k)]){ min = delta_T[HII_R_INDEX(i,j,k)];}
	ave += delta_T[HII_R_INDEX(i,j,k)];
      }
    }
  }
  ave /= (float)HII_TOT_NUM_PIXELS;
  fprintf(stderr, "Without velocities, max is %e, min is %e, ave is %e\n", max, min, ave);
  fprintf(LOG, "Without velocities, max is %e, min is %e, ave is %e\n", max, min, ave);

  if (USE_TS_IN_21CM){
    ave_Ts /= (double) HII_TOT_NUM_PIXELS;
    fprintf(stderr, "Ts, min = %e, max = %e, ave = %e\n", min_Ts, max_Ts, ave_Ts);
    fprintf(stderr, "corresponding to (1-trad/Ts of), min = %e, max = %e, ave = %e\n", 1-T_rad/min_Ts, 1-T_rad/max_Ts, 1-T_rad/ave_Ts);
  }


  // now write out the delta_T box
  if (!T_USE_VELOCITIES){
    switch(FIND_BUBBLE_ALGORITHM){
    case 2: 
      if(HII_EFF_FACTOR_SFR == 0){
        if (USE_HALO_FIELD)
  	    sprintf(filename, "../Boxes/delta_T_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
        else
  	    sprintf(filename, "../Boxes/delta_T_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX,ave, curr_Pop, HII_DIM, BOX_LEN);
  	  }
	  else{
        if (USE_HALO_FIELD)
	    sprintf(filename, "../Boxes/delta_T_SFR_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
        else
	    sprintf(filename, "../Boxes/delta_T_no_halos_SFR_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX,ave, curr_Pop, HII_DIM, BOX_LEN);
	  } 
      break;

    default:
	  if(HII_EFF_FACTOR_SFR == 0){
        if (USE_HALO_FIELD)
	    sprintf(filename, "../Boxes/sphere_delta_T_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX,ave, curr_Pop, HII_DIM, BOX_LEN);
        else
	    sprintf(filename, "../Boxes/sphere_delta_T_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
	  }
	  else{
        if (USE_HALO_FIELD)
	    sprintf(filename, "../Boxes/sphere_delta_T_SFR_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX,ave, curr_Pop, HII_DIM, BOX_LEN);
        else
	    sprintf(filename, "../Boxes/sphere_delta_T_no_halos_SFR_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
	  }
    }
    F = fopen(filename, "wb");
    fprintf(stderr, "\nWritting output delta_T box: %s\n", filename);
    if (mod_fwrite(delta_T, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "delta_T: Write error occured while writting delta_T box.\n");
    }
    fclose(F);
  }
  else{
    max = -1;
    min = 1e3;
    ave = 0;

  // let's take the derivative in k-space
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)v, (fftwf_complex *)v, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
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

	// take partial deriavative along the line of sight
	switch(VELOCITY_COMPONENT){
	case 1:
	  *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_x*I/(float)HII_TOT_NUM_PIXELS;
	  break;
	case 3:
	  *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_z*I/(float)HII_TOT_NUM_PIXELS;
	  break;
	default:
	  *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_y*I/(float)HII_TOT_NUM_PIXELS;
	}
      }
    }
  }
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)v, (float *)v, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();


  // now add the velocity correction to the delta_T maps
  max_v_deriv = fabs(MAX_DVDR*H);
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){

	dvdx = v[HII_R_FFT_INDEX(i,j,k)];

	// set maximum allowed gradient for this linear approximation
        if (fabs(dvdx) > max_v_deriv){
          if (dvdx < 0) dvdx = -max_v_deriv;
          else dvdx = max_v_deriv;
          nonlin_ct++;
        }

	delta_T[HII_R_INDEX(i,j,k)] /= (dvdx/H + 1.0);

	if (max < delta_T[HII_R_INDEX(i,j,k)]){
	  maxi = i;
	  maxj = j;
	  maxk = k;
	  maxdvdx = dvdx;
	  max = delta_T[HII_R_INDEX(i,j,k)];
	}
	if (min > delta_T[HII_R_INDEX(i,j,k)]){
	  mini = i;
	  minj = j;
	  mink = k;
	  mindvdx = dvdx;
	  min = delta_T[HII_R_INDEX(i,j,k)];
	}

	ave += delta_T[HII_R_INDEX(i,j,k)];
      }
    }
  }
  ave /= (HII_TOT_NUM_PIXELS+0.0);
  fprintf(LOG, "With velocities:\nMax is %e\t dvdx is %e, ave is %e\n", max, maxdvdx, ave);
  fprintf(stderr, "With velocities:\nMax is %e\t dvdx is %e, ave is %e\n", max, maxdvdx, ave);
  fprintf(LOG, "%llu out of %llu voxels (fraction=%e) exceeded max allowed velocity gradient\n", nonlin_ct, HII_TOT_NUM_PIXELS, nonlin_ct/(double)HII_TOT_NUM_PIXELS);
  fprintf(stderr, "%llu out of %llu voxels (fraction=%e) exceeded max allowed velocity gradient\n", nonlin_ct, HII_TOT_NUM_PIXELS, nonlin_ct/(double)HII_TOT_NUM_PIXELS);

  
  // now write out the delta_T box with velocity correction
  if(FIND_BUBBLE_ALGORITHM==2){
	if(HII_EFF_FACTOR_SFR == 0){
      if (USE_HALO_FIELD)
        sprintf(filename, "../Boxes/delta_T_v%i_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
      else
        sprintf(filename, "../Boxes/delta_T_v%i_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
	}
	else{
      if (USE_HALO_FIELD)
        sprintf(filename, "../Boxes/delta_T_SFR_v%i_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
      else
        sprintf(filename, "../Boxes/delta_T_SFR_v%i_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
	}
  }
  else{
	if(HII_EFF_FACTOR_SFR == 0){
      if (USE_HALO_FIELD)
        sprintf(filename, "../Boxes/sphere_delta_T_v%i_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX,  ave, curr_Pop, HII_DIM, BOX_LEN);
      else
        sprintf(filename, "../Boxes/sphere_delta_T_v%i_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX,  curr_alphaX, curr_TvirX, ave,  curr_Pop, HII_DIM, BOX_LEN);
	}
	else{
      if (USE_HALO_FIELD)
        sprintf(filename, "../Boxes/sphere_delta_T_SFR_v%i_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX,  ave, curr_Pop, HII_DIM, BOX_LEN);
      else
        sprintf(filename, "../Boxes/sphere_delta_T_SFR_v%i_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX,  curr_alphaX, curr_TvirX, ave,  curr_Pop, HII_DIM, BOX_LEN);
	}
  }
  F = fopen(filename, "wb");
  fprintf(stderr, "Writting output delta_T box: %s\n", filename);
  if (mod_fwrite(delta_T, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    fprintf(stderr, "delta_T: Write error occured while writting delta_T box.\n");
  }
  fclose(F);
}

// deallocate what we aren't using anymore 
 free(xH); free(deltax); free(v); if (USE_TS_IN_21CM){ free(Ts);}




  /******  PRINT OUT THE POWERSPECTRUM  *********/

 k_factor = 1.5;
  k_first_bin_ceil = DELTA_K;
  k_max = DELTA_K*HII_DIM;
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

  p_box =  (double *)malloc(sizeof(double)*NUM_BINS);
  k_ave =  (double *)malloc(sizeof(double)*NUM_BINS);
  in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*NUM_BINS);
  if (!p_box || !in_bin_ct || !k_ave){ // a bit sloppy, but whatever..
    fprintf(stderr, "delta_T.c: Error allocating memory.\nAborting...\n");
    fprintf(LOG, "delta_T.c: Error allocating memory.\nAborting...\n");
    free(delta_T); fclose(LOG);
    fftwf_cleanup_threads(); return -1;
  }
  for (ct=0; ct<NUM_BINS; ct++){
    p_box[ct] = k_ave[ct] = 0;
    in_bin_ct[ct] = 0;
  }

  deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!deldel_T){
    fprintf(stderr, "Unable to allocate memory for the deldel_T box!\n");
    fprintf(LOG, "Unable to allocate memory for the deldel_T box!\n");
    free(delta_T); fclose(LOG); free(p_box); free(k_ave); free(in_bin_ct);
    fftwf_cleanup_threads(); return -1;
  }

  // fill-up the real-space of the deldel box
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
	*((float *)deldel_T + HII_R_FFT_INDEX(i,j,k)) = (delta_T[HII_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);
	if (DIMENSIONAL_T_POWER_SPEC){
	  *((float *)deldel_T + HII_R_FFT_INDEX(i,j,k)) *= ave;
	}
	// Note: we include the V/N factor for the scaling after the fft
      }
    }
  }

  // transform to k-space
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T, (fftwf_complex *)deldel_T, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();

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
	    p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
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


  if (DIMENSIONAL_T_POWER_SPEC)
    sprintf(psoutputdir, "../Output_files/Deldel_T_power_spec");
  else
    sprintf(psoutputdir, "../Output_files/Deldel_T_power_spec/Dimensionless");
  sprintf(filename, "mkdir %s", psoutputdir);
  system(filename);

  // now lets print out the k bins
  if (T_USE_VELOCITIES){
	if(HII_EFF_FACTOR_SFR == 0){
      if (USE_HALO_FIELD)
        sprintf(filename, "%s/ps_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc_v%i", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN, VELOCITY_COMPONENT);
      else
        sprintf(filename, "%s/ps_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc_v%i", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN, VELOCITY_COMPONENT);
	}
	else{
      if (USE_HALO_FIELD)
        sprintf(filename, "%s/ps_SFR_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc_v%i", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN, VELOCITY_COMPONENT);
      else
        sprintf(filename, "%s/ps_no_halos_SFR_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc_v%i", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN, VELOCITY_COMPONENT);
	}
  }
  else{
	if(HII_EFF_FACTOR_SFR == 0){
      if (USE_HALO_FIELD)
        sprintf(filename, "%s/ps_nov_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
      else
        sprintf(filename, "%s/ps_nov_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
	}
	else{
      if (USE_HALO_FIELD)
        sprintf(filename, "%s/ps_SFR_nov_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
      else
        sprintf(filename, "%s/ps_SFR_nov_no_halos_z%06.2f_nf%f_useTs%i_zetaX%.1e_alphaX%.1f_TvirminX%.1e_aveTb%06.2f_Pop%i_%i_%.0fMpc", psoutputdir,REDSHIFT, nf, USE_TS_IN_21CM, curr_zetaX, curr_alphaX, curr_TvirX, ave, curr_Pop, HII_DIM, BOX_LEN);
	}
  }

  F = fopen(filename, "w");
  if (!F){
    fprintf(stderr, "delta_T.c: Couldn't open file %s for writting!\n", filename);
    fprintf(LOG, "delta_T.c: Couldn't open file %s for writting!\n", filename);
    free(delta_T); fclose(LOG); free(p_box); free(k_ave); free(in_bin_ct); fftwf_free(deldel_T);
  }
  for (ct=1; ct<NUM_BINS; ct++){
    if (in_bin_ct[ct]>0)
      fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
  }
  fclose(F); free(p_box); free(k_ave); free(in_bin_ct); fftwf_free(deldel_T);

  /****** END POWER SPECTRUM STUFF   ************/


  // deallocate
  free(delta_T); fclose(LOG);
  fftwf_cleanup_threads(); return 0;
}
