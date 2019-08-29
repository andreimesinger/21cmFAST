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
  float *xH, const_factor, *Ts, T_rad, pixel_Ts_factor, curr_alphaX, curr_MminX;
  double ave_Ts, min_Ts, max_Ts, temp, curr_zetaX;
  int ii;
  
  float d1_low, d1_high, d2_low, d2_high, gradient_component, min_gradient_component, subcell_width, x_val1, x_val2, subcell_displacement;
  float RSD_pos_new, RSD_pos_new_boundary_low,RSD_pos_new_boundary_high, fraction_within, fraction_outside, cell_distance;

  float *x_pos_offset, *x_pos, *delta_T_RSD_LOS;
  
  x_pos = calloc(N_RSD_STEPS,sizeof(float));
  x_pos_offset = calloc(N_RSD_STEPS,sizeof(float));
  delta_T_RSD_LOS = calloc(HII_DIM,sizeof(float));
  
  int HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 0;


  /************  BEGIN INITIALIZATION ****************************/
  if (!T_USE_VELOCITIES & SUBCELL_RSD){
    fprintf(stderr, "'SUBCELL_RSD' uses velocities. You MUST turn on 'T_USE_VELOCITIES' to use 'SUBCELL_RSD'.\nAborting!\n");
    return -1; 
  }

  if (SHARP_CUTOFF) HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 1;
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
  //strtok(filename, "f");
  token = strtok(filename, "f");
  nf = atof(strtok(NULL, "_"));

  
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
      if(SUBCELL_RSD) {
          // Converting the prefactors into the optical depth, tau. Factor of 1000 is the conversion of spin temperature from K to mK
          delta_T[HII_R_INDEX(i,j,k)] *= (1. + REDSHIFT)/(1000.*Ts[HII_R_INDEX(i,j,k)]);
      }
      else {
      
          pixel_Ts_factor = (1 - T_rad / Ts[HII_R_INDEX(i,j,k)]);
          delta_T[HII_R_INDEX(i,j,k)] *= pixel_Ts_factor;
      }
      
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

  x_val1 = 0.;
  x_val2 = 1.;
   
  subcell_width = (BOX_LEN/(float)HII_DIM)/(float)N_RSD_STEPS;
   
  float max_cell_distance;
   
  max_cell_distance = 0.;
  


    // check if we need to correct for velocities
  if (!T_USE_VELOCITIES){ //  we can stop here and print
    sprintf(filename, "../Boxes/delta_T_z%06.2f_nf%f_useTs%i_%i_%.0fMpc", REDSHIFT, nf, USE_TS_IN_21CM, HII_DIM, BOX_LEN);
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
  if(SUBCELL_RSD) {
    
        // now add the velocity correction to the delta_T maps
   min_gradient_component = 1.0;
   
   if(USE_TS_IN_21CM) {
       for (i=0; i<HII_DIM; i++){
           for (j=0; j<HII_DIM; j++){
               for (k=0; k<HII_DIM; k++){
    
                   gradient_component = fabs(v[HII_R_FFT_INDEX(i,j,k)]/H + 1.0);
    
                   // Calculate the brightness temperature, using the optical depth
                   if(gradient_component < FRACT_FLOAT_ERR) {
                       // Gradient component goes to zero, optical depth diverges. But, since we take exp(-tau), this goes to zero and (1 - exp(-tau)) goes to unity.
                       // Again, factors of 1000. are conversions from K to mK
                       delta_T[HII_R_INDEX(i,j,k)] = 1000.*(Ts[HII_R_INDEX(i,j,k)] - T_rad)/(1. + REDSHIFT);
                   }   
                   else {
                       delta_T[HII_R_INDEX(i,j,k)] = (1. - exp(- delta_T[HII_R_INDEX(i,j,k)]/gradient_component ))*1000.*(Ts[HII_R_INDEX(i,j,k)] - T_rad)/(1. + REDSHIFT);
                   }   
               }   
           }   
       }   
   }   
   else {
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

           }
         }
       }

   }
  
   // normalised units of cell length. 0 equals beginning of cell, 1 equals end of cell
   // These are the sub-cell central positions (x_pos_offset), and the corresponding normalised value (x_pos) between 0 and 1
   for(ii=0;ii<N_RSD_STEPS;ii++) {
       x_pos_offset[ii] = subcell_width*(float)ii + subcell_width/2.;
       x_pos[ii] = x_pos_offset[ii]/( BOX_LEN/(float)HII_DIM );
   }
   // Note to convert the velocity v, to a displacement in redshift space, convert from s -> r + (1+z)*v/H(z)
   // To convert the velocity within the array v to km/s, it is a*dD/dt*delta. Where the scale factor a comes from the continuity equation
   // The array v as defined in 21cmFAST is (ik/k^2)*dD/dt*delta, as it is defined as a comoving quantity (scale factor is implicit).
   // However, the conversion between real and redshift space also picks up a scale factor, therefore the scale factors drop out and therefore
   // the displacement of the sub-cells is purely determined from the array, v and the Hubble factor: v/H.
     
     for (i=0; i<HII_DIM; i++){
         for (j=0; j<HII_DIM; j++){
             
             // Generate the optical-depth for the specific line-of-sight with R.S.D
             for(k=0;k<HII_DIM;k++) {
                 delta_T_RSD_LOS[k] = 0.0;
             }
             
             for (k=0; k<HII_DIM; k++){
                 
                 if((fabs(delta_T[HII_R_INDEX(i,j,k)]) >= FRACT_FLOAT_ERR) && (xH[HII_R_INDEX(i,j,k)] >= FRACT_FLOAT_ERR)) {

                     if(k==0) {
                         d1_low = v[HII_R_FFT_INDEX(i,j,HII_DIM-1)]/H;
                         d2_low = v[HII_R_FFT_INDEX(i,j,k)]/H;
                     }
                     else {
                         d1_low = v[HII_R_FFT_INDEX(i,j,k-1)]/H;
                         d2_low = v[HII_R_FFT_INDEX(i,j,k)]/H;
                     }
                     // Displacements (converted from velocity) for the original cell centres straddling half of the sub-cells (cell after)
                     if(k==(HII_DIM-1)) {
                         d1_high = v[HII_R_FFT_INDEX(i,j,k)]/H;
                         d2_high = v[HII_R_FFT_INDEX(i,j,0)]/H;
                     }
                     else {
                         d1_high = v[HII_R_FFT_INDEX(i,j,k)]/H;
                         d2_high = v[HII_R_FFT_INDEX(i,j,k+1)]/H;
                     }
                    
                     for(ii=0;ii<N_RSD_STEPS;ii++) {
                         
                         // linearly interpolate the displacements to determine the corresponding displacements of the sub-cells
                         // Checking of 0.5 is for determining if we are left or right of the mid-point of the original cell (for the linear interpolation of the displacement)
                         // to use the appropriate cell
                         
                         if(x_pos[ii] <= 0.5) {
                             subcell_displacement = d1_low + ( (x_pos[ii] + 0.5 ) - x_val1)*( d2_low - d1_low )/( x_val2 - x_val1 );
                         }
                         else {
                             subcell_displacement = d1_high + ( (x_pos[ii] - 0.5 ) - x_val1)*( d2_high - d1_high )/( x_val2 - x_val1 );
                         }
                         
                         // The new centre of the sub-cell post R.S.D displacement. Normalised to units of cell width for determining it's displacement
                         RSD_pos_new = (x_pos_offset[ii] + subcell_displacement)/( BOX_LEN/(float)HII_DIM );
                         // The sub-cell boundaries of the sub-cell, for determining the fractional contribution of the sub-cell to neighbouring cells when
                         // the sub-cell straddles two cell positions
                         RSD_pos_new_boundary_low = RSD_pos_new - (subcell_width/2.)/( BOX_LEN/(float)HII_DIM );
                         RSD_pos_new_boundary_high = RSD_pos_new + (subcell_width/2.)/( BOX_LEN/(float)HII_DIM );
                         
                         if(RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_high < 1.0) {
                             // sub-cell has remained in the original cell (just add it back to the original cell)
                             
                             delta_T_RSD_LOS[k] += delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                         }
                         else if(RSD_pos_new_boundary_low < 0.0 && RSD_pos_new_boundary_high < 0.0) {
                             // sub-cell has moved completely into a new cell (toward the observer)
                             
                             // determine how far the sub-cell has moved in units of original cell boundary
                             cell_distance = ceil(fabs(RSD_pos_new_boundary_low))-1.;
                             
                             // Determine the location of the sub-cell relative to the original cell binning
                             if(fabs(RSD_pos_new_boundary_high) > cell_distance) {
                                 // sub-cell is entirely contained within the new cell (just add it to the new cell)
                                // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                                 if(k<((int)cell_distance+1)) {
                                     delta_T_RSD_LOS[k-((int)cell_distance+1) + HII_DIM] += delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 else {
                                     delta_T_RSD_LOS[k-((int)cell_distance+1)] += delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                             }
                             else {
                                 // sub-cell is partially contained within the cell
                                 
                                 // Determine the fraction of the sub-cell which is in either of the two original cells
                                 fraction_outside = (fabs(RSD_pos_new_boundary_low) - cell_distance)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                 fraction_within = 1. - fraction_outside;
                                 
                                 // Check if the first part of the sub-cell is at the box edge
                                 if(k<(((int)cell_distance))) {
                                     delta_T_RSD_LOS[k-((int)cell_distance) + HII_DIM] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 else {
                                     delta_T_RSD_LOS[k-((int)cell_distance)] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 // Check if the second part of the sub-cell is at the box edge
                                 if(k<(((int)cell_distance + 1))) {
                                     delta_T_RSD_LOS[k-((int)cell_distance+1) + HII_DIM] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 else {
                                     delta_T_RSD_LOS[k-((int)cell_distance+1)] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                             }
                         }
                        else if(RSD_pos_new_boundary_low < 0.0 && (RSD_pos_new_boundary_high > 0.0 && RSD_pos_new_boundary_high < 1.0)) {
                             // sub-cell has moved partially into a new cell (toward the observer)
                             
                             // Determine the fraction of the sub-cell which is in either of the two original cells
                             fraction_within = RSD_pos_new_boundary_high/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                             fraction_outside = 1. - fraction_within;
                             
                             // Check the periodic boundaries conditions and move the fraction of each sub-cell to the appropriate new cell
                             if(k==0) {
                                 delta_T_RSD_LOS[HII_DIM-1] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 delta_T_RSD_LOS[k] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                             }
                             else {
                                 delta_T_RSD_LOS[k-1] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 delta_T_RSD_LOS[k] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                             }
                         }
                         else if((RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_low < 1.0) && (RSD_pos_new_boundary_high >= 1.0)) {
                             // sub-cell has moved partially into a new cell (away from the observer)
                             
                             // Determine the fraction of the sub-cell which is in either of the two original cells
                             fraction_outside = (RSD_pos_new_boundary_high - 1.)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                             fraction_within = 1. - fraction_outside;
                             
                             // Check the periodic boundaries conditions and move the fraction of each sub-cell to the appropriate new cell
                             if(k==(HII_DIM-1)) {
                                 delta_T_RSD_LOS[k] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 delta_T_RSD_LOS[0] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                             }
                             else {
                                 delta_T_RSD_LOS[k] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 delta_T_RSD_LOS[k+1] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                             }
                         }
                        else {
                             // sub-cell has moved completely into a new cell (away from the observer)
                             
                             // determine how far the sub-cell has moved in units of original cell boundary
                             cell_distance = floor(fabs(RSD_pos_new_boundary_high));
                             
                             if(RSD_pos_new_boundary_low >= cell_distance) {
                                 // sub-cell is entirely contained within the new cell (just add it to the new cell)
                                 
                                 // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                                 if(k>(HII_DIM - 1 - (int)cell_distance)) {
                                     delta_T_RSD_LOS[k+(int)cell_distance - HII_DIM] += delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 else {
                                     delta_T_RSD_LOS[k+(int)cell_distance] += delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                             }
                             else {
                                 // sub-cell is partially contained within the cell
                                 
                                 // Determine the fraction of the sub-cell which is in either of the two original cells
                                 fraction_outside = (RSD_pos_new_boundary_high - cell_distance)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                 fraction_within = 1. - fraction_outside;
                                 
                                 // Check if the first part of the sub-cell is at the box edge
                                 if(k>(HII_DIM - 1 - ((int)cell_distance-1))) {
                                     delta_T_RSD_LOS[k+(int)cell_distance-1 - HII_DIM] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 else {
                                     delta_T_RSD_LOS[k+(int)cell_distance-1] += fraction_within*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 // Check if the second part of the sub-cell is at the box edge
                                 if(k>(HII_DIM - 1 - ((int)cell_distance))) {
                                     delta_T_RSD_LOS[k+(int)cell_distance - HII_DIM] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                                 else {
                                     delta_T_RSD_LOS[k+(int)cell_distance] += fraction_outside*delta_T[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                 }
                            }
                         }
                     }
                 }
             }
             
             for(k=0;k<HII_DIM;k++) {
                 delta_T[HII_R_INDEX(i,j,k)] = delta_T_RSD_LOS[k];
                 
                 ave += delta_T_RSD_LOS[k];
             }
         }
     }
    }
    else {
          
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
    }
      
     ave /= (HII_TOT_NUM_PIXELS+0.0);
      
      if(!SUBCELL_RSD) {
          fprintf(LOG, "With velocities:\nMax is %e\t dvdx is %e, ave is %e\n", max, maxdvdx, ave);
          fprintf(stderr, "With velocities:\nMax is %e\t dvdx is %e, ave is %e\n", max, maxdvdx, ave);
          fprintf(LOG, "%llu out of %llu voxels (fraction=%e) exceeded max allowed velocity gradient\n", nonlin_ct, HII_TOT_NUM_PIXELS, nonlin_ct/(double)HII_TOT_NUM_PIXELS);
          fprintf(stderr, "%llu out of %llu voxels (fraction=%e) exceeded max allowed velocity gradient\n", nonlin_ct, HII_TOT_NUM_PIXELS, nonlin_ct/(double)HII_TOT_NUM_PIXELS);
      }
      

  
  // now write out the delta_T box with velocity correction
  sprintf(filename, "../Boxes/delta_T_v%i_z%06.2f_nf%f_useTs%i_%i_%.0fMpc", VELOCITY_COMPONENT, REDSHIFT, nf, USE_TS_IN_21CM, HII_DIM, BOX_LEN);
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
    sprintf(filename, "%s/ps_z%06.2f_nf%f_useTs%i_aveTb%06.2f_%i_%.0fMpc_v%i", psoutputdir, REDSHIFT, nf, USE_TS_IN_21CM, ave, HII_DIM, BOX_LEN, VELOCITY_COMPONENT);
  }
  else{
    sprintf(filename, "%s/ps_z%06.2f_nf%f_useTs%i_aveTb%06.2f_%i_%.0fMpc", psoutputdir, REDSHIFT, nf, USE_TS_IN_21CM, ave, HII_DIM, BOX_LEN);
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
