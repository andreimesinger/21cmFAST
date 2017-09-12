#include "heating_helper_progs.c"


/*
  Program Ts calculates the spin temperature field, according to the perscription outlined in
  Mesinger, Furlanetto, Cen (2010).  The fluctuating component is sourced by the mean EPS collapsed fraction in spherical annyli surrounding each cell.

  Usage: Ts <REDSHIFT> [reload zp redshift]

  The last optional argument is the z' redshift output from which to reload intermediate
  evolution files in ../Boxes/Ts_evolution/

  Memory usage (in floats)~ (<NUMBER OF FILTER STEPS> + 3) x HII_DIM^3

  Author: Andrei Mesinger
  Date: 9.2.2009
*/


/* Maximum allowed value for the kinetic temperature.
   Useful to set to avoid some spurious behaviour 
   when the code is run with redshift poor resolution 
   and very high X-ray heating efficiency */
#define MAX_TK (float) 5e4 


int main(int argc, char ** argv){
  fftwf_complex *box, *unfiltered_box;
  fftwf_plan plan;
  unsigned long long ct, sample_ct;
  int R_ct,i,j,k, COMPUTE_Ts, x_e_ct;
  float REDSHIFT, growth_factor_z, R, R_factor, zp, mu_for_Ts, filling_factor_of_HI_zp;
  float *Tk_box, *x_e_box, *Ts, J_star_Lya, dzp, prev_zp, zpp, prev_zpp, prev_R;
  FILE *F, *GLOBAL_EVOL, *OUT;
  char filename[300];
  float dz, zeta_ion_eff, Tk_BC, xe_BC, nu, zprev, zcurr, curr_delNL0[NUM_FILTER_STEPS_FOR_Ts];
  double *evolve_ans, ans[2], dansdz[5], Tk_ave, J_alpha_ave, xalpha_ave, J_alpha_tot, Xheat_ave,
    Xion_ave;
double freq_int_heat_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts];
  int goodSteps,badSteps;
  int m_xHII_low, m_xHII_high, n_ct, zp_ct;
double freq_int_heat[NUM_FILTER_STEPS_FOR_Ts], freq_int_ion[NUM_FILTER_STEPS_FOR_Ts], freq_int_lya[NUM_FILTER_STEPS_FOR_Ts];
 double nuprime, fcoll_R, Ts_ave;
 float R_values[NUM_FILTER_STEPS_FOR_Ts], *delNL0[NUM_FILTER_STEPS_FOR_Ts], xHII_call, curr_xalpha;
 float z, Jalpha, TK, TS, xe, deltax;
 time_t start_time, curr_time;
 double J_alpha_threads[NUMCORES], xalpha_threads[NUMCORES], Xheat_threads[NUMCORES],
   Xion_threads[NUMCORES], lower_int_limit;
 float M_MIN_WDM =  M_J_WDM();
 int RESTART = 0;

 /*
 printf("Ratio of the cross-sections: %e\n Wieghted by the number density: %e\n", 
	HeI_ion_crosssec(HeI_NUIONIZATION)/HI_ion_crosssec(HeI_NUIONIZATION),
	He_No*HeI_ion_crosssec(HeI_NUIONIZATION)/(No*HI_ion_crosssec(HeI_NUIONIZATION)) 
);
	
 return 0;
 */
 /**********  BEGIN INITIALIZATION   **************************************/

 if (argc == 3){
   RESTART = 1;
   zp = atof(argv[2]);
 }
 else if (argc != 2){
   fprintf(stderr, "Usage: Ts <REDSHIFT>  [reload zp redshift]\nAborting...\n");
   return -1;
 }
 REDSHIFT = atof(argv[1]);
 system("mkdir ../Log_files");
 system("mkdir ../Output_files");
 system("mkdir ../Boxes/Ts_evolution/");
 system("mkdir ../Output_files/Ts_outs/");
 system("cp ../Parameter_files/* ../Output_files/Ts_outs/");
 system("cp ../Parameter_files/* ../Boxes/Ts_evolution/");
 init_ps();
 omp_set_num_threads(NUMCORES);
 growth_factor_z = dicke(REDSHIFT);
 if (X_RAY_Tvir_MIN < 9.99999e3) // neutral IGM
   mu_for_Ts = 1.22;
 else // ionized IGM
   mu_for_Ts = 0.6;
 
 //set the minimum ionizing source mass
 M_MIN_at_z = get_M_min_ion(REDSHIFT);
  
 // open log file
 if (!(LOG = fopen("../Log_files/Ts_log", "w") ) ){
   fprintf(stderr, "Unable to open log file for writting\nAborting...\n");
   return -1;
 }

 // Initialize some interpolation tables
 if (init_heat() < 0){
   fclose(LOG); fclose(GLOBAL_EVOL);
   return -1;
 }

 // check if we are in the really high z regime before the first stars; if so, simple
 if (REDSHIFT > Z_HEAT_MAX){
//(FgtrM(REDSHIFT, FMAX(TtoM(REDSHIFT, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM)) < 1e-15 ){
   xe = xion_RECFAST(REDSHIFT,0);
   TK = T_RECFAST(REDSHIFT,0);
   
   // open input
   sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", 
	   REDSHIFT, HII_DIM, BOX_LEN);
   F = fopen(filename, "rb");
   if ( !(F = fopen(filename, "rb") ) ){
     fprintf(stderr, "Error opening file %s for reading.\nAborting...\n", filename);
     fprintf(LOG, "Error opening file %s for reading.\nAborting...\n", filename);
     destruct_heat(); return -1;
   }
   fprintf(stderr, "Opened density file %s for reading\n", filename);
   fprintf(LOG, "Opened density file %s for reading\n", filename);

   // open output
   sprintf(filename, "../Boxes/Ts_z%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN); 
   if (!(OUT=fopen(filename, "wb"))){
     fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
     fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
     destruct_heat(); return -1;
   }
   fprintf(stderr, "Opened TS file %s for writting\n", filename);
   fprintf(LOG, "Opened TS file %s for writting\n", filename);

   // read file
   for (i=0; i<HII_DIM; i++){
     for (j=0; j<HII_DIM; j++){
       for (k=0; k<HII_DIM; k++){
	 if (fread(&deltax, sizeof(float), 1, F)!=1){
	  fprintf(stderr, "Error reading-in binary density file\nAborting...\n");
	  fprintf(LOG, "Error reading-in binary density file\nAborting...\n");
	  destruct_heat(); return -1;
	 }

	 // compute the spin temperature
	TS = get_Ts(REDSHIFT, deltax, TK, xe, 0, &curr_xalpha);

	// and print it out
	if (fwrite(&TS, sizeof(float), 1, OUT)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	  destruct_heat(); return -1;
	 }

       }
     }
   }

   destruct_heat(); fclose(F); fclose(OUT);
   return 0;
 }


 // open global evolution output file
 sprintf(filename, "../Output_files/Ts_outs/global_evolution_zetaIon%.2f_Nsteps%i_zprimestepfactor%.3f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_Pop%i_%i_%.0fMpc", HII_EFF_FACTOR, NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, Pop, HII_DIM, BOX_LEN);
 if (argc > 2) // restarting
   GLOBAL_EVOL = fopen(filename, "a");
 else
   GLOBAL_EVOL = fopen(filename, "w");
 if (!GLOBAL_EVOL){
   fprintf(stderr, "Unable to open global evolution file at %s\nAborting...\n",
	   filename);
   fprintf(LOG, "Unable to open global evolution file at %s\nAborting...\n",
	   filename);
   fclose(LOG);
   return -1;
 }
 


 // set boundary conditions for the evolution equations->  values of Tk and x_e at Z_HEAT_MAX
 if (XION_at_Z_HEAT_MAX > 0) // user has opted to use his/her own value
   xe_BC = XION_at_Z_HEAT_MAX;
 else// will use the results obtained from recfast
   xe_BC = xion_RECFAST(Z_HEAT_MAX,0);
 if (TK_at_Z_HEAT_MAX > 0)
   Tk_BC = TK_at_Z_HEAT_MAX;
 else
   Tk_BC = T_RECFAST(Z_HEAT_MAX,0);


  /******  Now allocate large arrays  ******/

  // allocate memory for the nonlinear density field and open file
  sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", 
	  REDSHIFT, HII_DIM, BOX_LEN);
  F = fopen(filename, "rb");
  if ( !(F = fopen(filename, "rb") ) ){
    fprintf(stderr, "Error opening file %s for reading.\nAborting...\n", filename);
    fprintf(LOG, "Error opening file %s for reading.\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL);
    destruct_heat();
    return -1;
  }
  if (!(box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for %s\nAborting...\n", filename);
    fprintf(LOG, "Error in memory allocation for %s\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL); fclose(F);   destruct_heat();
    return -1;
  }
  if (!(unfiltered_box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for %s\nAborting...\n", filename);
    fprintf(LOG, "Error in memory allocation for %s\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL);fclose(F);   destruct_heat(); fftwf_free(box);
    return -1;
  }
  fprintf(stderr, "Reading in deltax box\n");
  fprintf(LOG, "Reading in deltax box\n");
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
	if (fread((float *)unfiltered_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	  fprintf(stderr, "Error reading-in binary file %s\nAborting...\n", filename);
	  fprintf(LOG, "Error reading-in binary file %s\nAborting...\n", filename);
	  fftwf_free(box); fclose(GLOBAL_EVOL); fclose(F); fclose(LOG); fftwf_free(unfiltered_box);
	  destruct_heat();
	  return -1;
	}
      }
    }
  }
  fclose(F);


  /*** Transform unfiltered box to k-space to prepare for filtering ***/
  fprintf(stderr, "begin initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  fprintf(LOG, "begin initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)unfiltered_box, (fftwf_complex *)unfiltered_box, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
    unfiltered_box[ct] /= (float)HII_TOT_NUM_PIXELS;
  }
  fprintf(stderr, "end initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  fprintf(LOG, "end initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);


  /*** Create the z=0 non-linear density fields smoothed on scale R to be used in computing fcoll ***/
  R = L_FACTOR*BOX_LEN/(float)HII_DIM;
  R_factor = pow(R_XLy_MAX/R, 1/(float)NUM_FILTER_STEPS_FOR_Ts);
  //  R_factor = pow(E, log(HII_DIM)/(float)NUM_FILTER_STEPS_FOR_Ts);
  for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
    R_values[R_ct] = R;
    sigma_atR[R_ct] = sigma_z0(RtoM(R));
    fprintf(stderr, "Processing scale R= %06.2fMpc, time=%06.2f min\n", R, 
	    (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Processing scale R= %06.2fMpc, time=%06.2f min\n", R, 
	    (double)clock()/CLOCKS_PER_SEC/60.0);
    if (! (delNL0[R_ct] = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
      fprintf(stderr, "Error in memory allocation\nAborting...\n");
      fprintf(LOG, "Error in memory allocation\nAborting...\n");
      fclose(LOG); fclose(GLOBAL_EVOL);fftwf_free(box);  fftwf_free(unfiltered_box);
      for(ct=0; ct<R_ct; ct++)
	free(delNL0[ct]);
      destruct_heat();
      return -1;
    }

    // copy over unfiltered box
    memcpy(box, unfiltered_box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    if (R_ct > 0){ // don't filter on cell size
      HII_filter(box, HEAT_FILTER, R);
    }

    // now fft back to real space
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)box, (float *)box, FFTW_ESTIMATE);
    fftwf_execute(plan);

    // copy over the values
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  delNL0[R_ct][HII_R_INDEX(i,j,k)] = *((float *) box + HII_R_FFT_INDEX(i,j,k));
	  if (delNL0[R_ct][HII_R_INDEX(i,j,k)] < -1){ // correct for alliasing in the filtering step
	    delNL0[R_ct][HII_R_INDEX(i,j,k)] = -1+FRACT_FLOAT_ERR;
	  }

	  // and linearly extrapolate to z=0
	  delNL0[R_ct][HII_R_INDEX(i,j,k)] /= growth_factor_z; 
	}
      }
    }

    R *= R_factor;
  } //end for loop through the filter scales R
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  fftwf_free(box); fftwf_free(unfiltered_box);// we don't need this anymore

  // now lets allocate memory for our kinetic temperature and residual neutral fraction boxes
  if (!(Tk_box = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for Tk box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for Tk box\nAborting...\n");
    fclose(LOG);fclose(GLOBAL_EVOL);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }
  if (!(x_e_box = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for xe box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for xe box\nAborting...\n");
    fclose(LOG);  free(Tk_box);fclose(GLOBAL_EVOL);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }

  // and finally allocate memory for the spin temperature box
  if (!(Ts = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for Ts box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for Ts box\nAborting...\n");
    fclose(LOG);  fclose(GLOBAL_EVOL);free(Tk_box); free(x_e_box);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }


  // and initialize to the boundary values at Z_HEAT_END
  if (!RESTART){ // we are not restarting
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      Tk_box[ct] = Tk_BC;
      x_e_box[ct] = xe_BC;
    }
    x_e_ave = xe_BC;
    Tk_ave = Tk_BC;

    printf("Starting at at z_max=%f, Tk=%f, x_e=%e\n", Z_HEAT_MAX, Tk_ave, x_e_ave);
  }
  else{ // we need to load the evolution files from the intermediate output
    // first Tk
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
    if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "Ts.c: WARNING: Unable to open input file %s\nAborting\n", filename);
      fprintf(LOG, "Ts.c: WARNING: Unable to open input file %s\nAborting\n", filename);
      fclose(LOG); fclose(GLOBAL_EVOL); free(Tk_box); free(x_e_box); free(Ts);
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	free(delNL0[R_ct]);
      }
      destruct_heat();
      return -1;
    }
    else{
      if (mod_fread(Tk_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	fprintf(stderr, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
	fprintf(LOG, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
	fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts);
	for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	  free(delNL0[R_ct]);
	}
	destruct_heat();
      }
      fclose(F);
    }
    // then xe_neutral
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
      if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\nAborting\n", filename);
      fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\nAborting\n", filename);
      fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts);
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	free(delNL0[R_ct]);
      }
      destruct_heat();
    }
    else{
      if (mod_fread(x_e_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	fprintf(stderr, "Ts.c: Write error occured while reading xe box.\n");
	fprintf(LOG, "Ts.c: Write error occured while reading xe box.\n");
	fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts);
	for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	  free(delNL0[R_ct]);
	}
	destruct_heat();
      }
      fclose(F);
    }
    Tk_ave = x_e_ave = 0;
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      Tk_ave += Tk_box[box_ct];
      x_e_ave += x_e_box[box_ct];
    }
    Tk_ave /= (double) HII_TOT_NUM_PIXELS;
    x_e_ave /= (double) HII_TOT_NUM_PIXELS;
    fprintf(stderr, "Rebooting from z'=%f output. <Tk> = %f. <xe> = %e\n", zp, Tk_ave, x_e_ave);
    
  }

  /***************    END INITIALIZATION   *********************************/


  /*********  FOR DEBUGGING, set IGM to be homogeneous for testing purposes 
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++)
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS;box_ct++)
      delNL0[R_ct][box_ct] = 0;
    /*  *********/

  // main trapezoidal integral over z' (see eq. ? in Mesinger et al. 2009)
  if (!RESTART){
    zp = REDSHIFT*1.0001; //higher for rounding
    while (zp < Z_HEAT_MAX)
      zp = ((1+zp)*ZPRIME_STEP_FACTOR - 1);
    prev_zp = Z_HEAT_MAX;
  }
  else{
    prev_zp = zp;
  }
  zp = ((1+zp)/ ZPRIME_STEP_FACTOR - 1);
  dzp = zp - prev_zp;
  zp_ct=0;
  COMPUTE_Ts = 0;
  while (zp > REDSHIFT){
    // check if we will next compute the spin temperature (i.e. if this is the final zp step)
    if (Ts_verbose || (((1+zp) / ZPRIME_STEP_FACTOR) < (REDSHIFT+1)) )
      COMPUTE_Ts = 1;

    // check if we are in the really high z regime before the first stars..
    if (FgtrM(zp, FMAX(TtoM(zp, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM)) < 1e-15 )
      NO_LIGHT = 1;
    else
      NO_LIGHT = 0;

    M_MIN_at_zp = get_M_min_ion(zp);
    filling_factor_of_HI_zp = 1 - HII_EFF_FACTOR * FgtrM_st(zp, M_MIN_at_zp) / (1.0 - x_e_ave);
    if (filling_factor_of_HI_zp > 1) filling_factor_of_HI_zp=1;
    //    filling_factor_of_HI_zp = 1 - (1-nf_at_z)/FgtrM_st(REDSHIFT, M_MIN_at_z) * FgtrM_st(zp, M_MIN_at_zp);
    fprintf(stderr, "z'=%f, filling factor of HI is %e, without x_e would be %e, time=%06.2f min\n",
	    zp, filling_factor_of_HI_zp, 1 - HII_EFF_FACTOR * FgtrM_st(zp, M_MIN_at_zp), (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "z'=%f, filling factor of HI is %e, time=%06.2f min\n",
	    zp, filling_factor_of_HI_zp, (double)clock()/CLOCKS_PER_SEC/60.0);

    // let's initialize an array of redshifts (z'') corresponding to the 
    // far edge of the dz'' filtering shells
    // and the corresponding minimum halo scale, sigma_Tmin, 
    // as well as an array of the frequency integrals    
    fprintf(stderr, "Initializing look-up tables. Time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Initializing look-up tables. Time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
    time(&start_time);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){

      if (R_ct==0){
	prev_zpp = zp;
	prev_R = 0;
      }
      else{
	prev_zpp = zpp_edge[R_ct-1];
	prev_R = R_values[R_ct-1];
      }

      zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
      zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
      sigma_Tmin[R_ct] = sigma_z0(FMAX(TtoM(zpp, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM));

      // let's now normalize the total collapse fraction so that the mean is the
      // Sheth-Torman collapse fraction
      fcoll_R = 0;
      sample_ct=0;
      for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct+=(HII_TOT_NUM_PIXELS/1e5+1)){
	sample_ct++;
	fcoll_R += sigmaparam_FgtrM_bias(zpp, sigma_Tmin[R_ct], 
					 delNL0[R_ct][box_ct], sigma_atR[R_ct]);
      }
      fcoll_R /= (double) sample_ct;
      ST_over_PS[R_ct] = FgtrM_st(zpp, FMAX(TtoM(zpp, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM)) / fcoll_R;

      if (DEBUG_ON){
      printf("ST/PS=%g, mean_ST=%g, mean_ps=%g\n, ratios of mean=%g\n", ST_over_PS[R_ct], 
	     FgtrM_st(zpp, TtoM(zpp, X_RAY_Tvir_MIN, mu_for_Ts)), 
	     FgtrM(zpp, TtoM(zpp, X_RAY_Tvir_MIN, mu_for_Ts)),
	     FgtrM_st(zpp, TtoM(zpp, X_RAY_Tvir_MIN, mu_for_Ts))/FgtrM(zpp, TtoM(zpp, X_RAY_Tvir_MIN, mu_for_Ts))
	     );
      }


      lower_int_limit = FMAX(nu_tau_one(zp, zpp, x_e_ave, filling_factor_of_HI_zp), NU_X_THRESH);
      if (filling_factor_of_HI_zp < 0) filling_factor_of_HI_zp = 0; // for global evol; nu_tau_one above treats negative (post_reionization) inferred filling factors properly
/***************  PARALLELIZED LOOP ******************************************************************/
      // set up frequency integral table for later interpolation for the cell's x_e value
#pragma omp parallel shared(freq_int_heat_tbl, freq_int_ion_tbl, COMPUTE_Ts, freq_int_lya_tbl, zp, R_ct, x_e_ave, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, lower_int_limit) private(x_e_ct)
{
#pragma omp for
  for (x_e_ct = 0; x_e_ct < x_int_NXHII; x_e_ct++){
    freq_int_heat_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 0);
    freq_int_ion_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 1);
    if (COMPUTE_Ts)
      freq_int_lya_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 2);
  }
} // end omp declaration
/***************  END PARALLELIZED LOOP ******************************************************************/


      // and create the sum over Lya transitions from direct Lyn flux
      sum_lyn[R_ct] = 0;
      for (n_ct=NSPEC_MAX; n_ct>=2; n_ct--){
	if (zpp > zmax(zp, n_ct))
	  continue;

	nuprime = nu_n(n_ct)*(1+zpp)/(1.0+zp);
	sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0);
      }
    } // end loop over R_ct filter steps
    time(&curr_time);
    fprintf(stderr, "Finishing initializing look-up tables.  It took %06.2f min on the main thread. Time elapsed (total for all threads)=%06.2f\n", difftime(curr_time, start_time)/60.0, (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Finishing initializing look-up tables.  It took %06.2f min on the main thread. Total time elapsed (total for all threads)=%06.2f\n", difftime(curr_time, start_time)/60.0, (double)clock()/CLOCKS_PER_SEC/60.0);
    fflush(NULL);

    // scroll through each cell and update the temperature and residual ionization fraction
    growth_factor_zp = dicke(zp);
    dgrowth_factor_dzp = ddicke_dz(zp);
    dt_dzp = dtdz(zp);
    const_zp_prefactor = ZETA_X * X_RAY_SPEC_INDEX / NU_X_THRESH * C
      * F_STAR * OMb * RHOcrit * pow(CMperMPC, -3) * pow(1+zp, X_RAY_SPEC_INDEX+3);


    /********  LOOP THROUGH BOX *************/
    fprintf(stderr, "Looping through box at z'=%f, time elapsed  (total for all threads)= %06.2f min\n", zp, (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Looping through box at z'=%f, time elapsed  (total for all threads)= %06.2f min\n", zp, (double)clock()/CLOCKS_PER_SEC/60.0);
    fflush(NULL);
    time(&start_time);
    for (ct=0; ct<NUMCORES; ct++)
      J_alpha_threads[ct] = xalpha_threads[ct] = Xheat_threads[ct] = Xion_threads[ct] = 0;
    /***************  PARALLELIZED LOOP ******************************************************************/
#pragma omp parallel shared(COMPUTE_Ts, Tk_box, x_e_box, x_e_ave, delNL0, freq_int_heat_tbl, freq_int_ion_tbl, freq_int_lya_tbl, zp, dzp, Ts, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, growth_factor_zp, dgrowth_factor_dzp, NO_LIGHT, zpp_edge, sigma_atR, sigma_Tmin, ST_over_PS, sum_lyn, const_zp_prefactor, M_MIN_at_z, M_MIN_at_zp, dt_dzp, J_alpha_threads, xalpha_threads, Xheat_threads, Xion_threads) private(box_ct, ans, xHII_call, R_ct, curr_delNL0, m_xHII_low, m_xHII_high, freq_int_heat, freq_int_ion, freq_int_lya, dansdz, J_alpha_tot, curr_xalpha)
    {
#pragma omp for
      
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      if (!COMPUTE_Ts && (Tk_box[box_ct] > MAX_TK)) //just leave it alone and go to next value
	continue;

      // set to current values before updating
      ans[0] = x_e_box[box_ct];
      ans[1] = Tk_box[box_ct];

      /*
      if (DEBUG_ON){
	if (isnan(ans[0]))
	  fprintf(stderr, "Problem at cell %llu, x_e=%e\n", box_ct, ans[0]);
	if (isnan(ans[1]))
	  fprintf(stderr, "Problem at cell %llu, Tk=%e\n", box_ct, ans[1]);
      }
      */

      xHII_call = x_e_box[box_ct];

      // Check if ionized fraction is within boundaries; if not, adjust to be within
      if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
	xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
      } else if (xHII_call < x_int_XHII[0]) {
	xHII_call = 1.001*x_int_XHII[0];
      }
      //interpolate to correct nu integral value based on the cell's ionization state
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	curr_delNL0[R_ct] = delNL0[R_ct][box_ct];
	m_xHII_low = locate_xHII_index(xHII_call);
	m_xHII_high = m_xHII_low + 1;

	// heat
	freq_int_heat[R_ct] = (freq_int_heat_tbl[m_xHII_high][R_ct] - 
			       freq_int_heat_tbl[m_xHII_low][R_ct]) / 
	                  (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
	freq_int_heat[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
	freq_int_heat[R_ct] += freq_int_heat_tbl[m_xHII_low][R_ct];

	// ionization
	freq_int_ion[R_ct] = (freq_int_ion_tbl[m_xHII_high][R_ct] - 
			       freq_int_ion_tbl[m_xHII_low][R_ct]) / 
	                  (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
	freq_int_ion[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
	freq_int_ion[R_ct] += freq_int_ion_tbl[m_xHII_low][R_ct];

	// lya
	if (COMPUTE_Ts){
	  freq_int_lya[R_ct] = (freq_int_lya_tbl[m_xHII_high][R_ct] - 
				 freq_int_lya_tbl[m_xHII_low][R_ct]) / 
	    (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
	  freq_int_lya[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
	  freq_int_lya[R_ct] += freq_int_lya_tbl[m_xHII_low][R_ct];
	}
      }

      /********  finally compute the redshift derivatives *************/
      evolveInt(zp, curr_delNL0, freq_int_heat, freq_int_ion, freq_int_lya,
		COMPUTE_Ts, ans, dansdz);
 
      //update quantities
      x_e_box[box_ct] += dansdz[0] * dzp; // remember dzp is negative
      if (x_e_box[box_ct] > 1) // can do this late in evolution if dzp is too large
	x_e_box[box_ct] = 1 - FRACT_FLOAT_ERR;
      else if (x_e_box[box_ct] < 0)
	x_e_box[box_ct] = 0;
      if (Tk_box[box_ct] < MAX_TK)
	Tk_box[box_ct] += dansdz[1] * dzp;


      if (Tk_box[box_ct]<0){ // spurious bahaviour of the trapazoidalintegrator. generally overcooling in underdensities
	Tk_box[box_ct] = T_cmb*(1+zp);
      }
      if (COMPUTE_Ts){
	J_alpha_tot = dansdz[2]; //not really d/dz, but the lya flux
	Ts[box_ct] = get_Ts(zp, curr_delNL0[0]*growth_factor_zp,
			    Tk_box[box_ct], x_e_box[box_ct], J_alpha_tot, &curr_xalpha);
	J_alpha_threads[omp_get_thread_num()] += J_alpha_tot;
	xalpha_threads[omp_get_thread_num()] += curr_xalpha;
	Xheat_threads[omp_get_thread_num()] += dansdz[3];
	Xion_threads[omp_get_thread_num()] += dansdz[4];
      }
    }
    } // end parallelization pragma

/***************  END PARALLELIZED LOOP ******************************************************************/
    time(&curr_time);
    fprintf(stderr, "End scrolling through the box, which took %06.2 min\n", difftime(curr_time, start_time)/60.0);
    fprintf(LOG, "End scrolling through the box, which took %06.2 min\n", difftime(curr_time, start_time)/60.0);
    fflush(NULL);

    // compute new average values
    x_e_ave = 0; Tk_ave = 0; Ts_ave = 0; J_alpha_ave = 0; xalpha_ave = 0; Xheat_ave=0; Xion_ave=0;
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      x_e_ave += x_e_box[box_ct];
      Tk_ave += Tk_box[box_ct];
      if (COMPUTE_Ts)
	Ts_ave += Ts[box_ct];
    }
    for (ct=0; ct<NUMCORES; ct++){
      J_alpha_ave += J_alpha_threads[ct];
      xalpha_ave += xalpha_threads[ct];
      Xheat_ave += Xheat_threads[ct];
      Xion_ave += Xion_threads[ct];
    }
    Ts_ave /= (double)HII_TOT_NUM_PIXELS;
    x_e_ave /= (double)HII_TOT_NUM_PIXELS;
    Tk_ave /= (double)HII_TOT_NUM_PIXELS;
    J_alpha_ave /= (double)HII_TOT_NUM_PIXELS;
    xalpha_ave /= (double)HII_TOT_NUM_PIXELS;
    Xheat_ave /= (double)HII_TOT_NUM_PIXELS;
    Xion_ave /= (double)HII_TOT_NUM_PIXELS;
    // write to global evolution file
    fprintf(GLOBAL_EVOL, "%f\t%f\t%f\t%e\t%f\t%f\t%e\t%e\t%e\t%e\n", zp, filling_factor_of_HI_zp, Tk_ave, x_e_ave, Ts_ave, T_cmb*(1+zp), J_alpha_ave, xalpha_ave, Xheat_ave, Xion_ave);
    fflush(NULL);


    // output these intermediate boxes
    if ( Ts_verbose || (++zp_ct >= 10)){ // print every 10th z' evolution step, in case we need to restart
      zp_ct=0;
      fprintf(stderr, "Writting the intermediate output at zp = %.4f, <Tk>=%f, <x_e>=%e\n", zp, Tk_ave, x_e_ave);
      fprintf(LOG, "Writting the intermediate output at zp = %.4f, <Tk>=%f, <x_e>=%e\n", zp, Tk_ave, x_e_ave);
      fflush(NULL);

      // first Tk
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
      if (!(F=fopen(filename, "wb"))){
	fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
	fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
	if (mod_fwrite(Tk_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	}
	fclose(F);
      }
      // then xe_neutral
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
      if (!(F=fopen(filename, "wb"))){
	fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
	fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
	if (mod_fwrite(x_e_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	}
	fclose(F);
      }
    }

    // and the spin temperature if desired
  if ( COMPUTE_Ts ){
    sprintf(filename, "../Boxes/Ts_z%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN); 
      if (!(F=fopen(filename, "wb"))){
	fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
	fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
	if (mod_fwrite(Ts, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	}
	fclose(F);
      }
    }

    prev_zp = zp;
    zp = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
    dzp = zp - prev_zp;
  } // end main integral loop over z'



  //deallocate
  fclose(LOG); fclose(GLOBAL_EVOL); free(Tk_box); free(x_e_box); free(Ts);
  for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
    free(delNL0[R_ct]);
  }
  destruct_heat();
  return 0;
}


