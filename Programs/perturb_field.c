#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "bubble_helper_progs.c"

/*
  USAGE: perturb_field <REDSHIFT>

  PROGRAM PERTURB_FIELD uses the first-order Langragian displacement field
  to move the masses in the cells of the density field.
  The high-res density field is extrapolated to some high-redshift
  (INITIAL_REDSHIFT in ANAL_PARAMS.H), then uses the zeldovich approximation
  to move the grid "particles" onto the lower-res grid we use for the
  maps.  Then we recalculate the velocity fields on the perturbed grid.


  Output files:

  "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN
  -- This file contains the perturbed overdensity field, \delta, at <REDSHIFT>. The binary box has FFT padding

  "../Boxes/updated_vx_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN
  "../Boxes/updated_vy_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN
  "../Boxes/updated_vz_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN
  -- These files contain the velocity fields recalculated using the perturbed velocity fields at <REDSHIFT>.  The units are cMpc/s.  The boxes have FFT padding.
*/


int print_box_no_padding(float *box, int d, FILE *F){
  int i,j,k;
  /*
  printf("%e ", *((float *)box+HII_R_FFT_INDEX(0,0,0)));
  printf("%e ", *((float *)box+HII_R_FFT_INDEX(0,10,0)));
  printf("%e ",*((float *)box+HII_R_FFT_INDEX(0,0,10)));
  printf("%e ", *((float *)box+HII_R_FFT_INDEX(10,0,0)));
  printf("%e\n", *((float *)box+HII_R_FFT_INDEX(HII_DIM-1,HII_DIM-1,HII_DIM-1)));
  */
  for(i=0; i<d; i++){
    for(j=0; j<d; j++){
      for(k=0; k<d; k++){
	if (fwrite(box + HII_R_FFT_INDEX(i, j, k), sizeof(float), 1, F)!=1)
	  return -1;
      }
    }
  }
  return 0;
}


int process_velocity(fftwf_complex *updated, float dDdt_over_D, float REDSHIFT, int component){
  char filename[300];
  FILE *F;
  float k_x, k_y, k_z, k_sq;
  int n_x, n_y, n_z;
  fftwf_plan plan;

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
	
	k_sq = k_x*k_x + k_y*k_y + k_z*k_z;

	// now set the velocities
	if ((n_x==0) && (n_y==0) && (n_z==0)) // DC mode
	  updated[0] = 0;
	else{
	  if (component == 0) // x-component
	    updated[HII_C_INDEX(n_x,n_y,n_z)] *= dDdt_over_D*k_x*I/k_sq/(HII_TOT_NUM_PIXELS+0.0);
	  else if (component == 1)
	    updated[HII_C_INDEX(n_x,n_y,n_z)] *= dDdt_over_D*k_y*I/k_sq/(HII_TOT_NUM_PIXELS+0.0);
	  else
	    updated[HII_C_INDEX(n_x,n_y,n_z)] *= dDdt_over_D*k_z*I/k_sq/(HII_TOT_NUM_PIXELS+0.0);
	}
      }
    }
    //    fprintf(stderr, "%i ", n_x);
  }
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)updated, (float *)updated, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  if (component == 0)
    sprintf(filename, "../Boxes/updated_vx_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  else if (component == 1)
    sprintf(filename, "../Boxes/updated_vy_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  else
    sprintf(filename, "../Boxes/updated_vz_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  if (!(F=fopen(filename, "wb"))){
    fprintf(stderr, "Unable to open file %s to write to.\n", filename);
    return -1;
  }
  if (print_box_no_padding((float *)updated, HII_DIM, F) < 0){
    fprintf(stderr, "perturb_field: Write error occured writting deltax box!\n");
    fclose(F);
    return -1;
  }
  fclose(F);
  return 0;
}

int main (int argc, char ** argv){
  char filename[100];
  FILE *F;
  fftwf_complex *updated, *save_updated;
  fftwf_plan plan;
  float *vx, *vy, *vz, REDSHIFT, growth_factor, displacement_factor_2LPT, init_growth_factor, init_displacement_factor_2LPT, xf, yf, zf, *vx_2LPT, *vy_2LPT, *vz_2LPT;
  float *deltax, mass_factor, dDdt, f_pixel_factor;
  unsigned long long ct, HII_i, HII_j, HII_k;
  int i,j,k, xi, yi, zi;
  double ave_delta, new_ave_delta;
  /***************   BEGIN INITIALIZATION   **************************/

  // check usage
  if (argc != 2){
    fprintf(stderr, "USAGE: perturb_field <REDSHIFT>\nAborting...\n");
    return -1;
  }
  REDSHIFT = atof(argv[1]);
  // initialize and allocate thread info
  if (fftwf_init_threads()==0){
    fprintf(stderr, "perturb_field: ERROR: problem initializing fftwf threads\nAborting\n.");
    return -1;
  }
  omp_set_num_threads(NUMCORES);

  // perform a very rudimentary check to see if we are underresolved and not using the linear approx
  if ((BOX_LEN > DIM) && !EVOLVE_DENSITY_LINEARLY){
    fprintf(stderr, "perturb_field.c: WARNING: Resolution is likely too low for accurate evolved density fields\n It Is recommended that you either increase the resolution (DIM/Box_LEN) or set the EVOLVE_DENSITY_LINEARLY flag to 1\n");
  }

  // initialize power spectrum 
  init_ps();
  growth_factor = dicke(REDSHIFT);
  displacement_factor_2LPT = -(3.0/7.0) * growth_factor*growth_factor; // 2LPT eq. D8

  //  fprintf(stderr, "gf = %.2e\ndf = %.2e\n", growth_factor, displacement_factor_2LPT);


  dDdt = ddickedt(REDSHIFT); // time derivative of the growth factor (1/s)
  init_growth_factor = dicke(INITIAL_REDSHIFT);
  init_displacement_factor_2LPT = -(3.0/7.0) * init_growth_factor*init_growth_factor; // 2LPT eq. D8

  // allocate memory for the updated density, and initialize
  updated = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!updated){
    fprintf(stderr, "perturb_field: Error allocating memory for box.\nAborting...\n");
    free_ps(); return -1;
  }
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){ updated[ct]=0;}

  vx = (float *) fftwf_malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  if (!vx){
    fprintf(stderr, "perturb_field: Error allocating memory for velocity box\nAborting...\n");
    free_ps(); return -1;
  }

  // check if the linear evolution flag was set
  if (EVOLVE_DENSITY_LINEARLY){
    sprintf(filename, "../Boxes/smoothed_deltax_z0.00_%i_%.0fMpc", HII_DIM, BOX_LEN);
    if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "perturb_field.c: Unable to open file %s for reading.\nAborting\n", filename);
      fftwf_free(updated); fftwf_free(vx);
      free_ps(); return -1;
    }
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if (fread(((float *)updated + HII_R_FFT_INDEX(i,j,k)), sizeof(float), 1, F) != 1){
	    fprintf(stderr, "perturb_field.c: Error reading file %s.\nAborting\n", filename);
	    fftwf_free(updated); fclose(F); fftwf_free(vx);
	    free_ps(); return -1;
	  }

	  *((float *)updated + HII_R_FFT_INDEX(i,j,k)) *= growth_factor;
	}
      }
    }
    fclose(F);
   }

  // first order Zel'Dovich perturbation
  else{

    fprintf(stderr, "Openning velocity files\n");
    // allocate memory for the velocity boxes and read them in
    vy = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vy){
      fprintf(stderr, "perturb_field: Error allocating memory for velocity box\nAborting...\n");
      fftwf_free(vx); fftwf_free(updated);
      free_ps(); return -1;
    }
    vz = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vz){
      fprintf(stderr, "perturb_field: Error allocating memory for velocity box\nAborting...\n");
      fftwf_free(vx); fftwf_free(vy); fftwf_free(updated);
      free_ps(); return -1;
    }
    sprintf(filename, "../Boxes/vxoverddot_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
    if (mod_fread(vx, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "perturb_field: Read error occured while reading velocity box.\n");
      fftwf_free(vx);  fftwf_free(vy); fftwf_free(vz); fftwf_free(updated);
      free_ps(); return -1;
    }
    fclose(F);
    sprintf(filename, "../Boxes/vyoverddot_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
    if (mod_fread(vy, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "perturb_field: Read error occured while reading velocity box.\n");
      fftwf_free(vx);  fftwf_free(vy); fftwf_free(vz);fftwf_free(updated);
      free_ps(); return -1;
    }
    fclose(F);
    sprintf(filename, "../Boxes/vzoverddot_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
    if (mod_fread(vz, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "perturb_field: Read error occured while reading velocity box.\n");
      fftwf_free(vx);  fftwf_free(vy); fftwf_free(vz);fftwf_free(updated);
      free_ps(); return -1;
    }
    fclose(F);
    // now add the missing factor of D
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      vx[ct] *= (growth_factor-init_growth_factor) / BOX_LEN; // this is now comoving displacement in units of box size
      vy[ct] *= (growth_factor-init_growth_factor) / BOX_LEN; // this is now comoving displacement in units of box size
      vz[ct] *= (growth_factor-init_growth_factor) / BOX_LEN; // this is now comoving displacement in units of box size
    }

    
    // read in the linear density field
    deltax = (float *) fftwf_malloc(sizeof(float)*TOT_FFT_NUM_PIXELS);
    if (!deltax){
      fprintf(stderr, "perturb_field.c: Error allocating memory for box.\nAborting...\n");
      fftwf_free(vx);  fftwf_free(vy); fftwf_free(vz);fftwf_free(updated);
      free_ps(); return -1;
    }
    sprintf(filename, "../Boxes/deltax_z0.00_%i_%.0fMpc", DIM, BOX_LEN);
    F = fopen(filename, "rb");
    fprintf(stderr, "Reading in deltax box\n");
    if (mod_fread(deltax, sizeof(float)*TOT_FFT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "perturb_field: Read error occured while reading deltax box.\n");
      fclose(F); fftwf_free(vx);  fftwf_free(vy); fftwf_free(vz); fftwf_free(deltax);fftwf_free(updated);
      free_ps(); return -1;
    }
    fclose(F);


    // find factor of HII pixel size / deltax pixel size
    f_pixel_factor = DIM/(float)HII_DIM;
    mass_factor = pow(f_pixel_factor, 3);


/* ************************************************************************* *
 *                           BEGIN 2LPT PART                                 *
 * ************************************************************************* */
// reference: reference: Scoccimarro R., 1998, MNRAS, 299, 1097-1118 Appendix D
  if(SECOND_ORDER_LPT_CORRECTIONS){
      
    //fprintf(stderr, "Begin initialization 2LPT velocity field\nTotal elapsed time: %ds\n", time(NULL) - start_time);
    //last_time = time(NULL);
    // allocate memory for the velocity boxes and read them in
    vx_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vx_2LPT){
      fprintf(stderr, "perturb_field: Error allocating memory for 2LPT velocity box\nAborting...\n");
      free(vx);  free(vy); free(vz);
      return -1;
     }
    vy_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vy_2LPT){
      fprintf(stderr, "perturb_field: Error allocating memory for 2LPT velocity box\nAborting...\n");
      free(vx_2LPT);
      free(vx);  free(vy); free(vz);
      return -1;
    }
    vz_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vz_2LPT){
      fprintf(stderr, "perturb_field: Error allocating memory for 2LPT velocity box\nAborting...\n");
      free(vx_2LPT); free(vy_2LPT);
      free(vx);  free(vy); free(vz);
      return -1;
    }

    // read again velocities

    sprintf(filename, "../Boxes/vxoverddot_2LPT_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
    if (mod_fread(vx_2LPT, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "perturb_field: Read error occured while reading velocity 2LPT box.\n");
      free(vx);  free(vy); free(vz);
      free(vx_2LPT);  free(vy_2LPT); free(vz_2LPT);
      return -1;
    }
    fclose(F);
    //fprintf(stderr, "Read 2LPT vx velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
    //last_time = time(NULL);

    sprintf(filename, "../Boxes/vyoverddot_2LPT_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
    if (mod_fread(vy_2LPT, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "perturb_field: Read error occured while reading velocity 2LPT box.\n");
      free(vx);  free(vy); free(vz);
      free(vx_2LPT);  free(vy_2LPT); free(vz_2LPT);
      return -1;
    }
    fclose(F);
    //fprintf(stderr, "Read 2LPT vy velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
    //last_time = time(NULL);

    sprintf(filename, "../Boxes/vzoverddot_2LPT_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
   if (mod_fread(vz_2LPT, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "perturb_field: Read error occured while reading velocity box.\n");
      free(vx);  free(vy); free(vz);
      free(vx_2LPT);  free(vy_2LPT); free(vz_2LPT);
      return -1;
    }
    fclose(F);
    //fprintf(stderr, "Read 2LPT vz velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
    //last_time = time(NULL);

    // now add the missing factor in eq. D9
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      vx_2LPT[ct] *= (displacement_factor_2LPT - init_displacement_factor_2LPT) / BOX_LEN; // this is now comoving displacement in units of box size
      vy_2LPT[ct] *= (displacement_factor_2LPT - init_displacement_factor_2LPT) / BOX_LEN; // this is now comoving displacement in units of box size
      vz_2LPT[ct] *= (displacement_factor_2LPT - init_displacement_factor_2LPT) / BOX_LEN; // this is now comoving displacement in units of box size
    }
    //fprintf(stderr, "Read 2LPT velocity field\nTotal time: %ds\n", time(NULL) - start_time);
  }

/* ************************************************************************* *
 *                            END 2LPT PART                                  *
 * ************************************************************************* */


    /************  END INITIALIZATION ****************************/


    // go through the high-res box, mapping the mass onto the low-res (updated) box
    for (i=0; i<DIM;i++){
      for (j=0; j<DIM;j++){
	for (k=0; k<DIM;k++){

	  // map indeces to locations in units of box size
	  xf = (i+0.5)/(DIM+0.0);
	  yf = (j+0.5)/(DIM+0.0);
	  zf = (k+0.5)/(DIM+0.0);

	  // update locations
	  HII_i = (unsigned long long)(i/f_pixel_factor);
	  HII_j = (unsigned long long)(j/f_pixel_factor);
	  HII_k = (unsigned long long)(k/f_pixel_factor);
	  xf += vx[HII_R_INDEX(HII_i, HII_j, HII_k)];
	  yf += vy[HII_R_INDEX(HII_i, HII_j, HII_k)];
	  zf += vz[HII_R_INDEX(HII_i, HII_j, HII_k)];


    // 2LPT PART
    // add second order corrections
    if(SECOND_ORDER_LPT_CORRECTIONS){
      xf -= vx_2LPT[HII_R_INDEX(HII_i,HII_j,HII_k)];
      yf -= vy_2LPT[HII_R_INDEX(HII_i,HII_j,HII_k)];
      zf -= vz_2LPT[HII_R_INDEX(HII_i,HII_j,HII_k)];
    }

	  xf *= HII_DIM;
	  yf *= HII_DIM;
	  zf *= HII_DIM;
	  while (xf >= (float)HII_DIM){ xf -= HII_DIM;}
	  while (xf < 0){ xf += HII_DIM;}
	  while (yf >= (float)HII_DIM){ yf -= HII_DIM;}
	  while (yf < 0){ yf += HII_DIM;}
	  while (zf >= (float)HII_DIM){ zf -= HII_DIM;}
	  while (zf < 0){ zf += HII_DIM;}
	  xi = xf;
	  yi = yf;
	  zi = zf;
	  if (xi >= HII_DIM){ xi -= HII_DIM;}
	  if (xi < 0) {xi += HII_DIM;}
	  if (yi >= HII_DIM){ yi -= HII_DIM;}
	  if (yi < 0) {yi += HII_DIM;}
	  if (zi >= HII_DIM){ zi -= HII_DIM;}
	  if (zi < 0) {zi += HII_DIM;}

	  // now move the mass
	  *( (float *)updated + HII_R_FFT_INDEX(xi, yi, zi) ) +=
	    (1 + init_growth_factor*deltax[R_FFT_INDEX(i,j,k)]);
	}
      }
    }

    // renormalize to the new pixel size, and make into delta
    //    ave_delta = 0;
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  *((float *)updated + HII_R_FFT_INDEX(i,j,k) ) /= mass_factor;
	  *((float *)updated + HII_R_FFT_INDEX(i,j,k) ) -= 1;
	  //  ave_delta += *((float *)updated + HII_R_FFT_INDEX(i,j,k) );
	}
      }
    }
    //    ave_delta /= (double)HII_TOT_NUM_PIXELS;
    //    fprintf(stderr, "ave is %e\n", ave_delta);
 
    // deallocate
    fftwf_free(vy); fftwf_free(vz); fftwf_free(deltax);
  }


  /****  Print and convert to velocities *****/
  fprintf(stderr, "Done with PT. Printing density field and computing velocity components.\n");
  fftwf_plan_with_nthreads(NUMCORES); // use all processors for perturb_field
  save_updated = (fftwf_complex *) vx;
  sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  F=fopen(filename, "wb");
  if (EVOLVE_DENSITY_LINEARLY){
    if (print_box_no_padding((float *)updated, HII_DIM, F) < 0){
      fprintf(stderr, "perturb_field: Write error occured writting deltax box!\n");
      fftwf_free(updated); fftwf_free(vx); fclose(F);
      free_ps(); return -1;
    }

    // transform to k-space
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)updated, (fftwf_complex *)updated, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();

    // save a copy of the k-space density field
    memcpy(save_updated, updated, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  }
  else{
    // transform to k-space
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)updated, (fftwf_complex *)updated, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();

    //smooth the field
    if (!EVOLVE_DENSITY_LINEARLY && SMOOTH_EVOLVED_DENSITY_FIELD){
      HII_filter(updated, 2, R_smooth_density*BOX_LEN/(float)HII_DIM);
    }

    // save a copy of the k-space density field
    memcpy(save_updated, updated, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)updated, (float *)updated, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();

    // normalize after FFT
    for(i=0; i<HII_DIM; i++){
      for(j=0; j<HII_DIM; j++){
	for(k=0; k<HII_DIM; k++){
	  *((float *)updated + HII_R_FFT_INDEX(i,j,k)) /= (float)HII_TOT_NUM_PIXELS;
	  if (*((float *)updated + HII_R_FFT_INDEX(i,j,k)) < -1) // shouldn't happen
	    *((float *)updated + HII_R_FFT_INDEX(i,j,k)) = -1+FRACT_FLOAT_ERR;
	}
      }
    }

    if (print_box_no_padding((float *)updated, HII_DIM, F) < 0){
      fprintf(stderr, "perturb_field: Write error occured writting deltax box!\n");
      fftwf_free(updated); fftwf_free(vx); fclose(F);
      free_ps(); return -1;
    }
 
    memcpy(updated, save_updated, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  }
  fclose(F);

  // x-component
  fprintf(stderr, "Generate x-component\n");
  if (process_velocity(updated, dDdt/growth_factor, REDSHIFT, 0) < 0){
    fftwf_free(updated); fftwf_free(vx); fftwf_cleanup_threads(); free_ps(); return 0;
  }

  // y-component
  fprintf(stderr, "Generate y-component\n");
  memcpy(updated, save_updated, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (process_velocity(updated, dDdt/growth_factor, REDSHIFT, 1) < 0){
    fftwf_free(updated); fftwf_free(vx); fftwf_cleanup_threads(); free_ps(); return 0;
  }
  // z-component
  fprintf(stderr, "Generate z-component\n");
  memcpy(updated, save_updated, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  process_velocity(updated, dDdt/growth_factor, REDSHIFT, 2);


  // deallocate
  fftwf_free(updated);
  fftwf_free(vx);
  fftwf_cleanup_threads();
  free_ps(); return 0;
}
