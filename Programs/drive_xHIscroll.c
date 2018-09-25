#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"

/*
  Program DRIVE_xHIscroll.C scrolls through various values of ionizing
  efficiencies defined below, creating halo, evolved density, velocity,
  21cm fields, for various values of the global ionized fraction, at a given
  redshift.
  NOTE: this driver assumes that the IGM has already been heated to Ts>>Tcmb.
  If you wish to compute the spin temperature, use the other driver.

  USAGE: drive_xHIscroll [FLAG]
  setting the optional argument FLAG to 1 bypasses the ionizing efficiency
  parameters below, and instead tries to create fields at <xHI> = 0.1, 0.2,
  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
*/

#define Z (float) 10 // redshift

/*
  redshift scrolling parameters used by the drive_zscroll drive program
*/
#define EFF_START (float) 5 //inclusive
#define EFF_STOP (float) 15.1 // inclusive
#define EFF_STEP (float) 2.5


int main(int argc, char ** argv){
  float ion_eff, M, M_MIN, fcoll, x_i;
  char cmnd[1000];
  FILE *LOG;
  time_t start_time, curr_time;

  time(&start_time);

  // make appropriate directories
  system("mkdir ../Log_files");
  system("mkdir ../Boxes");
  system("mkdir ../Output_files");
  system("mkdir ../Output_files/DNDLNM_files");
  system("mkdir ../Output_files/FgtrM_files");
  system("mkdir ../Output_files/Halo_lists");
  system("mkdir ../Output_files/Size_distributions");
  system("mkdir ../Output_files/Deldel_T_power_spec");



  // open log file
  system("rm ../Log_files/*");
  LOG = log_open("../Log_files/drive_zscroll_noTs_log_file");
  if (!LOG){
    fprintf(stderr, "drive_zscroll_log_file.c: Unable to open log file\n Aborting...\n");
    return -1;
  }

  fprintf(stderr, "Calling init to set up the initial conditions\n");
  fprintf(LOG, "Calling init to set up the initial conditions\n");
  system("./init"); // you only need this call once per realization

  fprintf(stderr, "*************************************\n");

  //set the minimum source mass
  init_ps();
  M_MIN = get_M_min_ion(Z);
  fcoll = FgtrM_st(Z, M_MIN);
  fprintf(stderr, "Collapsed fraction above M_halo=%e is %e\n", M_MIN, fcoll);
  fprintf(LOG, "Collapsed fraction above M_halo=%e is %e\n", M_MIN, fcoll);

  // if USE_HALO_FIELD is turned on in ANAL_PARAMS.H, run the halo finder
  if (USE_HALO_FIELD){
    //  the following only depend on redshift, not ionization field
    // find halos
    sprintf(cmnd, "./find_halos %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);


    // shift halos accordig to their linear velocities
    sprintf(cmnd, "./update_halo_pos %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
  }

  // shift density field and update velocity field
  sprintf(cmnd, "./perturb_field %.2f", Z);
  time(&curr_time);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
  fflush(NULL);
  system(cmnd);
  // end of solely redshift dependent things, now do ionization stuff

  if ((argc==2) && (atoi(argv[1]) == 1)){
    x_i = 0.1;
    ion_eff = x_i/fcoll;
  }
  else{
    ion_eff = EFF_START;
  }
  while (1){
    // find bubbles
    sprintf(cmnd, "./find_HII_bubbles %.2f %f", Z, ion_eff);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);


    // generate size distributions, first ionized bubbles
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	  // New in v1.4
	  // These filenames are just a stopgap to avoid error. I have to modify this part.
	//sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff,  EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff,  HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
      else
	//sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	//sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/sphere_xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/sphere_xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
      else
	//sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/sphere_xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/sphere_xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);


    // generate size distributions, then neutral regions
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	  // New in v1.4
	  // These filenames are just a stopgap to avoid error. I have to modify this part.
	//sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
      else
	//sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	//sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/sphere_xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/sphere_xH_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
      else
	//sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/sphere_xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/sphere_xH_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_z%06.2f_%i_%.0fMpc", Z, ion_eff, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, Z, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);


    // do temperature map
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	  // New in v1.4
	  // These filenames are just a stopgap to avoid error. I have to modify this part.
	//sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc",  Z, Z, ion_eff, EFF_FACTOR_PL_INDEX,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc",  Z, Z, ion_eff, HII_EFF_FACTOR,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	//sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z, ion_eff, EFF_FACTOR_PL_INDEX,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z, ion_eff, HII_EFF_FACTOR,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	//sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z, ion_eff, EFF_FACTOR_PL_INDEX,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z, ion_eff, HII_EFF_FACTOR,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	//sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z, ion_eff, EFF_FACTOR_PL_INDEX,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z, ion_eff, HII_EFF_FACTOR,HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);

    fprintf(stderr, "*************************************\n");
    fflush(NULL);

    // go to next ioninzing efficiency
    if ((argc==2) && (atoi(argv[1]) == 1)){
      x_i += 0.1;
      ion_eff = x_i/fcoll;
      if (x_i > 0.99)
	break;
    }
    else{
      ion_eff += EFF_STEP;
      if (ion_eff > EFF_STOP)
	break;
    }
  }

  fclose(LOG);
  return 0;
}
