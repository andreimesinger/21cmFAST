#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <sys/wait.h>

#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"

/*
  Program DRIVE_ZSCROLL.C scrolls through the redshifts defined in ANAL_PARAMS.H creating halo, velocity, density, and ionization fields
*/

#define ZLOW (float) (6)
#define ZHIGH  Z_HEAT_MAX

int main(int argc, char ** argv){
  //float Z, M, M_MIN, nf;
  float Z, M, nf;
  char cmnd[1000];
  FILE *LOG;
  time_t start_time, curr_time;
  int status;



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
  system("mkdir ../Redshift_interpolate_filelists");
  //  system("mkdir ../Lighttravel_filelists");

  // remove some of the previous (astro) files which might conflict with current run
  system("rm ../Boxes/Ts_evolution/*");
  system("rm ../Boxes/Ts_*");
  system("rm ../Boxes/delta_T_*");
  system("rm ../Boxes/xH_*");
  system("rm ../Boxes/Nrec_*");
  system("rm ../Boxes/z_first*");
  system("rm ../Output_files/Deldel_T_power_spec/*");  

  init_ps();

  // open log file
  system("rm ../Log_files/*");
  LOG = log_open("../Log_files/drive_logzscroll_log_file");
  if (!LOG){
    fprintf(stderr, "drive_zscroll_log_file.c: Unable to open log file\n Aborting...\n");
    return -1;
  }


  fprintf(stderr, "Calling init to set up the initial conditions\n");
  fprintf(LOG, "Calling init to set up the initial conditions\n");
  system("./init"); // you only need this call once per realization

  Z = ZLOW*1.0001; // match rounding convention from Ts.c

   // call Ts on the lowest redshift
  if (USE_TS_IN_21CM){
    sprintf(cmnd, "./perturb_field %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
    
    sprintf(cmnd, "./Ts %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
  }


  // now go to highest redshift and step downwards
  while (Z < ZHIGH){
    Z = ((1+Z)*ZPRIME_STEP_FACTOR - 1);
  }
  Z = ((1+Z)/ ZPRIME_STEP_FACTOR - 1);
  while (Z >= ZLOW){

    //set the minimum source mass
    //M_MIN = get_M_min_ion(Z);

    // if USE_HALO_FIELD is turned on in ANAL_PARAMS.H, run the halo finder
    if (USE_HALO_FIELD){
      //  the following only depend on redshift, not ionization field
      // find halos
      sprintf(cmnd, "./find_halos %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fflush(NULL);
      system(cmnd);


      // shift halos accordig to their linear velocities
      sprintf(cmnd, "./update_halo_pos %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
    }

    // shift density field and update velocity field
    sprintf(cmnd, "./perturb_field %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
    // end of solely redshift dependent things, now do ionization stuff


    // if it is the lowest redshift, let's call Ts.c
    if (USE_TS_IN_21CM && (Z > Z_HEAT_MAX) ) { // NEW CONDITIONAL
      //    if (USE_TS_IN_21CM && (fabs(Z-ZLOW)/Z < 0.0002) ){
      sprintf(cmnd, "./Ts %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
    } // this will create all of the higher z Ts files in Boxes, provided Ts_verbose is turned on
    // in HEAT_PARAMS.H


    // find bubbles
    if (INHOMO_RECO)
      sprintf(cmnd, "./find_HII_bubbles %f %f", Z, (1+Z)*ZPRIME_STEP_FACTOR - 1 );
    else
      sprintf(cmnd, "./find_HII_bubbles %f", Z );
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    status = system(cmnd);
    nf = WEXITSTATUS(status) / 100.0;
    if (nf < 0){
      fprintf(stderr, "find_HII_bubbles exited...\nAborting run...\n");
      fprintf(LOG,  "find_HII_bubbles exited...\nAborting run...\n");
      return -1;
    }


    // do temperature map
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      else
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_nohalos_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      else
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      break;
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);

    fprintf(stderr, "*************************************\n");
    fflush(NULL);

    // update the redshift value according to the logarithmic stepping in the Ts.c routine
    Z = ((1+Z)/ZPRIME_STEP_FACTOR - 1);
  }

  
  // Create lightcone boxes from the coeval cubes
  sprintf(cmnd, "ls ../Boxes/xH_*%i_%.0fMpc > ../Redshift_interpolate_filelists/xH_%i_%.0fMpc", 
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);
  sprintf(cmnd, "./redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/xH_%i_%.0fMpc", HII_DIM, BOX_LEN);
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);

  sprintf(cmnd, "ls ../Boxes/delta_T_*%i_%.0fMpc > ../Redshift_interpolate_filelists/delta_T_%i_%.0fMpc", 
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);
  sprintf(cmnd, "./redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/delta_T_%i_%.0fMpc", HII_DIM, BOX_LEN);
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);

  if (INHOMO_RECO){
    sprintf(cmnd, "ls ../Boxes/Nrec_*%i_%.0fMpc > ../Redshift_interpolate_filelists/Nrec_%i_%.0fMpc", 
	    HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
    system(cmnd);
    sprintf(cmnd, "./redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/Nrec_%i_%.0fMpc", HII_DIM, BOX_LEN);
    system(cmnd);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
    fflush(NULL);
  }

  sprintf(cmnd, "./extract_delTps.pl 0.1 ../Output_files/Deldel_T_power_spec/ps_z0* > Power_k0.1");
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);

  fclose(LOG);
  return 0;
}
