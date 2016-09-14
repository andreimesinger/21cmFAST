#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"

/*
  Program DRIVE_ZSCROLL_REION_PARAM.C scrolls through the redshifts defined below,
  creating halo, evolved density, velocity, 21cm fields.
  NOTE: this driver assumes that the IGM has already been heated to Ts>>Tcmb.
  If you wish to compute the spin temperature, use the other driver.
*/

/*
  If KEEP_BOXES is set to 1, the code moves all ionization and 21cm boxes to a subdirectory in Boxes/
  If set to zero, the boxes are deleted.  If set to 2, then just the lightcone boxes are kept.
Note that HD space quickly gets used up if the
  boxes are kept and you go through many points in parameter space.
 */
#define KEEP_BOXES (int) 1

/*
  redshift scrolling parameters used by the drive_zscroll drive program
*/
#define ZSTART (float) 5.6 //inclusive
#define ZEND (float) 22.1 // inclusive
#define ZSTEP (float) 0.2
#define NUM_Z (int) ((ZEND-ZSTART)/(float)ZSTEP + 1.5)

/*  ZETA Parameter space */
#define ZETA_START (float) 10 //inclusive
#define ZETA_STOP (float) 100.1 // inclusive
#define ZETA_STEP (float) 20

/* Tvir min parameter space */
#define Tvir_START (float) 1e4 //inclusive (in K)
#define Tvir_STOP (float) 4e5 // inclusive (in K)
#define Tvir_FACTOR (float) 3.17

/* mfp parameter space */
#define MFP_START (float) 3.048
#define MFP_STOP (float) 60
#define MFP_FACTOR (float) 2.7


// returns the array index, i, such that arry[i-1] < value <= arry[i]
int find_index(float *arry, int len, float value){
  int i=1;
  if (value < arry[0])
    return 0;
  while ( (i < len) && (arry[i] < value) )
    i++;
  return i;
}

int main(int argc, char ** argv){
  float Z, M, M_MIN, TVIR, MFP, ZETA, zarry[NUM_Z], xHarry[NUM_Z], taue, z_re, del_z;
  float z1, z2, value, value1, value2,HII_RAM, TOT_RAM;
  char cmnd[1000], dirname[300], *token;
  FILE *LOG, *F;
  time_t start_time, curr_time;
  int status, i, index, index1, index2, astro_num_threads, th_index, tot_z, th_ct, proc_per_th;
  pid_t pid, pidarry[NUMCORES];


  /****** Initialization ******/
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
  system("mkdir ../Output_files/kSZ_power");
  system("mkdir ../Redshift_interpolate_filelists");
  system("mkdir ../Lighttravel_filelists");
  sprintf(cmnd, "rm ../Boxes/xH_*%i_%.0fMpc", HII_DIM, BOX_LEN);
  system(cmnd);

  // open log file
  // LOG = fopen("../Log_files/drive_zscroll_reion_param_log_file", "w");
  LOG = log_open("../Log_files/drive_zscroll_reion_param_log_file");
  if (!LOG){
    fprintf(stderr, "drive_zscroll_log_file.c: Unable to open log file\n Aborting...\n");
    return -1;
  }

  // calculate approximate RAM usage
  HII_RAM = sizeof(float)*pow(HII_DIM, 3)*5/1.0e9; // rough estimate of memory usage in GB
  TOT_RAM = HII_RAM*4/5.0 + sizeof(float)*pow(DIM, 3)/1.0e9; // rough estimate of total mem usage
  if (TOT_RAM > RAM){
    fprintf(stderr, "drive_zscroll_reion_param: WARNING: Your run will roughly use %.2fGB of RAM, and you claim to have %.2fGB.\n If this crashes, try lowering the resolution\n", TOT_RAM, RAM);
    fprintf(LOG, "drive_zscroll_reion_param: WARNING: Your run will roughly use %.2fGB of RAM, and you claim to have %.2fGB.\n If this crashes, try lowering the resolution\n", TOT_RAM, RAM);
  }

  fprintf(stderr, "Calling init to set up the initial conditions\n");
  fprintf(LOG, "Calling init to set up the initial conditions\n");
  system("init"); // you only need this call once per realization

  /**** First let's do the solely redshift-dependent code *****/

  Z = ZSTART;
  i=0;
  while (Z < (ZEND+0.0001)){
    // if USE_HALO_FIELD is turned on in ANAL_PARAMS.H, run the halo finder
    if (USE_HALO_FIELD){
      // find halos
      sprintf(cmnd, "find_halos %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fflush(NULL);
      system(cmnd);


      // shift halos accordig to their linear velocities
      sprintf(cmnd, "update_halo_pos %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fflush(NULL);
      system(cmnd);
    }

    // shift density field and update velocity field
    sprintf(cmnd, "perturb_field %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
    fflush(NULL);
    system(cmnd);

    // go to next redshift
    Z = ZSTART + (++i*ZSTEP);
  } // end of cosmological fields

  /**** We now have the PT density and velocity fields.  Let's create the redshift interpolated boxes for them.  ****/
  // create the lists of filenames
  sprintf(cmnd, "ls ../Boxes/updated_smoothed_deltax_*%i_%.0fMpc > ../Redshift_interpolate_filelists/delta_%i_%.0fMpc",
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);
  sprintf(cmnd, "ls ../Boxes/updated_vz_*%i_%.0fMpc > ../Redshift_interpolate_filelists/v_%i_%.0fMpc",
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);

  // first the density fields
  sprintf(cmnd, "redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/delta_%i_%.0fMpc", HII_DIM, BOX_LEN);
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);
  sprintf(cmnd, "ls ../Boxes/updated_smoothed_deltax_*%i_%.0fMpc_lighttravel > ../Lighttravel_filelists/delta_%i_%.0fMpc",
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);

  // then the velocity fields
  sprintf(cmnd, "redshift_interpolate_boxes 2 ../Redshift_interpolate_filelists/v_%i_%.0fMpc", HII_DIM, BOX_LEN);
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);
  sprintf(cmnd, "ls ../Boxes/updated_v*%i_%.0fMpc_lighttravel > ../Lighttravel_filelists/v_%i_%.0fMpc",
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);


  /**** start the astrophysical parameter loops  ****/
  // these we will thread depending on available memory
  astro_num_threads = FMIN(RAM/(float)HII_RAM, (float)NUMCORES) +0.5;
  if (astro_num_threads == 0) astro_num_threads=1;
  proc_per_th = NUMCORES/astro_num_threads + 1; //not really sure what this should be...
  fprintf(stderr, "Starting astro parameter runs.  Main loop will have %i threads, with each thread having max %i proc.\n", astro_num_threads, proc_per_th);
  fprintf(LOG, "Starting astro parameter runs.  Main loop will have %i threads, with each thread having max %i proc.\n", astro_num_threads, proc_per_th);

  MFP = MFP_START;
  while (MFP < MFP_STOP){

    ZETA = ZETA_START;
    while (ZETA < ZETA_STOP){

      TVIR = Tvir_START;
      while (TVIR < Tvir_STOP){


	// inner redhift loop
	Z = ZSTART;
	i=0;
	while (Z < (ZEND+0.0001)){
	  // thread the loop
	  th_index = 0;
	  while ((Z < (ZEND+0.0001)) && (th_index < astro_num_threads)){

	    switch (pid = fork()){

	      case -1: // fork has failed
		fprintf(stderr, "fork has failed\n");
		break;

	      case 0: // processed by child

		// compute M_MIN
		if (TVIR < 9.99999e3) // neutral IGM
		  M_MIN = TtoM(Z, TVIR, 1.22);
		else // ionized IGM
		  M_MIN = TtoM(Z, TVIR, 0.6);
		// check for WDM
		if (P_CUTOFF && ( M_MIN < M_J_WDM()))
		  M_MIN = M_J_WDM();

		// find bubbles
		sprintf(cmnd, "find_HII_bubbles -p %i %f %f %e %f", proc_per_th, Z,  ZETA, TVIR, MFP);

		printf("Thread %i now calling: %s\n", th_index, cmnd);
		system(cmnd);

		// do temperature map
		switch(FIND_BUBBLE_ALGORITHM){
		case 2:
		  if (USE_HALO_FIELD)
		    sprintf(cmnd, "delta_T -p %i %06.2f ../Boxes/xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", proc_per_th, Z, Z, ZETA, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
		  else
		    sprintf(cmnd, "delta_T -p %i %06.2f ../Boxes/xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", proc_per_th, Z, Z, ZETA, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
		  break;
		default:
		  if (USE_HALO_FIELD)
		    sprintf(cmnd, "delta_T -p %i %06.2f ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", proc_per_th, Z, Z, ZETA, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
		  else
		    sprintf(cmnd, "delta_T -p %i %06.2f ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", proc_per_th, Z, Z, ZETA, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
		}
		fprintf(stderr, "Thread %i now calling: %s\n", th_index, cmnd);
		system(cmnd);

		_exit(0);
		break;

	    default: // parent
	      pidarry[th_index] = pid;
	    }

	    th_index++;
	    Z = ZSTART + (++i*ZSTEP);
	  } // end of thread loop

	  // now wait for all of the children to finish
      	  for (th_ct=0; th_ct<th_index; th_ct++)
	    waitpid(pidarry[th_ct],&status,0);

	} // end main loop over redshift
	tot_z = i;

	/****************  Now let's create the redshift interpolated boxes for the xHI fields ****************/
	sprintf(cmnd, "ls ../Boxes/xH_*%i_%.0fMpc > ../Redshift_interpolate_filelists/xH_%i_%.0fMpc",
		HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
	system(cmnd);
	sprintf(cmnd, "redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/xH_%i_%.0fMpc", HII_DIM, BOX_LEN);
	system(cmnd);
	fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
	fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
	fflush(NULL);
	sprintf(cmnd, "ls ../Boxes/xH*%i_%.0fMpc_lighttravel > ../Lighttravel_filelists/xH_%i_%.0fMpc",
		HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
	system(cmnd);

	/*******  and redshift interpolate the 21cm boxes if u want... ******************/
	sprintf(cmnd, "ls ../Boxes/delta_T_*%i_%.0fMpc > ../Redshift_interpolate_filelists/delta_T_%i_%.0fMpc",
		HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
	system(cmnd);
	sprintf(cmnd, "redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/delta_T_%i_%.0fMpc", HII_DIM, BOX_LEN);
	system(cmnd);
	fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
	fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
	fflush(NULL);


	/***************** Get the secondary parameters ********************************************************/

	// first lets get a list of the xH files we just created
	sprintf(cmnd, "../Redshift_interpolate_filelists/xH_%i_%.0fMpc", HII_DIM, BOX_LEN);
	if (!(F=fopen(cmnd, "r"))){
	  fprintf(stderr, "drive_zscroll_reion_param: ERROR: unable to open %s for reading\nAborting\n.", cmnd);
	  fprintf(LOG, "drive_zscroll_reion_param: ERROR: unable to open %s for reading\nAborting\n.", cmnd);
	  fclose(LOG);
	  return -1;
	}

	// fill up the redshift and ionization arrays
	for (i=0; i < tot_z; i++){
	  fscanf(F, "%s\n", cmnd);
	  token = strtok(cmnd, "z");
	  token = strtok(NULL, "_");
	  zarry[i] = atof(token);
	  token = strtok(NULL, "f");
	  token = strtok(NULL, "_");
	  xHarry[i] = atof(token);
	}
	fclose(F);

	// now compute taue
	taue = tau_e(0, ZEND, zarry, xHarry, tot_z);

	// find z_re
	value = 0.5;
	index = find_index(xHarry, i, value);
	if ( (index == 0) || (index == i) ){ // outside of range
	  fprintf(stderr, "drive_zscroll_reion_param: WARNING: unable to compute z_re for parameter choice: (mfp, zeta, Tvir)=(%f, %f, %e)\nValue of xHI=0.5 is outside of redshift range.\nI will set z_re to -1. Consider narrowing parameter space or extending redshift range.\n", MFP, ZETA, TVIR);
	  fprintf(LOG, "drive_zscroll_reion_param: WARNING: unable to compute z_re for parameter choice: (mfp, zeta, Tvir)=(%f, %f, %e)\nValue of xHI=0.5 is outside of redshift range.\nI will set z_re to -1. Consider narrowing parameter space or extending redshift range.\n", MFP, ZETA, TVIR);
	  z_re = -1;
	}
	else{
	  z_re = zarry[index-1] + (value-xHarry[index-1])*
				   (zarry[index]-zarry[index-1]) / (xHarry[index] - xHarry[index-1]);
	  fprintf(stderr, "Interpolating between z(xH=%f)=%f and z(xH=%f)=%f to get z_re(xH=0.5)=%f\n",
		  xHarry[index-1], zarry[index-1], xHarry[index], zarry[index], z_re);
	  fprintf(LOG, "Interpolating between z(xH=%f)=%f and z(xH=%f)=%f to get z_re(xH=0.5)=%f\n",
		  xHarry[index-1], zarry[index-1], xHarry[index], zarry[index], z_re);
	}

	//delta z_re
	value1 = 0.25;
	index1 = find_index(xHarry, i, value1);
	value2 = 0.75;
	index2 = find_index(xHarry, i, value2);
	if ( (index1 == 0) || (index1 == i) || (index2 == 0) || (index2 == i)){ // outside of range
	  fprintf(stderr, "drive_zscroll_reion_param: WARNING: unable to compute  dz_re for parameter choice: (mfp, zeta, Tvir)=(%f, %f, %e)\nValues are outside of redshift range. I will set dz_re to -1. Consider narrowing parameter space or extending redshift range.\n", MFP, ZETA, TVIR);
	  fprintf(LOG, "drive_zscroll_reion_param: WARNING: unable to compute  dz_re for parameter choice: (mfp, zeta, Tvir)=(%f, %f, %e)\nValues are outside of redshift range.\nI will set dz_re to -1.\nConsider narrowing parameter space or extending redshift range.\n", MFP, ZETA, TVIR);
	  del_z = -1;
	}
	else{
	  z1 = zarry[index1-1] + (value1-xHarry[index1-1])*
				   (zarry[index1]-zarry[index1-1]) / (xHarry[index1] - xHarry[index1-1]);
	  z2 = zarry[index2-1] + (value2-xHarry[index2-1])*
				   (zarry[index2]-zarry[index2-1]) / (xHarry[index2] - xHarry[index2-1]);
	  del_z = z2-z1;
	}

	// and the associated directory filename
	sprintf(dirname, "Zeta%.1f_Tvir%.1e_mfp%.1f_Taue%.3f_zre%.3f_delz%.3f_%i_%.0fMpc",
		ZETA, TVIR, MFP, taue, z_re, del_z, HII_DIM, BOX_LEN);
	sprintf(cmnd, "mkdir ../Output_files/kSZ_power/%s", dirname);
	system(cmnd);
	sprintf(cmnd, "mkdir ../Output_files/Deldel_T_power_spec/%s", dirname);
	system(cmnd);
	sprintf(cmnd, "mv ../Output_files/Deldel_T_power_spec/ps* ../Output_files/Deldel_T_power_spec/%s", dirname);
	system(cmnd);


	/*********************  Create the kSZ power spectrum  ************************/

	sprintf(cmnd, "kSZ_power ../Lighttravel_filelists/delta_%i_%.0fMpc ../Lighttravel_filelists/xH_%i_%.0fMpc ../Lighttravel_filelists/v_%i_%.0fMpc ../Output_files/kSZ_power/%s/power", HII_DIM, BOX_LEN, HII_DIM, BOX_LEN, HII_DIM, BOX_LEN, dirname);
	fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
	fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
	fflush(NULL);
	system(cmnd);


	// print the reionization history
	sprintf(cmnd, "../Output_files/kSZ_power/%s/reion_hist", dirname);
	if (!(F=fopen(cmnd, "w"))){
	  fprintf(stderr, "drive_zscroll_reion_param: WARNING: unable to open %s for writing.\n", cmnd);
	  fprintf(LOG, "drive_zscroll_reion_param: WARNING: unable to open %s for writing.\n", cmnd);
	}
	for (i=0; i<tot_z; i++)
	  fprintf(F, "%.2f\t%f\n", zarry[i], xHarry[i]);
	fclose(F);
	sprintf(cmnd, "cp ../Output_files/kSZ_power/%s/reion_hist ../Output_files/Deldel_T_power_spec/%s", dirname, dirname);
	system(cmnd);


	// delete boxes
	if (KEEP_BOXES==0){
	  sprintf(cmnd, "rm ../Boxes/xH_*%i_%.0fMpc*", HII_DIM, BOX_LEN);
	  system(cmnd);
	  sprintf(cmnd, "rm ../Boxes/delta_T_*%i_%.0fMpc*", HII_DIM, BOX_LEN);
	  system(cmnd);
	}
	else if (KEEP_BOXES==2){
	  sprintf(cmnd, "rm ../Boxes/xH_*%i_%.0fMpc", HII_DIM, BOX_LEN);
	  system(cmnd);
	  sprintf(cmnd, "rm ../Boxes/delta_T_*%i_%.0fMpc", HII_DIM, BOX_LEN);
	  system(cmnd);
	}
	if (KEEP_BOXES != 0){
	  sprintf(cmnd, "mkdir ../Boxes/%s", dirname);
	  system(cmnd);
	  sprintf(cmnd, "mv ../Boxes/xH_*%i_%.0fMpc* ../Boxes/%s", HII_DIM, BOX_LEN, dirname);
	  system(cmnd);
	  sprintf(cmnd, "mv ../Boxes/delta_T_*%i_%.0fMpc* ../Boxes/%s", HII_DIM, BOX_LEN, dirname);
	  system(cmnd);
	}

	system("rm ../Log_files/HII*");
	system("rm ../Log_files/delta_T*");

	fprintf(stderr, "Finished processing %s\n", dirname);
	fprintf(LOG, "Finished processing %s\n", dirname);

	// move on in parameter space
	TVIR *= Tvir_FACTOR;
      } // end loop through Tvir
      ZETA += ZETA_STEP;
    } // end loop through zeta
    MFP *= MFP_FACTOR;
  } //end loop through mfp

  fclose(LOG);
  return 0;
}
