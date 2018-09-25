#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"

/********************************************************************
USAGE:  update_halo_pos <REDSHIFT>

Program UPDATE_HALO_POS reads in the linear velocity field, and uses
it to update halo locations with a corresponding displacement field
creating updated_halo_ in ../Output_files/Halo_lists/ directory
*********************************************************************/

float max(float a, float b){
  if (a > b) return a;
  return b;
}

int main(int argc, char ** argv){
  char filename[100];
  FILE *F, *OUT;
  float growth_factor, displacement_factor_2LPT, REDSHIFT, mass, xf, yf, zf, *vx, *vy, *vz, *vx_2LPT, *vy_2LPT, *vz_2LPT, z;
  int i,j,k, xi, yi, zi, DI;
  unsigned long long ct;
  float dz = 1e-10;
  time_t start_time, last_time;

  /******************   BEGIN INITIALIZATION     ********************************/
  // check arguments
  if (argc != 2){
    fprintf(stderr, "USAGE: update_halo_pos <redshift>\nAborting...\n");
    return -1;
  }
  REDSHIFT = atof(argv[1]);

  // initialize power spectrum 
  init_ps(0, 1e10);
  growth_factor = dicke(REDSHIFT);
  displacement_factor_2LPT = -(3.0/7.0) * growth_factor*growth_factor; // 2LPT eq. D8

  fprintf(stderr, "gf = %.2e\ndf = %.2e\n", growth_factor, displacement_factor_2LPT);

  fprintf(stderr, "Begin initialization velocity field\n");
  start_time = time(NULL);

  // allocate memory for the velocity boxes and read them in
  vx = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!vx){
    fprintf(stderr, "update_halo_pos: Error allocating memory for velocity box\nAborting...\n");
    return -1;
  }
  vy = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!vy){
    fprintf(stderr, "update_halo_pos: Error allocating memory for velocity box\nAborting...\n");
    free(vx);
    return -1;
  }
  vz = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!vz){
    fprintf(stderr, "update_halo_pos: Error allocating memory for velocity box\nAborting...\n");
    free(vx); free(vy);
    return -1;
  }
  sprintf(filename, "../Boxes/vxoverddot_%i_%.0fMpc", HII_DIM, BOX_LEN);
  F=fopen(filename, "rb");
  if (mod_fread(vx, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    fprintf(stderr, "update_halo_pos: Read error occured while reading velocity box.\n");
    free(vx);  free(vy); free(vz);
    return -1;
  }
  fclose(F);
  last_time = time(NULL);
  fprintf(stderr, "Read vx velocity field\nElapsed time: %ds\n", last_time - start_time);
  sprintf(filename, "../Boxes/vyoverddot_%i_%.0fMpc", HII_DIM, BOX_LEN);
  F=fopen(filename, "rb");
  if (mod_fread(vy, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    fprintf(stderr, "update_halo_pos: Read error occured while reading velocity box.\n");
    free(vx);  free(vy); free(vz);
    return -1;
  }
  fclose(F);
  fprintf(stderr, "Read vy velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
  last_time = time(NULL);
  sprintf(filename, "../Boxes/vzoverddot_%i_%.0fMpc", HII_DIM, BOX_LEN);
  F=fopen(filename, "rb");
  if (mod_fread(vz, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    fprintf(stderr, "update_halo_pos: Read error occured while reading velocity box.\n");
    free(vx);  free(vy); free(vz);
    return -1;
  }
  fclose(F);
  fprintf(stderr, "Read vz velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
  last_time = time(NULL);
  
  fprintf(stderr, "Read velocity field\n");
  // now add the missing factor of Ddot
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
    vx[ct] *= growth_factor / BOX_LEN; // this is now comoving displacement in units of box size
    vy[ct] *= growth_factor / BOX_LEN; // this is now comoving displacement in units of box size
    vz[ct] *= growth_factor / BOX_LEN; // this is now comoving displacement in units of box size
  }
  // open file to write to
  sprintf(filename, "../Output_files/Halo_lists/updated_halos_z%06.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
    OUT = fopen(filename, "w");
  if (!OUT){
    fprintf(stderr, "update_halo_pos: Error opening output file: %s\nAborting\n", filename);
    free(vx);  free(vy); free(vz);
    return -1;
  }

/* ************************************************************************* *
 *                           BEGIN 2LPT PART                                 *
 * ************************************************************************* */
// reference: reference: Scoccimarro R., 1998, MNRAS, 299, 1097-1118 Appendix D
  if(SECOND_ORDER_LPT_CORRECTIONS){
      
    fprintf(stderr, "Begin initialization 2LPT velocity field\nTotal elapsed time: %ds\n", time(NULL) - start_time);
    last_time = time(NULL);
    // allocate memory for the velocity boxes and read them in
    vx_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vx_2LPT){
      fprintf(stderr, "update_halo_pos: Error allocating memory for 2LPT velocity box\nAborting...\n");
      free(vx);  free(vy); free(vz);
      return -1;
     }
    vy_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vy_2LPT){
      fprintf(stderr, "update_halo_pos: Error allocating memory for 2LPT velocity box\nAborting...\n");
      free(vx_2LPT);
      free(vx);  free(vy); free(vz);
      return -1;
    }
    vz_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    if (!vz_2LPT){
      fprintf(stderr, "update_halo_pos: Error allocating memory for 2LPT velocity box\nAborting...\n");
      free(vx_2LPT); free(vy_2LPT);
      free(vx);  free(vy); free(vz);
      return -1;
    }

    // read again velocities

    sprintf(filename, "../Boxes/vxoverddot_2LPT_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
    if (mod_fread(vx_2LPT, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "update_halo_pos: Read error occured while reading velocity 2LPT box.\n");
      free(vx);  free(vy); free(vz);
      free(vx_2LPT);  free(vy_2LPT); free(vz_2LPT);
      return -1;
    }
    fclose(F);
    fprintf(stderr, "Read 2LPT vx velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
    last_time = time(NULL);

    sprintf(filename, "../Boxes/vyoverddot_2LPT_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
    if (mod_fread(vy_2LPT, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "update_halo_pos: Read error occured while reading velocity 2LPT box.\n");
      free(vx);  free(vy); free(vz);
      free(vx_2LPT);  free(vy_2LPT); free(vz_2LPT);
      return -1;
    }
    fclose(F);
    fprintf(stderr, "Read 2LPT vy velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
    last_time = time(NULL);

    sprintf(filename, "../Boxes/vzoverddot_2LPT_%i_%.0fMpc", HII_DIM, BOX_LEN);
    F=fopen(filename, "rb");
   if (mod_fread(vz_2LPT, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "update_halo_pos: Read error occured while reading velocity box.\n");
      free(vx);  free(vy); free(vz);
      free(vx_2LPT);  free(vy_2LPT); free(vz_2LPT);
      return -1;
    }
    fclose(F);
    fprintf(stderr, "Read 2LPT vz velocity field\nElapsed time: %ds\n", time(NULL) - last_time);
    last_time = time(NULL);

    // now add the missing factor in eq. D9
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      vx_2LPT[ct] *= displacement_factor_2LPT / BOX_LEN; // this is now comoving displacement in units of box size
      vy_2LPT[ct] *= displacement_factor_2LPT / BOX_LEN; // this is now comoving displacement in units of box size
      vz_2LPT[ct] *= displacement_factor_2LPT / BOX_LEN; // this is now comoving displacement in units of box size
    }
    fprintf(stderr, "Read 2LPT velocity field\nTotal time: %ds\n", time(NULL) - start_time);
  }

/* ************************************************************************* *
 *                            END 2LPT PART                                  *
 * ************************************************************************* */



  /******************   END INITIALIZATION     ********************************/

    float mean_correction = 0.0, mean_correction_2LPT = 0.0, mean_ratio = 0.0;
    float max_correction = 1e-10, max_correction_2LPT = 1e-10, max_ratio = 1e-10;
    int den = 0;

  last_time = time(NULL);
  // read in the halo list
  sprintf(filename, "../Output_files/Halo_lists/halos_z%.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
  F = fopen(filename, "r");
  if (!F){
    fprintf(stderr, "update_halo_pos: Error opening input file: %s\nAborting\n", filename);
    free(vx);  free(vy); free(vz);
    return -1;
  }
  // now read in all halos above our threshold into the smoothed halo field
  fprintf(stderr, "Updating halo positions\n");
  fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
  while (!feof(F)){
    i = xf*HII_DIM;
    j = yf*HII_DIM;
    k = zf*HII_DIM;

    // get new positions using linear velocity displacement from z=INITIAL
    xf += vx[HII_R_INDEX(i,j,k)];
    yf += vy[HII_R_INDEX(i,j,k)];
    zf += vz[HII_R_INDEX(i,j,k)];

    // 2LPT PART
    // add second order corrections
    if(SECOND_ORDER_LPT_CORRECTIONS){
      xf -= vx_2LPT[HII_R_INDEX(i,j,k)];
      yf -= vy_2LPT[HII_R_INDEX(i,j,k)];
      zf -= vz_2LPT[HII_R_INDEX(i,j,k)];

// DEBUG
      //fprintf(stderr, "Displacements ratio: %.2e\t%.2e\n%.2e\t%.2e\n%.2e\t%.2e\n",  vx[HII_R_INDEX(i,j,k)], vx_2LPT[HII_R_INDEX(i,j,k)], vy[HII_R_INDEX(i,j,k)], vy_2LPT[HII_R_INDEX(i,j,k)], vz[HII_R_INDEX(i,j,k)], vz_2LPT[HII_R_INDEX(i,j,k)]);
      if(fabs(mass - 1.72e10) < 0.1e10){
        den += 3;
        mean_correction += fabs(vx[HII_R_INDEX(i, j, k)]) + fabs(vy[HII_R_INDEX(i, j, k)]) + fabs(vz[HII_R_INDEX(i, j, k)]);
        mean_correction_2LPT += fabs(vx_2LPT[HII_R_INDEX(i, j, k)]) + fabs(vy_2LPT[HII_R_INDEX(i, j, k)]) + fabs(vz_2LPT[HII_R_INDEX(i, j, k)]);
        mean_ratio +=  fabs(vx_2LPT[HII_R_INDEX(i, j, k)]/vx[HII_R_INDEX(i, j, k)]) + fabs(vy_2LPT[HII_R_INDEX(i, j, k)]/vy[HII_R_INDEX(i, j, k)]) + fabs(vz_2LPT[HII_R_INDEX(i, j, k)]/vz[HII_R_INDEX(i, j, k)]); 

        max_correction = max(max_correction,  fabs(vx[HII_R_INDEX(i, j, k)]));
        max_correction = max(max_correction,  fabs(vy[HII_R_INDEX(i, j, k)]));
        max_correction = max(max_correction,  fabs(vz[HII_R_INDEX(i, j, k)]));

        max_correction_2LPT = max(max_correction_2LPT,  fabs(vx_2LPT[HII_R_INDEX(i, j, k)]));
        max_correction_2LPT = max(max_correction_2LPT,  fabs(vy_2LPT[HII_R_INDEX(i, j, k)]));
        max_correction_2LPT = max(max_correction_2LPT,  fabs(vz_2LPT[HII_R_INDEX(i, j, k)]));


        max_ratio =  max(max_ratio, fabs(vx_2LPT[HII_R_INDEX(i, j, k)]/vx[HII_R_INDEX(i, j, k)])); 
        max_ratio =  max(max_ratio, fabs(vy_2LPT[HII_R_INDEX(i, j, k)]/vy[HII_R_INDEX(i, j, k)])); 
        max_ratio =  max(max_ratio, fabs(vz_2LPT[HII_R_INDEX(i, j, k)]/vz[HII_R_INDEX(i, j, k)])); 
      }
// END


    }
    // check if we wrapped around, not the casting to ensure < 1.00000
    DI = 10000;
    xf = roundf(xf*DI);
    yf = roundf(yf*DI);
    zf = roundf(zf*DI);
    while (xf >= (float)DI){ xf -= DI;}
    while (xf < 0){ xf += DI;}
    while (yf >= (float)DI){ yf -= DI;}
    while (yf < 0){ yf += DI;}
    while (zf >= (float)DI){ zf -= DI;}
    while (zf < 0){ zf += DI;}
    xf = fabs(xf/(float)DI); // fabs gets rid of minus sign in -0.00000
    yf = fabs(yf/(float)DI);
    zf = fabs(zf/(float)DI);
    /*
    DI = 10000;
    xi = xf*DI;
    yi = yf*DI;
    zi = zf*DI;
    while (xi >= DI){ xi -= DI;}
    while (xi < 0){ xi += DI;}
    while (yi >= DI){ yi -= DI;}
    while (yi < 0){ yi += DI;}
    while (zi >= DI){ zi -= DI;}
    while (zi < 0){ zi += DI;}
    xf = ((float)xi) / ((float)DI);
    yf = ((float)yi) / ((float)DI);
    zf = ((float)zi) / ((float)DI);
    */
    // now write out the updated positions
    fprintf(OUT, "%e\t%f\t%f\t%f\n", mass, xf, yf, zf);

    fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
  }
  fprintf(stderr, "Done in %ds\nTotal elapsed time: %ds\n", time(NULL) - last_time, time(NULL) - start_time);

  mean_correction /= (float)den;
  mean_correction_2LPT /= (float)den;
  mean_ratio /= (float)den;
  fprintf(stderr, "mc = %.2e\tmc2lpt = %.2e\tmr = %.2e\n maxc = %.2e\tmaxc2lpt = %.2e\tmaxr = %.2e\n", mean_correction, mean_correction_2LPT, mean_ratio, max_correction, max_correction_2LPT, max_ratio);

  // deallocate
  free(vx_2LPT);  free(vy_2LPT); free(vz_2LPT); 
  free(vx);  free(vy); free(vz);  fclose(F); fclose(OUT);


  return 0;
}
