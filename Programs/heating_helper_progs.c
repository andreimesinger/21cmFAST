#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "bubble_helper_progs.c"
#include "elec_interp.c"

#define NSPEC_MAX (int) 23
#define RECFAST_NPTS (int) 501
#define KAPPA_10_NPTS (int) 27
#define KAPPA_10_elec_NPTS (int) 20
#define KAPPA_10_pH_NPTS (int) 17

/* Define some global variables; yeah i know it isn't "good practice" but doesn't matter */
double zpp_edge[NUM_FILTER_STEPS_FOR_Ts], sigma_atR[NUM_FILTER_STEPS_FOR_Ts], sigma_Tmin[NUM_FILTER_STEPS_FOR_Ts], ST_over_PS[NUM_FILTER_STEPS_FOR_Ts], sum_lyn[NUM_FILTER_STEPS_FOR_Ts], R_values[NUM_FILTER_STEPS_FOR_Ts];
unsigned long long box_ct;
double const_zp_prefactor, dt_dzp, x_e_ave;
double growth_factor_zp, dgrowth_factor_dzp, PS_ION_EFF;
int NO_LIGHT;
float M_MIN_at_z, M_MIN_at_zp;
int HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY; // New in v2
float F_STAR10,F_ESC10,ALPHA_STAR,ALPHA_ESC,M_TURN,T_AST,Mlim_Fstar,Mlim_Fesc,M_MIN,Splined_Fcoll; // New in v2
double X_LUMINOSITY;
float growth_zpp; // New in v2
static float determine_zpp_max, determine_zpp_min,zpp_bin_width; // new in v2
float *second_derivs_Nion_zpp[NUM_FILTER_STEPS_FOR_Ts]; // New
float *redshift_interp_table;
int Nsteps_zp; //New in v2 
float *zpp_interp_table; //New in v2
gsl_interp_accel *SFRDLow_zpp_spline_acc[NUM_FILTER_STEPS_FOR_Ts];
gsl_spline *SFRDLow_zpp_spline[NUM_FILTER_STEPS_FOR_Ts];

FILE *LOG;

/* initialization routine */
int init_heat();

/* destruction/deallocation routine */
void destruct_heat(); 

 /* returns the spectral emissity */
double spectral_emissivity(double nu_norm, int flag);

/* Ionization fraction from RECFAST. */
double xion_RECFAST(float z, int flag);

/* IGM temperature from RECFAST; includes Compton heating and adiabatic expansion only. */
double T_RECFAST(float z, int flag);

/* Main driver for evolution */
void evolveInt(float zp, float curr_delNL0[], double freq_int_heat[], 
	       double freq_int_ion[], double freq_int_lya[], 
	       int COMPUTE_Ts, double y[], double deriv[]);
		   //float Mturn, float ALPHA_STAR, float F_STAR10, float T_AST);

float dfcoll_dz(float z, float Tmin, float del_bias, float sig_bias);

/* Compton heating rate */
double dT_comp(double z, double TK, double xe);

/* Calculates the optical depth for a photon arriving at z = zp with frequency nu,
 emitted at z = zpp */
double tauX(double nu, double x_e, double zp, double zpp, double HI_filling_factor_zp); 

/* The total weighted HI + HeI + HeII  cross-section in pcm^-2 */
double species_weighted_x_ray_cross_section(double nu, double x_e); 

/* Returns the frequency threshold where \tau = 1 between zp and zpp,
 in the IGM with mean electron fraction x_e */
double nu_tau_one(double zp, double zpp, double x_e, double HI_filling_factor_zp); 

 /* Main integral driver for the frequency integral in the evolution equations */
double integrate_over_nu(double zp, double local_x_e, double lower_int_limit, int FLAG);

/* Returns the maximum redshift at which a Lyn transition contributes to Lya 
   flux at z */
float zmax(float z, int n);

/* Returns frequency of Lyman-n, in units of Lyman-alpha */
double nu_n(int n);

/* Returns recycling fraction (=fraction of photons converted into 
 * Lyalpha for Ly-n resonance */
double frecycle(int n);


 /* returns the spin temperature */
float get_Ts(float z, float delta, float TK, float xe, float Jalpha, float * curr_xalpha);

/* Spin Temperature helper functions */
double xcoll(double z, double TK, double delta, double xe);
double xcoll_HI(double z, double TK, double delta, double xe);
double xcoll_elec(double z, double TK, double delta, double xe);
double xcoll_prot(double z, double TK, double delta, double xe);
double kappa_10_pH(double T, int flag);
double kappa_10_elec(double T, int flag);
double kappa_10(double TK, int flag);
double xalpha_tilde(double z, double Jalpha, double TK, double TS,
		    double delta, double xe);
double taugp(double z, double delta, double xe);
double Salpha_tilde(double TK, double TS, double tauGP);
double Tc_eff(double TK, double TS);
/**********  END PROTOTYPE DEFINITIONS  *******************************/


int init_heat()
{
  kappa_10(1.0,1);
  if( kappa_10_elec(1.0,1) < 0)
    return -2;
  if( kappa_10_pH(1.0,1) < 0)
    return -3;
  if (T_RECFAST(100, 1) < 0)
    return -4;
  if (xion_RECFAST(100, 1) < 0)
    return -5;
  if (spectral_emissivity(0,1) < 0)
    return -6;

  initialize_interp_arrays();

  return 0;
}


void destruct_heat()
{
  kappa_10(1.0,2);
  kappa_10_elec(1.0,2);
  kappa_10_pH(1.0,2);
  T_RECFAST(100.0,2);
  xion_RECFAST(100.0,2);
  spectral_emissivity(0.0, 2);
}


/* Returns the maximum redshift at which a Lyn transition contributes to Lya 
   flux at z */
float zmax(float z, int n){
  double num, denom;
  num = 1 - pow(n+1, -2);
  denom = 1 - pow(n, -2);
  return (1+z)*num/denom - 1;
}


/* Returns frequency of Lyman-n, in units of Lyman-alpha */
double nu_n(int n)
{
  double ans;

  ans = 1.0 - pow(n, -2.0);
  ans /= 0.75;
  return ans;
}


/* Returns recycling fraction (=fraction of photons converted into 
 * Lyalpha for Ly-n resonance */
double frecycle(int n)
{
  switch (n){
  case 0:
    return 1;
  case 1:
    return 1;
  case 2:
    return 1;
  case 3:
    return 0;
  case 4:
    return 0.2609;
  case 5:
    return 0.3078;
  case 6:
    return 0.3259;
  case 7:
    return 0.3353;
  case 8:
    return 0.3410;
  case 9:
    return 0.3448;
  case 10:
    return 0.3476;
  case 11:
    return 0.3496;
  case 12:
    return 0.3512;
  case 13:
    return 0.3524;
  case 14:
    return 0.3535;
  case 15:
    return 0.3543;
  case 16:
    return 0.3550;
  case 17:
    return 0.3556;
  case 18:
    return 0.3561;
  case 19:
    return 0.3565;
  case 20:
    return 0.3569;
  case 21:
    return 0.3572;
  case 22:
    return 0.3575;
  case 23:
    return 0.3578;
  case 24:
    return 0.3580;
  case 25:
    return 0.3582;
  case 26:
    return 0.3584;
  case 27:
    return 0.3586;
  case 28:
    return 0.3587;
  case 29:
    return 0.3589;
  case 30:
    return 0.3590;
  default:
    return 0;
  }
}


/* Reads in and constructs table of the piecewise power-law fits to Pop 2 and 
 * Pop 3 stellar spectra, from Barkana */
double spectral_emissivity(double nu_norm, int flag)
{
  static int n[NSPEC_MAX];
  static float nu_n[NSPEC_MAX], alpha_S_2[NSPEC_MAX];
  static float alpha_S_3[NSPEC_MAX], N0_2[NSPEC_MAX], N0_3[NSPEC_MAX];
  double n0_fac;
  double ans, tot, lya;
  int i;
  FILE *F;

  if (flag == 1) {
    /* Read in the data */
    if (!(F = fopen(STELLAR_SPECTRA_FILENAME, "r"))){
      fprintf(stderr, "spectral_emissivity: Unable to open file: stellar_spectra.dat for reading\nAborting\n");
      fprintf(LOG, "spectral_emissivity: Unable to open file: stellar_spectra.dat for reading\nAborting\n");
      return -1;
    }

    for (i=1;i<NSPEC_MAX;i++) {
      fscanf(F, "%i %e %e %e %e", &n[i], &N0_2[i], &alpha_S_2[i], &N0_3[i], &alpha_S_3[i]);
      //      printf("%i\t%e\t%e\t%e\t%e\n", n[i], N0_2[i], alpha_S_2[i], N0_3[i], alpha_S_3[i]);
    }
    fclose(F);

    for (i=1;i<NSPEC_MAX;i++) {
      nu_n[i] = 4.0/3.0*(1.0-1.0/pow(n[i],2.0));
    }

    for (i=1;i<(NSPEC_MAX-1);i++) {
      n0_fac = (pow(nu_n[i+1],alpha_S_2[i]+1) - pow(nu_n[i],alpha_S_2[i]+1));
      N0_2[i] *= (alpha_S_2[i]+1)/n0_fac*Pop2_ion;
      n0_fac = (pow(nu_n[i+1],alpha_S_3[i]+1) - pow(nu_n[i],alpha_S_3[i]+1));
      N0_3[i] *= (alpha_S_3[i]+1)/n0_fac*Pop3_ion;
    }

    return 0.0;
  }

  ans = 0.0;
  for (i=1;i<(NSPEC_MAX-1);i++) {
    //    printf("checking between %e and %e\n", nu_n[i], nu_n[i+1]);
    if ((nu_norm >= nu_n[i]) && (nu_norm < nu_n[i+1])) {
      // We are in the correct spectral region
      if (Pop == 2)
	ans = N0_2[i]*pow(nu_norm,alpha_S_2[i]);
      else
	ans = N0_3[i]*pow(nu_norm,alpha_S_3[i]);
      //           printf("nu_norm=%e i=%i ans=%e\n", nu_norm, i, ans);
      return ans/Ly_alpha_HZ;
    }
  }

  i= NSPEC_MAX-1;
  if (Pop == 2)
    return  N0_2[i]*pow(nu_norm,alpha_S_2[i])/Ly_alpha_HZ;
  else
    return N0_3[i]*pow(nu_norm,alpha_S_3[i])/Ly_alpha_HZ;
  //  return 0;
}

/********************************************************************
 ************************** IGM Evolution ***************************
  This function creates the d/dz' integrands
*********************************************************************/
void evolveInt(float zp, float curr_delNL0[], double freq_int_heat[], 
	       double freq_int_ion[], double freq_int_lya[], 
	       int COMPUTE_Ts, double y[], double deriv[]){
  double  dfdzp, dadia_dzp, dcomp_dzp, dxheat_dt, ddz, dxion_source_dt, dxion_sink_dt;
  double zpp, dzpp, nu_temp;
  int zpp_ct,ithread;
  double T, x_e, dTdzp, dx_edzp, dfcoll, zpp_integrand;
  double dxe_dzp, n_b, dspec_dzp, dxheat_dzp, dxlya_dt, dstarlya_dt;
  // New in v2
  float growth_zpp,fcoll;

  x_e = y[0];
  T = y[1];
  n_b = N_b0 * pow(1+zp, 3) * (1+curr_delNL0[0]*growth_factor_zp);


  // First, let's do the trapazoidal integration over zpp
  ithread = omp_get_thread_num();

  dxheat_dt = 0;
  dxion_source_dt = 0;
  dxlya_dt = 0;
  dstarlya_dt = 0;
  if (!NO_LIGHT){
  for (zpp_ct = 0; zpp_ct < NUM_FILTER_STEPS_FOR_Ts; zpp_ct++){
    // set redshift of half annulus; dz'' is negative since we flipped limits of integral
    if (zpp_ct==0){
      zpp = (zpp_edge[0]+zp)*0.5;
      dzpp = zp - zpp_edge[0];
    }
    else{
      zpp = (zpp_edge[zpp_ct]+zpp_edge[zpp_ct-1])*0.5;
      dzpp = zpp_edge[zpp_ct-1] - zpp_edge[zpp_ct];
    }
	//New in v2
    if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
	  growth_zpp = dicke(zpp);
	  // Interpolate Fcoll -------------------------------------------------------------------------------------
	  if (curr_delNL0[zpp_ct]*growth_zpp < 1.5){
        if (curr_delNL0[zpp_ct]*growth_zpp < -1.) {
		  fcoll = 0;
        }
        else {
          fcoll = gsl_spline_eval(SFRDLow_zpp_spline[zpp_ct], log10(curr_delNL0[zpp_ct]*growth_zpp+1.), SFRDLow_zpp_spline_acc[zpp_ct]);
          fcoll= pow(10., fcoll);
        }
      }
      else {
        if (curr_delNL0[zpp_ct]*growth_zpp < 0.99*Deltac) {
          // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
          // However, such densities should always be collapsed, so just set f_coll to unity. 
          // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
          splint(Overdense_high_table-1,SFRD_z_high_table[zpp_ct]-1,second_derivs_Nion_zpp[zpp_ct]-1,NSFR_high,curr_delNL0[zpp_ct]*growth_zpp,&(fcoll));
        }
        else {
          fcoll = 1.;
        }
      }
      if (fcoll > 1.) fcoll = 1.;
	  // Find Fcoll end ----------------------------------------------------------------------------------

	  /* Instead of dfcoll/dz we compute fcoll/(T_AST*H(z)^-1)*(dt/dz), 
	    where T_AST is the typical star-formation timescale, in units of the Hubble time.
		This is the same parameter with 't_STAR' (defined in ANAL_PARAMS.H).
		If turn the new parametrization on, this is a free parameter.
		*/
	  dfcoll = ST_over_PS[zpp_ct]*(double)fcoll*hubble(zpp)/T_AST*fabs(dtdz(zpp))*fabs(dzpp);
	}
	else {
      dfcoll = dfcoll_dz(zpp, sigma_Tmin[zpp_ct], curr_delNL0[zpp_ct], sigma_atR[zpp_ct]);
      dfcoll *= ST_over_PS[zpp_ct] * dzpp; // this is now a positive quantity
	}
    zpp_integrand = dfcoll * (1+curr_delNL0[zpp_ct]*dicke(zpp)) * pow(1+zpp, -X_RAY_SPEC_INDEX);

    dxheat_dt += zpp_integrand * freq_int_heat[zpp_ct];
    dxion_source_dt += zpp_integrand * freq_int_ion[zpp_ct];
    if (COMPUTE_Ts){
      dxlya_dt += zpp_integrand * freq_int_lya[zpp_ct];
      dstarlya_dt += dfcoll * (1+curr_delNL0[zpp_ct]*dicke(zpp)) * pow(1+zp,2)*(1+zpp)
                     * sum_lyn[zpp_ct];
    }
  }

  // add prefactors
  dxheat_dt *= const_zp_prefactor;
  dxion_source_dt *= const_zp_prefactor;
  if (COMPUTE_Ts){
    dxlya_dt *= const_zp_prefactor*n_b;
    dstarlya_dt *= F_STAR10 * C * N_b0 / FOURPI;

    /*
    if ((dxlya_dt < 0) || (dstarlya_dt<0)){
         printf("***Jalpha_x=%e, Jalpha_star=%e\n", dxlya_dt, dstarlya_dt);
  for (zpp_ct = 0; zpp_ct < NUM_FILTER_STEPS_FOR_Ts; zpp_ct++)
    printf("freq_int=%e, sum_lyn=%e\n", freq_int_lya[zpp_ct], sum_lyn[zpp_ct]);
    }
    */
  }
  } // end NO_LIGHT if statement

  /**** Now we can solve the evolution equations  *****/
  /*** First let's do dxe_dzp ***/
  dxion_sink_dt = alpha_A(T) * CLUMPING_FACTOR * x_e*x_e * f_H * n_b;
  dxe_dzp = dt_dzp*(dxion_source_dt - dxion_sink_dt);
  //printf("In evolveInt: z=%.4f, dt_dzp=%.4e, dxion_source_dt=%.4e, dxion_sink_dt==%.4e\n",zpp,dt_dzp,dxion_source_dt,dxion_sink_dt);
  deriv[0] = dxe_dzp;
  /*
  if (box_ct==1000000){
    printf("%g\t%g\t%g\n", x_e, zpp_integrand);
    printf("(z', freq_int_ion[zpp_ct])\n");
    for (zpp_ct = 0; zpp_ct < NUM_FILTER_STEPS_FOR_Ts; zpp_ct++){
      printf("(%g, %g)\n", zpp_ct, freq_int_ion[zpp_ct]);
    }
    printf("source term is %e, sink term is %e\n", dt_dzp*dxion_source_dt, dt_dzp*dxion_sink_dt);
    printf("dxe_dzp=%e\n", dxe_dzp);
  }
  */

  /*** Next, let's get the temperature components ***/
  // first, adiabatic term
  dadia_dzp = 3/(1.0+zp);
  if (fabs(curr_delNL0[0]) > FRACT_FLOAT_ERR) // add adiabatic heating/cooling from structure formation
    dadia_dzp += dgrowth_factor_dzp/(1.0/curr_delNL0[0]+growth_factor_zp);
  dadia_dzp *= (2.0/3.0)*T;
  //  printf("dTstructure/dz=%e, total Adiabatic heat/dzp (structure + expansion) at zp=%.3f = %e\n", dgrowth_factor_dzp/(1.0/curr_delNL0[0]+growth_factor_zp), zp, dadia_dzp);

  // next heating due to the changing species
  dspec_dzp = - dxe_dzp * T / (1+x_e);
  //  printf("changing species, dT/dz'=%e\n", dspec_dzp);

  // next, Compton heating
  dcomp_dzp = dT_comp(zp, T, x_e);
  //  printf("Compton heat/dzp at zp=%.3f = %e\n", zp, dcomp_dzp);

  // lastly, X-ray heating
  dxheat_dzp = dxheat_dt * dt_dzp * 2.0 / 3.0 / k_B / (1.0+x_e);
  //  printf("x-ray heating, dT/dz=%e\n", dxheat_dzp);

  // summing them up...
  deriv[1] = dxheat_dzp + dcomp_dzp + dspec_dzp + dadia_dzp;
  //  printf("total dT/dz=%e\n", deriv[1]);

  /*** Finally, if we are at the last redshift step, Lya ***/
  deriv[2] = dxlya_dt + dstarlya_dt;

  // stuff for marcos
  deriv[3] = dxheat_dzp;
  deriv[4] = dt_dzp*dxion_source_dt;
}

/*
  Evaluates the frequency integral in the Tx evolution equation
  photons starting from zpp arive at zp, with mean IGM electron
  fraction of x_e (used to compute tau), and local electron 
  fraction local_x_e
  FLAG = 0 for heat integral
  FLAG = 1 for ionization integral
  FLAG = 2 for Lya integral
*/
double integrand_in_nu_heat_integral(double nu, void * params){
  double species_sum, fheat;
  float x_e = *(double *) params;

  // HI
  species_sum = interp_fheat((nu - NUIONIZATION)/NU_over_EV, x_e)
               * hplank*(nu - NUIONIZATION) * f_H * (1-x_e) * HI_ion_crosssec(nu);

  // HeI
  species_sum += interp_fheat((nu - HeI_NUIONIZATION)/NU_over_EV, x_e)
               * hplank*(nu - HeI_NUIONIZATION) * f_He * (1-x_e) * HeI_ion_crosssec(nu);

  // HeII
  species_sum += interp_fheat((nu - HeII_NUIONIZATION)/NU_over_EV, x_e)
               * hplank*(nu - HeII_NUIONIZATION) * f_He * x_e * HeII_ion_crosssec(nu);

  return species_sum * pow(nu/NU_X_THRESH, -X_RAY_SPEC_INDEX-1);
}
double integrand_in_nu_ion_integral(double nu, void * params){
  double species_sum, fheat, F_i;
  float x_e = *(double *) params;

  // photoionization of HI, prodicing e- of energy h*(nu - nu_HI)
  F_i = interp_nion_HI((nu - NUIONIZATION)/NU_over_EV, x_e) +
    interp_nion_HeI((nu - NUIONIZATION)/NU_over_EV, x_e) +
    interp_nion_HeII((nu - NUIONIZATION)/NU_over_EV, x_e) + 1;
  species_sum = F_i * f_H * (1-x_e) * HI_ion_crosssec(nu);

  // photoionization of HeI, prodicing e- of energy h*(nu - nu_HeI)
  F_i = interp_nion_HI((nu - HeI_NUIONIZATION)/NU_over_EV, x_e) +
    interp_nion_HeI((nu - HeI_NUIONIZATION)/NU_over_EV, x_e) +
    interp_nion_HeII((nu - HeI_NUIONIZATION)/NU_over_EV, x_e) + 1;
  species_sum += F_i * f_He * (1-x_e) * HeI_ion_crosssec(nu);

  // photoionization of HeII, prodicing e- of energy h*(nu - nu_HeII)
  F_i = interp_nion_HI((nu - HeII_NUIONIZATION)/NU_over_EV, x_e) +
    interp_nion_HeI((nu - HeII_NUIONIZATION)/NU_over_EV, x_e) +
    interp_nion_HeII((nu - HeII_NUIONIZATION)/NU_over_EV, x_e) + 1;
  species_sum += F_i * f_He * x_e * HeII_ion_crosssec(nu);

  return species_sum * pow(nu/NU_X_THRESH, -X_RAY_SPEC_INDEX-1);
}
double integrand_in_nu_lya_integral(double nu, void * params){
  double species_sum, fheat;
  float x_e = *(double *) params;

  // HI
  species_sum = interp_n_Lya((nu - NUIONIZATION)/NU_over_EV, x_e)
    * f_H * (double)(1-x_e) * HI_ion_crosssec(nu);

  // HeI
  species_sum += interp_n_Lya((nu - HeI_NUIONIZATION)/NU_over_EV, x_e)
    * f_He * (double)(1-x_e) * HeI_ion_crosssec(nu);

  // HeII
  species_sum += interp_n_Lya((nu - HeII_NUIONIZATION)/NU_over_EV, x_e)
    * f_He * (double)x_e * HeII_ion_crosssec(nu);

  return species_sum * pow(nu/NU_X_THRESH, -X_RAY_SPEC_INDEX-1);
}
double integrate_over_nu(double zp, double local_x_e, double lower_int_limit, int FLAG){
       double result, error;
       double rel_tol  = 0.01; //<- relative tolerance
       gsl_function F;
       gsl_integration_workspace * w 
	 = gsl_integration_workspace_alloc (1000);

       if (DEBUG_ON){
	 printf("integrate over nu, parameters: %f, %f, %e, %i, thread# %i\n", zp, local_x_e, lower_int_limit, FLAG, omp_get_thread_num());
       }
       /*
       if (DO_NOT_COMPARE_NUS)
	 lower_int_limit = NU_X_THRESH;
       else
	 lower_int_limit = FMAX(nu_tau_one(zp, zpp, x_e, HI_filling_factor_zp), NU_X_THRESH);
       */

       F.params = &local_x_e;

       //       printf("lower limit is %e\n", lower_int_limit / NU_over_EV);
       if (FLAG==0)
	 F.function = &integrand_in_nu_heat_integral;
       else if (FLAG==1)
	 F.function = &integrand_in_nu_ion_integral;
       else {
	 F.function = &integrand_in_nu_lya_integral;
       }

       gsl_integration_qag (&F, lower_int_limit, 100*lower_int_limit, 
			    0, rel_tol, 1000, GSL_INTEG_GAUSS61, w, &result, &error); 
       gsl_integration_workspace_free (w);


       // if it is the Lya integral, add prefactor
       if (FLAG == 2)
	 return result * C / FOURPI / Ly_alpha_HZ / hubble(zp);

       /*
       if (isnan(result))
	 fprintf(stderr, "We have a NaN in the intergrator with calling params: %g,%g,%g,%i\n", zp, local_x_e, lower_int_limit, FLAG);
       */
       return result;
}



/*
  The total weighted HI + HeI + HeII  cross-section in pcm^-2
  technically, the x_e should be local, line of sight (not global) here,
  but that would be very slow...
*/
double species_weighted_x_ray_cross_section(double nu, double x_e){
  double HI_factor, HeI_factor, HeII_factor;

  HI_factor = f_H * (1-x_e) * HI_ion_crosssec(nu);
  HeI_factor = f_He * (1-x_e) * HeI_ion_crosssec(nu);
  HeII_factor = f_He * x_e * HeII_ion_crosssec(nu);

  return HI_factor + HeI_factor + HeII_factor;
}



/*
  Calculates the optical depth for a photon arriving at z = zp with frequency nu,
  emitted at z = zpp.
  The filling factor of neutral IGM at zp is HI_filling_factor_zp.
*/
typedef struct{
  double nu_0, x_e, ion_eff;
} tauX_params;
double tauX_integrand(double zhat, void *params){
  double n, drpropdz, nuhat, HI_filling_factor_zhat, sigma_tilde, fcoll;
  tauX_params *p = (tauX_params *) params;
  float Splined_ans; 

  drpropdz = C * dtdz(zhat);
  n = N_b0 * pow(1+zhat, 3);
  nuhat = p->nu_0 * (1+zhat);
  // New in v2
  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
	Nion_ST_z(zhat,&(Splined_ans));
	fcoll = Splined_ans;
  }
  else {
    fcoll = FgtrM(zhat, M_MIN);
  }
  if (fcoll < 1e-20)
    HI_filling_factor_zhat = 1;
  else    
    HI_filling_factor_zhat = 1 - p->ion_eff * fcoll/(1.0 - x_e_ave); //simplification to use the <x_e> value at zp and not zhat.  should'nt matter much since the evolution in x_e_ave is slower than fcoll.  in principle should make an array to store past values of x_e_ave..
  if (HI_filling_factor_zhat < 1e-4) HI_filling_factor_zhat = 1e-4; //set a floor for post-reionization stability

  sigma_tilde = species_weighted_x_ray_cross_section(nuhat, p->x_e);

  return drpropdz * n * HI_filling_factor_zhat * sigma_tilde;
}
double tauX(double nu, double x_e, double zp, double zpp, double HI_filling_factor_zp){
//double tauX(double nu, double x_e, double zp, double zpp, double HI_filling_factor_zp, float M_TURN, float ALPHA_STAR, float F_STAR10){
  double result, error, fcoll;
       gsl_function F;
       double rel_tol  = 0.005; //<- relative tolerance
       gsl_integration_workspace * w 
	 = gsl_integration_workspace_alloc (1000);
       tauX_params p;
  float Splined_ans; 

       /*
       if (DEBUG_ON)
	 printf("in taux, parameters are: %e, %e, %f, %f, %e\n", nu, x_e, zp, zpp, HI_filling_factor_zp);
       */
       F.function = &tauX_integrand;
       p.nu_0 = nu/(1+zp);
       p.x_e = x_e;

       // New in v2
       if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
         p.ion_eff = N_GAMMA_UV*F_STAR10*F_ESC10;
       }    
       else {
         if (HI_filling_factor_zp > FRACT_FLOAT_ERR){
           // effective efficiency for the PS (not ST) mass function; quicker to compute...
           fcoll = FgtrM(zp, M_MIN);
           p.ion_eff = (1.0 - HI_filling_factor_zp) / fcoll * (1.0 - x_e_ave);
           PS_ION_EFF = p.ion_eff;
         }    
         else 
           p.ion_eff = PS_ION_EFF; // uses the previous one in post reionization regime
       }  

       F.params = &p;
       gsl_integration_qag (&F, zpp, zp, 0, rel_tol,
			    1000, GSL_INTEG_GAUSS61, w, &result, &error); 
       gsl_integration_workspace_free (w);

       /*
       if (DEBUG_ON)
	   printf("returning from tauX, return value=%e\n", result);
       */
       return result;
}



/*
  Returns the frequency threshold where \tau_X = 1, given parameter values of
  electron fraction in the IGM outside of HII regions, x_e, 
  recieved redshift, zp, and emitted redshift, zpp.
*/
typedef struct{
  double x_e, zp, zpp, HI_filling_factor_zp;
} nu_tau_one_params;
double nu_tau_one_helper(double nu, void * params){
  nu_tau_one_params *p = (nu_tau_one_params *) params;
  return tauX(nu, p->x_e, p->zp, p->zpp, p->HI_filling_factor_zp) - 1;
}
double nu_tau_one(double zp, double zpp, double x_e, double HI_filling_factor_zp){
  int status, iter, max_iter;
  const gsl_root_fsolver_type * T; 
  gsl_root_fsolver * s;
  gsl_function F;
  double x_lo, x_hi, r=0;
  double relative_error = 0.02;
  nu_tau_one_params p;

  if (DEBUG_ON){
    printf("in nu tau one, called with parameters: zp=%f, zpp=%f, x_e=%e, HI_filling_at_zp=%e\n",
	   zp, zpp, x_e, HI_filling_factor_zp);
  }

  // check if too ionized 
  if (x_e > 0.9999){
    fprintf(stderr,
     "Ts.c: WARNING: x_e value is too close to 1 for convergence in nu_tau_one\n");
    return -1;
  }

  // select solver and allocate memory
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T); // non-derivative based Brent method
  if (!s){
    fprintf(stderr, "Ts.c: Unable to allocate memory in function nu_tau_one!\n");
    return -1;
  }

  //check if lower bound has null
  if (tauX(HeI_NUIONIZATION, x_e, zp, zpp, HI_filling_factor_zp) < 1)
    return HeI_NUIONIZATION;

  // set frequency boundary values
  x_lo= HeI_NUIONIZATION;
  x_hi = 1e6 * NU_over_EV;

  // select function we wish to solve
  p.x_e = x_e;
  p.zp = zp;
  p.zpp = zpp;
  p.HI_filling_factor_zp = HI_filling_factor_zp;
  F.function = &nu_tau_one_helper;  
  F.params = &p;
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  // iterate until we guess close enough
  if (DEBUG_ON) printf ("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err(est)");
  iter = 0;
  max_iter = 100;
  do{
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      //      printf("iter%i, r=%e\n", iter, r);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, relative_error);
      if (DEBUG_ON){
	printf ("%5d [%.7e, %.7e] %.7e %.7e\n", iter, x_lo, x_hi, r, (x_hi - x_lo)/r);  
	fflush(NULL);
      }
}
  while (status == GSL_CONTINUE && iter < max_iter);

  // deallocate and return
  gsl_root_fsolver_free (s);
  if (DEBUG_ON) printf("Root found at %e eV", r/NU_over_EV);

  return r;
}


/*
  Redshift derivative of the conditional collapsed fraction
 */
float dfcoll_dz(float z, float sigma_min, float del_bias, float sig_bias)
{
  double dz,z1,z2;
  double mu, m1, m2;
  double fc1,fc2,ans;

  dz = 0.001;
  z1 = z + dz;
  z2 = z - dz;
  fc1 = sigmaparam_FgtrM_bias(z1, sigma_min, del_bias, sig_bias);
  fc2 = sigmaparam_FgtrM_bias(z2, sigma_min, del_bias, sig_bias);
  ans = (fc1 - fc2)/(2.0*dz);
  return ans;
}

/* Compton heating term */
double dT_comp(double z, double TK, double xe)
{
  double Trad,ans;

  Trad = T_cmb*(1.0+z);
  ans = (-1.51e-4) * (xe/(1.0+xe+f_He)) /(hubble(z)/Ho)/hlittle*pow(Trad,4.0)/(1.0+z);
  //fprintf(stderr, "%e\t%e\t%e\t%e\n", xe, pow(Trad,4.0), hubble(z)/Ho, ans);
  ans *= Trad - TK;
  return ans;
}



/********************************************************************
 ************************ RECFAST quantities ************************
 ********************************************************************/

/* IGM temperature from RECFAST; includes Compton heating and adiabatic
 * expansion only. */
double T_RECFAST(float z, int flag)
{
  double ans;
  static double zt[RECFAST_NPTS], TK[RECFAST_NPTS];
  static gsl_interp_accel *acc;
  static gsl_spline *spline;
  float currz, currTK, trash;
  int i;
  FILE *F;

  if (flag == 1) {
    // Read in the recfast data
    if ( !(F=fopen(RECFAST_FILENAME, "r")) ){
      fprintf(stderr, "T_RECFAST: Unable to open file: %s for reading\nAborting\n", RECFAST_FILENAME);
      fprintf(LOG, "T_RECFAST: Unable to open file: %s for reading\nAborting\n", RECFAST_FILENAME);
      return -1;
    }

    for (i=(RECFAST_NPTS-1);i>=0;i--) {
      fscanf(F, "%f %E %E %E", &currz, &trash, &trash, &currTK);
      zt[i] = currz;
      TK[i] = currTK;
    }
    fclose(F);

    // Set up spline table
    acc   = gsl_interp_accel_alloc ();
    spline  = gsl_spline_alloc (gsl_interp_cspline, RECFAST_NPTS);
    gsl_spline_init(spline, zt, TK, RECFAST_NPTS);

    return 0;
  }

  if (flag == 2) {
    // Free memory
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  if (z > zt[RECFAST_NPTS-1]) { // Called at z>500! Bail out
    fprintf(stderr, "Called xion_RECFAST with z=%f, bailing out!\n", z);
    fprintf(LOG, "Called xion_RECFAST with z=%f, bailing out!\n", z);
    return -1;
  }
  else { // Do spline
    ans = gsl_spline_eval (spline, z, acc);
  }
  return ans;
}


/* Ionization fraction from RECFAST. */
double xion_RECFAST(float z, int flag)
{
  static double zt[RECFAST_NPTS], xion[RECFAST_NPTS];
  static gsl_interp_accel *acc;
  static gsl_spline *spline;
  float trash, currz, currxion;
  double ans;
  int i;
  FILE *F;


  if (flag == 1) {
    // Initialize vectors
    if ( !(F=fopen(RECFAST_FILENAME, "r")) ){
      fprintf(stderr, "xion_RECFAST: Unable to open file: %s for reading\nAborting\n", RECFAST_FILENAME);
      fprintf(LOG, "xion_RECFAST: Unable to open file: %s for reading\nAborting\n", RECFAST_FILENAME);
      return -1;
    }

    for (i=(RECFAST_NPTS-1);i>=0;i--) {
      fscanf(F, "%f %E %E %E", &currz, &currxion, &trash, &trash);
      zt[i] = currz;
      xion[i] = currxion;
    }
    fclose(F);

    // Set up spline table
    acc   = gsl_interp_accel_alloc ();
    spline  = gsl_spline_alloc (gsl_interp_cspline, RECFAST_NPTS);
    gsl_spline_init(spline, zt, xion, RECFAST_NPTS);

    return 0;
  }

  if (flag == 2) {
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  if (z > zt[RECFAST_NPTS-1]) { // Called at z>500! Bail out
    fprintf(stderr, "Called xion_RECFAST with z=%f, bailing out!\n", z);
    fprintf(LOG, "Called xion_RECFAST with z=%f, bailing out!\n", z);
    return -1;
  }
  else { // Do spline
    ans = gsl_spline_eval (spline, z, acc);
  }
  return ans;
}


float get_Ts(float z, float delta, float TK, float xe, float Jalpha, float * curr_xalpha){
  double Trad,xc,xa_tilde;
  double TS,TSold,TSinv;
  double Tceff;

  Trad = T_cmb*(1.0+z);
  xc = xcoll(z,TK,delta,xe);
  if (Jalpha > 1.0e-20) { /* Must use WF effect */
    TS = Trad;
    TSold = 0.0;
    while (fabs(TS-TSold)/TS > 1.0e-3) {
      TSold = TS;
      xa_tilde = xalpha_tilde(z,Jalpha,TK,TS,delta,xe);
      Tceff = Tc_eff(TK,TS);
      TSinv = (1.0/Trad+xa_tilde/Tceff + xc/TK)/(1.0+xa_tilde+xc);
      TS = 1.0/TSinv;
    }
    *curr_xalpha = xa_tilde;
  } else { /* Collisions only */
    TSinv = (1.0/Trad + xc/TK)/(1.0 + xc);
    TS = 1.0/TSinv;
    *curr_xalpha = 0;
  }

  if((TS < 0.)&&(SUBCELL_RSD==1)) {
        // Found that in extreme cases it could produce a negative spin temperature
        // If negative, it is a very small number. Take the absolute value, the optical depth can deal with very large numbers, so ok to be small
        TS = fabs(TS);
    }
    

  return TS;
}


double xcoll(double z, double TK, double delta, double xe){
  return xcoll_HI(z,TK,delta,xe) + xcoll_elec(z,TK,delta,xe) + xcoll_prot(z,TK,delta,xe);
}

double xcoll_HI(double z, double TK, double delta, double xe)
{
  double krate,nH,Trad;
  double xcoll;

  Trad = T_cmb*(1.0+z);
  nH = (1.0-xe)*No*pow(1.0+z,3.0)*(1.0+delta);
  krate = kappa_10(TK,0);
  xcoll = T21/Trad*nH*krate/A10_HYPERFINE;
  //  printf("xcoll_HI=%f\n", xcoll);
  return xcoll;
}

/* Note that this assumes Helium ionized same as Hydrogen */
double xcoll_elec(double z, double TK, double delta, double xe)
{
  double krate,ne,Trad;
  double xcoll;

  Trad = T_cmb*(1.0+z);
  ne = xe*N_b0*pow(1.0+z,3.0)*(1.0+delta);
  krate = kappa_10_elec(TK,0);
  xcoll = T21/Trad*ne*krate/A10_HYPERFINE;
  return xcoll;
}

double xcoll_prot(double z, double TK, double delta, double xe)
{
  double krate,np,Trad;
  double xcoll;

  Trad = T_cmb*(1.0+z);
  np = xe*No*pow(1.0+z,3.0)*(1.0+delta);
  krate = kappa_10_pH(TK,0);
  xcoll = T21/Trad*np*krate/A10_HYPERFINE;
  return xcoll;
}


double kappa_10(double TK, int flag)
{
  int i;
  static double tkin[KAPPA_10_NPTS], kap[KAPPA_10_NPTS];
  static gsl_interp_accel *acc;
  static gsl_spline *spline;
  double ans;

  if (flag == 1) { /* Set up spline table */
    /* Initialize kappa from Zygelman (2005), Table 2, column 4 */
    tkin[0] = 1.0; kap[0] = 1.38e-13;
    tkin[1] = 2.0; kap[1] = 1.43e-13;
    tkin[2] = 4.0; kap[2] = 2.71e-13;
    tkin[3] = 6.0; kap[3] = 6.60e-13;
    tkin[4] = 8.0; kap[4] = 1.47e-12;
    tkin[5] = 10.0; kap[5] = 2.88e-12;
    tkin[6] = 15.0; kap[6] = 9.10e-12;
    tkin[7] = 20.0; kap[7] = 1.78e-11;
    tkin[8] = 25.0; kap[8] = 2.73e-11;
    tkin[9] = 30.0; kap[9] = 3.67e-11;
    tkin[10] = 40.0; kap[10] = 5.38e-11;
    tkin[11] = 50.0; kap[11] = 6.86e-11;
    tkin[12] = 60.0; kap[12] = 8.14e-11; 
    tkin[13] = 70.0; kap[13] = 9.25e-11;
    tkin[14] = 80.0; kap[14] = 1.02e-10;
    tkin[15] = 90.0; kap[15] = 1.11e-10;
    tkin[16] = 100.0; kap[16] = 1.19e-10;
    tkin[17] = 200.0; kap[17] = 1.75e-10;
    tkin[18] = 300.0; kap[18] = 2.09e-10;
    tkin[19] = 501.0; kap[19] = 2.565e-10;
    tkin[20] = 701.0; kap[20] = 2.91e-10;
    tkin[21] = 1000.0; kap[21] = 3.31e-10;
    tkin[22] = 2000.0; kap[22] = 4.27e-10;
    tkin[23] = 3000.0; kap[23] = 4.97e-10;
    tkin[24] = 5000.0; kap[24] = 6.03e-10;
    tkin[25] = 7000.0; kap[25] = 6.87e-10;
    tkin[26] = 10000.0; kap[26] = 7.87e-10;

    /* Convert to logs for interpolation */
    for (i=0;i<KAPPA_10_NPTS;i++) {
      tkin[i] = log(tkin[i]);
      kap[i] = log(kap[i]);
    }

    /* Set up spline table */
    acc   = gsl_interp_accel_alloc ();
    spline  = gsl_spline_alloc (gsl_interp_cspline, KAPPA_10_NPTS);
    gsl_spline_init(spline, tkin, kap, KAPPA_10_NPTS);
    return 0;
  } 

  if (flag == 2) { /* Clear memory */
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  if (log(TK) < tkin[0]) { /* Below 1 K, just use that value */
    ans = kap[0];
  } else if (log(TK) > tkin[KAPPA_10_NPTS-1]) { 
    /* Power law extrapolation */
    ans = log(exp(kap[KAPPA_10_NPTS-1])*pow(TK/exp(tkin[KAPPA_10_NPTS-1]),0.381));
  } else { /* Do spline */
    TK = log(TK);
    ans = gsl_spline_eval (spline, TK, acc);
  }
  return exp(ans);
}


/* Interpolate exact results for kappa_10^eH.  The table goes up to
 * 10^5 K, but values are only accurate for T<2x10^4 K.  From
 * Furlanetto & Furlanetto 2006 */
double kappa_10_elec(double T, int flag)
{
  static double TK[KAPPA_10_elec_NPTS], kappa[KAPPA_10_elec_NPTS];
  static gsl_interp_accel *acc;
  static gsl_spline *spline;
  double ans;
  int i;
  float curr_TK, curr_kappa;
  FILE *F;

  if (flag == 1) {
    if (!(F=fopen(KAPPA_EH_FILENAME, "r"))){
      fprintf(stderr, "Unable to open the kappa_10^eH file at %s\nAborting\n", KAPPA_EH_FILENAME);
      fprintf(LOG, "Unable to open the kappa_10^eH file at %s\nAborting\n", KAPPA_EH_FILENAME);      
      return 0;
    }

    for (i=0;i<KAPPA_10_elec_NPTS;i++) {
      fscanf(F, "%f %e", &curr_TK, &curr_kappa);
      TK[i] = curr_TK;
      kappa[i] = curr_kappa;
    }
    fclose(F);

    for (i=0;i<KAPPA_10_elec_NPTS;i++) {
      TK[i] = log(TK[i]);
      kappa[i] = log(kappa[i]);
    }

    /* Set up spline table */
    acc   = gsl_interp_accel_alloc ();
    spline  = gsl_spline_alloc (gsl_interp_cspline, KAPPA_10_elec_NPTS);
    gsl_spline_init(spline, TK, kappa, KAPPA_10_elec_NPTS);
    return 0;
  }

  if (flag == 2) {
    /* Free memory */
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  T = log(T);
  if (T < TK[0]) { /* Use TK=1 K value if called at lower temperature */
    ans = kappa[0];
  }
  else if (T > TK[KAPPA_10_elec_NPTS-1]) { 
    /* Power law extrapolation */
    ans  = kappa[KAPPA_10_elec_NPTS-1] + 
      ((kappa[KAPPA_10_elec_NPTS-1] - kappa[KAPPA_10_elec_NPTS-2]) / 
      (TK[KAPPA_10_elec_NPTS-1] - TK[KAPPA_10_elec_NPTS-2]) * 
       (T-TK[KAPPA_10_elec_NPTS-1]));
  }
  else { /* Do spline */
    ans = gsl_spline_eval (spline, T, acc);
  }
  return exp(ans);
}


/* Interpolate exact results for kappa_10^pH.  The table goes up to
 * 10^5 K, but values are only accurate for T<2x10^4 K.  From
 * Furlanetto & Furlanetto 2006 */
double kappa_10_pH(double T, int flag)
{
  static double TK[KAPPA_10_pH_NPTS], kappa[KAPPA_10_pH_NPTS];
  static gsl_interp_accel *acc;
  static gsl_spline *spline;
  double ans;
  int i;
  float curr_TK, curr_kappa;
  FILE *F;

  if (flag == 1) {
    if (!(F=fopen(KAPPA_PH_FILENAME, "r"))){
      fprintf(stderr, "Unable to open the kappa_10^pH file at %s\nAborting\n", KAPPA_PH_FILENAME);
      fprintf(LOG, "Unable to open the kappa_10^pH file at %s\nAborting\n", KAPPA_PH_FILENAME);      
      return 0;
    }


    for (i=0;i<KAPPA_10_pH_NPTS;i++) {
      fscanf(F, "%f %e", &curr_TK, &curr_kappa);
      TK[i] = curr_TK;
      kappa[i] = curr_kappa;
      //      fprintf(stderr, "scanning %f\t%e\n", TK[i], kappa[i]);
    }
    fclose(F);

    for (i=0;i<KAPPA_10_pH_NPTS;i++) {
      TK[i] = log(TK[i]);
      kappa[i] = log(kappa[i]);
    }

    /* Set up spline table */
    acc   = gsl_interp_accel_alloc ();
    spline  = gsl_spline_alloc (gsl_interp_cspline, KAPPA_10_pH_NPTS);
    gsl_spline_init(spline, TK, kappa, KAPPA_10_pH_NPTS);
    return 0;
  }

  if (flag == 2) {
    /* Free memory */
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  T = log(T);
  if (T < TK[0]){ /* Use TK=1 K value if called at lower temperature */
    ans = kappa[0];
  }
  else if (T > TK[KAPPA_10_pH_NPTS-1]) { 
    /* Power law extrapolation */
    ans  = kappa[KAPPA_10_pH_NPTS-1] + 
      ((kappa[KAPPA_10_pH_NPTS-1] - kappa[KAPPA_10_pH_NPTS-2]) / 
      (TK[KAPPA_10_pH_NPTS-1] - TK[KAPPA_10_pH_NPTS-2]) * 
       (T-TK[KAPPA_10_pH_NPTS-1]));
  } else { /* Do spline */
    ans = gsl_spline_eval (spline, T, acc);
  }
  ans = exp(ans);
  return ans;
}


/********************************************************************
 ********************* Wouthuysen-Field Coupling ********************
 ********************************************************************/

/* NOTE Jalpha is by number */
double xalpha_tilde(double z, double Jalpha, double TK, double TS,
		    double delta, double xe){
  double tgp,Stilde,x;

  tgp = taugp(z,delta,xe);
  Stilde = Salpha_tilde(TK,TS,tgp);
  x = 1.66e11/(1.0+z)*Stilde*Jalpha;
  return x;
}

// Compute the Gunn-Peterson optical depth.
double taugp(double z, double delta, double xe){
  return 1.342881e-7 / hubble(z)*No*pow(1+z,3) * (1.0+delta)*(1.0-xe);
}

double Salpha_tilde(double TK, double TS, double tauGP)
{
  double xi;
  double ans;

  xi = pow(1.0e-7*tauGP/TK/TK, 1.0/3.0);
  ans = 1.0 - 0.0631789/TK + 0.115995/TK/TK - 0.401403/TS/TK;
  ans += 0.336463/TS/TK/TK;
  ans /= 1.0 + 2.98394*xi + 1.53583*xi*xi + 3.85289*xi*xi*xi;
  return ans;
}

double Tc_eff(double TK, double TS)
{
  double ans;

  ans = 1.0/TK + 0.405535/TK*(1.0/TS - 1.0/TK);
  ans = 1.0/ans;
  return ans;
}
