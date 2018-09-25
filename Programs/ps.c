#ifndef _PS_
#define _PS_

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "../Parameter_files/COSMOLOGY.H"
#include "../Parameter_files/INIT_PARAMS.H"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "cosmo_progs.c"
#include "misc.c"

/* New in v1.1 */
#define ERFC_NPTS (int) 75
#define ERFC_PARAM_DELTA (float) 0.1
static double log_erfc_table[ERFC_NPTS], erfc_params[ERFC_NPTS];
static gsl_interp_accel *erfc_acc;
static gsl_spline *erfc_spline;

#define NR_END 1
#define FREE_ARG char*

#define MM 7
#define NSTACK 50

#define FUNC(x,y,z,xx,yy) ((*func)(x,y,z,xx,yy))
#define FUNC2(x1,x2,x3,x4,x5,x6) ((*func)(x1,x2,x3,x4,x5,x6))
#define EPS2 3.0e-11

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

#define NGaussLegendre 40  //defines the number of points in the Gauss-Legendre quadrature integration

#define SPLINE_NPTS (int) 250
#define NGLhigh 100
#define NGLlow 100

#define Nhigh 200
#define Nlow 100
#define NMass 200

/* New in v1.4 - part 1 of 3 */
#define NSFR_high 200
#define NSFR_low 250
#define NGL_SFR 100 

static double log_MFspline_table[SPLINE_NPTS], MFspline_params[SPLINE_NPTS];
static double log_MFspline_table_densgtr1[SPLINE_NPTS], MFspline_params_densgtr1[SPLINE_NPTS];
static gsl_interp_accel *MFspline_acc, *MFspline_densgtr1_acc;
static gsl_spline *MF_spline, *MF_spline_densgtr1;

static double Fcoll_spline_params[SPLINE_NPTS], log_Fcoll_spline_table[SPLINE_NPTS];
static gsl_interp_accel *Fcoll_spline_acc;
static gsl_spline *Fcoll_spline;

struct parameters_gsl_int_{
    double z_obs;
    double Mval;
    double M_Feed;
    double alpha_pl;
    double del_traj_1;
    double del_traj_2;
};

struct parameters_gsl_ST_int_{
    double z_obs;
    double M_Feed;
    double alpha_pl;
};

/* New in v1.4 - part 2 of 3: begin */
static double log10_overdense_spline_SFR[NSFR_low], log10_Fcoll_spline_SFR[NSFR_low];
static gsl_interp_accel *FcollLow_spline_acc;
static gsl_spline *FcollLow_spline;

void initialiseGL_FcollSFR(int n, float M_Min, float M_Max);
void FcollSpline_SFR(float Overdensity, float *splined_value);
void initialiseFcollSFR_spline(float z, float Mmin, float Mmax, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10);

struct parameters_gsl_SFR_int_{
    double z_obs;
    double Mdrop;
    double pl_star;
    double pl_esc;
	double frac_star;
	double frac_esc;
};

struct parameters_gsl_SFR_con_int_{
    double z_obs;
    double Mval;
    double delta1;
    double delta2;
    double Mdrop;
    double pl_star;
    double pl_esc;
	double frac_star;
	double frac_esc;
};
/* New in v1.4 - part 2 of 3: end */


/* DPL test 
   Once you verify the new parametrization works correctly, you don't need this functions */
double dFdlnM_st_DPL (double lnM, void *params);
double FgtrM_st_DPL(double z, double Mmin, double MFeedback, double M_i, double alpha_pl, double beta_pl);



unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);

float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

void gauleg(float x1, float x2, float x[], float w[], int n);

double FgtrlnM_general(double lnM, void *params);
double FgtrM_general(float z, float M1, float M_Max, float M2, float MFeedback, float alpha, float delta1, float delta2);

float FgtrConditionalM_second(float z, float M1, float M2, float MFeedback, float alpha, float delta1, float delta2);
float dNdM_conditional_second(float z, float M1, float M2, float delta1, float delta2);

float FgtrConditionallnM(float M1, struct parameters_gsl_int_ parameters_gsl_int);
float GaussLegengreQuad_Fcoll(int n, float z, float M2, float MFeedback, float alpha, float delta1, float delta2);

float *Overdense_spline_gsl,*Overdense_spline_GL_high,*Fcoll_spline_gsl,*Fcoll_spline_GL_high,*xi_low,*xi_high,*wi_high,*wi_low;
float *second_derivs_low_GL,*second_derivs_high_GL,*Overdense_spline_GL_low,*Fcoll_spline_GL_low;

float *Mass_Spline, *Sigma_Spline, *dSigmadm_Spline, *second_derivs_sigma, *second_derivs_dsigma;

void initialiseSplinedSigmaM(float M_Min, float M_Max);
void initialiseGL_Fcoll(int n_low, int n_high, float M_Min, float M_Max);
void initialiseGL_FcollDblPl(int n_low, int n_high, float M_Min, float M_feedback, float M_Max);
void initialiseFcoll_spline(float z, float Mmin, float Mmax, float Mval, float MFeedback, float alphapl);

double dFdlnM_st_PL (double lnM, void *params);
double FgtrM_st_PL(double z, double Mmin, double MFeedback, double alpha_pl);

/* New in v1.4 - part 3 of 3: start */
float *Overdense_spline_SFR,*Fcoll_spline_SFR,*second_derivs_SFR;
float *xi_SFR,*wi_SFR;
double Ion_Eff_Factor_SFR(double f_star10, double f_esc10);

float FgtrConditionallnM_GL_SFR(float M, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con);
float GaussLegendreQuad_FcollSFR(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10);
double dFgtrConditionallnM_SFR(double lnM, void *params);
double FgtrConditionalM_SFR(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10);

double dFdlnM_st_SFR (double lnM, void *params);
double FgtrM_st_SFR (double z, double M, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10);
/* New in v1.4 - part 3 of 3: end */

double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;


/*****     FUNCTION PROTOTYPES     *****/
double init_ps(); /* initialize global variables, MUST CALL THIS FIRST!!! returns R_CUTOFF */
void free_ps(); /* deallocates the gsl structures from init_ps */
double splined_erfc(double); /* returns erfc for x>=0, using cubic spline in logy-x space */
double deltolindel(float del, float z); /* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double lindeltodel(float lindel, float z); /* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double power_in_k(double k); /* Returns the value of the linear power spectrum density (i.e. <|delta_k|^2>/V) at a given k mode at z=0 */
double RtoM(double); /* R in Mpc, M in Msun */
double MtoR(double); /* R in Mpc, M in Msun */
double M_J_WDM(); /* returns the "effective Jeans mass" corresponding to the gas analog of WDM ; eq. 10 in BHO 2001 */
double sheth_delc(double del, double sig);
double dNdM_st(double z, double M);
double dNdM(double z, double M);
double dnbiasdM(double M, float z, double M_o, float del_o); /* dnbiasdM */
double FgtrM(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z
double FgtrM_st(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, with Sheth-Tormen correction
double FgtrM_bias(double z, double M, double del_bias, double sig_bias);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias);/* Uses sigma parameters instead of Mass for scale */

double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias); // as above, but this version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate integral version of eq. 2 in Barkana & Loeb 2008)
double dicke(double z); //calculates the dicke growth function at redshift z
double ddickedz(double z); /* Redshift derivative of the growth function at z */
double ddickedt(double z); /* Time derivative of the growth function at z */
double sigma_z0(double M); //calculates sigma at z=0 (no dicke)
double dsigmasqdm_z0(double M); //calculates d(sigma^2)/dm at z=0 (i.e. does not include dicke growth)
double TFmdm(double k); //Eisenstien & Hu power spectrum transfer function
void TFset_parameters();
float get_R_c();  // returns R_CUTOFF
double get_M_min_ion(float z);
/***************************************/



/* Returns the minimum source mass for ionizing sources, according to user specifications */
double get_M_min_ion(float z){
  double MMIN;

  MMIN = pow(10, LOG_MASS_TURNOVER);

  // check for WDM
  if (P_CUTOFF && ( MMIN < M_J_WDM()))
    MMIN = M_J_WDM();
  //  printf("Mmin is %e\n", MMIN);
  return MMIN;
}

/* Returns the minimum source mass for x-ray sources, according to user specifications */
double get_M_min_xray(float z){
  double MMIN;
  if (X_RAY_Tvir_MIN < 9.99999e3) //neutral IGM
    MMIN = TtoM(z, X_RAY_Tvir_MIN, 1.22);
  else // ionized IGM
    MMIN = TtoM(z, X_RAY_Tvir_MIN, 0.6);

  // check for WDM
  if (P_CUTOFF && ( MMIN < M_J_WDM()))
    MMIN = M_J_WDM();
  //  printf("Mmin is %e\n", MMIN);
  return MMIN;
}


/* returns the "effective Jeans mass" in Msun
   corresponding to the gas analog of WDM ; eq. 10 in Barkana+ 2001 */
double M_J_WDM(){
  double z_eq, fudge=60;
  if (!P_CUTOFF) 
    return 0;
  z_eq = 3600*(OMm-OMb)*hlittle*hlittle/0.15;
  return fudge*3.06e8 * (1.5/g_x) * sqrt((OMm-OMb)*hlittle*hlittle/0.15) * pow(M_WDM, -4) * pow(z_eq/3000.0, 1.5);
}


/* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double deltolindel(float del, float z){
  float onepdel = 1.0+del;
  return ( 1.68647 - 1.35*pow(onepdel,-2/3.0) + 0.78785*pow(onepdel,-0.58661) - 1.12431*pow(onepdel,-0.5) )/dicke(z);
}

/* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double lindeltodel(float lindel, float z){
  float prev_lindelguess, delcrit, delguess;
  float lindelguess, delmin, delmax, epsilon = 1.0e-7;

  // set the critical density corresponding to virialization
  // this will be maximum allowed del
  delcrit = Deltac_nonlinear(z)*rho_critz(z)/(OMm*RHOcrit*pow(1+z, 3)) - 1;

  delmin = -1;
  delmax = 500;
  prev_lindelguess = -1e10;
  while (1){
    delguess = 0.5*(delmax+delmin);
    lindelguess = deltolindel(delguess, z);
    //fprintf(stderr, "%e\t%e\n", delmin, delmax);
	      //    fprintf(stderr, "%e\t%e\t%e\n\n", delguess, lindelguess, lindel);
    if ((fabs((lindelguess-lindel)/lindel) < epsilon ) ||
	(fabs(lindelguess-lindel) < epsilon ) || 
	(fabs(prev_lindelguess - lindelguess) < TINY ))// close enough, or resolution loop
      return delguess;

    if (lindelguess > lindel)
      delmax = delguess;
    else
      delmin = delguess;

    // check if we are above delcrit (see above)
    if (delmin > delcrit){
      //      printf("exced max at lindel=%e\n", lindel);
      return delcrit;
    }

    prev_lindelguess = lindelguess;
  }
}


/* R in Mpc, M in Msun */
double RtoM(double R){
  // set M according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return (4.0/3.0)*PI*pow(R,3)*(OMm*RHOcrit);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow(2*PI, 1.5) * OMm*RHOcrit * pow(R, 3);
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}

/* R in Mpc, M in Msun */
double MtoR(double M){
  // set R according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return pow(3*M/(4*PI*OMm*RHOcrit), 1.0/3.0);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow( M/(pow(2*PI, 1.5) * OMm * RHOcrit), 1.0/3.0 );
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}


/* equation (5) from jenkis et al. (2001) */
double f_jenkins(float del, double sigsq){
  if (del < 0){  fprintf(stderr, "ERROR:  In function f_jenkins del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }

  //  fprintf(stderr, "%f\t%f\n", del, sqrt(sigsq));
  return sqrt(2/PI) * del/sqrt(sigsq) * pow(E, -0.5*del*del/sigsq);
}


float get_R_c(){
  return R_CUTOFF;
}

/* sheth correction to delta crit */
double sheth_delc(double del, double sig){
  return sqrt(SHETH_a)*del*(1 + SHETH_b*pow(sig*sig/(SHETH_a*del*del), SHETH_c));
}


/* dnbiasdM */
double dnbiasdM(double M, float z, double M_o, float del_o){
  double sigsq, del, sig_one, sig_o;

  if ((M_o-M) < TINY){
    fprintf(stderr, "WARNING:  In function dnbiasdM: M must be less than M_o!\nAborting...\n");
    return -1;
  }
  del = Deltac/dicke(z) - del_o;
  if (del < 0){  fprintf(stderr, "ERROR:  In function dnbiasdM: del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }
  sig_o = sigma_z0(M_o);
  sig_one = sigma_z0(M);
  sigsq = sig_one*sig_one - sig_o*sig_o;
  return -(RHOcrit*OMm)/M /sqrt(2*PI) *del*pow(sigsq,-1.5)*pow(E, -0.5*del*del/sigsq)*dsigmasqdm_z0(M);
}


/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function with Sheth-Torman correction for ellipsoidal collapse at 
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Sheth, Mo, Torman 2001
*/
double dNdM_st(double z, double M){
  double sigma, dsigmadm, nuhat, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * dicke_growth*dicke_growth/(2.0*sigma);
  nuhat = sqrt(SHETH_a) * Deltac / sigma;

  return (-OMm*RHOcrit/M) * (dsigmadm/sigma) * sqrt(2/PI)*SHETH_A * (1+ pow(nuhat, -2*SHETH_p)) * nuhat * pow(E, -nuhat*nuhat/2.0);
}



/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function at 
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Padmanabhan, pg. 214
*/
double dNdM(double z, double M){
  double sigma, dsigmadm, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * (dicke_growth*dicke_growth/(2*sigma));

  return (-OMm*RHOcrit/M) * sqrt(2/PI) * (Deltac/(sigma*sigma)) * dsigmadm * pow(E, -(Deltac*Deltac)/(2*sigma*sigma));
}


/*
  FUNCTION FgtrM_st(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
  Uses Sheth-Torman correction
*/
double dFdlnM_st (double lnM, void *params){
  double z = *(double *)params;
  double M = exp(lnM);
  return dNdM_st(z, M) * M * M;
}
double FgtrM_st(double z, double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = 0.001; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);

  F.function = &dFdlnM_st;
  F.params = &z;
  lower_limit = log(M);
  upper_limit = log(FMAX(1e16, M*100));
			   
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);
}


/*
  FUNCTION FgtrM(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
*/
double FgtrM(double z, double M){
  double del, sig;

  del = Deltac/dicke(z); //regular spherical collapse delta
  sig = sigma_z0(M);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
 */
double FgtrM_bias(double z, double M, double del_bias, double sig_bias){
  double del, sig, sigsmallR;

  sigsmallR = sigma_z0(M);

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

   return splined_erfc(del / (sqrt(2)*sig));
}


/* Uses sigma parameters instead of Mass for scale */
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias){
  double del, sig;

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "local_FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  Calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias.
  This version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate
  integral version of eq. 2 in Barkana & Loeb 2008)
 */
double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias){
  return FgtrM_st(z, M) / FgtrM(z, M) * FgtrM_bias(z, M, del_bias, sig_bias);
}


/*
  FUNCTION dicke(z)
  Computes the dicke growth function at redshift z, i.e. the z dependance part of sigma

  References: Peebles, "Large-Scale...", pg.53 (eq. 11.16). Includes omega<=1
  Nonzero Lambda case from Liddle et al, astro-ph/9512102, eqs. 6-8.
  and quintessence case from Wang et al, astro-ph/9804015

  Normalized to dicke(z=0)=1
*/
double dicke(double z){
  double omegaM_z, dick_z, dick_0, x, x_0;
  double tiny = 1e-4;

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return 1.0/(1.0+z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    dick_z = 2.5*omegaM_z / ( 1.0/70.0 + omegaM_z*(209-omegaM_z)/140.0 + pow(omegaM_z, 4.0/7.0) );
    dick_0 = 2.5*OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );
    return dick_z / (dick_0 * (1.0+z));
  }
  else if ( (OMtot < (1+tiny)) && (fabs(OMl) < tiny) ){ //open, zero lambda case (peebles, pg. 53)
    x_0 = 1.0/(OMm+0.0) - 1.0;
    dick_0 = 1 + 3.0/x_0 + 3*log(sqrt(1+x_0)-sqrt(x_0))*sqrt(1+x_0)/pow(x_0,1.5);
    x = fabs(1.0/(OMm+0.0) - 1.0) / (1+z);
    dick_z = 1 + 3.0/x + 3*log(sqrt(1+x)-sqrt(x))*sqrt(1+x)/pow(x,1.5);
    return dick_z/dick_0;
  }
  else if ( (OMl > (-tiny)) && (fabs(OMtot-1.0) < tiny) && (fabs(wl+1) > tiny) ){
    fprintf(stderr, "IN WANG\n");
    return -1;
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}


/* redshift derivative of the growth function at z */
double ddicke_dz(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;

   return (dicke(z+dz)-dicke(z))/dz;
}

/* Time derivative of the growth function at z */
double ddickedt(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;
  double tiny = 1e-4;

   return (dicke(z+dz)-dicke(z))/dz/dtdz(z); // lazy non-analytic form getting

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return -pow(1+z,-2)/dtdz(z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    domegaMdz = omegaM_z*3/(1+z) - OMm*pow(1+z,3)*pow(OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4), -2) * (3*OMm*(1+z)*(1+z) + 4*OMr*pow(1+z,3));
    dick_0 = OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );

    ddickdz = (domegaMdz/(1+z)) * (1.0/70.0*pow(omegaM_z,-2) + 1.0/140.0 + 3.0/7.0*pow(omegaM_z, -10.0/3.0)) * pow(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0) , -2);
    ddickdz -= pow(1+z,-2)/(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0));

    return ddickdz / dick_0 / dtdz(z);
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}

/*
  FUNCTION sigma_z0(M)
  Returns the standard deviation of the normalized, density excess (delta(x)) field,
  smoothed on the comoving scale of M (see filter definitions for M<->R conversion).
  The sigma is evaluated at z=0, with the time evolution contained in the dicke(z) factor,
  i.e. sigma(M,z) = sigma_z0(m) * dicke(z)

  normalized so that sigma_z0(M->8/h Mpc) = SIGMA8 in ../Parameter_files/COSMOLOGY.H

  NOTE: volume is normalized to = 1, so this is equvalent to the mass standard deviation

  M is in solar masses

  References: Padmanabhan, pg. 210, eq. 5.107
*/
double dsigma_dk(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, kR;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  // NOTE: only use top hat for SIGMA8 normalization
  kR = k*R;
  if ( (FILTER == 0) || (sigma_norm < 0) ){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0;} // w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  return k*k*p*w*w;
}
double sigma_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;

  R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return sigma_norm * sqrt(result);
}


/*
  Returns the value of the linear power spectrum DENSITY (i.e. <|delta_k|^2>/V)
  at a given k mode linearly extrapolated to z=0
*/
double power_in_k(double k){
  double p, T, gamma, q, aa, bb, cc;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  return p*TWOPI*PI*sigma_norm*sigma_norm;
}


/*
  FUNCTION dsigmasqdm_z0(M)
  returns  d/dm (sigma^2) (see function sigma), in units of Msun^-1
*/
double dsigmasq_dm(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, dwdr, drdm, kR;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu ApJ, 1999, 511, 5
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / (hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + pow(bb*k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  kR = k * R;
  if (FILTER == 0){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0; }// w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}

    // now do d(w^2)/dm = 2 w dw/dr dr/dm
    if ( (kR) < 1.0e-10 ){  dwdr = 0;}
    else{ dwdr = 9*cos(kR)*k/pow(kR,3) + 3*sin(kR)*(1 - 3/(kR*kR))/(kR*R);}
	    //3*k*( 3*cos(kR)/pow(kR,3) + sin(kR)*(-3*pow(kR, -4) + 1/(kR*kR)) );}
      //     dwdr = -1e8 * k / (R*1e3);
    drdm = 1.0 / (4.0*PI * OMm*RHOcrit * R*R);
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
    dwdr = - k*kR * w;
    drdm = 1.0 / (pow(2*PI, 1.5) * OMm*RHOcrit * 3*R*R);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  //  printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", k, R, p, w, dwdr, drdm, dsigmadk[1]);
  return k*k*p*2*w*dwdr*drdm * d2fact;
}
double dsigmasqdm_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;

  R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);
  d2fact = M*10000/sigma_z0(M);

  F.function = &dsigmasq_dm;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return sigma_norm * sigma_norm * result /d2fact;
}



/*
  FUNCTION TFmdm is the power spectrum transfer function from Eisenstein & Hu ApJ, 1999, 511, 5
*/
double TFmdm(double k){
  double q, gamma_eff, q_eff, TF_m, q_nu;

  q = k*pow(theta_cmb,2)/omhh;
  gamma_eff=sqrt(alpha_nu) + (1.0-sqrt(alpha_nu))/(1.0+pow(0.43*k*sound_horizon, 4));     
  q_eff = q/gamma_eff;
  TF_m= log(E+1.84*beta_c*sqrt(alpha_nu)*q_eff);
  TF_m /= TF_m + pow(q_eff,2) * (14.4 + 325.0/(1.0+60.5*pow(q_eff,1.11)));
  q_nu = 3.92*q/sqrt(f_nu/N_nu);
  TF_m *= 1.0 + (1.2*pow(f_nu,0.64)*pow(N_nu,0.3+0.6*f_nu)) / 
    (pow(q_nu,-1.6)+pow(q_nu,0.8));

  //   printf("%f  %e  %f  %f  %f  %f\n",omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu);
  // printf("%f  %f  %f  %f\n", beta_c,sound_horizon,theta_cmb,z_equality);
  //printf("%f  %e  %f  %f  %f\n\n",q, k, gamma_eff, q_nu, TF_m);


  return TF_m;
}


void TFset_parameters(){
  double z_drag, R_drag, R_equality, p_c, p_cb, f_c, f_cb, f_nub, k_equality;

  z_equality = 25000*omhh*pow(theta_cmb, -4) - 1.0;
  k_equality = 0.0746*omhh/(theta_cmb*theta_cmb);

  z_drag = 0.313*pow(omhh,-0.419) * (1 + 0.607*pow(omhh, 0.674));
  z_drag = 1 + z_drag*pow(OMb*hlittle*hlittle, 0.238*pow(omhh, 0.223));
  z_drag *= 1291 * pow(omhh, 0.251) / (1 + 0.659*pow(omhh, 0.828));

  y_d = (1 + z_equality) / (1.0 + z_drag);

  R_drag = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_drag);
  R_equality = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_equality);

  sound_horizon = 2.0/3.0/k_equality * sqrt(6.0/R_equality) * 
    log( (sqrt(1+R_drag) + sqrt(R_drag+R_equality)) / (1.0 + sqrt(R_equality)) );

  p_c = -(5 - sqrt(1 + 24*(1 - f_nu-f_baryon)))/4.0;
  p_cb = -(5 - sqrt(1 + 24*(1 - f_nu)))/4.0;
  f_c = 1 - f_nu - f_baryon;
  f_cb = 1 - f_nu;
  f_nub = f_nu+f_baryon;

  alpha_nu = (f_c/f_cb) * (2*(p_c+p_cb)+5)/(4*p_cb+5.0);
  alpha_nu *= 1 - 0.553*f_nub+0.126*pow(f_nub,3);
  alpha_nu /= 1-0.193*sqrt(f_nu)+0.169*f_nu;
  alpha_nu *= pow(1+y_d, p_c-p_cb);
  alpha_nu *= 1+ (p_cb-p_c)/2.0 * (1.0+1.0/(4.0*p_c+3.0)/(4.0*p_cb+7.0))/(1.0+y_d);
  beta_c = 1.0/(1.0-0.949*f_nub);
}




double init_ps(){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;
  int i;
  double x;

  // Set cuttoff scale for WDM (eq. 4 in Barkana et al. 2001) in comoving Mpc
  R_CUTOFF = 0.201*pow((OMm-OMb)*hlittle*hlittle/0.15, 0.15)*pow(g_x/1.5, -0.29)*pow(M_WDM, -1.15);

//  fprintf(stderr, "For M_DM = %.2e keV, R_CUTOFF is: %.2e comoving Mpc\n", M_WDM, R_CUTOFF);
  if (!P_CUTOFF)
//    fprintf(stderr, "But you have selected CDM, so this is ignored\n");

  omhh = OMm*hlittle*hlittle;
  theta_cmb = T_cmb / 2.7;

  // Translate Parameters into forms GLOBALVARIABLES form
  f_nu = OMn/OMm;
  f_baryon = OMb/OMm;
  if (f_nu < TINY) f_nu = 1e-10;
  if (f_baryon < TINY) f_baryon = 1e-10;


  TFset_parameters();

  sigma_norm = -1;

  R = 8.0/hlittle;
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);

  sigma_norm = SIGMA8/sqrt(result); //takes care of volume factor


  /* initialize the lookup table for erfc */
  /*
  for (i=0; i<=ERFC_NPTS; i++){
    erfc_params[i] = i*ERFC_PARAM_DELTA;
    log_erfc_table[i] = log(erfcc(erfc_params[i]));
  }
  // Set up spline table
  erfc_acc   = gsl_interp_accel_alloc ();
  erfc_spline  = gsl_spline_alloc (gsl_interp_cspline, ERFC_NPTS);
  gsl_spline_init(erfc_spline, erfc_params, log_erfc_table, ERFC_NPTS);
  */

  return R_CUTOFF;
}

void free_ps(){
  /*    gsl_spline_free (erfc_spline);
    gsl_interp_accel_free(erfc_acc);
  */
  return;
}

double splined_erfc(double x){
  if (x < 0){
    //    fprintf(stderr, "WARNING: Negative value %e passed to splined_erfc. Returning 1\n", x);
    return 1;
  }
  return erfcc(x); // the interpolation below doesn't seem to be stable in Ts.c
  if (x > ERFC_PARAM_DELTA*(ERFC_NPTS-1))
    return erfcc(x);
  else
    return exp(gsl_spline_eval(erfc_spline, x, erfc_acc));
}

float FgtrConditionalM_second(float z, float M1, float M2, float MFeedback, float alpha, float delta1, float delta2) {
    
    return exp(M1)*pow(exp(M1)/MFeedback,alpha)*dNdM_conditional_second(z,M1,M2,delta1,delta2)/sqrt(2.*PI);
}

float dNdM_conditional_second(float z, float M1, float M2, float delta1, float delta2){
    
    float sigma1, sigma2, dsigmadm, dicke_growth,dsigma_val;
    
    M1 = exp(M1);
    M2 = exp(M2);
    
    dicke_growth = dicke(z);
    
    splint(Mass_Spline-1,Sigma_Spline-1,second_derivs_sigma-1,(int)NMass,M1,&(sigma1));
    splint(Mass_Spline-1,Sigma_Spline-1,second_derivs_sigma-1,(int)NMass,M2,&(sigma2));
    
    sigma1 = sigma1*sigma1;
    sigma2 = sigma2*sigma2;

	//printf("sigma_min^2 = %.4e sigma^2 = %.4e\n",sigma1,sigma2);
    
    splint(Mass_Spline-1,dSigmadm_Spline-1,second_derivs_dsigma-1,(int)NMass,M1,&(dsigma_val));
    
    dsigmadm = -pow(10.,dsigma_val)/(2.0*sigma1); // This is actually sigma1^{2} as calculated above, however, it should just be sigma1. It cancels with the same factor below. Why I have decided to write it like that I don't know!
    
    if((sigma1 > sigma2)) {
        
        return -(( delta1 - delta2 )/dicke_growth)*( 2.*sigma1*dsigmadm )*( exp( - ( delta1 - delta2 )*( delta1 - delta2 )/( 2.*dicke_growth*dicke_growth*( sigma1 - sigma2 ) ) ) )/(pow( sigma1 - sigma2, 1.5));
    }
    else if(sigma1==sigma2) {
        
        return -(( delta1 - delta2 )/dicke_growth)*( 2.*sigma1*dsigmadm )*( exp( - ( delta1 - delta2 )*( delta1 - delta2 )/( 2.*dicke_growth*dicke_growth*( 1.e-6 ) ) ) )/(pow( 1.e-6, 1.5));
        
    }
    else {
        return 0.;
    }
}

void gauleg(float x1, float x2, float x[], float w[], int n)
//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns arrays x[1..n] and w[1..n] of length n,
//containing the abscissas and weights of the Gauss- Legendre n-point quadrature formula.
{
    
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;
    
    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++) {
        //High precision is a good idea for this routine.
        //The roots are symmetric in the interval, so we only have to find half of them.
        //Loop over the desired roots.
        
        z=cos(3.141592654*(i-0.25)/(n+0.5));
        
        //Starting with the above approximation to the ith root, we enter the main loop of refinement by Newton’s method.
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                //Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            //p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2,
            //the polynomial of one lower order.
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > EPS2);
        x[i]=xm-xl*z;
        x[n+1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i]=w[i];
    }
}

void nrerror(char error_text[])
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
    float *v;
    v = (float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if(!v) nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
/*Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
 x1 <x2 < :: : < xN, and given values yp1 and ypn for the first derivative of the interpolating
 function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
 the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
 ypn are equal to 1e30 or larger, the routine is signaled to set the corresponding boundary
 condition for a natural spline, with zero second derivative on that boundary.*/
{
    int i,k;
    float p,qn,sig,un,*u;
    int na,nb,check;
    u=vector(1,n-1);
    if (yp1 > 0.99e30)                     // The lower boundary condition is set either to be "natural"
        y2[1]=u[1]=0.0;
    else {                                 // or else to have a specified first derivative.
        y2[1] = -0.5;
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {                              //This is the decomposition loop of the tridiagonal algorithm.
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);                //y2 and u are used for temporary
        na = 1;
        nb = 1;
        check = 0;
        while(((float)(x[i+na*1]-x[i-nb*1])==(float)0.0)) {
            check = check + 1;
            if(check%2==0) {
                na = na + 1;
            }
            else {
                nb = nb + 1;
            }
            sig=(x[i]-x[i-1])/(x[i+na*1]-x[i-nb*1]);
        }
        p=sig*y2[i-1]+2.0;                                //storage of the decomposed
        y2[i]=(sig-1.0)/p;                                //  factors.
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
        
        if(((float)(x[i+1]-x[i])==(float)0.0) || ((float)(x[i]-x[i-1])==(float)0.0)) {
            na = 0;
            nb = 0;
            check = 0;
            while((float)(x[i+na*1]-x[i-nb])==(float)(0.0) || ((float)(x[i+na]-x[i-nb*1])==(float)0.0)) {
                check = check + 1;
                if(check%2==0) {
                    na = na + 1;
                }
                else {
                    nb = nb + 1;
                }
            }
            u[i]=(y[i+1]-y[i])/(x[i+na*1]-x[i-nb]) - (y[i]-y[i-1])/(x[i+na]-x[i-nb*1]);
            
            u[i]=(6.0*u[i]/(x[i+na*1]-x[i-nb*1])-sig*u[i-1])/p;
            
        }
    }
    if (ypn > 0.99e30)                        //The upper boundary condition is set either to be "natural"
        qn=un=0.0;
    else {                                    //or else to have a specified first derivative.
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    
    for (k=n-1;k>=1;k--) {                      //This is the backsubstitution loop of the tridiagonal
        y2[k]=y2[k]*y2[k+1]+u[k];               //algorithm.
    }
    free_vector(u,1,n-1);
}


void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
/*Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai's in order),
 and given the array y2a[1..n], which is the output from spline above, and given a value of
 x, this routine returns a cubic-spline interpolated value y.*/
{
    void nrerror(char error_text[]);
    int klo,khi,k;
    float h,b,a;
    klo=1;                                                  // We will find the right place in the table by means of
    khi=n;                                                  //bisection. This is optimal if sequential calls to this
    while (khi-klo > 1) {                                   //routine are at random values of x. If sequential calls
        k=(khi+klo) >> 1;                                     //are in order, and closely spaced, one would do better
        if (xa[k] > x) khi=k;                                 //to store previous values of klo and khi and test if
        else klo=k;                                           //they remain appropriate on the next call.
    }                                                           // klo and khi now bracket the input value of x.
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint");    //The xa's must be distinct.
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;                                            //Cubic spline polynomial is now evaluated.
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
    unsigned long *v;
    v = (unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
    if(!v) nrerror("allocation failure in lvector()");
    return v - nl + NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

double FgtrlnM_general(double lnM, void *params) {
    
    struct parameters_gsl_int_ vals = *(struct parameters_gsl_int_ *)params;
    
    float z = vals.z_obs;
    float M2 = vals.Mval;
    float MFeedback = vals.M_Feed;
    float alpha = vals.alpha_pl;
    float delta1 = vals.del_traj_1;
    float delta2 = vals.del_traj_2;
    
    return FgtrConditionalM_second(z,lnM,M2,MFeedback,alpha,delta1,delta2);
}
double FgtrM_general(float z, float M1, float M_Max, float M2, float MFeedback, float alpha, float delta1, float delta2) {
    
    double result, error, lower_limit, upper_limit;
    
    double rel_tol = 0.01;
    
    int size;
    size = 1000;
    
    //    printf("delta1 = %e Deltac = %e\n",delta1,Deltac);
    
    if((float)delta1==(float)Deltac) {
        
        gsl_function Fx;
        
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (size);
        
        Fx.function = &FgtrlnM_general;
        
        struct parameters_gsl_int_  parameters_gsl_int = {
            .z_obs = z,
            .Mval = M2,
            .M_Feed = MFeedback,
            .alpha_pl = alpha,
            .del_traj_1 = delta1,
            .del_traj_2 = delta2
        };
        
        Fx.params = &parameters_gsl_int;
        
        lower_limit = M1;
        upper_limit = M_Max;
        
        gsl_integration_qag (&Fx, lower_limit, upper_limit, 0, rel_tol, size, GSL_INTEG_GAUSS15, w, &result, &error);
        //        gsl_integration_qag (&Fx, lower_limit, upper_limit, 0, rel_tol, size, GSL_INTEG_GAUSS61, w, &result, &error);
        gsl_integration_workspace_free (w);
        
        if(delta2 > delta1) {
            return 1.;
        }
        else {
            return result;
        }
    }
}

float FgtrConditionallnM(float M1, struct parameters_gsl_int_ parameters_gsl_int) {
    
    float z = parameters_gsl_int.z_obs;
    float M2 = parameters_gsl_int.Mval;
    float MFeedback = parameters_gsl_int.M_Feed;
    float alpha = parameters_gsl_int.alpha_pl;
    float delta1 = parameters_gsl_int.del_traj_1;
    float delta2 = parameters_gsl_int.del_traj_2;
    
    return exp(M1)*pow(exp(M1)/MFeedback,alpha)*dNdM_conditional_second(z,M1,M2,delta1,delta2)/sqrt(2.*PI);
}

float GaussLegengreQuad_Fcoll(int n, float z, float M2, float MFeedback, float alpha, float delta1, float delta2)
{
    //Performs the Gauss-Legendre quadrature.
    int i;
    
    float integrand,x;
    integrand = 0.0;
    
    struct parameters_gsl_int_  parameters_gsl_int = {
        .z_obs = z,
        .Mval = M2,
        .M_Feed = MFeedback,
        .alpha_pl = alpha,
        .del_traj_1 = delta1,
        .del_traj_2 = delta2
    };
    
    if(delta2>delta1) {
        return 1.;
    }
    else {
        for(i=1;i<(n+1);i++) {
            x = xi_low[i];
            integrand += wi_low[i]*FgtrConditionallnM(x,parameters_gsl_int);
        }
        return integrand;
    }
}

/*
 FUNCTION FgtrM_st(z, M)
 Computes the fraction of mass contained in haloes with mass > M at redshift z
 Uses Sheth-Torman correction
 */
double dFdlnM_st_PL (double lnM, void *params){
    
    struct parameters_gsl_ST_int_ vals = *(struct parameters_gsl_ST_int_ *)params;
    
    double M = exp(lnM);
    float z = vals.z_obs;
    float MFeedback = vals.M_Feed;
    float alpha = vals.alpha_pl;
    
    return dNdM_st(z, M) * M * M * pow((M/MFeedback),alpha);
}
double FgtrM_st_PL(double z, double Mmin, double MFeedback, double alpha_pl){
    
    double result_lower, result_upper, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol  = 0.001; //<- relative tolerance
    gsl_integration_workspace * w_lower
    = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * w_upper
    = gsl_integration_workspace_alloc (1000);
    
    struct parameters_gsl_ST_int_  parameters_gsl_ST_lower = {
        .z_obs = z,
        .M_Feed = MFeedback,
        .alpha_pl = alpha_pl,
    };
    
    F.function = &dFdlnM_st_PL;
    F.params = &parameters_gsl_ST_lower;
    lower_limit = log(Mmin);
    upper_limit = log(1e16);
    
    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                         1000, GSL_INTEG_GAUSS61, w_lower, &result_lower, &error);
    gsl_integration_workspace_free (w_lower);
    
    return (result_lower) / (OMm*RHOcrit);
}

void initialiseSplinedSigmaM(float M_Min, float M_Max)
{
    int i;
    float Mass;
    
    Mass_Spline = calloc(NMass,sizeof(float));
    Sigma_Spline = calloc(NMass,sizeof(float));
    dSigmadm_Spline = calloc(NMass,sizeof(float));
    second_derivs_sigma = calloc(NMass,sizeof(float));
    second_derivs_dsigma = calloc(NMass,sizeof(float));
    
    for(i=0;i<NMass;i++) {
        Mass_Spline[i] = pow(10., log10(M_Min) + (float)i/(NMass-1)*( log10(M_Max) - log10(M_Min) ) );
        Sigma_Spline[i] = sigma_z0(Mass_Spline[i]);
        dSigmadm_Spline[i] = log10(-dsigmasqdm_z0(Mass_Spline[i]));
    }
    spline(Mass_Spline-1,Sigma_Spline-1,NMass,0,0,second_derivs_sigma-1);
    spline(Mass_Spline-1,dSigmadm_Spline-1,NMass,0,0,second_derivs_dsigma-1);
}

void initialiseGL_Fcoll(int n_low, int n_high, float M_Min, float M_Max)
{
    //calculates the weightings and the positions for Gauss-Legendre quadrature.
    
    gauleg(log(M_Min),log(M_Max),xi_low,wi_low,n_low);
    
    gauleg(log(M_Min),log(M_Max),xi_high,wi_high,n_high);
    
}

void initialiseFcoll_spline(float z, float Mmin, float Mmax, float Mval, float MFeedback, float alphapl)
{
    double overdense_val,overdense_small_low,overdense_small_high,overdense_large_low,overdense_large_high;
    int i;
    
    overdense_large_high = Deltac;
    overdense_large_low = 1.5;
    overdense_small_high = 1.5;
    overdense_small_low = -1. + 9.e-8;
    
    Fcoll_spline_acc   = gsl_interp_accel_alloc ();
    Fcoll_spline  = gsl_spline_alloc (gsl_interp_cspline, SPLINE_NPTS);
    
    for (i=0;i<SPLINE_NPTS;i++){
        overdense_val = log10(1.+overdense_small_low) + (float)i/(SPLINE_NPTS-1.)*(log10(1.+overdense_small_high) - log10(1.+overdense_small_low));
        
        log_Fcoll_spline_table[i] = log10(GaussLegengreQuad_Fcoll(NGLlow,z,log(Mval),MFeedback,alphapl,Deltac,pow(10.,overdense_val)-1.));
        Fcoll_spline_params[i] = overdense_val;
        
        if(log_Fcoll_spline_table[i]<-40.) {
            log_Fcoll_spline_table[i] = -40.;
        }
    }
    gsl_spline_init(Fcoll_spline, Fcoll_spline_params, log_Fcoll_spline_table, SPLINE_NPTS);
    
    for(i=0;i<Nhigh;i++) {
        Overdense_spline_GL_high[i] = overdense_large_low + (float)i/((float)Nhigh-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_GL_high[i] = FgtrM_general(z,log(Mmin),log(Mmax),log(Mval),MFeedback,alphapl,Deltac,Overdense_spline_GL_high[i]);
        
        if(Fcoll_spline_GL_high[i]<0.) {
            Fcoll_spline_GL_high[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_GL_high-1,Fcoll_spline_GL_high-1,Nhigh,0,0,second_derivs_high_GL-1);
}

void FcollSpline(float Overdensity, float *splined_value)
{
    int i;
    float returned_value;
    
    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(Fcoll_spline, log10(Overdensity+1.), Fcoll_spline_acc);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<Deltac) {
            splint(Overdense_spline_GL_high-1,Fcoll_spline_GL_high-1,second_derivs_high_GL-1,(int)Nhigh,Overdensity,&(returned_value));
        }
        else {
            returned_value = 1.;
        }
    }
    *splined_value = returned_value;
}

/* New in v1.4: begin */

double Ion_Eff_Factor_SFR(double f_star10, double f_esc10){
// Ionizing efficiency factor using the new parametrization

    double ionizing_photon_num = N_GAMMA_UV, baryon_frac = 1.;

    return ionizing_photon_num * baryon_frac * f_star10 * f_esc10;
}

/*
 FUNCTION FgtrM_st_SFR(z, Mdrop)
 Computes the fraction of mass contained in haloes with mass > M at redshift z
 Uses Sheth-Torman correction
*/
double dFdlnM_st_SFR(double lnM, void *params){
    struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
	double Fstar10 = vals.frac_star;
	double Fesc10 = vals.frac_esc;

	double Fstar = pow(M/1e10,Alpha_star);
	double Fesc = pow(M/1e10,Alpha_esc);

	if (Fstar*Fstar10 > 1.) Fstar = 1./Fstar10;
	else if (Fstar*Fstar10 < FRACT_FLOAT_ERR) Fstar = FRACT_FLOAT_ERR;

	if (Fesc*Fesc10 > 1.) Fesc = 1./Fesc10;
	else if (Fesc*Fstar10 < FRACT_FLOAT_ERR) Fesc = FRACT_FLOAT_ERR;

    //return dNdM_st(z,M) * M * M * exp(-MassTurnover/M) * pow(M/1e10,Alpha_star) * pow(M/1e10,Alpha_esc);
    return dNdM_st(z,M) * M * M * exp(-MassTurnover/M) * Fstar * Fesc;
}

double FgtrM_st_SFR(double z, double M, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10){
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
        .z_obs = z,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
		.frac_star = Fstar10,
		.frac_esc = Fesc10,
    };

    F.function = &dFdlnM_st_SFR;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M);
    upper_limit = log(FMAX(1e16, M*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    return result / (OMm*RHOcrit);

}

void initialiseGL_FcollSFR(int n, float M_Min, float M_Max){
    //calculates the weightings and the positions for Gauss-Legendre quadrature.
	//printf("%d %e %e\n",n,M_Min,M_Max);
    gauleg(log(M_Min),log(M_Max),xi_SFR,wi_SFR,n);

}

float FgtrConditionallnM_GL_SFR(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Fesc10 = parameters_gsl_SFR_con.frac_esc;

	float Fstar = pow(M/1e10,Alpha_star);
	float Fesc = pow(M/1e10,Alpha_esc);

	if (Fstar > 1.) Fstar = 1./Fstar10;
	else if (Fstar < 0.) Fstar = FRACT_FLOAT_ERR;

	if (Fesc > 1.) Fesc = 1./Fesc10;
	else if (Fstar < 0.) Fesc = FRACT_FLOAT_ERR;

    //return M*exp(-MassTurnover/M)*pow(M/1e10,Alpha_star)*pow(M/1e10,Alpha_esc)*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);

}

float GaussLegendreQuad_FcollSFR(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
		.frac_star = Fstar10,
		.frac_esc = Fesc10
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_SFR(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

double dFgtrConditionallnM_SFR(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
	double Fstar10 = vals.frac_star;
	double Fesc10 = vals.frac_esc;

	double Fstar = pow(M/1e10,Alpha_star);
	double Fesc = pow(M/1e10,Alpha_esc);

	if (Fstar > 1./Fstar10) Fstar = 1./Fstar10;
	else if (Fstar < 0.) Fstar = FRACT_FLOAT_ERR;

	if (Fesc > 1./Fstar10) Fesc = 1./Fesc10;
	else if (Fstar < 0.) Fesc = FRACT_FLOAT_ERR;

    //return M*exp(-MassTurnover/M)*pow(M/1e10,Alpha_star)*pow(M/1e10,Alpha_esc)*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
}

double FgtrConditionalM_SFR(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
    //printf("FgtrConditionalM_SFR start\n");

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
		.frac_star = Fstar10,
		.frac_esc = Fesc10
    };

    F.function = &dFgtrConditionallnM_SFR;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    //printf("FgtrConditionalM_SFR end\n");
    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

void initialiseFcollSFR_spline(float z, float Mmin, float Mmax, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10){
    //double overdense_val,overdense_small_low,overdense_small_high,overdense_large_low,overdense_large_high;
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;

    FcollLow_spline_acc = gsl_interp_accel_alloc ();
    FcollLow_spline = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);

    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_SFR[i] = overdense_val;
        log10_Fcoll_spline_SFR[i] = log10(GaussLegendreQuad_FcollSFR(NGL_SFR,z,log(Mmax),Deltac,pow(10.,overdense_val)-1.,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10));

        if(log10_Fcoll_spline_SFR[i] < -40.){
            log10_Fcoll_spline_SFR[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline, log10_overdense_spline_SFR, log10_Fcoll_spline_SFR, NSFR_low);


    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_SFR[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_SFR[i] = FgtrConditionalM_SFR(z,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10);

        if(Fcoll_spline_SFR[i]<0.) {
            Fcoll_spline_SFR[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_SFR-1,Fcoll_spline_SFR-1,NSFR_high,0,0,second_derivs_SFR-1);
}

void FcollSpline_SFR(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline, log10(Overdensity+1.), FcollLow_spline_acc);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<Deltac) {
            splint(Overdense_spline_SFR-1,Fcoll_spline_SFR-1,second_derivs_SFR-1,NSFR_high,Overdensity,&(returned_value));
        }
        else {
            returned_value = 1.;
        }
    }
	if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

/* New in v1.4: end */

/* DPL test 
   Once you verify the new parametrization works correctly, you don't need this functions */

double dFdlnM_st_DPL (double lnM, void *params){

    struct parameters_gsl_ST_int_ vals = *(struct parameters_gsl_ST_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    float MFeedback = vals.M_Feed;
    float alpha = vals.alpha_pl;
    
    return dNdM_st(z, M) * M * M * pow((M/MFeedback),alpha);
}

double FgtrM_st_DPL(double z, double Mmin, double MFeedback, double M_i, double alpha_pl, double beta_pl){

    double result_lower, result_upper, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol  = 0.001; //<- relative tolerance
    gsl_integration_workspace * w_lower
    = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * w_upper
    = gsl_integration_workspace_alloc (1000);

	if (M_i >= MFeedback) {
		result_lower = 0.;

	    struct parameters_gsl_ST_int_  parameters_gsl_ST_upper = {
	        .z_obs = z,
	        .M_Feed = MFeedback,
	        .alpha_pl = beta_pl,
	    };

	    F.function = &dFdlnM_st_DPL;
	    F.params = &parameters_gsl_ST_upper;
	    lower_limit = log(M_i);
	    upper_limit = log(1e16);

	    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                         1000, GSL_INTEG_GAUSS61, w_upper, &result_upper, &error);
	    gsl_integration_workspace_free (w_upper);
	}
	else if (M_i < MFeedback) {
	    struct parameters_gsl_ST_int_  parameters_gsl_ST_lower = {
	        .z_obs = z,
	        .M_Feed = MFeedback,
	        .alpha_pl = alpha_pl,
	    };

	    F.function = &dFdlnM_st_DPL;
	    F.params = &parameters_gsl_ST_lower;
	    lower_limit = log(M_i);//log(Mmin);
	    upper_limit = log(MFeedback);

	    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                         1000, GSL_INTEG_GAUSS61, w_lower, &result_lower, &error);
	    gsl_integration_workspace_free (w_lower);
	

	    struct parameters_gsl_ST_int_  parameters_gsl_ST_upper = {
	        .z_obs = z,
	        .M_Feed = MFeedback,
	        .alpha_pl = beta_pl,
	    };

	    F.function = &dFdlnM_st_DPL;
	    F.params = &parameters_gsl_ST_upper;
	    lower_limit = log(MFeedback);
	    upper_limit = log(1e16);

	    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                         1000, GSL_INTEG_GAUSS61, w_upper, &result_upper, &error);
	    gsl_integration_workspace_free (w_upper);
	}

    return (result_lower + result_upper) / (OMm*RHOcrit);
}




#endif
