#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "Variables.h"
#include "bubble_helper_progs.c"
#include "filter.c"
#include "gsl/gsl_sf_erf.h"

/* Throughout this and other 21cmMC drivers, the use of Deltac is not for checking against
    the z=0 collapsed fraction, but rather, it is used as a threshold for the validity of the 
    collapse fraction expression. This expression is only valid up to Deltac
 */

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
#define SPLINE_NPTS_DBLPL (int) 500
//#define SPLINE_NPTS_DBLPL (int) 350
#define NGLhigh 100
#define NGLlow 100
//#define NGLhigh 60
//#define NGLlow 60

//#define NGLhigh_DblPl 60
//#define NGLlow_DblPl 60

#define NGLhigh_DblPl 100
#define NGLlow_DblPl 100
//#define NGLhigh_DblPl 60
//#define NGLlow_DblPl 60
//#define NGLhigh_DblPl 40
//#define NGLlow_DblPl 40

#define Nhigh 200
//#define Nhigh 75
#define NhighDblpl 200
//#define NhighDblpl 75
#define Nlow 100
#define NMass 200

static double log_MFspline_table[SPLINE_NPTS], MFspline_params[SPLINE_NPTS];
static double log_MFspline_table_densgtr1[SPLINE_NPTS], MFspline_params_densgtr1[SPLINE_NPTS];
static gsl_interp_accel *MFspline_acc, *MFspline_densgtr1_acc;
static gsl_spline *MF_spline, *MF_spline_densgtr1;

static double Fcoll_spline_params[SPLINE_NPTS], log_Fcoll_spline_table[SPLINE_NPTS], Fcoll_spline_params_dblpl_alpha[SPLINE_NPTS_DBLPL], log_Fcoll_spline_table_dblpl_alpha[SPLINE_NPTS_DBLPL];
static double Fcoll_spline_params_dblpl_beta[SPLINE_NPTS_DBLPL], log_Fcoll_spline_table_dblpl_beta[SPLINE_NPTS_DBLPL];
static gsl_interp_accel *Fcoll_spline_acc, *Fcoll_splineDblpl_alpha_acc, *Fcoll_splineDblpl_beta_acc;
static gsl_spline *Fcoll_spline, *Fcoll_splineDblpl_alpha, *Fcoll_splineDblpl_beta;

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
float GaussLegengreQuad_Fcoll(int type, int n, float z, float M2, float MFeedback, float alpha, float delta1, float delta2);

float *Overdense_spline_gsl,*Overdense_spline_GL_high,*Fcoll_spline_gsl,*Fcoll_spline_GL_high,*xi_low,*xi_high,*wi_high,*wi_low;
float *xi_low_dblpl_alpha,*xi_high_dblpl_alpha,*wi_high_dblpl_alpha,*wi_low_dblpl_alpha,*xi_low_dblpl_beta,*xi_high_dblpl_beta,*wi_high_dblpl_beta,*wi_low_dblpl_beta;
float *second_derivs_low_GL,*second_derivs_high_GL,*Overdense_spline_GL_low,*Fcoll_spline_GL_low;
float *Fcoll_spline_GL_high_dblpl_alpha,*second_derivs_high_GL_dblpl_alpha;
float *Fcoll_spline_GL_high_dblpl_beta,*second_derivs_high_GL_dblpl_beta,*Overdense_spline_GL_high_dblpl;

float *Mass_Spline, *Sigma_Spline, *dSigmadm_Spline, *second_derivs_sigma, *second_derivs_dsigma;

void initialiseSplinedSigmaM(float M_Min, float M_Max);
void initialiseGL_Fcoll(int n_low, int n_high, float M_Min, float M_Max);
void initialiseGL_FcollDblPl(int n_low, int n_high, float M_Min, float M_feedback, float M_Max);
void initialiseFcoll_spline(float z, float Mmin, float Mmax, float Mval, float MFeedback, float alphapl);
void initialiseFcoll_spline_DblPl(float z, float Mmin, float Mmax, float Mval, float alphapl, float betapl, float Mfeedback);

double dFdlnM_st_DPL (double lnM, void *params);
double FgtrM_st_DPL(double z, double Mmin, double MFeedback, double alpha_pl, double beta_pl);

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

void initialiseGL_FcollDblPl(int n_low, int n_high, float M_Min, float M_feedback, float M_Max)
{
    //calculates the weightings and the positions for Gauss-Legendre quadrature.
    
    // For first power law index alpha
    gauleg(log(M_Min),log(M_feedback),xi_low_dblpl_alpha,wi_low_dblpl_alpha,n_low);
    
    gauleg(log(M_Min),log(M_feedback),xi_high_dblpl_alpha,wi_high_dblpl_alpha,n_high);
    
    // For second power law index beta
    gauleg(log(M_feedback),log(M_Max),xi_low_dblpl_beta,wi_low_dblpl_beta,n_low);
    
    gauleg(log(M_feedback),log(M_Max),xi_high_dblpl_beta,wi_high_dblpl_beta,n_high);
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

        log_Fcoll_spline_table[i] = log10(GaussLegengreQuad_Fcoll(1,NGLlow,z,log(Mval),MFeedback,alphapl,Deltac,pow(10.,overdense_val)-1.));
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

void initialiseFcoll_spline_DblPl(float z, float Mmin, float Mmax, float Mval, float alphapl, float betapl, float Mfeedback)
{
    double overdense_val,overdense_small_low,overdense_small_high,overdense_large_low,overdense_large_high;
    int i;
    
    overdense_large_high = Deltac;
    overdense_large_low = 1.5;
    overdense_small_high = 1.5;
    overdense_small_low = -1. + 9.e-8;
//    printf("Overdens boundary = %1.12e\n",overdense_small_low);
    
    // First set up spline for alpha
    
    Fcoll_splineDblpl_alpha_acc   = gsl_interp_accel_alloc ();
    Fcoll_splineDblpl_alpha  = gsl_spline_alloc (gsl_interp_cspline, SPLINE_NPTS_DBLPL);
    
    if(Mfeedback<=Mmin) {
        //This condition occurs if Mmin > Mfeedback, which should never be the case. Nevertheless, set to very small number in this case
        for (i=0;i<SPLINE_NPTS_DBLPL;i++){
            overdense_val = log10(1.+overdense_small_low) + (float)i/(SPLINE_NPTS_DBLPL-1.)*(log10(1.+overdense_small_high) - log10(1.+overdense_small_low));
            Fcoll_spline_params_dblpl_alpha[i] = overdense_val;
            
            log_Fcoll_spline_table_dblpl_alpha[i] = -45.;
        }
    }
    else {

        for (i=0;i<SPLINE_NPTS_DBLPL;i++){
            overdense_val = log10(1.+overdense_small_low) + (float)i/(SPLINE_NPTS_DBLPL-1.)*(log10(1.+overdense_small_high) - log10(1.+overdense_small_low));
            
            log_Fcoll_spline_table_dblpl_alpha[i] = log10(GaussLegengreQuad_Fcoll(3,NGLlow_DblPl,z,log(Mval),Mfeedback,alphapl,Deltac,pow(10.,overdense_val)-1.));
            Fcoll_spline_params_dblpl_alpha[i] = overdense_val;
            
            if(log_Fcoll_spline_table_dblpl_alpha[i]<-45.) {
                log_Fcoll_spline_table_dblpl_alpha[i] = -45.;
            }
        }
    }
    gsl_spline_init(Fcoll_splineDblpl_alpha, Fcoll_spline_params_dblpl_alpha, log_Fcoll_spline_table_dblpl_alpha, SPLINE_NPTS_DBLPL);

    for(i=0;i<NhighDblpl;i++) {
        Overdense_spline_GL_high_dblpl[i] = overdense_large_low + (float)i/((float)NhighDblpl-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_GL_high_dblpl_alpha[i] = FgtrM_general(z,log(Mmin),log(Mfeedback),log(Mval),Mfeedback,alphapl,Deltac,Overdense_spline_GL_high_dblpl[i]);
        
        //This condition occurs if Mmin > Mfeedback, which should never be the case. Nevertheless, set to very small number in this case
        if(Fcoll_spline_GL_high_dblpl_alpha[i]<0.) {
            Fcoll_spline_GL_high_dblpl_alpha[i]=pow(10.,-45.0);
        }
        
    }
    spline(Overdense_spline_GL_high_dblpl-1,Fcoll_spline_GL_high_dblpl_alpha-1,NhighDblpl,0,0,second_derivs_high_GL_dblpl_alpha-1);
    
    // Now set up spline for beta
    
    Fcoll_splineDblpl_beta_acc   = gsl_interp_accel_alloc ();
    Fcoll_splineDblpl_beta  = gsl_spline_alloc (gsl_interp_cspline, SPLINE_NPTS_DBLPL);
    
    for (i=0;i<SPLINE_NPTS_DBLPL;i++){
        overdense_val = log10(1.+overdense_small_low) + (float)i/(SPLINE_NPTS_DBLPL-1.)*(log10(1.+overdense_small_high) - log10(1.+overdense_small_low));
            
        log_Fcoll_spline_table_dblpl_beta[i] = log10(GaussLegengreQuad_Fcoll(5,NGLlow_DblPl,z,log(Mval),Mfeedback,betapl,Deltac,pow(10.,overdense_val)-1.));
        Fcoll_spline_params_dblpl_beta[i] = overdense_val;
            
        if(log_Fcoll_spline_table_dblpl_beta[i]<-45.) {
            log_Fcoll_spline_table_dblpl_beta[i] = -45.;
        }
    }
    gsl_spline_init(Fcoll_splineDblpl_beta, Fcoll_spline_params_dblpl_beta, log_Fcoll_spline_table_dblpl_beta, SPLINE_NPTS_DBLPL);
    
    for(i=0;i<NhighDblpl;i++) {
        Overdense_spline_GL_high_dblpl[i] = overdense_large_low + (float)i/((float)NhighDblpl-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_GL_high_dblpl_beta[i] = FgtrM_general(z,log(Mfeedback),log(Mmax),log(Mval),Mfeedback,betapl,Deltac,Overdense_spline_GL_high_dblpl[i]);
    }
    spline(Overdense_spline_GL_high_dblpl-1,Fcoll_spline_GL_high_dblpl_beta-1,NhighDblpl,0,0,second_derivs_high_GL_dblpl_beta-1);
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

void FcollSplineDblPl(float Overdensity, float Mfeedback, float M2, float *splined_value)
{
    int i;
    float returned_value,Fcoll_integrand;
    
//    printf("Overdensity = %1.12e\n",Overdensity);
//    printf("Density contrast = %1.12e\n",Overdensity+1.);
    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            
            if (M2 > Mfeedback) {

                returned_value = gsl_spline_eval(Fcoll_splineDblpl_alpha, log10(Overdensity+1.), Fcoll_splineDblpl_alpha_acc);
                returned_value = pow(10.,returned_value);
                returned_value = returned_value + pow(10., gsl_spline_eval(Fcoll_splineDblpl_beta, log10(Overdensity+1.), Fcoll_splineDblpl_beta_acc) );
            }
            else {
                returned_value = pow(10., gsl_spline_eval(Fcoll_splineDblpl_alpha, log10(Overdensity+1.), Fcoll_splineDblpl_alpha_acc) );
            }
        }
    }
    else {
        if(Overdensity<Deltac) {
            
            if (M2 > Mfeedback) {
                splint(Overdense_spline_GL_high_dblpl-1,Fcoll_spline_GL_high_dblpl_alpha-1,second_derivs_high_GL_dblpl_alpha-1,(int)NhighDblpl,Overdensity,&(Fcoll_integrand));
            
                returned_value = Fcoll_integrand;
                splint(Overdense_spline_GL_high_dblpl-1,Fcoll_spline_GL_high_dblpl_beta-1,second_derivs_high_GL_dblpl_beta-1,(int)NhighDblpl,Overdensity,&(Fcoll_integrand));
            
                returned_value = returned_value + Fcoll_integrand;
            }
            else {
                splint(Overdense_spline_GL_high_dblpl-1,Fcoll_spline_GL_high_dblpl_alpha-1,second_derivs_high_GL_dblpl_alpha-1,(int)NhighDblpl,Overdensity,&(Fcoll_integrand));
                
                returned_value = Fcoll_integrand;
            }
        }
        else {
            returned_value = 1.;
        }
    }
    
    *splined_value = returned_value;
}

/**** Arrays declared and used *****/

void init_21cmMC_arrays() {
    
    Overdense_spline_GL_low = calloc(Nlow,sizeof(float));
    Fcoll_spline_GL_low = calloc(Nlow,sizeof(float));
    second_derivs_low_GL = calloc(Nlow,sizeof(float));
    Overdense_spline_GL_high = calloc(Nhigh,sizeof(float));
    Fcoll_spline_GL_high = calloc(Nhigh,sizeof(float));
    second_derivs_high_GL = calloc(Nhigh,sizeof(float));
    Overdense_spline_GL_high_dblpl = calloc(NhighDblpl,sizeof(float));
    Fcoll_spline_GL_high_dblpl_alpha = calloc(NhighDblpl,sizeof(float));
    second_derivs_high_GL_dblpl_alpha= calloc(NhighDblpl,sizeof(float));
    Fcoll_spline_GL_high_dblpl_beta = calloc(NhighDblpl,sizeof(float));
    second_derivs_high_GL_dblpl_beta = calloc(NhighDblpl,sizeof(float));
    
    xH = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    Fcoll = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
    
    if(USE_TS_IN_21CM) {
        xe_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // Perhaps not necessary (within USE_TS_IN_21CM)
        xe_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // Perhaps not necessary (within USE_TS_IN_21CM)
    }
    if(USE_HALO_FIELD) {
        M_coll_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // Perhaps not necessary (within USE_HALO_FIELD)
        M_coll_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // Perhaps not necessary (within USE_HALO_FIELD)
    }
    
    deltax_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_unfiltered_original = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
    delta_T = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    v = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
    
    k_factor = 1.25;
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
    
    deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

    xi_low_dblpl_alpha = calloc((NGLlow_DblPl+1),sizeof(float));
    wi_low_dblpl_alpha = calloc((NGLlow_DblPl+1),sizeof(float));
    
    xi_high_dblpl_alpha = calloc((NGLhigh_DblPl+1),sizeof(float));
    wi_high_dblpl_alpha = calloc((NGLhigh_DblPl+1),sizeof(float));
    
    xi_low_dblpl_beta = calloc((NGLlow_DblPl+1),sizeof(float));
    wi_low_dblpl_beta = calloc((NGLlow_DblPl+1),sizeof(float));
    
    xi_high_dblpl_beta = calloc((NGLhigh_DblPl+1),sizeof(float));
    wi_high_dblpl_beta = calloc((NGLhigh_DblPl+1),sizeof(float));
    
    xi_low = calloc((NGLlow+1),sizeof(float));
    wi_low = calloc((NGLlow+1),sizeof(float));
    
    xi_high = calloc((NGLhigh+1),sizeof(float));
    wi_high = calloc((NGLhigh+1),sizeof(float));
    
}

void destroy_21cmMC_arrays() {
    
    fftwf_free(xH);
    
    if(USE_TS_IN_21CM) {
        fftwf_free(xe_unfiltered);
        fftwf_free(xe_filtered);
    }
    if(USE_HALO_FIELD) {
        fftwf_free(M_coll_unfiltered);
        fftwf_free(M_coll_filtered);
    }
    
    fftwf_free(deltax_unfiltered);
    fftwf_free(deltax_unfiltered_original);
    fftwf_free(deltax_filtered);
    free(deltax);
    free(Fcoll);
    free(delta_T);
    free(v);
        
    free(p_box);
    free(k_ave);
    free(in_bin_ct);
    
    fftwf_free(deldel_T);
    
    free(Overdense_spline_GL_low);
    free(Fcoll_spline_GL_low);
    free(second_derivs_low_GL);
    free(Overdense_spline_GL_high);
    free(Fcoll_spline_GL_high);
    free(second_derivs_high_GL);
    free(Overdense_spline_GL_high_dblpl);
    free(Fcoll_spline_GL_high_dblpl_alpha);
    free(second_derivs_high_GL_dblpl_alpha);
    free(Fcoll_spline_GL_high_dblpl_beta);
    free(second_derivs_high_GL_dblpl_beta);

    
    free(p_box);
    free(k_ave);
    free(in_bin_ct);
    
    free(deldel_T);
    
    free(xi_low_dblpl_alpha);
    free(wi_low_dblpl_alpha);
    
    free(xi_high_dblpl_alpha);
    free(wi_high_dblpl_alpha);
    
    free(xi_low_dblpl_beta);
    free(wi_low_dblpl_beta);
    
    free(xi_high_dblpl_beta);
    free(wi_high_dblpl_beta);
    
    free(xi_low);
    free(wi_low);
    
    free(xi_high);
    free(wi_high);
    
    free(Mass_Spline);
    free(Sigma_Spline);
    free(dSigmadm_Spline);
    free(second_derivs_sigma);
    free(second_derivs_dsigma);
    
}

int main(int argc, char ** argv){

    omp_set_num_threads(NUMCORES);

    char filename[500];
	FILE *F,*F1,*F2;
    fftwf_plan plan;
    
    // Various parameters to be used for the MCMC code
    float REDSHIFT, ION_EFF_FACTOR, MFP, TVIR_MIN, ALPHA, BETA, MFEEDBACK;
    int PERFORM_PS,short_completely_ionised;
    
    // Other parameters used in the code
    float growth_factor, init_growth_factor, xf, yf, zf;

    float mass, R, R_begin, pixel_mass, cell_length_factor;
	float ave_N_min_cell, M_MIN;
	int x,y,z, N_min_cell, LAST_FILTER_STEP;
	unsigned long long ion_ct;
	float f_coll_crit, pixel_volume, density_over_mean, erfc_num, erfc_denom, erfc_denom_cell, res_xH, Splined_Fcoll;
	float xHI_from_xrays, std_xrays;
	
    double global_xH, global_step_xH, ave_xHI_xrays, ave_den, ST_over_PS, mean_f_coll_st, mean_f_coll_ps, f_coll, ave_fcoll,TOT_NUM_NORM;
	const gsl_rng_type * T;
	gsl_rng * r;

    unsigned long long ct, aa, HII_i, HII_j, HII_k;
    
    int i,j,k,a, xi, yi, zi,ii;
    
    double ave_delta, new_ave_delta;
    
    float pixel_x_HI, pixel_deltax, H, dummy;
	int n_x, n_y, n_z, curr_Pop;
	double dvdx, ave, max_v_deriv;
	unsigned long long nonlin_ct, temp_ct;
	float nf, max, maxi, maxj, maxk, maxdvdx, min, mini, minj, mink, mindvdx;
	float k_x, k_y, k_z, k_mag, k_sq, k_mag_x, k_mag_y;
	float const_factor, T_rad, pixel_Ts_factor, curr_alphaX, curr_TvirX;
	double ave_Ts, min_Ts, max_Ts, temp, curr_zetaX;
    
    float collapsed_fraction,collapsed_fraction1,collapsed_fraction2,collapsed_fraction3;


    /**** Perform 'perturb_field.c' *****/
			
	/***************   BEGIN INITIALIZATION   **************************/
			 
    /*  Will remove all error checking etc. as it will be assumed that the MCMC will hand the correct parameter combinations etc. (also
        I will be removing most the file I/O therefore it won't be necessary to perform most of those checks */

// Redshift will remain constant through out the driver (handed to it by MCMC code)
	
	REDSHIFT = atof(argv[1]);
	ION_EFF_FACTOR = atof(argv[2]);
    MFP = atof(argv[3]);
//    TVIR_MIN = atof(argv[4]);
    TVIR_MIN = pow(10.,atof(argv[4]));
    
//  PERFORM PS flag decides whether or not to perform the PS computation (for example neutral fraction prior checking does not compute the PS)
    PERFORM_PS = atof(argv[5]);
    
    ALPHA = atof(argv[6]);
    BETA = atof(argv[7]);
    MFEEDBACK = pow(10.,atof(argv[8]));
    //    MFEEDBACK = pow(10.,atof(argv[12]));
    
	// perform a very rudimentary check to see if we are underresolved and not using the linear approx
	if ((BOX_LEN > DIM) && !EVOLVE_DENSITY_LINEARLY){
		printf("perturb_field.c: WARNING: Resolution is likely too low for accurate evolved density fields\n It Is recommended that you either increase the resolution (DIM/Box_LEN) or set the EVOLVE_DENSITY_LINEARLY flag to 1\n");
	}

	// initialize power spectrum 
    // set once, remove any later equivalent calls to these functions
	init_ps();
	growth_factor = dicke(REDSHIFT);
	   
    init_21cmMC_arrays();
    
    sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
    F = fopen(filename, "rb");
    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                if (fread((float *)deltax_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
                    printf("Read error occured while reading deltax box.\n");
                    return -1;
                }
            }
        }
    }
    fclose(F);

    switch(VELOCITY_COMPONENT){
        case 1:  sprintf(filename, "../Boxes/updated_vx_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
            break;
        case 3:  sprintf(filename, "../Boxes/updated_vz_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
            break;
        default: sprintf(filename, "../Boxes/updated_vy_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
    }
    F=fopen(filename, "rb");
    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                if (fread((float *)v + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
                    printf("Read error occured while reading velocity box.\n");
                    fclose(F);
                }
            }
        }
    }
    fclose(F);
    
    memcpy(deltax_unfiltered_original, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    
	/**** End of perform 'Perturb_Field.c' ******/
		
	/**** Perform 'find_HII_bubbles.c' ******/

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	pixel_volume = pow(BOX_LEN/(float)HII_DIM, 3);
	pixel_mass = RtoM(L_FACTOR*BOX_LEN/(float)HII_DIM); 
	f_coll_crit = 1/ION_EFF_FACTOR;
	cell_length_factor = L_FACTOR;

	// this parameter choice is sensitive to noise on the cell size, at least for the typical
	// cell sizes in RT simulations.  it probably doesn't matter for larger cell sizes.
	if (USE_HALO_FIELD && (FIND_BUBBLE_ALGORITHM==2)
		&& ((BOX_LEN/(float)HII_DIM) < 1)){ // fairly arbitrary length based on 2 runs i did
		cell_length_factor = 1; 
	}
	
	//set the minimum source mass
	if (TVIR_MIN > 0){ // use the virial temperature for Mmin
		if (TVIR_MIN < 9.99999e3) // neutral IGM
			M_MIN = TtoM(REDSHIFT, TVIR_MIN, 1.22);
		else // ionized IGM
			M_MIN = TtoM(REDSHIFT, TVIR_MIN, 0.6);
	}
	else if (TVIR_MIN < 0){ // use the mass
		M_MIN = ION_M_MIN;
	}
	// check for WDM
	if (P_CUTOFF && ( M_MIN < M_J_WDM())){
		printf( "The default Jeans mass of %e Msun is smaller than the scale supressed by the effective pressure of WDM.\n", M_MIN);
		M_MIN = M_J_WDM();
		printf( "Setting a new effective Jeans mass from WDM pressure supression of %e Msun\n", M_MIN);
	}
	
    if(MFEEDBACK > 0) {
        if(MFEEDBACK < 9.99999e3) { // neutral IGM
            MFEEDBACK = TtoM(REDSHIFT, MFEEDBACK, 1.22);
        }
        else {
            MFEEDBACK = TtoM(REDSHIFT, MFEEDBACK, 0.6);
        }
    }
    
    if(MFEEDBACK<M_MIN) {
        MFEEDBACK = M_MIN;
    }

	// initialize neutral fraction box
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
        xH[ct] = 1;
    }
	
	// lets check if we are going to bother with computing the inhmogeneous field at all...
    mean_f_coll_st = FgtrM_st_DPL(REDSHIFT,M_MIN,MFEEDBACK,ALPHA,BETA);
	if (mean_f_coll_st/f_coll_crit < HII_ROUND_ERR){ // way too small to ionize anything...
		printf( "The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n I will just declare everything to be neutral\n", mean_f_coll_st, f_coll_crit);
		
		if (USE_TS_IN_21CM){ // just use the x_e box
			sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
			if (!(F = fopen(filename, "rb"))){
				printf( "find_HII_bubbles: Unable to open x_e file at %s\nAborting...\n", filename);
				return -1;
			}
			for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
				if (fread(&xH[ct], sizeof(float), 1, F)!=1){
					printf( "find_HII_bubbles: Read error occured while reading xe box.\nAborting...\n");
                    return -1;
				}
				xH[ct] = 1-xH[ct]; // convert from x_e to xH
				global_xH += xH[ct];
			}
			fclose(F);
			global_xH /= (double)HII_TOT_NUM_PIXELS;
		}
		else{
			// find the neutral fraction
			global_xH = 1;
			        
            for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
                xH[ct] = 1;
            }
		}
		// print out the xH box
	}
	
	// see if we want to weigh by the number of neutral atoms from x-ray preheating
	if (USE_TS_IN_21CM){
		
		sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
		if (!(F = fopen(filename, "rb"))){
			printf( "find_HII_bubbles: Unable to open x_e file at %s\nAborting...\n", filename);
            return -1;
		}
		for (i=0; i<HII_DIM; i++){
			for (j=0; j<HII_DIM; j++){
				for (k=0; k<HII_DIM; k++){
					if (fread((float *)xe_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
						printf( "find_HII_bubbles: Read error occured while reading xe box.\nAborting...\n");
                        return -1;
					}
				}
			}
		}
		fclose(F);
	}
	
	if (USE_HALO_FIELD){

        // initialise the smoothed halo field
		for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS; ct++) {
            *((float *)M_coll_unfiltered + ct) = 0;
        }
		
		// read in the halo list
		sprintf(filename, "../Output_files/Halo_lists/updated_halos_z%06.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
		F = fopen(filename, "r");
		if (!F){
			printf( "find_HII_bubbles.c: Unable to open halo list file: %s\nAborting...\n", filename);
            return -1;
		}
		// now read in all halos above our threshold into the smoothed halo field
		fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
		while (!feof(F) && (mass>=M_MIN)){
			x = xf*HII_DIM;
			y = yf*HII_DIM;
			z = zf*HII_DIM;
			*((float *)M_coll_unfiltered + HII_R_FFT_INDEX(x, y, z)) += mass;
			fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
		}
		fclose(F);
	} // end of the USE_HALO_FIELD option
    
	if (USE_HALO_FIELD){
		plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)M_coll_unfiltered, (fftwf_complex *)M_coll_unfiltered, FFTW_ESTIMATE);
		fftwf_execute(plan);
	}
	if (USE_TS_IN_21CM){
		plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)xe_unfiltered, (fftwf_complex *)xe_unfiltered, FFTW_ESTIMATE);
		fftwf_execute(plan);    
	}
	plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
    
	// remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
	//  real space to k-space
	// Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
        
    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
        if (USE_HALO_FIELD){  M_coll_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);}
        if (USE_TS_IN_21CM){ xe_unfiltered[ct] /= (float)HII_TOT_NUM_PIXELS;}
        deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
    }
    
    int counter = 0;
    short_completely_ionised = 0;
	// loop through the filter radii (in Mpc)
	erfc_denom_cell=1; //dummy value
	R=fmin(MFP, L_FACTOR*BOX_LEN);
    R_begin = R;
	LAST_FILTER_STEP = 0;
    
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
    
    initialiseSplinedSigmaM(M_MIN,RtoM(MFP));
    
	while (!LAST_FILTER_STEP){//(R > (cell_length_factor*BOX_LEN/(HII_DIM+0.0))){
	       
        if ((R/DELTA_R_HII_FACTOR) <= (cell_length_factor*BOX_LEN/(float)HII_DIM)){
			LAST_FILTER_STEP = 1;
		}

		if (USE_HALO_FIELD){    memcpy(M_coll_filtered, M_coll_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);}
		if (USE_TS_IN_21CM){    memcpy(xe_filtered, xe_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);}
		memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
        
		if (USE_HALO_FIELD) { HII_filter(M_coll_filtered, HII_FILTER, R);}
		if (USE_TS_IN_21CM) { HII_filter(xe_filtered, HII_FILTER, R);}
		HII_filter(deltax_filtered, HII_FILTER, R);
        
		fftwf_execute(plan);
		// Check if this is the last filtering scale.  If so, we don't need deltax_unfiltered anymore.
		// We will re-read it to get the real-space field, which we will use to se the residual
		// neutral fraction
		ST_over_PS = 0;
		f_coll = 0;
		if (LAST_FILTER_STEP){
			
            memcpy(deltax_unfiltered, deltax_unfiltered_original, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
            
			if (USE_HALO_FIELD){
				for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS; ct++){    *((float *)M_coll_unfiltered + ct) = 0;  }
				sprintf(filename, "../Output_files/Halo_lists/updated_halos_z%06.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
				F = fopen(filename, "r");
				while (!feof(F) && (mass>=M_MIN)){
					x = xf*HII_DIM;
					y = yf*HII_DIM;
					z = zf*HII_DIM;
					*((float *)M_coll_unfiltered + HII_R_FFT_INDEX(x, y, z)) += mass;
					fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
				}
				fclose(F);
			} // end halo field option
			else{
                if(ALPHA==BETA) {
                    initialiseGL_Fcoll(NGLlow,NGLhigh,M_MIN,RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM));
                    initialiseFcoll_spline(REDSHIFT,M_MIN,RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM),RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM),MFEEDBACK,ALPHA);
                }
                else {
                    initialiseGL_FcollDblPl(NGLlow_DblPl,NGLhigh_DblPl,M_MIN,MFEEDBACK,RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM));
                    initialiseFcoll_spline_DblPl(REDSHIFT,M_MIN,RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM),RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM),ALPHA,BETA,MFEEDBACK);
                }
            				
                // renormalize the collapse fraction so that the mean matches ST,
                // since we are using the evolved (non-linear) density field
                TOT_NUM_NORM = 0.;
                f_coll = 0;
                for (x=0; x<HII_DIM; x++){
                    for (y=0; y<HII_DIM; y++){
                        for (z=0; z<HII_DIM; z++){
                            
                            if(*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)) < Deltac) {
                                        
                                if(ALPHA==BETA) {
                                    FcollSpline(*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)),&(Splined_Fcoll));
                                    Fcoll[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll;
                                }
                                else {
                                    FcollSplineDblPl(*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)),MFEEDBACK,RtoM(R),&(Splined_Fcoll));
                                    Fcoll[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll;
                                }
                                f_coll += Splined_Fcoll;
                            }
                            else {
                                f_coll += 1.;
                                TOT_NUM_NORM += 1.;
                            }
                        }
                    }
                }
                f_coll = (f_coll - TOT_NUM_NORM)/((double)HII_TOT_NUM_PIXELS - TOT_NUM_NORM);
                ST_over_PS = ((mean_f_coll_st*(double)HII_TOT_NUM_PIXELS - TOT_NUM_NORM)/((double)HII_TOT_NUM_PIXELS - TOT_NUM_NORM) )/f_coll;
            }
			
			// and the spin temperature option is set
			if (USE_TS_IN_21CM){
				sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
				F = fopen(filename, "rb");
				for (i=0; i<HII_DIM; i++){
					for (j=0; j<HII_DIM; j++){
						for (k=0; k<HII_DIM; k++){
							fread((float *)xe_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F);
						}
					}
				}
				fclose(F);
			} // end re-reading the x_e box
		} // end if last filter step conditional statement
        
		// not the last filter step, and we operating on the density field
		else if (!USE_HALO_FIELD){

            if(ALPHA==BETA) {
                initialiseGL_Fcoll(NGLlow,NGLhigh,M_MIN,RtoM(R));
                initialiseFcoll_spline(REDSHIFT,M_MIN,RtoM(R),RtoM(R),MFEEDBACK,ALPHA);
            }
            else {
                initialiseGL_FcollDblPl(NGLlow_DblPl,NGLhigh_DblPl,M_MIN,MFEEDBACK,RtoM(R));
                initialiseFcoll_spline_DblPl(REDSHIFT,M_MIN,RtoM(R),RtoM(R),ALPHA,BETA,MFEEDBACK);
            }

            // renormalize the collapse fraction so that the mean matches ST,
            // since we are using the evolved (non-linear) density field
            f_coll = 0.;
            TOT_NUM_NORM = 0.;
            for (x=0; x<HII_DIM; x++){
                for (y=0; y<HII_DIM; y++){
                    for (z=0; z<HII_DIM; z++){
                            
                        if(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) < Deltac) {
                            
                            if(ALPHA==BETA) {
                                FcollSpline(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)),&(Splined_Fcoll));
                                Fcoll[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll;
                            }
                            else {
                                FcollSplineDblPl(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)),MFEEDBACK,RtoM(R),&(Splined_Fcoll));
                                Fcoll[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll;
                            }
                            f_coll += Splined_Fcoll;
                        }
                        else {
                            f_coll += 1.0;
                            TOT_NUM_NORM += 1.0;
                        }
                    }
                }
            }
            f_coll = (f_coll - TOT_NUM_NORM)/((double)HII_TOT_NUM_PIXELS - TOT_NUM_NORM);
            ST_over_PS = ((mean_f_coll_st*(double)HII_TOT_NUM_PIXELS - TOT_NUM_NORM)/((double)HII_TOT_NUM_PIXELS - TOT_NUM_NORM) )/f_coll;
        }
		/************  MAIN LOOP THROUGH THE BOX **************/
            
        // now lets scroll through the filtered box
        ave_xHI_xrays = ave_den = ave_fcoll = std_xrays = 0;
            
        xHI_from_xrays = 1;
        
        ion_ct=0;
        for (x=0; x<HII_DIM; x++){
            for (y=0; y<HII_DIM; y++){
                for (z=0; z<HII_DIM; z++){
                    if (LAST_FILTER_STEP) {
                            
                        if(*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)) <= -1.) {
                                
                            if(ALPHA==BETA) {
                                FcollSpline((FRACT_FLOAT_ERR - 1.),&(Splined_Fcoll));
                            }
                            else {
                                FcollSplineDblPl((FRACT_FLOAT_ERR - 1.),MFEEDBACK,RtoM(R),&(Splined_Fcoll));
                            }
                            f_coll = ST_over_PS * Splined_Fcoll;
                                
                        }
                        else {
                            
                            if(*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)) < Deltac) {
                                    
                                f_coll = ST_over_PS * Fcoll[HII_R_FFT_INDEX(x,y,z)];
                            }
                            else {
                                f_coll = 1.;
                            }
                        }
                    }
                    else {
                        
                        if(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) <= -1.) {
                            
                            if(ALPHA==BETA) {
                                FcollSpline((FRACT_FLOAT_ERR - 1.),&(Splined_Fcoll));
                            }
                            else {
                                FcollSplineDblPl((FRACT_FLOAT_ERR - 1.),MFEEDBACK,RtoM(R),&(Splined_Fcoll));
                            }
                            f_coll = ST_over_PS * Splined_Fcoll;
                        }
                        else {
                                
                            if(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) < Deltac) {
                                
                                f_coll = ST_over_PS * Fcoll[HII_R_FFT_INDEX(x,y,z)];
                            }
                            else {
                                f_coll = 1.;
                            }
                        }
                    }
                        
                    // Removed two pieces of code here regarding the use of the Halo field and the calculation of the spin temperature. Neither condition is ever set for the determination of the 21cm PS
                        
                    // check if ionized!
                    if (f_coll > f_coll_crit){ // ionized
						
                        if (FIND_BUBBLE_ALGORITHM == 2) // center method
                            xH[HII_R_INDEX(x, y, z)] = 0;
						
                        else if (FIND_BUBBLE_ALGORITHM == 1) // sphere method
                            update_in_sphere(xH, HII_DIM, R/BOX_LEN, x/(HII_DIM+0.0), y/(HII_DIM+0.0), z/(HII_DIM+0.0));
						
                        else{
                            printf( "Incorrect choice of find bubble algorithm: %i\nAborting...", FIND_BUBBLE_ALGORITHM);
                            fflush(NULL);
                            z=HII_DIM;y=HII_DIM,x=HII_DIM;R=0;
                        }
                    }
					
                    // check if this is the last filtering step.
                    // if so, assign partial ionizations to those cells which aren't fully ionized
                    else if (LAST_FILTER_STEP && (xH[HII_R_INDEX(x, y, z)] > TINY)){

                        if(*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)) <= -1.) {
                                
                            if(ALPHA==BETA) {
                                FcollSpline((FRACT_FLOAT_ERR - 1.),&(Splined_Fcoll));
                            }
                            else {
                                FcollSplineDblPl((FRACT_FLOAT_ERR - 1.),MFEEDBACK,RtoM(R),&(Splined_Fcoll));
                            }
                            f_coll = ST_over_PS * Splined_Fcoll;
                                
//                              if (f_coll>1) f_coll=1;
                                
                            ave_N_min_cell = f_coll * pixel_mass*(FRACT_FLOAT_ERR) / M_MIN; // ave # of M_MIN halos in cell
                                
                            if (ave_N_min_cell < N_POISSON){
                                N_min_cell = (int) gsl_ran_poisson(r, ave_N_min_cell);
                                f_coll = N_min_cell * M_MIN / (pixel_mass*FRACT_FLOAT_ERR);
                            }
                        }
                            
                        else {
                            if(*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)) < Deltac) {
                                    
                                f_coll = ST_over_PS * Fcoll[HII_R_FFT_INDEX(x,y,z)];
                                    
                            }
                            else {
                                f_coll = 1.;
                            }
                                
//                            if (f_coll>1) f_coll=1;
                                
                            ave_N_min_cell = f_coll * pixel_mass*(1. + *((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z))) / M_MIN; // ave # of M_MIN halos in cell
                                
                            if (ave_N_min_cell < N_POISSON){
                                N_min_cell = (int) gsl_ran_poisson(r, ave_N_min_cell);
                                f_coll = N_min_cell * M_MIN / (pixel_mass*(1. + *((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z))));
                            }
                        }
                        //                                if (f_coll>1) f_coll=1;
                        res_xH = xHI_from_xrays - f_coll * ION_EFF_FACTOR;
                        // and make sure fraction doesn't blow up for underdense pixels
                        if (res_xH < 0)
                            res_xH = 0;
                        else if (res_xH > 1)
                            res_xH = 1;
                            
                        xH[HII_R_INDEX(x, y, z)] = res_xH;
                    } // end partial ionizations at last filtering step
					
                    // Debugging below
                    //	  ave_fcoll += f_coll;
                    // ave_xHI_xrays += xHI_from_xrays;
                    // std_xrays += pow(xHI_from_xrays - 0.205, 2);
                    // ave_den += density_over_mean;
                    // if (xH[HII_R_INDEX(x, y, z)] > 0.5)
                    // ion_ct++;
                } // k
            } // j
        } // i

		// Debugging
		//    printf( "Mean x-ray neutral fraction is %f, density is %f, %f of all pixels are neutral, <fcoll>=%f\n\n",
		// ave_xHI_xrays/(double)HII_TOT_NUM_PIXELS, ave_den/(double)HII_TOT_NUM_PIXELS, ion_ct/(double)HII_TOT_NUM_PIXELS, ave_fcoll/(double)HII_TOT_NUM_PIXELS);
		// fprintf(LOG, "end of main loop scroll, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
		// fflush(LOG);
        
//        printf("Neutral fraction each step\n");
        global_step_xH = 0;
        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
            global_step_xH += xH[ct];
        }
        global_step_xH /= (float)HII_TOT_NUM_PIXELS;
//        printf("Neutral fraction = %e\n",global_step_xH);
        
        if(global_step_xH==0.0) {
//            printf("Just end it now! R = %e R_begin = %e\n",R,R_begin);
            short_completely_ionised = 1;
            break;
        }
        
        R /= DELTA_R_HII_FACTOR;
	}
    
	// find the neutral fraction
	global_xH = 0;
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
        global_xH += xH[ct];
    }
	global_xH /= (float)HII_TOT_NUM_PIXELS;
//    printf("Neutral fraction = %e\n",global_xH);
    
	// deallocate
	gsl_rng_free (r);
    
	/**** End of perform 'find_HII_bubbles.c' *****/
	
	/**** Perform 'delta_T.c' ******/
	
    nf = global_xH;
  
    sprintf(filename, "NeutralFraction_%s_%s_%s_alpha%s_beta%s_lnTFeed%s.txt",argv[2],argv[3],argv[4],argv[6],argv[7],argv[8]);
//    sprintf(filename, "NeutralFraction.txt");
//    printf("%s\n",filename);
    F=fopen(filename, "wt");
    fprintf(F, "%lf\n",nf);
    fclose(F);
  
    
	/************  BEGIN INITIALIZATION ****************************/
    
    max = -1e3;
	min = 1e3;
	ave = 0;
	nonlin_ct=0;
	   
    if(PERFORM_PS==1) {
        
        T_rad = T_cmb*(1+REDSHIFT);
        H = hubble(REDSHIFT);
        const_factor = 27 * (OMb*hlittle*hlittle/0.023) *
        sqrt( (0.15/OMm/hlittle/hlittle) * (1+REDSHIFT)/10.0 );

        memcpy(deltax, deltax_unfiltered_original, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
		
        /************  END INITIALIZATION ****************************/
    
        // ok, lets fill the delta_T box; which will be the same size as the bubble box
        ave_Ts = max_Ts = 0;
        min_Ts = 1e5;
        temp=0;
//        temp_ct=0;
      
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
				
                    pixel_deltax = deltax[HII_R_FFT_INDEX(i,j,k)];
                    pixel_x_HI = xH[HII_R_INDEX(i,j,k)];
				
                    if (pixel_x_HI > TINY){
                        temp = pixel_deltax;
//                           temp_ct++;
                    }
				
                    delta_T[HII_R_INDEX(i,j,k)] = const_factor*pixel_x_HI*(1+pixel_deltax);
								
                    if (max < delta_T[HII_R_INDEX(i,j,k)]){ max = delta_T[HII_R_INDEX(i,j,k)];}
                    if (min > delta_T[HII_R_INDEX(i,j,k)]){ min = delta_T[HII_R_INDEX(i,j,k)];}
                    ave += delta_T[HII_R_INDEX(i,j,k)];
                }
            }
        }
        //	printf("%e\n", temp/(double)temp_ct);
        ave /= (float)HII_TOT_NUM_PIXELS;
//        printf( "Without velocities, max is %e, min is %e, ave is %e\n", max, min, ave);
	
        if (USE_TS_IN_21CM){
            ave_Ts /= (double) HII_TOT_NUM_PIXELS;
            printf( "Ts, min = %e, max = %e, ave = %e\n", min_Ts, max_Ts, ave_Ts);
            printf( "corresponding to (1-trad/Ts of), min = %e, max = %e, ave = %e\n", 1-T_rad/min_Ts, 1-T_rad/max_Ts, 1-T_rad/ave_Ts);
        }
		
        // now write out the delta_T box
        if (T_USE_VELOCITIES){
            max = -1;
            min = 1e3;
            ave = 0;
            
            // let's take the derivative in k-space
            plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)v, (fftwf_complex *)v, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);

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
//                          nonlin_ct++;
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
        }
        
        //        printf("average = %e\n",ave);
        
        /******  PRINT OUT THE POWERSPECTRUM  *********/
        
        for (ct=0; ct<NUM_BINS; ct++){
            p_box[ct] = k_ave[ct] = 0;
            in_bin_ct[ct] = 0;
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
//        fftwf_init_threads();
//        fftwf_plan_with_nthreads(NUMCORES);
        
        // transform to k-space
        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T, (fftwf_complex *)deldel_T, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

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
        
        // now lets print out the k bins
//	sprintf(filename, "delTps_estimate.txt");
        sprintf(filename, "delTps_estimate_%s_%s_%s_alpha%s_beta%s_lnTFeed%s.txt",argv[2],argv[3],argv[4],argv[6],argv[7],argv[8]);
//        printf("%s\n",filename);
        F=fopen(filename, "wt");
        for (ct=1; ct<NUM_BINS; ct++){
            if (in_bin_ct[ct]>0)
			fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
        }
        fclose(F);
    }
	/****** END POWER SPECTRUM STUFF   ************/
	
	/**** End of perform 'delta_T.c' *****/
	
//    destroy_21cmMC_arrays();

    return 0;
}

/*
 FUNCTIONS
 
 */

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

float GaussLegengreQuad_Fcoll(int type, int n, float z, float M2, float MFeedback, float alpha, float delta1, float delta2)
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
            if(type==1) {
                x = xi_low[i];
                integrand += wi_low[i]*FgtrConditionallnM(x,parameters_gsl_int);
            }
            if(type==2) {
                x = xi_high[i];
                integrand += wi_high[i]*FgtrConditionallnM(x,parameters_gsl_int);
            }
            if(type==3) {
                x = xi_low_dblpl_alpha[i];
                integrand += wi_low_dblpl_alpha[i]*FgtrConditionallnM(x,parameters_gsl_int);
            }
            if(type==4) {
                x = xi_high_dblpl_alpha[i];
                integrand += wi_high_dblpl_alpha[i]*FgtrConditionallnM(x,parameters_gsl_int);
            }
            if(type==5) {
                x = xi_low_dblpl_beta[i];
                integrand += wi_low_dblpl_beta[i]*FgtrConditionallnM(x,parameters_gsl_int);
            }
            if(type==6) {
                x = xi_high_dblpl_beta[i];
                integrand += wi_high_dblpl_beta[i]*FgtrConditionallnM(x,parameters_gsl_int);
            }
        }
        return integrand;
    }
}

/*
 FUNCTION FgtrM_st(z, M)
 Computes the fraction of mass contained in haloes with mass > M at redshift z
 Uses Sheth-Torman correction
 */
double dFdlnM_st_DPL (double lnM, void *params){

    struct parameters_gsl_ST_int_ vals = *(struct parameters_gsl_ST_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    float MFeedback = vals.M_Feed;
    float alpha = vals.alpha_pl;
    
    return dNdM_st(z, M) * M * M * pow((M/MFeedback),alpha);
}
double FgtrM_st_DPL(double z, double Mmin, double MFeedback, double alpha_pl, double beta_pl){

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
    
    F.function = &dFdlnM_st_DPL;
    F.params = &parameters_gsl_ST_lower;
    lower_limit = log(Mmin);
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
    
    return (result_lower + result_upper) / (OMm*RHOcrit);
}


