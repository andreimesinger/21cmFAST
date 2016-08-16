/************************************************************************
Run of the mill FFT functions.  
Written by: Matt McQuinn
Probably not so important to look in here
***********************************************************************/
#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"


void doFFT(float *R, fftwf_complex *Rf, long Nf, int dim);
void doInverseFFT(fftwf_complex *Rf, float *R, long Nf, int dim);
void initialize(float **K,float **Pk, int **count, int flag);
int knum(int i);
void makePk(float *K, float *Pk, int *Count, fftwf_complex *Rf);


/*Called to take Fourier transform of field*/
void doFFT(float *R, fftwf_complex *Rf, long Nf, int dim)
{
    int i;
    float dvol = pow(BOX_LEN/HII_DIM,dim);
    fftwf_plan p;
    if(dim == 2)
	p  = fftwf_plan_dft_r2c_2d(HII_DIM, HII_DIM, (float *)R, (fftwf_complex *)Rf, FFTW_ESTIMATE);
    else
	p  = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)R, (fftwf_complex *)Rf, FFTW_ESTIMATE);
    fftwf_execute(p);

    for(i = 0; i<Nf; i++)
    {
	Rf[i] *= dvol;  //d^3x differential
    }

    fftwf_destroy_plan(p);
}

/* Calculates inverse fourier transform of real field R*/
void doInverseFFT(fftwf_complex *Rf, float *R, long N, int dim)
{
    int i;
    fftwf_plan p;
    float dvol = pow(BOX_LEN,dim);

    if(dim == 2)
	p  = fftwf_plan_dft_c2r_2d(HII_DIM, HII_DIM, Rf, R, FFTW_ESTIMATE);
    else 
	p  = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, Rf, R, FFTW_ESTIMATE);
    fftwf_execute(p);
    for(i = 0; i<N; i++)
    {
	R[i] /= dvol;
    }

    fftwf_destroy_plan(p);
}


int knum(int i)
{
    if(i < HII_DIM/2 +1)
	return i;
    else 
	return -(HII_DIM - i);
} 



//allocates and initializes arrays for flag == 0, otherwise frees arrays
void initialize(float **K,float **Pk, int **Count, int flag)
{
    int i;

    if(flag == 0)
    {
	*K = (float *) malloc(sizeof(float)*(HII_DIM*sqrt(2)/2.+1));
	*Pk = (float *) malloc(sizeof(float)*(HII_DIM*sqrt(2)/2.+1));
	*Count = (int *) malloc(sizeof(int)*(HII_DIM*sqrt(2)/2. +1));  

//to be safe
	for(i =  0; i<= HII_DIM*sqrt(2)/2.; i++)
	{
	    Pk[0][i] = 0.;
	    Count[0][i] = 0;
	    K[0][i] = 0.;
	}
    }else if(flag==1){
	free(*Pk); 
	free(*Count);
	free(*K);
    }
    else{
	for(i =  0; i<= HII_DIM*sqrt(2)/2.; i++)
	{
	    Pk[0][i] = 0.;
	    Count[0][i] = 0;
	    K[0][i] = 0.;
	}
    }
}

/*******************************************************************
Makes power spectrum
*******************************************************************/
void makePk(float *K, float *Pk, int *Count, fftwf_complex *Rf)
{
    int i, k, l, a, n;
    float m;
    float dk = 2.*M_PI/BOX_LEN;
    float volume = pow(BOX_LEN, 2);

    for(n = 0; n< HII_DIM*sqrt(2)/2.; n++)
    {
	Count[n] = 0.;
	K[n] = 0.;
    }
  
    for(i =  0; i< HII_DIM; i++)
	    for(k =  0; k< HII_DIM/2 + 1; k++)
	    {
		a = 2;
		if(k == 0 || k == HII_DIM/2 + 1)
		    a = 1;
		l = i*(HII_DIM/2 +1) + k;
		m = sqrt(knum(i)*knum(i) + k*k);
		n = (int) (m + .5);
		K[n] += a*dk*m;
      
		Pk[n] += a*( crealf(Rf[l])*crealf(Rf[l]) + cimagf(Rf[l])*cimagf(Rf[l]) );

		Count[n]+=a;
	    }

  
    for(n = 0; n< HII_DIM*sqrt(2)/2.; n++)
    {
	K[n] /= Count[n];
	Pk[n] /= (Count[n]*volume);
    }
}
