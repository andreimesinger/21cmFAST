/*
  These include fundamental parameters of the realization.
  The cosmology is defined in COSMOLOGY.H
*/

#ifndef _INITS_PARAMETERS_H_
#define _INITS_PARAMETERS_H_

#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <string.h>
//#include <pthread.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include "COSMOLOGY.H"


/******** BEGIN USER CHANGABLE DEFINITIONS   **********/
#define RANDOM_SEED (long) (1) // seed for the random number generator
#define BOX_LEN (float) 300 // in Mpc
#define DIM (int) 800 // number of cells for the high-res box (sampling ICs) along a principal axis
// to avoid sampling issues, DIM should be at least 3 or 4 times HII_DIM, and an integer multiple)
#define HII_DIM (int) 200 // number of cells for the low-res box

/****** New in v1.1. ****** Threading parameters  ***/
#define NUMCORES (int) 2 // # of cores you wish to allocate (must be shared mem)
#define RAM (float) 8 // physical memory in GB available
/******** END USER CHANGABLE DEFINITIONS   **********/

#include "ANAL_PARAMS.H"
#include "HEAT_PARAMS.H"
#include "../Cosmo_c_files/misc.c"
#include "../Cosmo_c_files/cosmo_progs.c"
#include "../Cosmo_c_files/ps.c"
#include "../Cosmo_c_files/recombinations.c"

#define MIDDLE (DIM/2)
#define D (unsigned long long)DIM // the unsigned long long dimension
#define MID ((unsigned long long)MIDDLE)
#define VOLUME (BOX_LEN*BOX_LEN*BOX_LEN) // in Mpc^3
#define DELTA_K (TWOPI/BOX_LEN)
#define TOT_NUM_PIXELS ((unsigned long long)(D*D*D)) // no padding
#define TOT_FFT_NUM_PIXELS ((unsigned long long)(D*D*2llu*(MID+1llu)))
#define KSPACE_NUM_PIXELS ((unsigned long long)(D*D*(MID+1llu)))


/* Define some useful macros */
#define C_INDEX(x,y,z)((unsigned long long)((z)+(MID+1llu)*((y)+D*(x))))// for 3D complex array
#define R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(MID+1llu)*((y)+D*(x)))) // for 3D real array with the FFT padding
#define R_INDEX(x,y,z)((unsigned long long)((z)+D*((y)+D*(x)))) // for 3D real array with no padding

#endif
