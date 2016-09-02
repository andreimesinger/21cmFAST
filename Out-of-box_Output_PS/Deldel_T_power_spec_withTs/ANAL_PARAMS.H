#ifndef _ANAL_PARAM_H_
#define _ANAL_PARAM_H_


/******** BEGIN USER CHANGABLE DEFINITIONS   **********/
/* 
   Minimum halo mass counted as contributing ionizing photons in the bubble code.
   
*/
#define ION_M_MIN (double) (-1)
/* OR define it as a redshift-dependent halo virial temperature.. 
 NOTE: SET THE ONE YOU AREN'T USING TO A NEGATIVE VALUE.
 NOTE: If you are using the halo finder and studying redshift evolution over a 
       small redshift interval, it is probably best to use the ION_M_MIN option,
       since the discretization of the halo mass by the halo finder might lead
       to redshift discontinuities.
*/
#define ION_Tvir_MIN (double) (3e4)


/*
  HII efficiency factor, commonly refered to as \zeta: i.e. number of ionizing photons escaping into the IGM
  per baryon.

  NEW in v1.2:  The ionizing efficiency is allowed to have a power law dependence on halo mass:
  \zeta = HII_EFF_FACTOR * (M_halo / ION_M_MIN)^EFF_FACTOR_PL_INDEX
  see Nion_note.pdf in the root folder.
*/
#define HII_EFF_FACTOR (float) (20)


/*
  NEW IN V1.2: power law index for the ionizing efficiency scaling with halo mass:
  \zeta = HII_EFF_FACTOR * (M_halo / ION_M_MIN)^EFF_FACTOR_PL_INDEX

  The traditional constant mass to ionizing light ratio (as is assumed in previous versions,
  is obtained by setting EFF_FACTOR_PL_INDEX to zero.
 */
#define EFF_FACTOR_PL_INDEX (float) (0)


/*
  Maximum horizon for ionizing photons, in comoving Mpc.
  The bubble filtering procedure uses the minimum
  value of L_FACTOR*BOX_LEN and R_BUBBLE_MAX as the maximum allowed
  radius of a single, cohesive ionized region. (in practice, the ionization
  field is extreemly insensitive to this choice)
*/
#define R_BUBBLE_MAX (float) (20)


/*
  EVOLVE_DENSITY_LINEARLY = 1, evolve the density field with linear theory.
                If choosing this option, make sure that your cell size is
		in the linear regime at the redshift of interest

  EVOLVE_DENSITY_LINEARLY = 0, evolve the density field with 1st order perturbation theory.
                If choosing this option, make sure that you resolve small 
		enough scales, roughly we find BOX_LEN/DIM should be < 1Mpc
 */
#define EVOLVE_DENSITY_LINEARLY (int) (0)


/*
  If set to 1, the ZA density field is additionally smoothed (asside from the implicit
  boxcar smoothing performed when re-binning the ICs from DIM to HII_DIM) with a Gaussian
  filter of width 0.2*BOX_LEN/HII_DIM.  The implicit boxcar smoothing in perturb_field.c
  bins the density field on scale DIM/HII_DIM, similar to what Lagrangian codes do 
  when constructing Eulerian grids. In other words, the density field, \delta,
  is quantized into (DIM/HII_DIM)^3 values. If your usage requires smooth density fields,
  it is recommended to set SMOOTH_EVOLVED_FIELD to 1.  This also decreases the shot noise
  present in all grid based codes, though it overcompensates by an effective loss in resolution.
  New in v1.1
 */
#define SMOOTH_EVOLVED_DENSITY_FIELD (int) (1)

/*
  Allows one to set a flag allowing find_HII_bubbles to skip constructing the 
  ionization field if it is estimated that the mean neutral fraction, <xH>, is
  within HII_ROUND_ERR of 1. In other words, if <xH> > 1-HII_ROUND_ERR, 
  then find_HII_bubbles just prints a homogeneous xHI field
  of  1's.

  This is a new option in v1.1. Previous versions had a hardcoded value of 1e-15.
*/
#define HII_ROUND_ERR (float) (1e-3)


/*
  Choice of:
  1 - Mesinger & Furlanetto 2007 method of overlaping spheres:
  paint an ionized sphere with radius R, centered on pixel
  where R is filter radius
  This method, while somewhat more accurate, is slower than (2) especially
  in mostly ionized unverses, so only use for lower resolution boxes (HII_DIM<~400)
  2 - Center pixel only method (Zahn et al. 2007). this is faster.
*/
#define FIND_BUBBLE_ALGORITHM (int) (2)

/* 
   Minimum radius of an HII region in cMpc.  One can set this to 0, but should be careful with
   shot noise if the find_HII_bubble algorithm is run on a fine, non-linear density grid.
 */
#define R_BUBBLE_MIN (float) (1)


/*
  1 = Use the filtered and adjusted halo field to construct the ionization and 21-cm fields
  0 = Use the mean collapse fraction (eq. 2 in Barkana & Loeb 2008) to construct ...

  Note: option (1) is more accurate as the halo field captures the Poisson and clustering signature
  which dominates the collapse fraction field on small scales (i.e. where the mean halo number 
  is small).  However, option (0) is faster as it bypasses the halo finder, and it becomes
  increasingly accurate as the typical bubble size increases.
*/
#define USE_HALO_FIELD (int) (0)


/*
  If not using the halo field to generate HII regions, we provide the option of
  including Poisson scatter in the number of sources obtained through the conditional
  collapse fraction (which only gives the *mean* collapse fraction on a particular
  scale.  If the predicted mean collapse fraction is < N_POISSON * M_MIN,
  then Poisson scatter is added to mimic discrete halos on the subgrid scale
  (see Zahn+ 2010).

  NOTE: If you are interested in snapshots of the same realization at several redshifts,
  it is recommended to turn off this feature, as halos can stocastically
  "pop in and out of" existance from one redshift to the next...
*/
#define N_POISSON (int) (-1)


/*
  Parameter choice of whether to use velocity corrections in 21-cm fields
  1=use velocities in delta_T; 0=do not use velocities

  NOTE: The approximation used to include peculiar velocity effects works 
  only in the linear regime, so be careful using this (see Mesinger+2010)
*/
#define T_USE_VELOCITIES (int) (1)

/*                                                                                                    
  Maximum velocity gradient along the line of sight in units of the hubble parameter at z.
  This is only used in computing the 21cm fields.
  Note, setting this too high can add spurious 21cm power in the early stages, due to the
  1-e^-tau ~ tau approximation (see my 21cm intro paper and mao+2011).  However, this is still
  a good approximation at the <~10% level.  Future builds will include the mao+ model for
  redshift space distortions.
*/
#define MAX_DVDR (float) (0.2)

/*
  Component of the velocity to be used in 21-cm temperature maps
  1 = x
  2 = y
  3 = z
*/
#define VELOCITY_COMPONENT (int) (3)


/*
  0 = plot 21cm temperature power spectrum in non-dimensional units
  1 = plot 21cm...  in mK^2
*/
#define DIMENSIONAL_T_POWER_SPEC (int) (1)

#define DELTA_R_FACTOR (float) (1.1) // factor by which to scroll through filter radius for halos

#define DELTA_R_HII_FACTOR (float) (1.1) // factor by which to scroll through filter radius for bubbles

/*
Factor of the halo's radius, R, so that the effective radius
is R_eff = R_OVERLAP_FACTOR * R.  Halos whose centers are less than
R_eff away from another halo are not allowed.
R_OVERLAP_FACTOR = 1 is fully disjoint 
R_OVERLAP_FACTOR = 0 means that centers are allowed to luy on the edges of
neighboring halos
*/
#define R_OVERLAP_FACTOR (float) (1.0)

/*
0 = delta_crit is constant 1.68
1 = delta_crit is the sheth tormen ellipsoidal collapse
    correction to delta_crit (see ps.c)
*/
#define DELTA_CRIT_MODE (int) (1)


/*
Filter for the density field used to generate the halo field with EPS
0 = use real space top hat filter to smooth density field
1 = use k-space top hat filter
2 = use gaussian filter
*/
#define HALO_FILTER (int) (0)

/*
Filter for the Halo or density field used to generate ionization field
0 = use real space top hat filter
1 = use k-space top hat filter
2 = use gaussian filter
*/
#define HII_FILTER (int) (1)

/*
  0 = don't run the optimization code for zscroll (only useful at low-z)
  1 = do run ...
*/
#define OPTIMIZE (int) (0)
/*
  Minimum mass for which the optimization algorithm will be used
*/
#define OPTIMIZE_MIN_MASS (float) (1e11)

#define SIZE_RANDOM_SEED (-23456789) // seed for the size dist random number generator
#define LOS_RANDOM_SEED (-123456789) // seed for the extract LOS random number generator

#define INITIAL_REDSHIFT (float) 300 // used to perturb field

/************  END USER CHANGABLE DEFINITIONS  **************/

#define L_FACTOR (float) (0.620350491) // factor relating cube length to filter radius = (4PI/3)^(-1/3)

#define HII_D (unsigned long long) (HII_DIM)
#define HII_MIDDLE (HII_DIM/2)
#define HII_MID ((unsigned long long)HII_MIDDLE)

#define HII_TOT_NUM_PIXELS (unsigned long long)(HII_D*HII_D*HII_D)
#define HII_TOT_FFT_NUM_PIXELS ((unsigned long long)(HII_D*HII_D*2llu*(HII_MID+1llu)))
#define HII_KSPACE_NUM_PIXELS ((unsigned long long)(HII_D*HII_D*(HII_MID+1llu)))

/* INDEXING MACROS */
#define HII_C_INDEX(x,y,z)((unsigned long long)((z)+(HII_MID+1llu)*((y)+HII_D*(x))))// for 3D complex array
#define HII_R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(HII_MID+1llu)*((y)+HII_D*(x)))) // for 3D real array with the FFT padding
#define HII_R_INDEX(x,y,z)((unsigned long long)((z)+HII_D*((y)+HII_D*(x)))) // for 3D real array with no padding

#endif
