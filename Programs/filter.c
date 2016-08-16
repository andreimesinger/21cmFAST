#ifndef _FILTER_
#define _FILTER_

#include "../Parameter_files/INIT_PARAMS.H"

/*
  Function FILTER filters the k-space box, <box>, using filter type
  <filter_type> on a characteristic comoving scale <R> (in Mpc), where:
  0 = top-hat real space filter
  1 = top-hat k-space filter
  2 = gaussian

  Relavant box parameters are taken from INIT_PARAMS.H

  The function returns the filtered k field, <box>.
*/

void filter(fftwf_complex *box, int filter_type, float R){
  int n_x, n_z, n_y;
  float k_x, k_y, k_z, k_mag, kR;

  // loop through k-box
#pragma omp parallel shared(box, filter_type, R) private(k_x, k_y, k_z, k_mag, kR, n_x, n_z, n_y)
{

#pragma omp for 
  for (n_x=0; n_x<DIM; n_x++){
    if (n_x>MIDDLE) {k_x =(n_x-DIM) * DELTA_K;}
    else {k_x = n_x * DELTA_K;}

    for (n_y=0; n_y<DIM; n_y++){
      if (n_y>MIDDLE) {k_y =(n_y-DIM) * DELTA_K;}
      else {k_y = n_y * DELTA_K;}

      for (n_z=0; n_z<=MIDDLE; n_z++){ 
	k_z = n_z * DELTA_K;
	
	k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

	kR = k_mag*R; // real space top-hat
	if (filter_type == 0){ // real space top-hat
	  if (kR > 1e-4){
	    box[C_INDEX(n_x, n_y, n_z)] *= 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));
	  }
	}
	else if (filter_type == 1){ // k-space top hat
	  kR *= 0.413566994; // equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
	  if (kR > 1){
	    box[C_INDEX(n_x, n_y, n_z)] = 0;
	  }
	}
	else if (filter_type == 2){ // gaussian
	  kR *= 0.643; // equates integrated volume to the real space top-hat
	  box[C_INDEX(n_x, n_y, n_z)] *= pow(E, -kR*kR/2.0);
	}
	else{
	  if ( (n_x==0) && (n_y==0) && (n_z==0) )
	    fprintf(stderr, "filter.c: Warning, filter type %i is undefined\nBox is unfiltered\n", filter_type);
	}
      }
    }
  } // end looping through k box
  
}
  return;
}

#endif
