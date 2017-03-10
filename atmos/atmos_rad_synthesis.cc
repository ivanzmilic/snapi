#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "const.h"
#include "types.h"
#include "io.h"
#include "mathtools.h"
#include "mem.h"

#include "atmol/atmol.h"

#include "atmos.h"

// This source file contains following functions/methods:
// int int op_em_vector(fp_t,fp_t,fp_t*,int,fp_t***,fp_t**); - op and em for the 1D atmosphere, all wavelengths at once

int atmosphere::op_em_vector(fp_t *** Vlos, fp_t **** B, fp_t theta,fp_t phi,fp_t* lambda,int nlambda,fp_t****** op_vector,fp_t***** em_vector){
  // This function computes opacity and emissivity for atmosphere, all wavelengths at once. 
  // Returns 0 for success, 1 for error. 
  // Arguments:
  // fp_t *** Vlos - l.o.s. velocity
  // fp_t **** B - projected magnetic field in the reference frame of the ray
  // fp_t theta, phi - emergent angle for radiation
  // fp_t * lambda, int nlambda - wavelength array and number of wavelengths 
  // fp_t *** op_vector, NDx4x4 array of opacities, fp_t ** em_vector, ND x 4 array of emissivities

  // First we account for electron opacity and emissivity. Thomson scattering is very weakly dependent on wavelength
  // so we will just compute it for the middle wavelenght here:
  fp_t lambda_m = (lambda[nlambda] + lambda[1]) * 0.5;

  //memset input arrays to zero:
  memset(op_vector[1][x1l][x2l][x3l][1]+1,0,nlambda*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));
  memset(em_vector[1][x1l][x2l][x3l]+1,0,nlambda*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));

  for (int x1i=x1l;x1i<=x1h;++x1i)
	for (int x2i=x2l;x2i<=x2h;++x2i)
	  for (int x3i=x3l;x3i<=x3h;++x3i){
	    fp_t op = Ne[x1i][x2i][x3i] * 6.65E-25;
		fp_t em = op * Planck_f(lambda_m, T[x1i][x2i][x3i]);
		for (int l=1;l<=nlambda;++l){
		  em_vector[l][x1i][x2i][x3i][1] = em;
		  op_vector[l][x1i][x2i][x3i][1][1] = op;
		}
  }

  // Then add all the contributors from opacity and emissivity from atoms and molecules
  
  for (int a=0;a<=natm;++a){
  	atml[a]->op_em_vector(T,Ne,Vlos,Vt,B,theta,phi,lambda,nlambda,op_vector,em_vector);
  }

  // Before exiting, reorder the opacity, so absorption matrix is properly set-up:
  for (int l=1;l<=nlambda;++l)
  	for (int x1i=x1l;x1i<=x1h;++x1i)
  	  for (int x2i=x2l;x2i<=x2h;++x2i)
  	  	for (int x3i=x3l;x3i<=x3h;++x3i){
  	  	  op_vector[l][x1i][x2i][x3i][4][4] = op_vector[l][x1i][x2i][x3i][3][3] = op_vector[l][x1i][x2i][x3i][2][2] = op_vector[l][x1i][x2i][x3i][1][1];
  	  	  op_vector[l][x1i][x2i][x3i][2][1] = op_vector[l][x1i][x2i][x3i][1][2];
  	  	  op_vector[l][x1i][x2i][x3i][3][1] = op_vector[l][x1i][x2i][x3i][1][3];
  	  	  op_vector[l][x1i][x2i][x3i][4][1] = op_vector[l][x1i][x2i][x3i][1][4];
  	  	  op_vector[l][x1i][x2i][x3i][3][2] = -op_vector[l][x1i][x2i][x3i][2][3]; // Rho_V
  	  	  op_vector[l][x1i][x2i][x3i][2][4] = -op_vector[l][x1i][x2i][x3i][4][2]; // Rho_U
  	  	  op_vector[l][x1i][x2i][x3i][4][3] = -op_vector[l][x1i][x2i][x3i][3][4]; // Rho_Q
  }


  return 0;
}