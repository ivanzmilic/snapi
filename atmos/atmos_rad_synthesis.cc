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

// Synthesis refers to the fact that the opacity and emissivity are computed at all wavelengths simultaneously
// in order to maximise speed of execution. We then solve RTE and propagate perturbations wavelength-at-the time

fp_t **** atmosphere::opacity_vector_synth(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, 
	fp_t ****B,fp_t theta,fp_t phi,fp_t * lambda, int nlambda){
  
  	// This is the "synthesis" variant of opacity_vector method which computes wavelength and depth (1D) 
  	// dependent opacity, which is a 4x4 matrix. The thing that bothers me is that last two loop are 
  	// very short so I am not sure how will it work out.

	fp_t **** op_vector = ft4dim(1,nlambda,x3l,x3h,1,4,1,4);
	memset(op_vector[1][x3l][1]+1,0,nlambda*(x3h-x3l+1)*4*4*sizeof(fp_t));
	// Add thompson opacity to this:
	for (int l=1;l<=nlambda;++l)
		for (int x3i=x3l;x3i<=x3h;++x3i){
			fp_t op_thom = 6.65E-25 * Ne[x1l][x2l][x3i];
			for (int s=1;s<=4;++s)
				op_vector[l][x3i][s][s] = op_thom;
	}
	
	for (int a=0;a<natm;++a)
	op_vector = add(atml[a]->opacity_vector_synth(T_in,Ne_in,Vlos,Vt_in,B,theta,phi,lambda,nlambda),op_vector,
		1,nlambda,x3l,x3h,1,4,1,4);
	
	return op_vector;
}

fp_t *** atmosphere::emissivity_vector_synth(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, 
	fp_t ****B,fp_t theta,fp_t phi,fp_t * lambda, int nlambda){
  
 	// This is the "synthesis" variant of emissivity_vector method which computes wavelength and depth (1D) 
  	// dependent opacity, which is a 4x4 matrix. The thing that bothers me is that last two loop are 
  	// very short so I am not sure how will it work out. Same as above, except we also need to take into account emission
  	// which happens according to Planck function

	fp_t *** em_vector = ft3dim(1,nlambda,x3l,x3h,1,4);
	memset(em_vector[1][x3l]+1,0,nlambda*(x3h-x3l+1)*4*sizeof(fp_t));
	fp_t lambda_m = (lambda[nlambda] + lambda[1]) * 0.5;
	// Add thompson opacity to this:
	for (int x3i=x3l;x3i<=x3h;++x3i){
		fp_t local_Planck_function = Planck_f(lambda_m,T[x1l][x2l][x3i]);
		for (int l=1;l<=nlambda;++l){
			fp_t em_thom = 6.65E-25 * Ne[x1l][x2l][x3i] * local_Planck_function;
			for (int s=1;s<=4;++s)
				em_vector[l][x3i][s] = em_thom;
		}
	}
	
	for (int a=0;a<natm;++a)
	em_vector = add(atml[a]->emissivity_vector_synth(T_in,Ne_in,Vlos,Vt_in,B,theta,phi,lambda,nlambda),em_vector,
		1,nlambda,x3l,x3h,1,4);
	
	return em_vector;
}

fp_t ***** atmosphere::opacity_vector_pert_synth(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, 
	fp_t ****B,fp_t theta,fp_t phi,fp_t * lambda, int nlambda){
  
  /*fp_t ***** op = ft5dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  memset(op[x1l][x2l][x3l][1]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  fp_t *** thompson_scattering = thomson_sc(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        op[x1i][x2i][x3i][1][1] = op[x1i][x2i][x3i][2][2] = op[x1i][x2i][x3i][3][3] = op[x1i][x2i][x3i][4][4] = thompson_scattering[x1i][x2i][x3i];

  del_ft3dim(thompson_scattering,x1l,x1h,x2l,x3h,x3l,x3h);

  for (int a=0;a<natm;++a)
    add(atml[a]->opacity_vector(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);*/

  //return op;
}

fp_t **** atmosphere::emissivity_vector_pert_synth(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, 
	fp_t ****B,fp_t theta,fp_t phi,fp_t * lambda, int nlambda){
  
  /*fp_t ***** op = ft5dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  memset(op[x1l][x2l][x3l][1]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  fp_t *** thompson_scattering = thomson_sc(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        op[x1i][x2i][x3i][1][1] = op[x1i][x2i][x3i][2][2] = op[x1i][x2i][x3i][3][3] = op[x1i][x2i][x3i][4][4] = thompson_scattering[x1i][x2i][x3i];

  del_ft3dim(thompson_scattering,x1l,x1h,x2l,x3h,x3l,x3h);

  for (int a=0;a<natm;++a)
    add(atml[a]->opacity_vector(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);*/

  //return op;
}