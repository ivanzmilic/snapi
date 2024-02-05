#include <string.h>
#include <math.h>

#include "types.h"
#include "io.h"
#include "mem.h"

#include "atmos.h"

fp_t *atmosphere::anglesetup(fp_t *&th,fp_t *&ph,int &ntp)
{
  
  return 0;
}

// basic Cartesian 3D formal solver...
int atmosphere::formal(fp_t * spatial_grid, fp_t ***S,fp_t ***L,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary) // monochromatic, unidirectional formal solution...
{
// compute Stokes parameters on grid for specified direction
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        S[x1i][x2i][x3i]=0.0;
        if(L) L[x1i][x2i][x3i]=0.0;
      }
// compute opacity/emission
  io.msg(IOL_INFO,"atmosphere::formal: op[%d][%d][%d] = %E\n",x1l,x2l,x3l,op[x1l][x2l][x3l]);
  io.msg(IOL_INFO,"atmosphere::formal: em[%d][%d][%d] = %E\n",x1l,x2l,x3l,em[x1l][x2l][x3l]);
// 
// RAM management: we must compute the integrated intensity over each transition and all angles. 
// General approach: if we...
//
  return 0;

}

int atmosphere::optical_depth_scale(fp_t ***tau,fp_t ***op,fp_t t,fp_t p){

  return 0;
}

int atmosphere::compute_op_referent(){
  return 0;
}

int atmosphere::compute_tau_referent(){
  return 0;
}

int atmosphere::formal_with_lambda_operator(fp_t * spatial_grid, fp_t ***S,fp_t **lambda_full,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary) // monochromatic, unidirectional formal solution...
{
// compute Stokes parameters on grid for specified direction
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        S[x1i][x2i][x3i]=0.0;
      }
// compute opacity/emission
  io.msg(IOL_INFO,"atmosphere::formal: op[%d][%d][%d] = %E\n",x1l,x2l,x3l,op[x1l][x2l][x3l]);
  io.msg(IOL_INFO,"atmosphere::formal: em[%d][%d][%d] = %E\n",x1l,x2l,x3l,em[x1l][x2l][x3l]);
// 
// RAM management: we must compute the integrated intensity over each transition and all angles. 
// General approach: if we 
//
return 0;
}
int atmosphere::formal_with_responses(fp_t * spatial_grid, fp_t ***S,fp_t *** L, fp_t** response_to_op,fp_t** response_to_em,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary) // monochromatic, unidirectional formal solution...
{
  return 0;
}

// Full Stokes version
int atmosphere::formal(fp_t * spatial_grid, fp_t ****S,fp_t ***L,fp_t *****op,fp_t ****em,fp_t t,fp_t p, fp_t boundary){ // Formal for polarized
// compute Stokes parameters on grid for specified direction
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){ // check for Magneto-Optical Oscillations....
        fp_t **alpha=op[x1i][x2i][x3i];
        fp_t *eta=em[x1i][x2i][x3i];
        for(int s=0;s<=3;++s){
          S[x1i][x2i][x3i][s]=0.0;
          if(L) L[x1i][x2i][x3i]=0.0;
        }
      }
// compute opacity/emission
  io.msg(IOL_INFO,"atmosphere::formal: op[%d][%d][%d] = (%E,%E,%E,%E)\n",x1l,x2l,x3l,op[x1l][x2l][x3l][0][0],op[x1l][x2l][x3l][1][0],op[x1l][x2l][x3l][2][0],op[x1l][x2l][x3l][3][0]);
  io.msg(IOL_INFO,"atmosphere::formal: op[%d][%d][%d] = (%E,%E,%E,%E)\n",x1l,x2l,x3l,op[x1l][x2l][x3l][0][1],op[x1l][x2l][x3l][1][1],op[x1l][x2l][x3l][2][1],op[x1l][x2l][x3l][3][1]);
  io.msg(IOL_INFO,"atmosphere::formal: op[%d][%d][%d] = (%E,%E,%E,%E)\n",x1l,x2l,x3l,op[x1l][x2l][x3l][0][2],op[x1l][x2l][x3l][1][2],op[x1l][x2l][x3l][2][2],op[x1l][x2l][x3l][3][2]);
  io.msg(IOL_INFO,"atmosphere::formal: op[%d][%d][%d] = (%E,%E,%E,%E)\n",x1l,x2l,x3l,op[x1l][x2l][x3l][0][3],op[x1l][x2l][x3l][1][3],op[x1l][x2l][x3l][2][3],op[x1l][x2l][x3l][3][3]);
  io.msg(IOL_INFO,"atmosphere::formal: em[%d][%d][%d] = (%E,%E,%E,%E)\n",x1l,x2l,x3l,em[x1l][x2l][x3l][0],em[x1l][x2l][x3l][1],em[x1l][x2l][x3l][2],em[x1l][x2l][x3l][3]);
// 
// RAM management: we must compute the integrated intensity over each transition and all angles. 
// General approach: if we 
//

// integrate Unno-Rachowsky using DELO/HERM?
  return 0;
}

int atmosphere::formal_pert_numerical(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  return 0;
}

int atmosphere::formal_pert_analytical(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  return 0;
}

int atmosphere::formal_pert_jcdti(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  return 0;
}