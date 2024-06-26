#include <math.h>

#include "types.h"
#include "io.h"
#include "const.h"
#include "mem.h"

#include "atmos_pp.h"

// ADDED BY MILIC FOR DEBUGGING PURPOSES
#include <iostream>

fp_t *atmos_pp::anglesetup(fp_t *&th,fp_t *&ph,int &ntp)
//
// angle points (ang) and weights (wmu) for the standard double-Gauss angular quadrature
//
{
//  int nmu=8
//  fp_t ang[]={0.1996052399,0.4548356568,0.7032163503,0.9376088605,1.1502243207,1.3312789639,1.4689536057,1.5509399502};
//  fp_t wmu[]={0.0506142681,0.1111905172,0.1568533229,0.1813418917,0.1813418917,0.1568533229,0.1111905172,0.0506142681};
  //int nmu=5;
  //fp_t ang[]={0.3075109575,0.6931538261,1.0471975512,1.3379321422,1.5238690280};
  //fp_t wmu[]={0.1184634425,0.2393143352,0.2844444444,0.2393143352,0.1184634425};
  //int nmu=3;
  //fp_t ang[]={0.4793425352,1.0471975512,1.4578547042};
  //fp_t wmu[]={0.2777777778,0.4444444444,0.2777777778}; // ={1.0/3.6,4.0/9.0,1.0/3.6};

  int nmu=1;
  fp_t ang[]={1.0 / sqrt(3.0)};
  fp_t wmu[]={1.0};
//
  ntp=2*nmu;
  th=new fp_t [ntp]-1;
  for(int tp=1;tp<=nmu;++tp){
    th[tp]=ang[tp-1];
    th[ntp-tp+1]=pi-ang[tp-1];
  }
  ph=new fp_t [ntp]-1;
  for(int tp=1;tp<=ntp;++tp) ph[tp]=0.0; // value not really relevant in plane-parallel geometry
//
  fp_t *bin=new fp_t [ntp]-1;
  for(int tp=1;tp<=nmu;++tp){
    bin[tp]=wmu[tp-1]*2.0*pi;
    bin[ntp-tp+1]=wmu[tp-1]*2.0*pi;  
  }
  fp_t norm = 0;
  for (int tp=1;tp<=ntp;++tp)
    norm += bin[tp];
  return bin;
}

// plane parallel 1D formal solver
int atmos_pp::formal(fp_t * spatial_grid, fp_t ***S,fp_t ***L,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary) // monochromatic, unidirectional formal solution...
{
 
  fp_t *I=S[x1l][x2l],*oo=op[x1l][x2l],*ee=em[x1l][x2l]; 
  fp_t *l=(L)?L[x1l][x2l]:0;
// 
  fp_t ct=cos(t);
  int dz=(ct>0)?1:-1;
  int zl=(ct>0)?x3l:x3h,zh=(ct>0)?x3h:x3l;
//
  I[zl]=0.0; // I see a small problem there, as the lower boundary condition for the intensity is NOT always zero.
  // EDT by Milic. 
  // Lets try this:
  I[zl] = (ct > 0) ? 0 : ee[zl]/oo[zl]; // This sets intensity to zero if it is ingoing intensity but to the source function 
                                        // if it is the outgoing one.

  for(int z=zl+1;z<zh;z+=dz){ // This says: start from zl+1 until you reach zh, and increase by dz.... but this won't work! 
    fp_t dlf=-ct*(x3[z+dz]-x3[z]); // mu*dz always has the same sign
    fp_t dlb=-ct*(x3[z]-x3[z-dz]);
// will this work well? log(opacity) may be better...?
    fp_t dads=((oo[z]-oo[z-dz])*(dlf/dlb)-(oo[z]-oo[z+dz])*(dlb/dlf))/(dlf+dlb);
    fp_t da2ds2=((oo[z-dz]-oo[z])/dlb+(oo[z+dz]-oo[z])/dlf)/(dlf+dlb);
    fp_t tb=dlb*(oo[z]+dlb*(-dads/2.0+dlb*da2ds2/3.0)); // backward optical depth
    fp_t tf=dlf*(oo[z]+dlf*( dads/2.0+dlf*da2ds2/3.0)); // forward optical depth
//
    fp_t etf=exp(-tf);
    fp_t w0,w1,w2;       // integration weights
//      io.msg(IOL_INFO,"atmos_pp::formal: %d %E %E %E %E %E %E %E\n",z,tf,tb,dlf,dlb,x3[z+dz],x3[z],x3[z-dz]);
    
    if((tf>1.0E-14)){ // Milic: Happy to see that I am not the only one using this...
      if(tf<1.0E-5){ // Milic: This does not really have to be 1E-5, but ok, noone cares at the moment.
        w0=tf-(tf*tf/2.0)+(tf*tf*tf/6.0);
        w1=(tf*tf/2.0)-(tf*tf*tf/3.0);
        w2=(tf*tf*tf/3.0)-(tf*tf*tf*tf/6.0);
      }else{
        w0=1.0-etf;
        w1=w0-tf*etf;
        w2=2.0*w1-tf*tf*etf;
      }

// Calculate the contribution to I_mu
      fp_t sb=ee[z-dz]/oo[z-dz];
      fp_t so=ee[z]/oo[z];
      fp_t sf=ee[z+dz]/oo[z+dz];
      fp_t dsdt=((so-sb)*(tf/tb)-(so-sf)*(tb/tf))/(tf+tb);
      fp_t d2sdt2=((sb-so)/tb+(sf-so)/tf)/(tf+tb);
      if(l) l[z]=w0+(w1*(tf-tb)-w2)/(tf*tb);
      I[z]=I[z-dz]*etf+w0*so+w1*dsdt+w2*d2sdt2;
    }else{
      I[z]=I[z-dz];
      if(l) l[z]=0.0;
    }
  }
// final point
  fp_t dlf=ct*(x3[zh]-x3[zh-dz]); // mu*dz always has the same sign
  fp_t dlb=ct*(x3[zh-dz]-x3[zh-2*dz]);
// will this work out? log(opacity) may be better...?
  fp_t dads=((oo[zh-dz]-oo[zh-2*dz])*(dlf/dlb)-(oo[zh-dz]-oo[zh])*(dlb/dlf))/(dlf+dlb);
  fp_t da2ds2=((oo[zh-2*dz]-oo[zh-dz])/dlb+(oo[zh]-oo[zh-dz])/dlf)/(dlf+dlb);
  fp_t tb=dlb*(oo[zh-dz]+dlb*(-dads/2.0+dlb*da2ds2/3.0)); // backward optical depth
  fp_t tf=dlf*(oo[zh-dz]+dlf*( dads/2.0+dlf*da2ds2/3.0)); // forward optical depth
//
  fp_t etb=exp(-tb);
  fp_t w0,w1,w2;
  if(l) l[zh]=0.0;
  if((tb>1.0E-14)){
    if(tb<1.0E-5){
      w0=tb-(tb*tb/2.0)+(tb*tb*tb/6.0);
      w1=(tb*tb/2.0)-(tb*tb*tb/3.0);
      w2=(tb*tb*tb/3.0)-(tb*tb*tb*tb/6.0);
    }else{
      w0=1.0-etb;
      w1=w0-tb*etb;
      w2=2.0*w1-tb*tb*etb;
    }
// Calculate the formal integral
    fp_t s2=ee[zh]/oo[zh];
    fp_t s0=ee[zh-dz]/oo[zh-dz];
    fp_t s1=ee[zh-2*dz]/oo[zh-2*dz];

    fp_t dsdt=((s0-s2)*(tf/tb)-(s0-s1)*(tb/tf))/(tf+tb);
    fp_t d2sdt2=((s2-s0)/tb+(s1-s0)/tf)/(tf+tb);
    if(l) l[zh]=(w0*tb*(tf+tb)-w1*(tf+2.0*tb)+w2)/(tb*(tf+tb));
    I[zh]=I[zh-dz]*etb+w0*(s0-tb*(dsdt-d2sdt2*tb))+w1*(dsdt-2.0*d2sdt2*tb)+w2*d2sdt2; 
  }else{
    I[zh]=I[zh-dz];
    if(l) l[zh]=0.0;
  }
//
  return 0;
}
int atmos_pp::optical_depth_scale(fp_t ***op,fp_t ***tau,fp_t t,fp_t p) // monochromatic, unidirectional formal solution...
{
  return 0;
}

// plane parallel 1D formal solver
int atmos_pp::formal_with_lambda_operator(fp_t * spatial_grid, fp_t ***S,fp_t **L,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary) // monochromatic, unidirectional formal solution...
{
//
  return 0;
}

// Full Stokes version: Unno-Rachowsky using DELO?
int atmos_pp::formal(fp_t * spatial_grid, fp_t ****S,fp_t ***L,fp_t *****op,fp_t ****em,fp_t t,fp_t p, fp_t boundary) // Formal for polarized // monochromatic, unidirectional formal solution...
{
  if((x1l!=x1h)||(x2l!=x2h)) return io.msg(IOL_ERROR,"atmosphere::sc1dpp: 1D solver requested for region (%d:%d x %d:%d) \n",x1l,x1h,x2l,x2h);
  for(int x3i=x3l;x3i<=x3h;++x3i){
    fp_t **alpha=op[x1l][x2l][x3i];
    fp_t *eta=em[x1l][x2l][x3i];
    for(int s=0;s<=3;++s){
      S[x1l][x2l][x3i][s]=0.0;
      if(L) L[x1l][x2l][x3i]=0.0;
    }
  }
  return 0;
}

int atmos_pp::formal_pert_numerical(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){
  return 0;
}

int atmos_pp::formal_pert_analytical(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){
  return 0;
}

int atmos_pp::formal_pert_jcdti(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){
  return 0;
}
