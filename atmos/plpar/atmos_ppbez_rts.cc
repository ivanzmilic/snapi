#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "const.h"
#include "mem.h"
#include "mathtools.h"

#include "atmos_ppbez.h"

double const bezier_delimiter = 1E-10;
int const tau_switch_debug = 0;
double const expansion_delimiter = 1E-2;
bool const disable_monotonicity = false;

// Here we will use a Bezier formal solover in order to perform a monochromatic, unidirectional formal solution. We will try to use
// Bezier both for the opacity as well as for the formal solution. Milic_note_to_self: like in the rotating cylinders code. 

int atmos_ppbez::formal(fp_t * spatial_grid, fp_t ***S,fp_t ***L,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary){


  // Input quantities:
  // *** S = input Stokes (or scalar) intensity vector in which eventually the result of the computation is going to written
  // *** L = approximate lambda operator which is being computed in the same go
  // *** op = 3D opacity, for given freq and direction
  // *** em = 3D emissivity, for given freq and direction
  // t      = theta angle
  // p      = phi angle, relevant only in 3D case, but ok we can use it here as well in case of Zeeman or whatever
  // boundary = this specifies what kind of boundary condition are we using. Anything greater then zero means that intensity is specified and equal to boundary. 
  //            if it is -1 then boundary is equal to local source function.

  fp_t *I = S[x1l][x2l]; // "Current" value of intensity. I.e. the one for this "column" 
  fp_t *oo = op[x1l][x2l]; // Same for the opacity
  fp_t *ee = em[x1l][x2l]; // Same for emissivity
  fp_t *l = (L) ? L[x1l][x2l] : 0; // Now, if we pass zero, then do not compute the local operator, if we pass something we, surprise-surprise, do!
// 
  fp_t ct=cos(t); // Cosinus of the direction. If you are for first time looking @ this, do not be confused, the directions are the opposite: 
                  // t = 0 indicates ingoing ray while t = pi indicates the outgoing ray, i.e. angle is with respect to the direction
                  // observer - center of the star.
  
  int dz = (ct>0) ? 1 : -1; // This is increment along z which determines what do we do, i.e. if we increase or decrease the z
  int zl = (ct>0) ? x3l : x3h; // Where do we start from.
  int zh = (ct>0) ? x3h : x3l; // Where do we end.
//
  // Boundary condition for the intensity.
  I[zl]= (ct > 0) ? 0.0 : ee[zl] / oo[zl];

  if(ct < 0 && boundary >= 0) // lower boundary condition
      I[zl] = boundary;

  // And appropriate local operator.
  if(L && ct >= 0) l[zl] = 0.0; // Actually here is something that we might want to examine further; i.e. what happens with the local operator
                               // at the lower boundary condition, what makes sense is to set 1.0 at the lower one, no?
  if(L && ct <= 0) {

    if (boundary == -1)
      l[zl] = 1.0;
    else 
      l[zl] = 0.0;
  }

  // We will now compute the optical depth scale from here.
  // First derivatives of the opacity.
  int ND = x3h-x3l+1;
  fp_t * delta_tau = new fp_t [ND] - x3l;
  fp_t * opp_derivative = new fp_t [ND] - x3l;
  opp_derivative[zl] = (oo[zl+dz] - oo[zl]) / (spatial_grid[zl+dz] - spatial_grid[zl]) * (-ct);
  opp_derivative[zh] = (oo[zh] - oo[zh-dz]) / (spatial_grid[zh] - spatial_grid[zh-dz]) * (-ct);
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    opp_derivative[z] = 0.0; // default value
    fp_t h2 = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct); 
    fp_t h1 = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
    fp_t d2 = (oo[z+dz] - oo[z]) / h2;
    fp_t d1 = (oo[z] - oo[z-dz]) / h1;
    if (d1*d2 > 0.0 || disable_monotonicity){
      fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
      opp_derivative[z] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
    }
  }
  // Then from there the optical depth scale:
  delta_tau[zl] = 0.0; // By default.
  for (int z = zl+dz;(z-zh-dz)*dz<0;z+=dz){
    fp_t C0 = oo[z] - (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z];
    fp_t C1 = oo[z-dz] + (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z-dz];
    fp_t C = 0.5*(C0+C1);
    fp_t delta_t = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 3.0 * (oo[z-dz] + oo[z] + C);
    delta_tau[z] = delta_t;
    //if (delta_t <= 0.0){
    //  printf("Negative optical depth! %e \n", delta_t);
    //  exit(1);
    //}
  }

  // Compute the derivatives of the source function before advancing:
  fp_t * s_derivative = new fp_t[ND] - x3l;
  s_derivative[zl] = (ee[zl+dz]/oo[zl+dz] - ee[zl]/oo[zl]) / delta_tau[zl+dz];
  s_derivative[zh] = (ee[zh]/oo[zh] - ee[zh-dz]/oo[zh-dz]) / delta_tau[zh];
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    s_derivative[z] = 0.0; // default value
    fp_t h2 = delta_tau[z+dz];
    fp_t h1 = delta_tau[z];
    fp_t d2 = (ee[z+dz]/oo[z+dz] - ee[z]/oo[z]) / h2;
    fp_t d1 = (ee[z]/oo[z] - ee[z-dz]/oo[z-dz]) / h1;
    if(d1*d2>0.0 || disable_monotonicity){
      fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
      s_derivative[z] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
    }
  }

  // Now we need to formally solve RTE, layer by layer, EXCEPT THE LAST LAYER, which we do at the end.
  for(int z = zl+dz; (z-zh)*dz < 0; z+=dz){ // This one should work both ways, like Thirteen, and I have verified it works.
    
    fp_t tb = delta_tau[z];
    fp_t tf = delta_tau[z+dz];  
    // APPROXIMATE TO CHECK WHAT IS THE PROBLEM:
    if (tau_switch_debug){
      fp_t dlb = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
      fp_t dlf = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct);
      tb = dlb * (oo[z-dz] + oo[z]) * 0.5;
      tf = dlf * (oo[z] + oo[z+dz]) * 0.5;
    }
    
    // Now, when we have the upwind and dowind (backward and forward) optical paths, we can also compute the weights and control point. 
    // And, subsequently the contribution to the lambda operator.
    
    fp_t etb = (tb < expansion_delimiter) ? (1.0 - tb + tb*tb*0.5 - tb*tb*tb/6.0) : exp(-tb); // Exponent of the backward optical depth, will become useful.
    fp_t w0, w1, w2;       // Integration weights, for local, previous and control point.

    if (tf > bezier_delimiter){ // Bezier, ie. higher then 1st order
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb*(1.0/3.0 - tb*(1.0/4.0-tb/10.0));
        w1 = tb/3.0*(1.0-tb/4.0*(1.0-tb/5.0));
        w2 = tb*(1.0/3.0-tb*(1.0/6.0-tb/20.0));
      }
      else {
        w0 = (2.0 - etb * (tb * tb + 2.0 * tb + 2.0)) / tb / tb;
        //if (z == 29)
        //  printf("Formal  : %5.15e %5.15e %5.15e %5.15e \n",tb,etb, w0,(2.0 - etb * (tb * tb + 2.0 * tb + 2.0)) / tb / tb);
        w1 = 1.0 - 2.0 * (etb + tb - 1.0) / tb / tb; 
        w2 = 2.0 * (tb - 2.0  + etb * (tb + 2.0)) / tb / tb;
        
      }
    }
    else { // Linear
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb * (0.5 - tb * (1.0/3.0  - tb/8.0));
        w1 = 0.5 * tb * (1.0 - tb/3.0 * (1.0 - tb/4.0));
        w2 = 0.0;
      }
      else {
        w0 = 1.0 / tb - etb * (1.0 + 1.0 / tb);
        w1 = 1.0 - 1.0 / tb + etb / tb;
        w2 = 0.0;
      }
    }
    //if (z == 29)
      //printf("Formal : %5.15e %5.15e %5.15e \n",tf, tb, w0);

    // Before actually performing formal solution we have to compute the value of the control point, which is a bit slower when you use
    // BESSER but what can you do. (We can always revert o de la Cruz Rodriguez & Piskunov 2013)
    fp_t s0 = ee[z-dz] / oo[z-dz]; // Upwind source function.
    fp_t s1 = ee[z] / oo[z]; // Local source function.
    fp_t s2 = ee[z+dz] / oo [z+dz]; // Downwind source function.
    fp_t sp, C1, C0, C;

    if (tf > bezier_delimiter){
      
      C0 = s1 - tb/2.0*s_derivative[z];
      C1 = s0 + tb/2.0*s_derivative[z-dz];
      C = (C0+C1)*0.5;
    }
    else 
      C = 0.0;

    // Finally, calculate the contribution to I_mu
    I[z] = I[z-dz] * etb + w0 * s0 + w1 * s1 + w2 * C;

    //if (z == 29){
    //  printf("Formal. \n");
      //printf("%5.15e %5.15e %5.15e %5.15e %5.15e %5.15e %5.15e %5.15e \n", I[z-dz],etb,w0,s0,w1,s1,w2,C);
    //  printf("%5.15e \n", w0);
    //}

    //printf("%d %e %e %e %e %e %e  \n", z, etb, I[z-dz], w0, s0, w1, s1, w2, C); 
    
    // And compute the local operator
    if(L) l[z] = (w1+0.5*w2);
    //if (L){
    //  if (l[z] < 0){
    //    printf("Error! l[z] = %e @ %d. w1 = %e w2 = %e tb = %e \n", l[z],z,w1,w2,tb);
    //    exit(1);
    //  }
    //}
  }

  // Still the final point in the atmosphere is not covered here. That one we will compute in the linear approximation, both for 
  // the opacity and for the source function.
  if (zl != zh){
    int z = zh; // Just to be absolutely sure as I am a nub programmer sometimes (Milic)
    fp_t tb = delta_tau[z];
   
    fp_t etb = (tb < expansion_delimiter) ? 1.0 - tb + tb*tb*0.5 - tb*tb*tb/6.0 : exp(-tb);
    fp_t s0 = ee[z-dz] / oo[z-dz], s1 = ee[z] / oo[z];
    fp_t w0, w1;
    if (tb < expansion_delimiter){ // Expansion
      w0 = tb * (0.5 - tb * (1.0/3.0  - tb/8.0));;
      w1 = 0.5 * tb * (1.0 - tb/3.0 * (1.0 - tb/4.0));;
    }
    else {
      w0 = 1.0 / tb - etb * (1.0 + 1.0 / tb);
      w1 = 1.0 - 1.0 / tb + etb / tb;
    }
    // And finally compute the formal solution man:
    I[z] = I[z-dz] * etb + w0 * s0 + w1 * s1;
    // And compute the local operator:
    if(L) l[z] = w1;
  }

  //for (int x3i=x3l;x3i<=x3h;++x3i)
    //l[x3i] = 0.0;
    
  delete [](delta_tau+x3l);
  delete [](opp_derivative+x3l);
  delete [](s_derivative+x3l);
  return 0;
}

int atmos_ppbez::optical_depth_scale(fp_t ***tau,fp_t ***op,fp_t t,fp_t p){

  fp_t *tt = tau[x1l][x2l]; // "Current" value of intensity. I.e. the one for this "column" 
  fp_t *oo = op[x1l][x2l]; // Same for the opacity
   
  fp_t ct=cos(t); // Cosinus of the direction. If you are for first time looking @ this, do not be confused, the directions are the opposite: 
                  // t = 0 indicates ingoing ray while t = pi indicates the outgoing ray, i.e. angle is with respect to the direction
                  // observer - center of the star.
  
  int dz = (ct>0) ? 1 : -1; // This is increment along z which determines what do we do, i.e. if we increase or decrease the z
  int zl = (ct>0) ? x3l : x3h; // Where do we start from.
  int zh = (ct>0) ? x3h : x3l; // Where do we end.
//
  // We will now compute the optical depth scale from here.
  // First derivatives of the opacity.
  int ND = x3h-x3l+1;
  fp_t * opp_derivative = new fp_t [ND] - x3l;
  opp_derivative[zl] = (oo[zl+dz] - oo[zl]) / (x3[zl+dz] - x3[zl]) * (-ct);
  opp_derivative[zh] = (oo[zh] - oo[zh-dz]) / (x3[zh] - x3[zh-dz]) * (-ct);
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    opp_derivative[z] = 0.0; // default value
    fp_t h2 = (x3[z+dz] - x3[z]) / (-ct); 
    fp_t h1 = (x3[z] - x3[z-dz]) / (-ct);
    fp_t d2 = (oo[z+dz] - oo[z]) / h2;
    fp_t d1 = (oo[z] - oo[z-dz]) / h1;
    if (d1*d2 > 0 || disable_monotonicity){
      fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
      opp_derivative[z] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
    }
  }
  // Then from there the optical depth scale:
  tt[zl] = 1E-8; // By default.
  for (int z = zl+dz;(z-zh-dz)*dz<0;z+=dz){
    fp_t C0 = oo[z] - (x3[z] - x3[z-dz]) / (-ct) / 2.0 * opp_derivative[z];
    fp_t C1 = oo[z-dz] + (x3[z] - x3[z-dz]) / (-ct) / 2.0 * opp_derivative[z-dz];
    fp_t C = 0.5*(C0+C1);
    fp_t delta_t = (x3[z] - x3[z-dz]) / (-ct) / 3.0 * (oo[z-dz] + oo[z] + C);
    tt[z] = tt[z-dz] + delta_t;
  }

  delete [](opp_derivative+x3l);

  return 0;
}

int atmos_ppbez::compute_op_referent(){

  // First we have to compute the opacity.
  
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        fp_t T_local = T[x1i][x2i][x3i];
        chemeq(atml, natm, T_local, Nt[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        for (int a=0;a<natm;++a) atml[a]->lte(T_local, Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        op_referent[x1i][x2i][x3i] = opacity_continuum(T_local,Ne[x1i][x2i][x3i], 500E-7, x1i,x2i,x3i);
      }
  return 0;
}

int atmos_ppbez::compute_op_referent_derivative(){

  op_referent_derivative = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_referent_derivative[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){

        // Perturb temperature localy
        T[x1i][x2i][x3i] += delta_T*0.5;
        chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        for (int a=0;a<natm;++a) atml[a]->lte(T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        op_referent_derivative[1][x3i][x1i][x2i][x3i] = opacity_continuum(T[x1i][x2i][x3i],Ne[x1i][x2i][x3i], 500E-7, x1i,x2i,x3i);

        T[x1i][x2i][x3i] -= delta_T;
        chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        for (int a=0;a<natm;++a) atml[a]->lte(T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        op_referent_derivative[1][x3i][x1i][x2i][x3i] -= opacity_continuum(T[x1i][x2i][x3i],Ne[x1i][x2i][x3i], 500E-7, x1i,x2i,x3i);
        op_referent_derivative[1][x3i][x1i][x2i][x3i] /= delta_T;
        T[x1i][x2i][x3i] += delta_T*0.5;

        // Then the density
        fp_t delta_Nt = Nt[x1i][x2i][x3i] * delta_Nt_frac;
        Nt[x1i][x2i][x3i] += delta_Nt*0.5;
        chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        for (int a=0;a<natm;++a) atml[a]->lte(T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        op_referent_derivative[2][x3i][x1i][x2i][x3i] = opacity_continuum(T[x1i][x2i][x3i],Ne[x1i][x2i][x3i], 500E-7, x1i,x2i,x3i);

        Nt[x1i][x2i][x3i] -= delta_Nt;
        chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        for (int a=0;a<natm;++a) atml[a]->lte(T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        op_referent_derivative[2][x3i][x1i][x2i][x3i] -= opacity_continuum(T[x1i][x2i][x3i],Ne[x1i][x2i][x3i], 500E-7, x1i,x2i,x3i);
        op_referent_derivative[2][x3i][x1i][x2i][x3i] /= delta_Nt;
        Nt[x1i][x2i][x3i] += delta_Nt*0.5;

        chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
        for (int a=0;a<natm;++a) atml[a]->lte(T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], x1i, x2i, x3i);
  }
  return 0;
}


int atmos_ppbez::compute_tau_referent(){

  optical_depth_scale(tau_referent, op_referent, 0.0,0.0);

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        tau_referent[x1i][x2i][x3i] *= -1.0;

  return 0;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------

// Special, response function related version:

int atmos_ppbez::formal_with_lambda_operator(fp_t * spatial_grid, fp_t ***S,fp_t **lambda_full,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary){
  
  return 0;
}

int atmos_ppbez::formal_with_responses(fp_t * spatial_grid, fp_t *** S, fp_t *** L, fp_t ** response_to_op, fp_t ** response_to_em,fp_t ***op,fp_t ***em,fp_t t,fp_t p, fp_t boundary){
    
  // This is the same function as the one which performs the formal solution except it also returns the responses of all the intensities to all opacities
  // and emissivities
  // Input quantities:
  // *** S = input Stokes (or scalar) intensity vector in which eventually the result of the computation is going to written
  // *** L = approximate lambda operator which is being computed in the same go
  // ** response_to_op = Nx3 x Nx3 array which returns response of each intensity to each opacity
  // ** response_to_em = Nx3 x Nx3 array which returns response of each intensity to each emissivity, this should be very similar to lambda operator
  // *** op = 3D opacity, for given freq and direction
  // *** em = 3D emissivity, for given freq and direction
  // t      = theta angle
  // p      = phi angle, relevant only in 3D case, but ok we can use it here as well in case of Zeeman or whatever
  // boundary = this specifies what kind of boundary condition are we using. Anything greater then zero means that intensity is specified and equal to boundary. 
  //            if it is -1 then boundary is equal to local source function.
  int ND = x3h-x3l+1;
  fp_t *I = S[x1l][x2l]; // "Current" value of intensity. I.e. the one for this "column" 
  fp_t *oo = op[x1l][x2l]; // Same for the opacity
  fp_t *ee = em[x1l][x2l]; // Same for emissivity
  fp_t *l;
  if (L)
    l  = L[x1l][x2l]; // And the local operator
//
  fp_t ct=cos(t); // Cosinus of the direction. If you are for first time looking @ this, do not be confused, the directions are the opposite: 
                  // t = 0 indicates ingoing ray while t = pi indicates the outgoing ray, i.e. angle is with respect to the direction
                  // observer - center of the star.
  
  int dz = (ct>0) ? 1 : -1; // This is increment along z which determines what do we do, i.e. if we increase or decrease the z
  int zl = (ct>0) ? x3l : x3h; // Where do we start from.
  int zh = (ct>0) ? x3h : x3l; // Where do we end.

  memset(response_to_op[x3l]+x3l,0,ND*ND*sizeof(fp_t));
  memset(response_to_em[x3l]+x3l,0,ND*ND*sizeof(fp_t));
 
//
  // Boundary condition for the intensity.
  I[zl]= (ct > 0) ? 0.0 : ee[zl] / oo[zl];

  if(ct < 0 && boundary >= 0) // lower boundary condition
      I[zl] = boundary;

  // And appropriate local operator. From which also follow appropriate responses:
  if(ct >= 0){ 
    if (L)
      l[zl] = 0.0;
    for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
      response_to_op[zl][zp] = 0.0;
      response_to_em[zl][zp] = 0.0;
    }
  }

  if(ct <= 0) {
    if (L){
      if (boundary == -1)
        l[zl] = 1.0;
      else 
        l[zl] = 0.0;
    }

    // And sort out responses:
    if (boundary == -1){
      response_to_op[zl][zl] = -ee[zl]/oo[zl]/oo[zl];
      response_to_em[zl][zl] = 1.0/oo[zl];
      for (int zp=zl+dz;(zp-zh-dz)*dz<0;zp+=dz){
        response_to_op[zl][zp] = 0.0;
        response_to_em[zl][zp] = 0.0;
      } 
    }
    else {
      for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
        response_to_op[zl][zp] = 0.0;
        response_to_em[zl][zp] = 0.0;
      }
    }
  }

  // We will now compute the optical depth scale from here.
  // First derivatives of the opacity.

  fp_t * delta_tau = new fp_t [ND] - x3l;
  fp_t * opp_derivative = new fp_t [ND] - x3l;
  opp_derivative[zl] = (oo[zl+dz] - oo[zl]) / (spatial_grid[zl+dz] - spatial_grid[zl]) * (-ct);
  opp_derivative[zh] = (oo[zh] - oo[zh-dz]) / (spatial_grid[zh] - spatial_grid[zh-dz]) * (-ct);
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    opp_derivative[z] = 0.0; // default value
    fp_t h2 = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct); 
    fp_t h1 = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
    fp_t d2 = (oo[z+dz] - oo[z]) / h2;
    fp_t d1 = (oo[z] - oo[z-dz]) / h1;
    if(d1*d2 > 0.0 || disable_monotonicity){
      fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
      opp_derivative[z] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
    }
  }
  // Then from there the optical depth scale:
  delta_tau[zl] = 0.0; // By default.
  for (int z = zl+dz;(z-zh-dz)*dz<0;z+=dz){
    fp_t C0 = oo[z] - (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z];
    fp_t C1 = oo[z-dz] + (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z-dz];
    fp_t C = 0.5*(C0+C1);
    fp_t delta_t = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 3.0 * (oo[z-dz] + oo[z] + C);
    delta_tau[z] = delta_t;
  }
  //printf("tau[zh] = %e \n", tau[zh]);
  
  // Now we start from the beggining and start computing the responses. So, the first step is to find the derivative of the numerical
  // derivative with respect to the opacity. (This is not the error, we do find how the values of the derivative of the opacity, computed in the 
  // numerical way behave with the changes in the opacity)
  fp_t ** d_der_per_d_chi = ft2dim(x3l,x3h,x3l,x3h);
  
  memset(d_der_per_d_chi[x3l]+x3l,0,ND*ND*sizeof(fp_t));
  // First the first one:
  if (zh != zl){ 
    d_der_per_d_chi[zl][zl] = -1.0 / (spatial_grid[zl+dz]-spatial_grid[zl]) * (-ct);
    d_der_per_d_chi[zl][zl+dz] = 1.0 / (spatial_grid[zl+dz]-spatial_grid[zl]) * (-ct);
    // Then the last one:
    d_der_per_d_chi[zh][zh-dz] = -1.0 / (spatial_grid[zh]-spatial_grid[zh-dz]) * (-ct);
    d_der_per_d_chi[zh][zh] = 1.0 / (spatial_grid[zh]-spatial_grid[zh-dz]) * (-ct);
    // And then the ones in between:
    for (int i=zl+dz;(i-zh)*dz<0;i+=dz){
      fp_t h2 = (spatial_grid[i+dz] - spatial_grid[i]) / (-ct); 
      fp_t h1 = (spatial_grid[i] - spatial_grid[i-dz]) / (-ct);
      fp_t d2 = (oo[i+dz] - oo[i]) / h2;
      fp_t d1 = (oo[i] - oo[i-dz]) / h1;

      fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
      if (d1*d2 > 0.0 || disable_monotonicity){
        fp_t denominator = (alpha*d1 + (1.0 - alpha)*d2);
        // Depends on the previous one through d1:
        d_der_per_d_chi[i][i-dz] = (-1.0/h1*d2*denominator + alpha/h1*d1*d2)/denominator/denominator;
        // Depends on local one, both through d1 and d2:
        d_der_per_d_chi[i][i] = (1.0/h1*d2*denominator - 1.0/h2*d1*denominator - alpha/h1*d1*d2 + (1.0-alpha)/h2*d1*d2)/denominator/denominator;
        // And depends on the next one, through d2:
        d_der_per_d_chi[i][i+dz] = (1.0/h2*d1*denominator - (1.0-alpha)/h2*d1*d2)/denominator/denominator;
      }
    }
  }

  //for (int x3i=x3l;x3i<=x3h;++x3i)
    //for (int x3ii=x3l;x3ii<=x3h;++x3ii)
      //printf("%d %d %e \n", x3i, x3ii, d_der_per_d_chi[x3i][x3ii]);
  
  // After this one we compute the derivative of delta with respect to each opacity. It should be tridiagonal-ish matrix:
  fp_t ** d_delta_per_d_chi = ft2dim(x3l,x3h,x3l,x3h);
  memset(d_delta_per_d_chi[x3l]+x3l,0,ND*ND*sizeof(fp_t));
  // Start by adding trivial things:
  if (zl != zh){
    for (int i=zl+dz;(i-zh-dz)*dz<0;i+=dz){ // delta, we skip first as it is zero by default

      fp_t h_i = (spatial_grid[i] - spatial_grid[i-dz]) / (-ct);
      // Depends explicitily on local value of the opacity:
      d_delta_per_d_chi[i][i] = 1.0/2.0 * h_i;
      // Expllicitly on the previous opacity:
      d_delta_per_d_chi[i][i-dz] = 1.0/2.0 * h_i;
      // And then implicitly on the local and previous opacity:
      for (int j=zl;(j-zh-dz)*dz<0;j+=dz)
        d_delta_per_d_chi[i][j] += h_i*h_i/12.0 * d_der_per_d_chi[i-dz][j] - h_i*h_i/12.0 * d_der_per_d_chi[i][j];
    }
  }
  // Be aware that this is derivative of DELTA, that is optical depth spacing between to consecutive points, not the deviative of the optical depth
  // itself! 
  //for (int x3i=x3l;x3i<=x3h;++x3i)
    //printf("%d %e \n",x3i,d_delta_per_d_chi[x3i][x3l]);

// --------------------------------------------------------------------------------------------------------------------------------------------------

  // Compute the derivatives of the source function before advancing:
  fp_t * s_derivative = new fp_t[ND] - x3l;
  
  // We will then need two additional arrays, dS'/d delta and d S' / dS:
  fp_t ** d_sp_d_delta = ft2dim(x3l,x3h,x3l,x3h);
  fp_t ** d_sp_d_s = ft2dim(x3l,x3h,x3l,x3h);
  memset(d_sp_d_delta[x3l]+x3l,0,ND*ND*sizeof(fp_t));
  memset(d_sp_d_s[x3l]+x3l,0,ND*ND*sizeof(fp_t));
  if (zl != zh){
    // First point:
    s_derivative[zl] = (ee[zl+dz]/oo[zl+dz] - ee[zl]/oo[zl]) / delta_tau[zl+dz];
    
    d_sp_d_delta[zl][zl+dz] = -(ee[zl+dz]/oo[zl+dz] - ee[zl]/oo[zl]) / delta_tau[zl+dz] / delta_tau[zl+dz];
    
    d_sp_d_s[zl][zl] = -1.0 / delta_tau[zl+dz];
    d_sp_d_s[zl][zl+dz] = 1.0 / delta_tau[zl+dz];

    // Last point:
    s_derivative[zh] = (ee[zh]/oo[zh] - ee[zh-dz]/oo[zh-dz]) / delta_tau[zh];
    
    d_sp_d_delta[zh][zh] = -(ee[zh]/oo[zh] - ee[zh-dz]/oo[zh-dz]) / delta_tau[zh] / delta_tau[zh];
   
    d_sp_d_s[zh][zh] = 1.0 / delta_tau[zh];
    d_sp_d_s[zh][zh-dz] = -1.0 / delta_tau[zh];

    // Ones in between:
    for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
      s_derivative[z] = 0.0; // default value
      fp_t h2 = delta_tau[z+dz];
      fp_t h1 = delta_tau[z];
      fp_t d2 = (ee[z+dz]/oo[z+dz] - ee[z]/oo[z]) / h2;
      fp_t d1 = (ee[z]/oo[z] - ee[z-dz]/oo[z-dz]) / h1;

      if(d1*d2>0.0 || disable_monotonicity){
        fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
        s_derivative[z] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
        fp_t denominator = (alpha*d1 + (1.0-alpha)*d2);

        // And now there are are all individual derivatives, first with respect to source function.
        d_sp_d_s[z][z-dz] = -d2/h1*(denominator-alpha*d1) /denominator/denominator;
        d_sp_d_s[z][z] = ((d2/h1 - d1/h2) * denominator - (alpha/h1 - (1.0-alpha)/h2)*d1*d2)/denominator/denominator;
        d_sp_d_s[z][z+dz] = d1/h2 *(denominator-(1.0-alpha)*d2)/denominator/denominator;

        if (isnan(d_sp_d_s[z][z])){
          printf("%d %d %e %e %e %e \n", z, dz, d1, d2, denominator, alpha);
          printf("%e %e %e %e", s_derivative[z], ee[z+dz]/oo[z+dz], ee[z]/oo[z], h2);
          exit(1);
        }

        // Once with respect to the deltas are harder because alpha also depends on delta. So first compute that.
        fp_t d_alpha_i = h2/3.0/(h1+h2)/(h1+h2);
        fp_t d_alpha_ip = -h1/3.0/(h1+h2)/(h1+h2);

        d_sp_d_delta[z][z] = (-d1/h1*d2*denominator + alpha*d1/h1*d1*d2 - d_alpha_i*d1*d1*d2 + d_alpha_i*d2*d2*d1) /denominator/denominator;
        d_sp_d_delta[z][z+dz] =  (-d2/h2*d1*denominator + (1.0-alpha)*d2/h2*d1*d2 - d_alpha_ip*d1*d1*d2 +d_alpha_ip*d1*d2*d2)/denominator/denominator;
      }
    }
  }

  // Now we need to formally solve RTE, layer by layer, EXCEPT THE LAST LAYER, which we do at the end.
  for(int z = zl+dz; (z-zh)*dz < 0; z+=dz){ // This one should work both ways, like Thirteen, and I have verified it works.
    
    fp_t tb = delta_tau[z];
    fp_t tf = delta_tau[z+dz];  
    
    // APPROXIMATE TO CHECK WHAT IS THE PROBLEM:
    if (tau_switch_debug){
      fp_t dlb = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
      fp_t dlf = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct);
      tb = dlb * (oo[z-dz] + oo[z]) * 0.5;
      tf = dlf * (oo[z] + oo[z+dz]) * 0.5;
    }
    
    // Now, when we have the upwind and dowind (backward and forward) optical paths, we can also compute the weights and control point. 
    // And, subsequently the contribution to the lambda operator.
    
    fp_t etb = (tb < expansion_delimiter) ? (1.0 - tb + tb*tb*0.5 - tb*tb*tb/6.0) : exp(-tb); // Exponent of the backward optical depth, will become useful.
    fp_t w0, w1, w2;       // Integration weights, for local, previous and control point.

    // And then we also have to add the derivatives of all the stuff:
    fp_t d_etb, d_w0, d_w1, d_w2;

    d_etb = - etb;

    if (tf > bezier_delimiter){ // Bezier, ie. higher then 1st order
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb*(1.0/3.0 - tb*(1.0/4.0-tb/10.0));
        d_w0 = 1.0/3.0 - tb*(1.0/2.0 - 3.0*tb/10.0);
        w1 = tb/3.0*(1.0-tb/4.0*(1.0-tb/5.0));
        d_w1 = 1.0/3.0  - tb/2.0*(1.0/3.0 - tb/10.0);
        w2 = tb*(1.0/3.0-tb*(1.0/6.0-tb/20.0));
        d_w2 = 1.0/3.0 - tb*(1.0/3.0 - 3.0*tb/20.0);
      }
      else {
        w0 = (2.0 - etb * (tb * tb + 2.0 * tb + 2.0)) / tb / tb;
        //if (z == 29)
        //  printf("FormalR : %5.15e %5.15e %5.15e %5.15e \n",tb,etb, w0, (2.0 - etb * (tb * tb + 2.0 * tb + 2.0)) / tb / tb);
        d_w0 = -2.0*w0/tb + etb;
        w1 = 1.0 - 2.0 * (etb + tb - 1.0) / tb / tb;
        d_w1 = 2.0 *(1.0-w1) / tb - 2.0*(1.0-etb)/tb/tb; 
        w2 = 2.0 * (tb - 2.0  + etb * (tb + 2.0)) / tb / tb;
        d_w2 = -2.0*w2/tb + 2.0/tb/tb * (1.0-tb*etb-etb);
    
      }
      
    }
    else { // Linear
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb * (0.5 - tb * (1.0/3.0  - tb/8.0));
        d_w0 = 1.0/2.0 - tb*(2.0/3.0 - 3.0*tb/8.0);
        w1 = 0.5 * tb * (1.0 - tb/3.0 * (1.0 - tb/4.0));
        d_w1 = 0.5 - tb*(1.0/3.0 - tb/8.0);
        w2 = 0.0;

      }
      else {
        w0 = 1.0 / tb - etb * (1.0 + 1.0 / tb);
        d_w0 = -w0/tb + etb;
        w1 = 1.0 - 1.0 / tb + etb / tb;
        d_w1 = (1.0-w1)/tb - etb/tb;
        w2 = 0.0;
      }
      d_w2 = 0.0;
    }

    //if (z == 29)
      // /printf("FormalR : %5.15e %5.15e %5.15e \n",tf, tb, w0);

    // Before actually performing formal solution we have to compute the value of the control point, which is a bit slower when you use
    // BESSER but what can you do. (We can always revert o de la Cruz Rodriguez & Piskunov 2013)
    fp_t s0 = ee[z-dz] / oo[z-dz]; // Upwind source function.
    fp_t s1 = ee[z] / oo[z]; // Local source function.
    fp_t s2 = ee[z+dz] / oo [z+dz]; // Downwind source function.
    fp_t sp, C1, C0, C;

    if (tf > bezier_delimiter){
      
      C0 = s1 - tb/2.0*s_derivative[z];
      C1 = s0 + tb/2.0*s_derivative[z-dz];
      C = (C0+C1)*0.5;
    }
    else 
      C = 0.0;

    // Finally, calculate the contribution to I_mu
    I[z] = I[z-dz] * etb + w0 * s0 + w1 * s1 + w2 * C;

    //printf("Error at : %d \n", z);
    //if (z == 29){
    //  printf("Formal. \n");
      //printf("%5.15e %5.15e %5.15e %5.15e %5.15e %5.15e %5.15e %5.15e \n", I[z-dz],etb,w0,s0,w1,s1,w2,C);
    //  printf("%5.15e \n", w0);
    //}
    //exit(0);

  
    // And then we can also compute the derivative with respect to all opacities and emissivities. 
    // First the emissivities:
    // -----------------------------------------------------------------------------------------------
    // There are two which are, without a doubt local:
    response_to_em[z][z] = (w1 + 0.5 * w2)/oo[z];
    response_to_em[z][z-dz] = (w0+0.5*w2)/oo[z-dz];
    // Then some which are non local
    for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
      response_to_em[z][zp] += 0.25 * (-w2*tb*d_sp_d_s[z][zp] + w2*tb*d_sp_d_s[z-dz][zp]) / oo[zp];
      response_to_em[z][zp] += response_to_em[z-dz][zp] * etb;
    }
    // -----------------------------------------------------------------------------------------------
    // The story starts the same with the opacity:
    response_to_op[z][z] = -(w1 + 0.5 * w2)*ee[z]/oo[z]/oo[z];
    response_to_op[z][z-dz] = -(w0+0.5*w2)*ee[z-dz]/oo[z-dz]/oo[z-dz];
    // Then some which are non local
    for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
      response_to_op[z][zp] += 0.25 * (-w2*tb*d_sp_d_s[z][zp] + w2*tb*d_sp_d_s[z-dz][zp]) * (-ee[zp]/oo[zp]/oo[zp]);
      response_to_op[z][zp] += response_to_op[z-dz][zp] * etb;
    }
    // Then there is also dependance on delta:
    for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
      response_to_op[z][zp] += -I[z-dz] * etb * d_delta_per_d_chi[z][zp];
      response_to_op[z][zp] += (s0 * d_w0  + s1*d_w1 + C * d_w2) * d_delta_per_d_chi[z][zp];
      response_to_op[z][zp] += w2 * 0.25 * d_delta_per_d_chi[z][zp] * (s_derivative[z-dz] - s_derivative[z]);
      for (int zpp=zl;(zpp-zh-dz)*dz<0;zpp+=dz)
        response_to_op[z][zp] += w2 * 0.25 * tb * (d_sp_d_delta[z-dz][zpp] - d_sp_d_delta[z][zpp]) * d_delta_per_d_chi[zpp][zp];  
    }

    // And compute the local operator
    if(L){ 
      l[z] = (w1+0.5*w2);
    }    
  }

  // Still the final point in the atmosphere is not covered here. That one we will compute in the linear approximation, both for 
  // the opacity and for the source function.
  if (zh != zl){ // If it is not one and the same point:
    int z = zh; // Just to be absolutely sure as I am a nub programmer sometimes (Milic)
    fp_t tb = delta_tau[z];
   
    fp_t etb = (tb < expansion_delimiter) ? 1.0 - tb + tb*tb*0.5 - tb*tb*tb/6.0 : exp(-tb);
    fp_t s0 = ee[z-dz] / oo[z-dz], s1 = ee[z] / oo[z];
    fp_t w0, w1;
    fp_t d_w0, d_w1;
    if (tb < expansion_delimiter){ // Expansion
      w0 = tb * (0.5 - tb * (1.0/3.0  - tb/8.0));
      d_w0 = 1.0/2.0 - tb*(2.0/3.0 - 3.0*tb/8.0);
      w1 = 0.5 * tb * (1.0 - tb/3.0 * (1.0 - tb/4.0));
      d_w1 = 0.5 - tb*(1.0/3.0 - tb/8.0);
    }
    else {
      w0 = 1.0 / tb - etb * (1.0 + 1.0 / tb);
      d_w0 = -w0/tb + etb;
      w1 = 1.0 - 1.0 / tb + etb / tb;
      d_w1 = (1.0-w1)/tb - etb/tb;
    }
        
    // And finally compute the formal solution man:
    I[z] = I[z-dz] * etb + w0 * s0 + w1 * s1;
    // There are two which are, without a doubt local:
    response_to_em[z][z] = (w1)/oo[z];
    response_to_em[z][z-dz] = (w0)/oo[z-dz];
    // Then some which are non local
    for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
      response_to_em[z][zp] += response_to_em[z-dz][zp] * etb;
    }
    // -----------------------------------------------------------------------------------------------
    // The story starts the same with the opacity:
    response_to_op[z][z] = -(w1)*ee[z]/oo[z]/oo[z];
    response_to_op[z][z-dz] = -(w0)*ee[z-dz]/oo[z-dz]/oo[z-dz];
    // Then some which are non local
    for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
      response_to_op[z][zp] += response_to_op[z-dz][zp] * etb;
    }
    // Then there is also dependance on delta:
    for (int zp=zl;(zp-zh-dz)*dz<0;zp+=dz){
      response_to_op[z][zp] += -I[z-dz] * etb * d_delta_per_d_chi[z][zp];
      response_to_op[z][zp] += (s0 * d_w0  + s1*d_w1) * d_delta_per_d_chi[z][zp];  
    }
    // And compute the local operator:
    if(L) l[z] = w1;
  }

  /*for (int z=x3l;z<=x3h;++z)
    for (int zz=x3l;zz<=x3h;++zz){
      if (isnan(response_to_op[z][zz])){
        printf("Error: Response of local intensity to opacity for l = %d and ll = %d is NaN. Fix! \n", z, zz);
        printf("%e %e %e %e \n", d_der_per_d_chi[z][zz],d_delta_per_d_chi[z][zz], d_sp_d_delta[z][zz], d_sp_d_s[z][zz]) ;
        exit(1);

      }
      if (isnan(response_to_em[z][zz])){
        printf("Error: Response of local intensity to emissivity for l = %d and ll = %d is NaN. Fix! \n", z, zz);
        exit(1);
      }
  }*/
      
  delete [](delta_tau+x3l);
  delete [](opp_derivative+x3l);
  delete [](s_derivative+x3l);
  del_ft2dim(d_delta_per_d_chi,x3l,x3h,x3l,x3h);
  del_ft2dim(d_sp_d_delta,x3l,x3h,x3l,x3h);
  del_ft2dim(d_sp_d_s,x3l,x3h,x3l,x3h);
  del_ft2dim(d_der_per_d_chi,x3l,x3h,x3l,x3h);

  return 0;
}


int atmos_ppbez::formal_pert_numerical(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t ** lambda_full = ft2dim(x3l,x3h,x3l,x3h);
  fp_t *** Intensity = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t *** alo_temp = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    
  // Then we perform a formal solution:
  formal(x3, Intensity, alo_temp, op, em, theta, phi, boundary);
  memset(dS[x3l][x1l][x2l]+x3l,0,(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  // At the moment we do it in the numerical way:
  
  fp_t *** Intensity_perturbed = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t *** opacity_perturbed = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t *** emissivity_perturbed = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  
  
  for (int x3k=x3l;x3k<=x3h;++x3k){

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          opacity_perturbed[x1i][x2i][x3i] = 0.5*op_pert[x3k][x1i][x2i][x3i]+op[x1i][x2i][x3i];
          emissivity_perturbed[x1i][x2i][x3i] = 0.5*em_pert[x3k][x1i][x2i][x3i]+em[x1i][x2i][x3i];
        }
    formal(x3, Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          dS[x3k][x1i][x2i][x3i] = Intensity_perturbed[x1i][x2i][x3i];

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          opacity_perturbed[x1i][x2i][x3i] = -0.5*op_pert[x3k][x1i][x2i][x3i]+op[x1i][x2i][x3i];
          emissivity_perturbed[x1i][x2i][x3i] = -0.5*em_pert[x3k][x1i][x2i][x3i]+em[x1i][x2i][x3i];
        }
    formal(x3, Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          dS[x3k][x1i][x2i][x3i] -= Intensity_perturbed[x1i][x2i][x3i];
    

  }
  del_ft3dim(Intensity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(opacity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(emissivity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h);
  

  del_ft2dim(lambda_full, x3l,x3h,x3l,x3h);
  del_ft3dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(alo_temp, x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
  
}

int atmos_ppbez::formal_pert_numerical_taugrid(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, 
  fp_t*** op_ref, fp_t**** op_referent_pert, fp_t theta, fp_t phi, fp_t boundary){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t ** lambda_full = ft2dim(x3l,x3h,x3l,x3h);
  fp_t *** Intensity = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t *** alo_temp = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    
  // Then we perform a formal solution:
  formal(tau_referent[x1l][x2l], Intensity, alo_temp, op, em, theta, phi, boundary);
  memset(dS[x3l][x1l][x2l]+x3l,0,(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  // We need, in addition the referent intensity derivatives.
  
  // At the moment we do it in the numerical way:
  
  fp_t *** Intensity_perturbed = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t *** opacity_perturbed = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t *** emissivity_perturbed = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);

    
  for (int x3k=x3l;x3k<=x3h;++x3k){

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          opacity_perturbed[x1i][x2i][x3i] = 0.5*op_pert[x3k][x1i][x2i][x3i]+op[x1i][x2i][x3i];
          emissivity_perturbed[x1i][x2i][x3i] = 0.5*em_pert[x3k][x1i][x2i][x3i]+em[x1i][x2i][x3i];

          opacity_perturbed[x1i][x2i][x3i] /= (op_ref[x1i][x2i][x3i] + 0.5 * op_referent_pert[x3k][x1i][x2i][x3i]);
          emissivity_perturbed[x1i][x2i][x3i] /= (op_ref[x1i][x2i][x3i]+ 0.5 * op_referent_pert[x3k][x1i][x2i][x3i]);
        }
    formal(tau_referent[x1l][x2l], Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          dS[x3k][x1i][x2i][x3i] = Intensity_perturbed[x1i][x2i][x3i];

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          opacity_perturbed[x1i][x2i][x3i] = -0.5*op_pert[x3k][x1i][x2i][x3i]+op[x1i][x2i][x3i];
          emissivity_perturbed[x1i][x2i][x3i] = -0.5*em_pert[x3k][x1i][x2i][x3i]+em[x1i][x2i][x3i];

          opacity_perturbed[x1i][x2i][x3i] /= (op_ref[x1i][x2i][x3i] - 0.5 * op_referent_pert[x3k][x1i][x2i][x3i]);
          emissivity_perturbed[x1i][x2i][x3i] /= (op_ref[x1i][x2i][x3i] - 0.5 * op_referent_pert[x3k][x1i][x2i][x3i]);
        }
    formal(tau_referent[x1l][x2l], Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          dS[x3k][x1i][x2i][x3i] -= Intensity_perturbed[x1i][x2i][x3i];
    

  }
  del_ft3dim(Intensity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(opacity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(emissivity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h);
  

  del_ft2dim(lambda_full, x3l,x3h,x3l,x3h);
  del_ft3dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(alo_temp, x1l,x1h,x2l,x2h,x3l,x3h);


  return 0;
  
}

int atmos_ppbez::formal_pert_analytical(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t ** response_to_op = ft2dim(x3l,x3h,x3l,x3h);
  fp_t ** response_to_em = ft2dim(x3l,x3h,x3l,x3h);
  fp_t *** Intensity = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    
  // Then we perform a formal solution:
  memset(dS[x3l][x1l][x2l]+x3l,0,(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  formal_with_responses(x3, Intensity, 0, response_to_op, response_to_em, op, em, theta, phi, boundary);
 
  // For each perturbation:
  for (int x3k=x3l;x3k<=x3h;++x3k){

    // Perturbation of intensity at each point:
    for (int x3i=x3l;x3i<=x3h;++x3i)
      for (int x3ii=x3l;x3ii<=x3h;++x3ii){
        dS[x3k][x1l][x2l][x3i] += response_to_em[x3i][x3ii] * em_pert[x3k][x1l][x2l][x3ii] + response_to_op[x3i][x3ii] * op_pert[x3k][x1l][x2l][x3ii];
      }

  }

  del_ft2dim(response_to_op, x3l,x3h,x3l,x3h);
  del_ft2dim(response_to_em, x3l,x3h,x3l,x3h);
  del_ft3dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
  
}

int atmos_ppbez::formal_pert_analytical_taugrid(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t ** response_to_op = ft2dim(x3l,x3h,x3l,x3h);
  fp_t ** response_to_em = ft2dim(x3l,x3h,x3l,x3h);
  fp_t *** Intensity = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);

  memset(dS[x3l][x1l][x2l]+x3l,0,(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
    
  // Then we perform a formal solution:
  formal_with_responses(tau_referent[x1l][x2l], Intensity, 0, response_to_op, response_to_em, op, em, theta, phi, boundary);

  /*fp_t Intensity_ref = 0;
  for (int x3i=x3l;x3i<=x3h;++x3i){
    Intensity_ref = 0.0;
    for (int x3ii=x3l;x3ii<=x3h;++x3ii)
      Intensity_ref += response_to_em[x3i][x3ii] * em[x1l][x2l][x3ii];
    printf("%e %e \n", Intensity[x1l][x2l][x3i], Intensity_ref);
  }

  printf("Blingling\n");
  exit(1);*/
 
  // For each perturbation:
  for (int x3k=x3l;x3k<=x3h;++x3k){

    // Perturbation of intensity at each point:
    for (int x3i=x3l;x3i<=x3h;++x3i)
      for (int x3ii=x3l;x3ii<=x3h;++x3ii){
        dS[x3k][x1l][x2l][x3i] += response_to_em[x3i][x3ii] * em_pert[x3k][x1l][x2l][x3ii] + response_to_op[x3i][x3ii] * op_pert[x3k][x1l][x2l][x3ii];
      }

  }

  del_ft2dim(response_to_op, x3l,x3h,x3l,x3h);
  del_ft2dim(response_to_em, x3l,x3h,x3l,x3h);
  del_ft3dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
  
}

int atmos_ppbez::formal_pert_jcdti(fp_t **** dS, fp_t *** op, fp_t *** em, fp_t **** op_pert, fp_t **** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t ** response_to_op = ft2dim(x3l,x3h,x3l,x3h);
  fp_t ** response_to_em = ft2dim(x3l,x3h,x3l,x3h);
  fp_t *** Intensity = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t **** em_for_pert = ft4dim(x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
    
  // Then we perform a formal solution:
  memset(dS[x3l][x1l][x2l]+x3l,0,(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  formal_with_responses(x3, Intensity, 0, response_to_op, response_to_em, op, em, theta, phi, boundary);
 
  // For each perturbation:

  for (int x3k=x3l;x3k<=x3h;++x3k){

    // Perturbation of intensity at each point:
    // Compute emissivity for perturbation at each point:
    for (int x3i=x3l;x3i<=x3h;++x3i)
      em_for_pert[x3k][x1l][x2l][x3i] = em_pert[x3k][x1l][x2l][x3i] - op_pert[x3k][x1l][x2l][x3i] * Intensity[x1l][x2l][x3i];
    
    for (int x3i=x3l;x3i<=x3h;++x3i)
      for (int x3ii=x3l;x3ii<=x3h;++x3ii){
        dS[x3k][x1l][x2l][x3i] += response_to_em[x3i][x3ii] * em_for_pert[x3k][x1l][x2l][x3ii];
      
    }

  }

  del_ft2dim(response_to_op, x3l,x3h,x3l,x3h);
  del_ft2dim(response_to_em, x3l,x3h,x3l,x3h);
  del_ft3dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft4dim(em_for_pert,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
}


fp_t * atmos_ppbez::compute_full_lambda_operator(fp_t *** tau, int l, int dz){ // Non existant 

  return 0;    
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 

int atmos_ppbez::formal(fp_t * spatial_grid, fp_t ****S,fp_t ***L,fp_t *****op,fp_t ****em,fp_t t,fp_t p, fp_t boundary) // same as above but for polarized case
{
  // Input quantities:
  // **** S = input Stokes intensity vector in which eventually the result of the computation is going to written
  // *** L = approximate lambda operator which is being computed in the same go, most of the times this will be zero, as we rarely iterate on polarized quantities
  // ***** op = Opacity, for given freq and direction, it is a 5D quantity as it is a 3D array of 4x4 matrices
  // **** em = Emissivity, for given freq and direction, it is a 4D quantity
  // t      = theta angle
  // p      = phi angle, relevant only in 3D case, but ok we can use it here as well in case of Zeeman or whatever
  // boundary = this specifies what kind of boundary condition are we using. Anything greater then zero means that unpolarized intensity is specified and equal to boundary. 
  //            if it is -1 then boundary is equal to local source function. 

  fp_t **I = S[x1l][x2l]; // "Current" value of intensity. I.e. the one for this "column" 
  fp_t ***oo = op[x1l][x2l]; // Same for the opacity
  fp_t **ee = em[x1l][x2l]; // Same for emissivity
  fp_t *l = (L) ? L[x1l][x2l] : 0; // Now, if we pass zero, then do not compute the local operator, if we pass something we, surprise-surprise, do!
// 
  fp_t ct=cos(t); // Cosinus of the direction. If you are for first time looking @ this, do not be confused, the directions are the opposite: 
                  // t = 0 indicates ingoing ray while t = pi indicates the outgoing ray, i.e. angle is with respect to the direction
                  // observer - center of the star.
  
  int dz = (ct>0) ? 1 : -1; // This is increment along z which determines what do we do, i.e. if we increase or decrease the z
  int zl = (ct>0) ? x3l : x3h; // Where do we start from.
  int zh = (ct>0) ? x3h : x3l; // Where do we end.
//
  // Boundary condition for the intensity.
  I[zl][1]= (ct > 0) ? 0.0 : ee[zl][1] / oo[zl][1][1];

  if(ct < 0 && boundary >= 0) // lower boundary condition
      I[zl][1] = boundary;

  I[zl][2] = I[zl][3] = I[zl][4] = 0.0; // Is this correct

  // And appropriate local operator.
  if(L && ct >= 0) l[zl] = 0.0; // Actually here is something that we might want to examine further; i.e. what happens with the local operator
                               // at the lower boundary condition, what makes sense is to set 1.0 at the lower one, no?
  if(L && ct <= 0) {

    if (boundary == -1)
      l[zl] = 1.0;
    else 
      l[zl] = 0.0;
  }

  // We will now compute the optical depth scale from here.
  // First derivatives of the opacity.
  int ND = x3h-x3l+1;
  fp_t * delta_tau = new fp_t [ND] - x3l;
  fp_t * opp_derivative = new fp_t [ND] - x3l;  
  opp_derivative[zl] = (oo[zl+dz][1][1] - oo[zl][1][1]) / (spatial_grid[zl+dz] - spatial_grid[zl]) * (-ct);
  opp_derivative[zh] = (oo[zh][1][1] - oo[zh-dz][1][1]) / (spatial_grid[zh] - spatial_grid[zh-dz]) * (-ct);
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    opp_derivative[z] = 0.0; // default value
    fp_t h2 = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct); 
    fp_t h1 = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
    fp_t d2 = (oo[z+dz][1][1] - oo[z][1][1]) / h2;
    fp_t d1 = (oo[z][1][1] - oo[z-dz][1][1]) / h1;
    if (d1*d2 > 0 || disable_monotonicity){
      fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
      opp_derivative[z] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
    }
  }
  // Then from there the optical depth scale:
  delta_tau[zl] = 0.0; // By default.
  for (int z = zl+dz;(z-zh-dz)*dz<0;z+=dz){
    fp_t C0 = oo[z][1][1] - (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z];
    fp_t C1 = oo[z-dz][1][1] + (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z-dz];
    fp_t C = 0.5*(C0+C1);
    fp_t delta_t = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 3.0 * (oo[z-dz][1][1] + oo[z][1][1] + C);
    delta_tau[z] = delta_t;
  }

  // Compute the derivatives of the source function. (Following Piskunov & Jaime, 2013, this is j / \eta_I)
  fp_t ** s_derivative = ft2dim(x3l,x3h,1,4);
  fp_t ** source = ft2dim(x3l,x3h,1,4);
  for (int s=1;s<=4;++s){ // For each Stokes parameter
    s_derivative[zl][s] = (ee[zl+dz][s]/oo[zl+dz][1][1] - ee[zl][s]/oo[zl][1][1]) / delta_tau[zl+dz];
    s_derivative[zh][s] = (ee[zh][s]/oo[zh][1][1] - ee[zh-dz][s]/oo[zh-dz][1][1]) / delta_tau[zh];
  }
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    for (int s=1;s<=4;++s){
      s_derivative[z][s] = 0.0; // default value
      fp_t h2 = delta_tau[z+dz];
      fp_t h1 = delta_tau[z];
      fp_t d2 = (ee[z+dz][s]/oo[z+dz][1][1] - ee[z][s]/oo[z][1][1]) / h2;
      fp_t d1 = (ee[z][s]/oo[z][1][1] - ee[z-dz][s]/oo[z-dz][1][1]) / h1;
      if(d1*d2>0 || disable_monotonicity){
        fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
        s_derivative[z][s] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
      }
    }
  }

  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int s=1;s<=4;++s)
    source[x3i][s] = ee[x3i][s] / oo[x3i][1][1];
    

  // Finally, we need to compute the derivatives of the reduced absorption matrix, again from Piskunov & de la Cruz Rodriguez, 2013:

  fp_t *** K_reduced = ft3dim(x3l,x3h,1,4,1,4);
  fp_t *** K_reduced_derivative = ft3dim(x3l,x3h,1,4,1,4);

  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp)
        K_reduced[x3i][s][sp] = oo[x3i][s][sp] / oo[x3i][1][1];
      K_reduced[x3i][s][s] -= 1.0;
  }

  // And then you have to compute the derivative, in the same way as for the source function:

  for (int s=1;s<=4;++s)
    for (int sp=1;sp<=4;++sp){
      K_reduced_derivative[zl][s][sp] = (K_reduced[zl+dz][s][sp] - K_reduced[zl][s][sp]) / delta_tau[zl+dz];
      K_reduced_derivative[zh][s][sp] = (K_reduced[zh][s][sp] - K_reduced[zh-dz][s][sp]) / delta_tau[zh];
  }
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp){
        K_reduced_derivative[z][s][sp] = 0.0; // default value
        fp_t h2 = delta_tau[z+dz];
        fp_t h1 = delta_tau[z];
        fp_t d2 = (K_reduced[z+dz][s][sp] - K_reduced[z][s][sp]) / h2;
        fp_t d1 = (K_reduced[z][s][sp] - K_reduced[z-dz][s][sp]) / h1;
        if(d1*d2>0 || disable_monotonicity){
          fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
          K_reduced_derivative[z][s][sp] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
        }
      }
    }
  }

  fp_t * A_decomposed_by_LU = new fp_t [16];
  fp_t * I_solution = new fp_t [4];

  // Now we need to formally solve RTE, layer by layer, EXCEPT THE LAST LAYER, which we do at the end.
  for(int z = zl+dz; (z-zh)*dz < 0; z+=dz){ 
    
    fp_t tb = delta_tau[z];
    fp_t tf = delta_tau[z+dz];  
    // APPROXIMATE TO CHECK WHAT IS THE PROBLEM:
    if (tau_switch_debug){
      fp_t dlb = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
      fp_t dlf = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct);
      tb = dlb * (oo[z-dz][1][1] + oo[z][1][1]) * 0.5;
      tf = dlf * (oo[z][1][1] + oo[z+dz][1][1]) * 0.5;
    }
    
    // Now, when we have the upwind and dowind (backward and forward) optical paths, we can also compute the weights and control point. 
    // And, subsequently the contribution to the lambda operator.
    
    fp_t etb = (tb < expansion_delimiter) ? (1.0 - tb + tb*tb*0.5 - tb*tb*tb/6.0) : exp(-tb); // Exponent of the backward optical depth, will become useful.
    fp_t w0, w1, w2;       // Integration weights, for local, previous and control point.

    if (tf > bezier_delimiter){ // Bezier, ie. higher then 1st order
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb*(1.0/3.0 - tb*(1.0/4.0-tb/10.0));
        w1 = tb/3.0*(1.0-tb/4.0*(1.0-tb/5.0));
        w2 = tb*(1.0/3.0-tb*(1.0/6.0-tb/20.0));
      }
      else {
        w0 = (2.0 - etb * (tb * tb + 2.0 * tb + 2.0)) / tb / tb;
        w1 = 1.0 - 2.0 * (etb + tb - 1.0) / tb / tb; 
        w2 = 2.0 * (tb - 2.0  + etb * (tb + 2.0)) / tb / tb;
      }
    }
    else { // Linear
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb * (0.5 - tb * (1.0/3.0  - tb/8.0));
        w1 = 0.5 * tb * (1.0 - tb/3.0 * (1.0 - tb/4.0));
        w2 = 0.0;
      }
      else {
        w0 = 1.0 / tb - etb * (1.0 + 1.0 / tb);
        w1 = 1.0 - 1.0 / tb + etb / tb;
        w2 = 0.0;
      }
    }

    // Now this part is much different in the polarized case. 

    fp_t ** A = ft2dim(1,4,1,4);
    fp_t *  B = new fp_t [4]-1;

    // We need to fill in A and B, use equations 24 - 27 from page 3 of Jaimes paper. NOTE: For our notation, we need to change sign in front of
    // delta/2

    fp_t * C1k; // This is the thing in equation 25

    fp_t ** temp_matrix = multiply_square(K_reduced[z-dz], K_reduced[z-dz], 4);
    for (int s=1;s<=4;++s)
      for (int sp=1;sp<=4;++sp){
        temp_matrix[s][sp] += K_reduced[z-dz][s][sp] + K_reduced_derivative[z-dz][s][sp];
        temp_matrix[s][sp] *= tb/2.0;
        temp_matrix[s][sp] += K_reduced[z-dz][s][sp]; // QUESTIONABLE. WHY DID WE MISS THIS THE FIRST TIME?
      }
    fp_t * temp_vector;
    temp_vector = multiply_vector(temp_matrix, I[z-dz], 4);

    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp)
        temp_matrix[s][sp] = tb/2.0 * K_reduced[z-dz][s][sp];
      temp_matrix[s][s] += 1.0;
    }

    C1k = multiply_vector(temp_matrix, source[z-dz], 4);
    for (int s=1;s<=4;++s)
      C1k[s] += tb/2.0 * s_derivative[z-dz][s] - temp_vector[s];

    del_ft2dim(temp_matrix,1,4,1,4);
    delete [](temp_vector+1);

    // Done with C1k 

    fp_t * C0k; // Now equation (24)

    temp_matrix = ft2dim(1,4,1,4);
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp)
        temp_matrix[s][sp] = -tb/2.0 * K_reduced[z][s][sp];
      temp_matrix[s][s] += 1.0;
    }

    C0k = multiply_vector(temp_matrix, source[z],4);
    for (int s=1;s<=4;++s)
      C0k[s] -= tb/2.0 * s_derivative[z][s];

    del_ft2dim(temp_matrix,1,4,1,4);

    fp_t ** c0k; // This is the other part of equation (24)
    c0k = multiply_square(K_reduced[z],K_reduced[z],4);
    for (int s=1;s<=4;++s)
      for (int sp=1;sp<=4;++sp){
        c0k[s][sp] += K_reduced[z][s][sp] + K_reduced_derivative[z][s][sp];
        c0k[s][sp] *= tb/2.0;
        c0k[s][sp] -= K_reduced[z][s][sp];
      }

   // A is easy:
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp)
        A[s][sp] = w1 * K_reduced[z][s][sp] - w2/2.0 * c0k[s][sp];
      A[s][s] += 1.0;
    }

    // Now B:

    temp_vector = multiply_vector(K_reduced[z-dz], I[z-dz],4);
    for (int s=1;s<=4;++s)
      B[s] = -w0 * temp_vector[s] + w0 * source[z-dz][s] + w1 * source[z][s]  + I[z-dz][s] * etb + w2/2.0 * (C1k[s] + C0k[s]);


    delete [](temp_vector+1);
    delete [](C0k+1);
    delete [](C1k+1);
    del_ft2dim(c0k,1,4,1,4);

    
    // After this we solve:
    fp_t * A_to_solve = A[1]+1;
    fp_t * b = B+1;
    Crout(4,A_to_solve, A_decomposed_by_LU);
    solveCrout(4,A_decomposed_by_LU,b,I_solution);
    memcpy(I[z]+1,I_solution,4*sizeof(fp_t));

    /*for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp)
        printf("%e ", A[s][sp]);
      printf("   %e %e\n", B[s], I[z][s]);
    }
    printf("\n");*/
     
    // And compute the local operator
    if(L) l[z] = (w1+0.5*w2);

    del_ft2dim(A,1,4,1,4);
    delete [](B+1);
  }

  // Still the final point in the atmosphere is not covered here. That one we will compute in the linear approximation, both for 
  // the opacity and for the source function.
  if (zl != zh){
    int z = zh; // Just to be absolutely sure as I am a nub programmer sometimes (Milic)
    for (int s=1;s<=4;++s)
      I[z][s] = I[z-dz][s];
  }
    
  delete [](delta_tau+x3l);
  delete [](opp_derivative+x3l);
  del_ft2dim(source,x3l,x3h,1,4);
  del_ft2dim(s_derivative,x3l,x3h,1,4);
  del_ft3dim(K_reduced,x3l,x3h,1,4,1,4);
  del_ft3dim(K_reduced_derivative,x3l,x3h,1,4,1,4);
  delete []A_decomposed_by_LU;
  delete []I_solution;

  return 0;

}

int atmos_ppbez::formal_with_responses_jcdti(fp_t * spatial_grid, fp_t **** S,fp_t *** L,fp_t ***** op,fp_t **** em, 
  fp_t ***** op_pert, fp_t **** em_pert, fp_t **** S_pert, fp_t t,fp_t phi, fp_t boundary){return 0;}

int atmos_ppbez::formal_with_responses_full(fp_t * spatial_grid, fp_t ****S,fp_t ***L,fp_t *****op,fp_t ****em, 
fp_t ***** op_pert, fp_t **** em_pert, fp_t **** S_pert,fp_t t,fp_t p, fp_t boundary)
{
  // This is the function which follows the formal solution and then propagates the perturbations in op and em in order to compute emergent spectrum
  // Input quantities:
  // **** S = input Stokes intensity vector in which eventually the result of the computation is going to written
  // *** L = approximate lambda operator which is being computed in the same go, most of the times this will be zero, as we rarely iterate on polarized quantities
  // ***** op = Opacity, for given freq and direction, it is a 5D quantity as it is a 3D array of 4x4 matrices
  // **** em = Emissivity, for given freq and direction, it is a 4D quantity
  // ***** op_pert = perturbation (derivative) of the opacity
  // **** em_pert = perturbation (derivative) of the emissivity
  // **** I_pert = perturbation (derivative) of the intensity, everywhere. This is redundant as we ultimately only need the top one. But let it be.
  // t      = theta angle
  // p      = phi angle, relevant only in 3D case, but ok we can use it here as well in case of Zeeman or whatever
  // boundary = this specifies what kind of boundary condition are we using. Anything greater then zero means that unpolarized intensity is specified and equal to boundary. 
  //            if it is -1 then boundary is equal to local source function. 

  fp_t **I = S[x1l][x2l]; // "Current" value of intensity. I.e. the one for this "column" 
  fp_t **I_pert = S_pert[x1l][x2l];
  fp_t ***oo = op[x1l][x2l]; // Same for the opacity
  fp_t **ee = em[x1l][x2l]; // Same for emissivity
  fp_t *** oo_pert = op_pert[x1l][x2l];
  fp_t **ee_pert = em_pert[x1l][x2l];
  fp_t *l = (L) ? L[x1l][x2l] : 0; // Now, if we pass zero, then do not compute the local operator, if we pass something we, surprise-surprise, do!
// 
  fp_t ct=cos(t); // Cosinus of the direction. If you are for first time looking @ this, do not be confused, the directions are the opposite: 
                  // t = 0 indicates ingoing ray while t = pi indicates the outgoing ray, i.e. angle is with respect to the direction
                  // observer - center of the star.
  
  int dz = (ct>0) ? 1 : -1; // This is increment along z which determines what do we do, i.e. if we increase or decrease the z
  int zl = (ct>0) ? x3l : x3h; // Where do we start from.
  int zh = (ct>0) ? x3h : x3l; // Where do we end.
//
  // Boundary condition for the intensity.
  I[zl][1] = (ct > 0) ? 0.0 : ee[zl][1] / oo[zl][1][1];
  I_pert[zl][1] = (ct>0) ? 0.0 : (ee_pert[zl][1] * oo[zl][1][1] - oo_pert[zl][1][1] * ee[zl][1])/oo[zl][1][1]/oo[zl][1][1];

  if(ct < 0 && boundary >= 0){ // lower boundary condition
    I[zl][1] = boundary;
    I_pert[zl][1] = 0.0;;
  }

  I[zl][2] = I[zl][3] = I[zl][4] = 0.0; // Is this correct
  I_pert[zl][2] = I_pert[zl][3] = I_pert[zl][4] = 0.0;//

  // And appropriate local operator.
  if(L && ct >= 0) l[zl] = 0.0; // Actually here is something that we might want to examine further; i.e. what happens with the local operator
                               // at the lower boundary condition, what makes sense is to set 1.0 at the lower one, no?
  if(L && ct <= 0) {

    if (boundary == -1)
      l[zl] = 1.0;
    else 
      l[zl] = 0.0;
  }

  // We will now compute the optical depth scale from here.
  // First derivatives of the opacity.
  int ND = x3h-x3l+1;
  fp_t * delta_tau = new fp_t [ND] - x3l;
  fp_t * opp_derivative = new fp_t [ND] - x3l;
  fp_t * delta_tau_pert = new fp_t [ND] - x3l;
  fp_t * opp_derivative_pert = new fp_t [ND] - x3l;  
  
  opp_derivative[zl] = (oo[zl+dz][1][1] - oo[zl][1][1]) / (spatial_grid[zl+dz] - spatial_grid[zl]) * (-ct);
  opp_derivative[zh] = (oo[zh][1][1] - oo[zh-dz][1][1]) / (spatial_grid[zh] - spatial_grid[zh-dz]) * (-ct);
  opp_derivative_pert[zl] = (oo_pert[zl+dz][1][1] - oo_pert[zl][1][1]) / (spatial_grid[zl+dz] - spatial_grid[zl]) * (-ct);
  opp_derivative_pert[zh] = (oo_pert[zh][1][1] - oo_pert[zh-dz][1][1]) / (spatial_grid[zh] - spatial_grid[zh-dz]) * (-ct);
  
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    opp_derivative[z] = 0.0; // default value
    opp_derivative_pert[z] = 0.0; // default value
    fp_t h2 = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct); 
    fp_t h1 = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
    fp_t d2 = (oo[z+dz][1][1] - oo[z][1][1]) / h2;
    fp_t d1 = (oo[z][1][1] - oo[z-dz][1][1]) / h1;
    fp_t d2_pert = (oo_pert[z+dz][1][1] - oo_pert[z][1][1]) / h2;
    fp_t d1_pert = (oo_pert[z][1][1] - oo_pert[z-dz][1][1]) / h1;
    if (d1*d2 > 0 || disable_monotonicity){
      fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
      opp_derivative[z] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
      opp_derivative_pert[z] = opp_derivative[z] * (d1_pert/d1 + d2_pert/d2) - opp_derivative[z] * (alpha*d1_pert + (1.0-alpha) * d2_pert) / (alpha*d1 + (1.0-alpha)*d2);
    }
  }
  // Then from there the optical depth scale:
  delta_tau[zl] = 0.0; // By default.
  delta_tau_pert[zl] = 0.0;
  for (int z = zl+dz;(z-zh-dz)*dz<0;z+=dz){
    fp_t C0 = oo[z][1][1] - (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z];
    fp_t C1 = oo[z-dz][1][1] + (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative[z-dz];
    fp_t C = 0.5*(C0+C1);
    fp_t delta_t = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 3.0 * (oo[z-dz][1][1] + oo[z][1][1] + C);
    delta_tau[z] = delta_t;

    // Same for perturbation:
    fp_t C0_pert = oo_pert[z][1][1] - (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative_pert[z];
    fp_t C1_pert = oo_pert[z-dz][1][1] + (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 2.0 * opp_derivative_pert[z-dz];
    fp_t C_pert = 0.5*(C0_pert+C1_pert);
    fp_t delta_t_pert = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct) / 3.0 * (oo_pert[z-dz][1][1] + oo_pert[z][1][1] + C_pert);
    delta_tau_pert[z] = delta_t_pert;

  }

  // Compute the derivatives of the source function. (Following Piskunov & Jaime, 2013, this is j / \eta_I)
  fp_t ** s_derivative = ft2dim(x3l,x3h,1,4);
  fp_t ** source = ft2dim(x3l,x3h,1,4);
  
  fp_t ** s_derivative_pert = ft2dim(x3l,x3h,1,4);
  fp_t ** source_pert = ft2dim(x3l,x3h,1,4);

  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int s=1;s<=4;++s){
      source[x3i][s] = ee[x3i][s] / oo[x3i][1][1];
      source_pert[x3i][s] = (ee_pert[x3i][s] * oo[x3i][1][1] - oo_pert[x3i][1][1]*ee[x3i][s])/ oo[x3i][1][1]/ oo[x3i][1][1];
    }


  for (int s=1;s<=4;++s){ // For each Stokes parameter
    s_derivative[zl][s] = (source[zl+dz][s] - source[zl][s]) / delta_tau[zl+dz];
    s_derivative[zh][s] = (source[zh][s] - source[zh-dz][s]) / delta_tau[zh];
    s_derivative_pert[zl][s] = (source_pert[zl+dz][s] - source_pert[zl][s]) / delta_tau[zl+dz] - (source[zl+dz][s] - source[zl][s])*delta_tau_pert[zl+dz] / delta_tau[zl+dz]/delta_tau[zl+dz];
    s_derivative_pert[zh][s] = (source_pert[zh][s] - source_pert[zh-dz][s]) / delta_tau[zh] - (source[zh][s] - source[zh-dz][s])*delta_tau_pert[zh] / delta_tau[zh] / delta_tau[zh];
  }
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    for (int s=1;s<=4;++s){
      s_derivative[z][s] = 0.0; // default value
      s_derivative_pert[z][s] = 0.0;
      fp_t h2 = delta_tau[z+dz];
      fp_t h1 = delta_tau[z];
      fp_t h2_pert = delta_tau_pert[z+dz];
      fp_t h1_pert = delta_tau_pert[z];
      fp_t d2 = (source[z+dz][s] - source[z][s]) / h2;
      fp_t d1 = (source[z][s] - source[z-dz][s]) / h1;
      fp_t d2_pert = (source_pert[z+dz][s] - source_pert[z][s]) / h2 - (source[z+dz][s] - source[z][s])* h2_pert/h2/h2;
      fp_t d1_pert = (source_pert[z][s] - source_pert[z-dz][s]) / h1 - (source[z][s] - source[z-dz][s])* h1_pert/h1/h1;
      //if(s==1) printf("> %d %e %e %e %e \n",z, source_pert[z+dz][s],source_pert[z][s], h2,h2_pert);

      if(d1*d2>0 || disable_monotonicity){
        fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
        fp_t d_alpha = (h1_pert/(h1+h2) - h1*(h1_pert+h2_pert)/(h1+h2)/(h1+h2))/3.0;
        s_derivative[z][s] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
        s_derivative_pert[z][s] = s_derivative[z][s] * (d1_pert/d1 + d2_pert/d2) - s_derivative[z][s] * (alpha*d1_pert + (1.0-alpha)*d2_pert)/(alpha*d1 + (1.0-alpha)*d2);
        s_derivative_pert[z][s] -= s_derivative[z][s]/(alpha*d1 + (1.0-alpha)*d2)*(d1-d2)*d_alpha;
        //if(s==1) printf(">> %d %e %e %e %e \n",z, s_derivative[z][s], d1_pert, d2_pert, d_alpha);

      }
    }
    //printf("%d %e \n", z, s_derivative_pert[z][1]);
  }  
  
  // Finally, we need to compute the derivatives of the reduced absorption matrix, again from Piskunov & de la Cruz Rodriguez, 2013:

  fp_t *** K_reduced = ft3dim(x3l,x3h,1,4,1,4);
  fp_t *** K_reduced_derivative = ft3dim(x3l,x3h,1,4,1,4);

  fp_t *** K_reduced_pert = ft3dim(x3l,x3h,1,4,1,4);
  fp_t *** K_reduced_derivative_pert = ft3dim(x3l,x3h,1,4,1,4);

  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp){
        K_reduced[x3i][s][sp] = oo[x3i][s][sp]/oo[x3i][1][1];
        K_reduced_pert[x3i][s][sp] = oo_pert[x3i][s][sp]/oo[x3i][1][1] - oo[x3i][s][sp]*oo_pert[x3i][1][1]/oo[x3i][1][1]/oo[x3i][1][1];
      }
      K_reduced[x3i][s][s] -= 1.0;
  }

  // And then you have to compute the derivative, in the same way as for the source function:

  for (int s=1;s<=4;++s)
    for (int sp=1;sp<=4;++sp){
      K_reduced_derivative[zl][s][sp] = (K_reduced[zl+dz][s][sp] - K_reduced[zl][s][sp]) / delta_tau[zl+dz];
      K_reduced_derivative[zh][s][sp] = (K_reduced[zh][s][sp] - K_reduced[zh-dz][s][sp]) / delta_tau[zh];

      K_reduced_derivative_pert[zl][s][sp] = (K_reduced_pert[zl+dz][s][sp] - K_reduced_pert[zl][s][sp]) / delta_tau[zl+dz] -
        (K_reduced[zl+dz][s][sp] - K_reduced[zl][s][sp]) * delta_tau_pert[zl+dz] / delta_tau[zl+dz] / delta_tau[zl+dz];
      K_reduced_derivative_pert[zh][s][sp] = (K_reduced_pert[zh][s][sp] - K_reduced_pert[zh-dz][s][sp]) / delta_tau[zh] - 
        (K_reduced[zh][s][sp] - K_reduced[zh-dz][s][sp]) * delta_tau_pert[zh] / delta_tau[zh] / delta_tau[zh];
  }
  for (int z = zl+dz;(z-zh)*dz<0;z+=dz){
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp){
        K_reduced_derivative[z][s][sp] = 0.0; // default value
        K_reduced_derivative_pert[z][s][sp] = 0.0;
        fp_t h2 = delta_tau[z+dz];
        fp_t h1 = delta_tau[z];
        fp_t h2_pert = delta_tau_pert[z+dz];
        fp_t h1_pert = delta_tau_pert[z];
        fp_t d2 = (K_reduced[z+dz][s][sp] - K_reduced[z][s][sp]) / h2;
        fp_t d1 = (K_reduced[z][s][sp] - K_reduced[z-dz][s][sp]) / h1;
        fp_t d2_pert = (K_reduced_pert[z+dz][s][sp] - K_reduced_pert[z][s][sp]) / h2 - (K_reduced[z+dz][s][sp] - K_reduced[z][s][sp]) / h2 / h2 * h2_pert;
        fp_t d1_pert = (K_reduced_pert[z][s][sp] - K_reduced_pert[z-dz][s][sp]) / h1 - (K_reduced[z][s][sp] - K_reduced[z-dz][s][sp]) / h1 / h1 * h1_pert;
        if(d1*d2>0 || disable_monotonicity){
          fp_t alpha = (1.0 + h1/(h1+h2))/3.0;
          fp_t d_alpha = (h1_pert/(h1+h2) - h1*(h1_pert+h2_pert)/(h1+h2)/(h1+h2))/3.0;
          K_reduced_derivative[z][s][sp] = d1*d2 / (alpha*d1 + (1.0-alpha)*d2);
          K_reduced_derivative_pert[z][s][sp] = K_reduced_derivative[z][s][sp] * (d1_pert/d1 + d2_pert/d2) - K_reduced_derivative[z][s][sp] * (alpha*d1_pert + (1.0-alpha)*d2_pert)/(alpha*d1 + (1.0-alpha)*d2);
          K_reduced_derivative_pert[z][s][sp] -= K_reduced_derivative[z][s][sp]/(alpha*d1 + (1.0-alpha)*d2)*(d1-d2)*d_alpha;
        }
      }
    }
  }

  fp_t * A_decomposed_by_LU = new fp_t [16];
  fp_t * I_solution = new fp_t [4];

  // Now we need to formally solve RTE, layer by layer, EXCEPT THE LAST LAYER, which we do at the end.
  for(int z = zl+dz; (z-zh)*dz < 0; z+=dz){ 
    
    fp_t tb = delta_tau[z];
    fp_t tb_pert = delta_tau_pert[z];
    fp_t tf = delta_tau[z+dz];  
    // APPROXIMATE TO CHECK WHAT IS THE PROBLEM:
    if (tau_switch_debug){
      fp_t dlb = (spatial_grid[z] - spatial_grid[z-dz]) / (-ct);
      fp_t dlf = (spatial_grid[z+dz] - spatial_grid[z]) / (-ct);
      tb = dlb * (oo[z-dz][1][1] + oo[z][1][1]) * 0.5;
      tf = dlf * (oo[z][1][1] + oo[z+dz][1][1]) * 0.5;
    }
    
    // Now, when we have the upwind and dowind (backward and forward) optical paths, we can also compute the weights and control point. 
    // And, subsequently the contribution to the lambda operator.
    

    fp_t etb = (tb < expansion_delimiter) ? (1.0 - tb + tb*tb*0.5 - tb*tb*tb/6.0) : exp(-tb); // Exponent of the backward optical depth, will become useful.
    fp_t etb_pert = (tb < expansion_delimiter) ? (-tb_pert + tb*tb_pert - tb*tb*tb_pert/2.0) : -exp(-tb) * tb_pert;
    fp_t w0, w1, w2;       // Integration weights, for local, previous and control point.
    fp_t w0_pert, w1_pert, w2_pert;

    if (tf > bezier_delimiter){ // Bezier, ie. higher then 1st order
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb*(1.0/3.0 - tb*(1.0/4.0-tb/10.0));
        w0_pert = (1.0/3.0 - tb*(1.0/2.0 - 3.0*tb/10.0)) * tb_pert;
        w1 = tb/3.0*(1.0-tb/4.0*(1.0-tb/5.0));
        w1_pert = (1.0/3.0  - tb/2.0*(1.0/3.0 - tb/10.0)) * tb_pert;
        w2 = tb*(1.0/3.0-tb*(1.0/6.0-tb/20.0));
        w2_pert = (1.0/3.0 - tb/3.0 + 3.0*tb*tb/20.0) * tb_pert;
      }
      else {
        w0 = (2.0 - etb * (tb * tb + 2.0 * tb + 2.0)) / tb / tb;
        w0_pert = (-2.0*w0/tb + etb)*tb_pert;
        w1 = 1.0 - 2.0 * (etb + tb - 1.0) / tb / tb;
        w1_pert = (2.0 *(1.0-w1) / tb - 2.0*(1.0-etb)/tb/tb)*tb_pert; 
        w2 = 2.0 * (tb - 2.0  + etb * (tb + 2.0)) / tb / tb;
        w2_pert = (-2.0*w2/tb + 2.0/tb/tb * (1.0-tb*etb-etb))*tb_pert;
      }
      
    }
    else { // Linear
      if (tb < expansion_delimiter){ // Expansion
        w0 = tb * (0.5 - tb * (1.0/3.0  - tb/8.0));
        w0_pert = (1.0/2.0 - tb*(2.0/3.0 - 3.0*tb/8.0))*tb_pert;
        w1 = 0.5 * tb * (1.0 - tb/3.0 * (1.0 - tb/4.0));
        w1_pert = (0.5 - tb*(1.0/3.0 - tb/8.0))*tb_pert;
        w2 = 0.0;
        w2_pert = 0.0;
      }
      else {
        w0 = 1.0 / tb - etb * (1.0 + 1.0 / tb);
        w0_pert = (-w0/tb + etb)*tb_pert;
        w1 = 1.0 - 1.0 / tb + etb / tb;
        w1_pert = ((1.0-w1)/tb - etb/tb)*tb_pert;
        w2 = 0.0;
        w2_pert = 0.0;
      }
    }

    // Now this part is much different in the polarized case. 

    fp_t ** A = ft2dim(1,4,1,4);
    fp_t *  B = new fp_t [4]-1;
    fp_t ** A_pert = ft2dim(1,4,1,4);
    fp_t *  B_pert = new fp_t[4]-1;

    // We need to fill in A and B, use equations 24 - 27 from page 3 of Jaimes paper. NOTE: For our notation, we need to change sign in front of
    // delta/2.
    // Also, along the way we write down the perturbations for all these components. 

    fp_t * C1k; // This is the thing in equation 25
    fp_t * C1k_pert;

    fp_t ** temp_matrix = multiply_square(K_reduced[z-dz], K_reduced[z-dz], 4);
    fp_t ** temp_matrix_pert =ft2dim(1,4,1,4);
    memset(temp_matrix_pert[1]+1,0,16*sizeof(fp_t));
    for (int s=1;s<=4;++s)
      for (int sp=1;sp<=4;++sp)
        for (int kk=1;kk<=4;++kk)
          temp_matrix_pert[s][sp] += K_reduced[z-dz][s][kk] * K_reduced_pert[z-dz][kk][sp] + K_reduced_pert[z-dz][s][kk] * K_reduced[z-dz][kk][sp];

    for (int s=1;s<=4;++s)
      for (int sp=1;sp<=4;++sp){
        temp_matrix[s][sp] += K_reduced[z-dz][s][sp] + K_reduced_derivative[z-dz][s][sp];
        
        temp_matrix_pert[s][sp] += K_reduced_pert[z-dz][s][sp] + K_reduced_derivative_pert[z-dz][s][sp];
        temp_matrix_pert[s][sp] = temp_matrix_pert[s][sp] * tb/2.0 + tb_pert/2.0 * temp_matrix[s][sp];

        temp_matrix[s][sp] *= tb/2.0;
        temp_matrix[s][sp] += K_reduced[z-dz][s][sp]; // QUESTIONABLE. WHY DID WE MISS THIS THE FIRST TIME?
        temp_matrix_pert[s][sp] += K_reduced_pert[z-dz][s][sp]; // QUESTIONABLE. WHY DID WE MISS THIS THE FIRST TIME?

      }
    fp_t * temp_vector;
    fp_t * temp_vector_pert = new fp_t[4]-1;
    memset(temp_vector_pert+1,0,4*sizeof(fp_t));
    temp_vector = multiply_vector(temp_matrix, I[z-dz], 4);
    for(int s=1;s<=4;++s)
      for (int kk=1;kk<=4;++kk)
        temp_vector_pert[s] += temp_matrix[s][kk] * I_pert[z-dz][kk] + temp_matrix_pert[s][kk] * I[z-dz][kk];


    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp){
        temp_matrix[s][sp] = tb/2.0 * K_reduced[z-dz][s][sp];
        temp_matrix_pert[s][sp] = tb/2.0 * K_reduced_pert[z-dz][s][sp] + tb_pert/2.0 * K_reduced[z-dz][s][sp];
      }
      temp_matrix[s][s] += 1.0;
    }

    C1k = multiply_vector(temp_matrix, source[z-dz], 4);
    C1k_pert = new fp_t [4]-1;
    memset(C1k_pert+1,0,4*sizeof(fp_t));
    for(int s=1;s<=4;++s)
      for (int kk=1;kk<=4;++kk)
        C1k_pert[s] += temp_matrix[s][kk]*source_pert[z-dz][kk] + temp_matrix_pert[s][kk]*source[z-dz][kk];    
    for (int s=1;s<=4;++s){
      C1k[s] += tb/2.0 * s_derivative[z-dz][s] - temp_vector[s];
      C1k_pert[s] += tb/2.0*s_derivative_pert[z-dz][s] + tb_pert/2.0*s_derivative[z-dz][s] - temp_vector_pert[s];
    }

    del_ft2dim(temp_matrix,1,4,1,4);
    del_ft2dim(temp_matrix_pert,1,4,1,4);
    delete [](temp_vector+1);
    delete [](temp_vector_pert+1);

    // Done with C1k 

    fp_t * C0k; // Now equation (24)
    fp_t * C0k_pert;

    temp_matrix = ft2dim(1,4,1,4);
    temp_matrix_pert = ft2dim(1,4,1,4);
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp){
        temp_matrix[s][sp] = -tb/2.0 * K_reduced[z][s][sp];
        temp_matrix_pert[s][sp] = -tb_pert/2.0*K_reduced[z][s][sp] -tb/2.0 * K_reduced_pert[z][s][sp];
      }
      temp_matrix[s][s] += 1.0;
    }

    C0k = multiply_vector(temp_matrix, source[z],4);
    C0k_pert = new fp_t[4]-1;
    memset(C0k_pert+1,0,4*sizeof(fp_t));
    for(int s=1;s<=4;++s)
      for (int kk=1;kk<=4;++kk)
        C0k_pert[s] += temp_matrix[s][kk]*source_pert[z][kk] + temp_matrix_pert[s][kk]*source[z][kk];

    for (int s=1;s<=4;++s){
      C0k[s] -= tb/2.0 * s_derivative[z][s];
      C0k_pert[s] -= (tb_pert/2.0 * s_derivative[z][s] + tb/2.0 * s_derivative_pert[z][s]);
    }

    del_ft2dim(temp_matrix,1,4,1,4);
    del_ft2dim(temp_matrix_pert,1,4,1,4);

    fp_t ** c0k; // This is the other part of equation (24)
    c0k = multiply_square(K_reduced[z],K_reduced[z],4);
    fp_t ** c0k_pert;
    c0k_pert = ft2dim(1,4,1,4);
    memset(c0k_pert[1]+1,0,16*sizeof(fp_t));
    for (int s=1;s<=4;++s)
      for (int sp=1;sp<=4;++sp)
        for (int kk=1;kk<=4;++kk)
          c0k_pert[s][sp] += K_reduced[z][s][kk] * K_reduced_pert[z][kk][sp] + K_reduced_pert[z][s][kk] * K_reduced[z][kk][sp];


    for (int s=1;s<=4;++s)
      for (int sp=1;sp<=4;++sp){
        c0k[s][sp] += K_reduced[z][s][sp] + K_reduced_derivative[z][s][sp];
        c0k_pert[s][sp] += K_reduced_pert[z][s][sp] + K_reduced_derivative_pert[z][s][sp];
        c0k_pert[s][sp] = c0k_pert[s][sp] * tb/2.0 + c0k[s][sp] * tb_pert/2.0;
        c0k[s][sp] *= tb/2.0;
        c0k[s][sp] -= K_reduced[z][s][sp];
        c0k_pert[s][sp] -= K_reduced_pert[z][s][sp];
      }

   // A is easy:
    for (int s=1;s<=4;++s){
      for (int sp=1;sp<=4;++sp){
        A[s][sp] = w1 * K_reduced[z][s][sp] - w2/2.0 * c0k[s][sp];
        A_pert[s][sp] = w1_pert*K_reduced[z][s][sp] + w1*K_reduced_pert[z][s][sp] - w2_pert/2.0*c0k[s][sp] - w2/2.0*c0k_pert[s][sp];
      }
      A[s][s] += 1.0;
    }

    // Now B:

    temp_vector = multiply_vector(K_reduced[z-dz], I[z-dz],4);
    temp_vector_pert = new fp_t [4]-1;
    memset(temp_vector_pert+1,0,4*sizeof(fp_t));
    for(int s=1;s<=4;++s)
      for (int kk=1;kk<=4;++kk)
        temp_vector_pert[s] += K_reduced[z-dz][s][kk] * I_pert[z-dz][kk] + K_reduced_pert[z-dz][s][kk] * I[z-dz][kk];

    for (int s=1;s<=4;++s){
      B[s] = -w0 * temp_vector[s] + w0 * source[z-dz][s] + w1 * source[z][s]  + I[z-dz][s] * etb + w2/2.0 * (C1k[s] + C0k[s]);
      B_pert[s] = -w0_pert * temp_vector[s] -w0 * temp_vector_pert[s] + w0_pert * source[z-dz][s] + w0 * source_pert[z-dz][s] + 
        w1_pert * source[z][s]  + w1 * source_pert[z][s]  + I_pert[z-dz][s] * etb + I[z-dz][s] * etb_pert + w2_pert/2.0 * (C1k[s] + C0k[s]) + w2/2.0 * (C1k_pert[s] + C0k_pert[s]);
    }

    delete [](temp_vector+1);
    delete [](temp_vector_pert+1);
    delete [](C0k+1);
    delete [](C0k_pert+1);
    delete [](C1k+1);
    delete [](C1k_pert+1);
    del_ft2dim(c0k,1,4,1,4);
    del_ft2dim(c0k_pert,1,4,1,4);
    
    // After this we solve:
    fp_t * A_to_solve = A[1]+1;
    fp_t * b = B+1;
    Crout(4,A_to_solve, A_decomposed_by_LU);
    solveCrout(4,A_decomposed_by_LU,b,I_solution);
    memcpy(I[z]+1,I_solution,4*sizeof(fp_t));

    //printf(" %d %e\n", z, I[z][1]);
    //printf(" %d %e\n", z, I[z-dz][1]*etb + w0*source[z-dz][1] + w1*source[z][1] + w2/2.0*(source[z][1]+source[z-dz][1]+tb/2.0*(s_derivative[z-dz][1]-s_derivative[z][1])));

    //fp_t I_pert_by_hand = I_pert[z-dz][1]*etb + I[z-dz][1]*etb_pert + w0_pert*source[z-dz][1] + w0*source_pert[z-dz][1] + w1_pert*source[z][1] + w1*source_pert[z][1] + 
    //w2_pert/2.0*(source[z][1]+source[z-dz][1]+tb/2.0*(s_derivative[z-dz][1]-s_derivative[z][1])) + 
    //w2/2.0*(source_pert[z][1]+source_pert[z-dz][1]+tb/2.0*(s_derivative_pert[z-dz][1]-s_derivative_pert[z][1])) + 
    //w2/2.0*tb_pert/2.0*(s_derivative[z-dz][1]-s_derivative[z][1]);


    for (int s=1;s<=4;++s)
      for (int kk=1;kk<=4;++kk)
        B_pert[s] -= A_pert[s][kk] * I[z][kk];

    fp_t * b_pert = B_pert+1;

    solveCrout(4,A_decomposed_by_LU,b_pert,I_solution);
    memcpy(I_pert[z]+1,I_solution,4*sizeof(fp_t));

    //I_pert[z][1] = I_pert_by_hand; 
    // And compute the local operator
    if(L) l[z] = (w1+0.5*w2);

    del_ft2dim(A,1,4,1,4);
    delete [](B+1);
    del_ft2dim(A_pert,1,4,1,4);
    delete [](B_pert+1);
  }

  // Still the final point in the atmosphere is not covered here. That one we will compute in the linear approximation, both for 
  // the opacity and for the source function.
  if (zl != zh){
    int z = zh; // Just to be absolutely sure as I am a nub programmer sometimes (Milic)
    for (int s=1;s<=4;++s){
      I[z][s] = I[z-dz][s];
      I_pert[z][s] = I_pert[z-dz][s];
    }
  }
    
  delete [](delta_tau+x3l);
  delete [](opp_derivative+x3l);
  del_ft2dim(source,x3l,x3h,1,4);
  del_ft2dim(s_derivative,x3l,x3h,1,4);
  del_ft3dim(K_reduced,x3l,x3h,1,4,1,4);
  del_ft3dim(K_reduced_derivative,x3l,x3h,1,4,1,4);
  delete [](delta_tau_pert+x3l);
  delete [](opp_derivative_pert+x3l);
  del_ft2dim(source_pert,x3l,x3h,1,4);
  del_ft2dim(s_derivative_pert,x3l,x3h,1,4);
  del_ft3dim(K_reduced_pert,x3l,x3h,1,4,1,4);
  del_ft3dim(K_reduced_derivative_pert,x3l,x3h,1,4,1,4);
  
  delete []A_decomposed_by_LU;
  delete []I_solution;

  return 0;

}

// Numerical version for propagation of polarized intensity. Works fine for T and vt, for example. But for density not, because numbers are very slow,
// and hence, it simply fails, we will need to re-write it somehow.
int atmos_ppbez::formal_pert_numerical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t **** Intensity = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t *** alo_temp = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    
  // Then we perform a formal solution:
  formal(rt_grid, Intensity, alo_temp, op, em, theta, phi, boundary);
  memset(dS[x3l][x1l][x2l][x3l]+1,0,(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  
  // At the moment we do it in the numerical way:
  
  fp_t **** Intensity_perturbed = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t ***** opacity_perturbed = ft5dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  fp_t **** emissivity_perturbed = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  
  
  for (int x3k=x3l;x3k<=x3h;++x3k){

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s){
            emissivity_perturbed[x1i][x2i][x3i][s] = 0.5*em_pert[x3k][x1i][x2i][x3i][s]+em[x1i][x2i][x3i][s];
            for (int sp=1;sp<=4;++sp)
              opacity_perturbed[x1i][x2i][x3i][s][sp] = 0.5*op_pert[x3k][x1i][x2i][x3i][s][sp]+op[x1i][x2i][x3i][s][sp];
          }

    formal(rt_grid, Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    //printf("Up : %d %e \n", x3k, Intensity_perturbed[x1l][x2l][x3l][1]);

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s)
            dS[x3k][x1i][x2i][x3i][s] = Intensity_perturbed[x1i][x2i][x3i][s];

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s){
            emissivity_perturbed[x1i][x2i][x3i][s] = -0.5*em_pert[x3k][x1i][x2i][x3i][s]+em[x1i][x2i][x3i][s];
            for (int sp=1;sp<=4;++sp)
              opacity_perturbed[x1i][x2i][x3i][s][sp] = -0.5*op_pert[x3k][x1i][x2i][x3i][s][sp]+op[x1i][x2i][x3i][s][sp];
          }
          
    formal(rt_grid, Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    //printf("Low: %d %e \n", x3k, Intensity_perturbed[x1l][x2l][x3l][1]);


    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s)
            dS[x3k][x1i][x2i][x3i][s] -= Intensity_perturbed[x1i][x2i][x3i][s];
    

  }
  del_ft4dim(Intensity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  del_ft5dim(opacity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  del_ft4dim(emissivity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  
  del_ft4dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  del_ft3dim(alo_temp, x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
  
}

// Numerical version for propagation of polarized intensity. This is one for calculating responses to nodes directly.
// It is the same as the previous one, except it considers "parameters", that is "nodes" all together. 
int atmos_ppbez::formal_pert_numerical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, 
  fp_t theta, fp_t phi, fp_t boundary, int N_parameters){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t **** Intensity = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t *** alo_temp = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    
  // Then we perform a formal solution:
  formal(rt_grid, Intensity, alo_temp, op, em, theta, phi, boundary);
  memset(dS[x3l][x1l][x2l][x3l]+1,0,N_parameters*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  
  // At the moment we do it in the numerical way:
  
  fp_t **** Intensity_perturbed = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t ***** opacity_perturbed = ft5dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  fp_t **** emissivity_perturbed = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  
  
  for (int p=1;p<=N_parameters;++p){

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s){
            emissivity_perturbed[x1i][x2i][x3i][s] = 0.5*em_pert[p][x1i][x2i][x3i][s]+em[x1i][x2i][x3i][s];
            for (int sp=1;sp<=4;++sp)
              opacity_perturbed[x1i][x2i][x3i][s][sp] = 0.5*op_pert[p][x1i][x2i][x3i][s][sp]+op[x1i][x2i][x3i][s][sp];
          }

    formal(rt_grid, Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    //printf("Up : %d %e \n", x3k, Intensity_perturbed[x1l][x2l][x3l][1]);

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s)
            dS[p][x1i][x2i][x3i][s] = Intensity_perturbed[x1i][x2i][x3i][s];

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s){
            emissivity_perturbed[x1i][x2i][x3i][s] = -0.5*em_pert[p][x1i][x2i][x3i][s]+em[x1i][x2i][x3i][s];
            for (int sp=1;sp<=4;++sp)
              opacity_perturbed[x1i][x2i][x3i][s][sp] = -0.5*op_pert[p][x1i][x2i][x3i][s][sp]+op[x1i][x2i][x3i][s][sp];
          }
          
    formal(rt_grid, Intensity_perturbed, alo_temp, opacity_perturbed, emissivity_perturbed, theta, phi, boundary);

    //printf("Low: %d %e \n", x3k, Intensity_perturbed[x1l][x2l][x3l][1]);


    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s)
            dS[p][x1i][x2i][x3i][s] -= Intensity_perturbed[x1i][x2i][x3i][s];
    

  }
  del_ft4dim(Intensity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  del_ft5dim(opacity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  del_ft4dim(emissivity_perturbed, x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  
  del_ft4dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  del_ft3dim(alo_temp, x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
  
}

// Numerical version for propagation of 
int atmos_ppbez::formal_pert_analytical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary){

  // I sugggest we first do one formal solution with extraction of the full lambda operator.
  fp_t **** Intensity = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
    
  // Then we perform a formal solution:
  formal(rt_grid, Intensity, 0, op, em, theta, phi, boundary);
  memset(dS[x3l][x1l][x2l][x3l]+1,0,(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  
  // At the moment we do it in the numerical way:

  fp_t **** Intensity_perturbation = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
   
  // Then we just perform a formal solution at each depth:
  for (int x3k=x3l;x3k<=x3h;++x3k){
    formal_with_responses_full(rt_grid, Intensity,0,op,em,op_pert[x3k],em_pert[x3k],Intensity_perturbation,theta,phi,boundary);
  
    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          for (int s=1;s<=4;++s)
            dS[x3k][x1i][x2i][x3i][s] = Intensity_perturbation[x1i][x2i][x3i][s];
        }
  }
  del_ft4dim(Intensity_perturbation, x1l,x1h,x2l,x2h,x3l,x3h,1,4); 
  del_ft4dim(Intensity, x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  
  return 0;
  
}





