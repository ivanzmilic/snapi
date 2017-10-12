#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include "types.h"
#include "io.h"
#include "acfg.h"
#include "atmos.h"
#include "pack.h"
#include "const.h"

#include "atmos_pp.h"
#include "atmos_ppbez.h"
#include "mathtools.h"

#define TAU_MIN -5.0
#define TAU_MAX  1.0
#define TINY 1E-5

atmos_pp *atmos_pp_new(acfg *cfg,io_class &io_in)
{
//  io_in.msg(IOL_INFO,"atmos_pp_new: %s\n",cfg->rts);
  if(!strcmp(cfg->rts,"Bezier")) return new atmos_ppbez(cfg,io_in);
  if(!strcmp(cfg->rts,"QuadSC")) return new atmos_pp(cfg,io_in);
  return new atmos_pp(cfg,io_in);
}

atmos_pp *atmos_pp_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in)
{
  int08_t rtstype;
  offs+=::unpack(buf+offs,rtstype);
  switch(rtstype){
    case(ATMOS_RTS_BEZ): return new atmos_ppbez(buf,offs,do_swap,io_in);
    case(ATMOS_RTS_QSC): return new atmos_pp(buf,offs,do_swap,io_in);
  }
  return new atmos_pp(buf,offs,do_swap,io_in);
}

atmos_pp::atmos_pp(acfg *cfg,io_class &io_in):atmosphere(cfg,io_in)
{
  gtype=ATMOS_GEOM_PP;
  rtstype=ATMOS_RTS_QSC;
  if((x1l!=x1h)||(x2l!=x2h)) io.msg(IOL_ERROR,"atmos_pp::atmos_pp: plane parallel atmosphere may not extend over extended region (%d:%d x %d:%d) \n",x1l,x1h,x2l,x2h);
  io.msg(IOL_DEB1,"atmos_pp::atmos_pp: plane parallel atmosphere\n");
}

atmos_pp::atmos_pp(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):atmosphere(buf,offs,do_swap,io_in)
{
  gtype=ATMOS_GEOM_PP;
  rtstype=ATMOS_RTS_QSC;
  if((x1l!=x1h)||(x2l!=x2h)) io.msg(IOL_ERROR,"atmos_pp::atmos_pp: plane parallel atmosphere may not extend over extended region (%d:%d x %d:%d) \n",x1l,x1h,x2l,x2h);
  io.msg(IOL_DEB1,"atmos_pp::atmos_pp: plane parallel atmosphere\n");
}

atmos_pp::~atmos_pp(void)
{
}
// ================================================================================================

int atmos_pp::build_from_nodes(model * atmos_model){
  
  // We create the atmosphere on the basis of the previous atmosphere from given model. 
  // Model is the set of nodes for Temperature, microturbulent/systematic velocity and mag field. 
  // We enforce hydrostatic equilibrium by using the first order numerical scheme for differential eq. 
  // dp/dh = - \rho * g

  io.msg(IOL_INFO, "atmos_pp::build_from_nodes : building from nodes (%d..%d)\n",x3h,x3l);
  
  // Grid in tau, equidistantly spaced in log scale:
  int N_depths = x3h-x3l+1;
  fp_t * logtau = new fp_t [N_depths] - x3l; // This is tau500
  fp_t tau_min = TAU_MIN;
  fp_t tau_max = TAU_MAX;
  
  for (int x3i=x3l;x3i<=x3h;++x3i)
    logtau[x3i] = tau_min + (tau_max-tau_min) / (x3h-x3l) * (x3i-x3l);
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        tau_referent[x1i][x2i][x3i] = -pow(10.0,logtau[x3i]);


  // Interpolate from nodes all the quantities to this tau grid.

  // Temperature:
  int N_nodes_temp = atmos_model->get_N_nodes_temp();
  if (N_nodes_temp < 1){
    io.msg(IOL_ERROR, "atmos_pp::build_from_nodes : can't make atmosphere with less than one T node. Quitting...\n");
    return 1;  
  }
  fp_t * temp_nodes_tau = atmos_model->get_temp_nodes_tau();
  fp_t * temp_nodes_temp = atmos_model->get_temp_nodes_temp();
  atmospheric_interpolation(temp_nodes_tau, temp_nodes_temp, N_nodes_temp, logtau, T[x1l][x2l], x3l, x3h, 1);
  delete[](temp_nodes_tau+1);
  delete[](temp_nodes_temp+1);

  // Microturbulent velocity:
  int N_nodes_vt = atmos_model->get_N_nodes_vt();
  if (N_nodes_vt){
    fp_t * vt_nodes_tau = atmos_model->get_vt_nodes_tau();
    fp_t * vt_nodes_vt = atmos_model->get_vt_nodes_vt();
    atmospheric_interpolation(vt_nodes_tau, vt_nodes_vt, N_nodes_vt, logtau, Vt[x1l][x2l], x3l, x3h, 0);
    delete[](vt_nodes_tau+1);
    delete[](vt_nodes_vt+1);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i)
      Vt[x1l][x2l][x3i] = 0.0;
  
  // Systematic velocity: 
  int N_nodes_vs = atmos_model->get_N_nodes_vs();
  if (N_nodes_vs){
    fp_t * vs_nodes_tau = atmos_model->get_vs_nodes_tau();
    fp_t * vs_nodes_vs = atmos_model->get_vs_nodes_vs();
    atmospheric_interpolation(vs_nodes_tau, vs_nodes_vs, N_nodes_vs, logtau, Vz[x1l][x2l], x3l, x3h, 0);
    delete[](vs_nodes_tau+1);
    delete[](vs_nodes_vs+1);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i)
      Vz[x1l][x2l][x3i] = 0.0;

  // Magnetic field:

  int N_nodes_B = atmos_model->get_N_nodes_B();
  fp_t * B = new fp_t[x3h-x3l+1] - x3l;
    
  if (N_nodes_B){
    fp_t * B_nodes_tau = atmos_model->get_B_nodes_tau();
    fp_t * B_nodes_B = atmos_model->get_B_nodes_B();
    atmospheric_interpolation(B_nodes_tau,B_nodes_B, N_nodes_B, logtau, B, x3l, x3h, 0);
    delete[](B_nodes_tau+1);
    delete[](B_nodes_B+1);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i)
      B[x3i] = 0.0; 

  
  int N_nodes_theta = atmos_model->get_N_nodes_theta();
  fp_t * theta = new fp_t[x3h-x3l+1] - x3l;
    
  if (N_nodes_theta){
    fp_t * theta_nodes_tau = atmos_model->get_theta_nodes_tau();
    fp_t * theta_nodes_theta = atmos_model->get_theta_nodes_theta();
    atmospheric_interpolation(theta_nodes_tau,theta_nodes_theta, N_nodes_theta, logtau, theta, x3l, x3h, 0);
    delete[](theta_nodes_tau+1);
    delete[](theta_nodes_theta+1);
  }
  else
    for (int x3i=x3l;x3i<=x3h;++x3i)
      theta[x3i] = TINY; 


  int N_nodes_phi = atmos_model->get_N_nodes_phi();
  fp_t * phi = new fp_t[x3h-x3l+1] - x3l;

  if (N_nodes_phi){
    fp_t * phi_nodes_tau = atmos_model->get_phi_nodes_tau();
    fp_t * phi_nodes_phi = atmos_model->get_phi_nodes_phi();    
    atmospheric_interpolation(phi_nodes_tau,phi_nodes_phi, N_nodes_phi, logtau, phi, x3l, x3h, 0);
    delete[](phi_nodes_tau+1);
    delete[](phi_nodes_phi+1);
  }
  else
    for (int x3i=x3l;x3i<=x3h;++x3i)
      phi[x3i] = TINY; 
 
  for (int x3i=x3l;x3i<=x3h;++x3i){
    Bx[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * cos(phi[x3i]);
    By[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * sin(phi[x3i]);
    Bz[x1l][x2l][x3i] = B[x3i] * cos(theta[x3i]);
  }
  delete[](B+x3l);
  delete[](theta+x3l);
  delete[](phi+x3l);


// ------------------------------------------------------------------------------------------------
  
  // Density is not given, it is deduced from the temperature using hydrostatic equilibrium. 
  // Need populations to compute opacity correctly:
  popsetup();
  
  int MAX_ITER = 20;
  fp_t break_me = 1E-2;
  
  // Initial guess for particple density at the top. Assume pressure of 0.3 in CGS
  Nt[x1l][x2l][x3l] = 0.3/k/T[x1l][x2l][x3l];

  // Solve chemeq and ltepops to get the numbers and to compute density and opacity
  for (int iter=0;iter<MAX_ITER;++iter){
    chemeq(atml, natm, T[x1l][x2l][x3l], Nt[x1l][x2l][x3l], Ne[x1l][x2l][x3l], x1l, x2l, x3l);
    for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3l], Ne[x1l][x2l][x3l], x1l, x2l, x3l);
    op_referent[x1l][x2l][x3l] = opacity_continuum(T[x1l][x2l][x3l], Ne[x1l][x2l][x3l], lambda_referent, x1l,x2l,x3l);
    rho[x1l][x2l][x3l] = atml[0]->get_total_pop(x1l,x2l,x3l) * 1.4 * 1.67E-24; // in gram/cm^3, approximate
    // Correction to local particle density:
    fp_t dN = (pow(10.0, logtau[x3l]) * rho[x1l][x2l][x3l] * grav_acc / op_referent[x1l][x2l][x3l]) / k / T[x1l][x2l][x3l] - Nt[x1l][x2l][x3l];
    Nt[x1l][x2l][x3l] += dN;
    if (fabs(dN/Nt[x1l][x2l][x3l])<break_me)
      break;
  }
  // Go downward, enforcing differential form of HE equation:
  for (int x3i=x3l+1;x3i<=x3h;++x3i){ 

    // We have the eq: delta_p/delta_tau = g / kappa_mean
    // where delta_p is STEP (not correction) in pressure, delta_tau in optical depth and 
    // kappa_mean is sqrt(chi_i * chi_i+1 / rho_i / rho_i+1
    // we know everything except p_i+1 and chi/rho _i+1,so we iterate on these
    // delta_tau is known:
    fp_t delta_tau = (pow(10,logtau[x3i]) - pow(10,logtau[x3i-1])); 
    // inital guess for the step in p
    fp_t delta_p =  rho[x1l][x2l][x3i-1] * grav_acc / op_referent[x1l][x2l][x3i-1] * delta_tau; 
    // initial guess for the local Nt:
    Nt[x1l][x2l][x3i] = (Nt[x1l][x2l][x3i-1] * k * T[x1l][x2l][x3i-1] + delta_p) / k / T[x1l][x2l][x3i];

    // Now we iterate using HE equation and rho(N),chi(N)
    for (int iter=0;iter<MAX_ITER;++iter){
      
      chemeq(atml, natm, T[x1l][x2l][x3i], Nt[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
      for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
      // new opacity:
      op_referent[x1l][x2l][x3i] = opacity_continuum(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], lambda_referent, x1l,x2l,x3i); // New opacity
      // new mass density:
      rho[x1l][x2l][x3i] = atml[0]->get_total_pop(x1l,x2l,x3i) * 1.4 * 1.67E-24;
      fp_t kappa_mean = sqrt(op_referent[x1l][x2l][x3i] * op_referent[x1l][x2l][x3i-1] / rho[x1l][x2l][x3i] / rho[x1l][x2l][x3i-1]); // Geometrical mean of the value of kappa
      // follows new value of delta_p
      delta_p = grav_acc/kappa_mean * delta_tau; 
      fp_t dN = (Nt[x1l][x2l][x3i-1] * k * T[x1l][x2l][x3i-1] + delta_p) / k / T[x1l][x2l][x3i] -Nt[x1l][x2l][x3i];
      Nt[x1l][x2l][x3i] += dN;
      if (fabs(dN/Nt[x1l][x2l][x3l])<break_me)
        break;
    }
    // new height grid, rarely needed
    fp_t d_h = (pow(10,logtau[x3i]) - pow(10,logtau[x3i-1])) / sqrt(op_referent[x1l][x2l][x3i] * op_referent[x1l][x2l][x3i-1]);
  }

  for (int x3i=x3l;x3i<=x3h;++x3i){
    chemeq(atml, natm, T[x1l][x2l][x3i], Nt[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
    for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
    op_referent[x1l][x2l][x3i] = opacity_continuum(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], lambda_referent, x1l,x2l,x3i);
    rho[x1l][x2l][x3i] = atml[0]->get_total_pop(x1l,x2l,x3i) * 1.4 * 1.67E-24;
  }
  delete [](logtau+x3l);
  popclean();
  return 0;
}

// ================================================================================================

int atmos_pp::interpolate_from_nodes(model * atmos_model){
  
  // Similar to atmos_pp:build_from_nodes(model * atmos_model) 
  // except there is no HE enforcing part. We only interpolate quantities:

  io.msg(IOL_INFO, "atmos_pp::interpolating from nodes \n");
  
  // First you need to make a grid in tau. This is like this by default
  int N_depths = x3h-x3l+1;
  fp_t * logtau = new fp_t [N_depths] - x3l; // This is tau500
  
  fp_t tau_min = TAU_MIN;
  fp_t tau_max = TAU_MAX;
  
  for (int x3i=x3l;x3i<=x3h;++x3i)
    logtau[x3i] = tau_min + (tau_max-tau_min) / (x3h-x3l) * (x3i-x3l);

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        tau_referent[x1i][x2i][x3i] = -pow(10.0,logtau[x3i]);

  // Temperature:
  int N_nodes_temp = atmos_model->get_N_nodes_temp();
  if (N_nodes_temp < 1){
    io.msg(IOL_ERROR, "atmos_pp::build_from_nodes : can't make atmosphere with less than one T node. Quitting...\n");
    return 1;  
  }
  fp_t * temp_nodes_tau = atmos_model->get_temp_nodes_tau();
  fp_t * temp_nodes_temp = atmos_model->get_temp_nodes_temp();
  atmospheric_interpolation(temp_nodes_tau, temp_nodes_temp, N_nodes_temp, logtau, T[x1l][x2l], x3l, x3h, 1);
  delete[](temp_nodes_tau+1);
  delete[](temp_nodes_temp+1);

  // Microturbulent velocity:
  int N_nodes_vt = atmos_model->get_N_nodes_vt();
  if (N_nodes_vt){
    fp_t * vt_nodes_tau = atmos_model->get_vt_nodes_tau();
    fp_t * vt_nodes_vt = atmos_model->get_vt_nodes_vt();
    atmospheric_interpolation(vt_nodes_tau, vt_nodes_vt, N_nodes_vt, logtau, Vt[x1l][x2l], x3l, x3h, 0);
    delete[](vt_nodes_tau+1);
    delete[](vt_nodes_vt+1);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i)
      Vt[x1l][x2l][x3i] = 0.0;
  
  // Systematic velocity: 
  int N_nodes_vs = atmos_model->get_N_nodes_vs();
  if (N_nodes_vs){
    fp_t * vs_nodes_tau = atmos_model->get_vs_nodes_tau();
    fp_t * vs_nodes_vs = atmos_model->get_vs_nodes_vs();
    atmospheric_interpolation(vs_nodes_tau, vs_nodes_vs, N_nodes_vs, logtau, Vz[x1l][x2l], x3l, x3h, 0);
    delete[](vs_nodes_tau+1);
    delete[](vs_nodes_vs+1);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i)
      Vz[x1l][x2l][x3i] = 0.0;

  // Magnetic field:

  int N_nodes_B = atmos_model->get_N_nodes_B();
  fp_t * B = new fp_t[x3h-x3l+1] - x3l;
    
  if (N_nodes_B){
    fp_t * B_nodes_tau = atmos_model->get_B_nodes_tau();
    fp_t * B_nodes_B = atmos_model->get_B_nodes_B();
    atmospheric_interpolation(B_nodes_tau,B_nodes_B, N_nodes_B, logtau, B, x3l, x3h, 0);
    delete[](B_nodes_tau+1);
    delete[](B_nodes_B+1);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i)
      B[x3i] = 0.0; 

  
  int N_nodes_theta = atmos_model->get_N_nodes_theta();
  fp_t * theta = new fp_t[x3h-x3l+1] - x3l;
    
  if (N_nodes_theta){
    fp_t * theta_nodes_tau = atmos_model->get_theta_nodes_tau();
    fp_t * theta_nodes_theta = atmos_model->get_theta_nodes_theta();
    atmospheric_interpolation(theta_nodes_tau,theta_nodes_theta, N_nodes_theta, logtau, theta, x3l, x3h, 0);
    delete[](theta_nodes_tau+1);
    delete[](theta_nodes_theta+1);
  }
  else
    for (int x3i=x3l;x3i<=x3h;++x3i)
      theta[x3i] = TINY;


  int N_nodes_phi = atmos_model->get_N_nodes_phi();
  fp_t * phi = new fp_t[x3h-x3l+1] - x3l;

  if (N_nodes_phi){
    fp_t * phi_nodes_tau = atmos_model->get_phi_nodes_tau();
    fp_t * phi_nodes_phi = atmos_model->get_phi_nodes_phi();    
    atmospheric_interpolation(phi_nodes_tau,phi_nodes_phi, N_nodes_phi, logtau, phi, x3l, x3h, 0);
    delete[](phi_nodes_tau+1);
    delete[](phi_nodes_phi+1);
  }
  else
    for (int x3i=x3l;x3i<=x3h;++x3i)
      phi[x3i] = TINY;
 
  for (int x3i=x3l;x3i<=x3h;++x3i){
    Bx[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * cos(phi[x3i]);
    By[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * sin(phi[x3i]);
    Bz[x1l][x2l][x3i] = B[x3i] * cos(theta[x3i]);
  }
  delete[](B+x3l);
  delete[](theta+x3l);
  delete[](phi+x3l);

  delete[](logtau+x3l);
  return 0;
}

//=========================================================================================================================================

int atmos_pp::enforce_hequilibrium(){
  
  // We enforce hydrostatic equilibrium by using the first order numerical scheme for differential eq. 
  // dp/dh = - \rho * g
  // In the longer run, should be used to simplyfy build_from_nodes function
  // The code assumes that your taugrid is provided

  io.msg(IOL_INFO, "atmos_pp::enforcing hydrostatic equlibrium (%d..%d)\n",x3h,x3l);
  
// ------------------------------------------------------------------------------------------------
  
  // In essence what we do is compute density from given temperature run 
  // Need populations to compute opacity correctly:
  popsetup();
  
  int MAX_ITER = 20;
  fp_t break_me = 1E-2;
  
  // Initial guess for particple density at the top. Assume pressure of 0.3 in CGS
  Nt[x1l][x2l][x3l] = 0.3/k/T[x1l][x2l][x3l];

  // We need a logtau array to make things consistent with the above
  fp_t * logtau = new fp_t [x3h-x3l+1]-x3l;
  for (int x3i=x3l;x3i<=x3h;++x3i)
    logtau[x3i] = log10(-tau_referent[x1l][x2l][x3i]);
  
  // Solve chemeq and ltepops to get the numbers and to compute density and opacity
  for (int iter=0;iter<MAX_ITER;++iter){
    chemeq(atml, natm, T[x1l][x2l][x3l], Nt[x1l][x2l][x3l], Ne[x1l][x2l][x3l], x1l, x2l, x3l);
    for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3l], Ne[x1l][x2l][x3l], x1l, x2l, x3l);
    op_referent[x1l][x2l][x3l] = opacity_continuum(T[x1l][x2l][x3l], Ne[x1l][x2l][x3l], lambda_referent, x1l,x2l,x3l);
    rho[x1l][x2l][x3l] = atml[0]->get_total_pop(x1l,x2l,x3l) * 1.4 * 1.67E-24; // in gram/cm^3, approximate
    // Correction to local particle density:
    fp_t dN = (pow(10.0, logtau[x3l]) * rho[x1l][x2l][x3l] * grav_acc / op_referent[x1l][x2l][x3l]) / k / T[x1l][x2l][x3l] - Nt[x1l][x2l][x3l];
    Nt[x1l][x2l][x3l] += dN;
    if (fabs(dN/Nt[x1l][x2l][x3l])<break_me)
      break;
  }
  // Go downward, enforcing differential form of HE equation:
  for (int x3i=x3l+1;x3i<=x3h;++x3i){ 

    // We have the eq: delta_p/delta_tau = g / kappa_mean
    // where delta_p is STEP (not correction) in pressure, delta_tau in optical depth and 
    // kappa_mean is sqrt(chi_i * chi_i+1 / rho_i / rho_i+1
    // we know everything except p_i+1 and chi/rho _i+1,so we iterate on these
    // delta_tau is known:
    fp_t delta_tau = (pow(10,logtau[x3i]) - pow(10,logtau[x3i-1])); 
    // inital guess for the step in p
    fp_t delta_p =  rho[x1l][x2l][x3i-1] * grav_acc / op_referent[x1l][x2l][x3i-1] * delta_tau; 
    // initial guess for the local Nt:
    Nt[x1l][x2l][x3i] = (Nt[x1l][x2l][x3i-1] * k * T[x1l][x2l][x3i-1] + delta_p) / k / T[x1l][x2l][x3i];

    // Now we iterate using HE equation and rho(N),chi(N)
    for (int iter=0;iter<MAX_ITER;++iter){
      
      chemeq(atml, natm, T[x1l][x2l][x3i], Nt[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
      for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
      // new opacity:
      op_referent[x1l][x2l][x3i] = opacity_continuum(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], lambda_referent, x1l,x2l,x3i); // New opacity
      // new mass density:
      rho[x1l][x2l][x3i] = atml[0]->get_total_pop(x1l,x2l,x3i) * 1.4 * 1.67E-24;
      fp_t kappa_mean = sqrt(op_referent[x1l][x2l][x3i] * op_referent[x1l][x2l][x3i-1] / rho[x1l][x2l][x3i] / rho[x1l][x2l][x3i-1]); // Geometrical mean of the value of kappa
      // follows new value of delta_p
      delta_p = grav_acc/kappa_mean * delta_tau; 
      fp_t dN = (Nt[x1l][x2l][x3i-1] * k * T[x1l][x2l][x3i-1] + delta_p) / k / T[x1l][x2l][x3i] -Nt[x1l][x2l][x3i];
      Nt[x1l][x2l][x3i] += dN;
      if (fabs(dN/Nt[x1l][x2l][x3l])<break_me)
        break;
    }
    // new height grid, rarely needed
    fp_t d_h = (pow(10,logtau[x3i]) - pow(10,logtau[x3i-1])) / sqrt(op_referent[x1l][x2l][x3i] * op_referent[x1l][x2l][x3i-1]);
  }

  for (int x3i=x3l;x3i<=x3h;++x3i){
    chemeq(atml, natm, T[x1l][x2l][x3i], Nt[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
    for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
    op_referent[x1l][x2l][x3i] = opacity_continuum(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], lambda_referent, x1l,x2l,x3i);
    rho[x1l][x2l][x3i] = atml[0]->get_total_pop(x1l,x2l,x3i) * 1.4 * 1.67E-24;
  }
  delete [](logtau+x3l);
  popclean();
  return 0;
}

// ================================================================================================



