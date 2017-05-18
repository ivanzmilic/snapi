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
  io.msg(IOL_INFO,"atmos_pp::atmos_pp: plane parallel atmosphere\n");
}

atmos_pp::atmos_pp(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):atmosphere(buf,offs,do_swap,io_in)
{
  gtype=ATMOS_GEOM_PP;
  rtstype=ATMOS_RTS_QSC;
  if((x1l!=x1h)||(x2l!=x2h)) io.msg(IOL_ERROR,"atmos_pp::atmos_pp: plane parallel atmosphere may not extend over extended region (%d:%d x %d:%d) \n",x1l,x1h,x2l,x2h);
  io.msg(IOL_INFO,"atmos_pp::atmos_pp: plane parallel atmosphere\n");
}

atmos_pp::~atmos_pp(void)
{
}

int atmos_pp::build_from_nodes(model * atmos_model){
  
  // This is a working version which modifies the atmosphere with respect to given 'model'

  io.msg(IOL_INFO, "atmos_pp::building from nodes (%d..%d)\n",x3h,x3l);
  
  // First you need to make a grid in tau. This is like this by default
  int N_depths = x3h-x3l+1;
  fp_t * logtau = new fp_t [N_depths] - x3l; // This is tau500
  // Uniform in log_tau from the uppeermost to the lowermost node
  int N_nodes_temp = atmos_model->get_N_nodes_temp();
  fp_t * temp_nodes_tau = atmos_model->get_temp_nodes_tau();
  
  fp_t tau_min = temp_nodes_tau[1];
  fp_t tau_max = temp_nodes_tau[N_nodes_temp];
  
  for (int x3i=x3l;x3i<=x3h;++x3i)
    logtau[x3i] = tau_min + (tau_max-tau_min) / (x3h-x3l) * (x3i-x3l);

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        tau_referent[x1i][x2i][x3i] = -pow(10.0,logtau[x3i]);


  // Then we need to interpolate all the quantities to this tau grid:

  // Temperature:
  fp_t * temp_nodes_temp = atmos_model->get_temp_nodes_temp();
  atmospheric_interpolation(temp_nodes_tau, temp_nodes_temp, N_nodes_temp, logtau, T[x1l][x2l], x3l, x3h);
  delete[](temp_nodes_tau+1);
  delete[](temp_nodes_temp+1);
  

  // Microturbulent velocity:
  int N_nodes_vt = atmos_model->get_N_nodes_vt();
  if (N_nodes_vt){
    fp_t * vt_nodes_tau = atmos_model->get_vt_nodes_tau();
    fp_t * vt_nodes_vt = atmos_model->get_vt_nodes_vt();
    atmospheric_interpolation(vt_nodes_tau, vt_nodes_vt, N_nodes_vt, logtau, Vt[x1l][x2l], x3l, x3h);
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
    atmospheric_interpolation(vs_nodes_tau, vs_nodes_vs, N_nodes_vs, logtau, Vz[x1l][x2l], x3l, x3h);
    delete[](vs_nodes_tau+1);
    delete[](vs_nodes_vs+1);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i)
      Vz[x1l][x2l][x3i] = 0.0;

  // Magnetic field:

  int N_nodes_B = atmos_model->get_N_nodes_B();
  int N_nodes_theta = atmos_model->get_N_nodes_theta();
  int N_nodes_phi = atmos_model->get_N_nodes_phi();
  
  if (N_nodes_B && N_nodes_theta && N_nodes_phi){
    fp_t * B_nodes_tau = atmos_model->get_B_nodes_tau();
    fp_t * B_nodes_B = atmos_model->get_B_nodes_B();
    fp_t * theta_nodes_tau = atmos_model->get_theta_nodes_tau();
    fp_t * theta_nodes_theta = atmos_model->get_theta_nodes_theta();
    fp_t * phi_nodes_tau = atmos_model->get_phi_nodes_tau();
    fp_t * phi_nodes_phi = atmos_model->get_phi_nodes_phi();

    fp_t * B = new fp_t[x3h-x3l+1] - x3l;
    fp_t * theta = new fp_t[x3h-x3l+1] - x3l;
    fp_t * phi = new fp_t[x3h-x3l+1] - x3l;
  
    atmospheric_interpolation(B_nodes_tau,B_nodes_B, N_nodes_B, logtau, B, x3l, x3h);
    atmospheric_interpolation(theta_nodes_tau,theta_nodes_theta, N_nodes_theta, logtau, theta, x3l, x3h);
    atmospheric_interpolation(phi_nodes_tau,phi_nodes_phi, N_nodes_phi, logtau, phi, x3l, x3h);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      Bx[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * cos(phi[x3i]);
      By[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * sin(phi[x3i]);
      Bz[x1l][x2l][x3i] = B[x3i] * cos(theta[x3i]);
      //printf("%d %e %e %e \n", x3i, Bx[x1l][x2l][x3i], By[x1l][x2l][x3i], Bz[x1l][x2l][x3i]);
    }
    delete[](B_nodes_tau+1);
    delete[](B_nodes_B+1);
    delete[](theta_nodes_tau+1);
    delete[](theta_nodes_theta+1);
    delete[](phi_nodes_tau+1);
    delete[](phi_nodes_phi+1);
    delete[](B+x3l);
    delete[](theta+x3l);
    delete[](phi+x3l);
  }
  else 
    for (int x3i=x3l;x3i<=x3h;++x3i){ // WE NEED TO SORT THIS OUT
      Bx[x1l][x2l][x3i] = 0.001;
      By[x1l][x2l][x3i] = 0.001;
      Bz[x1l][x2l][x3i] = 0.001;
    }

  // We do not need the nodes any more.

  //FILE * out;
  //out = fopen("t_test.txt", "w");
  //for (int x3i=x3l;x3i<=x3h;++x3i)
  //  fprintf(out,"%d %e %e \n", x3i, logtau[x3i], T[x1l][x2l][x3i]);
  //fclose(out);

  
// -----------------------------------------------------------------------------------------------------------------------------------------
  
  // Now we need to iterate for the top pressure (density), initial guess can be the initial pressure on top:
  
  //io.msg(IOL_INFO, "atmos_pp::enforcing HD equilibrium. \n");
  popsetup();


  fp_t const grav_acc = 274.88E2; // cm/s^2

  fp_t lambda_reference = 500E-7; // Usually a referent lambda is 500 nm, but in general could be different
  int MAX_ITER = 20;
  fp_t break_me = 1E-2;
  
  Nt[x1l][x2l][x3l] = 0.3/k/T[x1l][x2l][x3l];

  // Solve chemeq and ltepops to get the numbers and to compute density and opacity
  

  for (int iter=0;iter<MAX_ITER;++iter){ // No more than 20 iterations
    
    chemeq(atml, natm, T[x1l][x2l][x3l], Nt[x1l][x2l][x3l], Ne[x1l][x2l][x3l], x1l, x2l, x3l);
    for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3l], Ne[x1l][x2l][x3l], x1l, x2l, x3l);
    // Now from this we need to get mass density and opacity
    op_referent[x1l][x2l][x3l] = opacity_continuum(T[x1l][x2l][x3l], Ne[x1l][x2l][x3l], lambda_reference, x1l,x2l,x3l);
    rho[x1l][x2l][x3l] = atml[0]->get_total_pop(x1l,x2l,x3l) * 1.4 * 1.67E-24; // in gram/cm^3, approximate.
    
    fp_t dN = (pow(10.0, logtau[x3l]) * rho[x1l][x2l][x3l] * grav_acc / op_referent[x1l][x2l][x3l]) / k / T[x1l][x2l][x3l] - Nt[x1l][x2l][x3l];
    Nt[x1l][x2l][x3l] += dN; 
    
  }

  //printf("Iterated top density to is : %e pressure %e \n", Nt[x1l][x2l][x3l], Nt[x1l][x2l][x3l] * k * T[x1l][x2l][x3l]);

  // Now we have converged the pressure at the top. We can now start integrating downward, the simplest way is to use linear approximation for the pressure and opacity derivative
  
  for (int x3i=x3l+1;x3i<=x3h;++x3i){

    fp_t delta_tau = (pow(10,logtau[x3i]) - pow(10,logtau[x3i-1]));
    fp_t delta_p =  rho[x1l][x2l][x3i-1] * grav_acc / op_referent[x1l][x2l][x3i-1] * delta_tau; // Initial value for delta_p
    
    for (int iter=0;iter<MAX_ITER;++iter){
      
      Nt[x1l][x2l][x3i] = (Nt[x1l][x2l][x3i-1] * k * T[x1l][x2l][x3i-1] + delta_p) / k / T[x1l][x2l][x3i];
      chemeq(atml, natm, T[x1l][x2l][x3i], Nt[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
      for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
      
      op_referent[x1l][x2l][x3i] = opacity_continuum(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], lambda_reference, x1l,x2l,x3i); // New opacity
      rho[x1l][x2l][x3i] = atml[0]->get_total_pop(x1l,x2l,x3i) * 1.4 * 1.67E-24; // New mass density
      fp_t kappa_mean = sqrt(op_referent[x1l][x2l][x3i] * op_referent[x1l][x2l][x3i-1] / rho[x1l][x2l][x3i] / rho[x1l][x2l][x3i-1]); // Geometrical mean of the value of kappa
      // New correction:
      fp_t P_local_correction = Nt[x1l][x2l][x3i-1] * k * T[x1l][x2l][x3i-1] + grav_acc / kappa_mean * delta_tau - Nt[x1l][x2l][x3i] * k * T[x1l][x2l][x3i];
    }

    fp_t d_h = (pow(10,logtau[x3i]) - pow(10,logtau[x3i-1])) / sqrt(op_referent[x1l][x2l][x3i] * op_referent[x1l][x2l][x3i-1]);
  }

  // Recompute, just in case referent opacity (DEBUG)
  for (int x3i=x3l;x3i<=x3h;++x3i)
    op_referent[x1l][x2l][x3i] = opacity_continuum(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], lambda_reference, x1l,x2l,x3i); // New opacity

  for (int x3i=x3l;x3i<=x3h;++x3i){
    chemeq(atml, natm, T[x1l][x2l][x3i], Nt[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
    for (int a=0;a<natm;++a) atml[a]->lte(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], x1l, x2l, x3i);
    op_referent[x1l][x2l][x3i] = opacity_continuum(T[x1l][x2l][x3i], Ne[x1l][x2l][x3i], lambda_reference, x1l,x2l,x3i);
    rho[x1l][x2l][x3i] = atml[0]->get_total_pop(x1l,x2l,x3i) * 1.4 * 1.67E-24;
  }

  delete [](logtau+x3l);

  popclean();
  return 0;
}

int atmos_pp::interpolate_from_nodes(model * atmos_model){
  
  // This is a working version which modifies the atmosphere with respect to given 'model'

  io.msg(IOL_INFO, "atmos_pp::interpolating from nodes \n");
  
  // First you need to make a grid in tau. This is like this by default
  int N_depths = x3h-x3l+1;
  fp_t * logtau = new fp_t [N_depths] - x3l; // This is tau500
  // Uniform in log_tau from -6 to 1.5 in tau500
  for (int x3i=x3l;x3i<=x3h;++x3i)
    logtau[x3i] = -5.0 + 5.5 / (x3h-x3l) * (x3i-x3l);

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        tau_referent[x1i][x2i][x3i] = -pow(10.0,logtau[x3i]);

  // Then we need to interpolate all the quantities to this tau grid:

  // Temperature:
  int N_nodes_temp = atmos_model->get_N_nodes_temp();
  fp_t * temp_nodes_tau = atmos_model->get_temp_nodes_tau();
  fp_t * temp_nodes_temp = atmos_model->get_temp_nodes_temp();
  atmospheric_interpolation(temp_nodes_tau, temp_nodes_temp, N_nodes_temp, logtau, T[x1l][x2l], x3l, x3h);

  // Microturbulent velocity:
  int N_nodes_vt = atmos_model->get_N_nodes_vt();
  fp_t * vt_nodes_tau = atmos_model->get_vt_nodes_tau();
  fp_t * vt_nodes_vt = atmos_model->get_vt_nodes_vt();
  atmospheric_interpolation(vt_nodes_tau, vt_nodes_vt, N_nodes_vt, logtau, Vt[x1l][x2l], x3l, x3h);
  
  // Systematic velocity: 
  int N_nodes_vs = atmos_model->get_N_nodes_vs();
  fp_t * vs_nodes_tau = atmos_model->get_vs_nodes_tau();
  fp_t * vs_nodes_vs = atmos_model->get_vs_nodes_vs();
  atmospheric_interpolation(vs_nodes_tau, vs_nodes_vs, N_nodes_vs, logtau, Vz[x1l][x2l], x3l, x3h);

  // Magnetic field:
  fp_t * B = new fp_t[x3h-x3l+1] - x3l;
  fp_t * theta = new fp_t[x3h-x3l+1] - x3l;
  fp_t * phi = new fp_t[x3h-x3l+1] - x3l;

  int N_nodes_B = atmos_model->get_N_nodes_B();
  fp_t * B_nodes_tau = atmos_model->get_B_nodes_tau();
  fp_t * B_nodes_B = atmos_model->get_B_nodes_B();
  int N_nodes_theta = atmos_model->get_N_nodes_theta();
  fp_t * theta_nodes_tau = atmos_model->get_theta_nodes_tau();
  fp_t * theta_nodes_theta = atmos_model->get_theta_nodes_theta();
  int N_nodes_phi = atmos_model->get_N_nodes_phi();
  fp_t * phi_nodes_tau = atmos_model->get_phi_nodes_tau();
  fp_t * phi_nodes_phi = atmos_model->get_phi_nodes_phi();
  
  atmospheric_interpolation(B_nodes_tau,B_nodes_B, N_nodes_B, logtau, B, x3l, x3h);
  atmospheric_interpolation(theta_nodes_tau,theta_nodes_theta, N_nodes_theta, logtau, theta, x3l, x3h);
  atmospheric_interpolation(phi_nodes_tau,phi_nodes_phi, N_nodes_phi, logtau, phi, x3l, x3h);

  if (N_nodes_B && N_nodes_theta && N_nodes_phi){
    for (int x3i=x3l;x3i<=x3h;++x3i){
      Bx[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * cos(phi[x3i]);
      By[x1l][x2l][x3i] = B[x3i] * sin(theta[x3i]) * sin(phi[x3i]);
      Bz[x1l][x2l][x3i] = B[x3i] * cos(theta[x3i]);
      //printf("%d %e %e %e \n", x3i, Bx[x1l][x2l][x3i], By[x1l][x2l][x3i], Bz[x1l][x2l][x3i]);
    }
  }

  // We do not need the nodes any more.
  delete[](temp_nodes_tau+1);
  delete[](temp_nodes_temp+1);
  delete[](vt_nodes_tau+1);
  delete[](vt_nodes_vt+1);
  delete[](vs_nodes_tau+1);
  delete[](vs_nodes_vs+1);
  delete[](B_nodes_tau+1);
  delete[](B_nodes_B+1);
  delete[](theta_nodes_tau+1);
  delete[](theta_nodes_theta+1);
  delete[](phi_nodes_tau+1);
  delete[](phi_nodes_phi+1);
  delete[](B+x3l);
  delete[](theta+x3l);
  delete[](phi+x3l);
  delete[](logtau+x3l);

  return 0;
}

