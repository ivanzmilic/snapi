#include <stdlib.h>
#include <string.h>
#include <ctime>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "const.h"
#include "math.h"

#include "atmos.h"
#include "profiles.h"
#include "obs.h"
#include "mathtools.h"

#define DELTA 1E-3


observable * atmosphere::stokes_lm_fit(observable * spectrum_to_fit, fp_t theta, fp_t phi, model * model_to_fit){

  // Perform LM fitting and return the best fitting spectrum
  
  // First extract the spectrum from the observation:
  fp_t ** S_to_fit = spectrum_to_fit->get_S_to_fit(1,1);
  int nlambda = spectrum_to_fit->get_n_lambda_to_fit();
  fp_t * lambda = spectrum_to_fit->get_lambda_to_fit();

  set_grid(1);
  
  // Set initial value of Levenberg-Marquardt parameter
  fp_t lm_parameter = 1E3;
  fp_t lm_multiplicator = 10.0;
  
  // Some fitting related parameters:
  fp_t metric = 0.0;
  int iter = 0;
  int MAX_ITER = 20;
  fp_t * chi_to_track = 0;
  int n_chi_to_track = 0;
  int corrected = 1;
  int to_break = 0;
  
  fp_t ws[4]; ws[0] = 1.0; ws[1] = ws[2] = 0.0; ws[3] = 0.0; // weights for Stokes parameters
  fp_t scattered_light = spectrum_to_fit->get_scattered_light();
  fp_t spectral_broadening = spectrum_to_fit->get_spectral_broadening();
  fp_t qs_level = spectrum_to_fit->get_synth_qs();
  int n_stokes_to_fit = 0; // A complicated piece of code to list what we want to fit
  for (int s=0;s<4;++s) 
    if (ws[s]) 
      ++n_stokes_to_fit;
  int stokes_to_fit[n_stokes_to_fit];
  int counter = 0;
  for (int s=0;s<4;++s) if (ws[s]){
    stokes_to_fit[counter] = s+1;
    ++counter;
  }
  
  fp_t *noise = new fp_t [nlambda]-1; // wavelength dependent noise
  for (int l=1;l<=nlambda;++l)
   noise[l] = sqrt(S_to_fit[1][l] * S_to_fit[1][1]) * 5E-3;
  
  observable * current_obs;
  fp_t *** derivatives_to_parameters;
  fp_t ** S_current;
  
  int N_parameters = model_to_fit->get_N_nodes_total();
  
  io.msg(IOL_INFO, "atmosphere::stokes_lm_fit : entering iterative procedure\n");
  
  for (iter=1;iter<=MAX_ITER;++iter){

    if (corrected){      
      // These quantities are only re-computed if the model has been modified:    
      derivatives_to_parameters = ft3dim(1,N_parameters,1,nlambda,1,4);     
      memset(derivatives_to_parameters[1][1]+1,0,N_parameters*nlambda*4*sizeof(fp_t));
      
      current_obs = obs_stokes_responses_to_nodes_new(model_to_fit, theta, phi, lambda, nlambda, derivatives_to_parameters); 
      current_obs->add_scattered_light(scattered_light,qs_level);
      if (spectral_broadening){
        current_obs->spectral_convolve(spectral_broadening,1,1);
        convolve_response_with_gauss(derivatives_to_parameters,lambda,N_parameters,nlambda,spectral_broadening);   
      }
      S_current = current_obs->get_S(1,1);
    }
    
    fp_t * residual = calc_residual(S_to_fit,S_current,nlambda,n_stokes_to_fit,stokes_to_fit);
    metric = calc_chisq(nlambda, n_stokes_to_fit, stokes_to_fit, residual, noise, ws);
    
    fprintf(stderr,"Iteration # %d metric = %e \n",iter,metric);
  
    fp_t ** J = ft2dim(1,n_stokes_to_fit*nlambda,1,N_parameters);
    for (int i=1;i<=N_parameters;++i) 
      for (int l=1;l<=nlambda;++l) 
        for (int s=1;s<=n_stokes_to_fit;++s){
          int stf = stokes_to_fit[s-1];
          J[(l-1)*n_stokes_to_fit+s][i] = derivatives_to_parameters[i][l][stf];
    }

    fp_t ** J_transpose = transpose(J,n_stokes_to_fit*nlambda,N_parameters);
    fp_t ** JTJ = multiply_with_transpose(J, n_stokes_to_fit*nlambda, N_parameters);
    
    for (int i=1;i<=N_parameters;++i) JTJ[i][i] *= (lm_parameter + 1.0);
    // Now correct
    fp_t * rhs = multiply_vector(J_transpose, residual, N_parameters, n_stokes_to_fit*nlambda);
    fp_t * correction = solve(JTJ, rhs, 1, N_parameters);
  
    // Apply the correction:
    model * test_model = clone(model_to_fit);
    test_model->correct(correction);
    build_from_nodes(test_model);
    // Compare again:
    observable *reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening); 
    fp_t ** S_reference = reference_obs->get_S(1,1);
    fp_t * residual_test = calc_residual(S_to_fit,S_reference,nlambda,n_stokes_to_fit,stokes_to_fit);
    fp_t metric_reference = calc_chisq(nlambda, n_stokes_to_fit, stokes_to_fit, residual_test, noise, ws);
    delete[](residual_test+1);
    
    fprintf(stderr,"Iteration # %d metric_reference = %e \n",iter,metric_reference);
    
    if (metric_reference < metric){
      // Everything is ok, and we can decrease lm_parameter:
      lm_parameter /= lm_multiplicator;
      look_for_best_lambda();
      model_to_fit->cpy_values_from(test_model);
      corrected=1;
      chi_to_track = add_to_1d_array(chi_to_track,n_chi_to_track,metric_reference);
      if (n_chi_to_track >=3)
        if ((chi_to_track[n_chi_to_track-2] - chi_to_track[n_chi_to_track-1]) / chi_to_track[n_chi_to_track-1] < DELTA &&
          (chi_to_track[n_chi_to_track-3] - chi_to_track[n_chi_to_track-2]) / chi_to_track[n_chi_to_track-2] < DELTA)
          to_break = 1;
    }
    else{
      lm_parameter *= lm_multiplicator;
      corrected = 0;
    }

    if(corrected || to_break || iter==MAX_ITER){
      
      del_ft3dim(derivatives_to_parameters,1,N_parameters,1,nlambda,1,4);
      delete current_obs;
      del_ft2dim(S_current,1,4,1,nlambda);
    }

    delete[](residual+1);
    delete test_model;
    del_ft2dim(J_transpose,1,N_parameters,1,n_stokes_to_fit*nlambda);
    del_ft2dim(JTJ,1,N_parameters,1,N_parameters);
    del_ft2dim(J,1,nlambda*n_stokes_to_fit,1,N_parameters);
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    delete [](rhs+1);
    delete [](correction+1);
    metric = 0.0;

    if (to_break)
      break;
  }

  io.msg(IOL_INFO, "fitting complete. Total number of iterations is : %d \n", iter-1);
  // Clean-up:
  del_ft2dim(S_to_fit,1,4,1,nlambda);
  delete[](lambda+1);
  delete[](noise+1);
  if (chi_to_track)
    delete[]chi_to_track;

  //model_to_fit->print();

  // Full version:
  lambda = spectrum_to_fit->get_lambda();
  nlambda = spectrum_to_fit->get_n_lambda();
  build_from_nodes(model_to_fit);
  observable *obs_to_return = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening);
      
  delete[](lambda+1);
  return obs_to_return;
}

fp_t * atmosphere::calc_residual(fp_t ** S_to_fit, fp_t ** S_current, int nlambda, int n_stokes_to_fit, int * stokes_to_fit){

  fp_t * residual = new fp_t[n_stokes_to_fit*nlambda]-1;
    for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=n_stokes_to_fit;++s){
        int stf=stokes_to_fit[s-1];
        residual[(l-1)*n_stokes_to_fit+stf] = S_to_fit[stf][l] - S_current[stf][l];
  }
  return residual;
}

fp_t atmosphere::calc_chisq(int nlambda, int n_stokes_to_fit, int * stokes_to_fit, fp_t * residual, fp_t * noise, fp_t * ws){
  fp_t chisq = 0.0;
  for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=n_stokes_to_fit;++s){
        int stf = stokes_to_fit[s-1];
        chisq += residual[(l-1)*n_stokes_to_fit+stf] * residual[(l-1)*n_stokes_to_fit+stf] 
          *ws[stf-1] / noise[l] / noise[l] / (n_stokes_to_fit*nlambda);
        residual[(l-1)*n_stokes_to_fit+stf] /= (noise[l]/noise[1]);
  }
  return chisq;
}

// ======================================================================================================================================================

observable * atmosphere::stokes_lm_nodeless_fit(observable * spectrum_to_fit, fp_t theta, fp_t phi){

  // Perform LM fitting and return the best fitting spectrum
  // Except here we are testing he possibility of doing a nodeless fit
  // Most of the function is based on the above atmosphere:stokes_lm_fit
  
  // First extract the spectrum from the observation:
  fp_t ** S_to_fit = spectrum_to_fit->get_S_to_fit(1,1);
  int nlambda = spectrum_to_fit->get_n_lambda_to_fit();
  fp_t * lambda = spectrum_to_fit->get_lambda_to_fit();

  set_grid(1); // we are still working in tau scale
  
  // Set initial value of Levenberg-Marquardt parameter
  fp_t lm_parameter = 1E3;
  fp_t lm_multiplicator = 5.0;
  
  // Some fitting related parameters:
  fp_t metric = 0.0;
  int iter = 0;
  int MAX_ITER = 50;
  int INITIAL_ITER = 15;
  fp_t * chi_to_track = 0;
  int n_chi_to_track = 0;
  int corrected = 1;
  int to_break = 0;
  
  fp_t ws[4]; ws[0] = 1.0; ws[1] = ws[2] = 0.0; ws[3] = 0.0; // weights for Stokes parameters
  fp_t scattered_light = spectrum_to_fit->get_scattered_light();
  fp_t spectral_broadening = spectrum_to_fit->get_spectral_broadening();
  fp_t qs_level = spectrum_to_fit->get_synth_qs();
  int n_stokes_to_fit = 0; // A complicated piece of code to list what we want to fit
  for (int s=0;s<4;++s) 
    if (ws[s]) 
      ++n_stokes_to_fit;
  int stokes_to_fit[n_stokes_to_fit];
  int counter = 0;
  for (int s=0;s<4;++s) if (ws[s]){
    stokes_to_fit[counter] = s+1;
    ++counter;
  }
  fp_t * logtau_grid = new fp_t [x3h+x3l+1]-x3l;
  for (int x3i=x3l;x3i<=x3h;++x3i)
    logtau_grid[x3i] = log10(-rt_grid[x3i]);
  
  fp_t *noise = new fp_t [nlambda]-1; // wavelength dependent noise
  for (int l=1;l<=nlambda;++l)
    noise[l] = sqrt(S_to_fit[1][1] * S_to_fit[1][l]) * 1E-2;

  
  fp_t **** full_stokes_responses;
  observable * current_obs;
  fp_t ** S_current;
  int N_parameters = 2 * (x3h-x3l+1);
  
  io.msg(IOL_INFO, "atmosphere::stokes_lm_nodeless_fit : entering iterative procedure\n");
  
  int NT=5; // Number of T nodes
  int NV=3; // Number of V nodes

  // Project the starting atmosphere on legendre basis:
  int ND=x3h-x3l+1;
  fp_t * alpha_T = project_on_legendre_basis(T[x1l][x2l]+x3l,logtau_grid+x3l,ND,NT);
  fp_t * T_r = reconstruct_from_legendre_basis(alpha_T,logtau_grid+x3l,ND,NT);
  fp_t * alpha_V = project_on_legendre_basis(Vz[x1l][x2l]+x3l,logtau_grid+x3l,ND,NV);
  fp_t * V_r = reconstruct_from_legendre_basis(alpha_V,logtau_grid+x3l,ND,NV);
  memcpy(T[x1l][x2l]+x3l,T_r,ND*sizeof(fp_t));
  memcpy(Vz[x1l][x3l]+x3l,V_r,ND*sizeof(fp_t));
  
  delete[]T_r;delete[]V_r;
  delete[]alpha_T;
  delete[]alpha_V;
  enforce_hequilibrium();


  for (iter=1;iter<=MAX_ITER;++iter){

    fprintf(stderr, "iteration #%d\n",iter);
    if (corrected){
      
      // These quantities are only re-computed if the model has been modified:    
      full_stokes_responses = ft4dim(1,7,x3l,x3h,1,nlambda,1,4);
      memset(full_stokes_responses[1][x3l][1]+1,0,7*(x3h-x3l+1)*nlambda*4*sizeof(fp_t));
      
      current_obs = obs_stokes_responses(theta, phi, lambda, nlambda,0,0,full_stokes_responses); 

      fp_t ** dN_dT = calculate_dN_dT();
        for (int x3k=x3l;x3k<=x3h;++x3k)
          for (int x3i=x3l;x3i<=x3h;++x3i)
            for (int l=1;l<=nlambda;++l)
              for (int s=1;s<=4;++s)
                full_stokes_responses[1][x3k][l][s] += full_stokes_responses[2][x3i][l][s] * dN_dT[x3i][x3k];
      del_ft2dim(dN_dT,x3l,x3h,x3l,x3h);

      current_obs->add_scattered_light(scattered_light,qs_level);
      if (spectral_broadening){
        current_obs->spectral_convolve(spectral_broadening,1,1);
        for (int p=1;p<=7;++p)
          convolve_response_with_gauss(full_stokes_responses[p]+x3l-1,lambda,(x3h-x3l+1),nlambda,spectral_broadening);   
      }
      S_current = current_obs->get_S(1,1);
    }
   
    fp_t * residual;
    residual = new fp_t [nlambda*n_stokes_to_fit]-1; 
    
    metric = 0.0;
  
    for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=n_stokes_to_fit;++s){
        int stf = stokes_to_fit[s-1];
        residual[(l-1)*n_stokes_to_fit+s] = S_to_fit[stf][l] - S_current[stf][l]; 
        metric += residual[(l-1)*n_stokes_to_fit+s] * residual[(l-1)*n_stokes_to_fit+s] 
          *ws[stf-1] / noise[l] / noise[l];
        residual[(l-1)*n_stokes_to_fit+s] /= (noise[l]/noise[1]);
    }

    printf("Starting metric = %f \n",metric);

    fp_t ** old_atmos = return_as_array();

    // What if I project corrections on legendre basis immediately
     
    observable * reference_obs;
    fp_t ** S_reference;

    fp_t BIC, BIC_TP, BIC_VP, BIC_TM, BIC_VM;
    fp_t amplify = 1.0;

    int NT_old = NT; int NV_old = NV;
    
    // For the proper one:
    fp_t ** corrections = calculate_legendre_corrections(full_stokes_responses,residual, lm_parameter, nlambda,NT,NV); 
    for (int x3i=x3l;x3i<=x3h;++x3i){
      T[x1l][x2l][x3i] += corrections[1][x3i-x3l+1];
      Vz[x1l][x2l][x3i] += corrections[2][x3i-x3l+1];
    }
    polish_extreme_values();    
    enforce_hequilibrium();
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening);
    S_reference = reference_obs->get_S(1,1);
    fp_t metric_reference = chi_sqr(S_to_fit,S_reference,noise,nlambda,stokes_to_fit,n_stokes_to_fit,ws);
    BIC = metric_reference + (NT+NV)*log(nlambda)*amplify;
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    
    // One more T mode:
    NT=NT_old+1;
    fp_t ** corrections_TP = calculate_legendre_corrections(full_stokes_responses,residual, lm_parameter, nlambda,NT,NV);
    copy_from_array(old_atmos);
    for (int x3i=x3l;x3i<=x3h;++x3i){
      T[x1l][x2l][x3i] += corrections_TP[1][x3i-x3l+1];
      Vz[x1l][x2l][x3i] += corrections_TP[2][x3i-x3l+1];
    }
    polish_extreme_values();    
    enforce_hequilibrium();
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening);
    S_reference = reference_obs->get_S(1,1);
    fp_t metric_reference_TP = chi_sqr(S_to_fit,S_reference,noise,nlambda,stokes_to_fit,n_stokes_to_fit,ws);
    BIC_TP = metric_reference_TP + (NT+NV)*log(nlambda)*amplify;
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    NT=NT_old;

    // One V mode more:
    NV=NV_old+1;
    fp_t ** corrections_VP = calculate_legendre_corrections(full_stokes_responses,residual, lm_parameter, nlambda,NT,NV);
    copy_from_array(old_atmos);
    for (int x3i=x3l;x3i<=x3h;++x3i){
      T[x1l][x2l][x3i] += corrections_VP[1][x3i-x3l+1];
      Vz[x1l][x2l][x3i] += corrections_VP[2][x3i-x3l+1];
    }
    polish_extreme_values();    
    enforce_hequilibrium();
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening);
    S_reference = reference_obs->get_S(1,1);
    fp_t metric_reference_VP = chi_sqr(S_to_fit,S_reference,noise,nlambda,stokes_to_fit,n_stokes_to_fit,ws);
    BIC_VP = metric_reference_VP + (NT+NV)*log(nlambda)*amplify;
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    NV=NV_old;

    // One less T mode:
    NT = (NT_old-1 > 0 ) ? NT_old-1 : NT_old;

    // Projecting part:
    copy_from_array(old_atmos);
    alpha_T = project_on_legendre_basis(T[x1l][x2l]+x3l,logtau_grid+x3l,ND,NT);
    T_r = reconstruct_from_legendre_basis(alpha_T,logtau_grid+x3l,ND,NT);
    memcpy(T[x1l][x2l]+x3l,T_r,ND*sizeof(fp_t));
    delete[]T_r;
    delete[]alpha_T;
    // -----------------------------------------------------------------------------

    fp_t ** corrections_TM = calculate_legendre_corrections(full_stokes_responses,residual, lm_parameter, nlambda,NT,NV); 
    for (int x3i=x3l;x3i<=x3h;++x3i){
      T[x1l][x2l][x3i] += corrections_TM[1][x3i-x3l+1];
      Vz[x1l][x2l][x3i] += corrections_TM[2][x3i-x3l+1];
    }
    polish_extreme_values();    
    enforce_hequilibrium();
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening);
    S_reference = reference_obs->get_S(1,1);
    fp_t metric_reference_TM = chi_sqr(S_to_fit,S_reference,noise,nlambda,stokes_to_fit,n_stokes_to_fit,ws);
    BIC_TM = metric_reference_TM + (NT+NV)*log(nlambda)*amplify;
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    NT = NT_old;

    // One less V mode:
    NV = (NV_old-1 > 0 ) ? NV_old-1 : NV_old;

    // Projecting part:
    copy_from_array(old_atmos);
    alpha_V = project_on_legendre_basis(Vz[x1l][x2l]+x3l,logtau_grid+x3l,ND,NV);
    V_r = reconstruct_from_legendre_basis(alpha_V,logtau_grid+x3l,ND,NV);
    memcpy(Vz[x1l][x2l]+x3l,V_r,ND*sizeof(fp_t));
    delete[]V_r;
    delete[]alpha_V;
    // -----------------------------------------------------------------------------
    
    fp_t ** corrections_VM = calculate_legendre_corrections(full_stokes_responses,residual, lm_parameter, nlambda,NT,NV);
    
    for (int x3i=x3l;x3i<=x3h;++x3i){
      T[x1l][x2l][x3i] += corrections_VM[1][x3i-x3l+1];
      Vz[x1l][x2l][x3i] += corrections_VM[2][x3i-x3l+1];
    }
    polish_extreme_values();    
    enforce_hequilibrium();
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening);
    S_reference = reference_obs->get_S(1,1);
    fp_t metric_reference_VM = chi_sqr(S_to_fit,S_reference,noise,nlambda,stokes_to_fit,n_stokes_to_fit,ws);
    BIC_VM = metric_reference_VM + (NT+NV)*log(nlambda)*amplify;
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    NV = NV_old;

    printf("Old metric = %f \n",metric);
    printf("%d %d %f %f %f %f %f \n", NT, NV, BIC,BIC_TP,BIC_VP,BIC_TM,BIC_VM);
    printf("%d %d %f %f %f %f %f \n", NT, NV, metric_reference,metric_reference_TP,
      metric_reference_VP,metric_reference_TM,metric_reference_VM);

    fp_t BICS[5] = {BIC,BIC_TP,BIC_VP,BIC_TM,BIC_VM};
    int indices[5] = {0,1,2,3,4};

    stupid_sort_indices_by_abs(BICS,indices,5);
    
    fp_t ** correctiones[] = {corrections,corrections_TP,corrections_VP,corrections_TM,corrections_VM};
    fp_t  metrics[] = {metric_reference,metric_reference_TP,metric_reference_VP,metric_reference_TM,metric_reference_VM};
  
    fp_t ** final_corrections;
    fp_t final_metric;

    int NT_test=NT, NV_test=NV;
    
    if (iter > INITIAL_ITER){

      int best = indices[4];
      if (best == 1)
        NT_test = NT_old+1;
      else if (best == 2)
        NV_test = NV_old+1;
      else if (best == 3)
        NT_test = (NT_old-1 > 0) ? NT_old-1 : NT_old;
      else if (best == 4)
        NV_test = (NV_old-1 > 0) ? NV_old-1 : NV_old;
      
      final_corrections = correctiones[indices[4]];
      final_metric = metrics[indices[4]];
    }
    else {
      final_corrections = corrections;
      final_metric = metric_reference;
    }
    // Apply corrections again
    copy_from_array(old_atmos);
    // Project again:
    printf("NT_test = %d NV_test = %d \n",NT_test,NV_test);
    if (NT_test < NT_old || NV_test < NV_old){ 
      alpha_T = project_on_legendre_basis(T[x1l][x2l]+x3l,logtau_grid+x3l,ND,NT_test);
      T_r = reconstruct_from_legendre_basis(alpha_T,logtau_grid+x3l,ND,NT_test);
      alpha_V = project_on_legendre_basis(Vz[x1l][x2l]+x3l,logtau_grid+x3l,ND,NV_test);
      V_r = reconstruct_from_legendre_basis(alpha_V,logtau_grid+x3l,ND,NV_test);
      memcpy(T[x1l][x2l]+x3l,T_r,ND*sizeof(fp_t));
      memcpy(Vz[x1l][x3l]+x3l,V_r,ND*sizeof(fp_t));
      delete[]T_r;delete[]V_r;
      delete[]alpha_T;
      delete[]alpha_V;
    }

    for (int x3i=x3l;x3i<=x3h;++x3i){
      T[x1l][x2l][x3i] += final_corrections[1][x3i-x3l+1];
      Vz[x1l][x2l][x3i] += final_corrections[2][x3i-x3l+1];
    }
    enforce_hequilibrium();

    del_ft2dim(corrections,1,2,x3l,x3h);
    del_ft2dim(corrections_TP,1,2,x3l,x3h);
    del_ft2dim(corrections_VP,1,2,x3l,x3h);
    del_ft2dim(corrections_TM,1,2,x3l,x3h);
    del_ft2dim(corrections_VM,1,2,x3l,x3h);
    
    // Compare again:
    
    if (final_metric < metric){
      // Everything is ok, and we can decrease lm_parameter:
      fprintf(stderr,"GOOD!\n");
      lm_parameter /= lm_multiplicator;
      corrected=1;
      chi_to_track = add_to_1d_array(chi_to_track,n_chi_to_track,metric_reference);
      del_ft2dim(old_atmos,1,12,x3l,x3h);
      if (n_chi_to_track >=3)
        if ((chi_to_track[n_chi_to_track-2] - chi_to_track[n_chi_to_track-1]) / chi_to_track[n_chi_to_track-1] < DELTA &&
          (chi_to_track[n_chi_to_track-3] - chi_to_track[n_chi_to_track-2]) / chi_to_track[n_chi_to_track-2] < DELTA)
          to_break = 1;
    }
    else{
      fprintf(stderr,"BAD!\n");
      copy_from_array(old_atmos);
      del_ft2dim(old_atmos,1,12,x3l,x3h);
      lm_parameter *= lm_multiplicator;
      corrected = 0;
    }

    if (iter > INITIAL_ITER && corrected){
      NT = NT_test;
      NV = NV_test;
    }

     // ===========================================================================
    //for (int x3i=x3l;x3i<=x3h;++x3i)
    //  fprintf(atmos_out,"%e %e %e \n",log10(-tau_referent[x1l][x2l][x3i]),T[x1l][x2l][x3i],Vz[x1l][x2l][x3i]);

    //for (int l=1;l<=nlambda;++l)
    //  fprintf(spectra_out, "%e %e\n",lambda[l],S_current[1][l]);
    // ===========================================================================


    if(corrected || to_break || iter==MAX_ITER){
      
      del_ft4dim(full_stokes_responses,1,7,x3l,x3h,1,nlambda,1,4);
      delete current_obs;
      del_ft2dim(S_current,1,4,1,nlambda);
    }
    
    delete[](residual+1);
    metric = 0.0;

    if (to_break)
      break;
  }

  io.msg(IOL_INFO, "fitting complete. Total number of iterations is : %d \n", iter-1);
  // Clean-up:
  del_ft2dim(S_to_fit,1,4,1,nlambda);
  delete[](lambda+1);
  delete[](noise+1);
  if (chi_to_track)
    delete[]chi_to_track;

  // Full version:
  lambda = spectrum_to_fit->get_lambda();
  nlambda = spectrum_to_fit->get_n_lambda();
  observable *obs_to_return = obs_stokes(theta, phi, lambda, nlambda);
  obs_to_return->add_scattered_light(scattered_light,qs_level);
  if (spectral_broadening)
    obs_to_return->spectral_convolve(spectral_broadening,1,1);

  //fclose(atmos_out);
  //fclose(spectra_out);
      
  delete[](lambda+1);
  delete[](logtau_grid+x3l);
  return obs_to_return;
}

// ========================================================================================================================================

observable* atmosphere::forward_evaluate(fp_t theta, fp_t phi, fp_t * lambda, int nlambda,
    fp_t scattered_light, fp_t qs, fp_t spectral_broadening){

  // proceduralized version which performs, in turn
  // 1) Stokes synthesis from current atmosphere
  // 2) Adds scattered light from the quiet Sun
  // 3) Convolves with instrumental profile

  observable *reference_obs = obs_stokes(theta, phi, lambda, nlambda);
  reference_obs->add_scattered_light(scattered_light,qs);
  if (spectral_broadening)
    reference_obs->spectral_convolve(spectral_broadening,1,1);  
  
  return reference_obs;
}

// ========================================================================================================================================

fp_t * atmosphere::calculate_svd_corrections(fp_t **** full_stokes_responses, fp_t * residual, fp_t lm_parameter, int nlambda, int param){
  
  // What we do here: 
  // We make a separate Hessian Matrix for each parameter we are fitting for.
  // Each Hessian Matrix we will decompose via SVD, we will keep some amount parameters and then 
  // we will use them to compute corrections. Hopefully this will work.
  // Smashing them all together seems wierd.

  fp_t * result;
  fp_t cutoff = 0.0;
  if (param == 1) cutoff = 1E-5;
  if (param == 4) cutoff = 1E-2;
  
  fprintf(stderr,"Calculating correction for parameter : %d \n",param);
  int ND = x3h-x3l+1;

  fp_t ** J = ft2dim(1,nlambda,1,ND);
  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int l=1;l<=nlambda;++l)
      J[l][x3i] = full_stokes_responses[param][x3i][l][1];

  fp_t ** JT = transpose(J,nlambda,ND);
  fp_t ** JTJ = multiply_with_transpose(J,nlambda,ND);
  for (int i=1;i<=ND;++i) JTJ[i][i] *= (lm_parameter + 1.0);
  fp_t * rhs = multiply_vector(JT, residual, ND, nlambda);
  // Now we have set everything up. Time to do SVD:
  fp_t ** U = JTJ;
  fp_t ** V = ft2dim(1,ND,1,ND);
  fp_t * w = new fp_t [ND]-1;

  if(svdcmp(U,ND,ND,w,V)<0) io.msg(IOL_WARN,"SVD: singular matrix?\n");
  
  // Now make pseudoinverse:
  fp_t ** U_T = transpose(U,ND,ND);
  fp_t ** V_T = transpose(V,ND,ND);
  fp_t ** W_inv = ft2dim(1,ND,1,ND);
  memset(W_inv[1]+1,0,ND*ND*sizeof(fp_t));
  for (int i=1;i<=ND;++i)
    if (w[i]/w[1] > cutoff)
      W_inv[i][i] = 1.0/w[i];

  fp_t ** temp = multiply_square(W_inv,U_T,ND);
  fp_t ** U_pseudoinverse = multiply_square(V,temp,ND);

  fp_t * correction = multiply_vector(U_pseudoinverse,rhs,ND);

  result = new fp_t[ND]-1;
  memcpy(result+1,correction+1,ND*sizeof(fp_t));

  del_ft2dim(J,1,nlambda,1,ND);
  del_ft2dim(JT,1,ND,1,nlambda);
  del_ft2dim(JTJ,1,ND,1,ND);
  del_ft2dim(V,1,ND,1,ND);
  del_ft2dim(V_T,1,ND,1,ND);
  del_ft2dim(U_T,1,ND,1,ND);
  del_ft2dim(W_inv,1,ND,1,ND);
  del_ft2dim(temp,1,ND,1,ND);
  del_ft2dim(U_pseudoinverse,1,ND,1,ND);
  delete[](w+1);
  delete[](correction+1);
  return result;
}

// ========================================================================================================================================

fp_t ** atmosphere::calculate_svd_corrections(fp_t **** full_stokes_responses, fp_t * residual, fp_t lm_parameter, int nlambda){
  
  // Overloaded version of the function from above. 
  // This one uses the approach that SIR and Andres are using. 
  // We decompose the matrix somehow.

  // Scale responses:
  fp_t T_scale = 5E3;
  fp_t V_scale = 1E5;
  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int l=1;l<=nlambda;++l){
      full_stokes_responses[1][x3i][l][1] *= T_scale;
      full_stokes_responses[4][x3i][l][1] *= V_scale;
  }

  int ND = x3h-x3l+1;
  fp_t ** J = ft2dim(1,nlambda,1,2*ND);
  for (int l=1;l<=nlambda;++l)
    for (int x3i=x3l;x3i<=x3h;++x3i){
      J[l][x3i-x3l+1] = full_stokes_responses[1][x3i][l][1];
      J[l][x3i-x3l+1+ND] = full_stokes_responses[4][x3i][l][1];
  }

  fp_t ** JT = transpose(J,nlambda,2*ND);
  fp_t ** JTJ = multiply_with_transpose(J,nlambda,2*ND);

  for (int i=1;i<=2*ND;++i) JTJ[i][i] *= (lm_parameter + 1.0);
  fp_t * rhs = multiply_vector(JT, residual, 2*ND, nlambda);

  // Separate the Hessian (JTJ) in two parts (in general case it 
  // will be in one part for each fitted parameter).
  fp_t ** Hessian_tmp = ft2dim(1,2*ND,1,2*ND);
  fp_t ** Hessian_vel = ft2dim(1,2*ND,1,2*ND);
  memset(Hessian_tmp[1]+1,0,4*ND*ND*sizeof(fp_t));
  memset(Hessian_vel[1]+1,0,4*ND*ND*sizeof(fp_t));
  for (int i=1;i<=ND;++i)
    for (int j=1;j<=2*ND;++j){
      Hessian_tmp[i][j] = JTJ[i][j];
      Hessian_vel[i+ND][j] = JTJ[i+ND][j];
  }

  fp_t * w;
  w = svd_treshold_square_matrix(Hessian_tmp,2*ND,1E-3,1);
  //for (int i=1;i<=2*ND;++i)
  //  fprintf(stderr,"%d %e \n",i,w[i]);
  delete[](w+1);
  w = svd_treshold_square_matrix(Hessian_vel,2*ND,1E-2,1);
  //for (int i=1;i<=2*ND;++i)
  //  fprintf(stderr,"%d %e \n",i,w[i]);
  delete[](w+1);
  
  for (int i=1;i<=2*ND;++i)
    for (int j=1;j<=2*ND;++j)
      JTJ[i][j] = Hessian_tmp[i][j] + Hessian_vel[i][j];

  w = svd_invert_square_matrix(JTJ,2*ND,1E-12,1);
  //for (int i=1;i<=2*ND;++i)
  //  fprintf(stderr,"%d %e \n",i,w[i]);
  delete[](w+1);

  fp_t * x = multiply_vector(JTJ,rhs,2*ND);
  fp_t ** corrections = ft2dim(1,2,1,ND);
  for (int i=1;i<=ND;++i){
    x[i] *= T_scale;
    x[i+ND] *= V_scale;
    corrections[1][i] = x[i];
    corrections[2][i] = x[i+ND];
  }
  delete[](x+1);

  // De-scale the responses:
  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int l=1;l<=nlambda;++l){
      full_stokes_responses[1][x3i][l][1] /= T_scale;
      full_stokes_responses[4][x3i][l][1] /= V_scale;
  }
  return corrections;
}

fp_t ** atmosphere::calculate_legendre_corrections(fp_t **** full_stokes_responses, fp_t * residual, fp_t lm_parameter, int nlambda,
  int NT, int NV){

  // Automatically asume that the logtau grid is equdistantly spaced and that the 
  // legendre polynomials are defined on the whole range
  fp_t * x = new fp_t[x3h-x3l+1]-x3l;
  for (int x3i=x3l;x3i<=x3h;++x3i){
    x[x3i] = -1.0 + (x3i-x3l)*2.0/(x3h-x3l);
  }

  fp_t scale_T = 5E3;
  fp_t scale_v = 1E5;

  // Scale response function:
  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int l=1;l<=nlambda;++l){
      full_stokes_responses[1][x3i][l][1] *= scale_T; // T
      full_stokes_responses[4][x3i][l][1] *= scale_v; // v
    }

  int N_components[2] = {NT,NV}; // # of components for T and V

  int N_parameters = N_components[0] + N_components[1]; // total number of parameters

  fp_t ** J = ft2dim(1,nlambda,1,N_parameters); // Jacobian
  memset(J[1]+1,0,nlambda*N_parameters*sizeof(fp_t));

  for (int i=1;i<=N_components[0];++i) // Calculate derivatives w.r.t. each Legendre polynomial
    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t P = Pn(i-1,x[x3i]); // legendre pol
      for (int l=1;l<=nlambda;++l)
        J[l][i] += full_stokes_responses[1][x3i][l][1] * P;
  }
  for (int i=1;i<=N_components[1];i++) // Calculate derivatives w.r.t. each Legendre polynomial
    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t P = Pn(i,x[x3i]); // legendre pol
      for (int l=1;l<=nlambda;++l)
        J[l][i+N_components[0]] += full_stokes_responses[4][x3i][l][1] * P;
  }
  
  fp_t ** JT = transpose(J,nlambda,N_parameters);
  fp_t ** JTJ = multiply_with_transpose(J,nlambda,N_parameters);

  for (int i=1;i<=N_parameters;++i) JTJ[i][i] *= (lm_parameter + 1.0);
  fp_t * rhs = multiply_vector(JT, residual, N_parameters, nlambda);

  del_ft2dim(JT,1,N_parameters,1,nlambda);
  del_ft2dim(J,1,nlambda,1,N_parameters);

  svd_invert_square_matrix(JTJ,N_parameters,0,0);

  fp_t * corrections_legendre_basis = multiply_vector(JTJ,rhs,N_parameters);

  del_ft2dim(JTJ,1,N_parameters,1,N_parameters);
  delete[](rhs+1);

  fp_t ** corrections = ft2dim(1,2,x3l,x3h);
  memset(corrections[1]+x3l,0,2*(x3h-x3l+1)*sizeof(fp_t));

  for (int i=1;i<=N_components[0];++i)
    for (int x3i=x3l;x3i<=x3h;++x3i)
      corrections[1][x3i] += Pn(i-1,x[x3i]) * corrections_legendre_basis[i]*scale_T;
  for (int i=1;i<=N_components[1];++i)
    for (int x3i=x3l;x3i<=x3h;++x3i)
      corrections[2][x3i] += Pn(i-1,x[x3i]) * corrections_legendre_basis[i+N_components[0]]*scale_v;

  delete[](corrections_legendre_basis+1);
  delete[](x+x3l);

  // Unscale response function:
  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int l=1;l<=nlambda;++l){
      full_stokes_responses[1][x3i][l][1] /= scale_T; // T
      full_stokes_responses[4][x3i][l][1] /= scale_v; // v
    }

  return corrections;
}


// ========================================================================================================================================

// FUNCTIONS THAT COMPUTE DERIVATIVES TO THE NODES:

 observable *atmosphere::obs_scalar_num_responses_to_nodes(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t ** response_to_parameters){

  // We first need to create the atmosphere from the model and compute the observable:

  build_from_nodes(atmos_model);
  observable *o = obs_scalar(theta, phi, lambda, nlambda);
  observable *o_temp;
  fp_t ** S_temp;

  io.msg(IOL_INFO, "atmosphere::computed observable from a model. \n");

  int N_parameters = atmos_model->get_N_nodes_total();
  
  io.msg(IOL_INFO, "atmosphere::numerically computing responses to total of %d model parameters. \n", N_parameters); 

  fp_t step = 5.0;

  fp_t * x3_old = new fp_t [x3h-x3l+1] - x3l;
  memcpy(x3_old+x3l,x3+x3l,(x3h-x3l+1)*sizeof(fp_t));

  for (int i=1;i<=N_parameters;++i){
    // First perturb to the one side
    atmos_model->perturb_node_value(i, step);

    build_from_nodes(atmos_model);

    o_temp = obs_scalar(theta, phi, lambda, nlambda);
    S_temp = o_temp->get_S(1,1);
    for (int l=1;l<=nlambda;++l)
      response_to_parameters[i][l] = S_temp[1][l];
    del_ft2dim(S_temp,1,1,1,nlambda);
    delete o_temp;
    
    atmos_model->perturb_node_value(i, -2.0 * step);
    build_from_nodes(atmos_model);
    o_temp = obs_scalar(theta, phi, lambda, nlambda);
    S_temp = o_temp->get_S(1,1);
    for (int l=1;l<=nlambda;++l){
      response_to_parameters[i][l] -= S_temp[1][l];
      response_to_parameters[i][l] /= (2.0 * step);
    }
    del_ft2dim(S_temp,1,1,1,nlambda);
    delete o_temp;

    // Finally return 
    atmos_model->perturb_node_value(i, 1.0 * step);
    build_from_nodes(atmos_model);
    
  }

  return o;
 }

  observable *atmosphere::obs_stokes_num_responses_to_nodes(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** response_to_parameters){

  // We first need to create the atmosphere from the model and compute the observable:

  build_from_nodes(atmos_model);
  observable *o = obs_stokes(theta, phi, lambda, nlambda);
  observable *o_temp;
  fp_t ** S_temp;

  io.msg(IOL_INFO, "atmosphere::computed polarized observable from a model. \n");

  int N_parameters = atmos_model->get_N_nodes_total();
  int N_nodes_T = atmos_model->get_N_nodes_temp();
  int N_nodes_vt = atmos_model->get_N_nodes_vt();
  int N_nodes_vs = atmos_model->get_N_nodes_vs();
  int N_nodes_B = atmos_model->get_N_nodes_B();
  
  io.msg(IOL_INFO, "atmosphere::numerically computing responses to total of %d model parameters. \n", N_parameters); 

  for (int i=1;i<=N_parameters;++i){

    fp_t step = 0.0;
    // First see what kind of step we need to use:
    if (i<=N_nodes_T)
      step = delta_T;
    else if (i<=N_nodes_T + N_nodes_vt)
      step = delta_vt;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs)
      step = delta_vr;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs+N_nodes_B)
      step = delta_B;    
    else 
      step = delta_angle;

    // First perturb to the one side
    atmos_model->perturb_node_value(i, step*0.5);

    build_from_nodes(atmos_model);

    o_temp = obs_stokes(theta, phi, lambda, nlambda);
    S_temp = o_temp->get_S(1,1);
    for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=4;++s)
      response_to_parameters[i][l][s] = S_temp[s][l];
    del_ft2dim(S_temp,1,4,1,nlambda);
    delete o_temp;
    
    atmos_model->perturb_node_value(i, -step);
    build_from_nodes(atmos_model);
    o_temp = obs_stokes(theta, phi, lambda, nlambda);
    S_temp = o_temp->get_S(1,1);
    for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=4;++s){
        response_to_parameters[i][l][s] -= S_temp[s][l];
        response_to_parameters[i][l][s] /= step;
    }
    del_ft2dim(S_temp,1,4,1,nlambda);
    delete o_temp;

    // Finally return 
    atmos_model->perturb_node_value(i, step*0.5);
    build_from_nodes(atmos_model);
    
  }
  return o;
 }

 observable *atmosphere::obs_scalar_num_responses_to_nodes_tau(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t ** response_to_parameters){

  // We first need to create the atmosphere from the model and compute the observable:

  build_from_nodes(atmos_model);
  observable *o = obs_scalar_tau(theta, phi, lambda, nlambda);
  observable *o_temp;
  fp_t ** S_temp;

  io.msg(IOL_INFO, "atmosphere::computed observable from a model. \n");

  int N_parameters = atmos_model->get_N_nodes_total();
  
  io.msg(IOL_INFO, "atmosphere::numerically computing responses to total of %d model parameters. \n", N_parameters); 

  fp_t step = 5.0;

  fp_t * x3_old = new fp_t [x3h-x3l+1] - x3l;
  memcpy(x3_old+x3l,x3+x3l,(x3h-x3l+1)*sizeof(fp_t));

  for (int i=1;i<=N_parameters;++i){
    // First perturb to the one side
    atmos_model->perturb_node_value(i, step);

    build_from_nodes(atmos_model);

    o_temp = obs_scalar_tau(theta, phi, lambda, nlambda);
    S_temp = o_temp->get_S(1,1);
    for (int l=1;l<=nlambda;++l)
      response_to_parameters[i][l] = S_temp[1][l];
    del_ft2dim(S_temp,1,1,1,nlambda);
    delete o_temp;
    
    atmos_model->perturb_node_value(i, -2.0 * step);
    build_from_nodes(atmos_model);
    o_temp = obs_scalar_tau(theta, phi, lambda, nlambda);
    S_temp = o_temp->get_S(1,1);
    for (int l=1;l<=nlambda;++l){
      response_to_parameters[i][l] -= S_temp[1][l];
      response_to_parameters[i][l] /= (2.0 * step);
    }
    del_ft2dim(S_temp,1,1,1,nlambda);
    delete o_temp;

    // Finally return 
    atmos_model->perturb_node_value(i, 1.0 * step);
    build_from_nodes(atmos_model);
    
  }

  return o;
 }


// Same as before but now using analytical approach which is supposed to be much much faster.
 observable *atmosphere::obs_scalar_responses_to_nodes(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t ** response_to_parameters){

  // We first need to create the atmosphere from the model and compute the observable:

  build_from_nodes(atmos_model);
  
  fp_t *** responses_depth_dependent;
  responses_depth_dependent = ft3dim(1,7,x3l,x3h,1,nlambda); 
  observable *o = obs_scalar_responses(theta, phi, lambda, nlambda, responses_depth_dependent);
  
  observable *o_temp;
  fp_t ** S_temp;

  io.msg(IOL_INFO, "atmosphere::computed observable from a model. \n");

  int N_parameters = atmos_model->get_N_nodes_total();
  memset(response_to_parameters[1]+1, 0, N_parameters*nlambda*sizeof(fp_t));

  io.msg(IOL_INFO, "atmosphere::analytically computing responses to total of %d model parameters. \n", N_parameters);

  // Lets plot some stuff here:
  FILE * output;
  output = fopen("perturbed_atmospheres.txt", "w");

  for (int x3i=x3l;x3i<=x3h;++x3i)
    fprintf(output,"%e %e %e \n", x3[x3i], T[x1l][x2l][x3i], Nt[x1l][x2l][x3i]); 

  fp_t step = delta_T; // We still need step because we want to numerically study the derivative of the interpolation scheme. It does not matter a lot.
  
  fp_t * depth = new fp_t [x3h-x3l+1];
  fp_t * original_depth = new fp_t [x3h-x3l+1];
  for (int x3i=0;x3i<x3h-x3l+1;++x3i)
    original_depth[x3i] = - x3[x3i+x3l];
     
  for (int i=1;i<=N_parameters;++i){
    
    // First perturb to the one side
  	fp_t * local_T_perturbations = new fp_t [x3h-x3l+1] - x3l;
  	fp_t * local_Nt_perturbations = new fp_t [x3h-x3l+1] - x3l;

    atmos_model->perturb_node_value(i, step * 0.5);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i)
      fprintf(output,"%e %e %e \n", x3[x3i], T[x1l][x2l][x3i], Nt[x1l][x2l][x3i]); 


    for (int x3i=0;x3i<x3h-x3l+1;++x3i)
      depth[x3i] = - x3[x3i+x3l];
  
    for (int x3i=x3l;x3i<=x3h;++x3i){
    	local_T_perturbations[x3i] = T[x1l][x2l][x3i];
      local_Nt_perturbations[x3i] = Nt[x1l][x2l][x3i];
      //local_T_perturbations[x3i] = interpol_1d(T[x1l][x2l]+x3l,depth,x3h-x3l+1,original_depth[x3i-x3l]);
      //local_Nt_perturbations[x3i] = interpol_1d(Nt[x1l][x2l]+x3l,depth,x3h-x3l+1,original_depth[x3i-x3l]);
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i)
      fprintf(output,"%e %e %e \n", x3[x3i], T[x1l][x2l][x3i], Nt[x1l][x2l][x3i]); 

    for (int x3i=0;x3i<x3h-x3l+1;++x3i)
      depth[x3i] = - x3[x3i+x3l];
    
    for (int x3i=x3l;x3i<=x3h;++x3i){
      local_T_perturbations[x3i] -= T[x1l][x2l][x3i];
      local_Nt_perturbations[x3i] -= Nt[x1l][x2l][x3i];
      //local_T_perturbations[x3i] -= interpol_1d(T[x1l][x2l]+x3l,depth,x3h-x3l+1,original_depth[x3i-x3l]);
      //local_Nt_perturbations[x3i] -= interpol_1d(Nt[x1l][x2l]+x3l,depth,x3h-x3l+1,original_depth[x3i-x3l]);
    	local_T_perturbations[x3i] /= (1.0 * step);
    	local_Nt_perturbations[x3i] /= (1.0 * step); 
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    build_from_nodes(atmos_model);
    
  
// ----------------------------------------------------------------------------------------------------------------------------------------------------

    // Now you need everything together:
    for (int x3i=x3l;x3i<=x3h;++x3i){
      //printf("%d %e %e \n", x3i, local_T_perturbations[x3i], local_Nt_perturbations[x3i]/Nt[x1l][x2l][x3i]);
    	for (int l=1;l<=nlambda;++l)
    	response_to_parameters[i][l] += responses_depth_dependent[1][x3i][l] * local_T_perturbations[x3i]
       + responses_depth_dependent[2][x3i][l] * local_Nt_perturbations[x3i];
    }

    delete [](local_T_perturbations+x3l);
    delete [](local_Nt_perturbations+x3l);
    
  }

  fclose(output);

  del_ft3dim(responses_depth_dependent,1,7,x3l,x3h,1,nlambda);
  delete []depth;
  delete []original_depth;
  return o;

 }

 // Same as before but now using analytical approach which is supposed to be much much faster.
 observable *atmosphere::obs_scalar_responses_to_nodes_tau(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t ** response_to_parameters){

  // We first need to create the atmosphere from the model and compute the observable:

  build_from_nodes(atmos_model);
  
  fp_t *** responses_depth_dependent;
  responses_depth_dependent = ft3dim(1,7,x3l,x3h,1,nlambda); 
  observable *o = obs_scalar_responses_tau(theta, phi, lambda, nlambda, responses_depth_dependent);

  io.msg(IOL_INFO, "atmosphere::computed observable from a model. \n");

  int N_parameters = atmos_model->get_N_nodes_total();
  int N_nodes_T = atmos_model->get_N_nodes_temp();
  int N_nodes_vt = atmos_model->get_N_nodes_vt();

  memset(response_to_parameters[1]+1, 0, N_parameters*nlambda*sizeof(fp_t));
  
  io.msg(IOL_INFO, "atmosphere::analytically computing responses to total of %d model parameters. \n", N_parameters);

  fp_t step = delta_T; // We still need step because we want to numerically study the derivative of the interpolation scheme. It does not matter a lot.
     
  for (int i=1;i<=N_nodes_T;++i){
    
    // First perturb to the one side
    fp_t * local_T_perturbations = new fp_t [x3h-x3l+1] - x3l;
    fp_t * local_Nt_perturbations = new fp_t [x3h-x3l+1] - x3l;

    atmos_model->perturb_node_value(i, step * 0.5);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      local_T_perturbations[x3i] = T[x1l][x2l][x3i];
      local_Nt_perturbations[x3i] = Nt[x1l][x2l][x3i];
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      local_T_perturbations[x3i] -= T[x1l][x2l][x3i];
      local_Nt_perturbations[x3i] -= Nt[x1l][x2l][x3i];
      local_T_perturbations[x3i] /= (1.0 * step);
      local_Nt_perturbations[x3i] /= (1.0 * step); 
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    build_from_nodes(atmos_model);
    

// ----------------------------------------------------------------------------------------------------------------------------------------------------

    // Now you need everything together:
    for (int x3i=x3l;x3i<=x3h;++x3i){
      //printf("%d %e %e \n", x3i, local_T_perturbations[x3i], local_Nt_perturbations[x3i]/Nt[x1l][x2l][x3i]);
      for (int l=1;l<=nlambda;++l)
      response_to_parameters[i][l] += responses_depth_dependent[1][x3i][l] * local_T_perturbations[x3i] + 
        responses_depth_dependent[2][x3i][l] * local_Nt_perturbations[x3i];
    }

    delete [](local_T_perturbations+x3l);
    delete [](local_Nt_perturbations+x3l);
    
  }

  // Then similar but for vt. However we do not have to build the atmosphere now

  step = delta_vt; // We still need step because we want to numerically study the derivative of the interpolation scheme. It does not matter a lot.
     
  for (int i=N_nodes_T+1;i<=N_nodes_T+N_nodes_vt;++i){
    
    // First perturb to the one side
    fp_t * local_vt_perturbations = new fp_t [x3h-x3l+1] - x3l;
    
    atmos_model->perturb_node_value(i, step * 0.5);
    build_from_nodes(atmos_model); // <------- Can be done better but cannot be bothered now. FIX! 

    for (int x3i=x3l;x3i<=x3h;++x3i)
      local_vt_perturbations[x3i] = Vt[x1l][x2l][x3i];
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      local_vt_perturbations[x3i] -= Vt[x1l][x2l][x3i];
      local_vt_perturbations[x3i] /= (1.0 * step);
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    build_from_nodes(atmos_model);
    
// ----------------------------------------------------------------------------------------------------------------------------------------------------

    // Now you need everything together:
    for (int x3i=x3l;x3i<=x3h;++x3i){
      //printf("%d %e %e \n", x3i, local_T_perturbations[x3i], local_Nt_perturbations[x3i]/Nt[x1l][x2l][x3i]);
      for (int l=1;l<=nlambda;++l)
      response_to_parameters[i][l] += responses_depth_dependent[3][x3i][l] * local_vt_perturbations[x3i];
        
    }

    delete [](local_vt_perturbations+x3l);
    
  }

  del_ft3dim(responses_depth_dependent,1,7,x3l,x3h,1,nlambda);
  return o;

 }

 // Same as the above except for polarized radiation (Zeeman mode so far):
 observable *atmosphere::obs_stokes_responses_to_nodes(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** response_to_parameters){

  // We first need to create the atmosphere from the model and compute the observable:
  build_from_nodes(atmos_model);
  
  // Now responses have an additional dimension
  fp_t **** responses_depth_dependent;
  responses_depth_dependent = ft4dim(1,7,x3l,x3h,1,nlambda,1,4); 
  observable *o = obs_stokes_responses(theta, phi, lambda, nlambda, responses_depth_dependent);

  io.msg(IOL_INFO, "atmosphere::computed polarized observable from a model. \n");

  int N_parameters = atmos_model->get_N_nodes_total();
  int N_nodes_T = atmos_model->get_N_nodes_temp();
  int N_nodes_vt = atmos_model->get_N_nodes_vt();
  int N_nodes_vs = atmos_model->get_N_nodes_vs();
  int N_nodes_B = atmos_model->get_N_nodes_B();
  int N_nodes_theta = atmos_model->get_N_nodes_theta();
  int N_nodes_phi = atmos_model->get_N_nodes_phi();

  memset(response_to_parameters[1][1]+1, 0, N_parameters*nlambda*4*sizeof(fp_t));
  
  io.msg(IOL_INFO, "atmosphere::analytically computing polarized responses to total of %d model parameters. \n", N_parameters);

  fp_t step = delta_T; // We still need step because we want to numerically study the derivative of the interpolation scheme. It does not matter a lot.
     
  for (int i=1;i<=N_nodes_T;++i){
    
    // First perturb to the one side
    fp_t * local_T_perturbations = new fp_t [x3h-x3l+1] - x3l;
    fp_t * local_Nt_perturbations = new fp_t [x3h-x3l+1] - x3l;

    atmos_model->perturb_node_value(i, step * 0.5);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      local_T_perturbations[x3i] = T[x1l][x2l][x3i];
      local_Nt_perturbations[x3i] = Nt[x1l][x2l][x3i];
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      local_T_perturbations[x3i] -= T[x1l][x2l][x3i];
      local_Nt_perturbations[x3i] -= Nt[x1l][x2l][x3i];
      local_T_perturbations[x3i] /= (1.0 * step);
      local_Nt_perturbations[x3i] /= (1.0 * step); 
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    build_from_nodes(atmos_model);
  
    // Now you need everything together:
    for (int x3i=x3l;x3i<=x3h;++x3i){
      //printf("%d %e %e \n", x3i, local_T_perturbations[x3i], local_Nt_perturbations[x3i]/Nt[x1l][x2l][x3i]);
      for (int l=1;l<=nlambda;++l)
        for (int s=1;s<=4;++s)
          response_to_parameters[i][l][s] += responses_depth_dependent[1][x3i][l][s] * local_T_perturbations[x3i] + 
            responses_depth_dependent[2][x3i][l][s] * local_Nt_perturbations[x3i];
    }

    delete [](local_T_perturbations+x3l);
    delete [](local_Nt_perturbations+x3l);
    
  }
  
  // Then for all other parameters it is much simpler as we do not have to build atmosphere, we can just interpolate.
  for (int i=N_nodes_T+1;i<=N_parameters;++i){

    if (i<=N_nodes_T+N_nodes_vt)
      step = delta_vt;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs)
      step = delta_vr;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs+N_nodes_B)
      step = delta_B;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs+N_nodes_B+N_nodes_theta)
      step = delta_angle;
    else 
      step = delta_angle;
    // First perturb to the one side
    fp_t * local_perturbations = new fp_t [x3h-x3l+1] - x3l;
    
    atmos_model->perturb_node_value(i, step * 0.5);
    interpolate_from_nodes(atmos_model); 

    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t B_mag = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i]+By[x1l][x2l][x3i]*By[x1l][x2l][x3i]+Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
      fp_t theta_B = acos(Bz[x1l][x2l][x3i]/B_mag);
      fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);
      local_perturbations[x3i] = Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B;
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    interpolate_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t B_mag = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i]+By[x1l][x2l][x3i]*By[x1l][x2l][x3i]+Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
      fp_t theta_B = acos(Bz[x1l][x2l][x3i]/B_mag);
      fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);
      local_perturbations[x3i] -= (Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B);
      local_perturbations[x3i] /= step;
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    interpolate_from_nodes(atmos_model);
  
    // Now you need everything together:
    int param = 0; // Which parameter are we accounting for? 
    if (i<=N_nodes_T+N_nodes_vt)
      param = 3;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs)
      param = 4;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs+N_nodes_B)
      param = 5;
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs+N_nodes_B+N_nodes_theta)
      param = 6;
    else 
      param = 7;
    for (int x3i=x3l;x3i<=x3h;++x3i){
      for (int l=1;l<=nlambda;++l)
        for (int s=1;s<=4;++s)
          response_to_parameters[i][l][s] += responses_depth_dependent[param][x3i][l][s] * local_perturbations[x3i];
    }
    delete [](local_perturbations+x3l);
  }

  del_ft4dim(responses_depth_dependent,1,7,x3l,x3h,1,nlambda,1,4);
  return o;

 }

 // Same as the above except for polarized radiation (Zeeman mode so far):
 observable *atmosphere::obs_stokes_responses_to_nodes_new(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** response_to_parameters){

  // We first need to create the atmosphere from the model and compute the observable:
  
  build_from_nodes(atmos_model);
  
  int N_parameters = atmos_model->get_N_nodes_total();
  int N_nodes_T = atmos_model->get_N_nodes_temp();
  int N_nodes_vt = atmos_model->get_N_nodes_vt();
  int N_nodes_vs = atmos_model->get_N_nodes_vs();
  int N_nodes_B = atmos_model->get_N_nodes_B();
  int N_nodes_theta = atmos_model->get_N_nodes_theta();
  int N_nodes_phi = atmos_model->get_N_nodes_phi();

  memset(response_to_parameters[1][1]+1, 0, N_parameters*nlambda*4*sizeof(fp_t));

  int N_depths = x3h-x3l+1;
  fp_t *** resp_atm_to_parameters = ft3dim(1,N_parameters,1,7,1,N_depths);
  memset(resp_atm_to_parameters[1][1]+1,0,N_parameters*7*N_depths*sizeof(fp_t));
  
  fp_t step = delta_T; // We still need step because we want to numerically study the derivative of the interpolation scheme. It does not matter a lot.
     
  for (int i=1;i<=N_nodes_T;++i){
    
    atmos_model->perturb_node_value(i, step * 0.5);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      resp_atm_to_parameters[i][1][x3i-x3l+1] = T[x1l][x2l][x3i];
      resp_atm_to_parameters[i][2][x3i-x3l+1] = Nt[x1l][x2l][x3i];
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    build_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      resp_atm_to_parameters[i][1][x3i-x3l+1] -= T[x1l][x2l][x3i];
      resp_atm_to_parameters[i][2][x3i-x3l+1] -= Nt[x1l][x2l][x3i];
      resp_atm_to_parameters[i][1][x3i-x3l+1] /= (1.0 * step);
      resp_atm_to_parameters[i][2][x3i-x3l+1] /= (1.0 * step); 
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    build_from_nodes(atmos_model);
  }
  
  int param = 0;
  // Then for all other parameters it is much simpler as we do not have to build atmosphere, we can just interpolate.
  for (int i=N_nodes_T+1;i<=N_parameters;++i){

    if (i<=N_nodes_T+N_nodes_vt){
      step = delta_vt;
      param = 3;
    }
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs){
      step = delta_vr;
      param = 4;
    }
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs+N_nodes_B){
      step = delta_B;
      param = 5;
    }
    else if (i<=N_nodes_T+N_nodes_vt+N_nodes_vs+N_nodes_B+N_nodes_theta){
      step = delta_angle;
      param = 6;
    }
    else {
      step = delta_angle;
      param = 7;
    }
    // First perturb to the one side
    fp_t * local_perturbations = new fp_t [x3h-x3l+1] - x3l;
    
    atmos_model->perturb_node_value(i, step * 0.5);
    interpolate_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t B_mag = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i]+By[x1l][x2l][x3i]*By[x1l][x2l][x3i]+Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
      fp_t theta_B = acos(Bz[x1l][x2l][x3i]/B_mag);
      fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);
      if (B_mag < 1E-3)
        B_mag = theta_B = phi_B = 0.0;
      resp_atm_to_parameters[i][param][x3i-x3l+1] = Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B;
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    interpolate_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t B_mag = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i]+By[x1l][x2l][x3i]*By[x1l][x2l][x3i]+Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
      fp_t theta_B = acos(Bz[x1l][x2l][x3i]/B_mag);
      fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);
      if (B_mag < 1E-3)
        B_mag = theta_B = phi_B = 0.0;
      resp_atm_to_parameters[i][param][x3i-x3l+1] -= (Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B);
      resp_atm_to_parameters[i][param][x3i-x3l+1] /= step;
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    interpolate_from_nodes(atmos_model);
    delete[](local_perturbations+x3l);
  
  }

  // Now we start the responses part.
  atmos_model->set_response_to_parameters(resp_atm_to_parameters,x3h-x3l+1);
  observable *o = obs_stokes_responses(theta, phi, lambda, nlambda, response_to_parameters,atmos_model,0);
    del_ft3dim(resp_atm_to_parameters,1,N_parameters,1,7,1,N_depths);
  
  return o;
  
 }

 // ===============================================================================================

int atmosphere::polish_extreme_values(){
  // This function just checks for bad values of physical parameters and puts them
  // back in the boundaries
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        if (T[x1i][x2i][x3i] < 3400.0) T[x1i][x2i][x3i] = 3400.0;
        if (T[x1i][x2i][x3i] > 12000.0) T[x1i][x2i][x3i] = 12000.0;
        if (Vz[x1i][x2i][x3i] < -2E6) Vz[x1i][x2i][x3i] = -2E6;
        if (Vz[x1i][x2i][x3i] > 2E6) Vz[x1i][x2i][x3i] = 2E6;

  }
  return 0;
}

 // ===============================================================================================

 /////////////////////
 // OLD:::::

 observable * atmosphere::scalar_lm_fit(observable * spectrum_to_fit, fp_t theta, fp_t phi, fp_t * lambda, int nlambda){};
