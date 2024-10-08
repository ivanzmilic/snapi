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

observable * atmosphere::stokes_lm_fit(observable * obs_to_fit, fp_t theta, fp_t phi, model * model_to_fit){

  // ------------------------------------------------------------------------------------------------------------------

  // Perform LM fitting and return the best fitting spectrum
  // observable * spectrum_to_fit : observable object to fit, contains stokes vector in physical units and the 
  //                                wavelenght grid 
  // fp_t theta, phi              : angles, theta for off-center observations, phi rarely used                                    
  // model * model_to_fit         : model object that will be used to construct the atmosphere, to store the solution:
  //                              : the values of parameters at the nodes will be packed. 
  // ------------------------------------------------------------------------------------------------------------------

  // First extract the spectrum, number of wavelengths and the wavelength grid 
  // from the observation:
  
  fp_t ** S_to_fit = obs_to_fit->get_S(1,1);
  int nlambda = obs_to_fit->get_n_lambda();
  fp_t * lambda = obs_to_fit->get_lambda();

  // Grid has to be set to be tau (in principle this can be changed):
  set_grid(1);
  
  // Set initial value of Levenberg-Marquardt parameter
  // Here we also hardcode how big changes we make in the LM jumps
  fp_t lm_parameter = obs_to_fit->get_start_lambda();
  fp_t lm_multiplicator = sqrt(10.0); // This is so technical, makes sense to have it hard-coded

  
  // Some fitting related parameters:
  fp_t metric = 0.0;
  // Iteration counter:
  int iter = 0;
  // Get what is maximum number of iterations: 
  int MAX_ITER = obs_to_fit->get_no_iterations();
    
  // At which chi-squared to we stop iterating:
  fp_t stopping_chisq = obs_to_fit->get_stopping_chisq();

  // Auxiliary variables that tell us when to terminate the iteration:
  fp_t * chi_to_track = 0;
  int n_chi_to_track = 0;
  int corrected = 1; // At the start we assume the state of the atmosphere has been corrected -> RF has to be recalculated
  int to_break = 0;

  // weights for Stokes parameters. They act like noise, linearly affect RFs and quadratically chisquared
  fp_t * ws = obs_to_fit->get_w_stokes();
  // wavelenght weights / mask 
  fp_t * wl = obs_to_fit->get_mask();

  // other fitting parameters
  fp_t spectral_broadening = obs_to_fit->get_spectral_broadening(); // This needs to be changed to that we can 
                                                                         // have different PSFs for different observables
  
  int n_spsf = obs_to_fit->get_n_spsf();
  fp_t * spsf = 0;
  if (n_spsf)
    spsf = obs_to_fit->get_spsf();
  
  fp_t qs_level = obs_to_fit->get_synth_qs();
 
  fp_t noise_level = 1E-3*qs_level; // The magnitude does not really matter. But keep it at something realistic.
  fp_t *noise_scaling = new fp_t [nlambda]-1; // wavelength dependent noise
  for (int l=1;l<=nlambda;++l)
    noise_scaling[l] = sqrt(S_to_fit[1][l]/S_to_fit[1][1]);
  fp_t *noise = new fp_t[nlambda]-1;
  for (int l=1;l<=nlambda;++l)
    noise[l] = noise_level * noise_scaling[l];
  
  observable * current_obs;
  fp_t *** derivatives_to_parameters; // [i, s, l], dS_{s,l} / d M_i
  fp_t ** S_current; // dS_{s,l}

  int N_parameters = model_to_fit->get_N_nodes_total();
  model_to_fit->bracket_parameter_values(); // This puts them in physical limits

  for (iter=1;iter<=MAX_ITER;++iter){

    if (corrected){      
      
      // These quantities are only re-computed if the model has been modified:    
      derivatives_to_parameters = ft3dim(1,N_parameters,1,nlambda,1,4);     
      memset(derivatives_to_parameters[1][1]+1,0,N_parameters*nlambda*4*sizeof(fp_t));

      // Calculate the spectrum and the responses:     
      current_obs = obs_stokes_responses_to_nodes(model_to_fit, theta, phi, lambda, nlambda, derivatives_to_parameters, 0); 
      
      // Apply spectral broadening if necessary:
      if (spectral_broadening){
        //fprintf(stderr,"atmosphere::stokes_lm_fit: convolving with gaussian psf");
        current_obs->spectral_convolve(spectral_broadening,1,1);
        convolve_response_with_gauss(derivatives_to_parameters,lambda,N_parameters,nlambda,spectral_broadening);   
      } 
      // Then check if you need F-P too (in principle both are allowed).
      else if (n_spsf){ 
        //fprintf(stderr,"atmosphere::stokes_lm_fit: convolving with given psf");
        current_obs->psf_convolve(n_spsf,spsf,1,1);
        convolve_response_with_psf(derivatives_to_parameters,lambda,N_parameters,nlambda,n_spsf,spsf);   
      }

      scale_rf(derivatives_to_parameters,model_to_fit,nlambda,N_parameters,ws,wl);
      S_current = current_obs->get_S(1,1);
    }

    fp_t * residual = calc_residual(S_to_fit,S_current,nlambda, ws, wl);
    metric = calc_chisq(S_to_fit, S_current, nlambda, ws, wl);
    //fprintf(stderr, "atmosphere::stokes_lm_fit: current chi-squared =  %e, \n", metric/4.0/nlambda/noise_level/noise_level);
    
    if (metric < stopping_chisq)
      to_break = 1;
    
    fp_t ** J = ft2dim(1,4*nlambda,1,N_parameters);
    for (int i=1;i<=N_parameters;++i) 
      for (int l=1;l<=nlambda;++l) 
        for (int s=1;s<=4;++s){
          J[(l-1)*4+s][i] = derivatives_to_parameters[i][l][s];
    }
    
    fp_t ** J_transpose = transpose(J,4*nlambda,N_parameters);
    fp_t ** JTJ = multiply_with_transpose(J, 4*nlambda, N_parameters);
    fp_t * rhs = multiply_vector(J_transpose, residual, N_parameters, 4*nlambda);    
    regularize_hessian(JTJ,rhs,model_to_fit);

    for (int i=1;i<=N_parameters;++i) JTJ[i][i] *= (lm_parameter + 1.0);
    // Now correct
    //fp_t * correction = solve(JTJ, rhs, 1, N_parameters); 
    fp_t ** correction = solve_svd(JTJ, &rhs-1, 1, N_parameters, 1);

    scale_corrections(correction[1],model_to_fit,N_parameters);

    // Reset the Hessian to it's original value in case it stays the same and our LM parameter changes:
    for (int i=1;i<=N_parameters;++i) JTJ[i][i] /= (lm_parameter + 1.0);

    // Apply the correction:
    model * test_model = clone(model_to_fit);
    test_model->correct(correction[1]);
    test_model->bracket_parameter_values(); // polish
    
    // Check if the corrected model is better:
    build_from_nodes(test_model);
    observable *reference_obs = forward_evaluate(theta,phi,lambda,nlambda,0,qs_level,spectral_broadening,n_spsf, spsf); 
    fp_t ** S_reference = reference_obs->get_S(1,1);
    fp_t metric_reference = calc_chisq(S_to_fit,S_reference,nlambda,ws,wl);
    test_model->print();
    fprintf(stderr, "atmosphere::stokes_lm_fit: test chi-squared = %e, \n", metric/4.0/nlambda/noise_level/noise_level);

    if (metric_reference < metric){ // If the solution is better:

      model_to_fit->cpy_values_from(test_model);
      lm_parameter /= lm_multiplicator;
      if (lm_parameter <= 1E-5) lm_parameter = 1E-5; // not sure if necessary
      corrected=1;
      
      // Tracked values of chisquared
      chi_to_track = add_to_1d_array(chi_to_track,n_chi_to_track,metric);
      if (n_chi_to_track >=3){
        fp_t change = fabs((chi_to_track[n_chi_to_track-2] - chi_to_track[n_chi_to_track-1]) / chi_to_track[n_chi_to_track-1]);
        if (change < DELTA) // This implies that our change is very small
          to_break = 1;
      }
    }
    else{ // Else the solution is worse and we are going toward the gradient decent:
      lm_parameter *= lm_multiplicator;
      corrected = 0;
    }

    //model_to_fit->print();

    if(corrected || to_break || iter==MAX_ITER){     
      del_ft3dim(derivatives_to_parameters,1,N_parameters,1,nlambda,1,4);
      delete current_obs;
      del_ft2dim(S_current,1,4,1,nlambda);
    }

    delete[](residual+1);
    delete test_model;
    del_ft2dim(J_transpose,1,N_parameters,1,4*nlambda);
    del_ft2dim(JTJ,1,N_parameters,1,N_parameters);
    del_ft2dim(J,1,nlambda*4,1,N_parameters);
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    delete [](rhs+1);
    del_ft2dim(correction,1,1,1,N_parameters);
    //delete [](correction+1);
    metric = 0.0;

    if (to_break)
      break;
  }
  
  // Clean-up these temporary ones
  del_ft2dim(S_to_fit,1,4,1,nlambda);
  delete[](lambda+1);
  if (chi_to_track)
    delete[]chi_to_track;

  lambda = obs_to_fit->get_lambda();
  nlambda = obs_to_fit->get_n_lambda();
  build_from_nodes(model_to_fit);
  observable *obs_to_return = forward_evaluate(theta,phi,lambda,nlambda,0,qs_level,spectral_broadening,n_spsf, spsf);
  model_to_fit->polish_angles();   
  
  delete[](lambda+1);
  delete[]ws;
  delete[](wl+1);
  delete[](noise+1);
  if (n_spsf)
    delete[](spsf+1);
  
  return obs_to_return;
}

fp_t * atmosphere::calc_residual(fp_t ** S_to_fit, fp_t ** S_current, int nlambda, fp_t * ws, fp_t * wl){

  fp_t * residual = new fp_t[4*nlambda]-1;
  for (int l=1;l<=nlambda;++l)
    for (int s=1;s<=4;++s){
      residual[(l-1)*4+s] = (S_to_fit[s][l] - S_current[s][l]) * ws[s-1] * wl[l];
  }
  return residual;
}

fp_t atmosphere::calc_chisq(fp_t ** S_to_fit, fp_t ** S_current, int nlambda, fp_t * ws, fp_t * wl){
  fp_t chisq = 0.0;
  for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=4;++s){
        fp_t residual_temp = (S_to_fit[s][l] - S_current[s][l]) * ws[s-1] * wl[l];
        chisq += residual_temp * residual_temp;
        
  }
  return chisq;
}

int atmosphere::scale_rf(fp_t *** derivatives_to_parameters, model* model_to_fit, int nlambda, int N_parameters, fp_t * w_stokes, fp_t * l_mask){

  fp_t scales[6] ={T_scale,vt_scale,vr_scale,B_scale,1.0,1.0};
  for (int i=1;i<=N_parameters;++i){
    int index=model_to_fit->which_parameter(i);
    for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=4;++s)
        derivatives_to_parameters[i][l][s] *= scales[index-1] * w_stokes[s-1] * l_mask[l];
  }
  return 0;
}

int atmosphere::scale_corrections(fp_t * corrections, model* model_to_fit, int N_parameters){
  fp_t scales[6] ={T_scale,vt_scale,vr_scale,B_scale,1.0,1.0};
  for (int i=1;i<=N_parameters;++i){
    int index=model_to_fit->which_parameter(i);
    corrections[i] *= scales[index-1];
  }
  return 0; 
}

int atmosphere::look_for_best_lambda(fp_t &lm_parameter, fp_t ** JTJ, int N_parameters,
  fp_t * rhs, model * model_to_fit, fp_t theta, fp_t phi, fp_t * lambda, int nlambda, fp_t scattered_light,
  fp_t qs_level, fp_t spectral_broadening, fp_t ** S_to_fit, int n_stokes_to_fit, int * stokes_to_fit,
  fp_t * ws, fp_t * noise, fp_t metric_old, fp_t * wl){
  
  fp_t lm_multiplicator=10.0;
  model * model_test;
  fp_t metric_prev = metric_old;

  int MAX_ITER = 10; // Maximum number of explorations.  
  int iter = 0;
  for (iter=1;iter<=MAX_ITER;++iter){
    
    lm_parameter /= lm_multiplicator;
    for (int i=1;i<=N_parameters;++i) JTJ[i][i] *= (1.0+lm_parameter);

    fp_t * correction = solve(JTJ, rhs, 1, N_parameters);
    scale_corrections(correction,model_to_fit,N_parameters);
    // Reset to original Hessian
    for (int i=1;i<=N_parameters;++i)
      JTJ[i][i] /= (1.0+lm_parameter);
  
    // Apply the correction:
    model_test = clone(model_to_fit); // We want to correct original model.
    model_test->correct(correction);
    model_test->bracket_parameter_values();
    build_from_nodes(model_test);
    // Compare again:
    observable *reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening,0,0); 
    fp_t ** S_reference = reference_obs->get_S(1,1);
    fp_t * residual_test = calc_residual(S_to_fit,S_reference,nlambda, ws, wl);
    fp_t metric_reference = calc_chisq(S_to_fit,S_reference,nlambda, ws, wl);
    delete[](residual_test+1);
    delete[](correction+1);
    del_ft2dim(S_reference,1,4,1,nlambda);
    delete reference_obs;
    //fprintf(stderr,"Correction attempt %d, m_ref = %e m_prev = %e, lambda = %e \n",iter, metric_reference,metric_prev, lm_parameter);
    if (metric_reference > metric_prev  || lm_parameter <= 1E-3/lm_multiplicator) { // If correction is bad, increase lm to previous value and break
      
      lm_parameter *= lm_multiplicator;
      // Delete the test model and update original model (again - awkward) according to the last value of lambda
      delete model_test;
      for (int i=1;i<=N_parameters;++i)
        JTJ[i][i] *= (1.0+lm_parameter);
      correction = solve(JTJ, rhs, 1, N_parameters);
      scale_corrections(correction,model_to_fit,N_parameters);
      model_to_fit->correct(correction);
      model_to_fit->bracket_parameter_values();
      delete[](correction+1);
      break; 
    }
    else if (iter==MAX_ITER){ //Else if we reached the last iteration leave lm_parameter as is 
                              // and copy the test model to the original model
      model_to_fit->cpy_values_from(model_test);
      delete model_test;
      break;
    }
    else { // Else everything is good, just update the metric reference and go to the next iteration
      metric_prev = metric_reference;
      delete model_test;
    }
  }
  //fprintf(stderr,"Found best lm = %f in %d attempts \n",lm_parameter,iter);
  return 0;
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
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening,0,0);
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
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening,0,0);
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
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening,0,0);
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
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening,0,0);
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
    reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening,0,0);
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
    fp_t scattered_light, fp_t qs, fp_t spectral_broadening, int n_spsf, fp_t * spsf){

  // proceduralized version which performs, in turn
  // 1) Stokes synthesis from current atmosphere
  // 2) Adds scattered light from the quiet Sun
  // 3) Convolves with instrumental profile

  observable *reference_obs = obs_stokes(theta, phi, lambda, nlambda);
  reference_obs->add_scattered_light(scattered_light,qs);
  if (spectral_broadening)
    reference_obs->spectral_convolve(spectral_broadening,1,1);
    reference_obs->set_n_spsf(0);
  if (n_spsf){
    reference_obs->set_n_spsf(n_spsf);
    reference_obs->set_spsf(spsf);
    reference_obs->psf_convolve(n_spsf, spsf, 1, 1);
  }
  
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
  fp_t V_scale = vt_scale;
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

  observable *atmosphere::obs_stokes_num_responses_to_nodes(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** response_to_parameters, fp_t filter_width){

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
   return 0; 
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

  return 0;

 }

// Same as the above except for polarized radiation (Zeeman mode so far):
observable *atmosphere::obs_stokes_responses_to_nodes(model * atmos_model, fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** response_to_parameters, fp_t filter_width){

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
  atmos_model->set_response_to_parameters(resp_atm_to_parameters,x3h-x3l+1);
  del_ft3dim(resp_atm_to_parameters,1,N_parameters,1,7,1,N_depths);
  
  // If there is a filter there is some gymnastics:
  // 1) Create a fine wavelenght grid and calculate spectrum and the responses on it
  // 2) Interpolate back to the original grid where exact filter positions are given. 
  // 3) Return the observable and responses on the original grid to the fitting method
  if (filter_width){
    int nlambda_fine = int((lambda[nlambda]+5.0*filter_width - lambda[1] - 5.0*filter_width)*2.0/filter_width);
    fp_t * lambda_fine = new fp_t [nlambda_fine]-1;
    for (int l=1;l<=nlambda_fine;++l){
      lambda_fine[l] = lambda[1]-5.0*filter_width + filter_width/2.0 * (l-1);
    }
    fp_t *** response_to_parameters_fine = ft3dim(1,N_parameters,1,nlambda_fine,1,4);
    observable * o_fine = obs_stokes_responses(theta, phi, lambda_fine, nlambda_fine, response_to_parameters_fine,atmos_model,0);
    o_fine->spectral_convolve(filter_width,1,1);
    convolve_response_with_gauss(response_to_parameters_fine,lambda_fine,N_parameters,nlambda_fine,filter_width);   
    fp_t ** S_fine = o_fine->get_S(1,1);
    fp_t **** S = ft4dim(1,1,1,1,1,4,1,nlambda);
    for (int s=1;s<=4;++s)
      for (int l=1;l<=nlambda;++l)
        S[1][1][s][l] = interpol_1d(S_fine[s],lambda_fine,nlambda_fine,lambda[l]);
    for (int i=1;i<=N_parameters;++i){
      fp_t ** response_temp = transpose(response_to_parameters_fine[i],nlambda,4);
      for (int s=1;s<=4;++s)
        for (int l=1;l<=nlambda;++l)
          response_to_parameters[i][l][s] = interpol_1d(response_temp[s],lambda_fine,nlambda_fine,lambda[l]);
      del_ft2dim(response_temp,1,4,1,nlambda);
    }
    del_ft2dim(S_fine,1,4,1,nlambda_fine);
    del_ft3dim(response_to_parameters_fine,1,N_parameters,1,nlambda_fine,1,4);
    delete[](lambda_fine+1);
    observable * o = new observable(1,1,4,nlambda); 
    o->set(S);
    o->set_lambda(lambda);
    del_ft4dim(S,1,1,1,1,1,4,1,nlambda);
    delete o_fine;
    return o;
  }
  else{
    observable *o = obs_stokes_responses(theta, phi, lambda, nlambda, response_to_parameters,atmos_model,0);
    return o;
  }
  return 0;
}

 // ===============================================================================================

int atmosphere::polish_extreme_values(){
  // This function just checks for bad values of physical parameters and puts them
  // back in the boundaries
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        if (T[x1i][x2i][x3i] < 3400.0) T[x1i][x2i][x3i] = 3400.0;
        if (T[x1i][x2i][x3i] > 20000.0) T[x1i][x2i][x3i] = 20000.0;
        if (Vz[x1i][x2i][x3i] < -2E6) Vz[x1i][x2i][x3i] = -2E6;
        if (Vz[x1i][x2i][x3i] > 2E6) Vz[x1i][x2i][x3i] = 2E6;

  }
  return 0;
}

void atmosphere::regularize_hessian(fp_t ** JTJ, fp_t * rhs, model * model_to_fit){
  for (int i=1;i<=6;++i)
    regularize_parameter(JTJ, rhs, model_to_fit,i);
}

void atmosphere::regularize_parameter(fp_t ** JTJ, fp_t * rhs, model * model_to_fit, int param){

    int reg_type = 0;
    switch(param){
      case 1 : reg_type = model_to_fit->get_temp_reg_type(); break;
      case 2 : reg_type = model_to_fit->get_vt_reg_type(); break;
      case 3 : reg_type = model_to_fit->get_vs_reg_type(); break;
      case 4 : reg_type = model_to_fit->get_B_reg_type(); break;
      case 5 : reg_type = model_to_fit->get_theta_reg_type(); break;
      case 6 : reg_type = model_to_fit->get_phi_reg_type(); break;
    }

    //printf("Param = %d Regtype = %d \n",param,reg_type);

    if (reg_type){

      // Load the strength of the regularization
      fp_t alpha = 0.0;
      int from=1; int to=0;
      switch(param){
        case 1 : alpha = model_to_fit->get_temp_reg_alpha(); break; 
        case 2 : alpha = model_to_fit->get_vt_reg_alpha(); break;
        case 3 : alpha = model_to_fit->get_vs_reg_alpha(); break;
        case 4 : alpha = model_to_fit->get_B_reg_alpha(); break;
        case 5 : alpha = model_to_fit->get_theta_reg_alpha(); break;
        case 6 : alpha = model_to_fit->get_phi_reg_alpha(); break;
      }
      int N_parameters = model_to_fit->get_N_nodes_total();
      for (int i=1;i<=N_parameters;++i)
        if (model_to_fit->which_parameter(i) == param){
          from = i;
          break;
      }
      for (int i=from;i<=N_parameters;++i)
        if (model_to_fit->which_parameter(i)!= param){
          to = i-1;
          break;
      }
      //fprintf(stderr,"Param = %d From = %d to = %d alpha = %e \n",param,from,to,alpha);
  
      // Regularization matrix. It's hessian is added to the hessian of the chisq
      fp_t ** Reg = ft2dim(1,N_parameters,1,N_parameters);
      memset(Reg[1]+1,0,N_parameters*N_parameters*sizeof(fp_t));
      
      // Here we want to add a regularization matrix:
      
      if (reg_type==2){ // Tikhonov
        for (int i=from+1;i<=to;++i){
          Reg[i][i] = 1.0;
          Reg[i][i-1] = -1.0;
        }
        Reg[from][from] = -1.0;
        Reg[from][from+1] = 1.0;
      }
      else if (reg_type==1){ // Departure from zero
        for (int i=from;i<=to;++i)
          Reg[i][i] = 1.0;
      }
      
      // Calculate Hessian of the regularization
      fp_t ** RegTReg = multiply_with_transpose(Reg,N_parameters,N_parameters);
      for (int i=1;i<=N_parameters;++i)
        for (int j=1;j<=N_parameters;++j){
          RegTReg[i][j] *= alpha*alpha;
          JTJ[i][j] += RegTReg[i][j];
      }
      del_ft2dim(RegTReg,1,N_parameters,1,N_parameters);
      
      // 'Residual' of the regularization:
      fp_t ** Reg_transpose = transpose(Reg,N_parameters,N_parameters);
      fp_t * reg_residual = new fp_t [N_parameters]-1;
      memset(reg_residual+1,0,N_parameters*sizeof(fp_t));

      fp_t * quantity;
      switch(param){
        case 1 : quantity = model_to_fit->get_temp_nodes_temp(); break;
        case 2 : quantity = model_to_fit->get_vt_nodes_vt(); break;
        case 3 : quantity = model_to_fit->get_vs_nodes_vs(); break;
        case 4 : quantity = model_to_fit->get_B_nodes_B(); break;
        case 5 : quantity = model_to_fit->get_theta_nodes_theta(); break; 
        case 6 : quantity = model_to_fit->get_phi_nodes_phi(); break;
      }
      //for (int i=from;i<=to;++i)
      //  printf("%d %e \n",i,quantity[i-from+1]);

      if (reg_type==2){ // Tikhonov
        for (int i=from+1;i<=to;++i)
          reg_residual[i] = (quantity[i+1-from]-quantity[i-from])*alpha;
        reg_residual[from] = (quantity[2] - quantity[1])*alpha;     
      }
      else if (reg_type==1){ // Departures from zero
        for (int i=from;i<=to;++i)
          reg_residual[i] = quantity[i+1-from]*alpha;
      }
      fp_t * reg_rhs = multiply_vector(Reg_transpose,reg_residual,N_parameters,N_parameters);
      //for (int i=1;i<=N_parameters;++i)
      //  rhs[i]-=reg_rhs[i];

      delete [](quantity+1);
      delete [](reg_residual+1);
      delete [](reg_rhs+1);
      del_ft2dim(Reg_transpose,1,N_parameters,1,N_parameters);
      del_ft2dim(Reg,1,N_parameters,1,N_parameters);
    }
}

// ===========================================================================================================
// ~~~ OBSOLETE ~~~

observable * atmosphere::stokes_lm_fit_old(observable * spectrum_to_fit, fp_t theta, fp_t phi, model * model_to_fit){
  return 0;
}
  /*
  //clock_t start = clock();

  // Perform LM fitting and return the best fitting spectrum
  
  // First extract the spectrum, number of wavelengths and the wavelength grid 
  // from the observation:
  fp_t ** S_to_fit = spectrum_to_fit->get_S_to_fit(1,1);
  int nlambda = spectrum_to_fit->get_n_lambda_to_fit();
  fp_t * lambda = spectrum_to_fit->get_lambda_to_fit();

  // Grid has to be set to be tau at the momen:
  set_grid(1);
  
  // Set initial value of Levenberg-Marquardt parameter
  // Here we also hardcode how big changes we make in the LM jumps
  fp_t lm_parameter = spectrum_to_fit->get_start_lambda();
  fp_t lm_multiplicator = 10.0;

  
  // Some fitting related parameters:
  fp_t metric = 0.0;
  // Iteration counter:
  int iter = 0;
  // Whether we use Jaime's method to look for optimum starting lambda:
  int search_for_optimum_lambda = 0;
  // Get what is maximum number of iterations: 
  int MAX_ITER = spectrum_to_fit->get_no_iterations();
  
  // If maximum number of iterations specified is negative, it is a crude command 
  // for asking the code to search for the optimal starting lambda:
  if (MAX_ITER < 0){
    MAX_ITER = -MAX_ITER;
    search_for_optimum_lambda = 1;
  }
  
  // At which chi-squared to we stop iterating:
  fp_t stopping_chisq = spectrum_to_fit->get_stopping_chisq();

  // Auxiliary variables that tell us when to terminate the iteration:
  fp_t * chi_to_track = 0;
  int n_chi_to_track = 0;
  int corrected = 1; // At the start we assume the state of the atmosphere has been corrected
  int to_break = 0;

  // weights for Stokes parameters. They enter like this in response scaling, and 
  // quadratically in chi_sq. basically they reduce the noise  
  fp_t * ws = spectrum_to_fit->get_w_stokes();

  // other fitting parameters
  fp_t scattered_light = spectrum_to_fit->get_scattered_light();
  fp_t spectral_broadening = spectrum_to_fit->get_spectral_broadening();
  fp_t qs_level = spectrum_to_fit->get_synth_qs();
  
  // A complicated piece of code to isolate what we want to fit
  int n_stokes_to_fit = 0; 
  for (int s=0;s<4;++s) 
    if (ws[s]) 
      ++n_stokes_to_fit;
  int stokes_to_fit[n_stokes_to_fit];
  int counter = 0;
  for (int s=0;s<4;++s) if (ws[s]){
    stokes_to_fit[counter] = s+1;
    ++counter;
  }

  fp_t noise_level = 1E-2*S_to_fit[1][1]; // The magnitude does not really matter.
  fp_t *noise_scaling = new fp_t [nlambda]-1; // wavelength dependent noise
  for (int l=1;l<=nlambda;++l)
    noise_scaling[l] = sqrt(S_to_fit[1][l]/S_to_fit[1][1]);
  fp_t *noise = new fp_t[nlambda]-1;
  for (int l=1;l<=nlambda;++l)
    noise[l] = noise_level * noise_scaling[l];

  observable * current_obs;
  fp_t *** derivatives_to_parameters;
  fp_t ** S_current;

  int N_parameters = model_to_fit->get_N_nodes_total();
  model_to_fit->bracket_parameter_values();

  io.msg(IOL_INFO, "atmosphere::stokes_lm_fit : entering iterative procedure\n");

  //clock_t point1 = clock();
  //fprintf(stderr,"Time needed to set-up: %e \n", (double)(point1-start)/CLOCKS_PER_SEC);

  for (iter=1;iter<=MAX_ITER;++iter){

    if (corrected){      
      
      // These quantities are only re-computed if the model has been modified:    
      derivatives_to_parameters = ft3dim(1,N_parameters,1,nlambda,1,4);     
      memset(derivatives_to_parameters[1][1]+1,0,N_parameters*nlambda*4*sizeof(fp_t));

      //fprintf(stderr,"We are starting from the following model: \n");
      //model_to_fit->print();
      
      // Calculate the spectrum and the responses and apply degradation to it:      
      current_obs = obs_stokes_responses_to_nodes(model_to_fit, theta, phi, lambda, nlambda, derivatives_to_parameters, 0); 
      
      current_obs->add_scattered_light(scattered_light,qs_level);

      //clock_t point2 = clock();

      //fprintf(stderr,"Time needed to calculate responses: %e \n", (double)(point2-point1)/CLOCKS_PER_SEC);      
      
      if (spectral_broadening){
        current_obs->spectral_convolve(spectral_broadening,1,1);
        convolve_response_with_gauss(derivatives_to_parameters,lambda,N_parameters,nlambda,spectral_broadening);   
      }

      //clock_t point3 = clock();

      //fprintf(stderr,"Time needed to convolve: %e \n", (double)(point3-point2)/CLOCKS_PER_SEC);      
      
      scale_rf(derivatives_to_parameters,model_to_fit,nlambda,N_parameters,ws,noise_scaling);
      S_current = current_obs->get_S(1,1);
    }

    //clock_t point3 = clock();
    
    fp_t * residual = calc_residual(S_to_fit,S_current,nlambda,n_stokes_to_fit,stokes_to_fit, ws);

    //for (int l=1; l<=nlambda; ++l)
    //  fprintf(stderr, "%e %e %e \n", lambda[l], S_to_fit[1][l], S_current[1][l]);
    metric = calc_chisq(nlambda, n_stokes_to_fit, stokes_to_fit, residual, noise, ws);
    if (metric < stopping_chisq)
      to_break = 1;
    fp_t ** J = ft2dim(1,n_stokes_to_fit*nlambda,1,N_parameters);
    for (int i=1;i<=N_parameters;++i) 
      for (int l=1;l<=nlambda;++l) 
        for (int s=1;s<=n_stokes_to_fit;++s){
          int stf = stokes_to_fit[s-1];
          J[(l-1)*n_stokes_to_fit+s][i] = derivatives_to_parameters[i][l][stf];
    }
    //for (int i=1; i<=N_parameters;++i){
    //  for (int ii=1; ii<=N_parameters; ++ii)
    //    fprintf(stderr, "%e ", J[i][ii]);
    //  fprintf(stderr,"\n");
    //}


    fp_t ** J_transpose = transpose(J,n_stokes_to_fit*nlambda,N_parameters);
    fp_t ** JTJ = multiply_with_transpose(J, n_stokes_to_fit*nlambda, N_parameters);
    fp_t * rhs = multiply_vector(J_transpose, residual, N_parameters, n_stokes_to_fit*nlambda);    
    regularize_hessian(JTJ,rhs,model_to_fit);

    for (int i=1;i<=N_parameters;++i) JTJ[i][i] *= (lm_parameter + 1.0);
    // Now correct
    fp_t * correction = solve(JTJ, rhs, 1, N_parameters); 
    scale_corrections(correction,model_to_fit,N_parameters);
    for (int i=1;i<=N_parameters;++i) JTJ[i][i] /= (lm_parameter + 1.0);
  
    // Apply the correction:
    model * test_model = clone(model_to_fit);
    test_model->correct(correction);
    test_model->bracket_parameter_values();

    //fprintf(stderr,"Model corrected. \n");
    //test_model->print();

    build_from_nodes(test_model);
    // Compare again:
    observable *reference_obs = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening); 
    fp_t ** S_reference = reference_obs->get_S(1,1);
    fp_t * residual_test = calc_residual(S_to_fit,S_reference,nlambda,n_stokes_to_fit,stokes_to_fit, ws);
    fp_t metric_reference = calc_chisq(nlambda, n_stokes_to_fit, stokes_to_fit, residual_test, noise, ws);
    delete[](residual_test+1);
    
    if (metric_reference < metric){

      //clock_t point3 = clock();
      
      // Everything is ok, and we can decrease lm_parameter:
      if (iter == 1 && search_for_optimum_lambda)
        look_for_best_lambda(lm_parameter, JTJ, N_parameters,
          rhs, model_to_fit, theta, phi, lambda, nlambda, scattered_light,
          qs_level, spectral_broadening, S_to_fit, n_stokes_to_fit, stokes_to_fit,
          ws, noise, metric_reference);
      
      else {
        model_to_fit->cpy_values_from(test_model);
        lm_parameter /= lm_multiplicator;
        if (lm_parameter <= 1E-5) lm_parameter = 1E-5;
      }
      
      corrected=1;
      chi_to_track = add_to_1d_array(chi_to_track,n_chi_to_track,metric);
      if (n_chi_to_track >=3){
        fp_t change = fabs((chi_to_track[n_chi_to_track-2] - chi_to_track[n_chi_to_track-1]) / chi_to_track[n_chi_to_track-1]);
        if (change < DELTA)
          to_break = 1;
      }
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

    //clock_t point4 = clock();

    //fprintf(stderr,"Time to do LM stuff: %e \n", (double)(point4-point3)/CLOCKS_PER_SEC); 
    //fprintf(stderr,"Iteration complete. \n");
    //model_to_fit->print();     
      
  }
  //io.msg(IOL_INFO, "fitting complete. Total number of iterations is : %d \n", iter-1);
  
  //clock_t point5 = clock();
  // Clean-up:
  del_ft2dim(S_to_fit,1,4,1,nlambda);
  delete[](lambda+1);
  delete[](noise+1);
  delete[](noise_scaling+1);
  if (chi_to_track)
    delete[]chi_to_track;

  // Full version:
  lambda = spectrum_to_fit->get_lambda();
  nlambda = spectrum_to_fit->get_n_lambda();
  build_from_nodes(model_to_fit);
  observable *obs_to_return = forward_evaluate(theta,phi,lambda,nlambda,scattered_light,qs_level,spectral_broadening);
   
  model_to_fit->polish_angles();   
  delete[](lambda+1);
  delete[]ws;

  //clock_t point6 = clock();

  //fprintf(stderr,"Clean -up time: %e \n", (double)(point6-point5)/CLOCKS_PER_SEC);      
  
  //fprintf(stderr,"Done! \n");
  //model_to_fit->print();     

  return obs_to_return;
} 

*/
