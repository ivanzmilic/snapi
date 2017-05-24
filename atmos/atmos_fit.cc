#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "const.h"
#include "math.h"

#include "atmos.h"
#include "profiles.h"
#include "obs.h"
#include "mathtools.h"

#include <ctime>


observable * atmosphere::stokes_lm_fit(observable * spectrum_to_fit, fp_t theta, fp_t phi, model * model_to_fit){

  
  // Perform LM fitting and return the best fitting spectrum
  // First extract the spectrum from the observation:

  fp_t ** S_to_fit = spectrum_to_fit->get_S_to_fit(1,1);
  int nlambda = spectrum_to_fit->get_n_lambda_to_fit();
  fp_t * lambda = spectrum_to_fit->get_lambda_to_fit();
  //fp_t ** S_to_fit = spectrum_to_fit->get_S(1,1);
  //int nlambda = spectrum_to_fit->get_n_lambda();
  //fp_t * lambda = spectrum_to_fit->get_lambda();

  /*FILE * fff = fopen("test_me.txt","w");
  for (int l=1;l<=nlambda;++l){
    fprintf(fff,"%e ",lambda[l]);
    for (int s=1;s<=4;++s) fprintf(fff,"%e ",S_to_fit[s][l]);
    fprintf(fff,"\n");
  }
  fclose(fff);
  exit(1);*/

  set_grid(1);
  
  // Set initial value of Levenberg-Marquardt parameter
  fp_t lm_parameter = 1E-3;
  
  // Some fitting related parameters:
  fp_t metric = 0.0;
  int iter = 0;
  int MAX_ITER = 30;
  fp_t * chi_to_track = 0;
  int n_chi_to_track = 0;
  int corrected = 0;
  fp_t ws[4]; ws[0] = 1.0; ws[1] = ws[2] = 0.0; ws[3] = 0.0; // weights for Stokes parameters
  fp_t *noise = new fp_t [nlambda]-1; // wavelength dependent noise
  for (int l=1;l<=nlambda;++l)
   noise[l] = sqrt(S_to_fit[1][l]*1E-6); 
  
  observable * current_obs;
  fp_t *** derivatives_to_parameters;
  fp_t ** S_current;
  
  int N_parameters = model_to_fit->get_N_nodes_total();
  
  io.msg(IOL_INFO, "atmosphere::stokes_lm_fit : entering iterative procedure\n");

  for (iter=1;iter<=MAX_ITER;++iter){

    int to_break = 0;
  
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
   
    fp_t * residual = new fp_t [nlambda*n_stokes_to_fit]-1;
      
    // Start by computing Chisq, and immediately the response of the current spectrum to the nodes
   
    if(corrected || iter == 1){
      derivatives_to_parameters = ft3dim(1,N_parameters,1,nlambda,1,4);
      memset(derivatives_to_parameters[1][1]+1,0,N_parameters*nlambda*4*sizeof(fp_t));
      current_obs = obs_stokes_responses_to_nodes_new(model_to_fit, theta, phi, lambda, nlambda, derivatives_to_parameters);    
      S_current = current_obs->get_S(1,1);
    }

    // This can also probably be absorbed, but for clarity lets not do it
    metric = 0.0;
    
    for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=n_stokes_to_fit;++s){
        int stf = stokes_to_fit[s-1];
        residual[(l-1)*n_stokes_to_fit+s] = S_to_fit[stf][l] - S_current[stf][l]; 
      
        metric += residual[(l-1)*n_stokes_to_fit+s] * residual[(l-1)*n_stokes_to_fit+s] 
        *ws[stf-1] / noise[l] / noise[l] / (n_stokes_to_fit*nlambda-N_parameters);
    }

    // Check if metric is maybe already close enough to interrupt:
    if (metric<1.0){
      to_break = 1;
      del_ft2dim(S_current,1,4,1,nlambda);
      del_ft3dim(derivatives_to_parameters,1,N_parameters,1,nlambda,1,4);
      delete[](residual+1);
      break;
    }

    // Correct either way
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
    observable *reference_obs = obs_stokes(theta, phi, lambda, nlambda);
    fp_t ** S_reference = reference_obs->get_S(1,1);

    // Compute new chi sq
    fp_t metric_reference = 0.0;
    for (int l=1;l<=nlambda;++l)
      for (int s=1;s<=n_stokes_to_fit;++s){
        int stf=stokes_to_fit[s-1];
        metric_reference += (S_to_fit[stf][l] - S_reference[stf][l]) * (S_to_fit[stf][l] - S_reference[stf][l])
         *ws[stf-1] / noise[l] / noise[l] / (n_stokes_to_fit*nlambda-N_parameters);
       }
    
    
    if (metric_reference < metric){
      // Everything is ok, and we can decrease lm_parameter:
      lm_parameter /= 10.0;
      corrected = 1;
      model_to_fit->cpy_values_from(test_model);
      chi_to_track = add_to_1d_array(chi_to_track,n_chi_to_track,metric_reference);
      // Check if we need to interrupt:
      if (n_chi_to_track >=3)
        if ((chi_to_track[n_chi_to_track-2] - chi_to_track[n_chi_to_track-1])/chi_to_track[n_chi_to_track-1]<1E-2
          && (chi_to_track[n_chi_to_track-3] - chi_to_track[n_chi_to_track-2])/chi_to_track[n_chi_to_track-2]<1E-2)
          to_break = 1;
    }
    else{
      lm_parameter *= 10.0;
      corrected = 0;
    }
    delete test_model;
    del_ft2dim(J_transpose,1,N_parameters,1,n_stokes_to_fit*nlambda);
    del_ft2dim(JTJ,1,N_parameters,1,N_parameters);
    del_ft2dim(J,1,nlambda*n_stokes_to_fit,1,N_parameters);
    del_ft2dim(S_reference,1,4,1,nlambda);
    if (corrected){
      delete current_obs;
      del_ft2dim(S_current,1,4,1,nlambda);
      del_ft3dim(derivatives_to_parameters,1,N_parameters,1,nlambda,1,4);
    }
    delete reference_obs;
    delete [](residual+1);
    delete [](rhs+1);
    delete [](correction+1);
    metric = 0.0;
    if(to_break)
      break;
  }

  io.msg(IOL_INFO, "fitting complete. Total number of iterations is : %d \n", iter-1);
  
  // Clean-up:
  del_ft2dim(S_to_fit,1,4,1,nlambda);
  delete[](lambda+1);
  delete[](noise+1);
  if (chi_to_track)
    delete[]chi_to_track;
  
  lambda = spectrum_to_fit->get_lambda();
  nlambda = spectrum_to_fit->get_n_lambda();
  observable *obs_to_return = obs_stokes(theta, phi, lambda, nlambda);
  delete[](lambda+1);
  return obs_to_return;
}

// ======================================================================================================================================================


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
    interpolate_from_nodes(atmos_model); // <------- Can be done better but cannot be bothered now. FIX! 

    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t B_mag = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i]+By[x1l][x2l][x3i]*By[x1l][x2l][x3i]+Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
      fp_t theta_B = acos(Bz[x1l][x2l][x3i]/B_mag);
      fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);
      //if (i>N_nodes_T+N_nodes_vt+N_nodes_vs)
      //  printf("%d %e %e %e",B_mag, theta_B, phi_B);
      local_perturbations[x3i] = Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B;
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    interpolate_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t B_mag = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i]+By[x1l][x2l][x3i]*By[x1l][x2l][x3i]+Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
      fp_t theta_B = acos(Bz[x1l][x2l][x3i]/B_mag);
      fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);
      //printf(B_mag, theta_B, phi_B);
      local_perturbations[x3i] -= (Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B);
      local_perturbations[x3i] /= step;
      //if (i>N_nodes_T+N_nodes_vt+N_nodes_vs)
      //  printf("%d %d %e \n", i, x3i, local_perturbations[x3i]);
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
      //printf("%d %e %e \n", x3i, local_T_perturbations[x3i], local_Nt_perturbations[x3i]/Nt[x1l][x2l][x3i]);
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
      //if (i>N_nodes_T+N_nodes_vt+N_nodes_vs)
      //  printf("%d %e %e %e",B_mag, theta_B, phi_B);
      resp_atm_to_parameters[i][param][x3i-x3l+1] = Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B;
    }
    
    atmos_model->perturb_node_value(i, -1.0 * step);
    interpolate_from_nodes(atmos_model);

    for (int x3i=x3l;x3i<=x3h;++x3i){
      fp_t B_mag = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i]+By[x1l][x2l][x3i]*By[x1l][x2l][x3i]+Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
      fp_t theta_B = acos(Bz[x1l][x2l][x3i]/B_mag);
      fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);
      resp_atm_to_parameters[i][param][x3i-x3l+1] -= (Vt[x1l][x2l][x3i] + Vz[x1l][x2l][x3i] + B_mag + theta_B + phi_B);
      resp_atm_to_parameters[i][param][x3i-x3l+1] /= step;
    }

    atmos_model->perturb_node_value(i, 0.5 * step);
    interpolate_from_nodes(atmos_model);
    delete[](local_perturbations+x3l);
  
  }

  // Now we start the responses part.

  atmos_model->set_response_to_parameters(resp_atm_to_parameters,x3h-x3l+1);

  observable *o = obs_stokes_responses(theta, phi, lambda, nlambda, response_to_parameters,atmos_model);
  
  del_ft3dim(resp_atm_to_parameters,1,N_parameters,1,7,1,N_depths);
  
  return o;
  
 }

 /////////////////////
 // OLD:::::

 observable * atmosphere::scalar_lm_fit(observable * spectrum_to_fit, fp_t theta, fp_t phi, fp_t * lambda, int nlambda){

  // Perform LM fitting and return the best fitting spectrum
  // First extract the spectrum from the observation:
  fp_t ** stokes_vector_to_fit = spectrum_to_fit->get_S(1,1);
  
  // This is a starting model
  int N_temp_nodes = 4;
  int N_vt_nodes = 1;
  int N_parameters = N_temp_nodes + N_vt_nodes;
  model * current_model = model_new(N_temp_nodes,N_vt_nodes,0,0);
    
  // Temperature nodes:
  fp_t * temp_nodes_tau = new fp_t [N_temp_nodes] -1;
  fp_t * temp_nodes_temp = new fp_t [N_temp_nodes] -1;
  temp_nodes_tau[1] = -4.0;
  temp_nodes_tau[2] = -2.5;
  temp_nodes_tau[3] = -1.0;
  temp_nodes_tau[4] = 0.5;
    
  temp_nodes_temp[1] = 6000.0;
  temp_nodes_temp[2] = 7000.0;
  temp_nodes_temp[3] = 8000.0;
  temp_nodes_temp[4] = 8500.0;
    
  current_model->set_temp_nodes(temp_nodes_tau, temp_nodes_temp);

  // Vt nodes:
  fp_t * vt_nodes_tau = new fp_t [N_vt_nodes] -1;
  fp_t * vt_nodes_vt = new fp_t [N_vt_nodes] -1;
  vt_nodes_tau[1] = 0;
  vt_nodes_vt[1] = 1.5E5;
  current_model->set_vt_nodes(vt_nodes_tau,vt_nodes_vt);

  fp_t metric = 0.0;
  fp_t lm_parameter = 1E-3;
  FILE * output;
  output = fopen("fitting_log.txt", "w");

  int iter = 0;
  int MAX_ITER = 1;

  io.msg(IOL_INFO, "atmosphere::scalar_lm_fit : entering iterative procedure\n");

  FILE * detailed_log;
  detailed_log = fopen("detailed_log.txt", "w");
  
  for (iter = 1; iter <=MAX_ITER; ++iter){
  
    fp_t ** derivatives_to_parameters_num = ft2dim(1,N_parameters,1,nlambda);
    memset(derivatives_to_parameters_num[1]+1,0,N_parameters*nlambda*sizeof(fp_t));
    fp_t ** derivatives_to_parameters = ft2dim(1,N_parameters,1,nlambda);
    memset(derivatives_to_parameters[1]+1,0,N_parameters*nlambda*sizeof(fp_t));
    fp_t * residual = new fp_t [nlambda] - 1;
      
    // Start by computing Chisq, and immediately the response of the current spectrum to the nodes

    observable *current_obs = obs_scalar_responses_to_nodes_tau(current_model, theta, phi, lambda, nlambda, derivatives_to_parameters);
    /*obs_scalar_num_responses_to_nodes_tau(current_model, theta, phi, lambda, nlambda, derivatives_to_parameters_num);
          
    for (int i=1;i<=N_parameters;++i){
      fprintf(detailed_log,"NODE# %d \n", i);
      for (int l=1;l<=nlambda;++l)
        fprintf(detailed_log,"%d %e %e %e %e \n", l, lambda[l-1]*1E8, derivatives_to_parameters_num[i][l], derivatives_to_parameters[i][l],
          (derivatives_to_parameters[i][l] - derivatives_to_parameters_num[i][l]) / derivatives_to_parameters_num[i][l]);
    }*/

    fp_t ** S = current_obs->get_S(1,1);

    for (int l=1;l<=nlambda;++l){
      residual[l] = stokes_vector_to_fit[1][l] - S[1][l]; 
      //printf("l = %d residual = %e \n", l, residual[l]);
      metric += residual[l] * residual[l] / 1E22;
    }
    fprintf(detailed_log, "Iteration # : %d Chisq : %e \n", iter, metric);
    // Now we need to correct:
    fp_t ** J = ft2dim(1,nlambda,1,N_parameters);
    for (int i=1;i<=N_parameters;++i) for (int l=1;l<=nlambda;++l) J[l][i] = derivatives_to_parameters[i][l];
    fp_t ** J_transpose = transpose(J,nlambda,N_parameters);
    fp_t ** JTJ = multiply_with_transpose(J, nlambda, N_parameters);
    for (int i=1;i<=N_parameters;++i) JTJ[i][i] *= (lm_parameter + 1.0);
    // Now correct
    fp_t * rhs = multiply_vector(J_transpose, residual, N_parameters, nlambda);
    fp_t * correction = solve(JTJ, rhs, 1, N_parameters);

    fprintf(detailed_log, "Jacobian, residual and finally, correction: \n");
    for (int i=1;i<=N_parameters;++i){
      for (int j=1;j<=N_parameters;++j)
        fprintf(detailed_log, " %e ", JTJ[i][j]);
      fprintf(detailed_log, " %e %e \n", residual[i], correction[i]);
    }

    for (int i=1;i<=N_parameters;++i)
      current_model->perturb_node_value(i, correction[i]);
    current_model->print();

    build_from_nodes(current_model);
    observable *reference_obs = obs_scalar_tau(theta, phi, lambda, nlambda);
    fp_t ** S_reference = reference_obs->get_S(1,1);

    // Compute new chi sq
    fp_t metric_reference = 0.0;
    for (int l=1;l<=nlambda;++l){
        metric_reference += (stokes_vector_to_fit[1][l] - S_reference[1][l]) * (stokes_vector_to_fit[1][l] - S_reference[1][l]) / 1E22;
      }
      fprintf(output, "%d %e \n", iter, metric);
      if (metric_reference < metric){
        // Everything is ok, and we can decrease lm_parameter:
        lm_parameter /= 10.0;
        //fprintf(output, "Good step! Decreasing lm parameter.\n");
      }
      else{
        // We are in the more non-linear regime, so we want linearization and hence, gradient descent
        lm_parameter *= 10.0;
        // And un-modify:
        for (int i=1;i<=N_parameters;++i)
        current_model->perturb_node_value(i, -1.0 * correction[i]);
      //fprintf(output, "Bad step! Undoing modification and increasing lm parameter.\n");
      }

    if (iter == MAX_ITER){
      FILE * result;
      result = fopen("spectrum_fitted.dat", "w");
      for (int l=1;l<=nlambda;++l)
        fprintf(result, "%e %e %e \n", lambda[l-1], stokes_vector_to_fit[1][l], S_reference[1][l]);
      fclose(result);
    }

    del_ft2dim(J_transpose,1,N_parameters,1,nlambda);
    del_ft2dim(JTJ,1,N_parameters,1,N_parameters);
    del_ft2dim(J,1,nlambda,1,N_parameters);
    del_ft2dim(derivatives_to_parameters,1,N_parameters,1,nlambda);
    del_ft2dim(derivatives_to_parameters_num,1,N_parameters,1,nlambda);
    del_ft2dim(S,1,1,1,nlambda);
    del_ft2dim(S_reference,1,1,1,nlambda);
    delete current_obs;
    delete reference_obs;

    delete [](residual+1);
    delete [](rhs+1);
    delete [](correction+1);
    metric = 0.0;
    // if (to_break < converged)
     //break;
  }

  fclose(output);
  fclose(detailed_log);

  io.msg(IOL_INFO, "fitting complete. Total number of iterations is : %d \n", iter-1);

  // Clean-up:
  del_ft2dim(stokes_vector_to_fit,1,1,1,nlambda);
  delete [](temp_nodes_tau+1);
  delete [](temp_nodes_temp+1);
  delete current_model;

  return 0;
}
