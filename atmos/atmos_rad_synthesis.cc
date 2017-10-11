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

// ================================================================================================

int atmosphere::op_em_vector(fp_t *** Vlos, fp_t **** B, fp_t theta,fp_t phi,fp_t* lambda,int nlambda,fp_t****** op_vector,fp_t***** em_vector){
  // This function computes opacity and emissivity for atmosphere, all wavelengths at once. 
  // Returns 0 for success, 1 for error. 
  // Arguments:
  // fp_t *** Vlos - l.o.s. velocity
  // fp_t **** B - projected magnetic field in the reference frame of the ray
  // fp_t theta, phi - emergent angle for radiation
  // fp_t * lambda, int nlambda - wavelength array and number of wavelengths 
  // fp_t ****** op_vector, NLxNXxNYxNZx4x4 array of opacities 
  // fp_t ***** em_vector, NLxNXxNYxNZx4 array of emissivities

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

  for (int a=0;a<natm;++a)
  	atml[a]->op_em_vector(T,Ne,Vlos,Vt,B,theta,phi,lambda,nlambda,op_vector,em_vector);
  
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

// ================================================================================================

int atmosphere::op_em_vector_plus_pert(fp_t *** Vlos, fp_t **** B, fp_t theta,fp_t phi,fp_t* lambda,int nlambda,
  fp_t****** op_vector,fp_t***** em_vector,fp_t******** op_vector_pert,fp_t******* em_vector_pert){
  // This function computes opacity and emissivity as well as the responses(perturbations) for atmosphere, all wavelengths at once. 
  // Returns 0 for success, 1 for error. 
  // Arguments:
  // fp_t *** Vlos - l.o.s. velocity
  // fp_t **** B - projected magnetic field in the reference frame of the ray
  // fp_t theta, phi - emergent angle for radiation
  // fp_t * lambda, int nlambda - wavelength array and number of wavelengths 
  // fp_t ****** op_vector, NLxNXxNYxNZx4x4 array of opacities 
  // fp_t ***** em_vector, NLxNXxNYxNZx4 array of emissivities
  // fp_t ******** op_vector_pert, NLxNPxNZxNXxNYxNZx4x4 array of opacity responses
  // fp_t ******* em_vector, NLxNPxNZxNXxNYxNZx4 array of emissivity responses

  // First we account for electron opacity and emissivity. Thomson scattering is very weakly dependent on wavelength
  // so we will just compute it for the middle wavelenght here:
  fp_t lambda_m = (lambda[nlambda] + lambda[1]) * 0.5;

  //memset input arrays to zero:
  memset(op_vector[1][x1l][x2l][x3l][1]+1,0,nlambda*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));
  memset(em_vector[1][x1l][x2l][x3l]+1,0,nlambda*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  memset(op_vector_pert[1][1][x3l][x1l][x2l][x3l][1]+1,0,nlambda*7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));
  memset(em_vector_pert[1][1][x3l][x1l][x2l][x3l]+1,0,nlambda*7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));


  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        fp_t op = Ne[x1i][x2i][x3i] * 6.65E-25;
        fp_t *op_pert = new fp_t[2]-1; // Only for temperature and density
        op_pert[1] = Ne_lte_der[1][x1i][x2i][x3i] * 6.65E-25;
        op_pert[2] = Ne_lte_der[2][x1i][x2i][x3i] * 6.65E-25;
        
        fp_t em = op * Planck_f(lambda_m, T[x1i][x2i][x3i]);
        fp_t *em_pert = new fp_t[2]-1; // Only for temperature and density
        em_pert[1] = em/op * op_pert[1] + Ne[x1i][x2i][x3i] * 6.65E-25 * Planck_f_derivative(lambda_m,T[x1i][x2i][x3i]);
        em_pert[2] = em/op * op_pert[2];
        
        for (int l=1;l<=nlambda;++l){
          op_vector[l][x1i][x2i][x3i][1][1] += op;
          em_vector[l][x1i][x2i][x3i][1] = +em;
          for (int p=1;p<=2;++p){
            op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][1] += op_pert[p];
            em_vector_pert[l][p][x3i][x1i][x2i][x3i][1] += em_pert[p];
          }
        }
        delete[](op_pert+1); delete[](em_pert+1);
  }

  // Then add all the contributors from opacity and emissivity from atoms and molecules
  
  for (int a=0;a<natm;++a)
    atml[a]->op_em_vector_plus_pert(T,Ne,Vlos,Vt,B,theta,phi,lambda,nlambda,op_vector,em_vector,op_vector_pert,em_vector_pert);
  
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
  // But also reorder the responses:
  for (int l=1;l<=nlambda;++l)
    for (int p=1;p<=7;++p)
      for (int x3k=x3l;x3k<=x3h;++x3k)
        for (int x1i=x1l;x1i<=x1h;++x1i)
          for (int x2i=x2l;x2i<=x2h;++x2i)
            for (int x3i=x3l;x3i<=x3h;++x3i){
              op_vector_pert[l][p][x3k][x1i][x2i][x3i][4][4] = op_vector_pert[l][p][x3k][x1i][x2i][x3i][3][3] = op_vector_pert[l][p][x3k][x1i][x2i][x3i][2][2] = op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][1];
              op_vector_pert[l][p][x3k][x1i][x2i][x3i][2][1] = op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][2];
              op_vector_pert[l][p][x3k][x1i][x2i][x3i][3][1] = op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][3];
              op_vector_pert[l][p][x3k][x1i][x2i][x3i][4][1] = op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][4];
              op_vector_pert[l][p][x3k][x1i][x2i][x3i][3][2] = -op_vector_pert[l][p][x3k][x1i][x2i][x3i][2][3]; // Rho_V
              op_vector_pert[l][p][x3k][x1i][x2i][x3i][2][4] = -op_vector_pert[l][p][x3k][x1i][x2i][x3i][4][2]; // Rho_U
              op_vector_pert[l][p][x3k][x1i][x2i][x3i][4][3] = -op_vector_pert[l][p][x3k][x1i][x2i][x3i][3][4]; // Rho_Q
  }
  return 0;
}

// ================================================================================================

observable *atmosphere::obs_stokes_responses(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** intensity_responses, model* input_model)
{
  // This function synthesises the spectrum, and computes the response functions of the 
  // unpolarized intensity. It can compute the derivatives of populations and number densities
  // analyticaly or numerically, but the final computation of perturbation is computed analyticaly.

  boundary_condition_for_rt = -1;
  popsetup(); // setup

  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);// all atoms now have pointer to this atmosphere

  if(tau_grid){ 
    compute_op_referent();
    compute_op_referent_derivative();
  }
  
  nltepops();
  respsetup();
  ne_lte_derivatives();
  
  //compute_nlte_population_responses_numerical(x3l,x3h);

  compute_nlte_population_responses(0);//

  io.msg(IOL_INFO,"atmosphere::obs: synthesizing observable and responses: theta=%f phi=%f\n",theta,phi);

  fp_t ***Vr=project(Vx,Vy,Vz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h);  // radial projection
  fp_t ****B=transform(Bx,By,Bz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h); // radial projection

  fp_t ****S=ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);   
  fp_t ***Lambda_approx = ft3dim(x1l, x1h, x2l, x2h, x3l, x3h);
  
  for (int a=0;a<natm;++a){
    atml[a]->rtsetup(x1l,x1h,x2l,x2h,x3l,x3h);
    atml[a]->zeeman_setup();
  }

  // Intensity perturbation:
  fp_t *****dS = ft5dim(x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t **** d_obs_a = ft4dim(1,7,x3l, x3h,1,nlambda,1,4);
  memset(d_obs_a[1][x3l][1]+1,0,7*(x3h-x3l+1)*4*nlambda*sizeof(fp_t));

  fp_t * lambda_vacuum = airtovac(lambda+1,nlambda);
  lambda_vacuum -=1;

  class observable *o=new observable(1,1,4,nlambda);

  fp_t****** op = ft6dim(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  fp_t***** em = ft5dim(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t******** op_pert = ft8dim(1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  fp_t******* em_pert = ft7dim(1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);

  op_em_vector_plus_pert(Vr,B,theta,phi,lambda_vacuum,nlambda,op,em,op_pert,em_pert);

  //for (int x3i=x3l;x3i<=x3h;++x3i)
    //fprintf(stderr,"%d %e %e \n",x3i,op_pert[1][3][x3i][x1l][x2l][x3i][1][1], em_pert[1][3][x3i][x1l][x2l][x3i][1]);

  // Normalize to referent opacity, for each wavelength:
  if (tau_grid)
    for (int l=1;l<=nlambda;++l){
      normalize_to_referent_opacity(op[l],em[l],op_pert[l],em_pert[l]);
      normalize_to_referent_opacity(op[l],em[l]);
  }

  fp_t *** atm_resp_to_parameters;
  int N_parameters = 0;

  if (input_model){
    atm_resp_to_parameters =  input_model->get_response_to_parameters();
    N_parameters = input_model->get_N_nodes_total();
  }
  
  for(int l=1;l<=nlambda;++l){

    // Now we formally solve for each wavelength:
    formal(rt_grid, S,0,op[l],em[l],theta,phi,boundary_condition_for_rt);
  
    // Sorting out the magnitude of perturbations so we can solve this properly:
    if (input_model==0){
      for (int x3k=x3l;x3k<=x3h;++x3k)
        for (int x1i=x1l;x1i<=x1h;++x1i)
          for (int x2i=x2l;x2i<=x2h;++x2i)
            for (int x3i=x3l;x3i<=x3h;++x3i)
              for (int s=1;s<=4;++s){
                for (int sp=1;sp<=4;++sp){
                  op_pert[l][2][x3k][x1i][x2i][x3i][s][sp] *= Nt[x1i][x2i][x3k]*delta_Nt_frac;
                  op_pert[l][6][x3k][x1i][x2i][x3i][s][sp] *= delta_angle;
                  op_pert[l][7][x3k][x1i][x2i][x3i][s][sp] *= delta_angle;
                }
                em_pert[l][2][x3k][x1i][x2i][x3i][s] *= Nt[x1i][x2i][x3k]*delta_Nt_frac;
                em_pert[l][6][x3k][x1i][x2i][x3i][s] *= delta_angle;
                em_pert[l][7][x3k][x1i][x2i][x3i][s] *= delta_angle;
      }
      for (int param=1;param<=7;++param){
        formal_pert_numerical(dS, op[l], em[l], op_pert[l][param], em_pert[l][param], theta, phi, boundary_condition_for_rt);
  
        if (param == 2)
          for (int x3i=x3l;x3i<=x3h;++x3i)
            for (int s=1;s<=4;++s)
              dS[x3i][x1l][x2l][x3l][s] /= Nt[x1l][x2l][x3i]*delta_Nt_frac;

        if (param == 6 || param == 7)
          for (int x3i=x3l;x3i<=x3h;++x3i)
            for (int s=1;s<=4;++s)
              dS[x3i][x1l][x2l][x3l][s] /= delta_angle;
      
        for (int x3i=x3l;x3i<=x3h;++x3i)
          for (int s=1;s<=4;++s)
            d_obs_a[param][x3i][l][s] = dS[x3i][x1l][x2l][x3l][s];
      }
    }
    else { // Else input model is provided and we need to extract  needed things from it.

      fp_t ****** op_pert_params = ft6dim(1,N_parameters,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
      fp_t ***** em_pert_params = ft5dim(1,N_parameters,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
      memset(op_pert_params[1][x1l][x2l][x3l][1]+1,0,N_parameters*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));
      memset(em_pert_params[1][x1l][x2l][x3l]+1,0,N_parameters*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));

      fp_t normalizer = 1.0;
      for (int p=1;p<=N_parameters;++p)
        for (int q=1;q<=7;++q){
          if (q==6 || q==7)
            normalizer = delta_angle;
          else 
            normalizer = 1.0;
          for (int x3k=x3l;x3k<=x3h;++x3k)
            for (int x1i=x1l;x1i<=x1h;++x1i)
              for (int x2i=x2l;x2i<=x2h;++x2i)
                for (int x3i=x3l;x3i<=x3h;++x3i)
                  for (int s=1;s<=4;++s){
                    em_pert_params[p][x1i][x2i][x3i][s] += em_pert[l][q][x3k][x1i][x2i][x3i][s] 
                      * atm_resp_to_parameters[p][q][x3k] * normalizer;
                    for (int sp=1;sp<=4;++sp)
                      op_pert_params[p][x1i][x2i][x3i][s][sp] += op_pert[l][q][x3k][x1i][x2i][x3i][s][sp]
                      * atm_resp_to_parameters[p][q][x3k] * normalizer;
        }
      }
      // But now when we have these "compressed responses" we need to propagate them further.
      formal_pert_numerical(dS, op[l], em[l], op_pert_params, em_pert_params, theta, phi, 
        boundary_condition_for_rt, N_parameters);

      del_ft6dim(op_pert_params,1,N_parameters,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
      del_ft5dim(em_pert_params,1,N_parameters,x1l,x1h,x2l,x2h,x3l,x3h,1,4);

      if (intensity_responses)
      for (int p=1;p<=N_parameters;++p){
        if (input_model->which_parameter(p)==5 || input_model->which_parameter(p)==6)
          normalizer = delta_angle;
        else 
          normalizer = 1.0; 
        for (int s=1;s<=4;++s)
          intensity_responses[p][l][s] = dS[p][x1l][x2l][x3l][s]/normalizer;
      }
    }

    // Add it to the observable
    o->set(S[x1l][x2l][x3l],lambda[l],1,1,l);    
  }
    
  // This should transform responses
  transform_responses(d_obs_a, theta, phi, 1, nlambda);

  // Write down the intensity perturbations:
  if (!intensity_responses && !input_model){
    FILE * out;
    out = fopen("stokes_intensity_responses_analytical.txt", "w");
      for (int param=1;param<=7;++param){
      for (int x3i=x3l;x3i<=x3h;++x3i){
        for (int l=1;l<=nlambda;++l){
          fp_t loc = 0;
          if (tau_grid) loc = log10(-tau_referent[x1l][x2l][x3i]); else loc = x3[x3i];
          fprintf(out,"%10.10e %10.10e", loc, lambda[l]);
          for (int s=1;s<=4;++s)
            fprintf(out," %10.10e", d_obs_a[param][x3i][l][s]);
          fprintf(out," \n");
        }
      }
    }
    fclose(out);
  }
  
  del_ft4dim(S,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  del_ft3dim(Lambda_approx,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft5dim(dS,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4); 
  del_ft4dim(d_obs_a,1,7,x3l,x3h,1,nlambda,1,4);
  del_ft6dim(op,1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  del_ft5dim(em,1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  del_ft8dim(op_pert,1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  del_ft7dim(em_pert,1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  delete [](lambda_vacuum+1);
  for(int a=0;a<natm;++a){
    atml[a]->zeeman_clear();
    atml[a]->rtclean(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h);
  }
  if (input_model)
    del_ft3dim(atm_resp_to_parameters,1,N_parameters,1,7,1,x3h-x3l+1);

  respclean();
  popclean(); // all done
  
  if (tau_grid) del_ft5dim(op_referent_derivative,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  
  io.msg(IOL_INFO,"atmosphere::obs_stokes_responses: observable and responses synthesized...\n");

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  return o;

}

