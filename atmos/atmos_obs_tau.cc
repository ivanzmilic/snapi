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

// Added on 24/06/2016
// This is the file which contains duplicates of basically all the functions give in atmos_obs.cc 
// except that they are all formulated in tau as an independent variable (except h)

observable *atmosphere::obs_scalar_tau(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda)
  // You have to write this down or you will forget. 
  // This function ONLY does the synthesis of the spectra. 
{

  
  boundary_condition_for_rt = -1;
  popsetup(); // setup
  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);// all atoms now have pointer to this atmosphere

  nltepops_taugrid();
  atml[0]->print_populations();

  io.msg(IOL_INFO,"atmosphere::obs: synthesizing observable: theta=%f phi=%f\n",theta,phi);

  fp_t ***Vr=project(Vx,Vy,Vz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h);  // radial projection
  fp_t ****B=transform(Bx,By,Bz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h); // radial projection
  
  fp_t ***S=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);   
  fp_t ***Lambda_approx = ft3dim(x1l, x1h, x2l, x2h, x3l, x3h);
  
  for (int a = 0; a<natm; ++a)
    atml[a]->prof_setup();

  class observable *o=new observable(1);

  for(int l=0;l<nlambda;++l){

    for (int a = 0; a<natm; ++a)
      atml[a]->prof_init();

    fp_t ***op=opacity(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);
    fp_t ***em=emissivity(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);

    // Normalize the opacity/emissivity
    fp_t *** op_rel = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    fp_t *** em_rel = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    
    for (int x3i=x3l;x3i<=x3h;++x3i){
      op_rel[x1l][x2l][x3i] = op[x1l][x2l][x3i] / op_referent[x1l][x2l][x3i];
      em_rel[x1l][x2l][x3i] = em[x1l][x2l][x3i] / op_referent[x1l][x2l][x3i];
    }

    formal(tau_referent[x1l][x2l], S,Lambda_approx,op_rel,em_rel,theta,phi, boundary_condition_for_rt);      // [Stokes] solution and approximate operator
    
    fp_t lambda_air = vactoair(lambda[l]);

    //printf("%1.15e %1.15e \n", lambda[l], S[x1l][x2l][x3l]);

    // Add it to the observable
    o->add(&(S[x1l][x2l][x3l])-1,lambda_air);
    // Delete all wavelength dependent quantities:            
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(op_rel,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(em_rel,x1l,x1h,x2l,x2h,x3l,x3h);
    
  }
  
  // Clear un-needed quantities
  for (int a = 0; a<natm; ++a){
    atml[a]->prof_clear();
  }
  
  del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(Lambda_approx,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
  popclean(); // all done

  io.msg(IOL_INFO,"atmosphere::obs: observable synthesized...\n");

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  return o;

}

observable *atmosphere::obs_scalar_num_responses_tau(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** intensity_responses){
  // This function computes the responses in completely numerical way. This means that we:
  // Perturb quantity in one direction (we use symmetric numerical derivatives), compute spectrum
  // Perturb quantity in other direction, compute spectrum
  // Subtract and divide by step. 
 
  // Compute observable:
  observable *o = obs_scalar_tau(theta, phi, lambda, nlambda);
  fp_t ** S = o->get_S(1,1);

  fp_t d_T = delta_T;
  int ns = 1;

  fp_t **** d_obs = ft4dim(1,7,x3l, x3h, 1, ns, 1, nlambda);
  memset(d_obs[1][x3l][1]+1,0,(x3h-x3l+1)*(ns)*(nlambda)*sizeof(fp_t));

  if (intensity_responses) // If provided, copy the intensity to the input array
    memset(intensity_responses[1][x3l]+1,0,(7*nlambda*(x3h-x3l+1)));

  for (int x3i = x3l;x3i<=x3h;++x3i){

    // Temperature finite-difference perturbations

    T[x1l][x2l][x3i] += d_T * 0.5;

    observable * temp_o = obs_scalar_tau(theta, phi, lambda, nlambda);

    fp_t ** Stemp = temp_o->get_S(1,1);

    memcpy(d_obs[1][x3i][1]+1, Stemp[1]+1, nlambda*sizeof(fp_t));

    delete temp_o;

    T[x1l][x2l][x3i] -= d_T;

    temp_o = obs_scalar_tau(theta, phi, lambda, nlambda);

    // Now we need to subtract, this is unfortunately more complicated:

    Stemp = temp_o->get_S(1,1);

    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[1][x3i][s][l] -= Stemp[s][l];
        d_obs[1][x3i][s][l] /= d_T;
      }

    delete temp_o;

    T[x1l][x2l][x3i] += 0.5 * d_T;

    // Density finite difference perturbations:

    fp_t d_Nt = 1E-2 * Nt[x1l][x2l][x3i];
    Nt[x1l][x2l][x3i] += 0.5 * d_Nt;

    temp_o = obs_scalar_tau(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);

    memcpy(d_obs[2][x3i][1]+1, Stemp[1]+1, nlambda*sizeof(fp_t));

    delete temp_o;

    Nt[x1l][x2l][x3i] -= d_Nt;

    temp_o = obs_scalar_tau(theta, phi, lambda, nlambda);

    // Now we need to subtract, this is unfortunately more complicated:

    Stemp = temp_o->get_S(1,1);

    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[2][x3i][s][l] -= Stemp[s][l];
        d_obs[2][x3i][s][l] /= d_Nt;
      }

    delete temp_o;
    Nt[x1l][x2l][x3i] += 0.5 * d_Nt;

    // Vt finite difference perturbations:

    Vt[x1l][x2l][x3i] += 0.5 * delta_vt;

    temp_o = obs_scalar_tau(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);

    memcpy(d_obs[3][x3i][1]+1, Stemp[1]+1, nlambda*sizeof(fp_t));

    delete temp_o;

    Vt[x1l][x2l][x3i] -= delta_vt;

    temp_o = obs_scalar_tau(theta, phi, lambda, nlambda);

    // Now we need to subtract, this is unfortunately more complicated:

    Stemp = temp_o->get_S(1,1);

    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[3][x3i][s][l] -= Stemp[s][l];
        d_obs[3][x3i][s][l] /= delta_vt;
      }

    delete temp_o;
    Vt[x1l][x2l][x3i] += 0.5 * delta_vt;

    io.msg(IOL_INFO, "atmosphere::obs_scalar_num_responses: computed perturbation with respect to x3i = %d \n", x3i);
      
  }

  // Please write it down:

  fp_t * lambda_air = new fp_t [nlambda];
  for (int l=0;l<nlambda;++l)
    lambda_air[l] = vactoair(lambda[l]);

  if (intensity_responses)
    // For the moment, just the temperature:
    for (int param=1;param<=7;++param)
      for (int l=1;l<=nlambda;++l)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          intensity_responses[param][x3i][l] = d_obs[param][x3i][1][l];


  FILE * out;
  out = fopen("intensity_responses_fin_diff.txt", "w");
  for (int param=1;param<=7;++param){
    fp_t norm = 1.0;
    for (int x3i=x3l;x3i<=x3h;++x3i){
      //norm = (param == 2) ? Nt[x1l][x2l][x3i] : 1.0;  
      for (int l=1;l<=nlambda;++l){
        fprintf(out,"%10.10e %10.10e", x3[x3i], lambda_air[l-1]);
        for (int s=1;s<=1;++s)
          fprintf(out," %10.10e", d_obs[param][x3i][1][l]*norm);
        fprintf(out," \n");
      }
    }
  }
  fclose(out);

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  io.msg(IOL_INFO, "atmosphere::obs_scalar_num_responses: computed responses of the intensity numericaly...\n");
  
  delete []lambda_air;
  del_ft4dim(d_obs,1,7,x3l,x3h,1,ns,1,nlambda);

  return o;

}

observable *atmosphere::obs_scalar_responses_tau(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** intensity_responses)
  // This function synthesises the spectrum, and computes the response functions of the 
  // unpolarized intensity. It can compute the derivatives of populations and number densities
  // analyticaly or numerically, but the final computation of perturbation is computed analyticaly.
  // Everything is done in tau grid.
{

  boundary_condition_for_rt = -1;
  popsetup(); // setup

  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);// all atoms now have pointer to this atmosphere

  compute_op_referent_derivative();
 
  nltepops_taugrid();
  respsetup();

  ne_lte_derivatives();
  //compute_nlte_population_responses_numerical_taugrid(x3l, x3h);
  compute_nlte_population_responses_taugrid(0);//


  io.msg(IOL_INFO,"atmosphere::obs: synthesizing observable and responses: theta=%f phi=%f\n",theta,phi);

  fp_t ***Vr=project(Vx,Vy,Vz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h);  // radial projection
  fp_t ****B=transform(Bx,By,Bz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h); // radial projection
  
  fp_t ***S=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);   
  fp_t ***Lambda_approx = ft3dim(x1l, x1h, x2l, x2h, x3l, x3h);
  
  for (int a = 0; a<natm; ++a)
    atml[a]->prof_setup();

  // Intensity perturbation:
  fp_t ****dS = ft4dim(x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t **** d_obs_a = ft4dim(1,7,x3l, x3h, 1, 1, 1, nlambda);

  fp_t * lambda_air = new fp_t [nlambda];

  class observable *o=new observable(1);

  if (intensity_responses) // If provided, copy the intensity to the input array
    memset(intensity_responses[1][x3l]+1,0,(7*nlambda*(x3h-x3l+1)));
  
  for(int l=0;l<nlambda;++l){

    for (int a = 0; a<natm; ++a)
      atml[a]->prof_init();

    fp_t ***op=opacity(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);
    fp_t ***em=emissivity(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);
    fp_t ***** op_pert = opacity_pert(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);
    fp_t ***** em_pert = emissivity_pert(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);

    
    for (int p=1;p<=7;++p)
   		for (int x3k=x3l;x3k<=x3h;++x3k)
   			for (int x1i=x1l;x1i<=x1h;++x1i)
   				for (int x2i=x2l;x2i<=x2h;++x2i)
   					for (int x3i=x3l;x3i<=x3h;++x3i){
   						 op_pert[p][x3k][x1i][x2i][x3i] = op_pert[p][x3k][x1i][x2i][x3i]/op_referent[x1i][x2i][x3i]
   							  -op[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i] * op_referent_derivative[p][x3k][x1i][x2i][x3i];
   						 em_pert[p][x3k][x1i][x2i][x3i] = em_pert[p][x3k][x1i][x2i][x3i]/op_referent[x1i][x2i][x3i] 
   							  -em[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i] * op_referent_derivative[p][x3k][x1i][x2i][x3i];		
              
   	}
   	for (int x1i=x1l;x1i<=x1h;++x1i)
   		for (int x2i=x2l;x2i<=x2h;++x2i)
   			for (int x3i=x3l;x3i<=x3h;++x3i){
   					op[x1i][x2i][x3i] /= op_referent[x1i][x2i][x3i];
   					em[x1i][x2i][x3i] /= op_referent[x1i][x2i][x3i];
   	}

    formal(tau_referent[x1l][x2l], S,Lambda_approx,op,em,theta,phi, boundary_condition_for_rt);      // [Stokes] solution and approximate operator
    for (int param=1;param<=7;++param){
      formal_pert_analytical_taugrid(dS, op, em, op_pert[param], em_pert[param], theta, phi, boundary_condition_for_rt);
      //formal_pert_numerical_taugrid(dS, op, em, op_pert[param], em_pert[param], op_referent, op_referent_derivative[param], theta, phi, boundary_condition_for_rt);

      for (int x3i=x3l;x3i<=x3h;++x3i){
        d_obs_a[param][x3i][1][l+1] = dS[x3i][x1l][x2l][x3l];
      }

      if (intensity_responses)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          intensity_responses[param][x3i][l+1] = dS[x3i][x1l][x2l][x3l];
    }

    // Go to lambda_air
    lambda_air[l] = vactoair(lambda[l]);

    //printf("%1.15e %1.15e \n", lambda[l], S[x1l][x2l][x3l]);

    // Add it to the observable
    o->add(&(S[x1l][x2l][x3l])-1,lambda_air[l]);

    
    // Delete all wavelength dependent quantities:            
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft5dim(em_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft5dim(op_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  }


  // Write down the intensity perturbations:
  FILE * out;
  out = fopen("intensity_responses_analytical.txt", "w");
  for (int param=1;param<=7;++param){
    fp_t norm = 1.0;
    for (int x3i=x3l;x3i<=x3h;++x3i){
      //norm = (param == 2) ? Nt[x1l][x2l][x3i] : 1.0;  
      for (int l=1;l<=nlambda;++l){
        fprintf(out,"%10.10e %10.10e", x3[x3i], lambda_air[l-1]);
        for (int s=1;s<=1;++s)
          fprintf(out," %10.10e", d_obs_a[param][x3i][1][l]*norm);
        fprintf(out," \n");
      }
    }
  }
  fclose(out);

  // Clear un-needed quantities
  for (int a = 0; a<natm; ++a){
    atml[a]->prof_clear();
  }

  del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(Lambda_approx,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
 
  del_ft4dim(dS,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h); 
  del_ft4dim(d_obs_a,1,7,x3l,x3h,1,1,1,nlambda);
  delete []lambda_air;

//
  popclean(); // all done
  respclean();

  delete_op_referent_derivative();

  io.msg(IOL_INFO,"atmosphere::obs_scalar_responses: observable and responses synthesized...\n");

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  return o;

}



