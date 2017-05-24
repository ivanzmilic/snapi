#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "const.h"
#include "math.h"

#include "atmos.h"
#include "profiles.h"
#include "obs.h"
#include "mathtools.h"

observable *atmosphere::obs_scalar(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda)
  // You have to write this down or you will forget. 
  // This function ONLY does the synthesis of the spectra. 
{

  
  boundary_condition_for_rt = -1;
  popsetup(); // setup
  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);// all atoms now have pointer to this atmosphere

  nltepops();
  atml[0]->print_populations(); // Print hydrogen populations. 

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
    
    //formal(tau_referent[x1l][x2l], S,Lambda_approx,op_rel,em_rel,theta,phi, boundary_condition_for_rt);      // [Stokes] solution and approximate operator
    formal(x3, S,Lambda_approx,op,em,theta,phi, boundary_condition_for_rt);

    fp_t lambda_air = vactoair(lambda[l]);

    // Add it to the observable
    o->add(&(S[x1l][x2l][x3l])-1,lambda_air);
    // Delete all wavelength dependent quantities:            
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
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

observable *atmosphere::obs_scalar_responses(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** intensity_responses)
  // This function synthesises the spectrum, and computes the response functions of the 
  // unpolarized intensity. It can compute the derivatives of populations and number densities
  // analyticaly or numerically, but the final computation of perturbation is computed analyticaly.
{

  boundary_condition_for_rt = -1;
  popsetup(); // setup

  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);// all atoms now have pointer to this atmosphere

  if(tau_grid) compute_op_referent_derivative();
 
  nltepops();
  respsetup();

  ne_lte_derivatives();
  //compute_nlte_population_responses_numerical(x3l, x3h);
  compute_nlte_population_responses(0);//

  io.msg(IOL_INFO,"atmosphere::obs: synthesizing observable and responses: theta=%f phi=%f\n",theta,phi);

  fp_t ***Vr=project(Vx,Vy,Vz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h);  // radial projection
  fp_t ****B=transform(Bx,By,Bz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h); // radial projection
  
  fp_t ***S=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);   
  fp_t ***Lambda_approx = ft3dim(x1l, x1h, x2l, x2h, x3l, x3h);
  
  for (int a = 0; a<natm; ++a)
    atml[a]->prof_setup();

  // Intensity perturbation:
  fp_t ****dS = ft4dim(x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t **** d_obs_a = ft4dim(1,7,x3l, x3h, 1, 2, 1, nlambda);
  memset(d_obs_a[1][x3l][1]+1,0,7*(x3h-x3l+1)*2*nlambda*sizeof(fp_t));

  fp_t * lambda_air = new fp_t [nlambda];

  class observable *o=new observable(1);

  if (intensity_responses) // If provided, copy the intensity to the input array
    memset(intensity_responses[1][x3l]+1,0,(7*nlambda*(x3h-x3l+1))*sizeof(fp_t));

  for(int l=0;l<nlambda;++l){

    for (int a = 0; a<natm; ++a)
      atml[a]->prof_init();

    fp_t ***op=opacity(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);
    fp_t ***em=emissivity(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);
    fp_t ***** op_pert = opacity_pert(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);
    fp_t ***** em_pert = emissivity_pert(T,Ne,Vr,Vt,B,theta,phi,lambda[l]);

    formal(rt_grid, S,Lambda_approx,op,em,theta,phi, boundary_condition_for_rt);      // [Stokes] solution and approximate operator
    for (int param=1;param<=7;++param){
      formal_pert_analytical(dS, op, em, op_pert[param], em_pert[param], theta, phi, boundary_condition_for_rt);
      //formal_pert_numerical(dS, op, em, op_pert[param], em_pert[param], theta, phi, boundary_condition_for_rt);

      for (int x3i=x3l;x3i<=x3h;++x3i)
        d_obs_a[param][x3i][1][l+1] = dS[x3i][x1l][x2l][x3l];

      if (intensity_responses)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          intensity_responses[param][x3i][l+1] = dS[x3i][x1l][x2l][x3l];
    }

    // Go to lambda_air
    lambda_air[l] = vactoair(lambda[l]);

    // Add it to the observable
    o->add(&(S[x1l][x2l][x3l])-1,lambda_air[l]);
    
    // Delete all wavelength dependent quantities:            
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft5dim(em_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft5dim(op_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  }

  FILE * out;
  out = fopen("intensity_responses_analytical.txt", "w");
  for (int param=1;param<=7;++param){
    for (int x3i=x3l;x3i<=x3h;++x3i){
      for (int l=1;l<=nlambda;++l){
        fprintf(out,"%10.10e %10.10e", x3[x3i], lambda_air[l-1]);
        fprintf(out," %10.10e", d_obs_a[param][x3i][1][l]);
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
//
  del_ft4dim(d_obs_a,1,7,x3l,x3h,1,2,1,nlambda);
  delete []lambda_air;

//
  //for(int a=0;a<natm;++a) atml[a]->radiation_moments_clean();
  respclean();
  popclean(); // all done

  io.msg(IOL_INFO,"atmosphere::obs_scalar_responses: observable and responses synthesized...\n");

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  return o;

}

observable *atmosphere::obs_scalar_num_responses(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t *** intensity_responses){
  // This function computes the responses in completely numerical way. This means that we:
  // Perturb quantity in one direction (we use symmetric numerical derivatives), compute spectrum
  // Perturb quantity in other direction, compute spectrum
  // Subtract and divide by step. 
 
  // Compute observable:
  observable *o = obs_scalar(theta, phi, lambda, nlambda);
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

    observable * temp_o = obs_scalar(theta, phi, lambda, nlambda);

    fp_t ** Stemp = temp_o->get_S(1,1);

    memcpy(d_obs[1][x3i][1]+1, Stemp[1]+1, nlambda*sizeof(fp_t));

    delete temp_o;

    T[x1l][x2l][x3i] -= d_T;

    temp_o = obs_scalar(theta, phi, lambda, nlambda);

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
    /*
    fp_t d_Nt = delta_Nt_frac * Nt[x1l][x2l][x3i];
    Nt[x1l][x2l][x3i] += 0.5 * d_Nt;
    
    temp_o = obs_scalar(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);

    memcpy(d_obs[2][x3i][1]+1, Stemp[1]+1, nlambda*sizeof(fp_t));

    delete temp_o;
   
    Nt[x1l][x2l][x3i] -= d_Nt;

    temp_o = obs_scalar(theta, phi, lambda, nlambda);

    // Now we need to subtract, this is unfortunately more complicated:

    Stemp = temp_o->get_S(1,1);

    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[2][x3i][s][l] -= Stemp[s][l];
        d_obs[2][x3i][s][l] /= d_Nt;
      }

    delete temp_o;
    Nt[x1l][x2l][x3i] += 0.5 * d_Nt;
    /*

     // Vt finite difference perturbations:

    Vt[x1l][x2l][x3i] += 0.5 * delta_vt;

    temp_o = obs_scalar(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);

    memcpy(d_obs[3][x3i][1]+1, Stemp[1]+1, nlambda*sizeof(fp_t));

    delete temp_o;

    Vt[x1l][x2l][x3i] -= delta_vt;

    temp_o = obs_scalar(theta, phi, lambda, nlambda);

    // Now we need to subtract, this is unfortunately more complicated:

    Stemp = temp_o->get_S(1,1);

    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[3][x3i][s][l] -= Stemp[s][l];
        d_obs[3][x3i][s][l] /= delta_vt;
      }

    delete temp_o;
    Vt[x1l][x2l][x3i] += 0.5 * delta_vt;
    */
    
    
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


fp_t *atmosphere::test_stokes(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda){}

// Same version as the above. This one however returns the observable, as intended

observable *atmosphere::obs_stokes(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda){

  // The similar one as the the atmosphere::obs:


  boundary_condition_for_rt = -1;
  popsetup(); // 
  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);

  nltepops();

  fp_t ***Vr=project(Vx,Vy,Vz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h);  // radial projection
  fp_t ****B=transform(Bx,By,Bz,theta,phi,x1l,x1h,x2l,x2h,x3l,x3h); // radial projection
  fp_t ****S=ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);

  for (int a=0;a<natm;++a){
    atml[a]->rtsetup(x1l,x1h,x2l,x2h,x3l,x3h);
    atml[a]->zeeman_setup();
  }

  class observable *o=new observable(1,1,4,nlambda);

  fp_t * lambda_vacuum = airtovac(lambda+1,nlambda);
  lambda_vacuum -=1;

  fp_t ****** op_vector = ft6dim(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  fp_t *****  em_vector = ft5dim(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  op_em_vector(Vr,B,theta,phi,lambda_vacuum,nlambda,op_vector,em_vector);

  for (int l = 1; l<=nlambda; ++l){

    if (tau_grid) normalize_to_referent_opacity(op_vector[l], em_vector[l]);

    formal(rt_grid, S,0,op_vector[l],em_vector[l],theta,phi,boundary_condition_for_rt); 

    o->set(S[x1l][x2l][x3l],lambda[l],1,1,l);
  }

  del_ft4dim(S,x1l, x1h, x2l, x2h, x3l, x3h, 1, 4);
  del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft6dim(op_vector,1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  del_ft5dim(em_vector,1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  delete[](lambda_vacuum+1);
  for (int a=0;a<natm;++a){
    atml[a]->rtclean(0,0,x1l,x1h,x2l,x2h,x3l,x3h);
    atml[a]->zeeman_clear();
  }
  popclean();
  
  io.msg(IOL_INFO,"atmosphere::obs: polarized observable synthesized...\n");

  // ------------------------------------------------------------------------------------------------------------

  return o;
}

observable *atmosphere::obs_stokes_responses(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t **** intensity_responses)
  // This function synthesises the spectrum, and computes the response functions of the 
  // unpolarized intensity. It can compute the derivatives of populations and number densities
  // analyticaly or numerically, but the final computation of perturbation is computed analyticaly, 
  // and propagated numerically or analytically, to get the stokes spctrum responses. .
{

  boundary_condition_for_rt = -1;
  popsetup(); // setup

  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);// all atoms now have pointer to this atmosphere

  if(tau_grid) compute_op_referent_derivative();

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
  
  for (int a=0;a<natm;++a)
    atml[a]->zeeman_setup();

  // Intensity perturbation:
  fp_t *****dS = ft5dim(x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t **** d_obs_a = ft4dim(1,7,x3l, x3h,1,nlambda,1,4);
  memset(d_obs_a[1][x3l][1]+1,0,7*(x3h-x3l+1)*4*nlambda*sizeof(fp_t));

  fp_t * lambda_air = new fp_t [nlambda];

  class observable *o=new observable(4);

  if (intensity_responses) // If provided, copy the intensity to the input array
    memset(intensity_responses[1][x3l][1]+1,0,(7*nlambda*(x3h-x3l+1))*4*sizeof(fp_t));

  clock_t begin = clock();
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  
  fp_t****** op = ft6dim(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  fp_t***** em = ft5dim(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  fp_t******** op_pert = ft8dim(1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  fp_t******* em_pert = ft7dim(1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  op_em_vector_plus_pert(Vr,B,theta,phi,lambda-1,nlambda,op,em,op_pert,em_pert);
  
  // Normalize to referent opacity, for each wavelength:
  if (tau_grid)
    for (int l=1;l<=nlambda;++l){
      normalize_to_referent_opacity(op[l],em[l],op_pert[l],em_pert[l]);
      normalize_to_referent_opacity(op[l],em[l]);
  }

  for(int l=0;l<nlambda;++l){

    // Now we formally solve for each wavelength:      
    formal(rt_grid, S,0,op[l+1],em[l+1],theta,phi,boundary_condition_for_rt);
  
    // Sorting out the magnitude of perturbations so we can solve this properly:
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i)
            for (int s=1;s<=4;++s){
              for (int sp=1;sp<=4;++sp){
                op_pert[l+1][2][x3k][x1i][x2i][x3i][s][sp] *= Nt[x1i][x2i][x3k]*delta_Nt_frac;
                op_pert[l+1][6][x3k][x1i][x2i][x3i][s][sp] *= delta_angle;
                op_pert[l+1][7][x3k][x1i][x2i][x3i][s][sp] *= delta_angle;
              }
              em_pert[l+1][2][x3k][x1i][x2i][x3i][s] *= Nt[x1i][x2i][x3k]*delta_Nt_frac;
              em_pert[l+1][6][x3k][x1i][x2i][x3i][s] *= delta_angle;
              em_pert[l+1][7][x3k][x1i][x2i][x3i][s] *= delta_angle;
    }
    for (int param=1;param<=7;++param){
      formal_pert_numerical(dS, op[l+1], em[l+1], op_pert[l+1][param], em_pert[l+1][param], theta, phi, boundary_condition_for_rt);
  
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
          d_obs_a[param][x3i][l+1][s] = dS[x3i][x1l][x2l][x3l][s];
    }
    // Go to lambda_air
    lambda_air[l] = vactoair(lambda[l]);
    // Add it to the observable
    o->add(S[x1l][x2l][x3l],lambda_air[l]);
  }
  
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time spent on op/em + pert = %f \n", time_spent);

  // This should transform responses
  transform_responses(d_obs_a, theta, phi, 1, nlambda);

  if (intensity_responses)
    memcpy(intensity_responses[1][x3l][1]+1,d_obs_a[1][x3l][1]+1,7*(x3h-x3l+1)*nlambda*4*sizeof(fp_t));
  
  // Write down the intensity perturbations only if we are doing forward modeling. If we fit, we don't want to lose 
  // time on this.
  if (!intensity_responses){
    FILE * out;
    out = fopen("stokes_intensity_responses_analytical.txt", "w");
    for (int param=1;param<=7;++param){
      for (int x3i=x3l;x3i<=x3h;++x3i){
        fp_t loc = 0;
        if (tau_grid) loc = log10(-tau_referent[x1l][x2l][x3i]); else loc = x3[x3i];  
        for (int l=1;l<=nlambda;++l){
          fprintf(out,"%10.10e %10.10e", loc, lambda_air[l-1]);
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
//
  del_ft4dim(d_obs_a,1,7,x3l,x3h,1,nlambda,1,4);

  del_ft6dim(op,1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  del_ft5dim(em,1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  del_ft8dim(op_pert,1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  del_ft7dim(em_pert,1,nlambda,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  delete []lambda_air;
  for(int a=0;a<natm;++a){
    atml[a]->zeeman_clear();
    atml[a]->rtclean(1,nlambda,x1l,x1h,x2l,x2h,x3l,x3h);
  }
  
  respclean();
  popclean(); // all done

  io.msg(IOL_INFO,"atmosphere::obs_stokes_responses: observable and responses synthesized...\n");

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  return o;

}

observable *atmosphere::obs_stokes_num_responses(fp_t theta,fp_t phi,fp_t *lambda,int32_t nlambda, fp_t **** intensity_responses){
  // This function computes the responses in completely numerical way. This means that we:
  // Perturb quantity in one direction (we use symmetric numerical derivatives), compute spectrum
  // Perturb quantity in other direction, compute spectrum
  // Subtract and divide by step. 
 
  // Compute observable:
  observable *o = obs_stokes(theta, phi, lambda, nlambda);
  
  fp_t ** S = o->get_S(1,1);
  int ns = 4;

  fp_t d_T = delta_T;
  
  fp_t **** d_obs = ft4dim(1,7,x3l, x3h, 1, ns, 1, nlambda);
  memset(d_obs[1][x3l][1]+1,0,(x3h-x3l+1)*(ns)*(nlambda)*sizeof(fp_t));

  if (intensity_responses) // If provided, copy the intensity to the input array
    memset(intensity_responses[1][x3l][1]+1,0,(7*nlambda*4*(x3h-x3l+1))*sizeof(fp_t));

  for (int x3i = x3l;x3i<=x3h;++x3i){

    // Temperature finite-difference perturbations
    
    T[x1l][x2l][x3i] += d_T * 0.5;

    observable * temp_o = obs_stokes(theta, phi, lambda, nlambda);

    fp_t ** Stemp = temp_o->get_S(1,1);
    memcpy(d_obs[1][x3i][1]+1, Stemp[1]+1, ns*nlambda*sizeof(fp_t));
    del_ft2dim(Stemp,1,ns,1,nlambda);
    
    delete temp_o;

    T[x1l][x2l][x3i] -= d_T;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);

    // Now we need to subtract, this is unfortunately more complicated:

    Stemp = temp_o->get_S(1,1);
    
    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[1][x3i][s][l] -= Stemp[s][l];
        d_obs[1][x3i][s][l] /= d_T;
      }
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;

    T[x1l][x2l][x3i] += 0.5 * d_T;
    
    // Density finite difference perturbations:
    
    fp_t d_Nt = delta_Nt_frac * Nt[x1l][x2l][x3i];
    Nt[x1l][x2l][x3i] += 0.5 * d_Nt;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);
    memcpy(d_obs[2][x3i][1]+1, Stemp[1]+1, ns*nlambda*sizeof(fp_t));
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;
   
    Nt[x1l][x2l][x3i] -= d_Nt;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);

    // Now we need to subtract, this is unfortunately more complicated:

    Stemp = temp_o->get_S(1,1);

    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[2][x3i][s][l] -= Stemp[s][l];
        d_obs[2][x3i][s][l] /= d_Nt;
      }
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;
    Nt[x1l][x2l][x3i] += 0.5 * d_Nt;

     // Vt finite difference perturbations:
    
    Vt[x1l][x2l][x3i] += 0.5 * delta_vt;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);
    memcpy(d_obs[3][x3i][1]+1, Stemp[1]+1, ns*nlambda*sizeof(fp_t));
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;

    Vt[x1l][x2l][x3i] -= delta_vt;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);
    // Now we need to subtract, this is unfortunately more complicated:
    Stemp = temp_o->get_S(1,1);
    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[3][x3i][s][l] -= Stemp[s][l];
        d_obs[3][x3i][s][l] /= delta_vt;
      }
    delete temp_o;
    del_ft2dim(Stemp,1,ns,1,nlambda);
    Vt[x1l][x2l][x3i] += 0.5 * delta_vt;
    
    
    // Then the radial, that is, macroscopic velocity

    Vz[x1l][x2l][x3i] += 0.5 * delta_vr;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);
    memcpy(d_obs[4][x3i][1]+1, Stemp[1]+1, ns*nlambda*sizeof(fp_t));
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;

    Vz[x1l][x2l][x3i] -= delta_vr;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);
    // Now we need to subtract, this is unfortunately more complicated:
    Stemp = temp_o->get_S(1,1);
    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[4][x3i][s][l] -= Stemp[s][l];
        d_obs[4][x3i][s][l] /= delta_vr;
      }
    delete temp_o;
    del_ft2dim(Stemp,1,ns,1,nlambda);
    Vz[x1l][x2l][x3i] += 0.5 * delta_vr;
    /*
    
    // Then Magnetic field, only z for the start
    fp_t B_total = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i] + By[x1l][x2l][x3i]*By[x1l][x2l][x3i] + Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
    fp_t scale_x = Bx[x1l][x2l][x3i]/B_total;
    fp_t scale_y = By[x1l][x2l][x3i]/B_total;
    fp_t scale_z = Bz[x1l][x2l][x3i]/B_total;
    if (B_total == 0.0) scale_x=scale_y=scale_z=1.0/sqrt(3.0); 
    
    Bx[x1l][x2l][x3i] += 0.5 * delta_B * scale_x;
    By[x1l][x2l][x3i] += 0.5 * delta_B * scale_y;
    Bz[x1l][x2l][x3i] += 0.5 * delta_B * scale_z;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);

    Stemp = temp_o->get_S(1,1);
    memcpy(d_obs[5][x3i][1]+1, Stemp[1]+1, ns*nlambda*sizeof(fp_t));
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;

    Bx[x1l][x2l][x3i] -= delta_B * scale_x;
    By[x1l][x2l][x3i] -= delta_B * scale_y;
    Bz[x1l][x2l][x3i] -= delta_B * scale_z;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);
    // Now we need to subtract, this is unfortunately more complicated:
    Stemp = temp_o->get_S(1,1);
    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[5][x3i][s][l] -= Stemp[s][l];
        d_obs[5][x3i][s][l] /= delta_B;
      }
    delete temp_o;
    del_ft2dim(Stemp,1,ns,1,nlambda);
    Bx[x1l][x2l][x3i] += 0.5 * delta_B * scale_x;
    By[x1l][x2l][x3i] += 0.5 * delta_B * scale_y;
    Bz[x1l][x2l][x3i] += 0.5 * delta_B * scale_z;

    // And then theta and phi
    fp_t theta_B = acos(Bz[x1l][x2l][x3i] / B_total);
    fp_t phi_B = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3i]);

    // First theta:

    fp_t d_Bx = cos(theta_B)*cos(phi_B)*B_total * delta_angle;
    fp_t d_By = cos(theta_B)*sin(phi_B)*B_total * delta_angle;
    fp_t d_Bz = -sin(theta_B)*B_total * delta_angle;
    Bx[x1l][x2l][x3i] += 0.5*d_Bx;
    By[x1l][x2l][x3i] += 0.5*d_By;
    Bz[x1l][x2l][x3i] += 0.5*d_Bz;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);
    Stemp = temp_o->get_S(1,1);
    memcpy(d_obs[6][x3i][1]+1, Stemp[1]+1, ns*nlambda*sizeof(fp_t));
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;

    Bx[x1l][x2l][x3i] -= d_Bx;
    By[x1l][x2l][x3i] -= d_By;
    Bz[x1l][x2l][x3i] -= d_Bz;
    temp_o = obs_stokes(theta, phi, lambda, nlambda);
    // Now we need to subtract, this is unfortunately more complicated:
    Stemp = temp_o->get_S(1,1);
    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[6][x3i][s][l] -= Stemp[s][l];
        d_obs[6][x3i][s][l] /= delta_angle;
      }
    delete temp_o;
    del_ft2dim(Stemp,1,ns,1,nlambda);
    Bx[x1l][x2l][x3i] += 0.5*d_Bx;
    By[x1l][x2l][x3i] += 0.5*d_By;
    Bz[x1l][x2l][x3i] += 0.5*d_Bz;

    // And then the same for phi:
    d_Bx = -sin(theta_B)*sin(phi_B)*B_total * delta_angle;
    d_By = sin(theta_B)*cos(phi_B)*B_total * delta_angle;
    d_Bz = 0.0;
    Bx[x1l][x2l][x3i] += 0.5*d_Bx;
    By[x1l][x2l][x3i] += 0.5*d_By;
    Bz[x1l][x2l][x3i] += 0.5*d_Bz;

    temp_o = obs_stokes(theta, phi, lambda, nlambda);
    Stemp = temp_o->get_S(1,1);
    memcpy(d_obs[7][x3i][1]+1, Stemp[1]+1, ns*nlambda*sizeof(fp_t));
    del_ft2dim(Stemp,1,ns,1,nlambda);
    delete temp_o;

    Bx[x1l][x2l][x3i] -= d_Bx;
    By[x1l][x2l][x3i] -= d_By;
    Bz[x1l][x2l][x3i] -= d_Bz;
    temp_o = obs_stokes(theta, phi, lambda, nlambda);
    // Now we need to subtract, this is unfortunately more complicated:
    Stemp = temp_o->get_S(1,1);
    for (int s=1;s<=ns;++s)
      for (int l=1;l<=nlambda;++l){
        d_obs[7][x3i][s][l] -= Stemp[s][l];
        d_obs[7][x3i][s][l] /= delta_angle;
      }
    delete temp_o;
    del_ft2dim(Stemp,1,ns,1,nlambda);
    Bx[x1l][x2l][x3i] += 0.5*d_Bx;
    By[x1l][x2l][x3i] += 0.5*d_By;
    Bz[x1l][x2l][x3i] += 0.5*d_Bz;
   */
    io.msg(IOL_INFO, "atmosphere::obs_stokes_num_responses: computed perturbation with respect to x3i = %d \n", x3i);
      
  }

  // Please write it down:

  fp_t * lambda_air = new fp_t [nlambda];
  for (int l=0;l<nlambda;++l)
    lambda_air[l] = vactoair(lambda[l]);

  if (intensity_responses){
    for (int param=1;param<=7;++param)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        for (int l=1;l<=nlambda;++l)
          for (int s=1;s<=ns;++s)
            intensity_responses[param][x3i][l][s] = d_obs[param][x3i][s][l];
  }
  
  
  FILE * out;
  out = fopen("stokes_intensity_responses_fin_diff.txt", "w");
  for (int param=1;param<=7;++param){
    fp_t norm = 1.0;
    for (int x3i=x3l;x3i<=x3h;++x3i){
      //norm = (param == 2) ? Nt[x1l][x2l][x3i] : 1.0;  
      for (int l=1;l<=nlambda;++l){
        fp_t loc = 0;
        if (tau_grid) loc = log10(-tau_referent[x1l][x2l][x3i]); else loc = x3[x3i];
        fprintf(out,"%10.10e %10.10e", loc, lambda_air[l-1]);
        for (int s=1;s<=ns;++s)
          fprintf(out," %10.10e", d_obs[param][x3i][s][l]*norm);
        fprintf(out," \n");
      }
    }
  }
  fclose(out);
  //exit(1);

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  io.msg(IOL_INFO, "atmosphere::obs_scalar_num_responses: computed responses of the intensity numericaly...\n");
  
  delete []lambda_air;
  del_ft4dim(d_obs,1,7,x3l,x3h,1,ns,1,nlambda);

  return o;

}