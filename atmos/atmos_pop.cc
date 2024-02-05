#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <fitsio.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "const.h"
#include <cmath>

#include "atmol/atmol.h"

#include "atmos.h"
#include "mathtools.h"


#include "../cmlib/fits_read.h"

void atmosphere::popsetup(void) // setup essential quantities
{
  for(int a=0;a<natm;++a) atml[a]->popsetup(x1l,x1h,x2l,x2h,x3l,x3h);
}

void atmosphere::popclean(void) // release memory
{
  for(int a=0;a<natm;++a) atml[a]->popclean(x1l,x1h,x2l,x2h,x3l,x3h);
}

void atmosphere::ltepops(void) // compute the populations in LTE
// **************************************************************************
// * Compute the LTE populations consistently with the chemical equilibrium *
// **************************************************************************
{
  for(int x1i=x1l;x1i<=x1h;++x1i) // Saha + Molecular 
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        chemeq(atml,natm,T[x1i][x2i][x3i],Nt[x1i][x2i][x3i],Ne[x1i][x2i][x3i],x1i,x2i,x3i); 
  }
  for(int a=0;a<natm;++a)
    atml[a]->lte(T,Ne); // Boltzmann

  io.msg(IOL_INFO,"atmosphere::ltepops:\n");  
  for(int a=0;a<natm;++a) atml[a]->info(); // Not sure this is needed but why not
}

fp_t atmosphere::ne_derivative(int x1i, int x2i, int x3i){

  // Compute the derivative of the electron density, using LTE chemeq

  fp_t derivative = 0;
  
  T[x1i][x2i][x3i] += delta_T * 0.5;
  chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i],Ne[x1i][x2i][x3i],x1i,x2i,x3i);
  derivative = Ne[x1i][x2i][x3i];
  T[x1i][x2i][x3i] -= delta_T;
  chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i],Ne[x1i][x2i][x3i],x1i,x2i,x3i);
  derivative -= Ne[x1i][x2i][x3i];
  derivative /= delta_T;
  T[x1i][x2i][x3i] += delta_T * 0.5;
  chemeq(atml, natm, T[x1i][x2i][x3i], Nt[x1i][x2i][x3i],Ne[x1i][x2i][x3i],x1i,x2i,x3i);
  
  return derivative;
}

int atmosphere::nltepops(void) // compute the NLTE populations (polarization free for now)
{


// Start by initializing everything to lte pops
  ltepops();
  int outcome = 0; // Returns 0 for all good, 1 for negative pops

  fp_t *lambda;
  int32_t nlambda=0;
  
  // Set-up transitions and get the lambda grid for the transitions
  for(int a=0;a<natm;++a) atml[a]->rtsetup(x1l,x1h,x2l,x2h,x3l,x3h); // initialize angular/wavelength redist/integration
  for(int a=0;a<natm;++a) lambda=atml[a]->getlambda(lambda,nlambda,T[x1l][x2l][x3h],Nt[x1l][x2l][x3h],Ne[x1l][x2l][x3h]); // compute wavelength grid for NLTE populations
  
  // If it turns out there are no wavelengths where RT needs to be done, we are done! 
  if (nlambda == 0){

    io.msg(IOL_INFO, "atmosphere::nltepops : no iterative procedure necessary, no nlte processes.\n");
    for(int a=0;a<natm;++a) atml[a]->rtclean(0,nlambda,x1l,x1h,x2l,x2h,x3l,x3h);
    return 0; 
  }

// radiative quantities
  fp_t ***S=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // Monochromatic intensity in given direction
  fp_t ***L=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // Monochromatic approximate operator for a given direction

// setup geometry specific angular quadrature, which is in 1D case gaussian quadrature, typically with 3 directions from -1,0 and viceversa.
  int ntp;
  fp_t *th,*ph;
  fp_t *bin=anglesetup(th,ph,ntp); // bin - integration weights
                                   // th - theta
                                   // ph - phi
                                   // ntp - number of angles, now as this is 1D, unpolarized case, it is the number of theta points. In general case it is N_theta x N_phi

  io.msg(IOL_INFO,"atmosphere::nltepops: lambda=%E .. %E nlambda = %d \n",lambda[0],lambda[nlambda-1], nlambda);
  
  // Set up integrating weights for the wavelength
  fp_t * lambda_w = new fp_t [nlambda];
  if (nlambda >= 2){
    lambda_w[0] = (lambda[1] - lambda[0]) * 0.5;
    lambda_w[nlambda - 1] = (lambda[nlambda-1] - lambda[nlambda-2]) * 0.5;
    for (int l = 1; l<nlambda-1; ++l)
      lambda_w[l] = (lambda[l+1] - lambda[l-1]) * 0.5;
  }
  else if (nlambda)
    lambda_w[0] = 1.0;
  
  io.msg(IOL_INFO, "atmosphere::nltepops: rt setup performed\n");

// NLTE population loop
  int32_t max_iter = 100;
  fp_t relative_change = 1.0;
  
  if (tau_grid) compute_op_referent();
  
  for (int a=0;a<natm;++a){
    atml[a]->compute_active_population(T, Ne);
  }
  
  int32_t iter = 0;

  // We treat the background opacity as constant, it is still wavelength and angle dependent:
  fp_t ***** op_background;
  fp_t ***** em_background;
  op_background = new fp_t****[ntp]-1;
  em_background = new fp_t****[ntp]-1;

  for (int tp=1;tp<=ntp;++tp){
    op_background[tp] = new fp_t***[nlambda];
    em_background[tp] = new fp_t***[nlambda];
  }

  for (iter = 0; iter<max_iter; ++iter){
      
    for(int a=0;a<natm;++a) atml[a]->rtinit();         // clear transition parameters for each atom
    
    for(int tp=1;tp<=ntp;++tp){                        // angular loop (loop over the rays)
        
        fp_t ***Vr=project(Vx,Vy,Vz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h);   // LOS projected velocity
        fp_t ****B=transform(Bx,By,Bz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h); // transform to (B,Inc,Az) representation

        for(int l=0;l<nlambda;++l){

          for (int a = 0; a<natm; ++a)
            atml[a]->prof_init();

          // If it is the first iteration calculate the background opacities, otherwise re-use old ones.
          if (iter==0){
            op_background[tp][l] = opacity_lte(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
            em_background[tp][l] = emissivity_lte(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
          }

          // This now implies that the opacity and emissivity are 1D only, and scalar:
          fp_t ***op = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
          fp_t ***em = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
          memcpy(op[x1l][x2l]+x3l,op_background[tp][l][x1l][x2l]+x3l,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
          memcpy(em[x1l][x2l]+x3l,em_background[tp][l][x1l][x2l]+x3l,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
          
          op=add(opacity_active_only(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]),op,x1l,x1h,x2l,x2h,x3l,x3h);
          em=add(emissivity_active_only(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]),em,x1l,x1h,x2l,x2h,x3l,x3h);

          if (tau_grid) normalize_to_referent_opacity(op,em);

          if(int rv=formal(rt_grid, S,L,op,em,th[tp],ph[tp], boundary_condition_for_rt)){ // solution and approximate operator
            io.msg(IOL_ERROR,"atmosphere::nltepops: for angle (%d), wavelength %d\n",tp,l);
            del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(L,x1l,x1h,x2l,x2h,x3l,x3h);
            return rv; // pass error to parent level
          }  
          if (tau_grid) de_normalize(op,em);  

          for(int a=0;a<natm;++a) atml[a]->add(S, L, op, lambda[l], lambda_w[l], bin[tp]); // give each species access to radiation field, that is, add the radiation field to the mean intensity
          del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h); // cannot be reused due to projections
          del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
        }
        del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
        del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);

      }

    // Update according to ALI
    relative_change = newpops(T,Nt,Ne,lambda,nlambda,1);
    io.msg(IOL_INFO, "\n atmosphere::nltepops : relative change after iteration %d is %.10e \n", iter, relative_change); 
    fprintf(stderr, "atmosphere::nltepops : relative change after iteration %d is %.5e \n", iter, relative_change);  
    
    if (relative_change < 0.0){
      // The populations are negative. Try normal lambda iteration:
      fprintf(stderr, "atmosphere::nltepops : found negative populations. trying lambda iteration...");
      relative_change = newpops(T,Nt,Ne,lambda,nlambda,0);
      if (relative_change < 0.0){ // If it is still bad:
        fprintf(stderr, "atmosphere::nltepops : found negative populations. could not fix. exiting...");
        outcome = -1;
        break;
      }
    }

    if (relative_change < 1E-2)
      break; 
  }
  io.msg(IOL_INFO, "atmosphere::nltepops : converged\n"); 
  
  for(int a=0;a<natm;++a) atml[a]->rtclean(ntp,nlambda,x1l,x1h,x2l,x2h,x3l,x3h); // uninitialize angular/wavelength redist/integration
  // cleanup background opacity
  for (int tp=1;tp<=ntp;++tp){
    for (int l=0;l<nlambda;++l){
      del_ft3dim(op_background[tp][l],x1l,x1h,x2l,x2h,x3l,x3h);
      del_ft3dim(em_background[tp][l],x1l,x1h,x2l,x2h,x3l,x3h);
    }
    delete[](op_background[tp]);
    delete[](em_background[tp]);
  }
  delete[](op_background+1);
  delete[](em_background+1);

// cleanup angular quadrature arrays
  delete[] (bin+1);
  delete[] (ph+1);
  delete[] (th+1);
  if (nlambda){
    delete[]lambda_w;
    delete[] lambda;
  }
// 
  del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(L,x1l,x1h,x2l,x2h,x3l,x3h);

  io.msg(IOL_INFO, "atmosphere::nltepops : solution converged in %6d iterations. Relative change is: %e \n", iter, relative_change); 
  return outcome;
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------

int atmosphere::atm_pop_setup(){


  // Find the total number of levels in the atmosphere:
  if (use_atm_lvls){
  
    //fprintf(stderr, "atmosphere::atm_pop_setup Total number of levels to consider is : %d \n", n_lvls);
    if (n_lvls){
      atm_lvl_pops = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls);
      use_atm_lvls = 2;
    }
    else
      atm_lvl_pops = 0;

  }
  //fprintf(stderr, "atmosphere::atm_pop_setup allocated \n");
  return 0;
}

int atmosphere::atm_pop_fill(){

  if (use_atm_lvls && atm_lvl_pops){
    for (int x1i=x1l; x1i<=x1h; ++x1i)
      for (int x2i=x2l; x2i<=x2h; ++x2i)
        for (int x3i=x3l; x3i<=x3h; ++x3i){
            int i = 1;
            for (int a=0;a<natm;++a)
              for (int z=0; z<atml[a]->get_no_ions(); ++z)
                for (int n=0; n<atml[a]->get_no_lvls(z); ++n){
                  atm_lvl_pops[x1i][x2i][x3i][i] = atml[a]->get_pop(x1i,x2i,x3i,z,n);
                  i++;
                }
            atm_lvl_pops[x1i][x2i][x3i][n_lvls] = Ne[x1i][x2i][x3i];
    }
  }
  return 0;
}

int atmosphere::atm_pop_invfill(){

  for (int x1i=x1l; x1i<=x1h; ++x1i)
    for (int x2i=x2l; x2i<=x2h; ++x2i)
      for (int x3i=x3l; x3i<=x3h; ++x3i){
        int i = 1; // Counts the levels
          for (int a=0;a<natm;++a)
            for (int z=0; z<atml[a]->get_no_ions(); ++z)
              for (int n=0; n<atml[a]->get_no_lvls(z); ++n){
                atml[a]->set_pop(x1i,x2i,x3i,z,n, atm_lvl_pops[x1i][x2i][x3i][i]);
                i++;
              }
  }
  return 0;
}

int atmosphere::atm_pop_output(){ // Posibbillity to output directly the data?

  if (use_atm_lvls == 2 && n_lvls){

    return 0;
  }
  return 0;
}

fp_t ** atmosphere::get_atm_pop(){

  if (use_atm_lvls == 2 && n_lvls){

    fp_t ** output = ft2dim(1,x3h-x3l+1,1,n_lvls);
    memcpy(output[1]+1,atm_lvl_pops[x1l][x2l][x3l]+1,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*n_lvls*sizeof(fp_t));
    return output;
  }
  return 0;
}

int atmosphere::atm_pop_clean(){

  if (use_atm_lvls == 2 && n_lvls){
    del_ft4dim(atm_lvl_pops, x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls);
    use_atm_lvls = 1;
  }
  return 0;
}

int atmosphere::get_total_atomic_levels(){
  return n_lvls;
}

int atmosphere::get_atm_pop_switch(){
  return use_atm_lvls;
}

fp_t atmosphere::newpops(fp_t ***T_in,fp_t ***Nt_in,fp_t ***Ne_in,fp_t *lambda,int32_t nlambda, int alo)

// *************************************************************************
// * solve the NLTE populations consistently with the chemical equilibrium *
// * This is necessary because the ionization fractions couple back to the *
// * particle number densities through the electron density and thus the   *
// * chemical equilibrium. This coupling may need to be explicitly treated *
// * to stabilize the solution or improve efficiency                       *
// *************************************************************************
// MvN has a good point here. But are still not there. We currently solve ONLY NLTE populations

{
  if (n_lvls){

    fprintf(stderr, "atmosphere:newpops: level of atoms that are considered is: %d \n", n_lvls);

    // Keeps relative changes - memory is cheap:
    fp_t **** relative_changes = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls);
    memset(relative_changes[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*n_lvls*sizeof(fp_t));

    // Using atm_lvl_pops to keep the proposed new values:
    atm_lvl_pops = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls);

    // This checks if we have negative populations.
    bool negative = false; 
  
    for(int x1i=x1l;x1i<=x1h;++x1i)
      for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i){ 
        
          int lvl_counter = 1;
          
          for(int a=0;a<natm;++a){
          
            // Get new populations for specific species:
            fp_t * new_pops_a = atml[a]->newpops(x1i, x2i, x3i, alo);
            
            for (int i=0; i<atml[a]->get_total_lvls();++i){
              // Propose that new pop
              atm_lvl_pops[x1i][x2i][x3i][lvl_counter] = new_pops_a[i];
              // Calculate the relative change:
              fp_t relative_change = new_pops_a[i] / atml[a]->get_mapped_pop(x1i,x2i,x3i,i) - 1.0;
              // Check for negative populations:
              if (relative_change < -1.0) negative = true;
              relative_changes[x1i][x2i][x3i][lvl_counter] = fabs(relative_change);
              lvl_counter ++;
            }
            delete []new_pops_a;
          }
    }
    fp_t max_rel_change = max_2d(relative_changes[x1l][x2l],x3l,x3h,1,n_lvls);
    fprintf(stderr, "atmos::newpops: max relative change is %e \n", max_rel_change);

    // Copy stuff:
    if (conserve_charge){
      //fprintf(stderr, "%e \n", Ne[x1l][x2l][44]);
      fprintf(stderr, "atmos::newpops: attempting to conserve charge...\n");
      enforce_conserve_charge();
      //fprintf(stderr, "%e \n", Ne[x1l][x2l][44]);
    }
    if (!negative){
      fprintf(stderr, "atmos::newpops: populations seem fine, updating...\n");
      atm_pop_invfill();
    }
    del_ft4dim(atm_lvl_pops,x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls);
    del_ft4dim(relative_changes,x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls);
    if (negative){ // If there are negative populations return error and tell the above the stop
      fprintf(stderr, "atmos::newpops: negative population. stopping...\n");
      return -1;
    }
    return max_rel_change;
  }
  return 0;
}

fp_t atmosphere::get_pop(int x1i, int x2i, int x3i, int species, int z, int i){

  return atml[species]->get_pop(x1i, x2i, x3i, z, i);
}

fp_t atmosphere::get_pop(int x1i, int x2i, int x3i, int species, int z){

  return atml[species]->get_pop(x1i, x2i, x3i, z);
}


fp_t atmosphere::get_Ne(int x1i, int x2i, int x3i){
  return Ne[x1i][x2i][x3i];
}

fp_t atmosphere::set_Ne(int x1i, int x2i, int x3i, fp_t input){
  Ne[x1i][x2i][x3i] = input;
  return 0;
}

int atmosphere::enforce_conserve_charge(){

  // Adjust electron populations according to the old and the new level populations

  for (int x1i=x1l; x1i<=x1h; ++x1i)
    for (int x2i=x2l; x2i<=x2h; ++x2i)
      for (int x3i=x3l; x3i<=x3h; ++x3i){
        int i = 1; // Counts the levels
          for (int a=0;a<natm;++a)
            for (int z=0; z<atml[a]->get_no_ions(); ++z) 
              for (int n=0; n<atml[a]->get_no_lvls(z); ++n){
                // atm_lvl_pops are the new ones, the other ones are the old ones
                fp_t ne_diff = z*(atm_lvl_pops[x1i][x2i][x3i][i] - atml[a]->get_pop(x1i,x2i,x3i,z,n));
                Ne[x1i][x2i][x3i] += ne_diff;
                //fprintf(stderr,"%d %d %d %d %d %e \n", x1i,x2i,x3i,z,n,ne_diff);
                i++;
              }

  }
  return 0;
}