#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "const.h"
#include <cmath>

#include "atmol/atmol.h"

#include "atmos.h"
#include "mathtools.h"

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
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        chemeq(atml,natm,T[x1i][x2i][x3i],Nt[x1i][x2i][x3i],Ne[x1i][x2i][x3i],x1i,x2i,x3i); 
  }
  for(int a=0;a<natm;++a)
    atml[a]->lte(T,Ne); // Boltzmann

  io.msg(IOL_INFO,"atmosphere::ltepops:\n");  
  for(int a=0;a<natm;++a) atml[a]->info();

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

  fp_t *lambda;
  int32_t nlambda=0;
  if (nlambda){
    lambda=new fp_t [nlambda];           // wavelengths for NLTE calculations, these are @ several pre-defined positions, for the continuum. 
    for (int l=0;l<nlambda;++l)
      lambda[l] = 500E-7 + (460E-7 - 50E-7) / (nlambda-1) * l;
  }
  if(nlambda == 1)lambda[0] = 500E-7;
  // Now we have some setup phase, where we pick wavelengths etc, that is, we try to cleverly sort out the wavelength grid required for our computations. For the moment it looks nice.
  for(int a=0;a<natm;++a) atml[a]->rtsetup(x1l,x1h,x2l,x2h,x3l,x3h); // initialize angular/wavelength redist/integration
  for(int a=0;a<natm;++a) lambda=atml[a]->getlambda(lambda,nlambda,T[x1l][x2l][x3h],Nt[x1l][x2l][x3h],Ne[x1l][x2l][x3h]); // compute wavelength grid for NLTE populations
 
  // If it turns out there are no wavelengths where RT needs to be done, we are done! 
  if (nlambda == 0){

    io.msg(IOL_INFO, "atmosphere::nltepops : no iterative procedure necessary, no nlte processes. we are done here.\n");
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
                                   // ntp - number of angles, now as this is 1D, unpolarized case, it is actually only number of theta points. In general case it is N_theta x N_phi

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
  int32_t max_iter = 200;
  fp_t relative_change = 1.0;
  
  if (tau_grid) compute_op_referent();

  for (int a=0;a<natm;++a){
    atml[a]->compute_active_population(T, Ne);
  }
  int32_t iter = 0;

  for (iter = 0; iter<max_iter; ++iter){
      
    for(int a=0;a<natm;++a) atml[a]->rtinit();         // clear transition parameters for each atom
    
    for(int tp=1;tp<=ntp;++tp){                        // angular loop (loop over the rays)
        
        fp_t ***Vr=project(Vx,Vy,Vz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h);   // LOS projected velocity
        fp_t ****B=transform(Bx,By,Bz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h); // transform to (B,Inc,Az) representation

        for(int l=0;l<nlambda;++l){

          for (int a = 0; a<natm; ++a)
            atml[a]->prof_init();

          fp_t ***op=opacity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
          fp_t ***em=emissivity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);

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

    relative_change = newpops(T,Nt,Ne,lambda,nlambda);

    //io.msg(IOL_INFO, "atmosphere::nltepops : relative change after iteration %d is %.10e \n", iter, relative_change); 
    //printf("atmosphere::nltepops : relative change after iteration %d is %.10e \n", iter, relative_change);  

    if (relative_change < 1E-3)
      break; 
  }
  io.msg(IOL_INFO, "atmosphere::nltepops : converged\n"); 
  //fclose(test_intensity);
  //fclose(test_opem);
  
  for(int a=0;a<natm;++a) atml[a]->rtclean(ntp,nlambda,x1l,x1h,x2l,x2h,x3l,x3h); // uninitialize angular/wavelength redist/integration

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
}

int atmosphere::nltepops_taugrid(void){
  
  // radiative quantities
  fp_t ***S=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // Monochromatic intensity in given direction
  fp_t ***L=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // Monochromatic approximate operator for a given direction

//
  int32_t nlambda=10;
  fp_t *lambda=new fp_t [nlambda];           // wavelengths for NLTE calculations, these are several pre-defined positions, for the continuum. 
  for (int l=0;l<nlambda;++l)
    lambda[l] = 300E-7 + (2000E-7 - 300E-7) / (nlambda-1) * l;

// setup geometry specific angular quadrature, which is in 1D case gaussian quadrature, typically with 3 directions from 0 to 1.
  int ntp;
  fp_t *th,*ph;
  fp_t *bin=anglesetup(th,ph,ntp); // bin - integration weights
                                   // th - theta
                                   // ph - phi
                                   // ntp - number of angles, now as this is 1D, unpolarized case, it is actually only number of theta points. In general case it is N_theta x N_phi

  io.msg(IOL_INFO,"atmosphere::nltepops_taugrid:\n");

  //Now we have some setup phase, where we pick wavelengths etc, that is, we try to cleverly sort out the wavelength grid required for our computations. For the moment it looks nice.
  for(int a=0;a<natm;++a) atml[a]->rtsetup(x1l,x1h,x2l,x2h,x3l,x3h); // initialize angular/wavelength redist/integration
  for(int a=0;a<natm;++a) lambda=atml[a]->getlambda(lambda,nlambda,T[x1l][x2l][x3h],Nt[x1l][x2l][x3h],Ne[x1l][x2l][x3h]); // compute wavelength grid for NLTE populations
 
  io.msg(IOL_INFO,"atmosphere::nltepops: lambda=%E .. %E nlambda = %d \n",lambda[0],lambda[nlambda-1], nlambda);
  
  // Set up integrating weights for the wavelength
  fp_t * lambda_w = new fp_t [nlambda];
  lambda_w[0] = (lambda[1] - lambda[0]) * 0.5;
  lambda_w[nlambda - 1] = (lambda[nlambda-1] - lambda[nlambda-2]) * 0.5;
  for (int l = 1; l<nlambda-1; ++l)
    lambda_w[l] = (lambda[l+1] - lambda[l-1]) * 0.5;
  
  io.msg(IOL_INFO, "atmosphere::nltepops: rt setup performed\n");

// NLTE population loop
  int32_t max_iter = 300;
  int32_t iter = 0;
  fp_t relative_change = 1.0;
  
  ltepops();
  compute_op_referent();

  for (int a=0;a<natm;++a){
    atml[a]->compute_active_population(T, Ne);
  }
  
  for (iter = 0; iter<max_iter; ++iter){
      
    for(int a=0;a<natm;++a) atml[a]->rtinit();         // clear transition parameters for each atom

    // You need to compute relevant opacity, but not the tau. We want tau to be independent
    
    for(int tp=1;tp<=ntp;++tp){                        // angular loop (loop over the "wavefronts")
        
        fp_t ***Vr=project(Vx,Vy,Vz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h);   // LOS projected velocity
        fp_t ****B=transform(Bx,By,Bz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h); // transform to (B,Inc,Az) representation

        for(int l=0;l<nlambda;++l){

          for (int a = 0; a<natm; ++a)
            atml[a]->prof_init();

          fp_t ***op=opacity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]); // angle dependent in the scalar case -> NO, at leasst not at the moment. 
          fp_t ***em=emissivity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);

          fp_t *** op_norm = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
          fp_t *** em_norm = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
          for (int x1i=x1l;x1i<=x1h;++x1i)
            for (int x2i=x2l;x2i<=x2h;++x2i)
              for (int x3i=x3l;x3i<=x3h;++x3i){
                op_norm[x1i][x2i][x3i] = op[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i];
                em_norm[x1i][x2i][x3i] = em[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i];
              }
          
          if(int rv=formal(tau_referent[x1l][x2l], S,L,op_norm, em_norm,th[tp],ph[tp], boundary_condition_for_rt)){ // solution and approximate operator
            io.msg(IOL_ERROR,"atmosphere::nltepops: for angle (%d), wavelength %d\n",tp,l);
            del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(op_norm,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(em_norm,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
            del_ft3dim(L,x1l,x1h,x2l,x2h,x3l,x3h);
            return rv; // pass error to parent level
          }   
          //for (int x3i=x3l;x3i<=x3h;++x3i)
            //printf("%d %e \n", x3i, S[x1l][x2l][x3i]);  

          for(int a=0;a<natm;++a) atml[a]->add(S, L, op, lambda[l], lambda_w[l], bin[tp]); // give each species access to radiation field, that is, add the radiation field to the mean intensity
          del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h); // cannot be reused due to projections
          del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
          del_ft3dim(op_norm,x1l,x1h,x2l,x2h,x3l,x3h);
          del_ft3dim(em_norm,x1l,x1h,x2l,x2h,x3l,x3h);
        }
        del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
        del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);

      }
    
    relative_change = newpops(T,Nt,Ne,lambda,nlambda);

    //io.msg(IOL_INFO, "atmosphere::nltepops : relative change after iteration %d is %.10e \n", iter, relative_change);   

    if (relative_change < 1E-5)
      break; 
  }
  
  for(int a=0;a<natm;++a) atml[a]->rtclean(ntp,nlambda,x1l,x1h,x2l,x2h,x3l,x3h); // uninitialize angular/wavelength redist/integration
  delete[] lambda;
// cleanup angular quadrature arrays
  delete[] (bin+1);
  delete[] (ph+1);
  delete[] (th+1);
  delete[]lambda_w;
// 
  del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(L,x1l,x1h,x2l,x2h,x3l,x3h);

  io.msg(IOL_INFO, "atmosphere::nltepops : solution converged in %6d iterations. Relative change is: %e \n", iter, relative_change); 

}

fp_t atmosphere::newpops(fp_t ***T_in,fp_t ***Nt_in,fp_t ***Ne_in,fp_t *lambda,int32_t nlambda)
// *************************************************************************
// * solve the NLTE populations consistently with the chemical equilibrium *
// * This is necessary because the ionization fractions couple back to the *
// * particle number densities through the electron density and thus the   *
// * chemical equilibrium. This coupling may need to be explicitly treated *
// * to stabilize the solution or improve efficiency                       *
// *************************************************************************

{

  //io.msg(IOL_INFO,"atmosphere::newpops \n");

  fp_t *** convergence = ft3dim(1, 1, 1, 1, x3l, x3h);
  int * rel_change_index = new int [x3h-x3l+1] - x3l;
  memset(convergence[1][1]+1, 0, (x3h - x3l+1)*sizeof(fp_t));

  uint08_t grv=0,rv;
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){ // the do-loop may be done as outer loop also, but this may be more efficient. Milic: It is probably more efficient because you do not "wait" for all the points to finish.
        
        do{   // NLTE <-> Chemical Equilibrium (see above)
          //rv=chemeq(atml,natm,T_in[x1i][x2i][x3i],Nt_in[x1i][x2i][x3i],Ne_in[x1i][x2i][x3i],x1i,x2i,x3i); // solve chemical and statistical equilibrium, using "known" ionization fractions and their derivatives

          fp_t * changes = new fp_t [natm];
          memset(changes, 0, natm*sizeof(fp_t));
          
          for(int a=0;a<natm;++a)
            changes[a] = atml[a]->pops(atml,natm,T_in[x1i][x2i][x3i],Ne_in[x1i][x2i][x3i],x1i,x2i,x3i); // solve the rate equations for each atom
          //printf("Change of atom %d is %e \n", 0, changes[0]);
          

          convergence[1][1][x3i] = max_1d(changes, 0, natm-1);
          rel_change_index[x3i] = x3i;//max_1d_index(changes, 0, natm-1);
          //printf("! %d %e %d \n", x3i, convergence[1][1][x3i], rel_change_index[x3i]);
          delete []changes;

          //rv=chemeq(atml,natm,T_in[x1i][x2i][x3i],Nt_in[x1i][x2i][x3i],Ne_in[x1i][x2i][x3i],x1i,x2i,x3i); // solve chemical and statistical equilibrium, using "known" ionization fractions and their derivatives
          rv = 0.0;
          
        }while(rv&EC_POP_ION_DEFECT); // chemical consistency when ionization fraction is consistent
      }

  fp_t r_c = max_1d(convergence[1][1], x3l, x3h);
  int max_index = max_1d_index(convergence[1][1], x3l, x3h);

  //printf ("Relative change in this iteration is: %e at point %d \n", r_c, max_index);

  del_ft3dim(convergence, 1, 1, 1, 1, x3l, x3h);
  delete [] (rel_change_index+x3l);
  //io.msg(IOL_INFO, "Populations updated! \n");
  //return grv; // return flags concerning any remaining defects
  return r_c;
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