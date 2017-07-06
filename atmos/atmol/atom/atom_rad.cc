#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>

#include "types.h"
#include "uts.h"
#include "io.h"
#include "pack.h"
#include "struts.h"
#include "mem.h"
#include "profiles.h"
#include "const.h"
#include "atomcfg.h"
#include "mathtools.h"

#include "H/hydrogen.h"
#include "He/helium.h"

#include "atomcfg.h"
#include "partf.h"  // partition functions
#include "bfcs.h"   // bound-free crossections
#include "colr.h"   // collisional rates
#include "atom.h"
// --------------------------------------------------------------------------------------------------------------------------------------------------

// Total opacity and emissivity functions. They both simply add three different contributions from three kinds of processes:

fp_t ***atom::opacity(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda)
{
  fp_t ***op=freefree_op(T,Ne,Vlos,lambda); // free-free opacity
  op=add(rayleigh_op(lambda),op,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  op=add(boundfree_op(Vlos,lambda),op,x1l,x1h,x2l,x2h,x3l,x3h);          // bound-free ionization
  op=add(boundbound_op(T,Ne,Vlos,Vt, B_vec, lambda),op,x1l,x1h,x2l,x2h,x3l,x3h); // bound-bound transitions
  
  return op;
}


fp_t ***atom::emissivity(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda)
{
  fp_t ***em=freefree_op(T,Ne,Vlos,lambda);
  em=add(rayleigh_op(lambda),em,x1l,x1h,x2l,x2h,x3l,x3h);
  
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i)
            em[x1i][x2i][x3i] *= Planck_f(lambda, T[x1i][x2i][x3i]);
  memset(em[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));  
  
  em=add(boundfree_em(Vlos,lambda),em,x1l,x1h,x2l,x2h,x3l,x3h);
  em=add(boundbound_em(T,Ne,Vlos,Vt, B_vec, lambda),em,x1l,x1h,x2l,x2h,x3l,x3h);

  return em;
}

fp_t ***atom::emissivity_polarized_dummy(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda, fp_t)
{
  fp_t ***em=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  for(int z=0;z<=Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t Aul = A[z][i][ii];
        for(int x1i=x1l;x1i<=x1h;++x1i)
          for(int x2i=x2l;x2i<=x2h;++x2i)
            for(int x3i=x3l;x3i<=x3h;++x3i){
              fp_t profile =  current_profile[x1i][x2i][x3i][tmap[z][i][ii]];
              em[x1i][x2i][x3i] += profile * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi * 
                (J_02[x1i][x2i][x3i][tmap[z][i][ii]]) * (1.0 - cos(theta) * cos(theta)) * sqrt(9./8.);
            }
      }       
  return em;
  
}

fp_t ***** atom::emissivity_polarized_perturbation_dummy(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda, fp_t){

  fp_t *****em_pert=ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  for(int z=0;z<=Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t Aul = A[z][i][ii];
        int upper_map = rmap[z][i];
        for(int x3k=x3l;x3k<=x3h;++x3k)
          for(int x1i=x1l;x1i<=x1h;++x1i)
            for(int x2i=x2l;x2i<=x2h;++x2i)
              for(int x3i=x3l;x3i<=x3h;++x3i){
                fp_t profile =  current_profile[x1i][x2i][x3i][tmap[z][i][ii]];
                em_pert[1][x3k][x1i][x2i][x3i] += profile *  level_responses[1][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * h * c / lam / 4.0 / pi;
                em_pert[2][x3k][x1i][x2i][x3i] += profile *  level_responses[2][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * h * c / lam / 4.0 / pi;
                em_pert[3][x3k][x1i][x2i][x3i] += profile *  level_responses[3][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * h * c / lam / 4.0 / pi;

                if (x3i == x3k){
                  
                  fp_t derivative = voigt_der_T(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  em_pert[1][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi;
                  derivative = voigt_der_density(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  em_pert[2][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi;
                  derivative = voigt_der_vt(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  em_pert[3][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi;
                }
                em_pert[1][x3k][x1i][x2i][x3i] *= J_02[x1i][x2i][x3i][tmap[z][i][ii]];
                em_pert[2][x3k][x1i][x2i][x3i] *= J_02[x1i][x2i][x3i][tmap[z][i][ii]];
                em_pert[3][x3k][x1i][x2i][x3i] *= J_02[x1i][x2i][x3i][tmap[z][i][ii]];

                em_pert[1][x3k][x1i][x2i][x3i] += pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi * profile * J_02_responses[1][x3k][x3i][tmap[z][i][ii]];
                em_pert[2][x3k][x1i][x2i][x3i] += pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi * profile * J_02_responses[2][x3k][x3i][tmap[z][i][ii]];
                em_pert[3][x3k][x1i][x2i][x3i] += pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi * profile * J_02_responses[3][x3k][x3i][tmap[z][i][ii]];

                em_pert[5][x3k][x1i][x2i][x3i] += pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi * profile * J_02_responses[5][x3k][x3i][tmap[z][i][ii]];
            }
      }       
  return em_pert;
  
  
}


// Perturbation variants:

fp_t ***** atom::opacity_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda)
{
  fp_t ***** op_pert = boundbound_op_pert(T,Ne,Vlos,Vt, B_vec, lambda);
  op_pert=add(boundfree_op_pert(Vlos,lambda),op_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  
  return op_pert;
}


fp_t ***** atom::emissivity_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda)
{
  fp_t ***** em_pert=boundbound_em_pert(T,Ne,Vlos,Vt, B_vec, lambda);
  em_pert=add(boundfree_em_pert(Vlos,lambda),em_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  return em_pert;

}

fp_t atom::opacity_continuum(fp_t T_in, fp_t Ne_in, fp_t lambda, int x1i, int x2i, int x3i){
  // this is the function which directly adds b-f and ff opacity to compute continuum opacity

  fp_t op_cont = 0.0;

  // Add F-F
  fp_t fct=(e*e*c*h*h*Ryd)/(6*pi*pi*pi*e0*me*c)*sqrt((2.0*pi)/(3.0*k*me*me*me));
  fp_t gff = 1.1;
  fp_t a = fct * gff * Ne_in /sqrt(T_in)/c/c/c;
  for (int z=1;z<=Z;++z)
    op_cont += a*sqr((fp_t)z)*pop[x1i][x2i][x3i].N[z];
  fp_t b=-(h*c)/(k*T_in);
  op_cont *= lambda*sqr(lambda)*(1.0-exp(b/lambda));

  op_cont = 0.0;
  // And then directly add B-F
  for (int z=0;z<Z;++z)
    for (int l=0;l<nl[z];++l){
      
      fp_t sigma = (bf[z][l]) ? bf[z][l]->U(lambda) : 0.0;
      if (sigma){
        fp_t pop_mod = pop[x1i][x2i][x3i].n[z+1][0] * Ne_in * saha_const * pow(T_in, -1.5) * fp_t(g[z][l]) / fp_t(g[z+1][0]) * exp((ip[z] - ee[z][l]) / k /  T_in);
        op_cont += sigma * (pop[x1i][x2i][x3i].n[z][l] - pop_mod * exp(- h * c / lambda / k / T_in));
      }   
  }

  // And then add Rayleigh scattering:
  if (strcmp(id,"H")==0){ // if H  
    // Using expression found in Lee & Kim (2002), equation 16.
    // This is valid far from resonances
    fp_t ratio = 141E-7 / lambda;
    fp_t cross_section = 8.41E-25  * pow(ratio,4) + 3.37E-24 * pow(ratio,6) + 4.71E-22 * pow(ratio,14);
    //op_cont += cross_section * pop[x1i][x2i][x3i].N[0];
  }

  return op_cont;
}

fp_t ** atom::opacity_continuum_pert(fp_t T_in, fp_t Ne_in, fp_t lambda, int x1i, int x2i, int x3i){

  return 0;
 
}



// --------------------------------------------------------------------------------------------------------------------------------------------------

// Now we have free-free functions

fp_t ***atom::freefree_op(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t lambda)
// *************************************************************************
// * Free-Free opacities:
// * default: Hydrogenic: SI/CGS invariant?
// * Hydrogenic expression from Mihalas, 2ed. p101, with charge dependence
// * and Rydberg constants retained.
// *************************************************************************
{
  fp_t ***op=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
//
  fp_t fct=(e*e*c*h*h*Ryd)/(6.0*pi*pi*pi*e0*me*c)*sqrt((2.0*pi)/(3.0*k*me*me*me));
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        fp_t gff=1.1; // see Berger, ApJ 1956 for expressions
        fp_t a=fct*gff*Ne[x1i][x2i][x3i]/(sqrt(T[x1i][x2i][x3i])*c*c*c);
        for(int z=1;z<=Z;++z) op[x1i][x2i][x3i]+=a*sqr((fp_t)z)*pop[x1i][x2i][x3i].N[z];
// wavelength dependence
        fp_t b=-(h*c)/(k*T[x1i][x2i][x3i]); // stimulated emission correction
        fp_t ll=lambda*(1.0+Vlos[x1i][x2i][x3i]/c);
        op[x1i][x2i][x3i]*=ll*sqr(ll)*(1.0-exp(b/ll));
      }
  return op;
}

fp_t ***** atom::freefree_op_pert(fp_t ***T, fp_t ***Ne, fp_t ***Vlos, fp_t lambda){
  fp_t ***** op_pert = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  return op_pert;
}

fp_t ***atom::freefree_em(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t lambda)
// *************************************************************************
// bremsstrahlung: B[Te/Ti]/alpha?
// *************************************************************************
{
  fp_t ***em=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
//
  fp_t fct=(e*e*c*h*h*Ryd)/(6*pi*pi*pi*e0*me*c)*sqrt((2.0*pi)/(3.0*k*me*me*me));
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        fp_t gff=1.1; // see Berger, ApJ 1956 for expressions
        fp_t a=fct*gff*Ne[x1i][x2i][x3i]/(sqrt(T[x1i][x2i][x3i])*c*c*c);
        for(int z=1;z<=Z;++z) em[x1i][x2i][x3i]+=a*sqr((fp_t)z)*pop[x1i][x2i][x3i].N[z];
// wavelength dependence
        fp_t b=-(h*c)/(k*T[x1i][x2i][x3i]); // stimulated emission correction
        fp_t ll=lambda*(1.0+Vlos[x1i][x2i][x3i]/c);
        em[x1i][x2i][x3i]*=ll*sqr(ll)*(1.0-exp(b/ll));
      }
  return em;
}

fp_t ***** atom::freefree_em_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t lambda){
  
 return 0;
}

// Rayleigh scattering:

fp_t *** atom::rayleigh_op(fp_t lambda){

	fp_t *** op = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
	memset(op[x1l][x2l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  if (Z == 1){ // if Hydrogen
   
    fp_t cross_section = 0;
    // Add all the possible resonances:
    int i = 0;
    int z = 0;
    fp_t lambda_0 = 1215.67E-8; // Ly alpha wavelength
    
    // This up there was literally from Stellar Atmospheres 3rd edition. It led to singularities. Singularities are bad. 
    // Now we will try approach from Lee, 2005 (found in references from Stellar Atmospheres)

    fp_t cp[5] = {1.26563, 3.73828125, 8.813130, 19.15379502, 39.92303};
    for (int ii=0;ii<=5;++ii)
      cross_section += cp[i] * pow(lambda_0/lambda, 2*ii);

    cross_section *= pow(lambda_0/lambda, 4);

    cross_section *= cross_section * 6.65E-25;
    for (int x1i=x3l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          op[x1i][x2i][x3i] = cross_section * pop[x1i][x2i][x3i].N[z];
  }
  return op;

}

fp_t *** atom::rayleigh_em(fp_t lambda){

  // This we will treat in LTE so far

  return 0;

}

// --------------------------------------------------------------------------------------------------------------------------------------------------

// Bound-free functions

fp_t ***atom::boundfree_op(fp_t ***Vlos,fp_t lambda) // absorption
// *************************************************************************
// * Bound-Free opacities as specified for each level...                   *
// *************************************************************************
{
  fp_t ***op=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
//
  for(int z=0;z<Z;++z) // skip final stage: it cannot be ionized
    for(int l=0;l<nl[z];++l)
      for(int x1i=x1l;x1i<=x1h;++x1i)
        for(int x2i=x2l;x2i<=x2h;++x2i)
          for(int x3i=x3l;x3i<=x3h;++x3i){
            fp_t sigma = (bf[z][l]) ? bf[z][l]->U(lambda * (1.0 + 0.0*Vlos[x1i][x2i][x3i]/c)) : 0.0;

            if (sigma){

              fp_t T = fetch_temperature(x1i, x2i, x3i);
    
            // But we also need some corrections now.
              fp_t pop_mod = pop[x1i][x2i][x3i].n[z+1][0] * fetch_Ne(x1i, x2i, x3i) * saha_const * pow(T, -1.5) 
              * fp_t(g[z][l]) / fp_t(g[z+1][0]) * exp((ip[z] - ee[z][l] - h * c / lambda)/k/T);
              op[x1i][x2i][x3i] += sigma * (pop[x1i][x2i][x3i].n[z][l] - pop_mod);
            }
          }
  return op;
}

fp_t ***** atom::boundfree_op_pert(fp_t ***Vlos, fp_t lambda){

  fp_t ***** op_pert = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  for(int z=0;z<=Z;++z)
    for(int l=0;l<nl[z];++l)
      for(int x1i=x1l;x1i<=x1h;++x1i)
        for(int x2i=x2l;x2i<=x2h;++x2i)
          for(int x3i=x3l;x3i<=x3h;++x3i)
            for (int x3k=x3l;x3k<=x3h;++x3k){

              op_pert[1][x3k][x1i][x2i][x3i] += op_derivative_to_level(x3i, rmap[z][l], lambda) * 
                level_responses[1][(x3i-x3l)*nmap+rmap[z][l]+1][x3k];
              op_pert[2][x3k][x1i][x2i][x3i] += op_derivative_to_level(x3i, rmap[z][l], lambda) * 
                level_responses[2][(x3i-x3l)*nmap+rmap[z][l]+1][x3k];
              op_pert[3][x3k][x1i][x2i][x3i] += op_derivative_to_level(x3i, rmap[z][l], lambda) * 
                level_responses[3][(x3i-x3l)*nmap+rmap[z][l]+1][x3k];
  }
  // Add explicit part:
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        fp_t * explicit_pert = op_derivative_explicit(x3i,0,0,lambda);
        op_pert[1][x3i][x1i][x2i][x3i] += explicit_pert[1];
        op_pert[2][x3i][x1i][x2i][x3i] += explicit_pert[2];
        op_pert[3][x3i][x1i][x2i][x3i] += explicit_pert[3];
        delete [](explicit_pert+1);
  }     
  return op_pert;
}


fp_t ***atom::boundfree_em(fp_t ***Vlos,fp_t lambda) // emissivity
// *************************************************************************
// * Bound-Free recombination...                                           *
// * Only recombination to particles in the ground state is considered.    *
// * This is considered reasonable since recombination with an excited     *
// * particle is usually followed by immediate autoionization (not always) *
// *************************************************************************
{
  fp_t ***em=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
//  
  for(int z=0;z<Z;++z) // recombination to stage z
    for(int l=0;l<nl[z];++l)
      for(int x1i=x1l;x1i<=x1h;++x1i)
        for(int x2i=x2l;x2i<=x2h;++x2i)
          for(int x3i=x3l;x3i<=x3h;++x3i){
            
            fp_t sigma = (bf[z][l]) ? bf[z][l]->U(lambda * (1.0 + 0.0*Vlos[x1i][x2i][x3i]/c)) : 0;
            if (sigma){
              fp_t pop_mod = pop[x1i][x2i][x3i].n[z+1][0] * fetch_Ne(x1i, x2i, x3i) * saha_const * pow(fetch_temperature(x1i, x2i, x3i), -1.5) * fp_t(g[z][l]) / fp_t(g[z+1][0])
              * exp((ip[z] - ee[z][l]) / k /  fetch_temperature(x1i, x2i, x3i));          
              em[x1i][x2i][x3i] += sigma * pop_mod * (1.0 - exp(- h * c / lambda / k / fetch_temperature(x1i, x2i, x3i))) * Planck_f(lambda, fetch_temperature(x1i,x2i,x3i));
            }
  }
  return em;
}

fp_t ***** atom::boundfree_em_pert(fp_t *** Vlos, fp_t lambda){
  fp_t ***** em_pert = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  for(int z=0;z<=Z;++z)
    for(int l=0;l<nl[z];++l)
      for(int x1i=x1l;x1i<=x1h;++x1i)
        for(int x2i=x2l;x2i<=x2h;++x2i)
          for(int x3i=x3l;x3i<=x3h;++x3i)
            for (int x3k=x3l;x3k<=x3h;++x3k){

              em_pert[1][x3k][x1i][x2i][x3i] += em_derivative_to_level(x3i, rmap[z][l], lambda) * 
                level_responses[1][(x3i-x3l)*nmap+rmap[z][l]+1][x3k];
              em_pert[2][x3k][x1i][x2i][x3i] += em_derivative_to_level(x3i, rmap[z][l], lambda) * 
                level_responses[2][(x3i-x3l)*nmap+rmap[z][l]+1][x3k];
              em_pert[3][x3k][x1i][x2i][x3i] += em_derivative_to_level(x3i, rmap[z][l], lambda) * 
                level_responses[3][(x3i-x3l)*nmap+rmap[z][l]+1][x3k];
  }
  // Add explicit part:
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        fp_t * explicit_pert = em_derivative_explicit(x3i,0,0,lambda);
        em_pert[1][x3i][x1i][x2i][x3i] += explicit_pert[1];
        em_pert[2][x3i][x1i][x2i][x3i] += explicit_pert[2];
        em_pert[3][x3i][x1i][x2i][x3i] += explicit_pert[3];
        delete [](explicit_pert+1);
  }
  return em_pert;
}

// --------------------------------------------------------------------------------------------------------------------------------------------------

// Bound-bound functions

fp_t ***atom::boundbound_em(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t lambda) // emission
{ 

  fp_t ***em=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  for(int z=0;z<=Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t Aul = A[z][i][ii];
        for(int x1i=x1l;x1i<=x1h;++x1i)
          for(int x2i=x2l;x2i<=x2h;++x2i)
            for(int x3i=x3l;x3i<=x3h;++x3i){
              fp_t profile =  current_profile[x1i][x2i][x3i][tmap[z][i][ii]];

              fp_t em_current_transition = profile * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lambda / 4.0 / pi;
              em[x1i][x2i][x3i] += em_current_transition;
            }
      }       
  return em;
}

fp_t ***** atom::boundbound_em_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t lambda){
  // The same as above, with additional loop for perturbation:
  fp_t *****em_pert=ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  for(int z=0;z<=Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t Aul = A[z][i][ii];
        int upper_map = rmap[z][i];
        for(int x3k=x3l;x3k<=x3h;++x3k)
          for(int x1i=x1l;x1i<=x1h;++x1i)
            for(int x2i=x2l;x2i<=x2h;++x2i)
              for(int x3i=x3l;x3i<=x3h;++x3i){
                fp_t profile =  current_profile[x1i][x2i][x3i][tmap[z][i][ii]];
                em_pert[1][x3k][x1i][x2i][x3i] += profile *  level_responses[1][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * h * c / lambda / 4.0 / pi;
                em_pert[2][x3k][x1i][x2i][x3i] += profile *  level_responses[2][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * h * c / lambda / 4.0 / pi;
                em_pert[3][x3k][x1i][x2i][x3i] += profile *  level_responses[3][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * h * c / lambda / 4.0 / pi;

                if (x3i == x3k){
                  
                  fp_t derivative = voigt_der_T(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  em_pert[1][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lambda / 4.0 / pi;
                  derivative = voigt_der_density(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  em_pert[2][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lambda / 4.0 / pi;
                  derivative = voigt_der_vt(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  em_pert[3][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lambda / 4.0 / pi;
                }
            }
      }       
  return em_pert;
}


// ----------------------------------------------------------------------------------------------------------------------------
// This is another function which is being inserted here by Milic on @ 10/10/2014: 
// I like the way fp_t *** atom::boundbound_em is written and I want another one just like that just for opacity.
// Overloaded version - accepts the density of collisional partners as well:

fp_t ***atom::boundbound_op(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t lambda) // absorption
{ 

  fp_t ***op=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  for(int z=0;z<=Z;++z) // // PAY ATENTION HERE < or <=  ??????
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t ar=damp_rad(A[z],ii,i);
        fp_t gc = 0.0;
        fp_t Blu=B[z][ii][i],Bul=B[z][i][ii];
        fp_t sf=(lam*lam)/(4.0*pi*c);


        for(int x1i=x1l;x1i<=x1h;++x1i)
          for(int x2i=x2l;x2i<=x2h;++x2i)
            for(int x3i=x3l;x3i<=x3h;++x3i){
              
              fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);

              gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
              //if (Z==2)
              //  printf("%d %e %e \n", x3i,parent_atm->get_x3(x3i),gc);
              fp_t a=sf*(ar+gc*turn_on_damping)/dld;
              fp_t x=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
              fp_t profile = 1.0 / dld * fvoigt(x, a);
              current_profile[x1i][x2i][x3i][tmap[z][i][ii]] = current_profile[x1i][x2i][x3i][tmap[z][ii][i]] = profile;

              fp_t op_current_transition = profile * pop[x1i][x2i][x3i].n[z][ii] * (Blu) * h * c / lambda / 4.0 / pi * 
                (1.0  - pop[x1i][x2i][x3i].n[z][i] / pop[x1i][x2i][x3i].n[z][ii] * Bul / Blu);
              
              op[x1i][x2i][x3i]+= op_current_transition;
            }
      }
  //printf("op_unpol = %e \n", op[x1l][x2l][25]);
  return op;
}

fp_t ***** atom::boundbound_op_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t lambda){
  fp_t *****op_pert=ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  for(int z=0;z<=Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t Bul = B[z][i][ii];
        fp_t Blu = B[z][ii][i];
        int  upper_map = rmap[z][i];
        int lower_map = rmap[z][ii];
        for(int x3k=x3l;x3k<=x3h;++x3k)
          for(int x1i=x1l;x1i<=x1h;++x1i)
            for(int x2i=x2l;x2i<=x2h;++x2i)
              for(int x3i=x3l;x3i<=x3h;++x3i){
                fp_t profile =  current_profile[x1i][x2i][x3i][tmap[z][i][ii]];
                op_pert[1][x3k][x1i][x2i][x3i] += profile * (level_responses[1][(x3i-x3l) * nmap+lower_map+1][x3k] * Blu - level_responses[1][(x3i-x3l)*nmap+upper_map+1][x3k] * Bul) 
                  * h * c / lambda / 4.0 / pi;

                op_pert[2][x3k][x1i][x2i][x3i] += profile * (level_responses[2][(x3i-x3l) * nmap+lower_map+1][x3k] * Blu - level_responses[2][(x3i-x3l)*nmap+upper_map+1][x3k] * Bul) 
                  * h * c / lambda / 4.0 / pi;

                op_pert[3][x3k][x1i][x2i][x3i] += profile * (level_responses[3][(x3i-x3l) * nmap+lower_map+1][x3k] * Blu - level_responses[3][(x3i-x3l)*nmap+upper_map+1][x3k] * Bul) 
                  * h * c / lambda / 4.0 / pi;


                if (x3i == x3k){
                  fp_t derivative = voigt_der_T(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  op_pert[1][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][ii] * (Blu) * h * c / lambda / 4.0 / pi * 
                    (1.0  - pop[x1i][x2i][x3i].n[z][i] / pop[x1i][x2i][x3i].n[z][ii] * Bul / Blu);

                  derivative = voigt_der_density(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  op_pert[2][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][ii] * (Blu) * h * c / lambda / 4.0 / pi * 
                    (1.0  - pop[x1i][x2i][x3i].n[z][i] / pop[x1i][x2i][x3i].n[z][ii] * Bul / Blu);

                  derivative = voigt_der_vt(x1i,x2i,x3i,z,i,ii,lambda, Vlos);
                  op_pert[3][x3k][x1i][x2i][x3i] += derivative * pop[x1i][x2i][x3i].n[z][ii] * (Blu) * h * c / lambda / 4.0 / pi * 
                    (1.0  - pop[x1i][x2i][x3i].n[z][i] / pop[x1i][x2i][x3i].n[z][ii] * Bul / Blu);
                }
            }
      }
  return op_pert;
}

fp_t atom::compute_x_scalar(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos){
  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);
  fp_t x = (lam-lambda*(1.0+vlos[x1i][x2i][x3i]/c))/dld;
  return x;
}

fp_t atom::compute_a_scalar(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos){
  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);
  fp_t dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  fp_t damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  fp_t a = damping / dld * (lam*lam)/(4.0*pi*c);
  return a;
}

int atom::compute_xa_der_scalar(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos, fp_t&x, fp_t &a, 
  fp_t * der_x, fp_t * der_a, fp_t * der_dld){

  // Memset derivatives to zero:
  memset(der_x+1,0,7*sizeof(fp_t));
  memset(der_a+1,0,7*sizeof(fp_t));
  memset(der_dld+1,0,7*sizeof(fp_t));
  
  // First x:
  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);
  x = (lam-lambda*(1.0+vlos[x1i][x2i][x3i]/c))/dld;
  
  // Then a:
  fp_t dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  fp_t damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  a = damping / dld * (lam*lam)/(4.0*pi*c);

  // Then derivative of x:
  fp_t dvd = sqrt(2.0 * k * Temp / mass + vt*vt);
  // Derivative of the doppler with:
  der_dld[1] = lam/(c*dvd) * k / mass;
  der_x[1] = -x / dld * der_dld[1];
  // Density does not influence
  der_x[2] = 0.0;
  // Microturbulent velocity
  der_dld[3] = lam/(c*dvd) * vt;
  der_x[3] = -x / dld * der_dld[3];
  // Systematic velocity
  der_x[4] = -lambda*(1.0/c)/dld;

  // Then derivative of a:
  // First the temperature
  // Derivative of the doppler with:
  der_a[1] = -a/dld * der_dld[1];
  der_a[1] += damp_col_der_T(x1i, x2i, x3i, z, i, ii, Temp, ne, lam)  / dld * (lam*lam)/(4.0*pi*c) * turn_on_damping;

  // Then the density:
  fp_t Nt = fetch_Nt(x1i,x2i,x3i);
  fp_t delta_Nt = delta_Nt_frac * Nt;
  fp_t delta_damp = dampcol * delta_Nt_frac;
  fp_t delta_a = delta_damp / dld * (lam*lam)/(4.0*pi*c);
  der_a[2] = delta_a / delta_Nt * turn_on_damping;

  // Then microturbulent:
  dvd = dld / lam * c;
  // Derivative of a is also non-zero because a depends on dld:
  der_a[3] = -a/dld * der_dld[3];

  return 0;
}
 
fp_t atom::voigt_der_T(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos){

  // Returns explicit derivative of voigt profile with respect to the temperature (so everything except dne/dT)
  // Actually for the moment, we will account for ne_derivative. This is to be debated

  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);
  fp_t x = (lam-lambda*(1.0+vlos[x1i][x2i][x3i]/c))/dld;
  
  // Compute the derivative :
  fp_t dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  fp_t damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  fp_t a = damping / dld * (lam*lam)/(4.0*pi*c);
  fp_t dvd = dld / lam * c;
  
  // Derivative of the doppler with:
  fp_t d_dld = lam/(c*dvd) * k / mass;
  // Derivative of x:
  fp_t dx_dT = -x / dld * d_dld;
  
  // Derivative of a is also non-zero because a depends on dld:
  fp_t da_dT = -a/dld * d_dld;
 
  da_dT += damp_col_der_T(x1i, x2i, x3i, z, i, ii, Temp, ne, lam)  / dld * (lam*lam)/(4.0*pi*c) * turn_on_damping;

  // Now we need to extract the derivatives:
  fp_t H,F,dH,dF;
  fvoigt(x,a,H,F,dH,dF);

  fp_t dH_dT = dH/dld * dx_dT - dF/dld * da_dT;

  fp_t derivative = dH_dT - H/dld/dld*d_dld;

  derivative = (fvoigt(x+dx_dT*0.5*delta_T,a+da_dT*0.5*delta_T) / (dld+0.5*d_dld*delta_T) - fvoigt(x-dx_dT*0.5*delta_T,a-da_dT*0.5*delta_T) / (dld-0.5*d_dld*delta_T))/delta_T;
  return derivative;
}
fp_t atom::voigt_der_density(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos){

  // Returns explicit derivative of voigt profile with respect to the temperature (so everything except dne/dT)
  // Actually for the moment, we will account for ne_derivative. This is to be debated

  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);
  fp_t x = (lam-lambda*(1.0+vlos[x1i][x2i][x3i]/c))/dld;
  // Compute the derivative :
  fp_t dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  fp_t damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  fp_t a = damping / dld * (lam*lam)/(4.0*pi*c);
  fp_t dvd = dld / lam * c;
  
  // Derivative of the doppler with:
  fp_t d_dld = lam/(c*dvd) * k / mass;
  // Derivative of x:
  fp_t dx = 0.0;
  // Derivative of a is also non-zero because a depends on dld:
  fp_t Nt = fetch_Nt(x1i,x2i,x3i);
  fp_t delta_Nt = delta_Nt_frac * Nt;
  fp_t delta_damp = dampcol * delta_Nt_frac * turn_on_damping;
  fp_t delta_a = delta_damp / dld * (lam*lam)/(4.0*pi*c);
  fp_t derivative = (fvoigt(x,a+delta_a*0.5) - fvoigt(x,a-delta_a*0.5)) / dld / delta_Nt;

  //return 0;
  
  // Numerical computation:
  /*
  fp_t Nt = fetch_Nt(x1i,x2i,x3i);
  fp_t delta_Nt = 1E-3 * Nt;
  parent_atm->set_Nt(x1i,x2i,x3i,Nt + delta_Nt/2.0);
  parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
  dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  a = damping / dld * (lam*lam)/(4.0*pi*c);
  fp_t upper = fvoigt(x,a)/dld;
  parent_atm->set_Nt(x1i,x2i,x3i,Nt - delta_Nt/2.0);
  parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
  dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  a = damping / dld * (lam*lam)/(4.0*pi*c);
  fp_t lower = fvoigt(x,a)/dld;
  fp_t derivative_numerical = (upper-lower) / delta_Nt;
  parent_atm->set_Nt(x1i,x2i,x3i,Nt);
  parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
  */
  
  //return derivative_numerical;
  return derivative * turn_on_damping;
}

fp_t atom::voigt_der_vt(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos){

  // Returns explicit derivative of voigt profile with respect to microturbulent velocity
  // This is supposed to be very straightforward

  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);
  fp_t x = (lam-lambda*(1.0+vlos[x1i][x2i][x3i]/c))/dld;
  // Compute the collisions :
  fp_t dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  fp_t damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  fp_t a = damping / dld * (lam*lam)/(4.0*pi*c);
  fp_t dvd = dld / lam * c;
  
  // Derivative of the doppler with:
  fp_t d_dld = lam/(c*dvd) * vt;
  // Derivative of x:
  fp_t dx_dT = -x / dld * d_dld;
  // Derivative of a is also non-zero because a depends on dld:
  fp_t da_dT = -a/dld * d_dld;

  // Now we need to extract the derivatives:
  fp_t H,F,dH,dF;
  fvoigtn(x,a,H,F,dH,dF);

  fp_t dH_dT = dH/dld * dx_dT - dF/dld * da_dT;
  
  fp_t derivative = dH_dT - H/dld/dld*d_dld;
  
  return derivative;
}

fp_t * atom::x_derivative(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos, fp_t **** B_vec, int trans_type, int m){
  
  fp_t * derivative = new fp_t [7]-1;
  memset(derivative+1,0,7*sizeof(fp_t));
  // Get all the necessary stuff and compute ingredients
  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);
  fp_t Bmag = B_vec[1][x1i][x2i][x3i];
  fp_t theta = B_vec[2][x1i][x2i][x3i] * pi /180.0;
  fp_t phi = B_vec[3][x1i][x2i][x3i] * pi/180.0;
  fp_t lambda_B = 4.67E-5 * lam * lam * Bmag; // in cm

  int tr = tmap[z][i][ii]; // Transition mapping 

  fp_t delta_lambda_zeeman = 0.0; // Zeeman shift, from the atomic set-up
  if (trans_type == 0) delta_lambda_zeeman = delta_lambda_p[tr][m];
  if (trans_type == 1) delta_lambda_zeeman = delta_lambda_b[tr][m];
  if (trans_type == 2) delta_lambda_zeeman = delta_lambda_r[tr][m];

  // First compute x itself: 
  fp_t x=(lam+delta_lambda_zeeman*lambda_B-lambda*(1.0+vlos[x1i][x2i][x3i]/c))/dld;

  // Sort out Temperature:
  fp_t dvd = sqrt(2.0 * k * Temp / mass + vt*vt);
  // Derivative of the doppler with:
  fp_t d_dld = lam/(c*dvd) * k / mass;
  derivative[1] = -x / dld * d_dld;
  // Density does not influence
  derivative[2] = 0.0;
  // Microturbulent velocity
  d_dld = lam/(c*dvd) * vt;
  derivative[3] = -x / dld * d_dld;
  // Systematic velocity
  derivative[4] = -lambda*(1.0/c)/dld;
  // Magnetic field magnitude:
  derivative[5] = (delta_lambda_zeeman*4.67E-5 * lam * lam)/dld;
  derivative[6] = derivative[7] = 0.0;

  return derivative;
}

fp_t * atom::a_derivative(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos, fp_t **** B_vec, int trans_type, int m){
  
  fp_t * derivative = new fp_t [7]-1;
  memset(derivative+1,0,7*sizeof(fp_t));
  
  fp_t lam=fabs(h*c/(ee[z][i]-ee[z][ii]));
  fp_t Temp = fetch_temperature(x1i,x2i,x3i);
  fp_t ne = fetch_Ne(x1i,x2i,x3i);
  fp_t vt = fetch_vt(x1i,x2i,x3i);
  fp_t dld = broad_dop(lam, Temp, vt);

  fp_t dampcol = damp_col(x1i, x2i, x3i, z, i, ii, Temp, ne, lam);
  fp_t damping =  damp_rad(A[z],ii,i)+dampcol*turn_on_damping;
  fp_t a = damping / dld * (lam*lam)/(4.0*pi*c);

  // First the temperature
  fp_t dvd = dld / lam * c;
  // Derivative of the doppler with:
  fp_t d_dld = lam/(c*dvd) * k / mass;
  derivative[1] = -a/dld * d_dld;
  derivative[1] += damp_col_der_T(x1i, x2i, x3i, z, i, ii, Temp, ne, lam)  / dld * (lam*lam)/(4.0*pi*c) * turn_on_damping;

  // Then the density:
  fp_t Nt = fetch_Nt(x1i,x2i,x3i);
  fp_t delta_Nt = delta_Nt_frac * Nt;
  fp_t delta_damp = dampcol * delta_Nt_frac;
  fp_t delta_a = delta_damp / dld * (lam*lam)/(4.0*pi*c);
  derivative[2] = delta_a / delta_Nt * turn_on_damping;

  // Then microturbulent:
  dvd = dld / lam * c;
  // Derivative of the doppler with:
  d_dld = lam/(c*dvd) * vt;
  // Derivative of a is also non-zero because a depends on dld:
  derivative[3] = -a/dld * d_dld;
  return derivative;
}

// ------------------------------------------------------------------------------------------------------------------------------------

fp_t atom::broad_dop(fp_t lambda_0,fp_t Temp_in,fp_t vturb_in){
  
  //return sqrt(2.0 * k*6000.0/mass)*lambda_0/c;
  return sqrt(2.0 * k*Temp_in/mass+ vturb_in*vturb_in)*lambda_0/c;
}

void atom::compute_profile_norm(fp_t theta, fp_t phi, fp_t * lambda, fp_t * wlambda, fp_t *** vlos, int nlambda){
  if (ntr){

    // We will also use this to compute the derivative of the profile norm:
    memset(norm_derivative[1][x1l][x2l][x3l]+1,0,7*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
    memset(norm[x1l][x2l][x3l]+1,0, (x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));

    for (int tr=1;tr<=ntr;++tr){    
      // z 
      int z = inverse_tmap[tr][2];
      // higher level
      int i = inverse_tmap[tr][4];
      // lower lvl
      int ii = inverse_tmap[tr][3];

      if (i<nl[z]){

        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t ar=damp_rad(A[z],ii,i);
        fp_t gc = 0.0;
        fp_t Blu=B[z][ii][i],Bul=B[z][i][ii];
        fp_t sf=(lam*lam)/(4.0*pi*c);
              
        for (int x1i=x1l;x1i<=x1h;++x1i)
          for (int x2i=x2l;x2i<=x2h;++x2i)
            for (int x3i=x1l;x3i<=x3h;++x3i){
                          
              // And if higher level is not the continuum, it is time to renormalize the profile
              for (int l=0;l<nlambda;++l){

                fp_t T = fetch_temperature(x1i,x2i,x3i);
                fp_t vturb = fetch_vt(x1i,x2i,x3i);
                fp_t ne = fetch_Ne(x1i,x2i,x3i);
                fp_t dld=broad_dop(lambda[l],T,vturb);
                gc = damp_col(x1i, x2i, x3i, z, i, ii, T, ne, lambda[l]);
                fp_t a=sf*(ar+gc*turn_on_damping)/dld;
                fp_t x=(lam-lambda[l]*(1.0+vlos[x1i][x2i][x3i]/c))/dld;
                fp_t profile = 1.0 / dld * fvoigt(x, a);
                fp_t der =  voigt_der_T(x1i, x2i, x3i, z, i, ii, lambda[l], vlos);
                norm[x1i][x2i][x3i][tr] += profile * wlambda[l];
                norm_derivative[1][x1i][x2i][x3i][tr] += der * wlambda[l]; 
                der =  voigt_der_density(x1i, x2i, x3i, z, i, ii, lambda[l], vlos);
                norm_derivative[2][x1i][x2i][x3i][tr] += der * wlambda[l]; 
                der =  voigt_der_vt(x1i, x2i, x3i, z, i, ii, lambda[l], vlos);
                norm_derivative[3][x1i][x2i][x3i][tr] += der * wlambda[l];                
              }
              //printf("%d %d %e %e \n", x3i, tr, norm[x1i][x2i][x3i][tr], norm_derivative[x1i][x2i][x3i][tr]);             
            }
          }
        }
    }
}
          
// ----------------------------------------------------------------------------------------------------------------------------

// VECTOR VARIANTS:

fp_t *****atom::opacity_vector(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t ****B_vec, fp_t theta,fp_t phi,fp_t lambda){

  // Now this will be interesting.

  fp_t ***** op;
  op = ft5dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, 4, 1, 4);

  fp_t *** op_bf = boundfree_op(Vlos,lambda);
  fp_t *** op_ff = freefree_op(T, Ne, Vlos, lambda);
  fp_t *** op_rayleigh = rayleigh_op(lambda);
  
  for (int x1i = x1l; x1i<=x1h; ++x1i)
    for (int x2i = x2l; x2i<=x2h; ++x2i)
      for (int x3i = x3l; x3i<=x3h; ++x3i)
        for (int s = 1; s<=4; ++s)
          for (int sp = 1; sp<=4; ++sp)
            op[x1i][x2i][x3i][s][sp] = (s == sp) ?  (1.0 * op_bf[x1i][x2i][x3i] + 0.0 * op_ff[x1i][x2i][x3i] + 0.0 * op_rayleigh[x1i][x2i][x3i]) : 0.0;

       
  del_ft3dim(op_bf, x1l, x1h, x2l, x2h, x3l, x3h);
  del_ft3dim(op_ff, x1l, x1h, x2l, x2h, x3l, x3h);
  del_ft3dim(op_rayleigh, x1l,x1h,x2l,x2h,x3l,x3h);

  add(boundbound_op_vector(T,Ne,Vlos,Vt,B_vec, lambda), op, x1l, x1h, x2l, x2h, x3l, x3h, 1, 4, 1, 4);
 
  return op;              
}

fp_t ****atom::emissivity_vector(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t ****B_vec, fp_t theta,fp_t phi,fp_t lambda){
  // This one can be postponed for later as it does not really have a function unless you want to deal with scattering polarization.
  // At the moment, we can just call scalar emissivity:

  fp_t **** em = ft4dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, 4);
  memset(em[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  
  fp_t ***em_scalar=freefree_op(T,Ne,Vlos,lambda);
  em_scalar=add(rayleigh_op(lambda),em_scalar,x1l,x1h,x2l,x2h,x3l,x3h);
  
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i)
            em_scalar[x1i][x2i][x3i] *= Planck_f(lambda, T[x1i][x2i][x3i]);
  
  memset(em_scalar[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));  
  
  em_scalar=add(boundfree_em(Vlos,lambda),em_scalar,x1l,x1h,x2l,x2h,x3l,x3h);

  //memset(em_scalar[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  // Al these above were unpolarized, so copy them to the polarized emission

  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i)
          em[x1i][x2i][x3i][1] = 1.0* em_scalar[x1i][x2i][x3i];

  
  
  // And then add b-b which can be polarized:
  em=add(boundbound_em_vector(T,Ne,Vlos,Vt, B_vec, lambda),em,x1l,x1h,x2l,x2h,x3l,x3h,1,4);

  del_ft3dim(em_scalar,x1l,x1h,x2l,x2h,x3l,x3h);
  return em;
}


fp_t ******* atom::opacity_vector_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda)
{

  fp_t ******* op_pert = ft7dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  memset(op_pert[1][x3l][x1l][x2l][x3l][1]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  fp_t ***** op_scalar_pert = boundfree_op_pert(Vlos,lambda);
  

  for (int p=1;p<=7;++p)
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i){
            op_pert[p][x3k][x1i][x2i][x3i][1][1] = op_pert[p][x3k][x1i][x2i][x3i][2][2] = op_pert[p][x3k][x1i][x2i][x3i][3][3] = 
              op_pert[p][x3k][x1i][x2i][x3i][4][4] = 1.0 * op_scalar_pert[p][x3k][x1i][x2i][x3i];
  }
  del_ft5dim(op_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  op_pert = add(boundbound_op_vector_pert(T,Ne,Vlos,Vt,B_vec,lambda),op_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  
  return op_pert;
}


fp_t ****** atom::emissivity_vector_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B_vec, fp_t theta,fp_t phi,fp_t lambda)
{
  
  fp_t ****** em_pert = ft6dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  memset(em_pert[1][x3l][x1l][x2l][x3l]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  
  fp_t ***** em_scalar_pert = boundfree_em_pert(Vlos,lambda);
  
  for (int p=1;p<=7;++p)
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i)
            em_pert[p][x3k][x1i][x2i][x3i][1] = 1.0 * em_scalar_pert[p][x3k][x1i][x2i][x3i];

  del_ft5dim(em_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);

  em_pert = add(boundbound_em_vector_pert(T,Ne,Vlos,Vt,B_vec,lambda),em_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);

  return em_pert;

}

fp_t ***** atom::boundbound_op_vector(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t ****B_vec, fp_t lambda){

 
  fp_t *****op=ft5dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  memset(op[x1l][x2l][x3l][1]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  for(int z=0;z<=Z;++z) //
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t ar=damp_rad(A[z],ii,i);
        fp_t gc = 0.0;
        fp_t Blu=B[z][ii][i],Bul=B[z][i][ii];
        fp_t sf=(lam*lam)/(4.0*pi*c);
        
        for(int x1i=x1l;x1i<=x1h;++x1i)
          for(int x2i=x2l;x2i<=x2h;++x2i)
            for(int x3i=x3l;x3i<=x3h;++x3i){

              
              fp_t Bmag = B_vec[1][x1i][x2i][x3i];
              //printf("%e \n", Bmag);
              fp_t theta = B_vec[2][x1i][x2i][x3i] * pi /180.0;
              fp_t phi = B_vec[3][x1i][x2i][x3i] * pi/180.0;
              if (Bmag < 0.001)  // very small
                theta = phi = 0.0;
              fp_t lambda_B = 4.67E-5 * lam * lam * Bmag; // in cm

              fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);
              gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
              fp_t a=sf*(ar+gc*turn_on_damping)/dld;

              fp_t x_0=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld; // Doppler shifted wavelength...
              fp_t x = 0.0;
              fp_t H_p, H_b, H_r; // Voigt Profiles
              fp_t F_p, F_b, F_r; // Faraday Voigt profiles
              H_p = H_b = H_r = 0.0;
              F_p = F_b = F_r = 0.0;
              int tr = tmap[z][i][ii]; // 

              fp_t constant_factor = pop[x1i][x2i][x3i].n[z][ii] * (Blu) * h * c / lambda / 4.0 / pi * (1.0  - pop[x1i][x2i][x3i].n[z][i] / pop[x1i][x2i][x3i].n[z][ii] * Bul / Blu);

              if (nm[tr][0] == 1 && nm[tr][1] == 0 && nm[tr][2] == 0){ // means J_u == J_l ==0, can only happen if J is not specified. Return scalar opacity

                x=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                fp_t profile_scalar = fvoigt(x,a) / dld;
                fp_t scalar_op = constant_factor * profile_scalar;
                //printf("%d Wow %d -> %d %e! \n",x3i, i,ii,scalar_op);
                for (int s=1;s<=4;++s)
                  op[x1i][x2i][x3i][s][s] += scalar_op;
              }
              else {
             
                // We now need to compute the total profiles:
                // PI component:
                for (int mi=0;mi<nm[tr][0];++mi){
                  x=(lam+delta_lambda_p[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t H_temp, F_temp, Hder,Fder;
                  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                  H_p += H_temp * S_p[tr][mi];
                  F_p += F_temp * S_p[tr][mi];
                }
                // Sigma B component:
                for (int mi=0;mi<nm[tr][1];++mi){
                  x=(lam+delta_lambda_b[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t H_temp, F_temp, Hder, Fder;
                  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                  H_b += H_temp * S_b[tr][mi];
                  F_b += F_temp * S_b[tr][mi];

                }
                // SIGMA R component:
                for (int mi=0;mi<nm[tr][2];++mi){
                  x=(lam+delta_lambda_r[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t H_temp, F_temp, Hder,Fder;
                  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                  H_r += H_temp * S_r[tr][mi];
                  F_r += F_temp * S_r[tr][mi];
                }
                // Normalize with dld
                H_p /= dld; H_b /= dld; H_r /= dld;
                F_p /= dld; F_b /= dld; F_r /= dld;

                // Eta_I
                op[x1i][x2i][x3i][1][1] += 0.5 * (H_p * sin(theta)*sin(theta) + 0.5 * (H_r + H_b) * (1.0+cos(theta)*cos(theta))) * constant_factor;
                op[x1i][x2i][x3i][2][2] = op[x1i][x2i][x3i][3][3] = op[x1i][x2i][x3i][4][4] = op[x1i][x2i][x3i][1][1];

                // Eta_Q 
                op[x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*cos(2.0*phi) * constant_factor;
                op[x1i][x2i][x3i][2][1] = op[x1i][x2i][x3i][1][2];

                // Eta_U
                op[x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*sin(2.0*phi) * constant_factor;
                op[x1i][x2i][x3i][3][1] = op[x1i][x2i][x3i][1][3];

                // Eta_V
                op[x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * cos(theta) * constant_factor;
                op[x1i][x2i][x3i][4][1] = op[x1i][x2i][x3i][1][4];

                // Rho_Q 
                op[x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*sin(theta)*sin(theta)*cos(2.0*phi) * constant_factor;
                op[x1i][x2i][x3i][4][3] = - op[x1i][x2i][x3i][3][4];

                //Rho_U
                op[x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*sin(theta)*sin(theta)*sin(2.0*phi) * constant_factor;
                op[x1i][x2i][x3i][2][4] = - op[x1i][x2i][x3i][4][2];

                //Rho_V
                op[x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*cos(theta) * constant_factor;
                op[x1i][x2i][x3i][3][2] = -op[x1i][x2i][x3i][2][3];
              }
            }
  }
  //printf("op_pol = %e %e %e %e \n", op[x1l][x2l][25][1][1], op[x1l][x2l][25][1][2],op[x1l][x2l][25][1][3],op[x1l][x2l][25][1][4]);
  return op;
}

fp_t **** atom::boundbound_em_vector(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t ****B_vec,  fp_t lambda){

  fp_t ****em=ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  memset(em[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));

  for(int z=0;z<=Z;++z) //
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t ar=damp_rad(A[z],ii,i);
        fp_t gc = 0.0;
        fp_t Aul = A[z][i][ii];
        fp_t sf=(lam*lam)/(4.0*pi*c);
       
       
        for(int x1i=x1l;x1i<=x1h;++x1i)
          for(int x2i=x2l;x2i<=x2h;++x2i)
            for(int x3i=x3l;x3i<=x3h;++x3i){

              fp_t Bmag = B_vec[1][x1i][x2i][x3i];
              fp_t theta = B_vec[2][x1i][x2i][x3i] * pi /180.0;
              fp_t phi = B_vec[3][x1i][x2i][x3i] * pi/180.0;
              if (Bmag < 0.001)  // very small
                theta = phi = 0.0;
              
              fp_t lambda_B = 4.67E-5 * lam * lam * Bmag; // in cm
              
              fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);
              gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
              fp_t a=sf*(ar+gc*turn_on_damping)/dld;

              fp_t x_0=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
              fp_t x = 0.0;
              fp_t H_p, H_b, H_r; // Voigt Profiles
              fp_t F_p, F_b, F_r; // Faraday Voigt profiles
              H_p = H_b = H_r = 0.0;
              F_p = F_b = F_r = 0.0;
              int tr = tmap[z][i][ii]; // 

              fp_t constant_factor = pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lambda / 4.0 / pi;

              if (nm[tr][0] == 1 && nm[tr][1] == 0 && nm[tr][2] == 0){ // means J_u == J_l == 0, can only happen if J is not specified. Return scalar opacity

                x=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                fp_t profile_scalar = fvoigt(x,a) / dld;
                fp_t scalar_em = constant_factor * profile_scalar;
                em[x1i][x2i][x3i][1] += scalar_em;
              }
              else {
             
                // We now need to compute the total profiles:
                // PI component:
                for (int mi=0;mi<nm[tr][0];++mi){
                  x=(lam+delta_lambda_p[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t H_temp, F_temp, Hder,Fder;
                  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                  H_p += H_temp * S_p[tr][mi];
                  F_p += F_temp * S_p[tr][mi];
                }
                // Sigma B component:
                for (int mi=0;mi<nm[tr][1];++mi){
                  x=(lam+delta_lambda_b[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t H_temp, F_temp, Hder, Fder;
                  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                  H_b += H_temp * S_b[tr][mi];
                  F_b += F_temp * S_b[tr][mi];
                }
                // SIGMA R component:
                for (int mi=0;mi<nm[tr][2];++mi){
                  x=(lam+delta_lambda_r[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t H_temp, F_temp, Hder,Fder;
                  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                  H_r += H_temp * S_r[tr][mi];
                  F_r += F_temp * S_r[tr][mi];
                }
                // Normalize with dld
                H_p /= dld; H_b /= dld; H_r /= dld;
                F_p /= dld; F_b /= dld; F_r /= dld;

                // j_I
                em[x1i][x2i][x3i][1] += 0.5 * (H_p * sin(theta)*sin(theta) + 0.5 * (H_r + H_b) * (1.0+cos(theta)*cos(theta))) * constant_factor;              
                // j_Q 
                em[x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*cos(2.0*phi) * constant_factor;
                // j_U
                em[x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*sin(2.0*phi) * constant_factor;
                // j_V
                em[x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * cos(theta) * constant_factor;
              }
            }
  }
  return em;
}

fp_t ******* atom::boundbound_op_vector_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t ****B_vec, fp_t lambda){
  fp_t *******op_pert=ft7dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  
  memset(op_pert[1][x3l][x1l][x2l][x3l][1]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  for(int z=0;z<=Z;++z) //
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        if (A[z][i][ii] > 1E0){
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t ar=damp_rad(A[z],ii,i);
        fp_t gc = 0.0;
        fp_t Blu=B[z][ii][i],Bul=B[z][i][ii];
        fp_t sf=(lam*lam)/(4.0*pi*c);
        
        for(int x3k=x3l;x3k<=x3h;++x3k)
          for(int x1i=x1l;x1i<=x1h;++x1i)
            for(int x2i=x2l;x2i<=x2h;++x2i)
              for(int x3i=x3l;x3i<=x3h;++x3i){

                fp_t Bmag = B_vec[1][x1i][x2i][x3i];
                fp_t theta = B_vec[2][x1i][x2i][x3i] * pi /180.0;
                fp_t phi = B_vec[3][x1i][x2i][x3i] * pi/180.0;
                if (Bmag < 0.001)  // very small
                  theta = phi = 0.0;
              
                fp_t lambda_B = 4.67E-5 * lam * lam * Bmag; // in cm

                fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);
                gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
                fp_t a=sf*(ar+gc*turn_on_damping)/dld;

                fp_t x_0=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld; // Doppler shifted wavelength...
                fp_t x = 0.0;
                fp_t H_p, H_b, H_r; // Voigt Profiles
                fp_t F_p, F_b, F_r; // Faraday Voigt profiles
                H_p = H_b = H_r = 0.0;
                F_p = F_b = F_r = 0.0;

                fp_t * H_p_der, * H_b_der, * H_r_der; // Voigt Profiles, derivatives
                fp_t * F_p_der, * F_b_der, * F_r_der; // Faraday Voigt profiles, derivatives
                H_p_der = new fp_t [7] - 1; H_b_der = new fp_t [7] - 1; H_r_der = new fp_t [7] - 1;
                F_p_der = new fp_t [7] - 1; F_b_der = new fp_t [7] - 1; F_r_der = new fp_t [7] - 1;
                memset(H_p_der+1,0,7*sizeof(fp_t));
                memset(H_b_der+1,0,7*sizeof(fp_t));
                memset(H_r_der+1,0,7*sizeof(fp_t));
                memset(F_p_der+1,0,7*sizeof(fp_t));
                memset(F_b_der+1,0,7*sizeof(fp_t));
                memset(F_r_der+1,0,7*sizeof(fp_t));

                int tr = tmap[z][i][ii]; // 
                int lower_map = rmap[z][ii];
                int upper_map = rmap[z][i];

                fp_t constant_factor = pop[x1i][x2i][x3i].n[z][ii] * (Blu) * h * c / lambda / 4.0 / pi * (1.0  - pop[x1i][x2i][x3i].n[z][i] / pop[x1i][x2i][x3i].n[z][ii] * Bul / Blu);

                if (nm[tr][0] == 1 && nm[tr][1] == 0 && nm[tr][2] == 0){ // means J_u == J_l ==0. Return scalar opacity

                  x=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t * x_der, * a_der, *dld_der;
                  x_der = new fp_t[7]-1;
                  a_der = new fp_t[7]-1;
                  dld_der = new fp_t[7]-1;
                  compute_xa_der_scalar(x1i, x2i, x3i, z, i, ii, lambda, Vlos, x, a, x_der, a_der, dld_der);
                  //x_der[4] = 0.0;
                  //a_der[4] = 0.0;
                  fp_t H,F,dH,dF;
                  fvoigtn(x,a,H,F,dH,dF);
                  for (int p=1;p<=7;++p){
                    fp_t scalar_op_pert = (level_responses[p][(x3i-x3l) * nmap+lower_map+1][x3k] * Blu - level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Bul)*
                        h * c / lambda / 4.0 / pi * H/dld;
                  
                    if (x3k==x3i && p == 1)
                      scalar_op_pert += constant_factor * (dH*x_der[1] - dF*a_der[1] - H/dld*dld_der[1])/dld;
                    else if (x3k==x3i && p == 2)
                      scalar_op_pert += constant_factor * (dH*x_der[2] - dF*a_der[2] - H/dld*dld_der[2])/dld;
                    else if (x3k==x3i && p == 3)
                      scalar_op_pert += constant_factor * (dH*x_der[3] - dF*a_der[3] - H/dld*dld_der[3])/dld;
                    else if (x3k==x3i && p == 4)
                      scalar_op_pert = constant_factor * (dH*x_der[4])/dld;
                    for (int s=1;s<=4;++s)
                      op_pert[p][x3k][x1i][x2i][x3i][s][s] += scalar_op_pert;
                  }
                  delete [](x_der+1); delete [](a_der+1); delete [](dld_der+1);
                }
                else {
                
                  fp_t * der_a;
                  if (x3k==x3i) 
                    der_a= a_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,0,0);

                  // We now need to compute the total profiles:
                  // PI component:
                  for (int mi=0;mi<nm[tr][0];++mi){
                    x=(lam+delta_lambda_p[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_p += H_temp * S_p[tr][mi];
                    F_p += F_temp * S_p[tr][mi];

                    if (x3k==x3i){
                      fp_t * der_x = x_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,0,mi);
                      //der_x[4] *= cos(pi-theta);

                      for (int p=1;p<=7;++p){
                        H_p_der[p] += (Hder * der_x[p] - Fder * der_a[p]) * S_p[tr][mi];
                        F_p_der[p] += (Fder * der_x[p] + Hder * der_a[p]) * S_p[tr][mi];
                      }
                      delete[](der_x+1);
                    }
                  }
                  // Sigma B component:
                  for (int mi=0;mi<nm[tr][1];++mi){
                    x=(lam+delta_lambda_b[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder, Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_b += H_temp * S_b[tr][mi];
                    F_b += F_temp * S_b[tr][mi];

                    if (x3k==x3i){
                      fp_t * der_x = x_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,1,mi);
                      //der_x[4] *= cos(pi-theta);

                      for (int p=1;p<=7;++p){
                        H_b_der[p] += (Hder * der_x[p] - Fder * der_a[p]) * S_b[tr][mi];
                        F_b_der[p] += (Fder * der_x[p] + Hder * der_a[p]) * S_b[tr][mi];
                      }
                      delete[](der_x+1);
                    }
                  }
                  // SIGMA R component:
                  for (int mi=0;mi<nm[tr][2];++mi){
                    x=(lam+delta_lambda_r[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_r += H_temp * S_r[tr][mi];
                    F_r += F_temp * S_r[tr][mi];

                    if (x3k==x3i){
                      fp_t * der_x = x_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,2,mi);
                      //der_x[4] *= cos(pi-theta);

                      for (int p=1;p<=7;++p){
                        H_r_der[p] += (Hder * der_x[p] - Fder * der_a[p]) * S_r[tr][mi];
                        F_r_der[p] += (Fder * der_x[p] + Hder * der_a[p]) * S_r[tr][mi];
                      }
                      delete[](der_x+1);
                    }
                  }
                  if (x3k==x3i) 
                    delete[](der_a+1);
                  // Normalize with dld
                  H_p /= dld; H_b /= dld; H_r /= dld;
                  F_p /= dld; F_b /= dld; F_r /= dld;

                  for (int p=1;p<=7;++p){
                    H_p_der[p] /= dld; H_b_der[p] /= dld; H_r_der[p] /= dld;
                    F_p_der[p] /= dld; F_b_der[p] /= dld; F_r_der[p] /= dld;
                  }
                  fp_t dvd = dld/lam * c;
                  fp_t d_dld_d_T = lam/c/dvd*k/mass;
                  fp_t d_dld_d_vt = lam/c/dvd*Vt[x1i][x2i][x3i];

                  if (x3k==x3i){
                    H_p_der[1] -= d_dld_d_T/dld * H_p; H_b_der[1] -= d_dld_d_T/dld * H_b; H_r_der[1] -= d_dld_d_T/dld * H_r;
                    H_p_der[3] -= d_dld_d_vt/dld * H_p; H_b_der[3] -= d_dld_d_vt/dld * H_b; H_r_der[3] -= d_dld_d_vt/dld * H_r;
                    F_p_der[1] -= d_dld_d_T/dld * F_p; F_b_der[1] -= d_dld_d_T/dld * F_b; F_r_der[1] -= d_dld_d_T/dld * F_r;
                    F_p_der[3] -= d_dld_d_vt/dld * F_p; F_b_der[3] -= d_dld_d_vt/dld * F_b; F_r_der[3] -= d_dld_d_vt/dld * F_r;
                  }

                  for (int p=1;p<=5;++p){

                    fp_t level_perturbation_factor = (level_responses[p][(x3i-x3l) * nmap+lower_map+1][x3k] * Blu - level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Bul)* h * c / lambda / 4.0 / pi;
                    // Eta_I
                    op_pert[p][x3k][x1i][x2i][x3i][1][1] += 0.5 * (H_p * sin(theta)*sin(theta) + 0.5 * (H_r + H_b) * (1.0+cos(theta)*cos(theta))) * level_perturbation_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][1][1] += 0.5 * (H_p_der[p] * sin(theta)*sin(theta) + 0.5 * (H_r_der[p] + H_b_der[p]) * (1.0+cos(theta)*cos(theta))) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][2] = op_pert[p][x3k][x1i][x2i][x3i][3][3] = op_pert[p][x3k][x1i][x2i][x3i][4][4] = op_pert[p][x3k][x1i][x2i][x3i][1][1];

                    // Eta_Q 
                    op_pert[p][x3k][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*cos(2.0*phi) * level_perturbation_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][1][2] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * sin(theta)*sin(theta)*cos(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][1] = op_pert[p][x3k][x1i][x2i][x3i][1][2];

                    // Eta_U
                    op_pert[p][x3k][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*sin(2.0*phi) * level_perturbation_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][1][3] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * sin(theta)*sin(theta)*sin(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][3][1] = op_pert[p][x3k][x1i][x2i][x3i][1][3];

                    // Eta_V
                    op_pert[p][x3k][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * cos(theta) * level_perturbation_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][1][4] += 0.5 * (H_r_der[p] - H_b_der[p]) * cos(theta) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][4][1] = op_pert[p][x3k][x1i][x2i][x3i][1][4];

                    // Rho_Q 
                    op_pert[p][x3k][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*sin(theta)*sin(theta)*cos(2.0*phi) * level_perturbation_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][3][4] += 0.5*(F_p_der[p]-0.5*(F_r_der[p]+F_b_der[p]))*sin(theta)*sin(theta)*cos(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][4][3] = - op_pert[p][x3k][x1i][x2i][x3i][3][4];

                    //Rho_U
                    op_pert[p][x3k][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*sin(theta)*sin(theta)*sin(2.0*phi) * level_perturbation_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][4][2] += 0.5*(F_p_der[p]-0.5*(F_r_der[p]+F_b_der[p]))*sin(theta)*sin(theta)*sin(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][4] = - op_pert[p][x3k][x1i][x2i][x3i][4][2];

                    //Rho_V
                    op_pert[p][x3k][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*cos(theta) * level_perturbation_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][3] += 0.5*(F_r_der[p]-F_b_der[p])*cos(theta) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][3][2] = -op_pert[p][x3k][x1i][x2i][x3i][2][3];
                  }
                  // But then theta and phi derivatives need to be separated because theta and phi are actually entering these equations:
                  // Theta:
                  if (x3k==x3i){
                    int p = 6;
                    op_pert[p][x3k][x1i][x2i][x3i][1][1] += 0.5 * (H_p *2.0*sin(theta)*cos(theta) + 0.5*(H_r+H_b)*(-2.0*cos(theta)*sin(theta))) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][2] = op_pert[p][x3k][x1i][x2i][x3i][3][3] = op_pert[p][x3k][x1i][x2i][x3i][4][4] = op_pert[p][x3k][x1i][x2i][x3i][1][1];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*sin(theta)*cos(theta)*cos(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][1] = op_pert[p][x3k][x1i][x2i][x3i][1][2];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*sin(theta)*cos(theta)*sin(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][3][1] = op_pert[p][x3k][x1i][x2i][x3i][1][3];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * (-sin(theta)) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][4][1] = op_pert[p][x3k][x1i][x2i][x3i][1][4];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*2.0*sin(theta)*cos(theta)*cos(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][4][3] = - op_pert[p][x3k][x1i][x2i][x3i][3][4];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*2.0*sin(theta)*cos(theta)*sin(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][4] = - op_pert[p][x3k][x1i][x2i][x3i][4][2];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*(-sin(theta)) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][3][2] = -op_pert[p][x3k][x1i][x2i][x3i][2][3];
                    // Phi:
                    p = 7; 
                    op_pert[p][x3k][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*(-2.0*sin(2.0*phi)) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][1] = op_pert[p][x3k][x1i][x2i][x3i][1][2];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*2.0*cos(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][3][1] = op_pert[p][x3k][x1i][x2i][x3i][1][3];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*sin(theta)*sin(theta)*(-2.0*sin(2.0*phi)) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][4][3] = - op_pert[p][x3k][x1i][x2i][x3i][3][4];
                    
                    op_pert[p][x3k][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*sin(theta)*sin(theta)*2.0*cos(2.0*phi) * constant_factor;
                    op_pert[p][x3k][x1i][x2i][x3i][2][4] = - op_pert[p][x3k][x1i][x2i][x3i][4][2];
                    
                  }
                }
                delete [](H_p_der+1); delete [](H_b_der+1); delete [](H_r_der+1);
                delete [](F_p_der+1); delete [](F_b_der+1); delete [](F_r_der+1);
              }
            }
  }
  //printf("op_pol = %e %e %e %e \n", op[x1l][x2l][25][1][1], op[x1l][x2l][25][1][2],op[x1l][x2l][25][1][3],op[x1l][x2l][25][1][4]);
  return op_pert;
}
fp_t ****** atom::boundbound_em_vector_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t ****B_vec, fp_t lambda){
  
  // Same as the above except for the emission
  fp_t ******em_pert=ft6dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  memset(em_pert[1][x3l][x1l][x2l][x3l]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));

  for(int z=0;z<=Z;++z) //
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t ar=damp_rad(A[z],ii,i);
        fp_t gc = 0.0;
        fp_t Aul = A[z][i][ii];
        fp_t sf=(lam*lam)/(4.0*pi*c);
        
        for(int x3k=x3l;x3k<=x3h;++x3k)
          for(int x1i=x1l;x1i<=x1h;++x1i)
            for(int x2i=x2l;x2i<=x2h;++x2i)
              for(int x3i=x3l;x3i<=x3h;++x3i){

                fp_t Bmag = B_vec[1][x1i][x2i][x3i];
                fp_t theta = B_vec[2][x1i][x2i][x3i] * pi /180.0;
                fp_t phi = B_vec[3][x1i][x2i][x3i] * pi/180.0;
                if (Bmag < 0.001)  // very small
                  theta = phi = 0.0;
                fp_t lambda_B = 4.67E-5 * lam * lam * Bmag; // in cm

                fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);
                gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
                fp_t a=sf*(ar+gc*turn_on_damping)/dld;

                fp_t x_0=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                fp_t H_p, H_b, H_r; // Voigt Profiles
                fp_t F_p, F_b, F_r; // Faraday Voigt profiles
                H_p = H_b = H_r = 0.0;
                F_p = F_b = F_r = 0.0;

                fp_t * H_p_der, * H_b_der, * H_r_der; // Voigt Profiles, derivatives
                fp_t * F_p_der, * F_b_der, * F_r_der; // Faraday Voigt profiles, derivatives
                H_p_der = new fp_t [7] - 1; H_b_der = new fp_t [7] - 1; H_r_der = new fp_t [7] - 1;
                F_p_der = new fp_t [7] - 1; F_b_der = new fp_t [7] - 1; F_r_der = new fp_t [7] - 1;
                memset(H_p_der+1,0,7*sizeof(fp_t));
                memset(H_b_der+1,0,7*sizeof(fp_t));
                memset(H_r_der+1,0,7*sizeof(fp_t));
                memset(F_p_der+1,0,7*sizeof(fp_t));
                memset(F_b_der+1,0,7*sizeof(fp_t));
                memset(F_r_der+1,0,7*sizeof(fp_t));

                int tr = tmap[z][i][ii]; // 
                int lower_map = rmap[z][ii];
                int upper_map = rmap[z][i];

                fp_t constant_factor = pop[x1i][x2i][x3i].n[z][i] * (Aul) * h * c / lambda / 4.0 / pi;

                if (nm[tr][0] == 1 && nm[tr][1] == 0 && nm[tr][2] == 0){ // means J_u == J_l ==0, can only happen if J is not specified. Return scalar opacity

                  fp_t x=(lam-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  
                  fp_t * x_der, * a_der, *dld_der;
                  x_der = new fp_t[7]-1;
                  a_der = new fp_t[7]-1;
                  dld_der = new fp_t[7]-1;
                  compute_xa_der_scalar(x1i, x2i, x3i, z, i, ii, lambda, Vlos, x, a, x_der, a_der, dld_der);
                  //x_der[4] = 0.0;
                  //a_der[4] = 0.0;
                  fp_t H,F,dH,dF;
                  fvoigtn(x,a,H,F,dH,dF);
                  for (int p=1;p<=7;++p)
                    em_pert[p][x3k][x1i][x2i][x3i][1] += (level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul)*h * c / lambda / 4.0 / pi * H/dld;
                    if (x3k==x3i){
                      em_pert[1][x3k][x1i][x2i][x3i][1] += constant_factor * (dH*x_der[1] - dF*a_der[1] - H/dld*dld_der[1])/dld;
                      em_pert[2][x3k][x1i][x2i][x3i][1] += constant_factor * (dH*x_der[2] - dF*a_der[2] - H/dld*dld_der[2])/dld;
                      em_pert[3][x3k][x1i][x2i][x3i][1] += constant_factor * (dH*x_der[3] - dF*a_der[3] - H/dld*dld_der[3])/dld;
                      em_pert[4][x3k][x1i][x2i][x3i][1] = constant_factor * (dH*x_der[4])/dld;
                    }
                  delete [](x_der+1); delete [](a_der+1); delete [](dld_der+1);
                }

                else {
                
                  fp_t * der_a;
                  if (x3k==x3i) 
                    der_a=a_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,0,0);

                  // We now need to compute the total profiles:
                  // PI component:
                  for (int mi=0;mi<nm[tr][0];++mi){
                    fp_t x=(lam+delta_lambda_p[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_p += H_temp * S_p[tr][mi];
                    F_p += F_temp * S_p[tr][mi];

                    if (x3k==x3i){
                      fp_t * der_x = x_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,0,mi);
                      //der_x[4] *= cos(pi-theta);

                      for (int p=1;p<=7;++p){
                        H_p_der[p] += (Hder * der_x[p] - Fder * der_a[p]) * S_p[tr][mi];
                        F_p_der[p] += (Fder * der_x[p] + Hder * der_a[p]) * S_p[tr][mi];
                      }
                      delete[](der_x+1);
                    }
                  }
                  // Sigma B component:
                  for (int mi=0;mi<nm[tr][1];++mi){
                    fp_t x=(lam+delta_lambda_b[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder, Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_b += H_temp * S_b[tr][mi];
                    F_b += F_temp * S_b[tr][mi];

                    if (x3k==x3i){
                      fp_t * der_x = x_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,1,mi);
                      //der_x[4] *= cos(pi-theta);

                      for (int p=1;p<=7;++p){
                        H_b_der[p] += (Hder * der_x[p] - Fder * der_a[p]) * S_b[tr][mi];
                        F_b_der[p] += (Fder * der_x[p] + Hder * der_a[p]) * S_b[tr][mi];
                      }
                      delete[](der_x+1);
                    }
                  }
                  // SIGMA R component:
                  for (int mi=0;mi<nm[tr][2];++mi){
                    fp_t x=(lam+delta_lambda_r[tr][mi]*lambda_B-lambda*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_r += H_temp * S_r[tr][mi];
                    F_r += F_temp * S_r[tr][mi];

                    if (x3k==x3i){
                      fp_t * der_x = x_derivative(x1i,x2i,x3i,z,i,ii,lambda,Vlos,B_vec,2,mi);
                      //der_x[4] *= cos(pi-theta);

                      for (int p=1;p<7;++p){
                        H_r_der[p] += (Hder * der_x[p] - Fder * der_a[p]) * S_r[tr][mi];
                        F_r_der[p] += (Fder * der_x[p] + Hder * der_a[p]) * S_r[tr][mi];
                      }
                      delete[](der_x+1);
                    }
                  }
                  if (x3k==x3i) 
                    delete[](der_a+1);
                  // Normalize with dld
                  H_p /= dld; H_b /= dld; H_r /= dld;
                  F_p /= dld; F_b /= dld; F_r /= dld;

                  for (int p=1;p<7;++p){
                    H_p_der[p] /= dld; H_b_der[p] /= dld; H_r_der[p] /= dld;
                    F_p_der[p] /= dld; F_b_der[p] /= dld; F_r_der[p] /= dld;
                  }
                  fp_t dvd = dld/lam * c;
                  fp_t d_dld_d_T = lam/c/dvd*k/mass;
                  fp_t d_dld_d_vt = lam/c/dvd*Vt[x1i][x2i][x3i];

                  if (x3k==x3i){
                    H_p_der[1] -= d_dld_d_T/dld * H_p; H_b_der[1] -= d_dld_d_T/dld * H_b; H_r_der[1] -= d_dld_d_T/dld * H_r;
                    H_p_der[3] -= d_dld_d_vt/dld * H_p; H_b_der[3] -= d_dld_d_vt/dld * H_b; H_r_der[3] -= d_dld_d_vt/dld * H_r;
                    F_p_der[1] -= d_dld_d_T/dld * F_p; F_b_der[1] -= d_dld_d_T/dld * F_b; F_r_der[1] -= d_dld_d_T/dld * F_r;
                    F_p_der[3] -= d_dld_d_vt/dld * F_p; F_b_der[3] -= d_dld_d_vt/dld * F_b; F_r_der[3] -= d_dld_d_vt/dld * F_r;
                  }

                  for (int p=1;p<=5;++p){

                    fp_t level_perturbation_factor = level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * h * c / lambda / 4.0 / pi;
                    // Eta_I
                    em_pert[p][x3k][x1i][x2i][x3i][1] += 0.5 * (H_p * sin(theta)*sin(theta) + 0.5 * (H_r + H_b) * (1.0+cos(theta)*cos(theta))) * level_perturbation_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][1] += 0.5 * (H_p_der[p] * sin(theta)*sin(theta) + 0.5 * (H_r_der[p] + H_b_der[p]) * (1.0+cos(theta)*cos(theta))) * constant_factor;
                    
                    // Eta_Q 
                    em_pert[p][x3k][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*cos(2.0*phi) * level_perturbation_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][2] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * sin(theta)*sin(theta)*cos(2.0*phi) * constant_factor;
                    
                    // Eta_U
                    em_pert[p][x3k][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*sin(2.0*phi) * level_perturbation_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][3] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * sin(theta)*sin(theta)*sin(2.0*phi) * constant_factor;
                    
                    // Eta_V
                    em_pert[p][x3k][x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * cos(theta) * level_perturbation_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][4] += 0.5 * (H_r_der[p] - H_b_der[p]) * cos(theta) * constant_factor;
                  }
                  if (x3k==x3i){
                    int p = 6; // Theta
                    em_pert[p][x3k][x1i][x2i][x3i][1] += 0.5 * (H_p *2.0*sin(theta)*cos(theta) + 0.5*(H_r+H_b)*(-2.0*cos(theta)*sin(theta))) * constant_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*sin(theta)*cos(theta)*cos(2.0*phi) * constant_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*sin(theta)*cos(theta)*sin(2.0*phi) * constant_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * (-sin(theta)) * constant_factor;
                    
                    p = 7; // Phi
                    em_pert[p][x3k][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*(-2.0*sin(2.0*phi)) * constant_factor;
                    em_pert[p][x3k][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * sin(theta)*sin(theta)*2.0*cos(2.0*phi) * constant_factor;
                  }
                }
                  
                delete [](H_p_der+1); delete [](H_b_der+1); delete [](H_r_der+1);
                delete [](F_p_der+1); delete [](F_b_der+1); delete [](F_r_der+1);
              }
  }
  //printf("op_pol = %e %e %e %e \n", op[x1l][x2l][25][1][1], op[x1l][x2l][25][1][2],op[x1l][x2l][25][1][3],op[x1l][x2l][25][1][4]);
  return em_pert;
}


// ----------------------------------------------------------------------------------------------------------------------------


// RATES METHODS: 


// ----------------------------------------------------------------------------------------------------------------------------
fp_t atom::R_ij(int z, int from, int to, fp_t JJ){

	return A[z][from][to] + B[z][from][to] * JJ;
}

fp_t atom::R_ij_local_ALO(int z, int from, int to, fp_t JJ, fp_t local_ALO_operator, fp_t from_pop, fp_t to_pop){

	// First compute line source funciton:

	if (to == from)
		return 0;

	int u, l;
	fp_t upper_pop, lower_pop;

	if (to > from){
		u = to;
		l = from;
		upper_pop = to_pop;
		lower_pop = from_pop;
	}
	else {
		u = from;
		l = to;
		upper_pop = from_pop;
		lower_pop = to_pop;
	}
	fp_t line_source_funtion = upper_pop * A[z][u][l] / (lower_pop * B[z][l][u] - upper_pop * B[z][u][l]);
  fp_t total_rate = A[z][from][to] * (1.0 - local_ALO_operator) + B[z][from][to] * (JJ - local_ALO_operator * line_source_funtion);
  //if (total_rate < 0){
    //printf("Radiative rate less than zero? %e %e %e %e \n", total_rate, local_ALO_operator, JJ, line_source_funtion);
    //exit(1);
    //total_rate = 0.0;
  //}
	return total_rate;

}

fp_t atom::C_ij_dummy (int z, int i, int ii, fp_t T){

	if (i > ii)
		return 1e5;
	if (ii > i)
		return 1e5 * exp((ee[z][i]-ee[z][ii]) / k / T) * fp_t(g[z][ii]) / fp_t(g[z][i]);
	return 0;
}

fp_t atom::C_ij(int z, int from, int to, fp_t T, fp_t Ne){

	// This is the default expression for computing collisional rates for given transition. By default we mean van Regemorter formula: 
	//Stel. Atm. 3rd edition and references therein, page 276

	fp_t oscillator_str = osc_str[z][from][to];
	fp_t en_difference = aps(ee[z][from] - ee[z][to]);
	fp_t u_0 = en_difference / k / T;
	fp_t Gamma_collisional;

	if (z > 0) {// There is one for ions:

		fp_t g_ = 0.2; // This is the one for levels of different principal quantum number; we will have to have another check once we include orbital quantum number as well
		fp_t alt = 0.276 * exp(u_0) * (-Exponential_Integral_Ei(-u_0));
		Gamma_collisional = (g_ > alt) ? g_ : alt;
	}
	else // Another one for neutrals
		Gamma_collisional = (u_0 > 14.0) ? 0.066 / sqrt(u_0) * (1.0 + 1.5 / u_0) : 0.276 * exp(u_0) * (-Exponential_Integral_Ei(-u_0)); 

	fp_t q_ij = 5.465E-11 * sqrt(T) * 14.5 * oscillator_str * 4.7468212E-22 / en_difference / en_difference * u_0 * exp(-u_0) * Gamma_collisional;

  //q_ij = 2E-15 * fp_t(g[z][to]) / fp_t(g[z][from]) * sqrt(T) * 1E6;

	q_ij = (from < to) ? q_ij : q_ij * exp(u_0) * fp_t(g[z][to]) / fp_t(g[z][from]);
  //q_ij = (from > to) ? q_ij : q_ij / exp(u_0) / fp_t(g[z][to]) * fp_t(g[z][from]);

	return (from != to) ? q_ij * Ne  : 0.0;
}

fp_t atom::C_ij_H(int z, int from, int to, fp_t T, fp_t N_partner){

	// Whatever the original collisional rate is, we will modify it for the collisions with neutral Hydrogen. 
	// For this we are using Drawin formula: 
	// Stel. Atm. 3rd edition, page 277

	fp_t oscillator_str = osc_str[z][from][to];
	fp_t en_difference = aps(ee[z][from] - ee[z][to]);
	fp_t u_0 = en_difference / k / T;
	
	fp_t a_0 = 5.29E-8;
	fp_t m_H = 1.6605402E-24;
	fp_t red_mass = mass * m_H / (mass + m_H);
	fp_t q_ij = 16.0 * pi * a_0 * a_0 * sqrt(2.0 * k * T / pi / red_mass) * 4.7468212E-22 / en_difference / en_difference * oscillator_str * mass / m_H / m_H * 9.10938291E-28 * 
		exp(-u_0) / (1.0 + 2.0 / u_0); 


	q_ij = (from < to) ? q_ij : q_ij * exp(u_0) * g[z][to] / g[z][from];

	return (from != to)? q_ij * N_partner : 0.0;
}

fp_t atom::R_i_cont(int z, int i, fp_t JJ, fp_t T){

	return JJ;
}

fp_t atom::R_cont_i(int z, int i, fp_t JJ, fp_t T, fp_t n_e){

	fp_t dE=ip[z]-ee[z][i];
  
  //if (i == 1)
    //printf("%e %e %e %e %e \n", fp_t(g[z][i]) / fp_t(g[z+1][0]), exp(dE / k / T), pow(T, -1.5), n_e, saha_const);
	
  return JJ * fp_t(g[z][i]) / fp_t(g[z+1][0]) * exp(dE / k / T) * pow(T, -1.5) * n_e * saha_const;

}
  

fp_t atom::C_i_cont(int z, int i, fp_t T, fp_t n_e){

	// For the moment we are using Seaton's formula.
	// Stellar Atmospheres 3rd edition, page 276

	fp_t dE=ip[z]-ee[z][i]; // Energy difference between the level and the next continuum state.
	fp_t u0 = dE / k / T;

 	fp_t g_i = 0.3;
	if (Z == 2)
		g_i = 0.2;
	else if (Z == 1)
		g_i = 0.1;

	// Now to obtain ionization cross_section for this particular energy. Here we use van Regemonter formula. And it is not completely clear how 
  // is this to be used, but as is it gives comparable results with RH equations.
	fp_t lambda = h * c / dE;
	fp_t sigma = (bf[z][i]) ? bf[z][i]->U(lambda) : -1.0;
  fp_t rate = 1.55E13 / sqrt(T) * g_i * sigma * exp(-u0) / u0 * n_e;

	return rate;

}

fp_t atom::C_i_cont_dummy(int z, int i, fp_t T, fp_t n_e){
  return 1E5;
}
fp_t atom::C_cont_i_dummy(int z, int i, fp_t T, fp_t n_e){
  fp_t dE=ip[z]-ee[z][i];
  fp_t rate = C_i_cont_dummy(z, i, T, n_e) * fp_t(g[z][i]) / fp_t(g[z+1][0]) * exp(dE / k / T) * pow(T, -1.5) * n_e * saha_const;

  return rate;
}

  

fp_t atom::C_cont_i(int z, int i, fp_t T, fp_t n_e){

	// This will follow from detailed balancing:
	fp_t dE=ip[z]-ee[z][i];
	fp_t rate = C_i_cont(z, i, T, n_e) * fp_t(g[z][i]) / fp_t(g[z+1][0]) * exp(dE / k / T) * pow(T, -1.5) * n_e * saha_const;

	return rate;
}

fp_t atom::boundfree_op(uint08_t z,uint16_t l,fp_t Vlos,fp_t lambda,struct pps &p) {};

