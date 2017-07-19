#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <math.h>
#include <ctime>

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

// A separate file versions of functions which compute opacity, emissivity, and their perturbation variants. 
// We compute opacity and emissivity separately 
// The file contains following methods: 
//
// 1) int atom::op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B, fp_t theta,fp_t phi,
//   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector);
// 
// 2) int atom::op_em_vector_plus_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** Bmag, fp_t theta,fp_t phi,
//   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector, fp_t ********op_vector_pert, fp_t*******em_vector_pert)


int atom::op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** Bmag, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  // Input quantities are usual ones, T, Ne, Vlos, Vt, B are atmospheric quantities
  // theta,phi,lambda and nlambda give us angles and wavelengths for which to compute opacity/emissivity
  // contributions from given atom are added to op_vector and em_vector.

  //rayleigh_op_em_vector(T,Ne,Vlos,theta,phi,lambda,nlambda,op_vector,em_vector); //not implemented
  //freefree_op_em_vector(T,Ne,Vlos,theta,phi,lambda,nlambda,op_vector,em_vector); //not implemented 
  boundfree_op_em_vector(T,Ne,Vlos,theta,phi,lambda,nlambda,op_vector,em_vector);
  boundbound_op_em_vector(T,Ne,Vlos,Vt,Bmag,theta,phi,lambda,nlambda,op_vector,em_vector);
  
  return 0;
}

int atom::op_em_vector_plus_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** Bmag, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector, fp_t ********op_vector_pert, fp_t*******em_vector_pert){

  // Input quantities are usual ones, T, Ne, Vlos, Vt, B are atmospheric quantities
  // theta,phi,lambda and nlambda give us angles and wavelengths for which to compute opacity/emissivity
  // contributions from given atom are added to op_vector and em_vector.
  boundfree_op_em_vector_plus_pert(T,Ne,Vlos,theta,phi,lambda,nlambda,op_vector,em_vector,op_vector_pert,em_vector_pert);
  
  //clock_t b=clock();
  boundbound_op_em_vector_plus_pert(T,Ne,Vlos,Vt,Bmag,theta,phi,lambda,nlambda,op_vector,em_vector,op_vector_pert,em_vector_pert);
  //clock_t a=clock();
  //printf("Time spent on BB op/em + pert = %f\n",double(a-b)/CLOCKS_PER_SEC);

  return 0;
}

int atom::freefree_op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos, fp_t theta,fp_t phi,
  fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  fp_t lambda_mean = 0.5*(lambda[1]+lambda[nlambda]);
  fp_t *** op = freefree_op(T,Ne,Vlos,lambda_mean);

  for (int l=1;l<=nlambda;++l)
    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          op_vector[l][x1i][x2i][x3i][1][1] += op[x1i][x2i][x3i];
          em_vector[l][x1i][x2i][x3i][1] += op[x1i][x2i][x3i]*Planck_f(lambda_mean,T[x1i][x2i][x3i]);
  }
  del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
}

int atom::rayleigh_op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos, fp_t theta,fp_t phi,
  fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  fp_t lambda_mean = 0.5*(lambda[1]+lambda[nlambda]);
  fp_t *** op = rayleigh_op(lambda_mean);

  for (int l=1;l<=nlambda;++l)
    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          op_vector[l][x1i][x2i][x3i][1][1] += op[x1i][x2i][x3i];
          em_vector[l][x1i][x2i][x3i][1] += op[x1i][x2i][x3i]*Planck_f(lambda_mean,T[x1i][x2i][x3i]);
  }
  del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
}

int atom::boundfree_op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  for (int z=0;z<Z;++z) // For each stage that can be ionized (i.e. not the last one)
    for (int i=0;i<nl[z];++i){ // for each level

    fp_t lambda_crit = (bf[z][i]) ? bf[z][i]->getlambda_crit() : 0.0;
    int nlambda_crit = 0;
    for (int l=1;l<=nlambda;++l){ // find out if there is 'edge'
    	if (lambda[l] > lambda_crit)
    		break;
    	++nlambda_crit;
    }
    fp_t lambda_mean = lambda[nlambda_crit];

    // then compute op/em for referent wavelength, and assign that value to all wavelengths 
    // that can ionize, i.e. to all the ones <= than lambda[nlambda_crit]
 	  if (nlambda_crit){
	    for (int x1i=x1l;x1i<=x1h;++x1i)
	      for (int x2i=x2l;x2i<=x2h;++x2i)
	        for (int x3i=x3l;x3i<=x3h;++x3i){
    		    fp_t sigma = (bf[z][i]) ? bf[z][i]->U(lambda_mean) : 0.0;
    		    
    		    fp_t pop_mod = pop[x1i][x2i][x3i].n[z+1][0] * fetch_Ne(x1i, x2i, x3i) * saha_const * pow(T[x1i][x2i][x3i], -1.5) 
              * fp_t(g[z][i]) / fp_t(g[z+1][0]) * exp((ip[z] - ee[z][i])/k/T[x1i][x2i][x3i]);

    		    fp_t op_loc = sigma * (pop[x1i][x2i][x3i].n[z][i] - pop_mod*exp(-h*c/lambda_mean/k/T[x1i][x2i][x3i]));
            fp_t em_loc = sigma * pop_mod * (1.0-exp(-h*c/lambda_mean/k/T[x1i][x2i][x3i])) * Planck_f(lambda_mean, T[x1i][x2i][x3i]);
    		    for (int l=1;l<=nlambda_crit;++l){
    			    op_vector[l][x1i][x2i][x3i][1][1] += op_loc;
    			    em_vector[l][x1i][x2i][x3i][1] += em_loc;
    		    } // wvl
          } // points in the atmosphere
 	  } // endif
  } // ionization stages, levels
  return 0;
}

int atom::boundfree_op_em_vector_plus_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector, fp_t********op_vector_pert,fp_t*******em_vector_pert){
  
  //memset(op_vector_pert[1][1][x3l][x1l][x2l][x3l][1]+1,0,7*(x3h-x3l+1)*(x3h-x3l+1)*16*sizeof(fp_t));;
  //memset(em_vector_pert[1][1][x3l][x1l][x2l][x3l]+1,0,7*(x3h-x3l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  //memset(op_vector[1][x1l][x2l][x3l][1]+1,0,(x3h-x3l+1)*16*sizeof(fp_t));;
  //memset(em_vector[1][x1l][x2l][x3l]+1,0,(x3h-x3l+1)*4*sizeof(fp_t));
  
  for (int z=0;z<Z;++z) // For each stage that can be ionized, except the last one
    for (int i=0;i<nl[z];++i){ // for each level

    fp_t lambda_crit = (bf[z][i]) ? bf[z][i]->getlambda_crit() : 0.0;
    int nlambda_crit = 0;
    for (int l=1;l<=nlambda;++l){ // find the 'edge' is there is one
    	if (lambda[l] > lambda_crit)
    		break;
    	++nlambda_crit;
    }
    fp_t lambda_mean = lambda[nlambda_crit];

    fp_t ***** op_pert_transition = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
    fp_t ***** em_pert_transition = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
    memset(op_pert_transition[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1));
    memset(em_pert_transition[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1));
 
 	  // We proceed the same way as above, except we need to loop over parameters and depths
    // which correspond to perturbation at some point.
 	  if (nlambda_crit){
 	    for (int x1i=x1l;x1i<=x1h;++x1i)
		    for (int x2i=x2l;x2i<=x2h;++x2i)
		      for (int x3i=x3l;x3i<=x3h;++x3i){

            fp_t sigma = (bf[z][i]) ? bf[z][i]->U(lambda_mean) : 0.0;
            // Then compute referent opacity and emissivity:
            fp_t pop_mod = pop[x1i][x2i][x3i].n[z+1][0] * fetch_Ne(x1i, x2i, x3i) * saha_const * pow(T[x1i][x2i][x3i], -1.5) 
              * fp_t(g[z][i]) / fp_t(g[z+1][0]) * exp((ip[z] - ee[z][i])/k/T[x1i][x2i][x3i]);

            fp_t op_loc = sigma * (pop[x1i][x2i][x3i].n[z][i] - pop_mod*exp(-h*c/lambda_mean/k/T[x1i][x2i][x3i]));
            fp_t em_loc = sigma * pop_mod * (1.0-exp(-h*c/lambda_mean/k/T[x1i][x2i][x3i])) * Planck_f(lambda_mean, T[x1i][x2i][x3i]);
            for (int l=1;l<=nlambda_crit;++l){
              op_vector[l][x1i][x2i][x3i][1][1] += op_loc;
              em_vector[l][x1i][x2i][x3i][1] += em_loc;
            }
            // ====================================================================================
            // Local, i.e. explicit dependence of opacity and emissivity on Temperature and density:
            // We will compute it for each term (f) separately by computing op_loc * f'/f

            fp_t op_loc_pert[3];
            fp_t em_loc_pert[3];
            
            // -- Opacity : -------------------------------
            // We need this temporary variable because there are two contributors, so we cannot simply 
            // scale the whole opacity with this. This is only the second one:
            fp_t pop_mod_factor = pop[x1i][x2i][x3i].n[z+1][0] * fetch_Ne(x1i, x2i, x3i) * saha_const * pow(T[x1i][x2i][x3i], -1.5) 
              * fp_t(g[z][i]) / fp_t(g[z+1][0]) * exp((ip[z] - ee[z][i] - h*c/lambda_mean)/k/T[x1i][x2i][x3i]);
            // Temperature dependence:      
            op_loc_pert[0] =  -pop_mod_factor * parent_atm->get_ne_lte_derivative(1,x1i,x2i,x3i)/fetch_Ne(x1i, x2i, x3i);
            op_loc_pert[0] += -pop_mod_factor * (-1.5)/T[x1i][x2i][x3i];
            op_loc_pert[0] +=  pop_mod_factor * (ip[z] - ee[z][i] - h*c/lambda_mean)/k/T[x1i][x2i][x3i]/T[x1i][x2i][x3i];
            // Density dependence:
            op_loc_pert[1] =  -pop_mod_factor * parent_atm->get_ne_lte_derivative(2,x1i,x2i,x3i)/fetch_Ne(x1i, x2i, x3i);
            
            // Add explicit dependencies:
            op_loc_pert[0] *= sigma;
            op_loc_pert[1] *= sigma;
            for (int p=1;p<=2;++p)
              op_pert_transition[p][x3i][x1i][x2i][x3i] = op_loc_pert[p-1];
              
            // -- Emissivity : ----------------------------
            // Temperature dependence:
            em_loc_pert[0] =  em_loc * parent_atm->get_ne_lte_derivative(1,x1i,x2i,x3i)/fetch_Ne(x1i, x2i, x3i);
            em_loc_pert[0] += em_loc * (-1.5)/T[x1i][x2i][x3i];
            em_loc_pert[0] += em_loc * -(ip[z] - ee[z][i])/k/T[x1i][x2i][x3i]/T[x1i][x2i][x3i];
            em_loc_pert[0] += em_loc * (-exp(-h*c/lambda_mean/k/T[x1i][x2i][x3i])*h*c/lambda_mean/k/T[x1i][x2i][x3i]/T[x1i][x2i][x3i]) / 
              (1.0-exp(-h*c/lambda_mean/k/T[x1i][x2i][x3i]));
            em_loc_pert[0] += em_loc * Planck_f_derivative(lambda_mean,T[x1i][x2i][x3i]) / Planck_f(lambda_mean,T[x1i][x2i][x3i]);
            // Density dependence:
            em_loc_pert[1] =  em_loc * parent_atm->get_ne_lte_derivative(2,x1i,x2i,x3i)/fetch_Ne(x1i, x2i, x3i);
            // Add explicit dependencies:
            for (int p=1;p<=2;++p)
              em_pert_transition[p][x3i][x1i][x2i][x3i] = em_loc_pert[p-1];
            
            // ====================================================================================
            // Implicit dependencies. I.e. factors which depend on the level perturbations:
            int x3l_p=x3l, x3h_p=x3h; // limits for computing perturbations
            if (!NLTE){
              x3l_p=x3h_p=x3i;
            }

            for (int x3k=x3l_p;x3k<=x3h_p;++x3k){
              op_loc_pert[0] =  -sigma * pop_mod_factor * level_responses[1][(x3i-x3l)*nmap+rmap[z+1][0]+1][x3k]/pop[x1i][x2i][x3i].n[z+1][0];
              op_loc_pert[0] +=  sigma * level_responses[1][(x3i-x3l)*nmap+rmap[z][i]+1][x3k];

              op_loc_pert[1] =  -sigma * pop_mod_factor * level_responses[2][(x3i-x3l)*nmap+rmap[z+1][0]+1][x3k]/pop[x1i][x2i][x3i].n[z+1][0];
              op_loc_pert[1] +=  sigma * level_responses[2][(x3i-x3l)*nmap+rmap[z][i]+1][x3k];

              op_loc_pert[2] =  -sigma * pop_mod_factor * level_responses[3][(x3i-x3l)*nmap+rmap[z+1][0]+1][x3k]/pop[x1i][x2i][x3i].n[z+1][0];
              op_loc_pert[2] +=  sigma * level_responses[3][(x3i-x3l)*nmap+rmap[z][i]+1][x3k];

              em_loc_pert[0] =   em_loc * level_responses[1][(x3i-x3l)*nmap+rmap[z+1][0]+1][x3k]/pop[x1i][x2i][x3i].n[z+1][0];
              em_loc_pert[1] =   em_loc * level_responses[2][(x3i-x3l)*nmap+rmap[z+1][0]+1][x3k]/pop[x1i][x2i][x3i].n[z+1][0];
              em_loc_pert[2] =   em_loc * level_responses[3][(x3i-x3l)*nmap+rmap[z+1][0]+1][x3k]/pop[x1i][x2i][x3i].n[z+1][0];

              for (int p=1;p<=3;++p){
                op_pert_transition[p][x3k][x1i][x2i][x3i] += op_loc_pert[p-1];
                em_pert_transition[p][x3k][x1i][x2i][x3i] += em_loc_pert[p-1];
              } // parameter
            } // depths of perturbation
            // ====================================================================================
      } // points in the atmosphere
 		} // endif
    for (int l=1;l<=nlambda_crit;++l)
      for (int p=1;p<=3;++p)
        for (int x3k=x3l;x3k<=x3h;++x3k)
          for (int x1i=x1l;x1i<=x1h;++x1i)
            for (int x2i=x2l;x2i<=x2h;++x2i)
              for (int x3i=x3l;x3i<=x3h;++x3i){
                op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][1] += op_pert_transition[p][x3k][x1i][x2i][x3i];
                em_vector_pert[l][p][x3k][x1i][x2i][x3i][1]    += em_pert_transition[p][x3k][x1i][x2i][x3i];
    }
  } // ionization stages & levels
  return 0;
}

int atom::boundbound_op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B_vec, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  for (int z=0;z<=Z;++z) // All ionization stages
    for (int i=1;i<nl[z];++i) // upper level
      for (int ii=0;ii<i;++ii){ // lower level
      	if (A[z][i][ii] > 1E1){ // If transition is important at all

      	  fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
          fp_t ar=damp_rad(A[z],ii,i);
          fp_t gc = 0.0;
          fp_t Blu=B[z][ii][i],Bul=B[z][i][ii]; fp_t Aul = A[z][i][ii];
          fp_t sf=(lam*lam)/(4.0*pi*c);
          int tr = tmap[z][i][ii];
          // Then spatially dependent stuff
          for (int x1i=x1l;x1i<=x1h;++x1i)
	    			for (int x2i=x2l;x2i<=x2h;++x2i)
	      	  	for (int x3i=x3l;x3i<=x3h;++x3i){
          		
          		fp_t Bmag = B_vec[1][x1i][x2i][x3i];
            	fp_t theta_B = B_vec[2][x1i][x2i][x3i] * pi /180.0;
            	fp_t phi_B = B_vec[3][x1i][x2i][x3i] * pi/180.0;
            	if (Bmag < 0.001)  // very small
                theta_B = phi_B = 0.0;
            	fp_t lambda_B = 4.67E-5 * lam * lam * Bmag; // in cm

            	fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);
            	gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
            	fp_t a=sf*(ar+gc*turn_on_damping)/dld;

            	fp_t op_loc = (pop[x1i][x2i][x3i].n[z][ii]*Blu - pop[x1i][x2i][x3i].n[z][i]*Bul)*h*c/lam/4.0/pi;
            	fp_t em_loc = pop[x1i][x2i][x3i].n[z][i]*Aul*h*c/lam/4.0/pi;

            	fp_t sp,st,cp,ct;
            	sp = sin(2.0*phi_B); cp = cos(2.0*phi_B);
            	st = sin(theta_B); ct = cos(theta_B); 

              // Check if the line is sensitive to zeeman effect
            	if (split[tr] == 0){ // If not sensitive:
            	  for (int l=1;l<=nlambda;++l){
            	  	fp_t x=(lam-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
            	  	fp_t profile = fvoigt(x,a) / dld;
            	  	op_vector[l][x1i][x2i][x3i][1][1] += op_loc * profile;
            	  	em_vector[l][x1i][x2i][x3i][1] += em_loc * profile;
            	  }
            	}
              else {
		          	for (int l=1;l<=nlambda;++l){
			          	fp_t x = 0.0;
			          	fp_t H_p, H_b, H_r; // Voigt Profiles
			          	fp_t F_p, F_b, F_r; // Faraday Voigt profiles
			          	H_p = H_b = H_r = 0.0; // Initial values
			          	F_p = F_b = F_r = 0.0;
			          	// Then add all the components together:
			          	// PI component:
			          	for (int mi=0;mi<nm[tr][0];++mi){
                  	x=(lam+delta_lambda_p[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                    H_p += H_temp * S_p[tr][mi];
                    F_p += F_temp * S_p[tr][mi];
                  }
                  // Sigma B component:
                	for (int mi=0;mi<nm[tr][1];++mi){
                	  x=(lam+delta_lambda_b[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                	  fp_t H_temp, F_temp, Hder, Fder;
                	  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                	  H_b += H_temp * S_b[tr][mi];
                	  F_b += F_temp * S_b[tr][mi];
                	}
                	// Sigma R component:
                	for (int mi=0;mi<nm[tr][2];++mi){
                	  x=(lam+delta_lambda_r[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                	  fp_t H_temp, F_temp, Hder,Fder;
                	  fvoigt(x,a,H_temp, F_temp, Hder, Fder);
                	  H_r += H_temp * S_r[tr][mi];
                	  F_r += F_temp * S_r[tr][mi];
                	}
                	// Normalize with dld
                	H_p /= dld; H_b /= dld; H_r /= dld;
                	F_p /= dld; F_b /= dld; F_r /= dld;
                	// Add to opacity:
                	op_vector[l][x1i][x2i][x3i][1][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * op_loc;	
                	op_vector[l][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * op_loc;
	                op_vector[l][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * op_loc;
	                op_vector[l][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * ct * op_loc;
	                op_vector[l][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*cp * op_loc	;
		              op_vector[l][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*sp * op_loc;
	                op_vector[l][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*ct * op_loc;
	                // And to emissivity:
	                em_vector[l][x1i][x2i][x3i][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * em_loc;              
                	em_vector[l][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * em_loc;
                	em_vector[l][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * em_loc;
                	em_vector[l][x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * ct * em_loc;
		           }
            	}
              
          }
      	}
  }
  return 0;
}

// ================================================================================================

/*int atom::boundbound_op_em_vector_plus_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B_vec, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector, fp_t******** op_vector_pert, fp_t******* em_vector_pert){

  // This function computes opacity, emissivity as well as responses (derivatives)
  // of both quantities, for the whole wavelength array. Input variables:
  // fp_t *** T,Ne,Vlos,Vt,B_vec - Temperature, electron density, velocities and magnetic field 
  // fp_t theta, phi, *lambda,int nlambda - angle for the ray, and wavelengt array + length
  // fp_t ****** op_vector, ***** em_vector - Stokes opacity and emissivity, dimensions are:
  //             lambda, x1,x2,x3,and then 4x4 for opacity and 4 for emissivity
  // fp_t ******** op_vector_pert, ******* em_vector_pert - Stokes opacity and emissivity derivatives
  //             lambda, parameters, x3k (perturbation), x1,x2,x3,and then again 4x4(or only 4)

  // wavelength-independent quantities
  fp_t *** op_total = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // line integrated opacity
  fp_t *** em_total = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // line integrated emissivity
  fp_t *** a        = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // total damping
  fp_t *** dld      = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // doppler width
  fp_t *** lambda_B = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // Zeeman splitting
  
  for (int z=0;z<=Z;++z) // All ionization stages
    for (int i=1;i<nl[z];++i) // upper level
      for (int ii=0;ii<i;++ii){ // lower level
      	if (A[z][i][ii] > 1E1){ // If transition is important at all

      	  fp_t lam=h*c/(ee[z][i]-ee[z][ii]); // linecenter lambda
          fp_t ar=damp_rad(A[z],ii,i);
          fp_t gc = 0.0;
          fp_t Blu=B[z][ii][i],Bul=B[z][i][ii]; fp_t Aul = A[z][i][ii];
          fp_t sf=(lam*lam)/(4.0*pi*c);
          int tr = tmap[z][i][ii];
          int lower_map = rmap[z][ii];
          int upper_map = rmap[z][i];
          
          // Then spatially dependent quantities:
          for (int x1i=x1l;x1i<=x1h;++x1i)
	    			for (int x2i=x2l;x2i<=x2h;++x2i)
	      	  	for (int x3i=x3l;x3i<=x3h;++x3i){
          		
          		fp_t lambda_B = 4.67E-5 * lam * lam * B_vec[1][x1i][x2i][x3i]; // in cm
            	fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);
            	gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
            	a[x1i][x2i][x3i]=sf*(ar+gc*turn_on_damping)/dld;

            	op_total[x1i][x2i][x3i] = (pop[x1i][x2i][x3i].n[z][ii]*Blu - pop[x1i][x2i][x3i].n[z][i]*Bul)*h*c/lam/4.0/pi;
            	em_total[x1i][x2i][x3i] = pop[x1i][x2i][x3i].n[z][i]*Aul*h*c/lam/4.0/pi;
          }// all atmosphere

          fp_t *x_der, *a_der, *dld_der;

          // Compute the wavelength dependent opacity/emissivity and the derivatives.
          // First compute local derivatives:
          // Check if the line is sensitive to zeeman effect
          if (split[tr]==0){ // If not sensitive
            for (int l=1;l<=nlambda;++l)
              for (int x1i=x1l;x1i<=x1h;++x1i)
                for (int x2i=x2l;x2i<=x2h;++x2i)
                  for (int x3i=x3l;x3i<=x3h;++x3i){
            	  	  fp_t x=(lam-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
            	  	fp_t H,F,dH,dF;
            	  	fvoigtn(x,a,H,F,dH,dF);
            	  	op_vector[l][x1i][x2i][x3i][1][1] += op_loc * H/dld;
            	  	em_vector[l][x1i][x2i][x3i][1] += em_loc * H/dld;

            	  	// And the derivatives of profiles:
            	  	x_der = new fp_t[7]-1;
                	a_der = new fp_t[7]-1;
                	dld_der = new fp_t[7]-1;
                	compute_xa_der_scalar(x1i, x2i, x3i, z, i, ii, lambda[l], Vlos, x, a, x_der, a_der, dld_der);
                	for (int p=1;p<=4;++p){
                	  op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][1] += op_loc * (dH*x_der[p] - dF*a_der[p] - H/dld*dld_der[p])/dld;
                	  em_vector_pert[l][p][x3i][x1i][x2i][x3i][1] += em_loc * (dH*x_der[p] - dF*a_der[p] - H/dld*dld_der[p])/dld;
                	}

                	// But then also depth dependent responses coming from level responses.
                	for (int p=1;p<=5;++p)
                		for (int x3k=x3l;x3k<=x3h;++x3k){
                		op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][1] += (level_responses[p][(x3i-x3l)*nmap+lower_map+1][x3k]*Blu - 
                			level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k]*Bul) * constant_factor * H/dld;
                		em_vector_pert[l][p][x3k][x1i][x2i][x3i][1] += level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k]*Aul*constant_factor * H/dld;
                	}
                	delete[](x_der+1);delete[](a_der+1);delete[](dld_der+1);
              	//}
              }
              else {
            	  // Derivative of a and allocating memory for all the derivatives of the profiles:
            	  a_der = a_derivative(x1i,x2i,x3i,z,i,ii,(lambda[nlambda]+lambda[1])*0.5,Vlos,B_vec,0,0);
            	  fp_t * H_p_der, * H_b_der, * H_r_der; // Voigt Profiles, derivatives
                fp_t * F_p_der, * F_b_der, * F_r_der; // Faraday Voigt profiles, derivatives
                H_p_der = new fp_t [7] - 1; H_b_der = new fp_t [7] - 1; H_r_der = new fp_t [7] - 1;
                F_p_der = new fp_t [7] - 1; F_b_der = new fp_t [7] - 1; F_r_der = new fp_t [7] - 1;
                
		          	//for (int l=1;l<=nlambda;++l){

		          		memset(H_p_der+1,0,7*sizeof(fp_t));
	                memset(H_b_der+1,0,7*sizeof(fp_t));
	                memset(H_r_der+1,0,7*sizeof(fp_t));
	                memset(F_p_der+1,0,7*sizeof(fp_t));
	                memset(F_b_der+1,0,7*sizeof(fp_t));
	                memset(F_r_der+1,0,7*sizeof(fp_t));
				          fp_t x = 0.0;
			          	fp_t H_p, H_b, H_r; // Voigt Profiles
			          	fp_t F_p, F_b, F_r; // Faraday Voigt profiles
			          	H_p = H_b = H_r = 0.0; // Initial values
			          	F_p = F_b = F_r = 0.0;
			          	
			          	// Then add all the components together:
			          	// PI component:
		          		for (int mi=0;mi<nm[tr][0];++mi){
                	  x=(lam+delta_lambda_p[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_p += H_temp * S_p[tr][mi];
                    F_p += F_temp * S_p[tr][mi];

                    // The derivative:
                    x_der = x_derivative(x1i,x2i,x3i,z,i,ii,lambda[l],Vlos,B_vec,0,mi);
                    for (int p=1;p<=7;++p){
                      H_p_der[p] += (Hder * x_der[p] - Fder * a_der[p]) * S_p[tr][mi];
                      F_p_der[p] += (Fder * x_der[p] + Hder * a_der[p]) * S_p[tr][mi];
                    }
                    
                    delete[](x_der+1);
                  }
                  // Sigma B component:
                	for (int mi=0;mi<nm[tr][1];++mi){
                	  x=(lam+delta_lambda_b[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                	  fp_t H_temp, F_temp, Hder, Fder;
                	  fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                	  H_b += H_temp * S_b[tr][mi];
                	  F_b += F_temp * S_b[tr][mi];

                	  // The derivative:
                    x_der = x_derivative(x1i,x2i,x3i,z,i,ii,lambda[l],Vlos,B_vec,1,mi);
                    for (int p=1;p<=7;++p){
                      H_b_der[p] += (Hder * x_der[p] - Fder * a_der[p]) * S_b[tr][mi];
                      F_b_der[p] += (Fder * x_der[p] + Hder * a_der[p]) * S_b[tr][mi];
                    }
                    delete[](x_der+1);
                	}
                	// Sigma R component:
                	for (int mi=0;mi<nm[tr][2];++mi){
                	  x=(lam+delta_lambda_r[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                	  fp_t H_temp, F_temp, Hder,Fder;
                	  fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                	  H_r += H_temp * S_r[tr][mi];
                	  F_r += F_temp * S_r[tr][mi];

                	  // The derivative:
                    x_der = x_derivative(x1i,x2i,x3i,z,i,ii,lambda[l],Vlos,B_vec,2,mi);
                    for (int p=1;p<=7;++p){
                      H_r_der[p] += (Hder * x_der[p] - Fder * a_der[p]) * S_r[tr][mi];
                      F_r_der[p] += (Fder * x_der[p] + Hder * a_der[p]) * S_r[tr][mi];
                    }
                    delete[](x_der+1);
                	}

                	//if (x3i==x3h){
                	//	printf("Atom %s wavelength %e ", id, lambda[l]);
                	//	printf("%e %e %e %e %e %e \n",H_p_der[1],H_b_der[1],H_r_der[1],F_p_der[1],F_b_der[1],F_r_der[1]);
                	//}

                	// Normalize with dld
                	H_p /= dld; H_b /= dld; H_r /= dld;
                	F_p /= dld; F_b /= dld; F_r /= dld;
                	// Add to opacity:
                	op_vector[l][x1i][x2i][x3i][1][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * op_loc;	
                	op_vector[l][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * op_loc;
	                op_vector[l][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * op_loc;
	                op_vector[l][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * ct * op_loc;
	                op_vector[l][x1i][x2i][x3i][3][4] += 0.5 * (F_p-0.5*(F_r+F_b))*st*st*cp * op_loc;
		            	op_vector[l][x1i][x2i][x3i][4][2] += 0.5 * (F_p-0.5*(F_r+F_b))*st*st*sp * op_loc;
	                op_vector[l][x1i][x2i][x3i][2][3] += 0.5 * (F_r-F_b)*ct * op_loc;
	                // And to emissivity:
	                em_vector[l][x1i][x2i][x3i][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * em_loc;              
                	em_vector[l][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * em_loc;
                	em_vector[l][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * em_loc;
                	em_vector[l][x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * ct * em_loc;

                	for (int p=1;p<=5;++p)
                		for (int x3k=x3l;x3k<=x3h;++x3k){
                			fp_t op_loc_pert = (level_responses[p][(x3i-x3l) * nmap+lower_map+1][x3k] * Blu - level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Bul) * constant_factor;
                			fp_t em_loc_pert = level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * constant_factor;

                			op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * op_loc_pert;	
		                	op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * op_loc_pert;
			                op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * op_loc_pert;
			                op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * ct * op_loc_pert;
			                op_vector_pert[l][p][x3k][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*cp * op_loc_pert;
				            	op_vector_pert[l][p][x3k][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*sp * op_loc_pert;
			                op_vector_pert[l][p][x3k][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*ct * op_loc_pert;
			                // And to emissivity:
			                em_vector_pert[l][p][x3k][x1i][x2i][x3i][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * em_loc_pert;              
		                	em_vector_pert[l][p][x3k][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * em_loc_pert;
		                	em_vector_pert[l][p][x3k][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * em_loc_pert;
		                	em_vector_pert[l][p][x3k][x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * ct * em_loc_pert;
                	}

                	// Now the same thing with the derivatives
                	// Normalize with dld:
                	for (int p=1;p<=7;++p){
                    H_p_der[p] /= dld; H_b_der[p] /= dld; H_r_der[p] /= dld;
                    F_p_der[p] /= dld; F_b_der[p] /= dld; F_r_der[p] /= dld;
                  }
                	fp_t dvd = dld/lam * c;
                	fp_t d_dld_d_T = lam/c/dvd*k/mass;
                	fp_t d_dld_d_vt = lam/c/dvd*Vt[x1i][x2i][x3i];
                	H_p_der[1] -= d_dld_d_T/dld * H_p; H_b_der[1] -= d_dld_d_T/dld * H_b; H_r_der[1] -= d_dld_d_T/dld * H_r;
                  H_p_der[3] -= d_dld_d_vt/dld * H_p; H_b_der[3] -= d_dld_d_vt/dld * H_b; H_r_der[3] -= d_dld_d_vt/dld * H_r;
                  F_p_der[1] -= d_dld_d_T/dld * F_p; F_b_der[1] -= d_dld_d_T/dld * F_b; F_r_der[1] -= d_dld_d_T/dld * F_r;
                  F_p_der[3] -= d_dld_d_vt/dld * F_p; F_b_der[3] -= d_dld_d_vt/dld * F_b; F_r_der[3] -= d_dld_d_vt/dld * F_r;
                  // Now these derivatives of profiles go into responses of opacity and emissivity:
                  for (int p=1;p<=5;++p){
                    // Add to opacity:
	                  op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][1] += 0.5 * (H_p_der[p] * st*st + 0.5 * (H_r_der[p] + H_b_der[p]) * (1.0+ct*ct)) * op_loc;	
	                  op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][2] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*cp * op_loc;
			              op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][3] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*sp * op_loc;
			              op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][4] += 0.5 * (H_r_der[p] - H_b_der[p]) * ct * op_loc;
			              op_vector_pert[l][p][x3i][x1i][x2i][x3i][3][4] += 0.5 * (F_p_der[p]-0.5*(F_r_der[p]+F_b_der[p]))*st*st*cp * op_loc;
				          	op_vector_pert[l][p][x3i][x1i][x2i][x3i][4][2] += 0.5 * (F_p_der[p]-0.5*(F_r_der[p]+F_b_der[p]))*st*st*sp * op_loc;
		              	op_vector_pert[l][p][x3i][x1i][x2i][x3i][2][3] += 0.5 * (F_r_der[p]-F_b_der[p])*ct * op_loc;
		              	// And to emissivity:
		              	em_vector_pert[l][p][x3i][x1i][x2i][x3i][1] += 0.5 * (H_p_der[p] * st*st + 0.5 * (H_r_der[p] + H_b_der[p]) * (1.0+ct*ct)) * em_loc;              
	                  em_vector_pert[l][p][x3i][x1i][x2i][x3i][2] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*cp * em_loc;
	                  em_vector_pert[l][p][x3i][x1i][x2i][x3i][3] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*sp * em_loc;
	                  em_vector_pert[l][p][x3i][x1i][x2i][x3i][4] += 0.5 * (H_r_der[p] - H_b_der[p]) * ct * em_loc;
                  }
                  // Then the derivative with respect to theta:
                  // Opacity:
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][1] += 0.5 * (H_p * 2.0*st*ct - 0.5 * (H_r + H_b) * (2.0*ct*st)) * op_loc;	
                	op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*cp * op_loc;
	                op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*sp * op_loc;
	                op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][4] += -0.5 * (H_r - H_b) * st * op_loc;
	                op_vector_pert[l][6][x3i][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*2.0*ct*st*cp * op_loc;
		            	op_vector_pert[l][6][x3i][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*2.0*ct*st*sp * op_loc;
	                op_vector_pert[l][6][x3i][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*ct * op_loc;
	                // And to emissivity:
	               	em_vector_pert[l][6][x3i][x1i][x2i][x3i][1] += 0.5 * (H_p * 2.0*st*ct - 0.5 * (H_r + H_b) * (2.0*ct*st)) * em_loc;              
                	em_vector_pert[l][6][x3i][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*cp * em_loc;
                	em_vector_pert[l][6][x3i][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*sp * em_loc;
                	em_vector_pert[l][6][x3i][x1i][x2i][x3i][4] += -0.5 * (H_r - H_b) * st * em_loc;
                	
                	// And then w.r.t to phi:
                	// Opacity:
                	op_vector_pert[l][7][x3i][x1i][x2i][x3i][1][2] += -0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*sp * op_loc;
	                op_vector_pert[l][7][x3i][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*cp * op_loc;
	                op_vector_pert[l][7][x3i][x1i][x2i][x3i][3][4] += -0.5*(F_p-0.5*(F_r+F_b))*st*st*2.0*sp * op_loc;
		            	op_vector_pert[l][7][x3i][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*2.0*cp * op_loc;
	                // And to emissivity:
	               	em_vector_pert[l][7][x3i][x1i][x2i][x3i][2] += -0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*sp * em_loc;
                	em_vector_pert[l][7][x3i][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*cp * em_loc;        
		          //}
		          delete[](H_p_der+1);delete[](H_r_der+1);delete[](H_b_der+1);
		          delete[](F_p_der+1);delete[](F_r_der+1);delete[](F_b_der+1);
		          delete[](a_der+1);
            }
              
          }
      	}
  }
  return 0;
}*/

int atom::boundbound_op_em_vector_plus_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B_vec, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector, fp_t******** op_vector_pert, fp_t******* em_vector_pert){

  for (int z=0;z<=Z;++z) // All ionization stages
    for (int i=1;i<nl[z];++i) // upper level
      for (int ii=0;ii<i;++ii){ // lower level
        if (A[z][i][ii] > 1E1){ // If transition is important at all

          fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
          fp_t ar=damp_rad(A[z],ii,i);
          fp_t gc = 0.0;
          fp_t Blu=B[z][ii][i],Bul=B[z][i][ii]; fp_t Aul = A[z][i][ii];
          fp_t sf=(lam*lam)/(4.0*pi*c);
          int tr = tmap[z][i][ii];
          int lower_map = rmap[z][ii];
          int upper_map = rmap[z][i];
          // Then spatially dependent stuff
          for (int l=1;l<=nlambda;++l)
          for (int x1i=x1l;x1i<=x1h;++x1i)
            for (int x2i=x2l;x2i<=x2h;++x2i)
              for (int x3i=x3l;x3i<=x3h;++x3i){
              
              fp_t Bmag = B_vec[1][x1i][x2i][x3i];
              fp_t theta_B = B_vec[2][x1i][x2i][x3i] * pi /180.0;
              fp_t phi_B = B_vec[3][x1i][x2i][x3i] * pi/180.0;
              if (Bmag < 0.001)  // very small
                theta_B = phi_B = 0.0;
              fp_t lambda_B = 4.67E-5 * lam * lam * Bmag; // in cm

              fp_t dld=broad_dop(lam,T[x1i][x2i][x3i],Vt[x1i][x2i][x3i]);
              gc = damp_col(x1i, x2i, x3i, z, i, ii, T[x1i][x2i][x3i], Ne[x1i][x2i][x3i], lam);
              fp_t a=sf*(ar+gc*turn_on_damping)/dld;

              fp_t op_loc = (pop[x1i][x2i][x3i].n[z][ii]*Blu - pop[x1i][x2i][x3i].n[z][i]*Bul)*h*c/lam/4.0/pi;
              fp_t em_loc = pop[x1i][x2i][x3i].n[z][i]*Aul*h*c/lam/4.0/pi;

              fp_t constant_factor = h*c/lam/4.0/pi;

              fp_t sp,st,cp,ct;
              sp = sin(2.0*phi_B); cp = cos(2.0*phi_B);
              st = sin(theta_B); ct = cos(theta_B); 

              fp_t *x_der, *a_der, *dld_der;

              // Check if the line is sensitive to zeeman effect
              if (split[tr]==0){ // If not sensitive
                //for (int l=1;l<=nlambda;++l){
                  fp_t x=(lam-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                  fp_t H,F,dH,dF;
                  fvoigtn(x,a,H,F,dH,dF);
                  op_vector[l][x1i][x2i][x3i][1][1] += op_loc * H/dld;
                  em_vector[l][x1i][x2i][x3i][1] += em_loc * H/dld;

                  // And the derivatives of profiles:
                  x_der = new fp_t[7]-1;
                  a_der = new fp_t[7]-1;
                  dld_der = new fp_t[7]-1;
                  compute_xa_der_scalar(x1i, x2i, x3i, z, i, ii, lambda[l], Vlos, x, a, x_der, a_der, dld_der);
                  for (int p=1;p<=4;++p){
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][1] += op_loc * (dH*x_der[p] - dF*a_der[p] - H/dld*dld_der[p])/dld;
                    em_vector_pert[l][p][x3i][x1i][x2i][x3i][1] += em_loc * (dH*x_der[p] - dF*a_der[p] - H/dld*dld_der[p])/dld;
                  }

                  // But then also depth dependent responses coming from level responses.
                  int x3l_p=x3l, x3h_p=x3h; // limits for computing perturbations
                  if (!NLTE){
                    x3l_p=x3h_p=x3i;
                  }
                  for (int p=1;p<=5;++p)
                    for (int x3k=x3l_p;x3k<=x3h_p;++x3k){
                    op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][1] += (level_responses[p][(x3i-x3l)*nmap+lower_map+1][x3k]*Blu - 
                      level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k]*Bul) * constant_factor * H/dld;
                    em_vector_pert[l][p][x3k][x1i][x2i][x3i][1] += level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k]*Aul*constant_factor * H/dld;
                  }
                  delete[](x_der+1);delete[](a_der+1);delete[](dld_der+1);
                //}
              }
              else {
                // Derivative of a and allocating memory for all the derivatives of the profiles:
                a_der = a_derivative(x1i,x2i,x3i,z,i,ii,(lambda[nlambda]+lambda[1])*0.5,Vlos,B_vec,0,0);
                fp_t * H_p_der, * H_b_der, * H_r_der; // Voigt Profiles, derivatives
                fp_t * F_p_der, * F_b_der, * F_r_der; // Faraday Voigt profiles, derivatives
                H_p_der = new fp_t [7] - 1; H_b_der = new fp_t [7] - 1; H_r_der = new fp_t [7] - 1;
                F_p_der = new fp_t [7] - 1; F_b_der = new fp_t [7] - 1; F_r_der = new fp_t [7] - 1;
                
                //for (int l=1;l<=nlambda;++l){

                  memset(H_p_der+1,0,7*sizeof(fp_t));
                  memset(H_b_der+1,0,7*sizeof(fp_t));
                  memset(H_r_der+1,0,7*sizeof(fp_t));
                  memset(F_p_der+1,0,7*sizeof(fp_t));
                  memset(F_b_der+1,0,7*sizeof(fp_t));
                  memset(F_r_der+1,0,7*sizeof(fp_t));
                  fp_t x = 0.0;
                  fp_t H_p, H_b, H_r; // Voigt Profiles
                  fp_t F_p, F_b, F_r; // Faraday Voigt profiles
                  H_p = H_b = H_r = 0.0; // Initial values
                  F_p = F_b = F_r = 0.0;
                  
                  // Then add all the components together:
                  // PI component:
                  for (int mi=0;mi<nm[tr][0];++mi){
                    x=(lam+delta_lambda_p[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_p += H_temp * S_p[tr][mi];
                    F_p += F_temp * S_p[tr][mi];

                    // The derivative:
                    x_der = x_derivative(x1i,x2i,x3i,z,i,ii,lambda[l],Vlos,B_vec,0,mi);
                    for (int p=1;p<=7;++p){
                      H_p_der[p] += (Hder * x_der[p] - Fder * a_der[p]) * S_p[tr][mi];
                      F_p_der[p] += (Fder * x_der[p] + Hder * a_der[p]) * S_p[tr][mi];
                    }
                    
                    delete[](x_der+1);
                  }
                  // Sigma B component:
                  for (int mi=0;mi<nm[tr][1];++mi){
                    x=(lam+delta_lambda_b[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder, Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_b += H_temp * S_b[tr][mi];
                    F_b += F_temp * S_b[tr][mi];

                    // The derivative:
                    x_der = x_derivative(x1i,x2i,x3i,z,i,ii,lambda[l],Vlos,B_vec,1,mi);
                    for (int p=1;p<=7;++p){
                      H_b_der[p] += (Hder * x_der[p] - Fder * a_der[p]) * S_b[tr][mi];
                      F_b_der[p] += (Fder * x_der[p] + Hder * a_der[p]) * S_b[tr][mi];
                    }
                    delete[](x_der+1);
                  }
                  // Sigma R component:
                  for (int mi=0;mi<nm[tr][2];++mi){
                    x=(lam+delta_lambda_r[tr][mi]*lambda_B-lambda[l]*(1.0+Vlos[x1i][x2i][x3i]/c))/dld;
                    fp_t H_temp, F_temp, Hder,Fder;
                    fvoigtn(x,a,H_temp, F_temp, Hder, Fder);
                    H_r += H_temp * S_r[tr][mi];
                    F_r += F_temp * S_r[tr][mi];

                    // The derivative:
                    x_der = x_derivative(x1i,x2i,x3i,z,i,ii,lambda[l],Vlos,B_vec,2,mi);
                    for (int p=1;p<=7;++p){
                      H_r_der[p] += (Hder * x_der[p] - Fder * a_der[p]) * S_r[tr][mi];
                      F_r_der[p] += (Fder * x_der[p] + Hder * a_der[p]) * S_r[tr][mi];
                    }
                    delete[](x_der+1);
                  }

                  // Normalize with dld
                  H_p /= dld; H_b /= dld; H_r /= dld;
                  F_p /= dld; F_b /= dld; F_r /= dld;
                  // Add to opacity:
                  op_vector[l][x1i][x2i][x3i][1][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * op_loc;  
                  op_vector[l][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * op_loc;
                  op_vector[l][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * op_loc;
                  op_vector[l][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * ct * op_loc;
                  op_vector[l][x1i][x2i][x3i][3][4] += 0.5 * (F_p-0.5*(F_r+F_b))*st*st*cp * op_loc;
                  op_vector[l][x1i][x2i][x3i][4][2] += 0.5 * (F_p-0.5*(F_r+F_b))*st*st*sp * op_loc;
                  op_vector[l][x1i][x2i][x3i][2][3] += 0.5 * (F_r-F_b)*ct * op_loc;
                  // And to emissivity:
                  em_vector[l][x1i][x2i][x3i][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * em_loc;              
                  em_vector[l][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * em_loc;
                  em_vector[l][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * em_loc;
                  em_vector[l][x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * ct * em_loc;

                  int x3l_p=x3l, x3h_p=x3h; // limits for computing perturbations
                  if (!NLTE){
                    x3l_p=x3h_p=x3i;
                  }

                  for (int p=1;p<=5;++p)
                    for (int x3k=x3l_p;x3k<=x3h_p;++x3k){
                      fp_t op_loc_pert = (level_responses[p][(x3i-x3l) * nmap+lower_map+1][x3k] * Blu - level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Bul) * constant_factor;
                      fp_t em_loc_pert = level_responses[p][(x3i-x3l)*nmap+upper_map+1][x3k] * Aul * constant_factor;

                      op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * op_loc_pert;  
                      op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * op_loc_pert;
                      op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * op_loc_pert;
                      op_vector_pert[l][p][x3k][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * ct * op_loc_pert;
                      op_vector_pert[l][p][x3k][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*cp * op_loc_pert;
                      op_vector_pert[l][p][x3k][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*sp * op_loc_pert;
                      op_vector_pert[l][p][x3k][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*ct * op_loc_pert;
                      // And to emissivity:
                      em_vector_pert[l][p][x3k][x1i][x2i][x3i][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*ct)) * em_loc_pert;              
                      em_vector_pert[l][p][x3k][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * em_loc_pert;
                      em_vector_pert[l][p][x3k][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * em_loc_pert;
                      em_vector_pert[l][p][x3k][x1i][x2i][x3i][4] += 0.5 * (H_r - H_b) * ct * em_loc_pert;
                  }

                  // Now the same thing with the derivatives
                  // Normalize with dld:
                  for (int p=1;p<=7;++p){
                    H_p_der[p] /= dld; H_b_der[p] /= dld; H_r_der[p] /= dld;
                    F_p_der[p] /= dld; F_b_der[p] /= dld; F_r_der[p] /= dld;
                  }
                  fp_t dvd = dld/lam * c;
                  fp_t d_dld_d_T = lam/c/dvd*k/mass;
                  fp_t d_dld_d_vt = lam/c/dvd*Vt[x1i][x2i][x3i];
                  H_p_der[1] -= d_dld_d_T/dld * H_p; H_b_der[1] -= d_dld_d_T/dld * H_b; H_r_der[1] -= d_dld_d_T/dld * H_r;
                  H_p_der[3] -= d_dld_d_vt/dld * H_p; H_b_der[3] -= d_dld_d_vt/dld * H_b; H_r_der[3] -= d_dld_d_vt/dld * H_r;
                  F_p_der[1] -= d_dld_d_T/dld * F_p; F_b_der[1] -= d_dld_d_T/dld * F_b; F_r_der[1] -= d_dld_d_T/dld * F_r;
                  F_p_der[3] -= d_dld_d_vt/dld * F_p; F_b_der[3] -= d_dld_d_vt/dld * F_b; F_r_der[3] -= d_dld_d_vt/dld * F_r;
                  // Now these derivatives of profiles go into responses of opacity and emissivity:
                  for (int p=1;p<=5;++p){
                    // Add to opacity:
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][1] += 0.5 * (H_p_der[p] * st*st + 0.5 * (H_r_der[p] + H_b_der[p]) * (1.0+ct*ct)) * op_loc;  
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][2] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*cp * op_loc;
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][3] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*sp * op_loc;
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][4] += 0.5 * (H_r_der[p] - H_b_der[p]) * ct * op_loc;
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][3][4] += 0.5 * (F_p_der[p]-0.5*(F_r_der[p]+F_b_der[p]))*st*st*cp * op_loc;
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][4][2] += 0.5 * (F_p_der[p]-0.5*(F_r_der[p]+F_b_der[p]))*st*st*sp * op_loc;
                    op_vector_pert[l][p][x3i][x1i][x2i][x3i][2][3] += 0.5 * (F_r_der[p]-F_b_der[p])*ct * op_loc;
                    // And to emissivity:
                    em_vector_pert[l][p][x3i][x1i][x2i][x3i][1] += 0.5 * (H_p_der[p] * st*st + 0.5 * (H_r_der[p] + H_b_der[p]) * (1.0+ct*ct)) * em_loc;              
                    em_vector_pert[l][p][x3i][x1i][x2i][x3i][2] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*cp * em_loc;
                    em_vector_pert[l][p][x3i][x1i][x2i][x3i][3] += 0.5 * (H_p_der[p]-0.5*(H_r_der[p]+H_b_der[p])) * st*st*sp * em_loc;
                    em_vector_pert[l][p][x3i][x1i][x2i][x3i][4] += 0.5 * (H_r_der[p] - H_b_der[p]) * ct * em_loc;
                  }
                  // Then the derivative with respect to theta:
                  // Opacity:
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][1] += 0.5 * (H_p * 2.0*st*ct - 0.5 * (H_r + H_b) * (2.0*ct*st)) * op_loc; 
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*cp * op_loc;
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*sp * op_loc;
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][1][4] += -0.5 * (H_r - H_b) * st * op_loc;
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*2.0*ct*st*cp * op_loc;
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*2.0*ct*st*sp * op_loc;
                  op_vector_pert[l][6][x3i][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*ct * op_loc;
                  // And to emissivity:
                  em_vector_pert[l][6][x3i][x1i][x2i][x3i][1] += 0.5 * (H_p * 2.0*st*ct - 0.5 * (H_r + H_b) * (2.0*ct*st)) * em_loc;              
                  em_vector_pert[l][6][x3i][x1i][x2i][x3i][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*cp * em_loc;
                  em_vector_pert[l][6][x3i][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * 2.0*st*ct*sp * em_loc;
                  em_vector_pert[l][6][x3i][x1i][x2i][x3i][4] += -0.5 * (H_r - H_b) * st * em_loc;
                  
                  // And then w.r.t to phi:
                  // Opacity:
                  op_vector_pert[l][7][x3i][x1i][x2i][x3i][1][2] += -0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*sp * op_loc;
                  op_vector_pert[l][7][x3i][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*cp * op_loc;
                  op_vector_pert[l][7][x3i][x1i][x2i][x3i][3][4] += -0.5*(F_p-0.5*(F_r+F_b))*st*st*2.0*sp * op_loc;
                  op_vector_pert[l][7][x3i][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*2.0*cp * op_loc;
                  // And to emissivity:
                  em_vector_pert[l][7][x3i][x1i][x2i][x3i][2] += -0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*sp * em_loc;
                  em_vector_pert[l][7][x3i][x1i][x2i][x3i][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*2.0*cp * em_loc;        
              //}
              delete[](H_p_der+1);delete[](H_r_der+1);delete[](H_b_der+1);
              delete[](F_p_der+1);delete[](F_r_der+1);delete[](F_b_der+1);
              delete[](a_der+1);
            }
              
          }
        }
  }
  return 0;
}


