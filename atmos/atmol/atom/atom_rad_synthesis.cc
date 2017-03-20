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

// A separate file versions of functions which compute opacity, emissivity, and their perturbation variants. 
// We compute opacity and emissivity separately 
// The file contains following methods: 
// int atom::op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B, fp_t theta,fp_t phi,
//   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector);

int atom::op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** Bmag, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  // Input quantities are usual ones, T, Ne, Vlos, Vt, B are atmospheric quantities
  // theta,phi,lambda and nlambda give us angles and wavelengths for which to compute opacity/emissivity
  // contributions from given atom are added to op_vector and em_vector.

  boundfree_op_em_vector(T,Ne,Vlos,theta,phi,lambda,nlambda,op_vector,em_vector);
  boundbound_op_em_vector(T,Ne,Vlos,Vt,Bmag,theta,phi,lambda,nlambda,op_vector,em_vector);

  return 0;
}


int atom::boundfree_op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  for (int z=0;z<Z;++z) // For each stage that can be ionized (i.e. not the last one)
    for (int i=0;i<nl[z];++i){ // for each level

    fp_t lambda_crit = (bf[z][i]) ? bf[z][i]->getlambda_crit() : 0.0;
    int nlambda_crit = 0;
    for (int l=1;l<=nlambda;++l){
    	if (lambda[l] > lambda_crit)
    		break;
    	++nlambda_crit;
    }
 
 	if (nlambda_crit){
	  for (int x1i=x1l;x1i<=x1h;++x1i)
	    for (int x2i=x2l;x2i<=x2h;++x2i)
	      for (int x3i=x3l;x3i<=x3h;++x3i){
    		fp_t sigma = (bf[z][i]) ? bf[z][i]->U(lambda_crit) : 0.0;
    		// Then compute referent opacity and emissivity:
    		fp_t pop_mod = pop[x1i][x2i][x3i].n[z+1][0] * fetch_Ne(x1i, x2i, x3i) * saha_const * pow(T[x1i][x2i][x3i], -1.5) 
              * fp_t(g[z][i]) / fp_t(g[z+1][0]) * exp((ip[z] - ee[z][i] - h * c / lambda_crit)/k/T[x1i][x2i][x3i]);

    		fp_t op_loc = sigma * (pop[x1i][x2i][x3i].n[z][i] - pop_mod);
    		fp_t em_loc = sigma * pop_mod * (1.0-exp(-h*c/lambda_crit/k/T[x1i][x2i][x2i])) * Planck_f(lambda_crit, T[x1i][x2i][x3i]);
    		for (int l=1;l<=nlambda_crit;++l){
    			op_vector[l][x1i][x2i][x3i][1][1] += op_loc;
    			em_vector[l][x1i][x2i][x3i][1] += em_loc;
    		}
      }
 	}
  } 

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
              	if (nm[tr][0] == 1 + nm[tr][1] == 0 + nm[tr][2] <= 1){
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
                	op_vector[l][x1i][x2i][x3i][1][1] += 0.5 * (H_p * st*st + 0.5 * (H_r + H_b) * (1.0+ct*st)) * op_loc;	
                	op_vector[l][x1i][x2i][x3i][1][2] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*cp * op_loc;
	                op_vector[l][x1i][x2i][x3i][1][3] += 0.5 * (H_p-0.5*(H_r+H_b)) * st*st*sp * op_loc;
	                op_vector[l][x1i][x2i][x3i][1][4] += 0.5 * (H_r - H_b) * ct * op_loc;
	                op_vector[l][x1i][x2i][x3i][3][4] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*cp * op_loc	;
		              op_vector[l][x1i][x2i][x3i][4][2] += 0.5*(F_p-0.5*(F_r+F_b))*st*st*sp * op_loc;
	                op_vector[l][x1i][x2i][x3i][2][3] += 0.5*(F_r-F_b)*st * op_loc;
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

