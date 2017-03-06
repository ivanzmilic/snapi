#include <math.h>
#include <stdlib.h>
#include <string.h>

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


fp_t atom::pops_broyden(atmol **atm,uint16_t natm,fp_t Temp,fp_t ne,int32_t ix1,int32_t ix2,int32_t ix3)

// Same as above, except it is done using Broyden matrix.
{

  int x3i_control = x3h;

  if(Jb){
    fp_t *J=Jb[ix1][ix2][ix3]; // integrated intensity
    fp_t *L=Ls[ix1][ix2][ix3]; // approximate lambda operator
    fp_t *nrm = norm[ix1][ix2][ix3]; // norm of the profile: <---- Who thought that THIS would be thing you need to save? Frequency-by-frequency approach is a bitch

    fp_t **M=ft2dim(1,nmap,1,nmap); // Allocate the rate matrix
    memset(M[1]+1,0,nmap*nmap*sizeof(fp_t)); // Set all elements to zero

    // Setup the rate equation:
    for(int i=0;i<nmap-1;++i){ // For all the 'levels' but the last one 
      uint16_t l=lmap[i]; // "current lvl
      uint08_t z=zmap[i]; // apropriate ionization stage
      
      // First account for b-b transitions within the particular ionization stage:
      for(int ll=0; ll<nl[z]; ++ll){ // For all the levels:
        
        fp_t JJ = (tmap[z][l][ll])?J[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] :-1.0;   // angular and frequency integrated intensity
        fp_t LL = (tmap[z][l][ll])?L[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] : 0.0;

       /*  if (ll != l && JJ > 0){
          fp_t lambda = h*c/aps(ee[z][l]-ee[z][ll]); 
          JJ = 2.d * h * c * c / lambda/ lambda / lambda / lambda/ lambda / (exp(1.439 / lambda / Temp) - 1.0);
        }*/
        
        /*if (ix3 == x3i_control && ll != l){
          fp_t lambda = h*c/aps(ee[z][l]-ee[z][ll]); 
          printf("map = %d z = %d l = %d ll = %d J = %e B = %e L = %e \n",tmap[z][l][ll], z, l, ll, JJ, 
            (2.d * h * c * c / lambda/ lambda / lambda / lambda/ lambda / (exp(1.439 / lambda / Temp) - 1.0)), LL);
        }*/

       
        // Transitions from this level

        //fp_t Radiative_rates = R_ij(z, l, ll, JJ);//A[z][l][ll] + B[z][l][ll] * JJ;
        fp_t Radiative_rates = R_ij_local_ALO(z, l, ll, JJ, LL, pop[ix1][ix2][ix3].n[z][l], pop[ix1][ix2][ix3].n[z][ll]);
        fp_t Collisional_rates = C_ij_dummy(z, l, ll, Temp);
        M[i+1][i+1] -= (Radiative_rates + Collisional_rates); 

        // Transitions to this level:
        //Radiative_rates = R_ij(z, ll, l, JJ);//A[z][ll][l] + B[z][ll][l] * JJ;
        Radiative_rates = R_ij_local_ALO(z, ll, l, JJ, LL, pop[ix1][ix2][ix3].n[z][ll], pop[ix1][ix2][ix3].n[z][l]);
        Collisional_rates = C_ij_dummy(z, ll, l, Temp);
        int dl = ll - l;
        M[i+1][i+1+dl] += (Radiative_rates + Collisional_rates);  
      }
      
      // Then account for b-f and f-b transitions:
      if(int dl=nl[z]-l){              // ground level of next ionization stage (if not this level)
        // ionization out of level z,l
        fp_t JJ=(tmap[z][l][nl[z]])?J[tmap[z][l][nl[z]]]:-1.0;     // angular and frequency integrated intensity
        /*
        // Transitions from this level:
        fp_t Radiative_rates = R_i_cont(z, l, JJ);
        fp_t Collisional_rates = C_i_cont(z, l, Temp);

        M[i+1][i+1] -= (Radiative_rates + Collisional_rates);

        Radiative_rates = R_cont_i (z, l, JJ);
        Collisional_rates = C_cont_i (z, l, Temp);

        M[i+1][nl[z]+1] += (Radiative_rates + Collisional_rates);*/
      }
        
        /*
        fp_t Rbf=rbf(z,l,JJ,Temp,ne,ix1,ix2,ix3);
        fp_t Cbf=cbf(z,l,atm,natm,Temp,ne,ix1,ix2,ix3);      // colissional ionization
        M[i+1][i+1]-=(Rbf+Cbf);        // b-f ionization from level l
        M[i+dl+1][i+1]+=(Rbf+Cbf);     // b-f rate contributes to next stage ground level
        // recombination to level z,l
        fp_t Rfb=rfb(z,l,JJ,Temp,ne,ix1,ix2,ix3);   // spontaneous recombination?
        fp_t Cfb=cfb(z,l,Cbf,Temp,ne,ix1,ix2,ix3);  // colissional recombination [negligible?]
        M[i+1][i+dl+1]+=(Rfb+Cfb);     // f-b recombination to level l
        M[i+dl+1][i+dl+1]-=(Rfb+Cfb);  // remove f-b from next stage ground level**/ 
    }
   for(int ii=1;ii<=nmap;++ii) M[nmap][ii]=1.0;  // final equation: sum(n)=cst (number of particles remains the same)


  // ------------------------------------------------------------------------------------------------------------------------------------------------
  // Broyden solver:
  
  fp_t ** Broyden_local = Broyden[ix1][ix2][ix3];

  fp_t * F = compute_F_broyden(ix1, ix2, ix3, Temp, J, nrm);

  // If Broyden matrix is not set_up, set it up:
  if (!Broyden_local[1][1]){


    //printf("about to initialize broyden matrix....\n");
    initialize_broyden_ALO(Broyden_local, ix1, ix2, ix3, Temp);

    /*for (int i = 1; i<=nmap-1; ++i)
      for (int ii = 1; ii<=nmap-1; ++ii)
        Broyden_local[i][ii] = M[i][ii];
    // But the last one has to be replaced with all 1s:
    for (int ii = 1; ii<=nmap-1; ++ii)
      Broyden_local[nmap-1][ii] = 1.0;*/
    
    invert(Broyden_local, Broyden_local, nmap-1, nmap-1);
    if (ix3 == x3i_control)
      printf("atom::pops_broyden : first iteration, initial Broyden's matrix set up using the approximate lambda matrix!\n");
  }

  // Then you update the populations: 
  
  fp_t * correction = multiply_vector(Broyden_local, F, nmap-1);
  
  fp_t delta = 0;

  for (int i = 0; i<nmap-1; ++i){
    fp_t rel_delta = abs(correction[i]) / pop[ix1][ix2][ix3].n[zmap[i-1]][lmap[i-1]];
    delta = (rel_delta > delta) ? rel_delta : delta;
  }
  
  for (int i = 0; i<nmap-1; ++i){
    correction[i+1] *= -1;
    pop[ix1][ix2][ix3].n[zmap[i]][lmap[i]] += correction[i+1];
    //printf("x3 = %d z = %d l = %d, n_l = %e \n", ix3, zmap[i], lmap[i], pop[ix1][ix2][ix3].n[zmap[i]][lmap[i]]);
  }



  //Store old F: // This actually we do not need, but keep it so far 
  //for (int i = 1; i<=nmap-1; ++i)
  //  F_old[ix1][ix2][ix3][i] = F[i];

  // Print F to debug:
  if (ix3 == x3i_control){
    for (int i = 1; i<=nmap-1; ++i)
      printf ("F[%d] = %e\n",i, F[i]);
  }
  delete [](F+1);
  
  F = compute_F_broyden(ix1, ix2, ix3, Temp, J, nrm);
  // Then compute the defect:
  fp_t * defect = new fp_t[nmap-1] -1;
  for (int i = 1; i<=nmap-1; ++i)
    defect[i] = F[i];// - F_old[ix1][ix2][ix3][i];

  // And finally, update the matrix:

  Broyden_update(Broyden_local, defect, correction, nmap-1);

  delete [](defect + 1);
  
  delete [](F+1);
  
  delete [](correction +1);
  
  return delta;
  
}

  return 0;
  
}

int atom::initialize_broyden_ALO(fp_t ** broyden_matrix, int ix1, int ix2, int ix3, fp_t T){

  // This is the function which will initialize the Broyden matrix:

  fp_t * J = Jb[ix1][ix2][ix3];
  fp_t * L = Ls[ix1][ix2][ix3];
  fp_t * nrm = norm[ix1][ix2][ix3];

  // Now let us go, transition by transition:

  for (int i = 0; i<nmap-1; ++i){

    // Get level index and ionization stage:
    int l = lmap[i];
    int z = zmap[i];


    // For each combination of levels:

    for (int ll = 0; ll < nl[z]; ++ll){

      fp_t JJ = (tmap[z][l][ll])?J[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] :-1.0;   // angular and frequency integrated intensity
      fp_t LL = (tmap[z][l][ll])?L[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] : 0.0;


      // First transitions from l to ll, which means negative
      fp_t Radiative_rates = R_ij(z, l, ll, JJ);
      fp_t Collisional_rates = C_ij_dummy(z, l, ll, T);
      broyden_matrix[i+1][i+1] -= (Radiative_rates + Collisional_rates); 

      // But now also the modifications due to the final derivative:

      fp_t s_der = 0;

      fp_t n_l = pop[ix1][ix2][ix3].n[z][l];
      fp_t n_ll = pop[ix1][ix2][ix3].n[z][ll];

      // Derivative of the source function S_l_ll with respect to n_l
      s_der = n_ll * A[z][l][ll] / B[z][ll][l] / (n_ll - n_l * B[z][l][ll] / B[z][ll][l]) / (n_ll - n_l * B[z][l][ll] / B[z][ll][l]);
      s_der -= n_ll * A[z][ll][l] / B[z][l][ll] / (n_l - n_ll * B[z][ll][l] / B[z][l][ll]) / (n_l - n_ll * B[z][ll][l] / B[z][l][ll]);

      s_der = (l == ll) ? 0 : s_der;

      //if (ix3 == x3l){
      //  printf("z = %d l = %d, S_der = %e \n", z, l, s_der);
      //}

      broyden_matrix[i+1][i+1] -= n_l * B[z][l][ll] * LL * s_der;
      broyden_matrix[i+1][i+1] += n_ll * B[z][ll][l] * LL * s_der;

      // Transitions to this level:
      Radiative_rates = R_ij(z, ll, l, JJ);
      Collisional_rates = C_ij_dummy(z, ll, l, T);
      int dl = ll - l;

      s_der = n_l * A[z][ll][l] / B[z][l][ll] / (n_l - n_ll * B[z][ll][l] / B[z][l][ll]) / (n_l - n_ll * B[z][ll][l] / B[z][l][ll]);
      s_der -= n_l * A[z][l][ll] / B[z][ll][l] / (n_ll - n_l * B[z][l][ll] / B[z][ll][l]) / (n_ll - n_l * B[z][l][ll] / B[z][ll][l]);

      s_der = (l == ll) ? 0 : s_der;


      //if (ix3 == x3l){
      //  printf("z = %d l = %d, S_der = %e \n", z, l, s_der);
      //}

      broyden_matrix[i+1][i+1+dl] += (Radiative_rates + Collisional_rates + n_ll * B[z][ll][l] * LL * s_der - n_l * B[z][l][ll] * LL * s_der); 

    }

  }
  for (int ii = 1; ii<=nmap-1; ++ii)
      broyden_matrix[nmap-1][ii] = 1.0;

  return 0;
}

fp_t * atom::compute_F_broyden(int ix1, int ix2, int ix3, fp_t Temp, fp_t * J, fp_t * nrm){

  // This is the function which computes F, i.e. acts on the level population with the 'rate' matrix in order to give us the value F

  // First get us the rate matrix:

  fp_t **M=ft2dim(1,nmap,1,nmap); // Allocate the rate matrix
  memset(M[1]+1,0,nmap*nmap*sizeof(fp_t)); // Set all elements to zero

  // Setup the rate equation:
  for(int i=0;i<nmap-1;++i){ // For all the 'levels' but the last one 
    uint16_t l=lmap[i]; // "current lvl
    uint08_t z=zmap[i]; // apropriate ionization stage
    
    // First account for b-b transitions within the particular ionization stage:
    for(int ll=0; ll<nl[z]; ++ll){ // For all the levels:
      fp_t JJ = (tmap[z][l][ll])?J[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] :-1.0;   // angular and frequency integrated intensity
      // Transitions from this level
      fp_t Radiative_rates = R_ij(z, l, ll, JJ);//A[z][l][ll] + B[z][l][ll] * JJ;
      fp_t Collisional_rates = C_ij_dummy(z, l, ll, Temp);
      M[i+1][i+1] -= (Radiative_rates + Collisional_rates); 

      // Transitions to this level:
      Radiative_rates = R_ij(z, ll, l, JJ);//A[z][ll][l] + B[z][ll][l] * JJ;
      Collisional_rates = C_ij_dummy(z, ll, l, Temp);
      int dl = ll - l;
      M[i+1][i+1+dl] += (Radiative_rates + Collisional_rates);  
    }
    
    // Then account for b-f and f-b transitions:
    if(int dl=nl[z]-l){              // ground level of next ionization stage (if not this level)
      // ionization out of level z,l
      fp_t JJ=(tmap[z][l][nl[z]])?J[tmap[z][l][nl[z]]]:-1.0;     // angular and frequency integrated intensity
      /*
      // Transitions from this level:
      fp_t Radiative_rates = R_i_cont(z, l, JJ);
      fp_t Collisional_rates = C_i_cont(z, l, Temp);

      M[i+1][i+1] -= (Radiative_rates + Collisional_rates);

      Radiative_rates = R_cont_i (z, l, JJ);
      Collisional_rates = C_cont_i (z, l, Temp);

      M[i+1][nl[z]+1] += (Radiative_rates + Collisional_rates);*/
    }
  }
 for(int ii=1;ii<=nmap;++ii) M[nmap][ii]=1.0;  // final equation: sum(n)=cst (number of particles remains the same)

  // This is 'classic' rate matrix, now apply it to the level populations:

  fp_t * F = new fp_t [nmap-1] - 1;

  for (int i = 1; i<=nmap-2; ++i){
    F[i] = 0.0;
    for (int ii = 1; ii<=nmap-1; ++ii)
      F[i] += M[i][ii] * pop[ix1][ix2][ix3].n[zmap[ii-1]][lmap[ii-1]];
  }
  F[nmap-1] = -get_pop(ix1, ix2, ix3, 0);
  for (int ii = 1; ii <= nmap-1; ++ii)
    F[nmap-1] += pop[ix1][ix2][ix3].n[zmap[ii-1]][lmap[ii-1]];
  
  del_ft2dim(M, 1, nmap, 1, nmap);

  return F;

}