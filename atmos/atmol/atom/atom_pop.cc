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

// -------------------------------------------------------------------------------------------------------------------------

fp_t * atom::newpops(int32_t x1i,int32_t x2i,int32_t x3i, int alo){

// Find and propose new values of the level populations given the location and the method

  if(Jb && NLTE){ // If it is an NLTE species with defined scattering integral

    fp_t Temp = parent_atm->get_T(x1i,x2i,x3i);
    fp_t ne   = parent_atm->get_Ne(x1i,x2i,x3i);

    fp_t *J_lu=Jb[x1i][x2i][x3i]; // scattering integral lower -> upper 
    fp_t *J_ul=Ju[x1i][x2i][x3i]; // scattering integral upper -> lower (same for the lines, different for b-f)
    fp_t *L=Ls[x1i][x2i][x3i]; // approximate lambda operator
    fp_t *nrm = norm[x1i][x2i][x3i]; // norm of the profile: 

    fp_t **M=ft2dim(1,nmap,1,nmap); // Allocate and zero the rate matrix
    memset(M[1]+1,0,nmap*nmap*sizeof(fp_t));

    fp_t* b = new fp_t[nmap]; // Allocate and zero r.h.s. of the SE equation
    memset(b,0,nmap*sizeof(fp_t));
  
    // Setup the rate equation:
    for(int i=0;i<nmap-1;++i){ // For all the 'levels' but the last one 
      
      uint16_t l=lmap[i]; // "current lvl
      uint08_t z=zmap[i]; // apropriate ionization stage

      // First account for b-b transitions within the particular ionization stage:
      for(int ll=0; ll<nl[z]; ++ll){ // For all the levels:

        fp_t JJ = (tmap[z][l][ll]) ? J_lu[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] :-1.0;   // angular and frequency integrated intensity
        fp_t LL = (tmap[z][l][ll]) ? L[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] : 0.0;

        // Transitions from this level
        fp_t Radiative_rates;
        if (alo == 1)
          Radiative_rates = R_ij_local_ALO(z, l, ll, JJ, LL, pop[x1i][x2i][x3i].n[z][l], pop[x1i][x2i][x3i].n[z][ll]);
        else
          Radiative_rates = R_ij(z,l,ll,JJ);
        fp_t Collisional_rates = C_ij(z, l, ll, Temp, ne);
        //if (Z==1) Collisional_rates += C_ij_H(z, l, ll, Temp, fetch_population(x1i, x2i, x3i, 0, 0)); // Modify for H collisions
        
        M[i+1][i+1] -= (Radiative_rates + Collisional_rates); 
        
        // Transitions to this level:
        if (alo)
          Radiative_rates = R_ij_local_ALO(z, ll, l, JJ, LL, pop[x1i][x2i][x3i].n[z][ll], pop[x1i][x2i][x3i].n[z][l]);
        else 
          Radiative_rates = R_ij(z,ll,l,JJ);
        Collisional_rates = C_ij(z, ll, l, Temp, ne);
        //if (Z==1) Collisional_rates += C_ij_H(z, ll, l, Temp, fetch_population(x1i, x2i, x3i, 0, 0));
        int dl = ll - l;

        M[i+1][i+1+dl] += (Radiative_rates + Collisional_rates);  
 
      }
      // Then account for b-f and f-b transitions:
      if(int dl=nl[z]-l){              // ground level of next ionization stage (if not this level)
        
        // ionization out of level z,l
        fp_t JJ=(tmap[z][l][nl[z]])?J_lu[tmap[z][l][nl[z]]]:-1.0;     // angular and frequency integrated intensity
        
        // Transitions from this level:
        fp_t Radiative_rates = R_i_cont(z, l, JJ, Temp);
        fp_t Collisional_rates = C_i_cont(z, l, Temp, ne);

        M[i+1][i+1] -= (Radiative_rates + Collisional_rates);
        // But now the same these rates should populate the continuum:
        M[i+1+dl][i+1] += (Radiative_rates + Collisional_rates);
        
        JJ=(tmap[z][l][nl[z]])?J_ul[tmap[z][l][nl[z]]]:-1.0;     // angular and frequency integrated intensity
        
        Radiative_rates = R_cont_i(z, l, JJ, Temp, ne);
        Collisional_rates = C_cont_i(z, l, Temp, ne);
        
        M[i+1][i+1+dl] += (Radiative_rates + Collisional_rates);
        // And these same rates should depopulate the continuum:
        M[i+1+dl][i+1+dl] -= (Radiative_rates + Collisional_rates);
      }
    }
  

    // Pick which level to replace with the conservation equation:
    // Default is the last one, but that is the poor choice, you want to replace one with the largest population.
    int level_to_replace = nmap-1; 
    fp_t maxpop = pop[x1i][x2i][x3i].n[zmap[level_to_replace]][lmap[level_to_replace]];
    
    for (int i=0;i<nmap-1;++i){
      if (pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] > maxpop){
        maxpop = pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]];
        level_to_replace = i;
      }
    }
  
    // Do the replacement
    for(int ii=1;ii<=nmap;++ii) M[level_to_replace+1][ii]=1.0;
      b[level_to_replace] = pop[x1i][x2i][x3i].Na;
 
    // Let us invert the matrix: 
    fp_t * M_to_solve = M[1] +1;
    fp_t * M_LU = new fp_t [nmap * nmap];
    fp_t * solution = new fp_t [nmap];  
    Crout(nmap,M_to_solve, M_LU);
    solveCrout(nmap,M_LU,b,solution);

    del_ft2dim(M,1,nmap,1,nmap);
    delete []M_LU;
    delete []b;
  
    return solution; // If not (J && NLTE)
  }
  else{ // Else it's a LTE thing, so it might have levels

    fp_t * old_pops = new fp_t [nmap]; //i indexes from zero
    for (int i=0; i<nmap;++i)
      old_pops[i] = pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]];

    return old_pops;
  }
}

// -----------------------------------------------------------------------------------------------------------------

fp_t atom::get_mapped_pop(int x1i, int x2i, int x3i, int i){

  if (nmap){
    int l=lmap[i]; // "current lvl
    int z=zmap[i]; // apropriate ionization stage

    return pop[x1i][x2i][x3i].n[z][l];
  }
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------

int atom::set_pop(int x1i, int x2i, int x3i, int z, int l, fp_t pop_in){

  pop[x1i][x2i][x3i].n[z][l] = pop_in;

  return 0; // Still can't make up my mind about void vs int
}

// -----------------------------------------------------------------------------------------------------------------

