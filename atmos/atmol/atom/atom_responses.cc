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

const int N_atm_parameters = 7;

void atom::compute_nlte_population_responses(){

  if(ntr){

    clock_t begin = clock();
    io.msg(IOL_INFO, "atom::now computing nlte population responses in analytical way \n ---- version in progress --- \n");

    int x1i = 1;
    int x2i = 1;
    int N_depths = x3h - x3l + 1;
    int N_responses = (x3h-x3l+1)*nmap; // Total number of responses
    
    // Lets us start with 'beta' coefficients. We stick to the temperature so far.
    for (int x3i = x3l; x3i<= x3h; ++x3i){ // For each depth point.
      // You are perturbing the temperature @ k and looking @ what is happening:
        for (int i=0; i<nmap; ++i){
          // Perturbations affect only locally, i.e. only beta at x3i is going to be non zero:
          for (int ii=0; ii<nmap; ++ii){
            beta_Temp[x3i][(x3i-x3l) * nmap + i+1] += (pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] * derivative_collisions_full_temp(x1i, x2i, x3i, i, ii) -
                  pop[x1i][x2i][x3i].n[zmap[ii]][lmap[ii]] * derivative_collisions_full_temp(x1i, x2i, x3i, ii, i));
          }
        }
    }
    // Modification, to acount for the last row of the SE:
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        if (x3k == x3i)
          beta_Temp[x3k][(x3k-x3l+1)*nmap] = derivative_active_population(x1l, x2l, x3k);
        else 
          beta_Temp[x3k][(x3i-x3l+1)*nmap] = 0.0;
      }
// ---------------------------------------------------------------------------------------------------------------------------------------------------------

    // The main change in beta_density should be due to the change in the collisions with e- and neutral hydrogen

    // Should be identical to the previous one:
    for (int x3i = x3l; x3i<= x3h; ++x3i){
        for (int i=0; i<nmap; ++i){
          for (int ii=0; ii<nmap; ++ii){
            beta_density[x3i][(x3i-x3l) * nmap + i+1] += (pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] * derivative_collisions_full_density(x1i, x2i, x3i, i, ii) -
                  pop[x1i][x2i][x3i].n[zmap[ii]][lmap[ii]] * derivative_collisions_full_density(x1i, x2i, x3i, ii, i));
          }
        }
    }
    // Modification, to acount for the last row of the SE:
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        if (x3k == x3i)
          beta_density[x3k][(x3k-x3l+1)*nmap] = derivative_active_population_density(x1l, x2l, x3k);
        else 
          beta_density[x3k][(x3i-x3l+1)*nmap] = 0.0;
      }

// ---------------------------------------------------------------------------------------------------------------------------------------------------------

   
    // Now we manipulate a bit the left hand side, what remains to be added, are the rates:
    fp_t * J_lu, * J_ul, * nrm;

    for (int l=x3l; l<=x3h; ++l){ // At each depth of the atmosphere.

      J_lu=Jb[x1l][x2l][l]; // integrated intensities in down-up and up-down directions:
      J_ul=Ju[x1l][x2l][l];
      nrm = norm[x1l][x2l][l];

      fp_t Temp = fetch_temperature(x1i, x2i, l);
      fp_t ne = fetch_Ne(x1i, x2i, l);
      
      // For each level
      for (int i = 1; i<=nmap-1; ++i){ // The same thing as in the solution of the populations, all the levels but the last one :)

        int l_i=lmap[i-1]; // "current lvl
        int z=zmap[i-1]; // apropriate ionization stage

        // First account for b-b transitions within the particular ionization stage:
        for(int l_ii=0; l_ii<nl[z]; ++l_ii){ // For all the levels:

          fp_t JJ = (tmap[z][l_i][l_ii]) ? J_lu[tmap[z][l_i][l_ii]] / nrm[tmap[z][l_i][l_ii]] :-1.0;   // angular and frequency integrated intensity
          
          // Transitions from this level
          fp_t Radiative_rates = 1.0 * R_ij(z, l_i, l_ii, JJ);  
          fp_t Collisional_rates = 1.0 * C_ij(z, l_i, l_ii, Temp, ne);
          //Collisional_rates += C_ij_H(z, l_i, l_ii, Temp, fetch_population(x1i, x2i, l, 0, 0)); // Modify for collisions with H
          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i] -= (Radiative_rates + Collisional_rates); 

          // Transitions to this level:
          Radiative_rates = 1.0 * R_ij(z, l_ii, l_i, JJ);
          Collisional_rates = 1.0 * C_ij(z, l_ii, l_i, Temp, ne);
          //Collisional_rates += C_ij_H(z, l_ii, l_i, Temp, fetch_population(x1l, x2l, l, 0, 0));
          int dl = l_ii - l_i;
          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i + dl] += (Radiative_rates + Collisional_rates);          
        }    
        // Then account for b-f and f-b transitions:
        if(int dl=nl[z]-l_i){              // ground level of next ionization stage (if not this level)
          // ionization out of level z,l
          fp_t JJ=(tmap[z][l_i][nl[z]])?J_lu[tmap[z][l_i][nl[z]]]:-1.0;     // angular and frequency integrated intensity 
          // Transitions from this level:
          fp_t Radiative_rates = 1.0 * R_i_cont(z, l_i, JJ, Temp);
          fp_t Collisional_rates = C_i_cont(z, l_i, Temp, ne);
          
          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i] -= (Radiative_rates + Collisional_rates); 
          response_matrix[(l-x3l) * nmap + i + dl][(l-x3l) * nmap + i] += (Radiative_rates + Collisional_rates); 
          
          JJ=(tmap[z][l_i][nl[z]])?J_ul[tmap[z][l_i][nl[z]]]:-1.0;     // angular and frequency integrated intensity
    
          Radiative_rates = 1.0 * R_cont_i(z, l_i, JJ, Temp, ne);
          Collisional_rates = C_cont_i(z, l_i, Temp, ne);

          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i+ dl] += (Radiative_rates + Collisional_rates);
          response_matrix[(l-x3l) * nmap + i + dl][(l-x3l) * nmap + i+ dl] -= (Radiative_rates + Collisional_rates);  
        }  
      }
    }
      // Polish the last line, it should be 1:
    for (int l = x3l; l<=x3h; ++l)
      for (int ll = x3l; ll<=x3h; ++ll)
        for (int ii = 1; ii <= nmap; ++ii){
          if (l == ll)
            response_matrix[(l-x3l) * nmap + nmap][(l-x3l) * nmap + ii] = 1.0;
          else 
            response_matrix[(l-x3l) * nmap + nmap][(ll-x3l) * nmap + ii] = 0.0;
    }
// ------------------------------------------------------------------------------------------------------------------------------------------------
  
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("Time spent on setting up linear system = %f \n", time_spent);

    // Finally we need to compute all the responses and that really means invert response matrix and then multiply by the right hand side:
   
    fp_t * M_to_solve = response_matrix[1] +1;
    fp_t * M_LU = new fp_t [N_depths * N_depths * nmap * nmap];
    fp_t * b; // Right hand side
    fp_t * solution = new fp_t [N_depths * nmap];
    
    //printf(">> Now solving via LU decomposition << \n");

    Crout(nmap * N_depths,M_to_solve, M_LU);

    begin = clock();
    for (int x3i = x3l; x3i <= x3h; ++x3i){

      b = beta_Temp[x3i] + 1;
      
      solveCrout(nmap * N_depths,M_LU,b,solution);

      // Write down the solution for Temp:
      for (int x3ii = x3l;x3ii<=x3h; ++x3ii)
        for (int i = 1; i<=nmap; ++i){
          level_responses[1][(x3ii - x3l) * nmap + i][x3i] = solution[(x3ii-x3l) * nmap + i -1];
        }

      // Now all the same again, but for density:
      b = beta_density[x3i] + 1;      
      solveCrout(nmap * N_depths,M_LU,b,solution);
      // Write down the solution for density:
      for (int x3ii = x3l;x3ii<=x3h; ++x3ii)
        for (int i = 1; i<=nmap; ++i){
          level_responses[2][(x3ii - x3l) * nmap + i][x3i] = solution[(x3ii-x3l) * nmap + i -1];
      }

      // Microturbulent velocity
      b = beta_v_micro[x3i] + 1;
      solveCrout(nmap * N_depths,M_LU,b,solution);
      // Write down the solution for microturbulent velocity
      for (int x3ii = x3l;x3ii<=x3h; ++x3ii)
        for (int i = 1; i<=nmap; ++i){
          level_responses[3][(x3ii - x3l) * nmap + i][x3i] = solution[(x3ii-x3l) * nmap + i -1];
      }
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("Time spent on solving linear system = %f \n", time_spent);
    io.msg(IOL_INFO, "responses computed for the atom with Z = %d using analytical approximation \n", Z);
    delete []M_LU;
    delete []solution;
  }

}

// Now the same function but much simpler, only taking into account collisional stuff:

void atom::compute_lte_population_responses(){

  if(ntr){

    io.msg(IOL_INFO, "atom::now computing lte population responses\n");

    // Small test, we put it here because why not:

    int x1i = 1;
    int x2i = 1;
    int N_depths = x3h - x3l + 1;
    int N_responses = (x3h-x3l+1)*nmap; // Total number of responses

    // LTE, so cancel everything:
    memset(beta_Temp[x3l]+1,0,N_depths*N_responses*sizeof(fp_t));
    memset(beta_density[x3l]+1,0,N_depths*N_responses*sizeof(fp_t));
    memset(beta_v_micro[x3l]+1,0,N_depths*N_responses*sizeof(fp_t));
     
    // Lets us start with 'beta' coefficients. We stick to the temperature so far.
    for (int x3i = x3l; x3i<= x3h; ++x3i){ // For each depth point.
      // You are perturbing the temperature @ k and looking @ what is happening:
        for (int i=0; i<nmap; ++i){
          // Perturbations affect only locally, i.e. only beta at x3i is going to be non zero:
          for (int ii=0; ii<nmap; ++ii){
            beta_Temp[x3i][(x3i-x3l) * nmap + i+1] += (pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] * derivative_collisions_full_temp(x1i, x2i, x3i, i, ii) -
                  pop[x1i][x2i][x3i].n[zmap[ii]][lmap[ii]] * derivative_collisions_full_temp(x1i, x2i, x3i, ii, i));
          }
        }
    }

    // Modification, to acount for the last row of the SE:
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        if (x3k == x3i)
          beta_Temp[x3k][(x3k-x3l+1)*nmap] = derivative_active_population(x1l, x2l, x3k);
        else 
          beta_Temp[x3k][(x3i-x3l+1)*nmap] = 0.0;
      }

// ---------------------------------------------------------------------------------------------------------------------------------------------------------

    // The main change in beta_density should be due to the change in the collisions with e- and neutral hydrogen

    // Should be identical to the previous one:
    for (int x3i = x3l; x3i<= x3h; ++x3i){
        for (int i=0; i<nmap; ++i){
          for (int ii=0; ii<nmap; ++ii){
            beta_density[x3i][(x3i-x3l) * nmap + i+1] += (pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] * derivative_collisions_full_density(x1i, x2i, x3i, i, ii) -
                  pop[x1i][x2i][x3i].n[zmap[ii]][lmap[ii]] * derivative_collisions_full_density(x1i, x2i, x3i, ii, i));
          }
        }
    }

    // Modification, to acount for the last row of the SE:
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        if (x3k == x3i)
          beta_density[x3k][(x3k-x3l+1)*nmap] = derivative_active_population_density(x1l, x2l, x3k);
        else 
          beta_density[x3k][(x3i-x3l+1)*nmap] = 0.0;
      }



// -------------------------------------------------------------------------------------------------------------------------------------------------
    memset(response_matrix[1]+1,0,N_responses*N_responses*sizeof(fp_t));

    for (int l=x3l; l<=x3h; ++l){ // At each depth of the atmosphere.

     
      fp_t Temp = fetch_temperature(x1i, x2i, l);
      fp_t ne = fetch_Ne(x1i, x2i, l);
      
      // For each level
      for (int i = 1; i<=nmap-1; ++i){ // The same thing as in the solution of the populations, all the levels but the last one :)

        int l_i=lmap[i-1]; // "current lvl
        int z=zmap[i-1]; // apropriate ionization stage

        // First account for b-b transitions within the particular ionization stage:
        for(int l_ii=0; l_ii<nl[z]; ++l_ii){ // For all the levels:
        
          // Transitions from this level
          fp_t Collisional_rates = 1.0 * C_ij(z, l_i, l_ii, Temp, ne);
          //if (l==x3l)
          //  printf("From level %d to level %d, rate is: %e \n", l_i, l_ii, Collisional_rates);
          //Collisional_rates += C_ij_H(z, l_i, l_ii, Temp, fetch_population(x1i, x2i, l, 0, 0)); // Modify for collisions with H
          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i] -= Collisional_rates; 

          // Transitions to this level:
          Collisional_rates = 1.0 * C_ij(z, l_ii, l_i, Temp, ne);
          //if (l==x3l)
          //  printf("From level %d to level %d, rate is: %e \n", l_ii, l_i, Collisional_rates);
          //Collisional_rates += C_ij_H(z, l_ii, l_i, Temp, fetch_population(x1l, x2l, l, 0, 0));
          int dl = l_ii - l_i;
          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i + dl] += Collisional_rates;          
        }    
        // Then account for b-f and f-b transitions:
        if(int dl=nl[z]-l_i){              // ground level of next ionization stage (if not this level)
          // ionization out of level z,l
          fp_t Collisional_rates = C_i_cont(z, l_i, Temp, ne);
          
          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i] -= Collisional_rates; 
          response_matrix[(l-x3l) * nmap + i + dl][(l-x3l) * nmap + i] += Collisional_rates; 
              
          Collisional_rates = C_cont_i(z, l_i, Temp, ne);

          response_matrix[(l-x3l) * nmap + i][(l-x3l) * nmap + i+ dl] += Collisional_rates;
          response_matrix[(l-x3l) * nmap + i + dl][(l-x3l) * nmap + i+ dl] -= Collisional_rates;  
        }  
      }
    }
      // Polish the last line, it should be 1:
    for (int l = x3l; l<=x3h; ++l)
      for (int ll = x3l; ll<=x3h; ++ll)
        for (int ii = 1; ii <= nmap; ++ii){
          if (l == ll)
            response_matrix[(l-x3l) * nmap + nmap][(l-x3l) * nmap + ii] = 1.0;
          else 
            response_matrix[(l-x3l) * nmap + nmap][(ll-x3l) * nmap + ii] = 0.0;
    }

    /*if (strcmp(id,"Fe") == 0){
      printf("I AM IRON!\n");
      FILE * rm = fopen("Fe_response_matrix.dat","w");
      for (int x3i=x3l;x3i<=x3h;++x3i)
        for (int i=1;i<=nmap;++i){
          for (int x3ii=x3l;x3ii<=x3h;++x3ii)
            for (int ii=1;ii<=nmap;++ii)
              fprintf(rm,"%e ",response_matrix[(x3i-x3l)*nmap+i][(x3ii-x3l)*nmap+ii]);
        fprintf(rm,"\n");
        }
      fclose(rm);
      FILE * rt = fopen("Fe_beta_T.dat","w");
      for (int x3i=x3l;x3i<=x3h;++x3i){
          for (int x3ii=x3l;x3ii<=x3h;++x3ii)
            for (int ii=1;ii<=nmap;++ii)
              fprintf(rt,"%e ",beta_Temp[x3i][(x3ii-x3l)*nmap+ii]);
        fprintf(rt,"\n");
        }
      fclose(rt);
      FILE * rd = fopen("Fe_beta_dens.dat","w");
      for (int x3i=x3l;x3i<=x3h;++x3i){
          for (int x3ii=x3l;x3ii<=x3h;++x3ii)
            for (int ii=1;ii<=nmap;++ii)
              fprintf(rd,"%e ",beta_density[x3i][(x3ii-x3l)*nmap+ii]);
        fprintf(rd,"\n");
        }
      fclose(rd);
    }*/

    // Finally we need to compute all the responses and that really means invert response matrix and then multiply by the right hand side:
   
    // DEBUG: --- this is the version where we invert matrix by matrix at each depth. 

    for (int x3i=x3l;x3i<=x3h;++x3i){

      fp_t * M_to_solve = new fp_t[nmap*nmap];
      for (int i=1;i<=nmap;++i)
        for (int ii=1;ii<=nmap;++ii)
          M_to_solve[(i-1)*nmap+ii-1] = response_matrix[(x3i-x3l)*nmap+i][(x3i-x3l)*nmap+ii];
      fp_t * M_LU = new fp_t[nmap*nmap];
      fp_t * b = new fp_t[nmap]; // rhs
      fp_t * solution = new fp_t[nmap];

      Crout(nmap,M_to_solve,M_LU);
      int x3k=x3i;
      memcpy(b,beta_Temp[x3k]+(x3i-x3l)*nmap+1,nmap*sizeof(fp_t));
      solveCrout(nmap,M_LU,b,solution);
      for (int i=1;i<=nmap;++i)
        level_responses[1][(x3i-x3l)*nmap+i][x3k] = solution[i-1];
      
      memcpy(b,beta_density[x3k]+(x3i-x3l)*nmap+1,nmap*sizeof(fp_t));
      solveCrout(nmap,M_LU,b,solution);
      for (int i=1;i<=nmap;++i)
        level_responses[2][(x3i-x3l)*nmap+i][x3k] = solution[i-1];
      
      memcpy(b,beta_v_micro[x3k]+(x3i-x3l)*nmap+1,nmap*sizeof(fp_t));
      solveCrout(nmap,M_LU,b,solution);
      for (int i=1;i<=nmap;++i)
        level_responses[3][(x3i-x3l)*nmap+i][x3k] = solution[i-1];

      delete[]M_to_solve;
      delete[]M_LU;
      delete[]b;
      delete[]solution;

    }
    /*fp_t * M_to_solve = response_matrix[1] +1;
    fp_t * M_LU = new fp_t [N_depths * N_depths * nmap * nmap];
    fp_t * b; // Right hand side
    fp_t * solution = new fp_t [N_depths * nmap];
    }
   
// ------------------------------------------------------------------------------------------------------------------------------------------------
    Crout(nmap * N_depths,M_to_solve, M_LU);

    for (int x3i = x3l; x3i <= x3h; ++x3i){

      b = beta_Temp[x3i] + 1;
      
      solveCrout(nmap * N_depths,M_LU,b,solution);

      // Write down the solution for Temp:
      for (int x3ii = x3l;x3ii<=x3h; ++x3ii)
        for (int i = 1; i<=nmap; ++i){
          level_responses[1][(x3ii - x3l) * nmap + i][x3i] = solution[(x3ii-x3l) * nmap + i -1];
        }

      // Now all the same again, but for density:
      b = beta_density[x3i] + 1;
      
      solveCrout(nmap * N_depths,M_LU,b,solution);

      // Write down the solution for Temp:
      for (int x3ii = x3l;x3ii<=x3h; ++x3ii)
        for (int i = 1; i<=nmap; ++i){
          level_responses[2][(x3ii - x3l) * nmap + i][x3i] = solution[(x3ii-x3l) * nmap + i -1];
      }

      // Microturbulent velocity
      b = beta_v_micro[x3i] + 1;
      
      solveCrout(nmap * N_depths,M_LU,b,solution);

      // Write down the solution for Temp:
      for (int x3ii = x3l;x3ii<=x3h; ++x3ii)
        for (int i = 1; i<=nmap; ++i){
          level_responses[3][(x3ii - x3l) * nmap + i][x3i] = solution[(x3ii-x3l) * nmap + i -1];
      }


    }
    */
    io.msg(IOL_INFO, "responses computed for the atom with Z = %d using analytical approximation \n", Z);

    //if (strcmp(id,"Fe")==0){
    //  print_population_responses("Fe_population_responses.dat",x3l,x3h);
    //  print_populations();
    //}
    //delete []M_LU;
    //delete []solution;
  }
}

fp_t atom::derivative_collisions_Temp(int x1i, int x2i, int x3i, int from, int to){
  return 0;
}

// Now fully numerical derivatives, also including the numerical derivative with respect to the electron density:

fp_t atom::derivative_collisions_full_temp(int x1i, int x2i, int x3i, int from , int to){

  fp_t derivative = 0.0;
  fp_t local_T = fetch_temperature(x1i, x2i, x3i);
  // Numerically compute the derivative of the rates for fixed Ne density:
  parent_atm->set_Temp(x1i, x2i, x3i, local_T + 0.5 * delta_T);
  derivative = collisional_rates(x1i, x2i, x3i, from, to);
  parent_atm->set_Temp(x1i, x2i, x3i, local_T - 0.5 * delta_T);
  derivative -= collisional_rates(x1i, x2i, x3i, from, to);
  derivative /= delta_T;
  parent_atm->set_Temp(x1i, x2i, x3i, local_T);
  fp_t rates = collisional_rates(x1i, x2i, x3i, from, to);
  derivative += parent_atm->get_ne_lte_derivative(1,x1i,x2i,x3i) * collisional_rates_der_ne(x1i, x2i, x3i, from, to, rates);
  
  return derivative;
}

fp_t atom::derivative_collisions_full_density(int x1i, int x2i, int x3i, int from , int to){

  fp_t derivative = 0.0;
  fp_t Nt = fetch_Nt(x1i, x2i, x3i);
  fp_t Nt_step = delta_Nt_frac*Nt;
  fp_t Ne = fetch_Ne(x1i,x2i,x3i);
  // Numerically compute the derivative of the rates for fixed Ne density:
  parent_atm->set_Nt(x1i, x2i, x3i, Nt + 0.5 * Nt_step);
  derivative = collisional_rates(x1i, x2i, x3i, from, to);
  parent_atm->set_Nt(x1i, x2i, x3i, Nt - 0.5 * Nt_step);
  derivative -= collisional_rates(x1i, x2i, x3i, from, to);
  derivative /= Nt_step;
  parent_atm->set_Nt(x1i, x2i, x3i, Nt);
  fp_t rates = collisional_rates(x1i, x2i, x3i, from, to);
  derivative += parent_atm->get_ne_lte_derivative(2,x1i,x2i,x3i) * collisional_rates_der_ne(x1i, x2i, x3i, from, to, rates);
  
  return derivative;
}

fp_t atom::collisional_rates(int x1i, int x2i, int x3i, int from, int to){
  // Find z and l of levels
  int z_from = zmap[from];
  int l_from = lmap[from];

  int z_to = zmap[to];
  int l_to = lmap[to];

  fp_t Temp = fetch_temperature(x1i, x2i, x3i);
  fp_t Ne = fetch_Ne(x1i, x2i, x3i);
  
  // If they are from the same ionization state:
  if (z_from == z_to)
    return C_ij(z_from, l_from, l_to, Temp, Ne);//+ C_ij_H(z_from, l_from, l_to, Temp, fetch_population(x1i, x2i, x3i, 0, 0));

  else if (z_from - z_to == -1 && l_to == 0)
    return C_i_cont(z_from, l_from, Temp, Ne);

  else if (z_from - z_to == 1 && l_from == 0)
    return C_cont_i(z_to, l_to, Temp, Ne);

  return 0;
}

fp_t atom::collisional_rates_der_ne(int x1i, int x2i, int x3i, int from, int to, fp_t rate_value){

  int z_from = zmap[from];
  int l_from = lmap[from];

  int z_to = zmap[to];
  int l_to = lmap[to];

  fp_t Ne = fetch_Ne(x1i, x2i, x3i);
  
  // If they are from the same ionization state:
  if (z_from == z_to)
    return rate_value/Ne;

  else if (z_from - z_to == -1 && l_to == 0)
    return rate_value/Ne;

  else if (z_from - z_to == 1 && l_from == 0)
    return rate_value*2.0/Ne;

  return 0;
}

int atom::responses_setup(){

  if (ntr){

    int N_depths = x3h - x3l + 1;
    int x1i = 1;
    int x2i = 1;
    beta_Temp = ft2dim(1, N_depths, 1, N_depths * nmap);
    beta_density = ft2dim(1, N_depths, 1, N_depths * nmap);
    beta_v_micro = ft2dim(1, N_depths, 1, N_depths * nmap);
    beta_v_r = ft2dim(1, N_depths, 1, N_depths * nmap);
    response_matrix = ft2dim(1, N_depths * nmap, 1, N_depths * nmap);
    level_responses = ft3dim(1,N_atm_parameters,1, N_depths * nmap, 1, N_depths);
    ionization_stage_responses = ft3dim(1,N_atm_parameters,x3l,x3h,0,Z);
    memset(ionization_stage_responses[1][x3l]+0,0,N_atm_parameters*(x3h-x3l+1)*(Z+1)*sizeof(fp_t));
    
  }

  return 0;
}

int atom::responses_init(){

  // What to use as the setup for responses? We assume that nlte solution has been done before, and therefore we will use ntr as a criterion for existing transitions
  if (ntr){
    int N_depths = x3h - x3l + 1;
    int x1i = 1;
    int x2i = 1;
    memset(beta_Temp[1] + 1,0, N_depths*N_depths*nmap*sizeof(fp_t));
    memset(beta_density[1] + 1,0, N_depths*N_depths*nmap*sizeof(fp_t));
    memset(beta_v_micro[1] + 1,0, N_depths*N_depths*nmap*sizeof(fp_t));
    memset(beta_v_r[1] + 1,0, N_depths*N_depths*nmap*sizeof(fp_t));
    memset(level_responses[1][1]+1, 0, N_atm_parameters*N_depths*nmap*N_depths*sizeof(fp_t));
    memset(response_matrix[1]+1,0,N_depths*N_depths*nmap*nmap*sizeof(fp_t));
    //memset(ionization_stage_responses[1][x3l]+0,0,N_atm_parameters*(x3h-x3l+1)*(Z+1)*sizeof(fp_t));
  }
  return 0;
}

int atom::responses_clear(){

  if (ntr){
    int N_depths = x3h - x3l + 1;
    int x1i = 1;
    int x2i = 1;
    del_ft3dim(level_responses,1,N_atm_parameters, 1, N_depths*nmap, 1, N_depths);
    del_ft2dim(response_matrix, 1, N_depths * nmap, 1, N_depths * nmap);
    del_ft2dim(beta_Temp, 1, N_depths, 1, N_depths * nmap);
    del_ft2dim(beta_density, 1, N_depths, 1, N_depths * nmap);
    del_ft2dim(beta_v_micro, 1, N_depths, 1, N_depths * nmap);
    del_ft2dim(beta_v_r, 1, N_depths, 1, N_depths * nmap);
    del_ft3dim(ionization_stage_responses,1,N_atm_parameters,x3l,x3h,0,Z);    
  }
  return 0;
}

void atom::add_to_ion_responses(int param,int x3i, fp_t sign){

  for (int z=0;z<=Z;++z)
    ionization_stage_responses[param][x3i][z] += pop[x1l][x2l][x3i].N[z] * sign;

}
void atom::divide_ion_responses(int param, int x3i, fp_t step){

  for (int z=0;z<=Z;++z)
    ionization_stage_responses[param][x3i][z] /= step;

}

int atom::add_response_contributions_new(fp_t *** I, fp_t ** response_to_op, fp_t ** response_to_em, fp_t *** opp, fp_t *** em, fp_t lambda, fp_t lambda_w, fp_t theta, fp_t phi, fp_t angular_weight, 
  fp_t *** vlos, fp_t ***** op_pert_lte, fp_t ***** em_pert_lte){

  int ncmp = 1; // total number of Stokes components
  int istaugrid = parent_atm->is_tau_grid();

  if (!Jb || !norm || !NLTE) // No reason to proceed.
    return 0;

  // Otherwise we start by doing the same as we do in NLTE solution itself. We compute radiative rates.
  // This costs very little, so we can keep it separate for clarity.
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){

        // Fetch them in the more convenient arrays:
        fp_t *Jm = Jb[x1i][x2i][x3i];
        fp_t *Jn = Ju[x1i][x2i][x3i];
        fp_t *nrm = norm[x1i][x2i][x3i];

        // For each transition:
        for (int tr=1; tr<=ntr; ++tr){

          int z_state = inverse_tmap[tr][2];
          int lower_level = inverse_tmap[tr][3];
          int upper_level = inverse_tmap[tr][4]; 

          if (upper_level < nl[z_state]){ // Line transitions:

            fp_t line_energy = ee[z_state][upper_level] - ee[z_state][lower_level];
            fp_t line_opacity = (pop[x1i][x2i][x3i].n[z_state][lower_level] * B[z_state][lower_level][upper_level] - pop[x1i][x2i][x3i].n[z_state][upper_level] * B[z_state][upper_level][lower_level]) 
              * current_profile[x1i][x2i][x3i][tr] * line_energy * 0.25 / pi;
            fp_t line_factor = line_opacity / opp[x1i][x2i][x3i];
            fp_t elementary_contribution = lambda_w * angular_weight * I[x1i][x2i][x3i] * current_profile[x1i][x2i][x3i][tr] * line_factor * 0.25 / pi;

            Jm[tr] += elementary_contribution;
            Jn[tr] += elementary_contribution;
          }
          // If it is a bound-free transition. There is no profile. Instead, compute rates directly, so:
          else if (upper_level == nl[z_state]){
            fp_t sigma = (bf[z_state][lower_level]) ? bf[z_state][lower_level]->U(lambda) : 0.0;
            Jm[tr] += lambda_w * angular_weight * sigma * lambda / h /c * I[x1i][x2i][x3i];
            Jn[tr] += lambda_w * angular_weight * sigma * lambda / h /c * (I[x1i][x2i][x3i] + 2.0 * h * c * c / pow(lambda, 5.0)) * exp(-h * c / lambda / k / fetch_temperature(x1i, x2i, x3i));
          }
        }
  }

  // We continue by pre-computing the line profile derivatives as they are used multiple times throughout the procedure.
  // This is also relatively computationaly cheap, and, furthermore, cannot be made much (if at all) faster.
  fp_t **** profile_derivative_T;
  profile_derivative_T = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  memset(profile_derivative_T[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
  fp_t **** profile_derivative_density;
  profile_derivative_density = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  memset(profile_derivative_density[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
  fp_t **** profile_derivative_vt;
  profile_derivative_vt = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  memset(profile_derivative_vt[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
  fp_t **** profile_derivative_vr;
  profile_derivative_vr = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  memset(profile_derivative_vr[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
    
  for (int x3i=x3l;x3i<=x3h;++x3i)
    for (int tr = 1;tr<=ntr;++tr){
      int z_state = inverse_tmap[tr][2];
      int lower_level = inverse_tmap[tr][3];
      int upper_level = inverse_tmap[tr][4]; 
      if (upper_level < nl[z_state]){
        int z_i = z_state;
        int l_i = upper_level;
        int l_ii = lower_level;
        fp_t x,a; x=a=0;
        fp_t * x_der;
        fp_t * a_der;
        fp_t * dld_der;
        x_der = new fp_t [7]-1;
        a_der = new fp_t [7]-1;
        dld_der = new fp_t [7]-1;
        memset(x_der+1,0,7*sizeof(fp_t));
        memset(a_der+1,0,7*sizeof(fp_t));
        memset(dld_der+1,0,7*sizeof(fp_t));
        if (A[z_i][l_i][l_ii] > 1.0){
          compute_xa_der_scalar(x1l, x2l, x3i, z_i, l_i, l_ii, lambda, vlos, x, a, x_der, a_der, dld_der);
          fp_t F,H,dF,dH;
          fvoigtn(x,a,H,F,dH,dF);
          fp_t lam=fabs(h*c/(ee[z_i][l_i]-ee[z_i][l_ii]));
          fp_t dld = broad_dop(lam, fetch_temperature(x1l,x2l,x3i),fetch_vt(x1l,x2l,x3i));

          profile_derivative_T[x1l][x2l][x3i][tr] = dH * x_der[1] / dld - dF * a_der[1] / dld - H/dld/dld*dld_der[1];
          profile_derivative_density[x1l][x2l][x3i][tr] = dH * x_der[2] / dld - dF * a_der[2] / dld - H/dld/dld*dld_der[2];
          profile_derivative_vt[x1l][x2l][x3i][tr] = dH * x_der[3] / dld - dF * a_der[3] / dld - H/dld/dld*dld_der[3];
          profile_derivative_vr[x1l][x2l][x3i][tr] = (dH * x_der[4] / dld - dF * a_der[4] / dld - H/dld/dld*dld_der[4])*0.0*(cos(theta));
        }
        delete [](x_der+1);
        delete [](a_der+1);
        delete [](dld_der+1);
      }
  }
  fp_t *** profile_derivatives = ft3dim(x3l,x3h,1,7,1,ntr);
  memset(profile_derivatives[x3l][1]+1,0,(x3h-x3l+1)*7*ntr*sizeof(fp_t));
  
  for (int l=x3l;l<=x3h;++l)
    for (int tr=1;tr<=ntr;++tr){
      profile_derivatives[l][1][tr] = profile_derivative_T[x1l][x2l][l][tr];
      profile_derivatives[l][2][tr] = profile_derivative_density[x1l][x2l][l][tr];
      profile_derivatives[l][3][tr] = profile_derivative_vt[x1l][x2l][l][tr];
  }

  fp_t ** bf_op_derivative_ex;
  fp_t ** bf_em_derivative_ex;
  fp_t ** bb_op_derivative_ex;
  fp_t ** bb_em_derivative_ex;
  bf_op_derivative_ex = new fp_t* [x3h-x3l+1]-x3l;
  bf_em_derivative_ex = new fp_t* [x3h-x3l+1]-x3l;
  bb_op_derivative_ex = new fp_t* [x3h-x3l+1]-x3l;
  bb_em_derivative_ex = new fp_t* [x3h-x3l+1]-x3l;


  for (int ll=x3l;ll<=x3h;++ll){
    bf_op_derivative_ex[ll] = op_derivative_explicit(ll,0,0,lambda);
    bf_em_derivative_ex[ll] = em_derivative_explicit(ll,0,0,lambda);
    bb_op_derivative_ex[ll] = bb_op_derivative_explicit(x1l,x2l,ll,lambda,profile_derivatives[ll]);
    bb_em_derivative_ex[ll] = bb_em_derivative_explicit(x1l,x2l,ll,lambda,profile_derivatives[ll]);
  }

  // And then derivatives of op and em w.r.t. level populations, again pre-compute:
  fp_t ** d_op_d_ni = ft2dim(x3l,x3h,0,nmap-1);
  fp_t ** d_em_d_ni = ft2dim(x3l,x3h,0,nmap-1);
  fp_t ** d_op_d_ni_bb = ft2dim(x3l,x3h,0,nmap-1);
  fp_t ** d_em_d_ni_bb = ft2dim(x3l,x3h,0,nmap-1);
  for (int ll=x3l;ll<=x3h;++ll)
    for (int i=0;i<nmap;++i){
      d_op_d_ni[ll][i] = op_derivative_to_level(ll,i,lambda);
      d_em_d_ni[ll][i] = em_derivative_to_level(ll,i,lambda);
      d_op_d_ni_bb[ll][i] = bb_op_derivative_to_level(x1l,x2l,ll,i,lambda);
      d_em_d_ni_bb[ll][i] = bb_em_derivative_to_level(x1l,x2l,ll,i,lambda);
  }

  // Now we need an effective method to sum-up and distribute responses to the l.h.s and r.h.s.
  // Let's first try BOOM-ing it straight from the other function. We compute r.h.s separately

   // Lets deal with the "response matrix" separately:
  for (int l=x3l;l<=x3h;++l) // For each depth point
    for (int tr=1;tr<=ntr;++tr){ // For each posible transition

      // Look up what levels we are talking about:
      int z = inverse_tmap[tr][2];
      int ii = inverse_tmap[tr][3];
      int i = inverse_tmap[tr][4];
      int type = 0; // 0 don't account, 1 b-f, 2-bb
      if (i == nl[z]){
        if (lambda <= bf[z][ii]->getlambda_crit())
          type = 1;
      }
      else {
        if (A[z][i][ii] > 1.0)
          type = 2;
      }
      // So if it is a line:
      if (type == 2){ // If this is a line transition

        // First we add local, explicit perturbations which go to the rhs:

        // This does not change whatsoever, regardless of the type of the derivative
        fp_t constant_factor =  lambda_w * angular_weight * 0.25 / pi;
        fp_t line_energy = ee[z][i] - ee[z][ii];
        fp_t lam = h * c / line_energy;
        fp_t line_opacity = (pop[x1l][x2l][l].n[z][ii]*B[z][ii][i] - pop[x1l][x2l][l].n[z][i]*B[z][i][ii]) 
          * line_energy * 0.25 / pi;
        fp_t line_factor = line_opacity * current_profile[x1l][x2l][l][tr] / opp[x1l][x2l][l];      
        fp_t J_chunk = (B[z][ii][i]*pop[x1l][x2l][l].n[z][ii]-B[z][i][ii]*pop[x1l][x2l][l].n[z][i])*
        I[x1l][x2l][l] * lambda_w * angular_weight * current_profile[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr]
          * line_factor * 0.25 / pi;     

        // TEMPERATURE profile derivative, influencing locally:
        fp_t derivative = profile_derivative_T[x1l][x2l][l][tr];
        fp_t derived_part = 0.0; // This is a part which multiplies J_chunk to get the derivative of J_chunk.
        // It has the form profile' / profile + norm' / norm .... etc.   
        // Profile derivative:
        derived_part += derivative / current_profile[x1l][x2l][l][tr];
        // Norm derivative:
        derived_part -= norm_derivative[1][x1l][x2l][l][tr]/ norm[x1l][x2l][l][tr];
        // Profile derivative influencing line factor:
        derived_part +=  (1.0 - line_factor)/line_factor * line_opacity / opp[x1l][x2l][l] * derivative;
        // Background opacity derivative influencing the line factor:
        derived_part -= (op_pert_lte[1][l][x1l][x2l][l]+bf_op_derivative_ex[l][1]) / opp[x1l][x2l][l];
        // Then add everything to the appropriate beta
        beta_Temp[l][(l-x3l)*nmap + rmap[z][i]+1] -= J_chunk * derived_part;
        beta_Temp[l][(l-x3l)*nmap + rmap[z][ii]+1] += J_chunk * derived_part;

        // DENSITY profile derivative:
        derivative = profile_derivative_density[x1l][x2l][l][tr];
        derived_part = 0.0;
        // Profile derivative:
        derived_part += derivative / current_profile[x1l][x2l][l][tr];
        // Norm derivative:
        derived_part -= norm_derivative[2][x1l][x2l][l][tr]/ norm[x1l][x2l][l][tr];
        // Profile derivative influencing line factor:
        derived_part +=  (1.0 - line_factor)/line_factor * line_opacity / opp[x1l][x2l][l] * derivative;
        // Background opacity derivative influencing the line factor:
        derived_part -= (op_pert_lte[2][l][x1l][x2l][l]+bf_op_derivative_ex[l][2]) / opp[x1l][x2l][l];
        // Then add everything to the appropriate beta
        beta_density[l][(l-x3l)*nmap + rmap[z][i]+1] -= J_chunk * derived_part;
        beta_density[l][(l-x3l)*nmap + rmap[z][ii]+1] += J_chunk * derived_part;

        // V_micro profile derivative:
        derivative = profile_derivative_vt[x1l][x2l][l][tr];
        derived_part = 0.0;
        // Profile derivative:
        derived_part += derivative / current_profile[x1l][x2l][l][tr];
        // Norm derivative:
        derived_part -= norm_derivative[3][x1l][x2l][l][tr]/ norm[x1l][x2l][l][tr];
        // Profile derivative influencing line factor:
        derived_part +=  (1.0 - line_factor)/line_factor * line_opacity / opp[x1l][x2l][l] * derivative;
        // Then add everything to the appropriate beta
        beta_v_micro[l][(l-x3l)*nmap + rmap[z][i]+1] -= J_chunk * derived_part;
        beta_v_micro[l][(l-x3l)*nmap + rmap[z][ii]+1] += J_chunk * derived_part;

        // -----------------------------------------------------------------------------------------------------------------

        // Now we can pre-compute dR_dI, to shorten the notation:
        fp_t dR_dI = (B[z][ii][i]*pop[x1l][x2l][l].n[z][ii]-B[z][i][ii]*pop[x1l][x2l][l].n[z][i])*
          lambda_w * angular_weight * current_profile[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * line_factor * 0.25 / pi; 
        // Now time for the non-local explicit dependencies:
        for (int ll=x3l;ll<=x3h;++ll){

          fp_t op_ref = 1.0;
          if (istaugrid)
            op_ref = parent_atm->get_op_referent(x1l,x2l,ll);

          fp_t d_op_d_profile = (pop[x1l][x2l][ll].n[z][ii]*B[z][ii][i]-pop[x1l][x2l][ll].n[z][i]*B[z][i][ii])*h*c/lam*0.25/pi;
          fp_t d_em_d_profile = pop[x1l][x2l][ll].n[z][i]*A[z][i][ii]*h*c/lam*0.25/pi;

          // TEMPERATURE:
          // First main contributor is the change of the profile in other points:
          derivative = profile_derivative_T[x1l][x2l][ll][tr];
          fp_t d_op_d_qk = derivative*d_op_d_profile;
          fp_t d_em_d_qk = derivative*d_em_d_profile;
          // Add b-f contributions, and background contributors:
          d_op_d_qk += bf_op_derivative_ex[ll][1] + op_pert_lte[1][ll][x1l][x2l][ll];
          d_em_d_qk += bf_em_derivative_ex[ll][1] + em_pert_lte[1][ll][x1l][x2l][ll];

          fp_t remote_response = dR_dI * (response_to_op[l][ll] * d_op_d_qk + response_to_em[l][ll] * d_em_d_qk);
          beta_Temp[ll][(l-x3l) * nmap +rmap[z][i] +1] -= remote_response/op_ref;
          beta_Temp[ll][(l-x3l) * nmap +rmap[z][ii] +1] += remote_response/op_ref;

          // Then, an additional term if we are using taugrid as a primary grid:
          if (istaugrid){      
            fp_t op_ref_der = parent_atm->get_op_referent_der(1,ll,x1l,x2l,ll);
            beta_Temp[ll][(l-x3l) * nmap +rmap[z][i] +1] += dR_dI *  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + 
              response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            beta_Temp[ll][(l-x3l) * nmap +rmap[z][ii] +1] -= dR_dI *  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + 
              response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
          }
          
          // DENSITY:
          // First main contributor is the change of the profile i other points:
          derivative = profile_derivative_density[x1l][x2l][ll][tr];
          d_op_d_qk = derivative*d_op_d_profile;
          d_em_d_qk = derivative*d_em_d_profile;
          // Add b-f contributions + background contributors
          d_op_d_qk += bf_op_derivative_ex[ll][2] + op_pert_lte[2][ll][x1l][x2l][ll];;
          d_em_d_qk += bf_em_derivative_ex[ll][2] + em_pert_lte[2][ll][x1l][x2l][ll];;
          remote_response = dR_dI * (response_to_op[l][ll] * d_op_d_qk + response_to_em[l][ll] * d_em_d_qk);
          beta_density[ll][(l-x3l) * nmap +rmap[z][i] +1] -= remote_response/op_ref;
          beta_density[ll][(l-x3l) * nmap +rmap[z][ii] +1] += remote_response/op_ref;

          // Then, an additional term if we are using taugrid as a primary grid:
          if (istaugrid){      
            fp_t op_ref_der = parent_atm->get_op_referent_der(2,ll,x1l,x2l,ll);
            beta_density[ll][(l-x3l) * nmap +rmap[z][i] +1] += dR_dI *  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + 
              response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            beta_density[ll][(l-x3l) * nmap +rmap[z][ii] +1] -= dR_dI *  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + 
              response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
          }

          // MICROTURBULENT VELOCITY:

          // First main contributor is the change of the profile i other points:
          derivative = profile_derivative_vt[x1l][x2l][ll][tr];
          d_op_d_qk = derivative*d_op_d_profile;
          d_em_d_qk = derivative*d_em_d_profile;
          // + Background sources
          d_op_d_qk += op_pert_lte[3][ll][x1l][x2l][ll];
          d_em_d_qk += em_pert_lte[3][ll][x1l][x2l][ll];
          remote_response = dR_dI * (response_to_op[l][ll] * d_op_d_qk + response_to_em[l][ll] * d_em_d_qk);
          beta_v_micro[ll][(l-x3l) * nmap +rmap[z][i] +1] -= remote_response/op_ref;
          beta_v_micro[ll][(l-x3l) * nmap +rmap[z][ii] +1] += remote_response/op_ref;

          // Then, an additional term if we are using taugrid as a primary grid:
          if (istaugrid){      
            fp_t op_ref_der = parent_atm->get_op_referent_der(3,ll,x1l,x2l,ll);
            beta_v_micro[ll][(l-x3l) * nmap +rmap[z][i] +1] += dR_dI *  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + 
              response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            beta_v_micro[ll][(l-x3l) * nmap +rmap[z][ii] +1] -= dR_dI *  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + 
              response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
          }

          // Now we need to sort out left hand side:
          fp_t der_i = line_energy * 0.25 / pi * current_profile[x1l][x2l][ll][tr] * (A[z][i][ii] * response_to_em[l][ll] - B[z][i][ii] * response_to_op[l][ll]); 
          fp_t der_ii = line_energy * 0.25 / pi * current_profile[x1l][x2l][ll][tr] * (B[z][ii][i] * response_to_op[l][ll]);
          der_i /= op_ref;
          der_ii /= op_ref;

          if (l==ll){
            der_i -= (1.0 - line_factor) / opp[x1l][x2l][l] * B[z][i][ii] * I[x1l][x2l][l] * line_energy * 0.25 / pi * current_profile[x1l][x2l][l][tr] / line_factor;
            der_ii += (1.0 - line_factor) / opp[x1l][x2l][l] * B[z][ii][i] * I[x1l][x2l][l] * line_energy * 0.25 / pi * current_profile[x1l][x2l][l][tr] / line_factor;
          }
              
          fp_t nd = der_ii * dR_dI;
          fp_t dd = der_i * dR_dI;

          response_matrix[(l-x3l) * nmap + rmap[z][i]+1][(ll-x3l) * nmap + rmap[z][ii]+1] += nd;
          response_matrix[(l-x3l) * nmap + rmap[z][i]+1][(ll-x3l) * nmap + rmap[z][i] +1] += dd;
          response_matrix[(l-x3l) * nmap + rmap[z][ii]+1][(ll-x3l) * nmap + rmap[z][i] +1] -= dd;
          response_matrix[(l-x3l) * nmap + rmap[z][ii]+1][(ll-x3l) * nmap + rmap[z][ii]+1] -= nd;

          for (int iii=0;iii<nmap;++iii){ // Now note that this notation here is using >>map<<
            response_matrix[(l-x3l) * nmap + rmap[z][i]+1][(ll-x3l) * nmap + iii+1] += dR_dI / op_ref * 
              (response_to_op[l][ll] * d_op_d_ni[ll][iii] + response_to_em[l][ll] * d_em_d_ni[ll][iii]);
            response_matrix[(l-x3l) * nmap + rmap[z][ii]+1][(ll-x3l) * nmap + iii+1] -= dR_dI / op_ref * 
              (response_to_op[l][ll] * d_op_d_ni[ll][iii] + response_to_em[l][ll] * d_em_d_ni[ll][iii]);
            
            if (l==ll){
              response_matrix[(l-x3l) * nmap + i+1][(ll-x3l) * nmap + iii+1] -= J_chunk * d_op_d_ni[l][iii] / opp[x1l][x2l][l];
              response_matrix[(l-x3l) * nmap + ii+1][(ll-x3l) * nmap + iii+1] += J_chunk * d_op_d_ni[l][iii] / opp[x1l][x2l][l];         
            }
          }
          
        }
      }

      else if (type == 1){ // Else it is the f-b transition (i->ii):

        // Keep in mind that ii -> lower level 
        // i is equal to nl[z], z is z corresponding to ii

        fp_t T = fetch_temperature(x1l,x2l,l);
        fp_t n_e = fetch_Ne(x1l,x2l,l);
        fp_t ddR_i_cont_dI = dR_i_cont_dI(z, ii, T, n_e, lambda);
        fp_t ddR_cont_i_dI = dR_cont_i_dI(z, ii, T, n_e, lambda);

        // Now, explicit dependency of radiative rates on atmospheric parameters. Only applies of local:
        fp_t * ddR_i_cont_dqk = dR_i_cont_dqk(z,ii,x1l,x2l,l,T,n_e,I[x1l][x2l][l],lambda);
        fp_t * ddR_cont_i_dqk = dR_cont_i_dqk(z,ii,x1l,x2l,l,T,n_e,I[x1l][x2l][l],lambda);

        beta_Temp[l][(l-x3l)*nmap+rmap[z+1][0]+1] += pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dqk[1]*angular_weight*lambda_w;
        beta_Temp[l][(l-x3l)*nmap+rmap[z][ii]+1] -= pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dqk[1]*angular_weight*lambda_w;
        beta_Temp[l][(l-x3l)*nmap+rmap[z+1][0]+1] -= pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dqk[1]*angular_weight*lambda_w;
        beta_Temp[l][(l-x3l)*nmap+rmap[z][ii]+1] += pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dqk[1]*angular_weight*lambda_w;
        
        beta_density[l][(l-x3l)*nmap+rmap[z+1][0]+1] += pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dqk[2]*angular_weight*lambda_w;
        beta_density[l][(l-x3l)*nmap+rmap[z][ii]+1] -= pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dqk[2]*angular_weight*lambda_w;
        beta_density[l][(l-x3l)*nmap+rmap[z+1][0]+1] -= pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dqk[2]*angular_weight*lambda_w;
        beta_density[l][(l-x3l)*nmap+rmap[z][ii]+1] += pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dqk[2]*angular_weight*lambda_w;

        delete [](ddR_i_cont_dqk+1);
        delete [](ddR_cont_i_dqk+1);
        
        for (int ll=x3l;ll<=x3h;++ll){

          fp_t op_ref = 1.0;
          if (istaugrid)
            op_ref = parent_atm->get_op_referent(x1l,x2l,ll);

          // Separate known perturbations to intensity. These come from known perturbations of chi and eta
          fp_t * dI_dqk_lte = new fp_t [7]-1;
          memset(dI_dqk_lte+1,0,7*sizeof(fp_t));
          dI_dqk_lte[1] = response_to_op[l][ll] * op_pert_lte[1][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[1][ll][x1l][x2l][ll];
          dI_dqk_lte[2] = response_to_op[l][ll] * op_pert_lte[2][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[2][ll][x1l][x2l][ll];

          // Additional known perturbations to the intensity come from known perturbations to opacity and emissivity. Here we go step-by-step
          // so we will add first the contributions from known perturbations to opacity and emissivity
      
          dI_dqk_lte[1] += response_to_op[l][ll] * bb_op_derivative_ex[ll][1] + response_to_em[l][ll] * bb_em_derivative_ex[ll][1];
          dI_dqk_lte[2] += response_to_op[l][ll] * bb_op_derivative_ex[ll][2] + response_to_em[l][ll] * bb_em_derivative_ex[ll][2];
          dI_dqk_lte[1] += response_to_op[l][ll] * bf_op_derivative_ex[ll][1] + response_to_em[l][ll] * bf_em_derivative_ex[ll][1];
          dI_dqk_lte[2] += response_to_op[l][ll] * bf_op_derivative_ex[ll][2] + response_to_em[l][ll] * bf_em_derivative_ex[ll][2];

          // This is contribution of known Intensity perturbations to the free-bound rates
          beta_Temp[ll][(l-x3l)*nmap+rmap[z+1][0]+1] += pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*dI_dqk_lte[1]*angular_weight*lambda_w/op_ref;
          beta_Temp[ll][(l-x3l)*nmap+rmap[z][ii]+1] -= pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*dI_dqk_lte[1]*angular_weight*lambda_w/op_ref;
          beta_density[ll][(l-x3l)*nmap+rmap[z+1][0]+1] += pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*dI_dqk_lte[2]*angular_weight*lambda_w/op_ref;
          beta_density[ll][(l-x3l)*nmap+rmap[z][ii]+1] -= pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*dI_dqk_lte[2]*angular_weight*lambda_w/op_ref;

          // This is the contribution of known intensity perturbations to the bound-free rates
          beta_Temp[ll][(l-x3l)*nmap+rmap[z+1][0]+1] -= pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dI*dI_dqk_lte[1]*angular_weight*lambda_w/op_ref;
          beta_Temp[ll][(l-x3l)*nmap+rmap[z][ii]+1] += pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dI*dI_dqk_lte[1]*angular_weight*lambda_w/op_ref;
          beta_density[ll][(l-x3l)*nmap+rmap[z+1][0]+1] -= pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dI*dI_dqk_lte[2]*angular_weight*lambda_w/op_ref;
          beta_density[ll][(l-x3l)*nmap+rmap[z][ii]+1] += pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dI*dI_dqk_lte[2]*angular_weight*lambda_w/op_ref;

          if (istaugrid){ // if tau grid

            fp_t op_ref_der = parent_atm->get_op_referent_der(1,ll,x1l,x2l,ll);
            beta_Temp[ll][(l-x3l)*nmap+rmap[z+1][0]+1] -= pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*angular_weight*lambda_w* 
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            beta_Temp[ll][(l-x3l)*nmap+rmap[z][ii]+1] += pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*angular_weight*lambda_w* 
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            
            beta_Temp[ll][(l-x3l)*nmap+rmap[z+1][0]+1] += pop[x1l][x2l][l].n[z][ii] * ddR_i_cont_dI * angular_weight * lambda_w *
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            beta_Temp[ll][(l-x3l)*nmap+rmap[z][ii]+1] -= pop[x1l][x2l][l].n[z][ii] * ddR_i_cont_dI * angular_weight * lambda_w *
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;

            op_ref_der = parent_atm->get_op_referent_der(2,ll,x1l,x2l,ll);
            beta_density[ll][(l-x3l)*nmap+rmap[z+1][0]+1] -= pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*angular_weight*lambda_w* 
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            beta_density[ll][(l-x3l)*nmap+rmap[z][ii]+1] += pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI*angular_weight*lambda_w* 
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            
            beta_density[ll][(l-x3l)*nmap+rmap[z+1][0]+1] += pop[x1l][x2l][l].n[z][ii] * ddR_i_cont_dI * angular_weight * lambda_w *
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
            beta_density[ll][(l-x3l)*nmap+rmap[z][ii]+1] -= pop[x1l][x2l][l].n[z][ii] * ddR_i_cont_dI * angular_weight * lambda_w *
              (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
          }

          delete [](dI_dqk_lte+1);
                
          // And finally, modification of the coupling matrix because of the contribution of all the level populations of other levels to 
          // this transition:
          for (int iii=0;iii<nmap;++iii){ // iii is now notation in rmap
            fp_t d_op_d_n = d_op_d_ni_bb[ll][iii] + d_op_d_ni[ll][iii];
            fp_t d_em_d_n = d_em_d_ni_bb[ll][iii] + d_em_d_ni[ll][iii];

            response_matrix[(l-x3l)*nmap+rmap[z+1][0]+1][(ll-x3l)*nmap+iii+1] += pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dI* 
              (response_to_op[l][ll]*d_op_d_n + response_to_em[l][ll]*d_em_d_n)*lambda_w*angular_weight/op_ref;
            response_matrix[(l-x3l)*nmap+rmap[z][ii]+1][(ll-x3l)*nmap+iii+1] -= pop[x1l][x2l][l].n[z][ii]*ddR_i_cont_dI* 
              (response_to_op[l][ll]*d_op_d_n + response_to_em[l][ll]*d_em_d_n)*lambda_w*angular_weight/op_ref;

            response_matrix[(l-x3l)*nmap+rmap[z+1][0]+1][(ll-x3l)*nmap+iii+1] -= pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI* 
              (response_to_op[l][ll]*d_op_d_n + response_to_em[l][ll]*d_em_d_n)*lambda_w*angular_weight/op_ref;
            response_matrix[(l-x3l)*nmap+rmap[z][ii]+1][(ll-x3l)*nmap+iii+1] += pop[x1l][x2l][l].n[z+1][0]*ddR_cont_i_dI* 
              (response_to_op[l][ll]*d_op_d_n + response_to_em[l][ll]*d_em_d_n)*lambda_w*angular_weight/op_ref;
          }
        }

      }
  }
  del_ft4dim(profile_derivative_T,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  del_ft4dim(profile_derivative_density,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  del_ft4dim(profile_derivative_vt,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  del_ft4dim(profile_derivative_vr,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  for (int ll=x3l;ll<=x3h;++ll){
    delete [](bf_op_derivative_ex[ll]+1);
    delete [](bf_em_derivative_ex[ll]+1);
    delete [](bb_op_derivative_ex[ll]+1);
    delete [](bb_em_derivative_ex[ll]+1);
  }
  delete[](bf_op_derivative_ex+x3l);
  delete[](bf_em_derivative_ex+x3l);
  delete[](bb_op_derivative_ex+x3l);
  delete[](bb_em_derivative_ex+x3l);
  del_ft2dim(d_op_d_ni,x3l,x3h,0,nmap-1);
  del_ft2dim(d_em_d_ni,x3l,x3h,0,nmap-1);
  del_ft2dim(d_op_d_ni_bb,x3l,x3h,0,nmap-1);
  del_ft2dim(d_em_d_ni_bb,x3l,x3h,0,nmap-1);

  del_ft3dim(profile_derivatives,x3l,x3h,1,7,1,ntr);

  return 0; 
}

int atom::add_response_contributions(fp_t *** I, fp_t ** response_to_op, fp_t ** response_to_em, fp_t *** opp, fp_t *** em, fp_t lambda, fp_t lambda_w, fp_t theta, fp_t phi, fp_t angular_weight, 
  fp_t *** vlos, fp_t ***** op_pert_lte, fp_t ***** em_pert_lte){

  // This function requires some explanation. 
  // The same way we are handling the NLTE problem itself, we are going to try to handle the problem of finding response functions.
  // In principle you should compute all the quantities which are needed for the computation of radiative rates, but then also some other stuff.

  int ncmp = 1; // total number of Stokes components
  int istaugrid = parent_atm->is_tau_grid();

  // If this is the proper transition, i.e. if it has mean intensiy, approximate operator and `norm' 
  if (Jb && norm && NLTE){
    for(int x1i=x1l;x1i<=x1h;++x1i)
      for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i){

          // Fetch them in the more convenient arrays:
          fp_t *Jm = Jb[x1i][x2i][x3i];
          fp_t *Jn = Ju[x1i][x2i][x3i];
          fp_t *nrm = norm[x1i][x2i][x3i];

          // For each Stokes component and for each transition.

          for (int tr=1; tr<=ntr; ++tr){

            // Now, for the general case which involves the presence of the continuum and the overlapping lines, the approximate lambda operator for the transition 'tr'
            // is actually first multiplied by the ratio of the opacity due to the transition in question to the total opacity. We need to extract the transition.

            int z_state = inverse_tmap[tr][2];
            int lower_level = inverse_tmap[tr][3];
            int upper_level = inverse_tmap[tr][4]; 

            if (upper_level < nl[z_state]){ // Line transitions:

              fp_t line_energy = ee[z_state][upper_level] - ee[z_state][lower_level];

              fp_t line_opacity = (pop[x1i][x2i][x3i].n[z_state][lower_level] * B[z_state][lower_level][upper_level] - pop[x1i][x2i][x3i].n[z_state][upper_level] * B[z_state][upper_level][lower_level]) 
                * current_profile[x1i][x2i][x3i][tr] * line_energy * 0.25 / pi;
              fp_t line_factor = line_opacity / opp[x1i][x2i][x3i];
              //line_factor = 1.0;
              
              fp_t elementary_contribution = lambda_w * angular_weight * I[x1i][x2i][x3i] * current_profile[x1i][x2i][x3i][tr] * line_factor * 0.25 / pi;

              Jm[tr] += elementary_contribution;
              Jn[tr] += elementary_contribution;
            }

            // If it is a bound-free transition. There is no profile. Instead, compute rates directly, so:

            else if (upper_level == nl[z_state]){
              fp_t sigma = (bf[z_state][lower_level]) ? bf[z_state][lower_level]->U(lambda) : 0.0;
              Jm[tr] += lambda_w * angular_weight * sigma * lambda / h /c * I[x1i][x2i][x3i];
              Jn[tr] += lambda_w * angular_weight * sigma * lambda / h /c * (I[x1i][x2i][x3i] + 2.0 * h * c * c / pow(lambda, 5.0)) * exp(-h * c / lambda / k / fetch_temperature(x1i, x2i, x3i));
            }
          }
        }
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Before solving response matrix it is advisable to compute derivatives of profiles to relevant quantities as they are used multiple times in the code:
    fp_t **** profile_derivative_T;
    profile_derivative_T = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
    memset(profile_derivative_T[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
    fp_t **** profile_derivative_density;
    profile_derivative_density = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
    memset(profile_derivative_density[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
    fp_t **** profile_derivative_vt;
    profile_derivative_vt = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
    memset(profile_derivative_vt[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
    fp_t **** profile_derivative_vr;
    profile_derivative_vr = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
    memset(profile_derivative_vr[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
    
    for (int x3i=x3l;x3i<=x3h;++x3i)
      for (int tr = 1;tr<=ntr;++tr){
          int z_state = inverse_tmap[tr][2];
          int lower_level = inverse_tmap[tr][3];
          int upper_level = inverse_tmap[tr][4]; 
          if (upper_level < nl[z_state]){
            int z_i = z_state;
            int l_i = upper_level;
            int l_ii = lower_level;
            fp_t x,a; x=a=0;
            fp_t * x_der;
            fp_t * a_der;
            fp_t * dld_der;
            x_der = new fp_t [7]-1;
            a_der = new fp_t [7]-1;
            dld_der = new fp_t [7]-1;
            memset(x_der+1,0,7*sizeof(fp_t));
            memset(a_der+1,0,7*sizeof(fp_t));
            memset(dld_der+1,0,7*sizeof(fp_t));
            if (A[z_i][l_i][l_ii] > 1.0){
              compute_xa_der_scalar(x1l, x2l, x3i, z_i, l_i, l_ii, lambda, vlos, x, a, x_der, a_der, dld_der);
              fp_t F,H,dF,dH;
              fvoigtn(x,a,H,F,dH,dF);
              fp_t lam=fabs(h*c/(ee[z_i][l_i]-ee[z_i][l_ii]));
              fp_t dld = broad_dop(lam, fetch_temperature(x1l,x2l,x3i),fetch_vt(x1l,x2l,x3i));

              profile_derivative_T[x1l][x2l][x3i][tr] = dH * x_der[1] / dld - dF * a_der[1] / dld - H/dld/dld*dld_der[1];
              profile_derivative_density[x1l][x2l][x3i][tr] = dH * x_der[2] / dld - dF * a_der[2] / dld - H/dld/dld*dld_der[2];
              profile_derivative_vt[x1l][x2l][x3i][tr] = dH * x_der[3] / dld - dF * a_der[3] / dld - H/dld/dld*dld_der[3];
              profile_derivative_vr[x1l][x2l][x3i][tr] = (dH * x_der[4] / dld - dF * a_der[4] / dld - H/dld/dld*dld_der[4])*0.0*(cos(theta));
            }
            delete [](x_der+1);
            delete [](a_der+1);
            delete [](dld_der+1);
          }
    }
   
    for (int l = x3l; l<=x3h; ++l){
      for (int ll = x3l; ll <= x3h; ++ll){

          for (int i=0; i<nmap; ++i)
            for (int ii = 0; ii<nmap; ++ii){

              // Find the upper and lower level and the z:
              // This is super complicated:
              int z_i = zmap[i];
              int l_i = lmap[i];
              int z_ii = zmap[ii];
              int l_ii = lmap[ii];

              int is_bf = 0; // 0 if it is not, 1 if it is bound->free and -1 if it is free->bound

              if (z_i == z_ii-1 && l_ii == 0)
                if (lambda <= bf[z_i][l_i]->getlambda_crit()) // We only care if it can actually ionize
                  is_bf = 1; // bound-free
              else if (z_i == z_ii+1 && l_i == 0)
                if (lambda <= bf[z_ii][l_ii]->getlambda_crit()) // Same here
                  is_bf = -1; // free-bound

              fp_t op_ref = 1.0;
              if (parent_atm->is_tau_grid())
                op_ref = parent_atm->get_op_referent(x1l,x2l,ll);

              // If it is a b-b transition (line):
              if ((z_i == z_ii) && (l_i != l_ii)){

                int tr = tmap[z_i][l_i][l_ii];
                int z_state = inverse_tmap[tr][2];
                int lower_level = inverse_tmap[tr][3];
                int upper_level = inverse_tmap[tr][4];
                if (A[z_state][upper_level][lower_level] > 1.0){ 

                // Lets first write it in the column where upper is i and lower is ii
            
                fp_t line_energy = ee[z_state][upper_level] - ee[z_state][lower_level];
                fp_t lam = h * c / line_energy;

                fp_t line_opacity = (pop[x1l][x2l][l].n[z_state][lower_level] * B[z_state][lower_level][upper_level] - pop[x1l][x2l][l].n[z_state][upper_level] * B[z_state][upper_level][lower_level]) 
                   * line_energy * 0.25 / pi;
                fp_t line_factor = line_opacity * current_profile[x1l][x2l][l][tr] / opp[x1l][x2l][l];
                //line_factor = 1.0;
                
                fp_t elementary_contribution = lambda_w * angular_weight * current_profile[x1l][x2l][l][tr] * line_factor * 0.25 / pi;

                fp_t der_l_i=0.0;
                fp_t der_l_ii=0.0;
                if(l_i > l_ii){
                  der_l_i = line_energy * 0.25 / pi * current_profile[x1l][x2l][ll][tr] * (A[z_i][l_i][l_ii] * response_to_em[l][ll] - B[z_i][l_i][l_ii] * response_to_op[l][ll]); 
                  der_l_ii = line_energy * 0.25 / pi * current_profile[x1l][x2l][ll][tr] * (B[z_i][l_ii][l_i] * response_to_op[l][ll]);
                  der_l_i /= op_ref;
                  der_l_ii /= op_ref;

                  if (l==ll){
                    der_l_i -= (1.0 - line_factor) / opp[x1l][x2l][l] * B[z_i][l_i][l_ii] * I[x1l][x2l][l] * line_energy * 0.25 / pi * current_profile[x1l][x2l][l][tr] / line_factor;
                    der_l_ii += (1.0 - line_factor) / opp[x1l][x2l][l] * B[z_i][l_ii][l_i] * I[x1l][x2l][l] * line_energy * 0.25 / pi * current_profile[x1l][x2l][l][tr] / line_factor;
                  }
	  						}
                else if (l_ii > l_i){
                  der_l_ii = line_energy * 0.25 / pi * current_profile[x1l][x2l][ll][tr] * (A[z_i][l_ii][l_i] * response_to_em[l][ll] - B[z_i][l_ii][l_i] * response_to_op[l][ll]); 
                  der_l_i = line_energy * 0.25 / pi * current_profile[x1l][x2l][ll][tr] * (B[z_i][l_i][l_ii] * response_to_op[l][ll]);
                  der_l_i /= op_ref;
                  der_l_ii /= op_ref;

                  if (l==ll){
                    der_l_i += (1.0 - line_factor) / opp[x1l][x2l][l] * B[z_i][l_i][l_ii] * I[x1l][x2l][l] * line_energy * 0.25 / pi * current_profile[x1l][x2l][l][tr] / line_factor;
                    der_l_ii -= (1.0 - line_factor) / opp[x1l][x2l][l] * B[z_i][l_ii][l_i] * I[x1l][x2l][l] * line_energy * 0.25 / pi * current_profile[x1l][x2l][l][tr] / line_factor;
                  }
                }

                fp_t nd = (pop[x1l][x2l][l].n[z_i][l_ii] * B[z_i][l_ii][l_i] - pop[x1l][x2l][l].n[z_i][l_i] * B[z_i][l_i][l_ii]) *  
                  der_l_ii * elementary_contribution;
                fp_t dd = (pop[x1l][x2l][l].n[z_i][l_ii] * B[z_i][l_ii][l_i] - pop[x1l][x2l][l].n[z_i][l_i] * B[z_i][l_i][l_ii]) *  
                  der_l_i * elementary_contribution;

                nd /= norm[x1l][x2l][l][tr];
                dd /= norm[x1l][x2l][l][tr];

             // -----------------------------------------------------------------------------------------------------------------------------------------------------------
                response_matrix[(l-x3l) * nmap + i+1][(ll-x3l) * nmap + ii+1] += nd;
                response_matrix[(l-x3l) * nmap + i+1][(ll-x3l) * nmap + i +1] += dd;

             // Now we also need to take into account b-f contributions to these b-b transitions:
                
                for (int iii=0;iii<nmap;++iii){
                  response_matrix[(l-x3l) * nmap + i+1][(ll-x3l) * nmap + iii+1] += 
                    (pop[x1l][x2l][l].n[z_i][l_ii] * B[z_i][l_ii][l_i] - pop[x1l][x2l][l].n[z_i][l_i] * B[z_i][l_i][l_ii]) * elementary_contribution * 
                    (response_to_op[l][ll] * op_derivative_to_level(ll,iii,lambda) + response_to_em[l][ll] * em_derivative_to_level(ll,iii,lambda)) / op_ref / norm[x1l][x2l][l][tr];
                  if (l==ll)
                    response_matrix[(l-x3l) * nmap + i+1][(ll-x3l) * nmap + iii+1] -= 
                      (pop[x1l][x2l][l].n[z_i][l_ii] * B[z_i][l_ii][l_i] - pop[x1l][x2l][l].n[z_i][l_i] * B[z_i][l_i][l_ii]) * elementary_contribution * 
                      op_derivative_to_level(ll,iii,lambda) / opp[x1l][x2l][l] * I[x1l][x2l][l] / norm[x1l][x2l][l][tr];
                }
                

              // Now the derivative of the line profile, which needs to be taken into account on r.h.s.
                if (l == ll){  

                  fp_t constant_factor = (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * 
                   lambda_w * angular_weight * 0.25 / pi;     
                  
                  // --------- TEMPERATURE ----------------------------------------------------------------------------------------------//

                  fp_t derivative = profile_derivative_T[x1l][x2l][l][tr];

                  fp_t derived_part = 0.0;

                  // Profile derivative
                  derived_part += I[x1l][x2l][l] * derivative * line_factor / norm[x1l][x2l][l][tr];
                  // Norm derivative
                  derived_part -= I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] * line_factor 
                    / norm[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * norm_derivative[1][x1l][x2l][l][tr];
                  // Profile derivative influencing line factor:
                  derived_part += I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * (1.0 - line_factor) * line_opacity / opp[x1l][x2l][l]
                    * derivative;
                  // Background opacity derivative influencing the line factor:
                  derived_part -= I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * line_factor * op_pert_lte[1][l][x1l][x2l][l] / opp[x1l][x2l][l];

                  beta_Temp[l][(l-x3l)*nmap + i+1] += constant_factor * derived_part;
                  
                  // ------------ DENSITY -----------------------------------------------------------------------------------------------

                  derivative = profile_derivative_density[x1l][x2l][l][tr];
                  
                  derived_part = 0.0;

                  // Profile derivative
                  derived_part += I[x1l][x2l][l] * derivative * line_factor / norm[x1l][x2l][l][tr];
                  // Norm derivative
                  derived_part -= I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] * line_factor 
                    / norm[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * norm_derivative[2][x1l][x2l][l][tr];
                  // Profile derivative influencing line factor:
                  derived_part += I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * (1.0 - line_factor) * line_opacity / opp[x1l][x2l][l]
                    * derivative;
                  // Background opacity derivative influencing the line factor:
                  derived_part -= I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * line_factor * op_pert_lte[2][l][x1l][x2l][l] / opp[x1l][x2l][l];

                  beta_density[l][(l-x3l)*nmap + i+1] += constant_factor * derived_part;

                  // ---------- MICROTURBULENT VELOCITY -----------------------------------------------------------------------------------

                  derivative = profile_derivative_vt[x1l][x2l][l][tr];

                  derived_part = 0.0;
                  derived_part += I[x1l][x2l][l] * derivative * line_factor / norm[x1l][x2l][l][tr];
                  derived_part -= I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] * line_factor 
                    / norm[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * norm_derivative[3][x1l][x2l][l][tr];
                  derived_part += I[x1l][x2l][l] * current_profile[x1l][x2l][l][tr] / norm[x1l][x2l][l][tr] * (1.0 - line_factor) * line_opacity / opp[x1l][x2l][l]
                    * derivative;
                  beta_v_micro[l][(l-x3l)*nmap + i+1] += constant_factor * derived_part;

                }
                // ========================================================================================================================

                // NON-LOCAL PART ON THE R.H.S:

                fp_t derivative = profile_derivative_T[x1l][x2l][ll][tr];

                // I am really running out of these names
                fp_t J_chunk = derivative * (response_to_op[l][ll] * (pop[x1l][x2l][ll].n[z_i][lower_level] * B[z_i][lower_level][upper_level] - pop[x1l][x2l][ll].n[z_i][upper_level] * B[z_i][upper_level][lower_level])   
                  + response_to_em[l][ll] * pop[x1l][x2l][ll].n[z_i][upper_level] * A[z_i][upper_level][lower_level]) * line_energy * 0.25 / pi * elementary_contribution 
                  / norm[x1l][x2l][l][tr];
                  
                J_chunk /= op_ref;

                beta_Temp[ll][(l-x3l) * nmap +i + 1] += (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * J_chunk;

                // Finally there is additional contribution from other species that we know:
                fp_t chunk = (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) 
                  * elementary_contribution * (response_to_op[l][ll] * op_pert_lte[1][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[1][ll][x1l][x2l][ll]) / op_ref/ norm[x1l][x2l][l][tr];
                beta_Temp[ll][(l-x3l) * nmap +i +1] += chunk;

                //if (isnan(em_pert_lte[1][ll][x1l][x2l][ll])){
                //  printf("AAAAAAAAAAA!\n");
                  //exit(1);
                //}
                /*if (isnan(chunk)){
                    printf("%e %e %e %e %e %e \n", elementary_contribution,response_to_op[l][ll],op_pert_lte[1][ll][x1l][x2l][ll],response_to_em[l][ll],em_pert_lte[1][ll][x1l][x2l][ll],norm[x1l][x2l][l][tr]);
                    printf("!3! %d %d %d \n", ll,l,i);
                    exit(1);
                  }*/
                // Then, an additional term if we are using taugrid as a primary grid:
                if (parent_atm->is_tau_grid()){
                  
                  fp_t op_ref_der = parent_atm->get_op_referent_der(1,ll,x1l,x2l,ll);
                  
                  beta_Temp[ll][(l-x3l) * nmap +i +1] -= (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) 
                  * elementary_contribution * (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] 
                  + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref / norm[x1l][x2l][l][tr];;
                }
                // ------------------ DENSITY -----------------------------------------------------------------------------------------------------------------
                fp_t derivative_d = profile_derivative_density[x1l][x2l][ll][tr];

                J_chunk = derivative_d * (response_to_op[l][ll] * (pop[x1l][x2l][ll].n[z_i][lower_level] * B[z_i][lower_level][upper_level] - pop[x1l][x2l][ll].n[z_i][upper_level] * B[z_i][upper_level][lower_level])   
                  + response_to_em[l][ll] * pop[x1l][x2l][ll].n[z_i][upper_level] * A[z_i][upper_level][lower_level]) * line_energy * 0.25 / pi * elementary_contribution 
                  / norm[x1l][x2l][l][tr];

                J_chunk /= op_ref;

                beta_density[ll][(l-x3l) * nmap +i + 1] += (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * J_chunk;

                // Finally there is additional contribution from other spieces that we know:
                beta_density[ll][(l-x3l) * nmap +i +1] += (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) 
                  * elementary_contribution * (response_to_op[l][ll] * op_pert_lte[2][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[2][ll][x1l][x2l][ll]) / op_ref / norm[x1l][x2l][l][tr];;
                
                if (parent_atm->is_tau_grid()){
                  fp_t op_ref_der = parent_atm->get_op_referent_der(2,ll,x1l,x2l,ll);
                  beta_density[ll][(l-x3l) * nmap +i +1] -= (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) 
                  * elementary_contribution * (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] 
                  + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref / norm[x1l][x2l][l][tr];
                }

                // ---------------- MICROTURBULENT VELOCITY ---------------------------------------------------------------------------------------------------

                fp_t derivative_vt = profile_derivative_vt[x1l][x2l][ll][tr];

                J_chunk = derivative_vt * (response_to_op[l][ll] * (pop[x1l][x2l][ll].n[z_i][lower_level] * B[z_i][lower_level][upper_level] - pop[x1l][x2l][ll].n[z_i][upper_level] * B[z_i][upper_level][lower_level])   
                  + response_to_em[l][ll] * pop[x1l][x2l][ll].n[z_i][upper_level] * A[z_i][upper_level][lower_level]) * line_energy * 0.25 / pi * elementary_contribution 
                  / norm[x1l][x2l][l][tr];

                J_chunk /= op_ref;

                beta_v_micro[ll][(l-x3l) * nmap +i + 1] += (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * J_chunk;

                // Some additional, hopefully clean(r)er formulations:
                
                fp_t * bf_op_derivative_ex = op_derivative_explicit(ll,0,0,lambda);
                fp_t * bf_em_derivative_ex = em_derivative_explicit(ll,0,0,lambda);

                beta_Temp[ll][(l-x3l) * nmap +i + 1] += (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * 
                  (response_to_op[l][ll] * bf_op_derivative_ex[1] + response_to_em[l][ll] * bf_em_derivative_ex[1]) * elementary_contribution / op_ref / norm[x1l][x2l][l][tr];
                beta_density[ll][(l-x3l) * nmap +i + 1] += (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * 
                  (response_to_op[l][ll] * bf_op_derivative_ex[2] + response_to_em[l][ll] * bf_em_derivative_ex[2]) * elementary_contribution / op_ref / norm[x1l][x2l][l][tr];

                // And additional contributions from the 'line factor opacity'
                if (l==ll){
                  beta_Temp[ll][(l-x3l) * nmap +i + 1] -= (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * 
                  bf_op_derivative_ex[1] / opp[x1l][x2l][l] * I[x1l][x2l][l] / norm[x1l][x2l][l][tr] * elementary_contribution;
                  beta_density[ll][(l-x3l) * nmap +i + 1] -= (B[z_i][l_i][l_ii] * pop[x1l][x2l][l].n[z_i][l_i] - B[z_i][l_ii][l_i] * pop[x1l][x2l][l].n[z_i][l_ii]) * 
                  bf_op_derivative_ex[2] / opp[x1l][x2l][l] * I[x1l][x2l][l] / norm[x1l][x2l][l][tr] * elementary_contribution;
                }
                delete[](bf_op_derivative_ex+1);
                delete[](bf_em_derivative_ex+1);
              }  
            }

            // Otherwise it has to be a b-f transition, f-b transition:
            else if (is_bf == 1){ // Such elegant very wow
        
              // Separate known perturbations to intensity. These come from known perturbations of chi and eta
              fp_t * dI_dqk_lte = new fp_t [7]-1;
              memset(dI_dqk_lte+1,0,7*sizeof(fp_t));
              dI_dqk_lte[1] = response_to_op[l][ll] * op_pert_lte[1][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[1][ll][x1l][x2l][ll];
              dI_dqk_lte[2] = response_to_op[l][ll] * op_pert_lte[2][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[2][ll][x1l][x2l][ll];

              // Additional known perturbations to the intensity come from known perturbations to opacity and emissivity. Here we go step-by-step
              // so we will add first the contributions from known perturbations to opacity and emissivity
              fp_t ** profile_derivatives = ft2dim(1,7,1,ntr);
              memset(profile_derivatives[1]+1,0,7*ntr*sizeof(fp_t));
              for (int tr=1;tr<=ntr;++tr){
                profile_derivatives[1][tr] = profile_derivative_T[x1l][x2l][ll][tr];
                profile_derivatives[2][tr] = profile_derivative_density[x1l][x2l][ll][tr];
                profile_derivatives[3][tr] = profile_derivative_vt[x1l][x2l][ll][tr];
              }

              fp_t * bb_op_derivative_ex = bb_op_derivative_explicit(x1l,x2l,ll,lambda,profile_derivatives);
              fp_t * bb_em_derivative_ex = bb_em_derivative_explicit(x1l,x2l,ll,lambda,profile_derivatives);
              fp_t * bf_op_derivative_ex = op_derivative_explicit(ll,0,0,lambda);
              fp_t * bf_em_derivative_ex = em_derivative_explicit(ll,0,0,lambda);

              dI_dqk_lte[1] += response_to_op[l][ll] * bb_op_derivative_ex[1] + response_to_em[l][ll] * bb_em_derivative_ex[1];
              dI_dqk_lte[2] += response_to_op[l][ll] * bb_op_derivative_ex[2] + response_to_em[l][ll] * bb_em_derivative_ex[2];
              dI_dqk_lte[1] += response_to_op[l][ll] * bf_op_derivative_ex[1] + response_to_em[l][ll] * bf_em_derivative_ex[1];
              dI_dqk_lte[2] += response_to_op[l][ll] * bf_op_derivative_ex[2] + response_to_em[l][ll] * bf_em_derivative_ex[2];

              del_ft2dim(profile_derivatives,1,7,1,ntr);
              delete[](bb_op_derivative_ex+1);
              delete[](bb_em_derivative_ex+1);
              delete[](bf_op_derivative_ex+1);
              delete[](bf_em_derivative_ex+1);

              fp_t T = fetch_temperature(x1l,x2l,l);
              fp_t n_e = fetch_Ne(x1l,x2l,l);

              fp_t ddR_i_cont_dI = dR_i_cont_dI(z_i, l_i, T, n_e, lambda);
              fp_t ddR_cont_i_dI = dR_cont_i_dI(z_i, l_i, T, n_e, lambda);

              // This is contribution of known Intensity perturbations to the bound-free rates: 
              beta_Temp[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_i_cont_dI  * dI_dqk_lte[1] * angular_weight * lambda_w/op_ref;
              beta_density[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_i_cont_dI  * dI_dqk_lte[2] * angular_weight * lambda_w/op_ref;

              // This is the contribution of known intensity perturbations to the free-bound rates, this comes with the opposite sign from the above:
              beta_Temp[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_cont_i_dI * dI_dqk_lte[1] * angular_weight * lambda_w/op_ref;
              beta_density[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_cont_i_dI * dI_dqk_lte[2] * angular_weight * lambda_w/op_ref;

              if (parent_atm->is_tau_grid()){ // if tau grid

                // First sort out temperature
                fp_t op_ref_der = parent_atm->get_op_referent_der(1,ll,x1l,x2l,ll);
                beta_Temp[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_i][l_i] * ddR_i_cont_dI * angular_weight * lambda_w * 
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
                beta_Temp[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_cont_i_dI * angular_weight * lambda_w *
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
                // And then the good old density
                op_ref_der = parent_atm->get_op_referent_der(2,ll,x1l,x2l,ll);
                beta_density[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_i][l_i] * ddR_i_cont_dI * angular_weight * lambda_w * 
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
                beta_density[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_cont_i_dI * angular_weight * lambda_w *
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
              }

              delete [](dI_dqk_lte+1);
                
              // Now, explicit dependency of radiative rates on atmospheric parameters. Only applies of local:
              if (l==ll){   
                fp_t * ddR_i_cont_dqk = dR_i_cont_dqk(z_i,l_i,x1l,x2l,l,T,n_e,I[x1l][x2l][l],lambda);
                fp_t * ddR_cont_i_dqk = dR_cont_i_dqk(z_i,l_i,x1l,x2l,l,T,n_e,I[x1l][x2l][l],lambda);

                beta_Temp[l][(l-x3l)*nmap+i+1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_i_cont_dqk[1] * angular_weight * lambda_w;
                beta_Temp[l][(l-x3l)*nmap+i+1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_cont_i_dqk[1] * angular_weight * lambda_w;
                beta_density[l][(l-x3l)*nmap+i+1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_i_cont_dqk[2] * angular_weight * lambda_w;
                beta_density[l][(l-x3l)*nmap+i+1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_cont_i_dqk[2] * angular_weight * lambda_w;

                delete [](ddR_i_cont_dqk+1);
                delete [](ddR_cont_i_dqk+1);
              }

              // And finally, modification of the coupling matrix because of the contribution of all the level populations of other levels to 
              // this transition:

              for (int iii=0;iii<nmap;++iii){
                fp_t d_op_d_n = bb_op_derivative_to_level(x1l,x2l,ll,iii,lambda);
                fp_t d_em_d_n = bb_em_derivative_to_level(x1l,x2l,ll,iii,lambda);

                d_op_d_n += op_derivative_to_level(ll,iii,lambda); 
                d_em_d_n += em_derivative_to_level(ll,iii,lambda); 

                response_matrix[(l-x3l)*nmap+i+1][(ll-x3l)*nmap+iii+1] += pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_cont_i_dI * 
                  (response_to_op[l][ll] * d_op_d_n + response_to_em[l][ll] * d_em_d_n) * lambda_w * angular_weight/op_ref;

                response_matrix[(l-x3l)*nmap+i+1][(ll-x3l)*nmap+iii+1] -= pop[x1l][x2l][l].n[z_i][l_i] * ddR_i_cont_dI * 
                  (response_to_op[l][ll] * d_op_d_n + response_to_em[l][ll] * d_em_d_n) * lambda_w * angular_weight/op_ref;
              }
            }
            else if (is_bf == -1){ // Even more elegant even more wow

             // This is all practically the same as in the case of b-f transitions, except that we have to take the opposite level
             // That is l_i is now continuum while l_ii is bound.
             // Separate known perturbations to intensity. These come from known perturbations of chi and eta.
              fp_t * dI_dqk_lte = new fp_t [7]-1;
              memset(dI_dqk_lte+1,0,7*sizeof(fp_t));
              dI_dqk_lte[1] = response_to_op[l][ll] * op_pert_lte[1][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[1][ll][x1l][x2l][ll];
              dI_dqk_lte[2] = response_to_op[l][ll] * op_pert_lte[2][ll][x1l][x2l][ll] + response_to_em[l][ll] * em_pert_lte[2][ll][x1l][x2l][ll];

              // Additional known perturbations to the intensity come from known perturbations to opacity and emissivity. Here we go step-by-step (ooo baby)
              // so we will add first the contributions from known perturbations to opacity and emissivity
              fp_t ** profile_derivatives = ft2dim(1,7,1,ntr);
              memset(profile_derivatives[1]+1,0,7*ntr*sizeof(fp_t));
              for (int tr=1;tr<=ntr;++tr){
                profile_derivatives[1][tr] = profile_derivative_T[x1l][x2l][ll][tr];
                profile_derivatives[2][tr] = profile_derivative_density[x1l][x2l][ll][tr];
                profile_derivatives[3][tr] = profile_derivative_vt[x1l][x2l][ll][tr];
              }

              fp_t * bb_op_derivative_ex = bb_op_derivative_explicit(x1l,x2l,ll,lambda,profile_derivatives);
              fp_t * bb_em_derivative_ex = bb_em_derivative_explicit(x1l,x2l,ll,lambda,profile_derivatives);
              fp_t * bf_op_derivative_ex = op_derivative_explicit(ll,0,0,lambda);
              fp_t * bf_em_derivative_ex = em_derivative_explicit(ll,0,0,lambda);

              dI_dqk_lte[1] += response_to_op[l][ll] * bb_op_derivative_ex[1] + response_to_em[l][ll] * bb_em_derivative_ex[1];
              dI_dqk_lte[2] += response_to_op[l][ll] * bb_op_derivative_ex[2] + response_to_em[l][ll] * bb_em_derivative_ex[2];
              dI_dqk_lte[1] += response_to_op[l][ll] * bf_op_derivative_ex[1] + response_to_em[l][ll] * bf_em_derivative_ex[1];
              dI_dqk_lte[2] += response_to_op[l][ll] * bf_op_derivative_ex[2] + response_to_em[l][ll] * bf_em_derivative_ex[2];

              del_ft2dim(profile_derivatives,1,7,1,ntr);
              delete[](bb_op_derivative_ex+1);
              delete[](bb_em_derivative_ex+1);
              delete[](bf_op_derivative_ex+1);
              delete[](bf_em_derivative_ex+1);

      
              fp_t T = fetch_temperature(x1l,x2l,l);
              fp_t n_e = fetch_Ne(x1l,x2l,l);

              fp_t ddR_i_cont_dI = dR_i_cont_dI(z_ii, l_ii, T, n_e, lambda);
              fp_t ddR_cont_i_dI = dR_cont_i_dI(z_ii, l_ii, T, n_e, lambda);

              // This is contribution of known Intensity perturbations to the free-bound rates
              beta_Temp[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_cont_i_dI  * dI_dqk_lte[1] * angular_weight * lambda_w/op_ref;
              beta_density[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_cont_i_dI  * dI_dqk_lte[2] * angular_weight * lambda_w/op_ref;

              // This is the contribution of known intensity perturbations to the bound-free rates
              beta_Temp[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_i_cont_dI * dI_dqk_lte[1] * angular_weight * lambda_w/op_ref;
              beta_density[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_i_cont_dI * dI_dqk_lte[2] * angular_weight * lambda_w/op_ref;

              if (parent_atm->is_tau_grid()){ // if tau grid

                op_ref = parent_atm->get_op_referent(x1l,x2l,ll);
                // First sort out temperature
                fp_t op_ref_der = parent_atm->get_op_referent_der(1,ll,x1l,x2l,ll);
                beta_Temp[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_i][l_i] * ddR_cont_i_dI * angular_weight * lambda_w * 
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
                beta_Temp[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_i_cont_dI * angular_weight * lambda_w *
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
                // And then the good old density
                op_ref_der = parent_atm->get_op_referent_der(2,ll,x1l,x2l,ll);
                beta_density[ll][(l-x3l) * nmap +i + 1] -= pop[x1l][x2l][l].n[z_i][l_i] * ddR_cont_i_dI * angular_weight * lambda_w * 
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
                beta_density[ll][(l-x3l) * nmap +i + 1] += pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_i_cont_dI * angular_weight * lambda_w *
                  (response_to_op[l][ll] * op_ref_der * opp[x1l][x2l][ll] + response_to_em[l][ll] * op_ref_der * em[x1l][x2l][ll]) / op_ref/op_ref;
              }

              delete [](dI_dqk_lte+1);
                
              // Now, explicit dependency of radiative rates on atmospheric parameters. Only applies of local:
              if (l==ll){   
                fp_t * ddR_i_cont_dqk = dR_i_cont_dqk(z_ii,l_ii,x1l,x2l,l,T,n_e,I[x1l][x2l][l],lambda);
                fp_t * ddR_cont_i_dqk = dR_cont_i_dqk(z_ii,l_ii,x1l,x2l,l,T,n_e,I[x1l][x2l][l],lambda);

                beta_Temp[l][(l-x3l)*nmap+i+1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_cont_i_dqk[1] * angular_weight * lambda_w;
                beta_Temp[l][(l-x3l)*nmap+i+1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_i_cont_dqk[1] * angular_weight * lambda_w;
                beta_density[l][(l-x3l)*nmap+i+1] += pop[x1l][x2l][l].n[z_i][l_i] * ddR_cont_i_dqk[2] * angular_weight * lambda_w;
                beta_density[l][(l-x3l)*nmap+i+1] -= pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_i_cont_dqk[2] * angular_weight * lambda_w;

                delete [](ddR_i_cont_dqk+1);
                delete [](ddR_cont_i_dqk+1);
              }

              // And finally, modification of the coupling matrix because of the contribution of all the level populations of other levels to 
              // this transition:

              for (int iii=0;iii<nmap;++iii){
                fp_t d_op_d_n = bb_op_derivative_to_level(x1l,x2l,ll,iii,lambda);
                fp_t d_em_d_n = bb_em_derivative_to_level(x1l,x2l,ll,iii,lambda);

                d_op_d_n += op_derivative_to_level(ll,iii,lambda); 
                d_em_d_n += em_derivative_to_level(ll,iii,lambda); 

                response_matrix[(l-x3l)*nmap+i+1][(ll-x3l)*nmap+iii+1] += pop[x1l][x2l][l].n[z_ii][l_ii] * ddR_i_cont_dI * 
                  (response_to_op[l][ll] * d_op_d_n + response_to_em[l][ll] * d_em_d_n) * lambda_w * angular_weight/op_ref;

                response_matrix[(l-x3l)*nmap+i+1][(ll-x3l)*nmap+iii+1] -= pop[x1l][x2l][l].n[z_i][l_i] * ddR_cont_i_dI * 
                  (response_to_op[l][ll] * d_op_d_n + response_to_em[l][ll] * d_em_d_n) * lambda_w * angular_weight/op_ref;
              }
            }
          }
        }
                
    } 
    del_ft4dim(profile_derivative_T,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr); 
    del_ft4dim(profile_derivative_density,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr); 
    del_ft4dim(profile_derivative_vt,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr); 
    del_ft4dim(profile_derivative_vr,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);          
  }
  return 0;
}

// Same version but now written to work @ taugrid: REDUNDANT NOW

int atom::add_response_contributions_taugrid(fp_t *** I, fp_t *** op_ref, fp_t ***** op_ref_der, fp_t ** response_to_op, fp_t ** response_to_em, fp_t *** opp, fp_t *** em, fp_t lambda, fp_t lambda_w, fp_t theta, fp_t phi, fp_t angular_weight, 
  fp_t *** vlos, fp_t ***** op_pert_lte, fp_t ***** em_pert_lte){};

int atom::add_pops_to_response(int depth_of_perturbation, int parameter_no){

  if(ntr){

    for (int x3i = x3l;x3i<=x3h; ++x3i)
      for (int i = 1; i<=nmap; ++i)
        level_responses[parameter_no][(x3i - x3l) * nmap + i][depth_of_perturbation] = pop[x1l][x2l][x3i].n[zmap[i-1]][lmap[i-1]]; 
  }

  return 0;
}
int atom::subtract_pops_from_response(int depth_of_perturbation, int parameter_no){

  if(ntr){

    for (int x3i = x3l;x3i<=x3h; ++x3i)
      for (int i = 1; i<=nmap; ++i)
        level_responses[parameter_no][(x3i - x3l) * nmap + i][depth_of_perturbation] -= pop[x1l][x2l][x3i].n[zmap[i-1]][lmap[i-1]]; 
  }

  return 0;
}
int atom::divide_responses_by_step(int depth_of_perturbation, int parameter_no,  fp_t step){

  if(ntr){
  for (int x3i = x3l;x3i<=x3h; ++x3i)
    for (int i = 1; i<=nmap; ++i)
      level_responses[parameter_no][(x3i - x3l) * nmap + i][depth_of_perturbation] /= step;
  }

  return 0;
}

int atom::print_population_responses(const char* filename, int depth_of_perturbation){

  FILE * output;
  output = fopen(filename, "w");

  for (int param = 1; param<=N_atm_parameters;++param)
    for (int x3i = x3l; x3i<=x3h; ++x3i){
      for (int i=1; i<=nmap; ++i)
        fprintf(output, "%15.15e ", level_responses[param][(x3i - x3l) * nmap + i][depth_of_perturbation]);// / pop[x1l][x2l][x3i].n[zmap[i-1]][lmap[i-1]]);
      fprintf(output, "\n");
    }

  fclose(output);
  return 0;
}

int atom::print_population_responses(const char* filename, int from, int to){

  if (!ntr){
    io.msg(IOL_INFO, "atom::responses : the atom has no levels. Moving on... \n");
    return 0;
  }

  FILE * output;
  output = fopen(filename, "w");

  for (int param = 1; param<=N_atm_parameters;++param)  
    for (int depth = from; depth <= to; ++ depth){
      for (int x3i = x3l; x3i<=x3h; ++x3i){
      //int x3i = depth;
        fprintf(output, "%d %d ", depth, x3i);
        for (int i=1; i<=nmap; ++i)
          fprintf(output, "%15.15e ", level_responses[param][(x3i - x3l) * nmap + i][depth]);// / pop[x1l][x2l][x3i].n[zmap[i-1]][lmap[i-1]]);
          fprintf(output, "\n");
        }
    }

  fclose(output);
  return 0;
}



fp_t atom::op_derivative_to_level(int depth, int i, fp_t lambda){

  // Basically here we can take into account all the processes which happen in the atom. b-b, b-f, f-f, Raighley, etc...
  // For the moment, we will only take into account b-f, but later we will restructure the code so it is used for all the 
  // transitions in the rate matrix. So, let's start:

  int z = zmap[i];
  int l = lmap[i];

  fp_t derivative = 0.0;
  fp_t sigma = (bf[z][l]) ? bf[z][l]->U(lambda) : 0.0;
  derivative += sigma;
  
  fp_t T = fetch_temperature(x1l,x2l,depth);

  if (z>0 && l == 0)
    for (int ll=0;ll<nl[z-1];++ll){
      fp_t sigmaa = (bf[z-1][ll]) ? bf[z-1][ll]->U(lambda) : 0.0;
      fp_t delta_E = ip[z-1] - ee[z-1][ll] - h *c / lambda;
      derivative -= sigmaa * fetch_Ne(x1l,x2l,depth) * saha_const * pow(T,-1.5) * fp_t(g[z-1][ll]) / fp_t(g[z][l]) * exp(delta_E / k /T);
    }
  return derivative;
}

fp_t atom::em_derivative_to_level(int depth, int i, fp_t lambda){
  int z = zmap[i];
  int l = lmap[i];

  fp_t derivative = 0;
  fp_t T = fetch_temperature(x1l,x2l,depth);

  if (z>0 && l == 0)
    for (int ll=0;ll<nl[z-1];++ll){
      fp_t sigmaa = (bf[z-1][ll]) ? bf[z-1][ll]->U(lambda) : 0.0;
      fp_t delta_E = ip[z-1] - ee[z-1][ll];
      derivative += sigmaa * fetch_Ne(x1l,x2l,depth) * saha_const * pow(T,-1.5) * fp_t(g[z-1][ll]) / fp_t(g[z][l]) * exp(delta_E / k /T) *
        (1.0 - exp(-h*c/lambda/k/T)) * Planck_f(lambda, T);
    }
  return derivative;
}


fp_t * atom::op_derivative_explicit(int depth, int i, int depthp, fp_t lambda){ // These return arrays, because we have derivatives w.r.t different atmospheric parameters
  
  // Now these are very similar to the ones above but here actually compute the known perturbation to the opacity. 
  // Again we will only involve the b-f opacity, for the moment.

  // At this interpretation i  and depth do not do anything! 

  fp_t * derivative = new fp_t [7]-1;
  memset(derivative+1,0,7*sizeof(fp_t));
  //return derivative;
  fp_t T = fetch_temperature(x1l, x2l, depth);
  fp_t ne = fetch_Ne(x1l,x2l,depth);
  
  for(int z=0;z<Z;++z) // skip final stage: it cannot be ionized
    for(int l=0;l<nl[z];++l){
      
      fp_t sigma = (bf[z][l]) ? bf[z][l]->U(lambda * (1.0)) : 0.0;
      if (sigma){
        // But we also need some corrections now.
        fp_t delta_E = (ip[z] - ee[z][l] - h*c/lambda);
        
        fp_t constant_part = -pop[x1l][x2l][depth].n[z+1][0] * saha_const * fp_t(g[z][l]) / fp_t(g[z+1][0]) * sigma;

        derivative[1] += constant_part * (parent_atm->get_ne_lte_derivative(1,x1l,x2l,depth) * pow(T,-1.5) * exp(delta_E/k/T));
        
        derivative[1] += constant_part * ne * (-1.5) * pow(T,-2.5) * exp(delta_E/k/T);
        derivative[1] += constant_part * ne * pow(T,-1.5) * exp(delta_E/k/T) * (-delta_E/k/T/T);

        derivative[2] += constant_part * (parent_atm->get_ne_lte_derivative(2,x1l,x2l,depth) * pow(T,-1.5) * exp(delta_E/k/T));
      }
  }
  return derivative;
}

fp_t * atom::em_derivative_explicit(int depth, int, int, fp_t lambda){

  // In this implementation, depth is enough 

  fp_t * derivative = new fp_t [7]-1;
  memset(derivative+1,0,7*sizeof(fp_t));
  //return derivative;
  fp_t T = fetch_temperature(x1l, x2l, depth);
  fp_t ne = fetch_Ne(x1l,x2l,depth);
  
  for(int z=0;z<Z;++z) // skip final stage: it cannot be ionized
    for(int l=0;l<nl[z];++l){
      
      fp_t sigma = (bf[z][l]) ? bf[z][l]->U(lambda * (1.0)) : 0.0;
      if (sigma){

        fp_t B_planck = Planck_f(lambda, T);

        // But we also need some corrections now.
        fp_t delta_E = (ip[z] - ee[z][l]);
        fp_t pop_mod = pop[x1l][x2l][depth].n[z+1][0] * ne * 2.07E-16 * pow(T, -1.5) * fp_t(g[z][l]) / fp_t(g[z+1][0]) * exp(delta_E / k /  T);
        
        derivative[1] += sigma *  (-pop_mod * exp(- h * c / lambda / k / T) * h*c/lambda/k/T/T) * B_planck;
        fp_t constant_part = pop[x1l][x2l][depth].n[z+1][0] * 2.07E-16 * fp_t(g[z][l]) / fp_t(g[z+1][0]) * B_planck * sigma * (1.0 - exp(-h*c/lambda/k/T));
        
        derivative[1] += constant_part * (parent_atm->get_ne_lte_derivative(1,x1l,x2l,depth) * pow(T,-1.5) * exp(delta_E/k/T) +
          ne * (-1.5 * pow(T,-2.5)) * exp(delta_E/k/T) + ne * pow(T,-1.5) * exp(delta_E/k/T) * -delta_E/k/T/T);
        // Finally, additional part because of Planck function:
        derivative[1] += pop_mod * sigma * (1.0 - exp(-h*c/lambda/k/T)) * Planck_f_derivative(lambda, T);
          
        derivative[2] += constant_part * (parent_atm->get_ne_lte_derivative(2,x1l,x2l,depth) * pow(T,-1.5) * exp(delta_E/k/T));
      }

  }
  return derivative;
  
}

fp_t atom::bb_op_derivative_to_level(int x1i, int x2i, int x3i, int i, fp_t lambda){

  int z = zmap[i];
  int l = lmap[i];
  
  fp_t derivative = 0.0;

  // We need all the transitions where this level is lower, it contributes through absorption:
  for (int ll=l+1;ll<nl[z];++ll)
    derivative += B[z][l][ll] * h * c / lambda * 0.25 / pi * current_profile[x1i][x2i][x3i][tmap[z][l][ll]];

  // But there it is higher, it contributes throught stimulated emission:
  for (int ll=0;ll<l;++ll)
    derivative -= B[z][l][ll] * h * c / lambda * 0.25 / pi * current_profile[x1i][x2i][x3i][tmap[z][l][ll]];

  return derivative;
}

fp_t atom::bb_em_derivative_to_level(int x1i, int x2i, int x3i, int i, fp_t lambda){

  int z = zmap[i];
  int l = lmap[i];
  
  fp_t derivative = 0.0;

  // For the emissions, only transitions which matter are the transitions where this level is the upper one: 

  for (int ll=0;ll<l;++ll)
    derivative += A[z][l][ll] * h * c / lambda * 0.25 / pi * current_profile[x1i][x2i][x3i][tmap[z][l][ll]];

  return derivative;
}
  

fp_t * atom::bb_op_derivative_explicit(int x1i, int x2i, int x3i, fp_t lambda, fp_t ** profile_derivatives){
  
  fp_t * derivative = new fp_t[7]-1;
  memset(derivative+1,0,7*sizeof(fp_t));

  for(int z=0;z<=Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t profile =  current_profile[x1i][x2i][x3i][tmap[z][i][ii]];
          
        derivative[1] += profile_derivatives[1][tmap[z][i][ii]] * (pop[x1i][x2i][x3i].n[z][ii] * B[z][ii][i] - pop[x1i][x2i][x3i].n[z][i] * B[z][i][ii]) * h * c / lam / 4.0 / pi;
        derivative[2] += profile_derivatives[2][tmap[z][i][ii]] * (pop[x1i][x2i][x3i].n[z][ii] * B[z][ii][i] - pop[x1i][x2i][x3i].n[z][i] * B[z][i][ii]) * h * c / lam / 4.0 / pi;
        derivative[3] += profile_derivatives[3][tmap[z][i][ii]] * (pop[x1i][x2i][x3i].n[z][ii] * B[z][ii][i] - pop[x1i][x2i][x3i].n[z][i] * B[z][i][ii]) * h * c / lam / 4.0 / pi;
  }       

  return derivative;
}

fp_t * atom::bb_em_derivative_explicit(int x1i, int x2i, int x3i, fp_t lambda, fp_t ** profile_derivatives){

  fp_t * derivative = new fp_t[7]-1;
  memset(derivative+1,0,7*sizeof(fp_t));

  for(int z=0;z<=Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t Aul = A[z][i][ii];
        int upper_map = rmap[z][i];
        fp_t profile =  current_profile[x1i][x2i][x3i][tmap[z][i][ii]];
          
        derivative[1] += profile_derivatives[1][tmap[z][i][ii]] * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi;
        derivative[2] += profile_derivatives[2][tmap[z][i][ii]] * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi;
        derivative[3] += profile_derivatives[3][tmap[z][i][ii]] * pop[x1i][x2i][x3i].n[z][i] * Aul * h * c / lam / 4.0 / pi;
  }       

  return derivative;
}


// And then even more. These are really important and maybe we should think about putting all this in a separate file:

fp_t * atom::dR_i_cont_dqk(int z, int i, int x1i, int x2i, int x3i, fp_t Temp, fp_t Ne, fp_t I, fp_t lambda){
  
  // Explicit derivative of bound-free rates from level i of ionization stage z with respect to atmospheric parameters
  // Bound-free radiative transitions do not have any aspect which explicitly depends on the temperature/pressure/velocity/magneticfield

  fp_t * derivative = new fp_t [7] - 1;
  memset(derivative+1,0,7*sizeof(fp_t));
  return derivative;
}

fp_t * atom::dR_cont_i_dqk(int z, int i, int x1i, int x2i, int x3i, fp_t Temp, fp_t Ne, fp_t I, fp_t lambda){

  // // Explicit derivative of free-bound rates from level i of ionization stage z with respect to atmospheric parameters
  // Here there are some 

  fp_t * derivative = new fp_t [7] - 1;
  memset(derivative+1,0,7*sizeof(fp_t));

  fp_t sigma = (bf[z][i]) ? bf[z][i]->U(lambda * (1.0)) : 0.0;

  fp_t constant_part = (I+ 2.0 * h * c*c / pow(lambda, 5.0)) * sigma * lambda / h / c * fp_t(g[z][i]) / fp_t(g[z+1][0]) * saha_const;
  
  fp_t delta_E = ip[z] - ee[z][i];

  derivative[1] = exp((delta_E - h*c/lambda)/k/Temp) * pow(Temp,-1.5) * parent_atm->get_ne_lte_derivative(1,x1i,x2i,x3i);
  derivative[1] += exp((delta_E - h*c/lambda)/k/Temp) * (-1.5) * pow(Temp,-2.5) * Ne;
  derivative[1] += exp((delta_E - h*c/lambda)/k/Temp) * (-delta_E + h*c/lambda)/k/Temp/Temp * pow(Temp,-1.5) * Ne;
  derivative[1] *= constant_part;

  derivative[2] = constant_part * exp((delta_E - h*c/lambda)/k/Temp) * pow(Temp,-1.5) * parent_atm->get_ne_lte_derivative(2,x1i,x2i,x3i);
  return derivative;
}

fp_t atom::dR_i_cont_dI(int z, int i, fp_t Temp, fp_t Ne, fp_t lambda){
  
  fp_t sigma = (bf[z][i]) ? bf[z][i]->U(lambda * (1.0)) : 0.0;
  return sigma * lambda / h /c;  
}
fp_t atom::dR_cont_i_dI(int z, int i, fp_t Temp, fp_t Ne, fp_t lambda){
  
  fp_t sigma = (bf[z][i]) ? bf[z][i]->U(lambda * (1.0)) : 0.0;
  fp_t delta_E = ip[z] - ee[z][i];
  fp_t derivative = sigma * lambda / h / c * exp((delta_E - h*c/lambda)/k/Temp) * fp_t(g[z][i]) / fp_t (g[z+1][0]) * pow(Temp,-1.5) * Ne * saha_const;
  return derivative;

}
