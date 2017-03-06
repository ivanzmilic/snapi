#ifndef  __ATOM_H__ // __ATOM_H__
#define  __ATOM_H__

#include "types.h"
#include "io.h"
#include "atmol.h"
#include "partf.h"
#include "bfcs.h"
#include "colr.h"

struct pps{
  fp_t **n; // level populations [zl..zh][0..nl]: sum(nl[z])=N[z]?
  fp_t *N;  // total number density for each ionization stage: sum(N)=Nt; // But then we have conflict with Na, what is the issue here
  fp_t Nt;  // total number density for this element
  fp_t Na;  // total number density of active levels 
//
  pps(void){
    n=0;
    N=0;
  };
  ~pps(void){
    if(n) delete[] n[0];
    delete[] n;
    delete[] N;
  };
};

// level flags

#define FL_NONE 0x01
#define FL_LTE  0x02

class atom:public atmol{
protected:
  uint08_t Z;          // 0..Z
  fp_t abund;          // abundance
  uint08_t NLTE;            // is the atom in NLTE (if not, levels are assumed to be in LTE but lines are included)
//
  fp_t **ee;           // level energy
  fp_t ***A;      // spontaneous transition coefficient
  fp_t ***B;      // stimulated transition coefficient
 
  fp_t ***osc_str; // oscillator strength, easier to just pre-compute and then use
  fp_t *** col_dam_cross_section; // collisional damping cross-section, also easier to pre-compute
  fp_t *** alpha_col_dam; // collisional damping coefficient alpha


  uint16_t *nl;        // number of levels/ion
  uint16_t nmap,*lmap; // mapping of level index -> z,l
  uint16_t **rmap;     // reverse mapping (z,l -> level index)
  uint08_t *zmap;      // mapping of level index -> z,l
  uint08_t **j,**g;
  uint08_t **flags;

  // Additional quantities which characterize a level
  fp_t ** g_lande; // lande factor of the level in question
  int ** nm; // Number of 'magnetic' sub transitions for each transition
  fp_t ** S_p; // Relative strengths of p transitions
  fp_t ** S_b; // Relative strengths of b transitions
  fp_t ** S_r; // Relative strengths of r transitions
  fp_t ** delta_lambda_p; // Zeeman shifts of p transitions
  fp_t ** delta_lambda_b; // Zeeman shifts of b transitions
  fp_t ** delta_lambda_r; // Zeeman shifts of r transitions

  uint08_t ** j_qn; // Quantum number j, total angular momentum
  uint08_t ** l_qn; // Quantum number l, orbital quantum number

  fp_t *ip;            // ionization potentials
//
  class pf **partf;
  class bfcs ***bf;    // bound-free cross-section
  class colr ****cr;   // collisional rates
//
// space variant variables
// The rate matrix scales as n^2, for 50 levels at 288*288*200 points this is approximately 1GB
// This is a major RAM consumer... 
// It gets worse: if the density matrix formalism is used, the populations become coupled (quantum coherencies). 
// This leads to off-diagonal elements in the density matrix (uses lots of memory)
// We want to describe this somehow in a blocked form, to avoid a big matrix with many zero elements
//
  int32_t x1l,x1h,x2l,x2h,x3l,x3h;
  
// Populations:
  struct pps ***pop;

// List of transitions is needed? It is needed, in a way.
// This is transition related stuff:
  uint32_t ***tmap,ntr;
  // This is inverse mapping:
  uint32_t ** inverse_tmap;


  fp_t ****Jb; // angular and frequency integrated intensity up->down
  fp_t ****Ju; // angular and frequency integrated intensity down-up
  // Random thoughts in the morning (02/06/2015): These are not "frequency and angle integrated intensities." It is better to think of them as terms which contain 
  // rates, or simple as radiative rates themselves. It turns out that in the expressions this really has the dimension of the intensity, but it is already 
  // annoying me with number of missunderstandings. The other think is that the rates themselves 

  // Oom estimates for scattering polarization.

  fp_t **** J_02; // J_02 component of radiation tensor, which most directly influences emergent scattering polarization.
  fp_t **** J_02_responses; // Responses, just 1D in space, 6D array is too much
  fp_t **** J_00_responses; // Responses for J
  
  fp_t ****Ls; // approximate lambda operator
  fp_t ****current_profile; // This saves the value of the current profile function for each transition @ each point
  fp_t ****norm; // This is intended to normalize the wavelength integration
  fp_t *****norm_derivative; // This is the derivative of the 

  fp_t *****Broyden; // Inverse broyden matrix <--- At the moment un-necesary.

  fp_t ** beta_Temp; // The 'local' or 'pseudo local' part of the response to the temperature. 
  fp_t ** beta_density;
  fp_t ** beta_v_micro; // The same for the microturbulent velocity. 
  fp_t ** beta_v_r;
  // In principle we might add more, but these are the most important at the moment...

  // This is the big response matrix
  fp_t ** response_matrix;

  // These are the reponses themselves:
  fp_t *** level_responses;
  fp_t *** ionization_stage_responses; 

  fp_t* add(fp_t*,fp_t*,int32_t);
  fp_t*** add(fp_t***,fp_t***,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  fp_t ****add(fp_t ****, fp_t ****, int32_t,int32_t,int32_t,int32_t,int32_t,int32_t, int32_t, int32_t);
  fp_t***** add(fp_t*****, fp_t*****, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
  fp_t****** add(fp_t******, fp_t******, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
  fp_t******* add(fp_t*******, fp_t*******, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
//
  fp_t *rayleigh_em(fp_t *lambda,int32_t nlambda);
  fp_t *rayleigh_op(fp_t *lambda,int32_t nlambda);
  virtual fp_t *freefree_op(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
  virtual fp_t *freefree_em(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
  virtual fp_t *boundfree_op(fp_t*,int32_t,int32_t,int32_t,int32_t);
  virtual fp_t *boundfree_em(fp_t*,int32_t,int32_t,int32_t,int32_t);
  fp_t *boundbound_em(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
  fp_t *boundbound_op(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
//
  virtual fp_t ***freefree_op(fp_t***,fp_t***,fp_t***,fp_t);
  virtual fp_t ***freefree_em(fp_t***,fp_t***,fp_t***,fp_t);
  virtual fp_t ***boundfree_op(fp_t***,fp_t);
  virtual fp_t ***boundfree_em(fp_t***,fp_t);
  virtual fp_t ***rayleigh_op(fp_t lambda);
  virtual fp_t ***rayleigh_em(fp_t lambda);

  virtual fp_t ***** freefree_op_pert(fp_t***,fp_t***,fp_t***,fp_t);
  virtual fp_t ***** freefree_em_pert(fp_t***,fp_t***,fp_t***,fp_t);
  virtual fp_t ***** boundfree_op_pert(fp_t***,fp_t);
  virtual fp_t ***** boundfree_em_pert(fp_t***,fp_t);
  
  fp_t ***boundbound_em(fp_t***,fp_t***,fp_t***,fp_t***,fp_t);
  fp_t ***boundbound_op(fp_t***,fp_t***,fp_t***,fp_t***,fp_t);

//
  // Overloaded versions of b-b functions, which also take magnetic field. 
  // Oringally they took the concentration of the collisional partner but we have dumped them. 
  fp_t ***boundbound_em(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t);                       // |
  fp_t ***boundbound_op(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****, fp_t);
  fp_t ***** boundbound_em_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t);                       // |
  fp_t ***** boundbound_op_pert(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****, fp_t);

  // Then there are functions which give us the partial response of the opacity and emissivity to specific level
  // Level is in this case enumerated with 0.... nmap-1, these assume 1D! 

  fp_t op_derivative_to_level(int, int, fp_t lambda);
  fp_t em_derivative_to_level(int, int, fp_t lambda);

  fp_t * op_derivative_explicit(int, int, int, fp_t lambda); // These return arrays, because we have derivatives w.r.t different atmospheric parameters
  fp_t * em_derivative_explicit(int, int, int, fp_t lambda);

  // Do they really have to assume 1D? Not really

  fp_t bb_op_derivative_to_level(int x1i, int x2i, int x3i, int i, fp_t lambda);
  fp_t bb_em_derivative_to_level(int x1i, int x2i, int x3i, int i, fp_t lambda);
  fp_t * bb_op_derivative_explicit(int x1i, int x2i, int x3i, fp_t lambda, fp_t ** profile_derivatives);
  fp_t * bb_em_derivative_explicit(int x1i, int x2i, int x3i, fp_t lambda, fp_t ** profile_derivatives);
                          // |
// ---------------------------------------------------------------------------------------------- |
//
// -----------------------------------------------------------------------------------------------|
// Some functions for the computation of the weights. First of all, the function which returns the|
// profile for given transition and we use it to compute the weights if not more                  |
//
  int store_profile (int x1i, int x2i, int x3i, int transition_number, fp_t lambda);
  fp_t get_profile (int x1i, int x2i, int x3i, int transition_number, fp_t lambda);
  fp_t compute_profile(int x1i, int x2i, int x3i, int transition_number, fp_t lambda);                                                                                                
  int store_current_profile(int x1i, int x2i, int x3i, int transition_number, fp_t current_profile_in);
  fp_t get_current_profile(int x1i, int x2i, int x3i, int transition_number);

  // Profile derivatives with respect to various quantities:
  fp_t voigt_der_T(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos);
  fp_t voigt_der_density(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos);
  fp_t voigt_der_vt(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos);
  fp_t voigt_der_v(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos){
    return 0;
  }
  // Alternative approach to computing profile responses. Computing derivatives of x and a to atmospheric parameters in one
  // procedure. And then subsituting back in the derivative. Might end up actually faster in the end. And is also somewhat more compact. 
  fp_t compute_x_scalar(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos);
  fp_t compute_a_scalar(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos);
  int compute_xa_der_scalar(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos, fp_t&, fp_t&, fp_t*, fp_t*, fp_t*);

  fp_t * x_derivative(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos, fp_t **** B_vec, int trans_type, int m);
  fp_t * a_derivative(int x1i, int x2i, int x3i, int z, int i, int ii, fp_t lambda, fp_t *** vlos, fp_t **** B_vec, int trans_type, int m);
 // -----------------------------------------------------------------------------------------------|
  virtual fp_t damp_rad(fp_t **,uint16_t,uint16_t);
  virtual fp_t damp_col(fp_t,fp_t,fp_t,fp_t,fp_t,int08_t);
  virtual fp_t damp_col(fp_t,fp_t,fp_t,fp_t,fp_t, fp_t, int08_t); // Overloaded version, last fp_t number is the density of the collisional partner.
  
  // New ones, ones that we are actually using
  virtual fp_t damp_col(int ix1, int ix2, int ix3, int z, int i_from, int i_to, fp_t Temp, fp_t Ne, fp_t lambda_0); // This version is the one we are going to use from now on.
  virtual fp_t damp_col_der_T(int ix1, int ix2, int ix3, int z, int i_from, int i_to, fp_t Temp, fp_t Ne, fp_t lambda_0); // This version is the one we are going to use from now on.
  

  virtual int compute_damp_col(int, int, int);
  virtual int interpolate_col_damp_sp(fp_t, fp_t, int, int, int);
  virtual int interpolate_col_damp_pd(fp_t, fp_t, int, int, int);
  virtual int interpolate_col_damp_df(fp_t, fp_t, int, int, int);

  virtual fp_t broad_dop(fp_t,fp_t,fp_t);
 // -----------------------------------------------------------------------------------------------

  // New transitions:

  // bound - bound`

  virtual fp_t R_ij(int, int, int, fp_t);
  virtual fp_t C_ij_dummy(int, int, int, fp_t); 
  virtual fp_t C_ij(int, int, int, fp_t, fp_t);
  virtual fp_t C_ij_H(int, int, int, fp_t, fp_t);
  

  virtual fp_t R_i_cont(int z, int i, fp_t JJ, fp_t T);
  virtual fp_t R_cont_i(int z, int i, fp_t JJ, fp_t T, fp_t n_e);
  virtual fp_t C_i_cont(int z, int i, fp_t T, fp_t n_e);
  virtual fp_t C_cont_i(int z, int i, fp_t T, fp_t n_e);

  // Then we need to modify radiative rates with alo:
  virtual fp_t R_ij_local_ALO(int, int, int, fp_t, fp_t, fp_t, fp_t); 

  // And it will turn out that is is good to have a methods for computing the derivatives of collisions with respect to the temperature.

  virtual fp_t derivative_collisions_Temp(int, int, int, int, int);
  virtual fp_t collisional_rates(int, int, int, int, int);
  virtual fp_t derivative_collisions_full_temp(int, int, int, int, int);
  virtual fp_t derivative_collisions_full_density(int, int, int, int, int);

  // More derivative related functions, but now related to collisional rates:

  virtual fp_t * dR_i_cont_dqk(int z, int i, int x1i, int x2i, int x3i, fp_t Temp, fp_t Ne, fp_t I, fp_t lambda);
  virtual fp_t * dR_cont_i_dqk(int z, int i, int x1i, int x2i, int x3i, fp_t Temp, fp_t Ne, fp_t I, fp_t lambda);

  virtual fp_t dR_i_cont_dI(int z, int i, fp_t Temp, fp_t Ne, fp_t lambda);
  virtual fp_t dR_cont_i_dI(int z, int i, fp_t Temp, fp_t Ne, fp_t lambda);

//
  fp_t boundbound_op(uint08_t,uint16_t,uint16_t,fp_t,fp_t,fp_t,fp_t,fp_t,struct pps&);
  fp_t boundfree_op(uint08_t,uint16_t,fp_t,fp_t,struct pps&); // absorption
  fp_t **weight(fp_t,fp_t,fp_t,fp_t,fp_t,struct pps&,fp_t,fp_t,fp_t);
public:
  atom(atmcfg*,io_class&);
  atom(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~atom(void);
// transport/sync
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
//
// radiative properties
//
  virtual fp_t *opacity(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
  virtual fp_t *emissivity(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
  virtual fp_t ***opacity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  
  virtual fp_t opacity_continuum(fp_t, fp_t, fp_t, int, int, int);
  virtual fp_t ** opacity_continuum_pert(fp_t, fp_t, fp_t, int, int, int);
  
  virtual fp_t ***emissivity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  
  virtual fp_t ***emissivity_polarized_dummy(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t, fp_t);
  virtual fp_t *****emissivity_polarized_perturbation_dummy(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t, fp_t);
  
  virtual fp_t ***** opacity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ***** emissivity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
//  virtual fp_t *****opacity(fp_t***,fp_t***,fp_t***,fp_t***,fp_t,fp_t,fp_t);
//  virtual fp_t ****emissivity(fp_t***,fp_t***,fp_t***,fp_t***,fp_t,fp_t,fp_t);
//
// vector quantities:
  virtual fp_t *****opacity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ****emissivity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  // Then the perturbations:
  virtual fp_t *******opacity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ******emissivity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  
  // The only relevant ones for the vector case are b-b opacity and emissivity:
  fp_t ****boundbound_em_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t ****, fp_t);
  fp_t *****boundbound_op_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t ****, fp_t);
  fp_t ******boundbound_em_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t ****, fp_t);
  fp_t *******boundbound_op_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t ****, fp_t);

// perturbations of vector quantities:
  //virtual fp_t ******opacity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t){return 0;};
  //virtual fp_t ****emissivity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t){return 0;};
// 
  virtual void popsetup(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  virtual void popclean(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
//
  virtual void lte(fp_t***,fp_t***);
  virtual void lte(fp_t, fp_t, int, int, int);
  virtual void compute_active_population(fp_t ***, fp_t ***);
  virtual fp_t derivative_active_population(int, int, int);
  virtual fp_t derivative_active_population_density(int, int, int);


  virtual fp_t *getlambda(fp_t*,int32_t&,fp_t,fp_t,fp_t);
//  virtual void pupsetup(fp_t*,int32_t,fp_t*,uint16_t,fp_t*,uint16_t);
 
  virtual uint08_t rtsetup(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  virtual uint08_t rtclean(int,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  virtual void rtinit(void);
  virtual void prof_init(void);

  virtual void prof_setup(void);
  virtual void prof_clear(void);

  virtual void zeeman_setup(void);
  virtual void zeeman_clear(void);

  virtual void radiation_moments_setup();
  virtual void radiation_moments_init();
  virtual void radiation_moments_clean();

  virtual void compute_profile_norm(fp_t, fp_t, fp_t *, fp_t *, fp_t ***, int);
 
  virtual fp_t pops(atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t);
  virtual fp_t pops_broyden(atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t);
  virtual int initialize_broyden_ALO(fp_t **, int, int, int, fp_t);
  virtual fp_t * compute_F_broyden(int, int, int, fp_t, fp_t *, fp_t *);

  virtual void add(fp_t***,fp_t***,fp_t***,fp_t***,fp_t***,fp_t***,fp_t,fp_t,fp_t,fp_t);
  // New one:
  virtual void add(fp_t***, fp_t***, fp_t ***, fp_t, fp_t, fp_t);
  virtual void add_to_radiation_tensor(fp_t***, fp_t***, fp_t ***, fp_t, fp_t, fp_t, fp_t);
  virtual void add_to_radiation_tensor_perturbation(fp_t***, fp_t**,fp_t **, fp_t ***, fp_t ***, fp_t, fp_t, fp_t, fp_t, fp_t, fp_t ***, fp_t *****, fp_t *****);
// *************************************************************************************
// * routines needed for the new combined chemical/ionization equilibrium calculations *
// *************************************************************************************
  virtual void ionfrc(fp_t,fp_t,fp_t*&,fp_t*&,int&);
  virtual fp_t abundance(void){ return abund; };
  virtual void pupdate(fp_t,fp_t*,int,int32_t,int32_t,int32_t);
  virtual void randomize_populations(fp_t);
  virtual void print_populations();
//
  virtual void info(void);

// New ones: Milic, February 2015
  virtual fp_t get_total_pop(int x1i, int x2i, int x3i){
    return pop[x1i][x2i][x3i].Nt;
  }

  virtual fp_t get_pop(int x1i, int x2i, int x3i, int index_ion);
  virtual fp_t get_pop(int x1i, int x2i, int x3i, int index_ion, int index_lvl);
  virtual fp_t get_active_pop(int x1i, int x2i, int x3i);

  virtual fp_t get_J(int x1i, int x2i, int x3i, int transition);
  virtual fp_t get_L(int x1i, int x2i, int x3i, int transition);
  virtual fp_t get_norm(int x1i, int x2i, int x3i, int transition);

  virtual void print_radiation_field_tensor();

  // -----------------------------------------------------------------------------------------------------------------------------------
  // Response fuction related stuff:
  virtual int responses_setup();
  virtual int responses_init();
  virtual int responses_clear();
  virtual int add_response_contributions(fp_t***, fp_t**,fp_t **, fp_t ***, fp_t ***, fp_t, fp_t, fp_t, fp_t, fp_t, fp_t ***, fp_t *****, fp_t *****);
  virtual int add_response_contributions_taugrid(fp_t***, fp_t ***, fp_t *****, fp_t**,fp_t **, fp_t ***, fp_t ***, fp_t, fp_t, fp_t, fp_t, fp_t, fp_t ***, fp_t *****, fp_t *****);
  virtual int add_pops_to_response(int , int);
  virtual int subtract_pops_from_response(int, int);
  virtual int divide_responses_by_step(int,int, fp_t);
  virtual int print_population_responses(const char *, int);
  virtual int print_population_responses(const char *, int, int);
  virtual void compute_nlte_population_responses();
  virtual void compute_lte_population_responses();
  virtual void add_to_ion_responses(int param,int x3i, fp_t sign);
  virtual void divide_ion_responses(int param,int x3i, fp_t step);
  virtual fp_t get_population_response(int parameter, int x3k, int x1i, int x2i, int x3i, int z);
  
  virtual int check_if_nlte();
  
};

atmol *atom_new(atmcfg*,io_class&);
atmol *atom_new(uint64_t,uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in);

#endif              // __ATOM_H__
