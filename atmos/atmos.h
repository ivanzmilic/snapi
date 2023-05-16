#ifndef __ATMOS_H__  // __ATMOS_H__
#define __ATMOS_H__

class model;

#include <string.h>
#include "types.h"
#include "io.h"
#include "flags.h"
#include "obs.h"

#include "../grid/grid.h"

class atmosphere;

#include "atmol/atmol.h"

#include "acfg.h"
#include "../jsub/cfg.h"


#define ATMOS_FLAG_DEF       0x01UL
#define ATMOS_FLAG_INIT      0x02UL

#define ATMOS_FLAG_MASK      (ATMOS_FLAG_DEF|ATMOS_FLAG_INIT)

struct typemap{
  const char *type_name;
  uint08_t type_code;
  uint08_t equals(const char *type_in){ return (strcmp(type_in,type_name))?0:type_code; };
};

#define ATMOS_TYPE_SPINOR    1
#define ATMOS_TYPE_SPINOR3D  4
#define ATMOS_TYPE_MHD       2
#define ATMOS_TYPE_MURAM     3

#define ATMOS_GEOM_PP      1
#define ATMOS_GEOM_2D      2
#define ATMOS_GEOM_3D      3

#define ATMOS_RTS_QSC      1
#define ATMOS_RTS_BEZ      2

#define FILETYPES {{"SPINOR",ATMOS_TYPE_SPINOR},{"SPINOR3D",ATMOS_TYPE_SPINOR3D},{"MHD",ATMOS_TYPE_MHD},{"MuRAM",ATMOS_TYPE_MURAM},{0,0}}


class atmosphere:public grid{
protected:
  class reg flags;
  int08_t rtstype; // RT solver type
//
//  class param *parm; // parameterization: external?
//
  class atmol **atml;
  uint16_t natm;
//
  char *id,*fname;
  uint08_t ftype,gtype;
//
  fp_t ***T;              // temperature
  fp_t ***rho;            // mass density
  fp_t ***Nt,***Ne;       // electron pressure
  fp_t **** Ne_lte_der; // response of electron density in lte to all the perturbations
  fp_t ***Bx,***By,***Bz; // magnetic field
  fp_t ***Vx,***Vy,***Vz; // velocity field
  fp_t ***Vt;             // turbulent velocity
  fp_t ***tau_referent; // referent optical depth scale, from the top of the atmosphere down. 
  
  fp_t ***op_referent; // referent opacity, usually @ 500 nm, but can be any sort of mean opacity, this is by definition in LTE
  fp_t *****op_referent_derivative; // derivative (strictry local), of the referent opacity. It is strictly local because it is in LTE

  fp_t boundary_condition_for_rt; // I am not sure that this is the best solution for an atmosphere, but this i easiest to work with at the moment.
  int tau_grid; // This is a 'switch' which determines whether we use tau grid as a coordinate grid. If we do, it is a grid in continuum optical
                 // depth at lambda=500. If we do not, we use h which is given in the input atmosphere. 
  fp_t * rt_grid; // This is grid we use to perform radiative transfer solution. It will either be geometrical grid (x3) or it will be optical depth grid
                  // (tau_referent)

  // Opacity fudge values:
  fp_t * lambda_of;
  fp_t * value_of;
  int N_of;
  int conserve_charge; // whether we conserve charge or not

  // Atomic level populations to be stored and input / output when neeed, a la RH
  // Not sure if this is a good idea, but let's go.
  int use_atm_lvls; // Whether we are using atomic levels as a part of the atmospheric model
  int n_lvls; // How many atomic levels are we considering, default is 0;
  fp_t **** atm_lvl_pops; // Level populations 
//
//
  virtual int08_t resize(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
//
//
  fp_t* add(fp_t*,fp_t*,int32_t);
  fp_t** add(fp_t**,fp_t**,int32_t,int32_t,int32_t,int32_t);
  fp_t*** add(fp_t***,fp_t***,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  fp_t ****add(fp_t ****, fp_t ****, int32_t,int32_t,int32_t,int32_t,int32_t,int32_t, int32_t, int32_t);
  fp_t***** add(fp_t*****, fp_t*****, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
  fp_t****** add(fp_t******, fp_t******, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
  fp_t******* add(fp_t*******, fp_t*******, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);


// general utilities
  fp_t ***project(fp_t***,fp_t***,fp_t***,fp_t,fp_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  fp_t ****transform(fp_t***,fp_t***,fp_t***,fp_t,fp_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  void transform_responses(fp_t ****, fp_t, fp_t, int, int);
// atmos_rad.cc: radiative quantities
// frequency-by-frequency
  fp_t ***opacity(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ***emissivity(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);

  fp_t ***opacity_lte(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ***emissivity_lte(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);

  fp_t ***opacity_active_only(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ***emissivity_active_only(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);

  // Active opacity/emissivity. All wavelengths simultaneously.
  // Written taking into account that this is a method, so we don't have
  // to pass everything.
  
  fp_t opacity_continuum(fp_t, fp_t, fp_t, int, int, int);
  fp_t ** opacity_continuum_derivative(fp_t, fp_t, fp_t, int, int, int);

  
  fp_t ***thomson_sc(fp_t***,fp_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  fp_t ***thomson_em(fp_t***,fp_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
// And then their corresponding sister functions:
  fp_t ***** opacity_pert(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ***** emissivity_pert(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ***** opacity_pert_lte(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ***** emissivity_pert_lte(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ***** thomson_sc_pert(fp_t***,fp_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  fp_t ***** thomson_em_pert(fp_t***,fp_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);

  
  int op_em_pert_numerical_scalar(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t,fp_t *****, fp_t *****);
  
// Vector versions, frequency by frequency:
  fp_t *****opacity_vector(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ****emissivity_vector(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ******* opacity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);
  fp_t ****** emissivity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***,fp_t****,fp_t,fp_t,fp_t);

// ======================================================================================================================
  // Special versions of functions which compute opacity and emissivity for full wavelength grid in one go.
  // We still account for the atmosphere as if it were 3D, although we might not need it like that.

  // Calculates opacity and emissivity all wavelengths at once, op and em together:
  int op_em_scalar_active_only(fp_t ***, fp_t ****, fp_t, fp_t, fp_t*,int32_t,fp_t ****,fp_t ****);
  int op_em_vector(fp_t ***, fp_t ****, fp_t,fp_t,fp_t*,int,fp_t******,fp_t*****); // With less arguments as most are already contained in the atmosphere
  int op_em_vector_pert(fp_t ***, fp_t ****, fp_t,fp_t,fp_t*,int,fp_t********,fp_t*******){return 0;}; // Basically the same as the one above, except taking different arguments
  int op_em_vector_plus_pert(fp_t ***, fp_t ****, fp_t,fp_t,fp_t*,int,fp_t******,fp_t*****,fp_t********,fp_t*******);

// =======================================================================================================================

// Referent optical depth scale related functions:
  virtual int compute_op_referent();
  virtual int compute_tau_referent();
  virtual int compute_op_referent_derivative(){return 0;};
  virtual void delete_op_referent_derivative();
  virtual void normalize_to_referent_opacity(fp_t ***, fp_t ***);
  virtual void normalize_to_referent_opacity(fp_t *****, fp_t ****); // Overloaded vector vesion
  virtual void normalize_to_referent_opacity(fp_t *****, fp_t ****, fp_t *******, fp_t ******); // Overloaded for perturbations of vector quantities
  virtual void de_normalize(fp_t ***, fp_t ***); // Basically undo what the function above does
  virtual void de_normalize(fp_t *****, fp_t ****); // Overloaded vector version

// point-by-point
  //fp_t *opacity(fp_t*,int32_t,int32_t,int32_t,int32_t);
  //fp_t *emissivity(fp_t*,int32_t,int32_t,int32_t,int32_t);
  //fp_t *thomson_op(fp_t,fp_t*,int32_t);
  //fp_t *thomson_em(fp_t,fp_t,fp_t*,int32_t);
// atmos_chm.cc: chemistry
  uint08_t chemeq(class atmol**,int,fp_t,fp_t,fp_t&,int32_t,int32_t,int32_t);  
// atmos_pop.cc: populations
  
  fp_t newpops(fp_t***,fp_t***,fp_t***,fp_t*,int32_t);
  // ^
  //changed so it returns the greatest relative difference between new and old pops
  int nltepops(void);
  int nltepops_taugrid(void);
  void ltepops(void);
  void popsetup(void);
  void popclean(void);
  void respsetup(void);
  void respclean(void);
  fp_t ne_derivative(int, int, int);
  void ne_lte_derivatives();
  void clear_ne_lte_derivatives();
  int atm_pop_setup(void);
  int atm_pop_fill(void);
  int atm_pop_clean(void);


// atmos_rts.cc: radiative transfer solver[s]
  virtual fp_t *anglesetup(fp_t*&,fp_t*&,int&);
  virtual int formal(fp_t*, fp_t***,fp_t***,fp_t***,fp_t***,fp_t,fp_t, fp_t);
  virtual int formal_with_lambda_operator(fp_t*, fp_t***,fp_t**,fp_t***,fp_t***,fp_t,fp_t, fp_t);
  virtual int formal_with_responses(fp_t*, fp_t ***, fp_t ***, fp_t **, fp_t **, fp_t ***, fp_t ***, fp_t,fp_t,fp_t);
  // polarized:
  virtual int formal(fp_t *, fp_t ****,fp_t ***,fp_t *****,fp_t ****,fp_t ,fp_t , fp_t ); // Formal for polarized
  // radiative transfer solver for perturbations:
  virtual int formal_pert_numerical(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  virtual int formal_pert_analytical(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  virtual int formal_pert_jcdti(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  // polarized radiative transfer solver for perturbations
  virtual int formal_pert_numerical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary)
  {return 0;};
  virtual int formal_pert_numerical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary, int N_parameters)
  {return 0;};
   virtual int formal_pert_analytical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary)
  {return 0;};
  virtual int formal_with_responses_jcdti(fp_t *, fp_t ****,fp_t ***,fp_t *****,fp_t ****, fp_t *****, fp_t ****, fp_t ****, fp_t ,fp_t , fp_t )
  {return 0;};
  virtual int formal_with_responses_full(fp_t *, fp_t ****,fp_t ***,fp_t *****,fp_t ****, fp_t *****, fp_t ****, fp_t ****, fp_t ,fp_t , fp_t )
  {return 0;};

  virtual int formal_pert_analytical_taugrid(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  virtual int formal_pert_numerical_taugrid(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****,fp_t ***, fp_t ****, fp_t, fp_t, fp_t){return 0;};
// atmos_fio.cc: file I/O
  int08_t read_atmos(const char*,const char*,uint08_t,io_class*);
  int08_t read_spinor(const char*,const char*,io_class*);
  int08_t read_spinor3d(const char*,const char*,io_class*);
  int08_t read_mhd(const char*,const char*,io_class*);
  int08_t read_muram(const char*,const char*,io_class*);
public:
  atmosphere(acfg*,io_class&);
  atmosphere(uint08_t*,int32_t&,uint08_t,io_class&);

  virtual atmosphere* extract(int i, int j, io_class&);

  virtual int build_from_nodes(model *);
  virtual int interpolate_from_nodes(model *); // the same as above except it does not re-evaluate HE
  virtual int enforce_hequilibrium(){return 0;}; // enforces hydrostatic equilibrium
  virtual fp_t ** calculate_dN_dT(){return 0;}; 

  virtual ~atmosphere(void);
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
//
  virtual int08_t init(const char*,io_class*);
//  obs observable(lamda,prm);
//  void restruct(obs *difference,fp_t *meritfunc());
  virtual void set_grid(int);
  virtual int get_grid(){
    return tau_grid; 
  };
//
// atmos_obs.cc: generate observable
  virtual observable *obs_scalar(fp_t,fp_t,fp_t*,int32_t);
  virtual observable *obs_scalar_responses(fp_t,fp_t,fp_t*,int32_t, fp_t ***);
  virtual observable *obs_scalar_num_responses(fp_t,fp_t,fp_t*,int32_t, fp_t ***);
// These are now responses to nodes.
  virtual observable *obs_scalar_num_responses_to_nodes(model *, fp_t, fp_t, fp_t *, int32_t, fp_t **); // The same as before,except it will also return the derivatives to nodes (written in the last argument)
  virtual observable *obs_scalar_responses_to_nodes(model *, fp_t, fp_t, fp_t *, int32_t, fp_t **);

// Duplicated, where tau is an independent variable:
  virtual observable *obs_scalar_tau(fp_t,fp_t,fp_t*,int32_t);
  virtual observable *obs_scalar_num_responses_tau(fp_t,fp_t,fp_t*,int32_t, fp_t ***);
  virtual observable *obs_scalar_responses_tau(fp_t,fp_t,fp_t*,int32_t, fp_t ***);
// Responses to nodes:
  virtual observable *obs_scalar_num_responses_to_nodes_tau(model *, fp_t, fp_t, fp_t *, int32_t, fp_t **); // The same as before,except it will also return the derivatives to nodes (written in the last argument)
  virtual observable *obs_scalar_responses_to_nodes_tau(model *, fp_t, fp_t, fp_t *, int32_t, fp_t **);

// Vector case. Here we have already generalized. So no need to split between "tau" and "geometrical" functions. This should be cleaned up in 
// the 'final' version of the code
  virtual observable *obs_stokes_responses_to_nodes(model *, fp_t, fp_t, fp_t *, int32_t, fp_t ***, fp_t filter_width); // 'New' version
                                                                                                         // hence the 3D vector in last argument
  virtual observable *obs_stokes_num_responses_to_nodes(model *, fp_t, fp_t, fp_t *, int32_t, fp_t ***, fp_t filter_width); // The same as above, except it is for vector case.
  virtual int scale_rf(fp_t ***, model*, int, int, fp_t *, fp_t *);
  virtual int scale_corrections(fp_t *, model*, int);                                                                                                      // hence the 3D vector in last argument
//
  virtual fp_t *test_stokes(fp_t,fp_t,fp_t*,int32_t); // Keep this for debugging purposes. In the end you can delete it.
  virtual observable *obs_stokes(fp_t,fp_t,fp_t*,int32_t); // Same as the obs_scalar
  virtual observable *obs_stokes_responses(fp_t,fp_t,fp_t*,int32_t, fp_t ****); // Same as obs_scalar_responses, except it works for full Stokes vector
  virtual observable *forward_evaluate(fp_t theta, fp_t phi, fp_t * lambda, int nlambda,fp_t scattered_light, fp_t qs, fp_t spectral_broadening);
  virtual fp_t * calc_residual(fp_t **, fp_t **, int, int, int*, fp_t *);
  virtual fp_t calc_chisq(int, int, int*, fp_t *, fp_t *, fp_t *);
  virtual int look_for_best_lambda(fp_t &lm_parameter, fp_t ** JTJ, int N_parameters,
  fp_t * rhs, model * model_current, fp_t theta, fp_t phi, fp_t * lambda, int nlambda, fp_t scattered_light,
  fp_t qs_level, fp_t spectral_broadening, fp_t ** S_to_fit, int n_stokes_to_fit, int * stokes_to_fit,
  fp_t * ws, fp_t * noise, fp_t metric_old);
  
  virtual observable *obs_stokes_responses(fp_t,fp_t,fp_t*,int32_t, fp_t ***, model*, fp_t ****); // Same as obs_scalar_responses, except it works for full Stokes vector
  virtual observable *obs_stokes_num_responses(fp_t,fp_t,fp_t*,int32_t, fp_t ****); // 
  
  // Debug:
  virtual int optical_depth_scale(fp_t ***, fp_t ***, fp_t, fp_t);
  virtual fp_t * compute_full_lambda_operator(fp_t ***, int, int){
    return 0;
  };

// atmos_fit.cc various fitting examples, routines and testing:
  //virtual observable *scalar_lm_fit(observable *, fp_t, fp_t, fp_t *, int); // Function which performs a levenberg-marquard fit
  virtual observable *stokes_lm_fit(observable *, fp_t, fp_t, model *); // Function which performs a levenberg-marquard fit
  virtual observable *stokes_lm_nodeless_fit(observable *, fp_t, fp_t); // LM fits trying out nodeless inversion
  virtual fp_t * calculate_svd_corrections(fp_t ****,fp_t *, fp_t, int, int);
  virtual fp_t ** calculate_svd_corrections(fp_t ****, fp_t *, fp_t, int);
  virtual fp_t ** calculate_legendre_corrections(fp_t ****, fp_t *, fp_t, int, int, int);
  virtual void regularize_hessian(fp_t **, fp_t *, model*);
  virtual void regularize_parameter(fp_t **, fp_t *, model *, int);

  virtual int polish_extreme_values();
  
  fp_t get_pop(int, int, int, int, int, int);
  fp_t get_pop(int, int, int, int, int);
  fp_t get_T(int, int, int);
  fp_t get_Ne(int, int, int);
  fp_t set_Ne(int, int, int,fp_t);
  fp_t get_Nt(int, int, int);
  fp_t get_ne_lte_derivative(int, int, int, int);
  fp_t get_neutral_H_derivative_lte(int,int, int, int);
  fp_t get_partf(int who, int z, fp_t T_in, fp_t Ne_in);
  fp_t get_vt(int, int, int);
  fp_t * get_magnetic_field(int, int, int);
  int execute_chemeq_for_point(int, int, int); // Just to execture chemeq in point x1i, x2i, x3i, populations will be automatically updated
  void set_Temp(int, int, int, fp_t); // We might have merged this with the previous one, but it might be handy to "remotely" change the atmospheric parameters
  void set_Nt(int, int, int, fp_t);
  fp_t get_x1(int index){
    return x1[index];
  }
  fp_t get_x2(int index){
    return x2[index];
  }
  fp_t get_x3(int index){
    return x3[index];
  }
  fp_t get_op_referent(int x1i, int x2i, int x3i){
    return op_referent[x1i][x2i][x3i];
  };
  fp_t get_op_referent_der(int p, int x3k, int x1i, int x2i, int x3i){
    return op_referent_derivative[p][x3k][x1i][x2i][x3i];
  };
  int get_N_depths();
  fp_t ** return_as_array();
  int copy_from_array(fp_t **);
  virtual fp_t get_opacity_fudge(fp_t lambda);
  void set_conserve_charge(int);
  int get_conserve_charge();
//
// ----------------------------------------------------------------------------------------------------------------
// -------------- RESPONSE FUNCTIONS RELATED FUNCTIONS --------------------------------------------------------------

void compute_responses(){

} // Compute complete set of nlte responses
void compute_lte_population_responses(){

}
int compute_nlte_population_responses(int lvl_of_approximation);
void compute_nlte_population_responses_numerical(int from, int to);
void compute_anisotropy_responses();

void compute_nlte_population_responses_taugrid(int lvl_of_approximation);
void compute_nlte_population_responses_numerical_taugrid(int from, int to);


fp_t **** compute_intensity_response_numerical(int from, int to, fp_t theta, fp_t phi, fp_t lambda);
};

atmosphere *atmos_new(acfg*,io_class&);
atmosphere *atmos_new(uint08_t*,int32_t&,uint08_t,io_class&);

// ---------------------------------------------------------------------------------------------------------------

class model{

private:
  int N_nodes_temp;
  int N_nodes_vt;
  int N_nodes_vs;
  int N_nodes_B; 
  int N_nodes_theta;
  int N_nodes_phi;

  int N_depths; // This depends on the specific atmosphere we are using:
  fp_t tau_min;
  fp_t tau_max;
  
  // Does it make sense to have separate number of nodes for different components of velocity/magnetic field? 
  // To me it does not. But maybe we want to allow for possibillity for that in the code?
  // Let's leave this as a starting point. If we get things working we will code new one.

  // Lets define nodes. For each "couple" of arrays, first ones are x coordinates of nodes, second ones are values of nodes.

  // Temperature nodes
  fp_t * temp_nodes_tau;
  fp_t * temp_nodes_temp;
  int temp_reg_type;
  fp_t temp_reg_alpha;

  // Microturbulent velocity nodes
  fp_t * vt_nodes_tau;
  fp_t * vt_nodes_vt;
  int vt_reg_type;
  fp_t vt_reg_alpha;

  // Systematic velocity nodes

  fp_t * vs_nodes_tau;
  fp_t * vs_nodes_vs;
  int vs_reg_type;
  fp_t vs_reg_alpha;

  // Magnetic field nodes

  fp_t * B_nodes_tau;
  fp_t * B_nodes_B;
  int B_reg_type;
  fp_t B_reg_alpha;

  // Theta nodes 

  fp_t * theta_nodes_tau;
  fp_t * theta_nodes_theta;
  int theta_reg_type;
  fp_t theta_reg_alpha;

  // Phi nodes: 

  fp_t * phi_nodes_tau;
  fp_t * phi_nodes_phi;
  int phi_reg_type;
  fp_t phi_reg_alpha;

  // Mapping matrix which describes how is every parameter mapped onto perturbations of other ones:
  fp_t *** response_to_parameters;

public:

  // We will keep these virtual so we can make some derived class from this

  model(); // Creates default values with zero free parameters
  model(int, int, int, int); // Creates model with values of nodes 
  model(int, int, int, int, int, int); // Creates model with nodes for all 6 parameters
  ~model(void);
  model(mcfg*,io_class&); // Create from config file
  model(uint08_t*,int32_t&,uint08_t,io_class&); // Unpack from buffer

  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);

  // Now methods which set nodes
  int set_temp_nodes(fp_t *, fp_t *);
  int set_vt_nodes(fp_t *, fp_t *);
  int set_vs_nodes(fp_t *, fp_t *);
  int set_B_nodes(fp_t *, fp_t *);
  int set_theta_nodes(fp_t *, fp_t *);
  int set_phi_nodes(fp_t *, fp_t *);

  int get_N_nodes_temp();
  fp_t * get_temp_nodes_tau();
  fp_t * get_temp_nodes_temp();
  int get_temp_reg_type(){return temp_reg_type;};
  fp_t get_temp_reg_alpha(){return temp_reg_alpha;};
  void set_temp_reg(int,fp_t);

  int get_N_nodes_vt();
  fp_t * get_vt_nodes_tau();
  fp_t * get_vt_nodes_vt();
  int get_vt_reg_type(){return vt_reg_type;};
  fp_t get_vt_reg_alpha(){return vt_reg_alpha;};
  void set_vt_reg(int,fp_t);

  int get_N_nodes_vs();
  fp_t * get_vs_nodes_tau();
  fp_t * get_vs_nodes_vs();
  int get_vs_reg_type(){return vs_reg_type;};
  fp_t get_vs_reg_alpha(){return vs_reg_alpha;};
  void set_vs_reg(int,fp_t);

  int get_N_nodes_B();
  fp_t * get_B_nodes_tau();
  fp_t * get_B_nodes_B();
  int get_B_reg_type(){return B_reg_type;};
  fp_t get_B_reg_alpha(){return B_reg_alpha;};
  void set_B_reg(int,fp_t);

  int get_N_nodes_theta();
  fp_t * get_theta_nodes_tau();
  fp_t * get_theta_nodes_theta();
  int get_theta_reg_type(){return theta_reg_type;};
  fp_t get_theta_reg_alpha(){return theta_reg_alpha;};
  void set_theta_reg(int,fp_t);

  int get_N_nodes_phi();
  fp_t * get_phi_nodes_tau();
  fp_t * get_phi_nodes_phi();
  int get_phi_reg_type(){return phi_reg_type;};
  fp_t get_phi_reg_alpha(){return phi_reg_alpha;};
  void set_phi_reg(int,fp_t);

  int get_N_nodes_total();
  int perturb_node_value(int, fp_t);
  fp_t get_parameter(int);
  int which_parameter(int);

  int print();
  int print_to_file(FILE *);

  int set_response_to_parameters(fp_t ***, int);
  fp_t *** get_response_to_parameters();

  int correct(fp_t * correction);
  int set_parameters(fp_t * par_input);
  int bracket_parameter_values();
  int polish_angles();
  
  int cpy_values_from(model *);

  fp_t get_tau_min(){
    return tau_min;
  }
  fp_t get_tau_max(){
    return tau_max;
  }
  void set_tau_min(fp_t input){tau_min=input;};
  void set_tau_max(fp_t input){tau_max=input;};
};

model * clone(model *);
model * model_new(int, int, int, int);
model * model_new(int, int, int, int, int, int);
model * model_new(mcfg*,io_class&);
model * model_new(uint08_t*,int32_t&,uint08_t,io_class&);

class modelcube{
private:
  int nx,ny;
  int N_parameters;
  int N_nodes_temp, N_nodes_vt, N_nodes_vs, N_nodes_B, N_nodes_theta, N_nodes_phi;
  fp_t * temp_nodes_tau;
  fp_t * vt_nodes_tau;
  fp_t * vs_nodes_tau;
  fp_t * B_nodes_tau;
  fp_t * theta_nodes_tau;
  fp_t * phi_nodes_tau;

  fp_t *** data;
public:
  modelcube();
  modelcube(model*, int nx, int ny);
  ~modelcube();
  void add_model(model*,int i, int j);
  fp_t *** get_data(int &nx_in, int &ny_in, int &np_in);
  void set_data(fp_t ***);
  void simple_print(const char*);
};

#endif              // __ATMOS_H__
