#ifndef  __ATMOL_H__ // __ATMOL_H__
#define  __ATMOL_H__

#include <string.h>

#include "types.h"
#include "io.h"

#include "../../grid/grid.h"
#include "../atmos.h"

#define EC_OK             0x00
#define EC_DIM_MISMATCH   0x01
#define EC_POP_RAD_DEFECT 0x02
#define EC_POP_ION_DEFECT 0x04

class atmol;

struct ion{
  atmol *atom; // atom
  int08_t z;  // ionization stage
};

class atmol{
protected:
  io_class &io;
  uint64_t numid;
  char *name,*id;
  fp_t mass;

  atmosphere * parent_atm; // Pointer to the parent atmosphere which is used to compute all the 

public:
  atmol(const char*,const char*,io_class&);
  atmol(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~atmol(void);
// transport
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
  const char *get_frm(void){ return id; };
  uint64_t get_id(void){ return numid; };
  int08_t has_id(uint64_t numid_in){ return numid==numid_in; };
  int08_t has_id(const char *id_in){ return !strcmp(id,id_in); };
// opacity sources
  virtual fp_t *opacity(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
  virtual fp_t *emissivity(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t);
// we are using these two:
  virtual fp_t ***opacity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t opacity_continuum(fp_t, fp_t, fp_t, int, int, int){
    return 0;
  };
  virtual fp_t ** opacity_continuum_pert(fp_t, fp_t, fp_t, int, int, int){
    return 0;
  };
  virtual fp_t ***emissivity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ***emissivity_polarized_dummy(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t, fp_t){
    return 0;
  };
  virtual fp_t *****emissivity_polarized_perturbation_dummy(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t, fp_t){
    return 0;
  };
// and then their perturbations
  virtual fp_t ***** opacity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ***** emissivity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
// vector quantities:
  virtual fp_t *****opacity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ****emissivity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  // perturbations of vector quantities:
  virtual fp_t ******* opacity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t){return 0;};
  virtual fp_t ******  emissivity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t){return 0;};
  
  // Merged version with all wavelength points at once:
  virtual int op_em_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t*,int,fp_t ******, fp_t *****){return 0;};
  virtual int op_em_vector_plus_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t*,int,
    fp_t ******, fp_t *****, fp_t ********, fp_t *******){return 0;};
 
  
// compute populations
  virtual void popsetup(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t){};
  virtual void popclean(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t){};
  virtual void lte(fp_t***,fp_t***){};
  virtual void lte(fp_t, fp_t, int, int, int){};
  virtual void compute_active_population(fp_t ***, fp_t ***){};
  virtual fp_t derivative_active_population(int, int, int){
    return 0;
  }

  virtual void compute_profile_norm(fp_t, fp_t, fp_t *, fp_t *, fp_t ***, int){};

  virtual fp_t *getlambda(fp_t *lambda,int32_t&,fp_t,fp_t,fp_t){ return lambda; };
  virtual uint08_t rtsetup(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  virtual uint08_t rtclean(int,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  virtual void rtinit(void);
  virtual void prof_init(void);
  virtual void zeeman_setup(void){};
  virtual void zeeman_clear(void){};
  virtual void radiation_moments_setup(){};
  virtual void radiation_moments_init(){};
  virtual void radiation_moments_clean(){};

  virtual int responses_setup(){
    return 0;
  }
   virtual int responses_init(){
    return 0;
  }
  virtual int responses_clear(){
    return 0;
  }
  virtual int add_response_contributions(fp_t***, fp_t**, fp_t **, fp_t ***, fp_t ***, fp_t, fp_t, fp_t, fp_t, fp_t, fp_t ***, fp_t *****, fp_t *****){
    return 0;
  }
   virtual int add_response_contributions_new(fp_t***, fp_t**, fp_t **, fp_t ***, fp_t ***, fp_t, fp_t, fp_t, fp_t, fp_t, fp_t ***, fp_t *****, fp_t *****){
    return 0;
  }

  virtual int add_response_contributions_taugrid(fp_t***, fp_t ***, fp_t *****, fp_t**, fp_t **, fp_t ***, fp_t ***, fp_t, fp_t, fp_t, fp_t, fp_t, fp_t ***, fp_t *****, fp_t *****){
    return 0;
  }

  virtual int add_pops_to_response(int, int){
    return 0;
  }
  virtual int subtract_pops_from_response(int, int){
    return 0;
  }
  virtual int divide_responses_by_step(int, int, fp_t){
    return 0;
  }
  virtual int print_population_responses(const char*, int){
    return 0;
  }
   virtual int print_population_responses(const char*, int, int){
    return 0;
  }

  virtual void add_to_ion_responses(int,int, fp_t sign){};
  virtual void divide_ion_responses(int,int,fp_t){};

  virtual fp_t get_population_response(int param, int x3k, int x1i, int x2i, int x3i, int z){
    return 0;
  }
  
  virtual void prof_setup(void){
  }
  virtual void prof_clear(void){
    
  }
  virtual fp_t pops(atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t);
  virtual fp_t pops_broyden(atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t);
  virtual void add(fp_t***,fp_t***,fp_t***,fp_t***,fp_t***,fp_t***,fp_t,fp_t,fp_t,fp_t);
  virtual void add(fp_t***, fp_t***, fp_t ***, fp_t, fp_t, fp_t){
  }
  virtual void add_to_radiation_tensor(fp_t***, fp_t***, fp_t ***, fp_t, fp_t, fp_t, fp_t){
  }
  virtual void add_to_radiation_tensor_perturbation(fp_t***, fp_t**,fp_t **, fp_t ***, fp_t ***, fp_t, fp_t, fp_t, fp_t, fp_t, fp_t ***, fp_t *****, fp_t *****){
  }

  int set_parent_atmosphere(atmosphere * atm_in);
  int clear_parent_atmosphere();
  fp_t fetch_population(int x1i, int x2i, int x3i, int species, int z, int i);
  fp_t fetch_population(int x1i, int x2i, int x3i, int species, int z);
  fp_t fetch_temperature (int x1i, int x2i, int x3i);
  fp_t fetch_Ne(int x1i, int x2i, int x3i);
  fp_t fetch_vt(int x1i, int x2i, int x3i);
  fp_t * fetch_magnetic_field(int x1i, int x2i, int x3i);
  fp_t fetch_Nt(int x1i, int x2i, int x3i);

// *************************************************************************************
// * routines needed for the new combined chemical/ionization equilibrium calculations *
// *************************************************************************************
// atomic properties
  virtual fp_t abundance(void){ return -1.0; };
  virtual void ionfrc(fp_t,fp_t,fp_t*&,fp_t*&,int&){};
// molecular properties
  virtual fp_t dissoc(fp_t,fp_t){ return 0.0; } // default: not a molecule
  virtual ion *components(int &nc_o){ return 0; }
  virtual void pupdate(fp_t,fp_t*,int,int32_t,int32_t,int32_t){
  };
  virtual void randomize_populations(fp_t){
  };
//
  virtual void print_populations(){
  };
  virtual void info(void){};


  // New functions: for debugging and some other new stuff. (Milic 2014/15)
  virtual fp_t get_total_pop(int x1i, int x2i, int x3i){
    //printf("This is only atmol, so not total pop.\n");
    return 0;
  }
  virtual fp_t get_pop(int x1i, int x2i, int x3i, int){
    return 0; // As this is the atmol, not atom, nor molecule (nor H- molecule as a matter of fact), do not return anything
  }
  virtual fp_t get_pop(int x1i, int x2i, int x3i, int, int){
    return 0; // As this is the atmol, not atom, nor molecule (nor H- molecule as a matter of fact), do not return anything
  }

  virtual fp_t get_J(int, int, int, int){
    return 0;
  }

  virtual void print_radiation_field_tensor(){

  }

  virtual fp_t get_L(int, int, int, int){
    return 0;
  }
  virtual fp_t get_norm(int, int, int, int){
    return 0;
  }

  virtual fp_t get_active_pop(int, int, int){
    return 0;
  }
  virtual fp_t get_partf(int z, fp_t T, fp_t){
    return 1.0;
  };

  // ------------------------------------------------------------------------------------------------------------------------------------------------------
  // Some response function related stuff:
  virtual void compute_nlte_population_responses(){}
  virtual void compute_lte_population_responses(){}
  virtual void compute_lte_population_responses_analytical(fp_t ***, fp_t ***){}

  virtual int check_if_nlte(){
    return 0; // By default it is not nlte, because it is not an atom.
  }

};

#include "mol/molcfg.h"
#include "atom/atomcfg.h"

atmol *atmol_new(atmcfg*,io_class&);
atmol *atmol_new(molcfg*,atmol**,int,io_class&);
atmol *atmol_new(uint08_t*,int32_t&,uint08_t,atmol**,int,io_class&);

atmol **append(atmol *q,atmol **p,int &len);

#endif              // __ATMOL_H__
