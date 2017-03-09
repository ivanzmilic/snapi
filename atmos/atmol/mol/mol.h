#ifndef  __MOL_H__ // __MOL_H__
#define  __MOL_H__

#include "types.h"
#include "io.h"

#include "molcfg.h"

#include "atmol.h"

#include "moleqc.h"

class mol:public atmol{
protected:
  struct ion *comp; // references to the constituents
  uint08_t nc;
//
  eqc **K;
  uint08_t nk;
//
  fp_t ***N; // I have slightly changed this. This will ONLY be allocated if the molecule is of type h_minus_mol
public:
  mol(molcfg*,atmol**,int,io_class&);
  mol(uint08_t*,int32_t&,uint08_t,atmol**,int,io_class&);
  virtual ~mol(void);
// transport
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,atmol**,int,io_class&);
//  virtual fp_t *opacity(fp_t,fp_t,fp_t*,int32_t);
//  virtual fp_t *emissivity(fp_t,fp_t,fp_t*,int32_t);
  virtual fp_t dissoc(fp_t,fp_t);
  virtual ion *components(int &nc_o){ nc_o=nc; return comp; };

  virtual void popsetup(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t){};

  virtual void pupdate(fp_t,fp_t*,int); // <----- This is bad, different parameter list overloads, does not override the function
  virtual void pupdate(fp_t N_in,fp_t*,int,int32_t x1i,int32_t x2i,int32_t x3i){ // <------------- This one should override.

  }

  virtual void print_total_pop(int x1i, int x2i, int x3i){
    printf("This is mol, we still did not implement the function. \n");
  }

};

atmol *mol_new(molcfg*,atmol**,int,io_class&);
atmol *mol_new(uint32_t,uint08_t *buf,int32_t &offs,uint08_t do_swap,atmol **atm,int na,io_class &io_in);

// Below is edit by Milic, on 20/10/2014. We are adding a derived class, almost identical to molecule
// Essentially EVERYTHING is the same except that we will have to have different:
// - Constructors/destructors 
// - popsetup (remember that mol pcatically does not have it.)
// - pupdate
// - also opacity computation

class h_minus_mol:public mol{

protected:
  fp_t N_t; // total number of H- molecules
  fp_t **** dN; // Response of the total population to all the parameters

  int32_t x1l,x1h,x2l,x2h,x3l,x3h; // These are the same as in the story with atoms.

  // Constructors/destuctors, identical to ones from the mol class:

public:
  h_minus_mol(molcfg*,atmol**,int,io_class&);
  h_minus_mol(uint08_t*,int32_t&,uint08_t,atmol**,int,io_class&);
  
  // Now, other ones which we will have to override here:

  virtual void popsetup(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in);
  virtual void popclean(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in);

  // Some other ones, more of the technical nature:

  fp_t ***add(fp_t***,fp_t***,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  fp_t *****add(fp_t*****,fp_t*****,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);

  
  virtual void print_total_pop(int x1i, int x2i, int x3i){
    printf("This is h_minus_mol, population is: %e \n", N[x1i][x2i][x3i]);
  }

  void pupdate(fp_t N_in,fp_t*,int,int32_t x1i,int32_t x2i,int32_t x3i);

  // Now the ones which deal with opacity and emissivity:
  virtual fp_t ***opacity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t opacity_continuum(fp_t, fp_t, fp_t, int, int, int);
  virtual fp_t ** opacity_continuum_pert(fp_t, fp_t, fp_t, int, int, int);
  virtual fp_t ***emissivity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);

  virtual fp_t ***** opacity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t *****  emissivity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);

  virtual fp_t ***** freefree_op_pert(fp_t***,fp_t***,fp_t***,fp_t);
  virtual fp_t ***** boundfree_op_pert(fp_t***,fp_t);


  // --------------------------------------------------------------------------------------------
  virtual fp_t *****opacity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ****emissivity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t *******opacity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  virtual fp_t ******emissivity_vector_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t);
  
  // And specific kinds of opacity and emissivity:
  virtual fp_t ***freefree_op(fp_t***,fp_t***,fp_t***,fp_t);
  virtual fp_t ***freefree_em(fp_t***,fp_t***,fp_t***,fp_t);
  virtual fp_t ***boundfree_op(fp_t***,fp_t);
  virtual fp_t ***boundfree_em(fp_t***,fp_t);

  // -------------------------------------------------------------------------------------------
  // And there have to be some related to computing the responses:
  virtual void compute_lte_population_responses();
  virtual void compute_nlte_population_responses();
  virtual int responses_clear();
  virtual fp_t get_population_response(int parameter, int x3k, int x1i, int x2i, int x3i, int z);

};

// Now two identical ones as above, except they will return h_minus_mol, instead of mol

atmol *h_minus_mol_new(molcfg*,atmol**,int,io_class&);
atmol *h_minus_mol_new(uint32_t,uint08_t *buf,int32_t &offs,uint08_t do_swap,atmol **atm,int na,io_class &io_in);


#endif              // __MOL_H__
