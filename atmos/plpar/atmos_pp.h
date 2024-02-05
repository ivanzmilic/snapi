#ifndef __ATMOS_PP_H__  // __ATMOS_PP_H__
#define __ATMOS_PP_H__

//
// ***********************************
// * Plane parallel 1D atmosphere    *
// ***********************************
//

#include "types.h"
#include "io.h"
#include "acfg.h"
#include "atmos.h"

class atmos_pp:public atmosphere{
protected:
  virtual fp_t *anglesetup(fp_t*&,fp_t*&,int&);
  virtual int formal(fp_t *, fp_t***,fp_t***,fp_t***,fp_t***,fp_t,fp_t, fp_t);
  virtual int formal_with_lambda_operator(fp_t*, fp_t***,fp_t**,fp_t***,fp_t***,fp_t,fp_t, fp_t);
  virtual int formal (fp_t*, fp_t ****S,fp_t ***L,fp_t *****op,fp_t ****em,fp_t t,fp_t p, fp_t boundary); // Formal for polarized
  
  virtual int formal_pert_numerical(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  virtual int formal_pert_analytical(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  virtual int formal_pert_jcdti(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
public:
  atmos_pp(acfg*,io_class&);
  atmos_pp(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~atmos_pp(void);

  virtual fp_t * compute_full_lambda_operator(fp_t ***, int, int){
    return 0;
  };
  virtual int optical_depth_scale(fp_t***,fp_t***,fp_t,fp_t);
  virtual int build_from_nodes(model*);
  virtual int interpolate_from_nodes(model*);
  virtual int enforce_hequilibrium();
  virtual fp_t ** calculate_dN_dT(); 
//
//  virtual int08_t init(const char*,io_class*);
//  virtual fp_t *obs(fp_t,fp_t,fp_t*,int32_t);
};

atmos_pp *atmos_pp_new(acfg *cfg,io_class &io_in);
atmos_pp *atmos_pp_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in);

#endif // __ATMOS_PP_H__
