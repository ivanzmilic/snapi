#ifndef __ATMOS_PPBEZ_H__  // __ATMOS_PP_H__
#define __ATMOS_PPBEZ_H__

//
// ***********************************
// * Plane parallel 1D atmosphere    *
// ***********************************
//

#include "types.h"
#include "io.h"
#include "acfg.h"
#include "atmos.h"
#include "atmos_pp.h"

class atmos_ppbez:public atmos_pp{
protected:
  virtual int formal(fp_t*, fp_t***,fp_t***,fp_t***,fp_t***,fp_t,fp_t, fp_t);
  virtual int formal_with_lambda_operator(fp_t *, fp_t***,fp_t**,fp_t***,fp_t***,fp_t,fp_t, fp_t);
  virtual int formal_with_responses(fp_t *, fp_t***, fp_t ***, fp_t**, fp_t**,fp_t***,fp_t***,fp_t,fp_t, fp_t);
  
  
  virtual int formal_pert_numerical(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  virtual int formal_pert_analytical(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);
  virtual int formal_pert_jcdti(fp_t ****, fp_t ***, fp_t ***, fp_t ****, fp_t ****, fp_t, fp_t, fp_t);

  // Polarized versions:
  virtual int formal (fp_t*, fp_t ****S,fp_t ***L,fp_t *****op,fp_t ****em,fp_t t,fp_t p, fp_t boundary); // Formal for polarized
  virtual int formal_pert_numerical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary);
  virtual int formal_pert_numerical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary, int N_parameters);
  virtual int formal_pert_analytical(fp_t ***** dS, fp_t ***** op, fp_t **** em, fp_t ****** op_pert, fp_t ***** em_pert, fp_t theta, fp_t phi, fp_t boundary);
  virtual int formal_with_responses_jcdti(fp_t *, fp_t ****,fp_t ***,fp_t *****,fp_t ****, fp_t *****, fp_t ****, fp_t ****, fp_t ,fp_t , fp_t );
  virtual int formal_with_responses_full(fp_t *, fp_t ****,fp_t ***,fp_t *****,fp_t ****, fp_t *****, fp_t ****, fp_t ****, fp_t ,fp_t , fp_t );
  
  virtual int compute_op_referent();
  virtual int compute_op_referent_derivative();
  virtual int compute_tau_referent();

public:
  atmos_ppbez(acfg *cfg,io_class &io_in);
  atmos_ppbez(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in);
  virtual fp_t * compute_full_lambda_operator(fp_t ***, int, int);
  virtual int optical_depth_scale(fp_t***,fp_t***,fp_t,fp_t);
  virtual ~atmos_ppbez(void);
};

#endif // __ATMOS_PPBEZ_H__
