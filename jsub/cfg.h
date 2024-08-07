#ifndef __CFG_H__ // __CFG_H__
#define __CFG_H__ 

#include "types.h"
#include "io.h"
#include "cmdln.h"

#include "acfg.h"

struct gcfg;

struct ocfg{                  // configuration per object ("color")
  
  char *id;
  
  fp_t az,el;
  
  fp_t *lambda;
  int nlambda;
  
  fp_t spectral_broadening; // If you want a simple gaussian broadening, just use fixed width
  fp_t *spsf; // Spectral psf, presumed given on the equidistant grid with the spacing same as the wavelength
  int n_spsf; // number of wavelength points for the spectral psf 
  
  char *name;
  
  int to_invert;
  int xl,xh,yl,yh,ll,lh;
  fp_t * weight; // what are the weights for the fitting
  int return_model; // After inverting, should we return model? 
  int return_atmos; // After inverting, should we return full atmosphere?
  fp_t scattered_light; // percentage of gray, scattered light in the image.
  fp_t obs_qs; // observed quiet sun intensity (in observed units)
  fp_t synth_qs; // quiet sun resulting from falc (or other) model
  int  no_iterations; // maximum number of iterations for the inversion
  fp_t starting_lambda; // starting value of lambda parameter
  fp_t stopping_chisq; // stop fitting when this chisq is reached (good for multiple cycles)
  fp_t * w_stokes; // weights for different stokes parameters

//
  ocfg(char *odata,struct gcfg &gdata,io_class&);
  ~ocfg(void);
};

struct parcfg{ // configuration for each of the parameter, nodes and values
  char *id;
  int n; // number of nodes
  fp_t *tau;
  fp_t *value;
  int reg_type; // type of regularization
  fp_t reg_alpha; // coefficient of regularization

  parcfg(char *pardata,io_class&);
  ~parcfg(void);

};

struct mcfg{                  // configuration for model, which consists from nodes
  char *id;
  parcfg **par;
  int np;
  int read_from_file; // whetwher to read cube of models from file. Needed for multi-step approach.
  char *filename;
  fp_t tau_min; // Upper tau (typically -4 to -8)
  fp_t tau_max; // Lower tau (typically 0 to 2)
  mcfg(char *mdata,struct gcfg &gdata,io_class&);
  ~mcfg(void);
};

struct gcfg{                  // global configuration
  struct acfg **atm;
  struct ocfg **obs;
  struct mcfg **mod;
  int no,na,nm;
//
  gcfg(cmdln&,io_class&);
  ~gcfg(void);
};

#endif            // __CFG_H__
