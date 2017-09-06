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
  char *name;
  int nlambda;
  int to_invert;
  int xl,xh,yl,yh,ll,lh;
  fp_t * weight; // what are the weights for the fitting
  int return_model; // After inverting, should we return model? 
  int return_atmos; // After inverting, should we return full atmosphere?
  fp_t scattered_light; // percentage of gray, scattered light in the image.
  fp_t spectral_broadening; // spectral broadening due to the instrument
  fp_t obs_qs; // observed quiet sun intensity (in observed units)
  fp_t synth_qs; // quiet sun resulting from falc (or other) model

//
  ocfg(char *odata,struct gcfg &gdata,io_class&);
  ~ocfg(void);
};

struct parcfg{ // configuration for each of the parameter, nodes and values
  char *id;
  int n; // number of nodes
  fp_t *tau;
  fp_t *value;

  parcfg(char *pardata,io_class&);
  ~parcfg(void);

};

struct mcfg{                  // configuration for model, which consists from nodes
  char *id;
  parcfg **par;
  int np;
  int read_from_file; // whetwher to read cube of models from file. Needed for multi-step approach.
  char *filename;
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
