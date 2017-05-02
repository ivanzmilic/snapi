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
//
  ocfg(char *odata,struct gcfg &gdata,io_class&);
  ~ocfg(void);
};

struct parcfg{ // configuration for each of the parameter, nodes and values
  char *id;
  int n;
  fp_t *tau;
  fp_t *value;

  parcfg(char *pardata,io_class&);
  ~parcfg(void);

};

struct mcfg{                  // configuration for model, which consists from nodes
  char *id;
  parcfg **par;
  int np;
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
