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
//
  ocfg(char *odata,struct gcfg &gdata,io_class&);
  ~ocfg(void);
};

struct gcfg{                  // global configuration
  struct acfg **atm;
  struct ocfg **obs;
  int no,na;
//
  gcfg(cmdln&,io_class&);
  ~gcfg(void);
};

#endif            // __CFG_H__
