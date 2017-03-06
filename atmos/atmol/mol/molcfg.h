#ifndef  __MOLCFG_H__  // __MOLCFG_H__
#define  __MOLCFG_H__

#include "types.h"
#include "io.h"

struct eqcfg{
 char *type;
// constant
 fp_t K;
// polynomial
 fp_t *cfs;
 int08_t ncfs;
// Saha
// fp_t ediss;
public:
  eqcfg(char*,const char*,const char*,io_class&);
  ~eqcfg(void);
};

struct molcfg{
  char *id,*name;
  char **cid;
  char **zid;
//
  eqcfg **K;
  int nk;
public:
  molcfg(char*,io_class&);
  ~molcfg(void);
};

#endif                 // __MOLCFG_H__
