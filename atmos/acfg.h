#ifndef __ACFG_H__ // __ACFG_H__
#define __ACFG_H__ 

#include "types.h"
#include "io.h"

#include "atmol/atom/atomcfg.h"
#include "atmol/mol/molcfg.h"

struct acfg{                  // atmosphere config
  char *geom,*rts;
//
  atmcfg **atm;
  int natm;
//
  molcfg **mol;
  int nmol;
//
  char *id;
  char *filename,*filetype;
  int nx,ny,nz;               // dimensions
  fp_t dx,dy,dz;              // grid parameters
//
  acfg(void);
  acfg(char *adata,io_class&);
  ~acfg(void);
};

#endif            // __ACFG_H__
