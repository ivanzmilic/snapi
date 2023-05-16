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

  int of; //opacity fudge yes/no
  char *of_filename; // opacity fudge filename
  int N_of; // number of wavelenghts for opacity fudge
  int conserve_charge; // How to treat ne (0 - LTE, 1 - simple charge conservation...)
  int tau_grid; // RT grid (0 - height (default), 1 - logtau_500)
  int use_atm_lvls; // Whether to save / load atomic levels. Only saving supported so far.
//
  acfg(void);
  acfg(char *adata,io_class&);
  ~acfg(void);
};

#endif            // __ACFG_H__
