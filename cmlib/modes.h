#ifndef __MODES_H__ // __MODES_H__
#define __MODES_H__ 

#include "types.h"
#include "io.h"

struct klmc{
  fp_t *c,v;
  int *m,nm;
};

struct modes{
  float64_t ****mode,**pupil,area;
  float64_t **covar;
  int32_t nph,zmin,zmax,kmin,kmax,mmin,mmax;
  modes(struct klmc*,fp_t,fp_t,int,int,int,int*,int,int*,int**,int**,int,int,fp_t,fp_t,fp_t**,io_class&);
  modes(byte*,byte,io_class&);
  ~modes(void);
  byte *compress(int&,int,byte,io_class &io);
};

struct klmc *kl_cfg(int kl_min_mode,int kl_max_mode);
void kl_uncfg(struct klmc *cfs,int kl_min_mode,int kl_max_mode);

#endif              // __MODES_H__
