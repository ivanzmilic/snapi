#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "types.h"
#include "io.h"
#include "pack.h"

#include "molcfg.h"

#include "const.h"

#include "moleqc.h"

class eqc *eqc_new(struct eqcfg *cfg,io_class &io)
{
  if(!cfg) return new eqc;
  if(!strcmp(cfg->type,"CONST")) return new eqc_cnst(cfg->K,io);
  if(!strcmp(cfg->type,"POLY")) return new eqc_poly(cfg->cfs,cfg->ncfs,io);
  if(!strcmp(cfg->type,"HMIN")) return new eqc_hmin(io);
//  if(!strcmp(cfg->eqctype,"SAHA")) return new eqc_saha(cfg->K[0]);
  return new eqc;
}

class eqc *eqc_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  uint08_t type;
  ::unpack(buf+offs,type);
  switch(type){
    case(EQC_TYPE_CNST): return new eqc_cnst(buf,offs,do_swap,io);
    case(EQC_TYPE_POLY): return new eqc_poly(buf,offs,do_swap,io);
    case(EQC_TYPE_HMIN): return new eqc_hmin(buf,offs,do_swap,io);
//    case(EQC_TYPE_SAHA): return new eqc_saha(buf,offs,do_swap,io);
  }
  return new eqc(buf,offs,do_swap,io);
}

//
// chemical equilibrium constant class
//

eqc::eqc(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

int32_t eqc::size(io_class &io)
{
  return sizeof(uint08_t);
}

int32_t eqc::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::pack(buf,type);
  return offs;
}

int32_t eqc::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::unpack(buf,type);
  return offs;
}

//
// chemical equilibrium constant class: constant
//

eqc_cnst::eqc_cnst(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):eqc(buf,offs,do_swap,io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

int32_t eqc_cnst::size(io_class &io)
{
  int32_t sz=eqc::size(io);
  sz+=sizeof(fp_t);
  return sz;
}

int32_t eqc_cnst::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=eqc::pack(buf,do_swap,io);
  offs+=::pack(buf+offs,coef,do_swap);
  return offs;
}

int32_t eqc_cnst::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only local stuff
  int32_t offs=::unpack(buf,coef,do_swap);
  return offs;
}

//
// chemical equilibrium constant class: polynomial description
//

eqc_poly::eqc_poly(fp_t *cfs_in,uint08_t nc_in,io_class &io):eqc(EQC_TYPE_POLY)
{
  cfs=new fp_t [nc=nc_in];
  memcpy(cfs,cfs_in,nc*sizeof(fp_t));
}

eqc_poly::eqc_poly(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):eqc(buf,offs,do_swap,io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

eqc_poly::~eqc_poly(void)
{
  if(cfs) delete[] cfs;
}

fp_t eqc_poly::K(fp_t T,fp_t)
{
  fp_t x=5040.0/T;
  fp_t sum=cfs[nc-1];
  for(int i=nc-2;i>=0;--i) sum=cfs[i]+x*sum;
  return exp(2.302585093*sum); // 10^sum
}

int32_t eqc_poly::size(io_class &io)
{
  int32_t sz=eqc::size(io);
  sz+=sizeof(uint08_t);
  sz+=nc*sizeof(fp_t);
  return sz;
}

int32_t eqc_poly::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=eqc::pack(buf,do_swap,io);
  offs+=::pack(buf+offs,nc);
  offs+=::pack(buf+offs,cfs,0,nc-1,do_swap);
  return offs;
}

int32_t eqc_poly::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only local stuff
  int32_t offs=::unpack(buf,nc);
  offs+=::unpack(buf+offs,cfs=new fp_t [nc],0,nc-1,do_swap);
  return offs;
}

eqc_hmin::eqc_hmin(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):eqc(buf,offs,do_swap,io)
{
}

fp_t eqc_hmin::K(fp_t T,fp_t ne)
{
/*  fp_t x=5040.0/T;
  fp_t dsum=2.0*ne;                  // Hydrogen only approximation
  fp_t D=sqrt((e0*k*T)/(e*e*dsum));  // SI version of Debeye length, see also Mihalas 2ed. 293
  fp_t dxi=e*e/D;                    // Mihalas 2ed. 295
//
  fp_t xx=1.202E9/(x*x*sqrt(x));
  return 2.0*xx*exp(-2.302585*(0.754-2.0*dxi)*x); // H-
*/
  // Note by Ivan: 07/07/2017: This number is actually not so important since we 
  // will be computing opacity of H- from number densities of H and e- anyway! 
  fp_t x=5040.0/T;
  fp_t sum=0.1249+(2.5*log(T)/2.302585)-0.747*x;

  fp_t old_value = exp(2.302585093*sum);

  return exp(2.302585093*sum); // 10^sum
}
