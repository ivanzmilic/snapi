#include <math.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "pack.h"
#include "uts.h"
#include "const.h"

#include "atomcfg.h"

#include "colr.h"

class colr *cr_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  int08_t type;
  ::unpack(buf+offs,type);
  switch(type){
    case(COLR_TYPE_TAB): return new tbcr(buf,offs,do_swap,io);
  }
  return new colr(buf,offs,do_swap,io);
}

class colr *cr_new(struct crcfg *cfg,io_class &io_in)
{
  if(!cfg) return new colr(COLR_TYPE_NONE);
  if(!strcmp(cfg->crtype,"TAB")) return new tbcr(cfg->t,cfg->v,cfg->n,io_in);
  return new colr(COLR_TYPE_NONE);
}

colr::colr(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

colr::colr(int08_t type_in)
{
  type=type_in;
}

int32_t colr::size(io_class &io)
{
  return sizeof(int08_t);
}

int32_t colr::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::pack(buf,type);
  return offs;
}

int32_t colr::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::unpack(buf,type);
  return offs;
}

// *******************************
// * Tabulated collisional rates *
// *******************************

tbcr::tbcr(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):colr(buf,offs,do_swap,io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

tbcr::tbcr(fp_t *t_in,fp_t *v_in,int n_in,io_class &io):colr(COLR_TYPE_TAB)
{
  n=n_in;
  t=new fp_t [n];
  v=new fp_t [n];
  memcpy(t,t_in,n*sizeof(fp_t));
  memcpy(v,v_in,n*sizeof(fp_t));
//
}

tbcr::~tbcr(void)
{
  if(t) delete[] t;
  if(v) delete[] v;
}

int32_t tbcr::size(io_class &io)
{
  int sz=colr::size(io);
  sz+=sizeof(int32_t);  // n
  sz+=2*n*sizeof(fp_t); // l,v
  return sz;
}

int32_t tbcr::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=colr::pack(buf,do_swap,io);
// local stuff
  offs+=::pack(buf+offs,n,do_swap);
  offs+=::pack(buf+offs,t,0,n-1,do_swap);
  offs+=::pack(buf+offs,v,0,n-1,do_swap);
  return offs;
}

int32_t tbcr::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only unpack local stuff
  int32_t offs=::unpack(buf,n,do_swap);
  offs+=::unpack(buf+offs,t=new fp_t [n],0,n-1,do_swap);
  offs+=::unpack(buf+offs,v=new fp_t [n],0,n-1,do_swap);
//
  return offs;
}

fp_t tbcr::interpol(fp_t tmp)
{ // simple linear interpolation
  for(int i=1;i<n;++i) if(t[i]>tmp) return (v[i]*(tmp-t[i-1])+v[i-1]*(t[i]-tmp))/(t[i]-t[i-1]);
  return v[n-1];
}

fp_t tbcr::C(fp_t Temp,fp_t ne)
{
  return 0.0;
}
