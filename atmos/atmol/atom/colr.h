#ifndef  __COLR_H__ // __COLR_H__
#define  __COLR_H__

#include "types.h"
#include "io.h"

#include "atomcfg.h"

#define COLR_TYPE_NONE   0
#define COLR_TYPE_TAB    1


class colr{     // generic collision rate from tables
protected:
  int08_t type;
public:
  colr(uint08_t*,int32_t&,uint08_t,io_class&);
  colr(int08_t);
  colr(void);
  virtual ~colr(void){};
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
//
  virtual fp_t C(fp_t,fp_t){ return 0.0; }
  int08_t get_type(){return type;};
};


class tbcr:public colr{ // tabulated collision rates
  int n;
  fp_t *t,*v;
  fp_t interpol(fp_t);
public:
  tbcr(uint08_t*,int32_t&,uint08_t,io_class&);
  tbcr(fp_t*,fp_t*,int,io_class&);
  virtual ~tbcr(void);
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
//
  virtual fp_t C(fp_t,fp_t);
};


class colr *cr_new(uint08_t*,int32_t&,uint08_t,io_class&);
class colr *cr_new(fp_t,fp_t,fp_t,fp_t,struct crcfg*,io_class&);

#endif              // __COLR_H__
