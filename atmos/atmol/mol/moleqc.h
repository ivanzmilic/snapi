#ifndef  __MOLEQC_H__ // __MOLEQC_H__
#define  __MOLEQC_H__

#include "types.h"
#include "io.h"

#include "molcfg.h"

#define EQC_TYPE_NONE   0
#define EQC_TYPE_CNST   1
#define EQC_TYPE_POLY   2
#define EQC_TYPE_HMIN   3
//#define EQC_TYPE_SAHA   4

class eqc{                 // dissociation coefficient
protected:
  uint08_t type;
public:
  eqc(uint08_t*,int32_t&,uint08_t,io_class&);
  eqc(uint08_t type_in){ type=type_in; }
  eqc(void){ type=EQC_TYPE_NONE; }
  virtual ~eqc(void){}
  virtual fp_t K(fp_t,fp_t){ return 0.0; }
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
};

class eqc_cnst:public eqc{ // constant dissociation coefficient
protected:
  fp_t coef;
public:
  eqc_cnst(uint08_t*,int32_t&,uint08_t,io_class&);
  eqc_cnst(fp_t coef_in,io_class&):eqc(EQC_TYPE_CNST){ coef=coef_in; }
  virtual ~eqc_cnst(void){}
  virtual fp_t K(fp_t,fp_t){ return coef; }
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
};

class eqc_poly:public eqc{ // polynomial dissociation coefficient
protected:
  fp_t *cfs;
  int08_t nc;
public:
  eqc_poly(fp_t*,uint08_t,io_class&);
  eqc_poly(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~eqc_poly(void);
  virtual fp_t K(fp_t,fp_t);
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
};

class eqc_hmin:public eqc{ // H-minus dissociation coefficient
public:
  eqc_hmin(uint08_t*,int32_t&,uint08_t,io_class&);
  eqc_hmin(io_class&):eqc(EQC_TYPE_HMIN){}
  virtual ~eqc_hmin(void){}
  virtual fp_t K(fp_t,fp_t);
};

/*
class eqc_saha:public eqc{ // saha type dissociation coefficient

};
*/
class eqc *eqc_new(struct eqcfg *cfg,io_class &io);
class eqc *eqc_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io);

#endif  // __MOLEQC_H__
