#ifndef  __PARTF_H__ // __PARTF_H__
#define  __PARTF_H__

#include "types.h"
#include "io.h"

#include "atomcfg.h"

#define PF_TYPE_NONE   0
#define PF_TYPE_CONST  1
#define PF_TYPE_TRAV   2
#define PF_TYPE_IRWIN  3 

//
// MULTI OPTIONS:
//   PARTITION FUNCTIONS TO BE INTERPOLATED IN T
//   PARTITION FUNCTIONS FOLLOWING TRAVING ET AL., ABH. HAMB. VIII, 1 (1966) (QTRAV)
//   THE PARTITION FUNCTION IS CONSTANT
//   THE PARTITION FUNCTION IS NLTE
//

class pf{ // partition function
protected:
  int08_t type;
public:
  pf(uint08_t*,int32_t&,uint08_t,io_class&);
  pf(int08_t);
//  pf(void);
  virtual ~pf(void){};
  virtual fp_t U(fp_t,fp_t,io_class&){ return 1.0; }
  virtual fp_t dU(fp_t,fp_t,io_class&){ return 0.0; }
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
};

class cpf:public pf{ // partition function
  fp_t val;
public:
  cpf(uint08_t*,int32_t&,uint08_t,io_class&);
  cpf(fp_t);
  virtual ~cpf(void){};
  virtual fp_t U(fp_t,fp_t,io_class&){ return val; };
  virtual fp_t dU(fp_t,fp_t,io_class&){ return 0.0; };
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
};

class trav:public pf{ // Traving et. al. partition function
  fp_t z,g0,*g2,*ee,*ll;
  fp_t **a,**g;
  int n,*m;
  int08_t *asym;
  fp_t qas_b(fp_t,fp_t,fp_t);
  fp_t qas_f(fp_t,fp_t,fp_t);
  fp_t dqas_f(fp_t l,fp_t p,fp_t Dz);
public:
  trav(uint08_t*,int32_t&,uint08_t,io_class&);
  trav(int08_t,fp_t,tpfcfg**,int);
  ~trav(void);
  virtual fp_t U(fp_t,fp_t,io_class&);
  virtual fp_t dU(fp_t,fp_t,io_class&);
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
};

class irwin:public pf{ // Traving et. al. partition function
  fp_t z;
  fp_t *a;
  int n; // number of coefficients
public:
  irwin(uint08_t*,int32_t&,uint08_t,io_class&);
  irwin(int08_t,iwpfcfg*);
  ~irwin(void);
  virtual fp_t U(fp_t,fp_t,io_class&);
  virtual fp_t dU(fp_t,fp_t,io_class&);
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
};

class pf *pf_new(uint08_t*,int32_t&,uint08_t,io_class&);
class pf *pf_new(int08_t,struct pcfg*);

#endif              // __PARTF_H__
