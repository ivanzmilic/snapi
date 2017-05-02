#ifndef  __BFCS_H__ // __BFCS_H__
#define  __BFCS_H__

#include "types.h"
#include "io.h"

#include "atomcfg.h"

#define BFCS_TYPE_NONE   0
#define BFCS_TYPE_HYD    1
#define BFCS_TYPE_HEL    2
#define BFCS_TYPE_TAB    3
#define BFCS_TYPE_PHY    4


class bfcs{ // bound-free crossection
protected:
  int08_t type;
public:
  bfcs(uint08_t*,int32_t&,uint08_t,io_class&);
  bfcs(int08_t);
  bfcs(void);
  virtual ~bfcs(void){};
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
//
  virtual fp_t *U(fp_t*,int32_t){ return 0; }
  virtual fp_t U(fp_t){ return 0.0; }
//
  virtual fp_t *getlambda(int&){ return 0; }
  virtual fp_t getlambda_crit(){return 0;}
};

class hybf:public bfcs{ // hydrogenic
  fp_t *g,*l,l_b,a_b;
  int n;
public:
  hybf(uint08_t*,int32_t&,uint08_t,io_class&);
  hybf(fp_t,fp_t,fp_t,fp_t,fp_t*,fp_t*,int,io_class&);
  virtual ~hybf(void);
//
  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
//
  virtual fp_t *U(fp_t*,int32_t);
  virtual fp_t U(fp_t);
//
  virtual fp_t *getlambda(int&);
  virtual fp_t getlambda_crit();
};

class tabbf:public bfcs{ // tabulated
  fp_t * lambda_tabulated; 
  fp_t * bf_tabulated; // We have to deviate here a bit from previous formalism, but what can we do :/
  int n; // number of points for interpolateion
public:
  tabbf(uint08_t*,int32_t&,uint08_t,io_class&); // We can decipher these two from what is done
  tabbf(fp_t,fp_t,fp_t,fp_t,fp_t*,fp_t*,int,io_class&);
  virtual ~tabbf(void);

  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);

  virtual fp_t U(fp_t);
  virtual fp_t *getlambda(int&);
  virtual fp_t getlambda_crit();
};

class phbf:public bfcs{ // pseudo-hydrogenic, only critical cross-section is given
  fp_t lambda_crit;
  fp_t bf_crit; 
public:
  phbf(uint08_t*,int32_t&,uint08_t,io_class&); // We can decipher these two from what is done
  phbf(fp_t,fp_t,fp_t,fp_t,fp_t*,fp_t*,int,io_class&);
  
  virtual ~phbf(void);

  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);

  virtual fp_t U(fp_t);
  virtual fp_t *getlambda(int&);
  virtual fp_t getlambda_crit();
};

class bfcs *bf_new(uint08_t*,int32_t&,uint08_t,io_class&);
class bfcs *bf_new(fp_t,fp_t,fp_t,fp_t,struct bfcfg*,io_class&);

#endif              // __BFCS_H__
