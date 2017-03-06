#ifndef  __HYDROGEN_H__ // __HYDROGEN_H__ 
#define  __HYDROGEN_H__

#include "types.h"
#include "io.h"

#include "atomcfg.h"
#include "atom.h"

#include "hydrogen.h"

class hydrogen:public atom{
protected:
  virtual fp_t cbb(uint08_t,uint16_t,uint16_t,atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t);
  virtual fp_t cbf(uint08_t,uint16_t,atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t);
//  virtual fp_t cfb(uint08_t,uint16_t,atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t);
public:
  hydrogen(atmcfg*,io_class&);
  hydrogen(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~hydrogen(void);
};

#endif                  // __HYDROGEN_H__
