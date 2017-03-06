#ifndef  __HELIUM_H__ // __HELIUM_H__
#define  __HELIUM_H__

#include "types.h"
#include "io.h"

#include "atomcfg.h"
#include "atom.h"

#include "helium.h"

class helium:public atom{
protected:
  virtual fp_t *boundfree(fp_t*,int32_t);
public:
  helium(atmcfg*,io_class&);
  helium(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~helium(void);
};

#endif              // __HELIUM_H__
