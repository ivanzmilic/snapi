#include <math.h>
#include "const.h"
#include "types.h"
#include "io.h"

#include "atomcfg.h"

#include "helium.h"

helium::helium(atmcfg *cfg,io_class &io_in):atom(cfg,io_in)
{
//  io_in.msg(IOL_INFO,"helium::helium: %s[%s]\n",cfg->name,cfg->id);
}

helium::helium(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):atom(buf,offs,do_swap,io_in)
{
//  io_in.msg(IOL_INFO,"helium::helium: %lX %d\n",buf,offs);
}

helium::~helium(void)
{
}

fp_t *helium::boundfree(fp_t *lambda,int32_t nlambda)
{
  return 0;
}
