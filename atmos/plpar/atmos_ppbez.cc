#include <math.h>
#include "types.h"
#include "io.h"
#include "acfg.h"
#include "atmos.h"

#include "atmos_pp.h"
#include "atmos_ppbez.h"

atmos_ppbez::atmos_ppbez(acfg *cfg,io_class &io_in):atmos_pp(cfg,io_in)
{
  rtstype=ATMOS_RTS_BEZ;
  io.msg(IOL_DEB1,"atmos_ppbez::atmos_ppbez: using Bezier Spline interpolation\n");
}

atmos_ppbez::atmos_ppbez(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):atmos_pp(buf,offs,do_swap,io_in)
{
  rtstype=ATMOS_RTS_BEZ;
  io.msg(IOL_DEB1,"atmos_ppbez::atmos_ppbez: using Bezier Spline interpolation\n");
}

atmos_ppbez::~atmos_ppbez(void)
{
}

