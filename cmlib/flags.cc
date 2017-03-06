#include <string.h>

#include "types.h"
#include "io.h"
#include "pack.h"

#include "flags.h"

reg::reg(uint64_t mask_in)
{
  data=0;
  mask=mask_in;
}

int32_t reg::size(void)
{
  return 2*sizeof(uint64_t);
}

int32_t reg::pack(uint08_t *p,io_class &io)
{
  int offs=::pack(p,(uint08_t*)&data,0,sizeof(uint64_t)-1);
  offs+=::pack(p+offs,(uint08_t*)&mask,0,sizeof(uint64_t)-1);
  return offs;
}

int32_t reg::unpack(uint08_t *p,io_class &io)
{
  int offs=::unpack(p,(uint08_t*)&data,0,sizeof(uint64_t)-1);
  offs+=::unpack(p+offs,(uint08_t*)&mask,0,sizeof(uint64_t)-1);
  return offs;
}

void reg::set_mask(uint64_t mask_in)
{
  mask=mask_in;
}

uint64_t reg::get_mask(void)
{
  return mask;
}

void reg::clear(void)
{
  data&=~mask;
}

void reg::clear(uint64_t fl)
{
  data&=~fl;
}

uint08_t reg::is_clr(uint64_t msk)
{
  return (data&msk)==0;
}

void reg::set(void)
{
  data|=mask;
}

void reg::set(uint64_t fl)
{
  data|=fl;
}

uint08_t reg::is_set(uint64_t msk)
{
  return (data&msk)==msk;
}
