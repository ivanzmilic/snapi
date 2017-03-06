#ifndef __FLAGS_H__   // __FLAGS_H__
#define __FLAGS_H__ 

#include "types.h"
#include "io.h"

class reg{
  uint64_t data,mask;
public:
  reg(uint64_t);
//
  int32_t size(void);
  int32_t pack(uint08_t*,io_class&);
  int32_t unpack(uint08_t*,io_class&);
//
  void set_mask(uint64_t);
  uint64_t get_mask(void);
//
  void clear(void);
  void clear(uint64_t);
  uint08_t is_clr(uint64_t msk);
  void set(void);
  void set(uint64_t);
  uint08_t is_set(uint64_t msk);
};

#endif                // __FLAGS_H__
