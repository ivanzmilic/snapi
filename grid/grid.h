#ifndef __GRID_H__  // __GRID_H__
#define __GRID_H__

#include "types.h"
#include "io.h"

#include "node.h"

class grid:public node{
protected:
  int32_t x1l,x1h,x2l,x2h,x3l,x3h;  // dimensions
  fp_t *x1,*x2,*x3;                 // coordinates
//  int32_t *dim;   // size of each dimension
//  fp_t **x;       // coordinates
  virtual int08_t resize(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t);
  virtual int08_t resize(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,fp_t,fp_t,fp_t,fp_t,fp_t,fp_t);
public:
  grid(io_class&);
  grid(uint08_t,int32_t*,io_class&);
  grid(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~grid(void);
//
  virtual int32_t size(void);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t *,uint08_t,io_class&);
  virtual grid* extract_grid(int i, int j, io_class&);
//
};

#endif              // __GRID_H__
