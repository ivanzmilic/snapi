#include <string.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "pack.h"

#include "grid.h"

grid::grid(io_class &io_in):node(0,io_in)
{
  x1l=x1h=x2l=x2h=x3l=x3h=0;
//
  (x1=new fp_t [x1h-x1l+1]-x1l)[x1l]=0.0;
  (x2=new fp_t [x2h-x2l+1]-x2l)[x2l]=0.0;
  (x3=new fp_t [x3h-x3l+1]-x3l)[x3l]=0.0;
}

grid::grid(uint08_t ndim_in,int32_t *dim_in,io_class &io_in):node(ndim_in,io_in)
{
//  memcpy((dim=new int32_t [ndim]-1)+1,dim_in+1,ndim*sizeof(int32_t));
//  x=ft2dim(1,ndim,1,dim);
}

grid::grid(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):node(buf,offs,do_swap,io_in)
{
  printf("grid::grid : can I even make a grid from a buffer?\n");
  offs+=unpack(buf+offs,do_swap,io_in);
  printf("grid::grid : yes I can!\n");
}

grid::~grid(void)
{
  if(x1) delete[] (x1+x1l);
  if(x2) delete[] (x2+x2l);
  if(x3) delete[] (x3+x3l);
}


int08_t grid::resize(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
{
  if(x1) delete[] (x1+x1l);
  if(x2) delete[] (x2+x2l);
  if(x3) delete[] (x3+x3l);
//
  x1l=x1l_in;
  x1h=x1h_in;
  x2l=x2l_in;
  x2h=x2h_in;
  x3l=x3l_in;
  x3h=x3h_in;
//
  x1=new fp_t [x1h-x1l+1]-x1l;
  x2=new fp_t [x2h-x2l+1]-x2l;
  x3=new fp_t [x3h-x3l+1]-x3l;
  for(int i=x1l;i<=x1h;++i) x1[i]=0.0;
  for(int i=x2l;i<=x2h;++i) x2[i]=0.0;
  for(int i=x3l;i<=x3h;++i) x3[i]=0.0;
//
  return 0;
}

int08_t grid::resize(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in,fp_t minx1,fp_t maxx1,fp_t minx2,fp_t maxx2,fp_t minx3,fp_t maxx3)
{
  if(x1) delete[] (x1+x1l);
  if(x2) delete[] (x2+x2l);
  if(x3) delete[] (x3+x3l);
//
  x1l=x1l_in;
  x1h=x1h_in;
  x2l=x2l_in;
  x2h=x2h_in;
  x3l=x3l_in;
  x3h=x3h_in;
//
  x1=new fp_t [x1h-x1l+1]-x1l;
  x2=new fp_t [x2h-x2l+1]-x2l;
  x3=new fp_t [x3h-x3l+1]-x3l;
  for(int i=x1l;i<=x1h;++i) x1[i]=minx1+(maxx1-minx1)*(fp_t)(i-x1l)/(fp_t)(x1h-x1l);
  for(int i=x2l;i<=x2h;++i) x2[i]=minx2+(maxx2-minx2)*(fp_t)(i-x2l)/(fp_t)(x2h-x2l);
  for(int i=x3l;i<=x3h;++i) x3[i]=minx3+(maxx3-minx3)*(fp_t)(i-x3l)/(fp_t)(x3h-x3l);
//
  return 0;
}

int32_t grid::size(void)
{
  int32_t sz=node::size();
// add local size
  sz+=6*sizeof(int32_t); // x?l,x?h
  sz+=(x1h-x1l+1)*sizeof(fp_t); // x1
  sz+=(x2h-x2l+1)*sizeof(fp_t); // x2
  sz+=(x3h-x3l+1)*sizeof(fp_t); // x3
// 
  return sz;
}

int32_t grid::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
  int32_t offs=node::pack(buf,do_swap,io_in);
// pack local stuff
  offs+=::pack(buf+offs,x1l,do_swap);
  offs+=::pack(buf+offs,x1h,do_swap);
  offs+=::pack(buf+offs,x2l,do_swap);
  offs+=::pack(buf+offs,x2h,do_swap);
  offs+=::pack(buf+offs,x3l,do_swap);
  offs+=::pack(buf+offs,x3h,do_swap);
//
  offs+=::pack(buf+offs,x1,x1l,x1h,do_swap);
  offs+=::pack(buf+offs,x2,x2l,x2h,do_swap);
  offs+=::pack(buf+offs,x3,x3l,x3h,do_swap);
//
  return offs;
}

int32_t grid::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
// only unpack local stuff  - nodes unpacked automatically in the contructor - this makes a problem somehow?

  int32_t offs=::unpack(buf,x1l,do_swap);
  offs+=::unpack(buf+offs,x1h,do_swap);
  offs+=::unpack(buf+offs,x2l,do_swap);
  offs+=::unpack(buf+offs,x2h,do_swap);
  offs+=::unpack(buf+offs,x3l,do_swap);
  offs+=::unpack(buf+offs,x3h,do_swap);
  printf("Dimensions = %d %d %d %d %d %d\n",x1l,x1h,x2l,x2h,x3l,x3h);
//
  offs+=::unpack(buf+offs,x1=new fp_t [x1h-x1l+1]-x1l,x1l,x1h,do_swap);
  offs+=::unpack(buf+offs,x2=new fp_t [x2h-x2l+1]-x2l,x2l,x2h,do_swap);
  offs+=::unpack(buf+offs,x3=new fp_t [x3h-x3l+1]-x3l,x3l,x3h,do_swap);
//
  return offs;
}
