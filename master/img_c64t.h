#ifndef __IMG_C64T_H__  // __IMG_C64T_H__
#define __IMG_C64T_H__

#include "types.h"
#include "io.h"

class image_c64t;

#include "img.h"
#include "fftw3.h"

class image_c64t:public image_t{
  complex_t **pic;
  int nx,ny;
  fftw_plan fplan,bplan;
public:
  image_c64t(complex_t**,int,int,fftw_plan,fftw_plan,io_class&);
  virtual ~image_c64t(void);
  virtual image_t *clone(io_class&);
  fp_t **convolve(fp_t**);
};

#endif                   // __IMG_C64T_H__
