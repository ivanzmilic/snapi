#ifndef __IMG_MAT_H__  // __IMG_MAT_H__
#define __IMG_MAT_H__

#include "types.h"
#include "io.h"

class image_mat;

#include "img.h"

class image_mat:public image_t{
  int nz;
  float32_t ***pic;
public:
  image_mat(const char*,const char*,io_class&);
  virtual ~image_mat(void);
  image_t *modgain(image_t*,fp_t*,int,int);
  virtual fp_t **mult(fp_t**);
};

#endif                   // __IMG_C64T_H__
