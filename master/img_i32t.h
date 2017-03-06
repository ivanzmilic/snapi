#ifndef __IMG_I32T_H__  // __IMG_I32T_H__
#define __IMG_I32T_H__

#include "types.h"
#include "io.h"
#include "jnfo.h"

class image_i32t;

#include "img.h"
#include "img_i32t.h"

class image_i32t:public image_t{
  int32_t **pic;
public:
  image_i32t(int,int);
  image_i32t(const char*,io_class&);
  image_i32t(image_t*,io_class&);
//  image_i32t(const char*,const char*,io_class&);
  image_i32t(const char *dir,const char *fntmpl,int num,io_class &io);
  virtual ~image_i32t(void);
  virtual fp_t **subimage(int,int,int,int);
  virtual void paste(fp_t**,int,int,int*,int*,int*,int*,int,int,fp_t,io_class &io);
  virtual void write(char *,struct jnfo&,int,io_class*);
  virtual void clip(io_class*);
  virtual int status(void){ return !pic; }
  virtual fp_t mean(int);
  virtual void mult(fp_t);
  virtual fp_t **safediv(fp_t**);
  virtual fp_t **mult(fp_t**);
  virtual fp_t **sub(fp_t**);
  virtual fp_t **add(fp_t**);
  virtual fp_t **get(void);
  virtual void put(fp_t**);
//  virtual void flatfield(image_t*,image_t*,image_t*,image_t*,byte,io_class&);
  virtual void clip(int16_t&,int16_t&,int16_t&,int16_t&,io_class*);
  virtual fp_t noise_sigma(int);
  virtual image_t *clone(io_class&);
};

#endif                   // __IMG_I32T_H__
