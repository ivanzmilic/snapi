#ifndef __IMG_I16T_H__  // __IMG_I16T_H__
#define __IMG_I16T_H__

#include "types.h"
#include "io.h"

class image_i16t;

#include "img.h"
//#include "img_f32t.h"
#include "jnfo.h"

class image_i16t:public image_t{
  int16_t **pic;
public:
  image_i16t(int,int);
  image_i16t(image_t*,io_class &io);
  image_i16t(const char *path,io_class &io);
//  image_i16t(const char *dir,const char *fntmpl,io_class &io);
  image_i16t(const char *dir,const char *fntmpl,int num,io_class &io);
  virtual ~image_i16t(void);
  virtual int status(void){ return !pic; }
  virtual fp_t **subimage(int,int,int,int);
  virtual void paste(fp_t**,int,int,int*,int*,int*,int*,int,int,fp_t,io_class &io);
  virtual void clip(io_class*);
  virtual void write(char *,struct jnfo&,int,io_class*);
  void clip(int,int,int,int);
  virtual void clip(int16_t&,int16_t&,int16_t&,int16_t&,io_class*);
  fp_t noise_sigma(int);
  virtual fp_t mean(int);
  virtual void mult(fp_t);
  virtual fp_t **mult(fp_t**);
  virtual fp_t **sub(fp_t**);
  virtual fp_t **add(fp_t**);
  virtual fp_t **get(void);
  virtual void put(fp_t**);
//  virtual void flatfield(image_t*,image_t*,image_t*,image_t*,byte,io_class&);
  virtual image_t *clone(io_class&);
};

#endif                   // __IMG_I16T_H__
