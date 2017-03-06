#ifndef __IMG__H__  // __IMG_H__
#define __IMG__H__

#include "types.h"
#include "io.h"
#include "jnfo.h"

class image_t;

image_t *readimage(const char *dir,const char *fntmpl,int num,io_class &io);
char *imagepath(const char *dir,const char *fntmpl,int num,io_class &io);

class image_t{
protected:
  int nx,ny;
  char *header,*path;
  void descatter(fp_t**,image_t*,image_t*,io_class&);
  void fft_reorder(fp_t**,int,int);
public:
  image_t(void);
  virtual ~image_t(void);
  virtual fp_t **subimage(int,int,int,int){ return 0; };
  virtual int status(void){ return 1; }
  virtual void paste(fp_t**,int,int,int*,int*,int*,int*,int,int,fp_t,io_class &io){};
  virtual void clip(io_class*){};
  virtual void write(char *,struct jnfo&,int,io_class*){};
  virtual fp_t mean(int){ return -1.0; };
  virtual void mult(fp_t){};
  virtual fp_t **mult(fp_t**){ return 0; };
  virtual fp_t **sub(fp_t**){ return 0; };
  virtual fp_t **add(fp_t**){ return 0; };
  virtual fp_t **get(void){ return 0; };
  virtual void put(fp_t**){ };
  virtual void flatfield(image_t*,image_t*,image_t*,image_t*,image_t*,byte,io_class&);
  virtual void clip(int16_t&,int16_t&,int16_t&,int16_t&,io_class*){};
  virtual fp_t noise_sigma(int){ return -1; };
  virtual image_t *clone(io_class&){ return 0; };
  void info(int &nx_out,int &ny_out,char *&header_out)
  {
    nx_out=nx;
    ny_out=ny;
    header_out=header;
  }
  image_t *FT(io_class&);
};

void image_preproc(int&,int&,image_t****,fp_t***,fp_t**,int,int*,int**,char *&,int,io_class&);
#endif                   // __IMG_H__
