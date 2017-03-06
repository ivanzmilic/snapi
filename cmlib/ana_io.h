#ifndef __ANA_READ_H__      // __ANA_READ_H__
#define __ANA_READ_H__      // __ANA_READ_H__

#define INT8	        0
#define INT16	        1
#define INT32	        2
#define FLOAT32	        3
#define FLOAT64	        4
#define INT64	        5

#define ANA_VAR_SZ {1,2,4,4,8,8}

#define MFBD_LITTLE_ENDIAN    0
#define MFBD_BIG_ENDIAN       1

#include "types.h"
#include "io.h"

int ana_data_type(const char*,io_class&);
byte *ana_fzread(char *file_name,int *&ds,int &nd,char *&header,int &type,io_class &io);
void ana_fcwrite(byte *data,char *file_name,int *ds,int nd,char *header,int type,io_class &io);
void ana_fzwrite(byte *data,char *file_name,int *ds,int nd,char *header,int type,io_class &io);

struct fzhead{                    // first block for fz files
  int synch_pattern;
  byte subf;
  byte source;
  byte nhb,datyp,ndim,file_class;
  byte cbytes[4];	      // can't do as int because of %4 rule
  byte free[178];
  int dim[16];
  char txt[256];
};

struct compresshead{
  int tsize,nblocks,bsize;
  byte slice_size,type;
};

struct ana_file{
  io_class &io;
  FILE *f;
  int type,tsize,t_endian;
  int fsize,felem;
  int size,elem;
//
  ana_file(const char *file_name,int *ds,int nd,const char *header,int type_in,io_class &io);
  ~ana_file(void);
  int append(void *data,int n,int type);
};

#endif                      // __ANA_READ_H__
