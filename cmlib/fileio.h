#ifndef __FILEIO_H__   // __FILEIO_H__
#define __FILEIO_H__

#include "types.h"
#include "io.h"

#define MFBD_I08T   0
#define MFBD_I16T   1
#define MFBD_I32T   2
#define MFBD_I64T   3
#define MFBD_F32T   4
#define MFBD_F64T   5

#define MFBD_TYPE_SIZES {1,2,4,8,4,8}
#define MFBD_TYPE_NAMES {"byte","int16","int32","int64","float32","float64"}

int data_type(const char*,io_class&);
fp_t ****read_file(char*,int&,int&,int&,int&,io_class&);
fp_t ***read_file(char*,int&,int&,int&,io_class&);
fp_t **read_file(char*,int&,int&,io_class&);
fp_t *read_file(char*,int&,io_class&);
byte *read_image(char*,int&,int&,int&,char*&,io_class&);
byte *read_file(char*,int&,int&,int&,int&,char*&,io_class&);
int16_t **read_file(char*,int&,int&,char*&,io_class&);
void write_file(char*,fp_t***,int,int,int,io_class&);
void write_file(char*,fp_t****,int,int,int*,int,io_class&);
void write_file(char*,float32_t****,int,int,int,int,io_class&);
void write_file(char*,fp_t****,int,int,int,int,io_class&);

#endif                 // __FILEIO_H__
