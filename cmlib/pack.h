#ifndef __PACK_H__    // __PACK_H__
#define __PACK_H__

#include <sys/types.h>
#include <time.h>
#include "types.h"

int pack(byte *dest,int08_t &x);
int unpack(byte *src,int08_t &x);
int pack(byte *dest,uint08_t &x);
int unpack(byte *src,uint08_t &x);
int pack(byte *dest,int16_t &x,byte);
int unpack(byte *src,int16_t &x,byte);
int pack(byte *dest,uint16_t &x,byte);
int unpack(byte *src,uint16_t &x,byte);
int pack(byte *dest,int32_t &x,byte);
int unpack(byte *src,int32_t &x,byte);
int pack(byte *dest,uint32_t &x,byte);
int unpack(byte *src,uint32_t &x,byte);
int pack(byte *dest,int64_t &x,byte);
int unpack(byte *src,int64_t &x,byte);
int pack(byte *dest,uint64_t &x,byte);
int unpack(byte *src,uint64_t &x,byte);
//int pack(byte *dest,uid_t &x,byte);
//int unpack(byte *src,uid_t &x,byte);
int pack(byte *dest,float64_t &x,byte);
int unpack(byte *src,float64_t &x,byte);
int pack(byte *dest,float32_t &x,byte);
int unpack(byte *src,float32_t &x,byte);
int pack(byte *dest,int08_t *x,int xl,int xh);
int unpack(byte *src,int08_t *x,int xl,int xh);
int pack(byte *dest,uint08_t *x,int xl,int xh);
int unpack(byte *src,uint08_t *x,int xl,int xh);
int pack(byte *dest,int16_t *x,int xl,int xh,byte);
int unpack(byte *src,int16_t *x,int xl,int xh,byte);
int pack(byte *dest,uint16_t *x,int xl,int xh,byte);
int unpack(byte *src,uint16_t *x,int xl,int xh,byte);
int pack(byte *dest,int32_t *x,int xl,int xh,byte);
int unpack(byte *src,int32_t *x,int xl,int xh,byte);
int pack(byte *dest,uint32_t *x,int xl,int xh,byte);
int unpack(byte *src,uint32_t *x,int xl,int xh,byte);
int pack(byte *dest,float64_t *x,int xl,int xh,byte);
int unpack(byte *src,float64_t *x,int xl,int xh,byte);
int pack(byte *dest,float32_t *x,int xl,int xh,byte);
int unpack(byte *src,float32_t *x,int xl,int xh,byte);
int pack(byte *dest,int16_t **x,int xl,int xh,int yl,int yh,byte);
int unpack(byte *dest,int16_t **x,int xl,int xh,int yl,int yh,byte);
int pack(byte *dest,int16_t **x,int xl,int xh,int yl,int *yh,byte);
int unpack(byte *dest,int16_t **x,int xl,int xh,int yl,int *yh,byte);
int pack(byte *dest,int32_t **x,int xl,int xh,int yl,int *yh,byte);
int unpack(byte *dest,int32_t **x,int xl,int xh,int yl,int *yh,byte);
int pack(byte *dest,fp_t **x,int xl,int xh,int yl,int yh,byte);
int unpack(byte *dest,fp_t **x,int xl,int xh,int yl,int yh,byte);
int pack(byte *dest,fp_t **x,int xl,int xh,int yl,int *yh,byte);
int unpack(byte *dest,fp_t **x,int xl,int xh,int yl,int *yh,byte);
int pack(byte *dest,byte **x,int xl,int xh,int yl,int *yh);
int unpack(byte *dest,byte **x,int xl,int xh,int yl,int *yh);
int pack(byte *dest,int32_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte);
int unpack(byte *dest,int32_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte);
int pack(byte *dest,fp_t ***x,int xl,int xh,int yl,int yh,int zl,int zh,byte);
int unpack(byte *dest,fp_t ***x,int xl,int xh,int yl,int yh,int zl,int zh,byte);
int pack(byte *dest,fp_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte);
int unpack(byte *dest,fp_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte);
int pack(byte *dest,fp_t ****x,int wl,int wh,int xl,int xh,int yl,int *yh,int zl,int zh,byte);
int unpack(byte *dest,fp_t ****x,int wl,int wh,int xl,int xh,int yl,int *yh,int zl,int zh,byte);
int pack(byte *dest,char *s);
int unpack(byte *src,char *&s);
int pack(byte *dest,float32_t **x,int xl,int xh,int yl,int yh,byte);
int unpack(byte *dest,float32_t **x,int xl,int xh,int yl,int yh,byte);

#endif                // __PACK_H__
