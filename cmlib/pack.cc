#include <string.h>
#include <sys/types.h>
#include <time.h>

#include "uts.h"
#include "pack.h"

int pack(byte *dest,int08_t &x)
{
  memcpy(dest,&x,sizeof(x));
  return sizeof(x);
}

int unpack(byte *src,int08_t &x)
{
  memcpy(&x,src,sizeof(x));
  return sizeof(x);
}

int pack(byte *dest,uint08_t &x)
{
  memcpy(dest,&x,sizeof(x));
  return sizeof(x);
}

int unpack(byte *src,uint08_t &x)
{
  memcpy(&x,src,sizeof(x));
  return sizeof(x);
}

int pack(byte *dest,int16_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,int16_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,uint16_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,uint16_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,int32_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,int32_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,uint32_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,uint32_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,int64_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,int64_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,uint64_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,uint64_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

/*
int unpack(byte *src,uid_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,uid_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}
*/

int pack(byte *dest,float64_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,float64_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,float32_t &x,byte swap_endian)
{
  memcpy(dest,&x,sizeof(x));
  if(swap_endian) swap(dest,sizeof(x),1);
  return sizeof(x);
}

int unpack(byte *src,float32_t &x,byte swap_endian)
{
  memcpy(&x,src,sizeof(x));
  if(swap_endian) swap(&x,sizeof(x),1);
  return sizeof(x);
}

int pack(byte *dest,int08_t *x,int xl,int xh)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(int08_t));
  return sz*sizeof(int08_t);
}

int unpack(byte *src,int08_t *x,int xl,int xh)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(int08_t));
  return sz*sizeof(int08_t);
}

int pack(byte *dest,uint08_t *x,int xl,int xh)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(uint08_t));
  return sz*sizeof(uint08_t);
}

int unpack(byte *src,uint08_t *x,int xl,int xh)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(uint08_t));
  return sz*sizeof(uint08_t);
}

int pack(byte *dest,int16_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(int16_t));
  if(swap_endian) swap(dest,sizeof(int16_t),sz);
  return sz*sizeof(int16_t);
}

int unpack(byte *src,int16_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(int16_t));
  if(swap_endian) swap(x+xl,sizeof(int16_t),sz);
  return sz*sizeof(int16_t);
}

int pack(byte *dest,uint16_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(uint16_t));
  if(swap_endian) swap(dest,sizeof(uint16_t),sz);
  return sz*sizeof(uint16_t);
}

int unpack(byte *src,uint16_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(uint16_t));
  if(swap_endian) swap(x+xl,sizeof(uint16_t),sz);
  return sz*sizeof(uint16_t);
}

int pack(byte *dest,int32_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(int32_t));
  if(swap_endian) swap(dest,sizeof(int32_t),sz);
  return sz*sizeof(int32_t);
}

int unpack(byte *src,int32_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(int32_t));
  if(swap_endian) swap(x+xl,sizeof(int32_t),sz);
  return sz*sizeof(int32_t);
}

int pack(byte *dest,uint32_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(uint32_t));
  if(swap_endian) swap(dest,sizeof(uint32_t),sz);
  return sz*sizeof(uint32_t);
}

int unpack(byte *src,uint32_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(uint32_t));
  if(swap_endian) swap(x+xl,sizeof(uint32_t),sz);
  return sz*sizeof(uint32_t);
}

int pack(byte *dest,float64_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(float64_t));
  if(swap_endian) swap(dest,sizeof(float64_t),sz);
  return sz*sizeof(float64_t);
}

int unpack(byte *src,float64_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(float64_t));
  if(swap_endian) swap(x+xl,sizeof(float64_t),sz);
  return sz*sizeof(float64_t);
}

int pack(byte *dest,float32_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(dest,x+xl,sz*sizeof(float32_t));
  if(swap_endian) swap(dest,sizeof(float32_t),sz);
  return sz*sizeof(float32_t);
}

int unpack(byte *src,float32_t *x,int xl,int xh,byte swap_endian)
{
  int sz=(xh-xl+1);
  memcpy(x+xl,src,sz*sizeof(float32_t));
  if(swap_endian) swap(x+xl,sizeof(float32_t),sz);
  return sz*sizeof(float32_t);
}

int pack(byte *dest,int16_t **x,int xl,int xh,int yl,int yh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1);
  memcpy(dest,x[xl]+yl,sz*sizeof(int16_t));
  if(swap_endian) swap(dest,sizeof(int16_t),sz);
  return sz*sizeof(int16_t);
}

int unpack(byte *src,int16_t **x,int xl,int xh,int yl,int yh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1);
  memcpy(x[xl]+yl,src,sz*sizeof(int16_t));
  if(swap_endian) swap(x[xl]+yl,sizeof(int16_t),sz);
  return sz*sizeof(int16_t);
}

int pack(byte *dest,int16_t **x,int xl,int xh,int yl,int *yh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(dest,x[xl]+yl,sz*sizeof(int16_t));
  if(swap_endian) swap(dest,sizeof(int16_t),sz);
  return sz*sizeof(int16_t);
}

int unpack(byte *src,int16_t **x,int xl,int xh,int yl,int *yh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(x[xl]+yl,src,sz*sizeof(int16_t));
  if(swap_endian) swap(x[xl]+yl,sizeof(int16_t),sz);
  return sz*sizeof(int16_t);
}

int pack(byte *dest,int32_t **x,int xl,int xh,int yl,int *yh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(dest,x[xl]+yl,sz*sizeof(int32_t));
  if(swap_endian) swap(dest,sizeof(int32_t),sz);
  return sz*sizeof(int32_t);
}

int unpack(byte *src,int32_t **x,int xl,int xh,int yl,int *yh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(x[xl]+yl,src,sz*sizeof(int32_t));
  if(swap_endian) swap(x[xl]+yl,sizeof(int32_t),sz);
  return sz*sizeof(int32_t);
}

int pack(byte *dest,fp_t **x,int xl,int xh,int yl,int yh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1);
  memcpy(dest,x[xl]+yl,sz*sizeof(fp_t));
  if(swap_endian) swap(dest,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int unpack(byte *src,fp_t **x,int xl,int xh,int yl,int yh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1);
  memcpy(x[xl]+yl,src,sz*sizeof(fp_t));
  if(swap_endian) swap(x[xl]+yl,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int pack(byte *dest,fp_t **x,int xl,int xh,int yl,int *yh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(dest,x[xl]+yl,sz*sizeof(fp_t));
  if(swap_endian) swap(dest,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int unpack(byte *src,fp_t **x,int xl,int xh,int yl,int *yh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(x[xl]+yl,src,sz*sizeof(fp_t));
  if(swap_endian) swap(x[xl]+yl,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int pack(byte *dest,byte **x,int xl,int xh,int yl,int *yh)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(dest,x[xl]+yl,sz*sizeof(byte));
  return sz*sizeof(byte);
}

int unpack(byte *src,byte **x,int xl,int xh,int yl,int *yh)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  memcpy(x[xl]+yl,src,sz*sizeof(byte));
  return sz*sizeof(byte);
}

int pack(byte *dest,int32_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix)
    for(int iy=yl;iy<=yh[ix];++iy) sz+=(zh[ix][iy]-zl+1);
  memcpy(dest,x[xl][yl]+zl,sz*sizeof(int32_t));
  if(swap_endian) swap(dest,sizeof(int32_t),sz);
  return sz*sizeof(int32_t);
}

int unpack(byte *src,int32_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte swap_endian)
{
  int nx=(xh-xl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix)
    for(int iy=yl;iy<=yh[ix];++iy) sz+=(zh[ix][iy]-zl+1);
  memcpy(x[xl][yl]+zl,src,sz*sizeof(int32_t));
  if(swap_endian) swap(x[xl][yl]+zl,sizeof(int32_t),sz);
  return sz*sizeof(int32_t);
}

int pack(byte *dest,fp_t ***x,int xl,int xh,int yl,int yh,int zl,int zh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1)*(zh-zl+1);
  memcpy(dest,x[xl][yl]+zl,sz*sizeof(fp_t));
  if(swap_endian) swap(dest,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int unpack(byte *src,fp_t ***x,int xl,int xh,int yl,int yh,int zl,int zh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1)*(zh-zl+1);
  memcpy(x[xl][yl]+zl,src,sz*sizeof(fp_t));
  if(swap_endian) swap(x[xl][yl]+zl,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int pack(byte *dest,fp_t ****x,int xl,int xh,int yl,int yh,int zl,int zh,int wl,int wh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1)*(zh-zl+1)*(wh-wl+1);
  memcpy(dest,x[xl][yl][zl]+wl,sz*sizeof(fp_t));
  if(swap_endian) swap(dest,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int unpack(byte *src,fp_t ****x,int xl,int xh,int yl,int yh,int zl,int zh,int wl,int wh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1)*(zh-zl+1)*(wh-wl+1);
  memcpy(x[xl][yl][zl]+wl,src,sz*sizeof(fp_t));
  if(swap_endian) swap(x[xl][yl][zl]+wl,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}


int pack(byte *dest,fp_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte swap_endian)
{
  int sz=0;
  for(int ix=xl;ix<=xh;++ix)
    for(int iy=yl;iy<=yh[ix];++iy) sz+=(zh[ix][iy]-zl+1);
  memcpy(dest,x[xl][yl]+zl,sz*sizeof(fp_t));
  if(swap_endian) swap(dest,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int unpack(byte *src,fp_t ***x,int xl,int xh,int yl,int *yh,int zl,int **zh,byte swap_endian)
{
  int sz=0;
  for(int ix=xl;ix<=xh;++ix)
    for(int iy=yl;iy<=yh[ix];++iy) sz+=(zh[ix][iy]-zl+1);
  memcpy(x[xl][yl]+zl,src,sz*sizeof(fp_t));
  if(swap_endian) swap(x[xl][yl]+zl,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int pack(byte *dest,fp_t ****x,int wl,int wh,int xl,int xh,int yl,int *yh,int zl,int zh,byte swap_endian)
{
  int nw=(wh-wl+1),nx=(xh-xl+1),nz=(zh-zl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  sz*=nw*nz;
  memcpy(dest,x[wl][xl][yl]+zl,sz*sizeof(fp_t));
  if(swap_endian) swap(dest,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int unpack(byte *src,fp_t ****x,int wl,int wh,int xl,int xh,int yl,int *yh,int zl,int zh,byte swap_endian)
{
  int nw=(wh-wl+1),nx=(xh-xl+1),nz=(zh-zl+1),sz=0;
  for(int ix=xl;ix<=xh;++ix) sz+=(yh[ix]-yl+1);
  sz*=nw*nz;
  memcpy(x[wl][xl][yl]+zl,src,sz*sizeof(fp_t));
  if(swap_endian) swap(x[wl][xl][yl]+zl,sizeof(fp_t),sz);
  return sz*sizeof(fp_t);
}

int pack(byte *dest,char *s)
{
  int sz=strlen(s)+1;
  memcpy(dest,s,sz);
  return sz;
}

int unpack(byte *src,char *&s)
{
  int sz=strlen((char*)src)+1;
  memcpy(s=new char [sz],src,sz);
  return sz;
}

int pack(byte *dest,float32_t **x,int xl,int xh,int yl,int yh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1);
  memcpy(dest,x[xl]+yl,sz*sizeof(float32_t));
  if(swap_endian) swap(dest,sizeof(float32_t),sz);
  return sz*sizeof(float32_t);
}

int unpack(byte *src,float32_t **x,int xl,int xh,int yl,int yh,byte swap_endian)
{
  int sz=(xh-xl+1)*(yh-yl+1);
  memcpy(x[xl]+yl,src,sz*sizeof(float32_t));
  if(swap_endian) swap(x[xl]+yl,sizeof(float32_t),sz);
  return sz*sizeof(float32_t);
}

