#ifndef __UTS_H__  // __UTS_H__
#define __UTS_H__

#include <string.h>
#include <math.h>
#include "types.h"

void delay(fp_t t);
int contains(int *in,int nn,int num);
int where(int *in,int nn,int num);

static __inline int pack_int(char *p,int &i)
{
  memcpy(p,&i,sizeof(int));
  return sizeof(int);
}

static __inline int unpack_int(char *p,int &i)
{
  memcpy(&i,p,sizeof(int));
  return sizeof(int);
}

static __inline int pack_chr(char *p,char *q,int n)
{
  memcpy(p,q,n);
  return n;
}

static __inline int unpack_chr(char *p,char *q,int n)
{
  memcpy(q,p,n);
  return n;
}

static __inline fp_t sqr(fp_t x)
{
  return x*x;
}

static __inline int sqr(int x)
{
  return x*x;
}

static __inline void swap(int16_t &a,int16_t &b)
{
  int16_t c=a;
  a=b;
  b=c;
}

static __inline void swap(uint16_t &a,uint16_t &b)
{
  uint16_t c=a;
  a=b;
  b=c;
}

static __inline void swap(int32_t &a,int32_t &b)
{
  int32_t c=a;
  a=b;
  b=c;
}

static __inline void swap(uint32_t &a,uint32_t &b)
{
  uint32_t c=a;
  a=b;
  b=c;
}

static __inline void swap(int64_t &a,int64_t &b)
{
  int64_t c=a;
  a=b;
  b=c;
}

static __inline void swap(float32_t &a,float32_t &b)
{
  float32_t c=a;
  a=b;
  b=c;
}

static __inline void swap(float64_t &a,float64_t &b)
{
  float64_t c=a;
  a=b;
  b=c;
}

static __inline int sign(int f)
{
  return (f>=0)?1:-1;
}

static __inline fp_t sign(fp_t f)
{
  return (f>=0.0)?1.0:-1.0;
}

static __inline fp_t sign(fp_t a,fp_t b)
{
  return (b>=0.0)?fabs(a):-fabs(a);
}

static __inline int i_round(fp_t x)
{
  return (int)((fabs(x)+0.5)*sign(x));
}

static __inline fp_t rest(fp_t x)
{
  return x-i_round(x);
}

static __inline fp_t max(fp_t a,fp_t b)
{
  return (a>b)?a:b;
}

static __inline fp_t min(fp_t a,fp_t b)
{
  return (a>b)?b:a;
}

static __inline int max(int a,int b)
{
  return (a>b)?a:b;
}

static __inline int min(int a,int b)
{
  return (a>b)?b:a;
}

static __inline fp_t dsh(fp_t prv,fp_t nu_0)
{
#define C 3E10
  return nu_0*(1.0+prv/C); // Approximate but close enough!
#undef C
}

void swap(char &a,char &b);
void swap(int64_t *x,int n);
void swap(int32_t *x,int n);
void swap(int16_t *x,int n);
void swap(void *x,int sz,int n);

void sort(int *x,int nx);
void sort(fp_t *x,int nx);
void sort(int *x,int *y,int nx);
void sort(fp_t *x,int *y,int nx);

int get_numbers(const char *numbers,fp_t *&nout,int &ni,fp_t dn);
int get_numbers(const char *numbers,int *&img_num,int &ni);
int get_numbers(const char *numbers,fp_t *&data,int &nd);
int get_number(const char*,int&);
int get_number(const char *number,fp_t &fpnum);

void in_place_transpose(float32_t *data,int L,int M);
void in_place_transpose(int32_t *data,int L,int M);
void in_place_transpose(int16_t *data,int L,int M);

#endif             // __UTS_H__
