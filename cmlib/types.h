#ifndef __TYPES_H__  // __TYPES_H__
#define __TYPES_H__

#include "../config.h"

typedef unsigned char byte;
typedef unsigned char uint08_t;
typedef signed char int08_t;
typedef short int int16_t;
typedef unsigned short int uint16_t;
typedef int int32_t;
typedef unsigned int uint32_t;
#ifdef LONG_INT64_T
  typedef long int int64_t;
  typedef unsigned long int uint64_t;
  #define I64FMT "ld"
#else
  typedef long long int int64_t;
  typedef unsigned long long int uint64_t;
  #define I64FMT "lld"
#endif
typedef float float32_t;
typedef double float64_t;

typedef float64_t fp_t;  // the default floating point type

struct complex{
  fp_t re;
  fp_t im;
//
  complex(void){};
//  complex(void):re(0.0),im(0.0){};
  complex(fp_t a,fp_t b):re(a),im(b){};
//
  complex operator+(complex a)
  {
    complex x(re+a.re,im+a.im);
    return x;
  }
  complex operator-(complex a)
  {
    complex x(re-a.re,im-a.im);
    return x;
  }
  complex operator*(complex a)
  {
    complex x(re*a.re-im*a.im,re*a.im+im*a.re);
    return x;
  }
  complex operator/(complex a)
  {
    fp_t f=a.re*a.re+a.im*a.im;
    complex x((re*a.re+im*a.im)/f,(im*a.re-re*a.im)/f);
    return x;
  }
//
  complex operator+(fp_t a)
  {
    complex x(re+a,im);
    return x;
  }
  complex operator-(fp_t a)
  {
    complex x(re-a,im);
    return x;
  }
  complex operator*(fp_t a)
  {
    complex x(re*a,im*a);
    return x;
  }
  complex operator/(fp_t a)
  {
    complex x(re/a,im/a);
    return x;
  }
//
};

static __inline complex operator+(fp_t a,complex b)
{
  complex x(b.re+a,b.im);
  return x;
};

static __inline complex operator-(fp_t a,complex b)
{
  complex x(a-b.re,-b.im);
  return x;
};

static __inline complex operator*(fp_t a,complex b)
{
  complex x(b.re*a,b.im*a);
  return x;
};

static __inline complex operator/(fp_t a,complex b)
{
  complex x(b.re/a,b.im/a);
  return x;
};

typedef struct complex complex_t;

#endif               // __TYPES_H__
