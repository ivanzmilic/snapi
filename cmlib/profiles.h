#ifndef __PROFILES_H__      //  __PROFILES_H__
#define __PROFILES_H__

#include "types.h"

void fvoigt(fp_t x,fp_t a,fp_t &hh,fp_t &ff,fp_t &dh,fp_t &df);
void fvoigtn(fp_t x,fp_t a,fp_t &hh,fp_t &ff,fp_t &dh,fp_t &df);
fp_t fvoigt(fp_t v,fp_t a);
void voigt_num_der(fp_t x, fp_t a,fp_t &hh,fp_t &ff,fp_t &dh,fp_t &df);

#endif                      //  __PROFILES_H__
 
