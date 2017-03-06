#ifndef __NR_FFT_H__  // __NR_FFT_H__
#define __NR_FFT_H__

#include "types.h"

void fft(complex_t*,int,int);       // 1D
void fft_init(int,int);
void fft_2d(complex_t**&,int,int,int);  // 2D
void fft_done(int,int);
void fft_reorder(complex_t **f,int np);
void fft_reorder_2d(complex_t **f,int nx,int ny);

#endif                // __NR_FFT_H__
