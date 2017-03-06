#ifndef __MEM_H__  // __MEM_H__
#define __MEM_H__

float *f1dim(int,int);
void del_f1dim(float*,int,int);
float **f2dim(int,int,int,int);
void del_f2dim(float**,int,int,int,int);
float ***f3dim(int,int,int,int,int,int);
void del_f3dim(float***,int,int,int,int,int,int);

double *d1dim(int,int);
void del_d1dim(double*,int,int);
double **d2dim(int,int,int,int);
void del_d2dim(double**,int,int,int,int);
double ***d3dim(int,int,int,int,int,int);
void del_d3dim(double***,int,int,int,int,int,int);
double ****d4dim(int,int,int,int,int,int,int,int);
void del_d4dim(double****,int,int,int,int,int,int,int,int);
double *****d5dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h);
void del_d5dim(double *****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h);
double ******d6dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h);
void del_d6dim(double ******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h);

int *i1dim(int,int);
void del_i1dim(int*,int,int);
int **i2dim(int,int,int,int);
void del_i2dim(int**,int,int,int,int);
int ***i3dim(int,int,int,int,int,int);
void del_i3dim(int***,int,int,int,int,int,int);
int ****i4dim(int,int,int,int,int,int,int,int);
void del_i4dim(int****,int,int,int,int,int,int,int,int);
//
int **i2dim(int,int,int,int*);
void del_i2dim(int**,int,int,int,int*);
int ***i3dim(int,int,int,int*,int,int**);
void del_i3dim(int***,int,int,int,int*,int,int**);
//
int ***i3dim(int,int,int,int,int,int*);
void del_i3dim(int***,int,int,int,int,int,int*);

void ***v2dim(int,int,int,int);
void del_v2dim(void***,int,int,int,int);
void ****v3dim(int,int,int,int,int,int);
void del_v3dim(void****,int,int,int,int,int,int);
//
void ****v3dim(int,int,int,int*,int,int**);
void del_v3dim(void****,int,int,int,int*,int,int**);
//
void ***v2dim(int,int,int,int*);
void del_v2dim(void***,int,int,int,int*);
void ****v3dim(int,int,int,int,int,int*);
void del_v3dim(void****,int,int,int,int,int,int*);

#include "types.h"

fp_t **ft2dim(int,int,int,int);
void del_ft2dim(fp_t**,int,int,int,int);
fp_t ***ft3dim(int,int,int,int,int,int);
void del_ft3dim(fp_t***,int,int,int,int,int,int);
fp_t ****ft4dim(int,int,int,int,int,int,int,int);
void del_ft4dim(fp_t****,int,int,int,int,int,int,int,int);
fp_t *****ft5dim(int,int,int,int,int,int,int,int, int, int);
void del_ft5dim(fp_t*****,int,int,int,int,int,int,int,int, int, int);
fp_t ******ft6dim(int,int,int,int,int,int,int,int,int,int,int,int);
void del_ft6dim(fp_t******,int,int,int,int,int,int,int,int,int,int,int,int);
fp_t *******ft7dim(int,int,int,int,int,int,int,int,int,int,int,int,int,int);
void del_ft7dim(fp_t*******,int,int,int,int,int,int,int,int,int,int,int,int,int,int);

// variable dimension arrays
fp_t **ft2dim(int,int,int,int*);
void del_ft2dim(fp_t**,int,int,int,int*);
fp_t ***ft3dim(int,int,int,int*,int,int**);
void del_ft3dim(fp_t***,int,int,int,int*,int,int**);
//
fp_t ****ft4dim(int,int,int,int,int,int*,int,int);
void del_ft4dim(fp_t****,int,int,int,int,int,int*,int,int);
fp_t ***ft3dim(int,int,int,int*,int,int);
void del_ft3dim(fp_t***,int,int,int,int*,int,int);
//
fp_t ***ft3dim(int,int,int,int*,int,int*);
void del_ft3dim(fp_t***,int,int,int,int*,int,int*);
fp_t ***ft3dim(int,int,int,int,int,int*);
void del_ft3dim(fp_t***,int,int,int,int,int,int*);

complex_t **ct2dim(int,int,int,int);
void del_ct2dim(complex_t**,int,int,int,int);
complex_t ***ct3dim(int,int,int,int,int,int);
void del_ct3dim(complex_t***,int,int,int,int,int,int);

int08_t **i08t2dim(int,int,int,int);
void del_i08t2dim(int08_t**,int,int,int,int);

int16_t **i16t2dim(int,int,int,int);
int16_t **i16t2dim(int16_t*,int,int,int,int);
void del_i16t2dim(int16_t**,int,int,int,int);
int16_t ***i16t3dim(int,int,int,int,int,int);
int16_t ***i16t3dim(int16_t*,int,int,int,int,int,int);
void del_i16t3dim(int16_t***,int,int,int,int,int,int);
int16_t **i16t2dim(int,int,int,int*);
void del_i16t2dim(int16_t**,int,int,int,int*);

int32_t **i32t2dim(int,int,int,int);
int32_t **i32t2dim(int32_t*,int,int,int,int);
void del_i32t2dim(int32_t**,int,int,int,int);
int32_t ***i32t3dim(int,int,int,int,int,int);
int32_t ***i32t3dim(int32_t*,int,int,int,int,int,int);
void del_i32t3dim(int32_t***,int,int,int,int,int,int);

uint32_t **ui32t2dim(int,int,int,int);
uint32_t **ui32t2dim(uint32_t*,int,int,int,int);
void del_ui32t2dim(uint32_t**,int,int,int,int);
uint32_t ***ui32t3dim(int,int,int,int,int,int);
uint32_t ***ui32t3dim(uint32_t*,int,int,int,int,int,int);
void del_ui32t3dim(uint32_t***,int,int,int,int,int,int);

int64_t **i64t2dim(int,int,int,int);
int64_t **i64t2dim(int64_t*,int,int,int,int);
void del_i64t2dim(int64_t**,int,int,int,int);
int64_t ***i64t3dim(int,int,int,int,int,int);
int64_t ***i64t3dim(int64_t*,int,int,int,int,int,int);
void del_i64t3dim(int64_t***,int,int,int,int,int,int);

float32_t **f32t2dim(int,int,int,int*);
void del_f32t2dim(float32_t**,int,int,int,int*);

float32_t **f32t2dim(int,int,int,int);
float32_t **f32t2dim(float32_t*,int,int,int,int);
void del_f32t2dim(float32_t**,int,int,int,int);
float32_t ***f32t3dim(int,int,int,int,int,int);
float32_t ***f32t3dim(float32_t*,int,int,int,int,int,int);
void del_f32t3dim(float32_t***,int,int,int,int,int,int);
float32_t ****f32t4dim(int,int,int,int,int,int,int,int);
void del_f32t4dim(float32_t****,int,int,int,int,int,int,int,int);

float64_t **f64t2dim(int,int,int,int);
float64_t **f64t2dim(float64_t*,int,int,int,int);
void del_f64t2dim(float64_t**,int,int,int,int);
float64_t ***f64t3dim(int,int,int,int,int,int);
float64_t ***f64t3dim(float64_t*,int,int,int,int,int,int);
void del_f64t3dim(float64_t***,int,int,int,int,int,int);

fp_t **ft2dim(fp_t*,int,int,int,int);
fp_t ***ft3dim(fp_t*,int,int,int,int,int,int);
fp_t ****ft4dim(fp_t*,int,int,int,int,int,int*,int,int);

byte **b2dim(int,int,int,int);
void del_b2dim(byte**,int,int,int,int);
byte ***b3dim(int,int,int,int,int,int);
void del_b3dim(byte***,int,int,int,int,int,int);
byte **b2dim(int,int,int,int*);
void del_b2dim(byte**,int,int,int,int*);
//

#endif             // __MEM_H__
