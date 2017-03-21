#include "mem.h"
#include <stdio.h>

float *f1dim(int x1l,int x1h)
{
  return new float [x1h-x1l+1] - x1l;
}

void del_f1dim(float *p,int x1l,int x1h)
{
  delete[] (p+x1l);
}

float **f2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  float **p;
  p=new float* [nx1] - x1l;
  p[x1l]=new float [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_f2dim(float **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

float ***f3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  float ***p;
  p=new float** [nx1] - x1l;
  p[x1l]=new float* [nx1*nx2] - x2l;
  p[x1l][x2l]=new float [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_f3dim(float ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double *d1dim(int x1l,int x1h)
{
  return new double [x1h-x1l+1] - x1l;
}

void del_d1dim(double *p,int x1l,int x1h)
{
  delete[] (p+x1l);
}

double **d2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  double **p;
  p=new double* [nx1] - x1l;
  p[x1l]=new double [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_d2dim(double **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double ***d3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  double ***p;
  p=new double** [nx1] - x1l;
  p[x1l]=new double* [nx1*nx2] - x2l;
  p[x1l][x2l]=new double [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_d3dim(double ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double ****d4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1;
  double ****p;
  p=new double*** [nx1] - x1l;
  p[x1l]=new double** [nx1*nx2] - x2l;
  p[x1l][x2l]=new double* [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new double [nx1*nx2*nx3*nx4] - x4l;
  for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  return p;
}

void del_d4dim(double ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double *****d5dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1;
  double *****p;
  p=new double**** [nx1] - x1l;
  p[x1l]=new double*** [nx1*nx2] - x2l;
  p[x1l][x2l]=new double** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new double* [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new double [nx1*nx2*nx3*nx4*nx5] - x5l;
  for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;    
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5;
  }
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;
    }
  }	
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;
      for(int x3=x3l+1;x3<=x3h;++x3){
        p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
        p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
        for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;
      }
    }
  }
  return p;
}

void del_d5dim(double *****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h)
{
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double ******d6dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1,nx6=x6h-x6l+1;
  double ******p;
  p=new double***** [nx1] - x1l;
  p[x1l]=new double**** [nx1*nx2] - x2l;
  p[x1l][x2l]=new double*** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new double** [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new double* [nx1*nx2*nx3*nx4*nx5] - x5l;
  p[x1l][x2l][x3l][x4l][x5l]=new double [nx1*nx2*nx3*nx4*nx5*nx6] - x6l;
  for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4l][x5]=p[x1l][x2l][x3l][x4l][x5-1]+nx6;
  for(int x4=x4l+1;x4<=x4h;++x4){
    p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;    
    p[x1l][x2l][x3l][x4][x5l]=p[x1l][x2l][x3l][x4-1][x5l]+nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4][x5]=p[x1l][x2l][x3l][x4][x5-1]+nx6;
  }
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    p[x1l][x2l][x3][x4l][x5l]=p[x1l][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4l][x5]=p[x1l][x2l][x3][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5;    
      p[x1l][x2l][x3][x4][x5l]=p[x1l][x2l][x3][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4][x5]=p[x1l][x2l][x3][x4][x5-1]+nx6;
    }
  }
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    p[x1l][x2][x3l][x4l][x5l]=p[x1l][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4l][x5]=p[x1l][x2][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;    
      p[x1l][x2][x3l][x4][x5l]=p[x1l][x2][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4][x5]=p[x1l][x2][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      p[x1l][x2][x3][x4l][x5l]=p[x1l][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4l][x5]=p[x1l][x2][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;    
        p[x1l][x2][x3][x4][x5l]=p[x1l][x2][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4][x5]=p[x1l][x2][x3][x4][x5-1]+nx6;
      }
    }
  }
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    p[x1][x2l][x3l][x4l][x5l]=p[x1-1][x2l][x3l][x4l][x5l]+nx2*nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4l][x5]=p[x1][x2l][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;    
      p[x1][x2l][x3l][x4][x5l]=p[x1][x2l][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4][x5]=p[x1][x2l][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      p[x1][x2l][x3][x4l][x5l]=p[x1][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4l][x5]=p[x1][x2l][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;    
        p[x1][x2l][x3][x4][x5l]=p[x1][x2l][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4][x5]=p[x1][x2l][x3][x4][x5-1]+nx6;
      }
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      p[x1][x2][x3l][x4l][x5l]=p[x1][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4l][x5]=p[x1][x2][x3l][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;    
        p[x1][x2][x3l][x4][x5l]=p[x1][x2][x3l][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4][x5]=p[x1][x2][x3l][x4][x5-1]+nx6;
      }
      for(int x3=x3l+1;x3<=x3h;++x3){
        p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
        p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
        p[x1][x2][x3][x4l][x5l]=p[x1][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4l][x5]=p[x1][x2][x3][x4l][x5-1]+nx6;
        for(int x4=x4l+1;x4<=x4h;++x4){
          p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;    
          p[x1][x2][x3][x4][x5l]=p[x1][x2][x3][x4-1][x5l]+nx5*nx6;
          for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4][x5]=p[x1][x2][x3][x4][x5-1]+nx6;
        }
      }
    }
  }
  return p;
}

void del_d6dim(double ******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  delete[] (p[x1l][x2l][x3l][x4l][x5l]+x6l);
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int *i1dim(int x1l,int x1h)
{
  return new int [x1h-x1l+1] - x1l;
}

void del_i1dim(int *p,int x1l,int x1h)
{
  delete[] (p+x1l);
}

int **i2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int **p;
  p=new int* [nx1] - x1l;
  p[x1l]=new int [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_i2dim(int **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ***i3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int ***p;
  p=new int** [nx1] - x1l;
  p[x1l]=new int* [nx1*nx2] - x2l;
  p[x1l][x2l]=new int [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_i3dim(int ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ****i4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1;
  int ****p;
  p=new int*** [nx1] - x1l;
  p[x1l]=new int** [nx1*nx2] - x2l;
  p[x1l][x2l]=new int* [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new int [nx1*nx2*nx3*nx4] - x4l;
  for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  return p;
}

void del_i4dim(int ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int **i2dim(int x1l,int x1h,int x2l,int *x2h)
{
  int nx1=x1h-x1l+1,nx12=0;
  int **p;
  p=new int* [nx1] - x1l;
  for(int x1=x1l;x1<=x1h;++x1) nx12+=(x2h[x1]-x2l+1);
  p[x1l]=new int [nx12] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+(x2h[x1-1]-x2l+1);
  return p;
}

void del_i2dim(int **p,int x1l,int x1h,int x2l,int *x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ***i3dim(int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  int nx1=x1h-x1l+1;
  int *nx2=new int [nx1]-x1l,nx12=0;
  int **nx3=i2dim(x1l,x1h,x2l,x2h),nx123=0;
  for(int x1=x1l;x1<=x1h;++x1){
    nx12+=(nx2[x1]=x2h[x1]-x2l+1);
    for(int x2=x2l;x2<=x2h[x1];++x2) nx123+=(nx3[x1][x2]=x3h[x1][x2]-x3l+1);
  }
//
  int ***p;
  p=new int** [nx1]-x1l;
  p[x1l]=new int* [nx12]-x2l;
  p[x1l][x2l]=new int [nx123]-x3l;
  for(int x2=x2l+1;x2<=x2h[x1l];++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l][x2-1];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2[x1-1];
    p[x1][x2l]=p[x1-1][x2h[x1-1]]+nx3[x1-1][x2h[x1-1]];
    for(int x2=x2l+1;x2<=x2h[x1];++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1][x2-1];
  }
//
  del_i2dim(nx3,x1l,x1h,x2l,x2h);
  delete[] (nx2+x1l);
//
  return p;
}

void del_i3dim(int ***p,int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ***i3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nxt=0;
  int ***p;
  p=new int** [nx1] - x1l;
  p[x1l]=new int* [nx1*nx2] - x2l;
  int *nx3=new int [nx1] -x1l;
  for(int x1=x1l;x1<=x1h;++x1) nxt+=(nx3[x1]=x3h[x1]-x3l+1);
  p[x1l][x2l]=new int [nx2*nxt] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3[x1-1];
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1];
  }
  delete[] (nx3+x1l);
  return p;
}

void del_i3dim(int ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ***v2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  void ***p;
  p=new void** [nx1] - x1l;
  p[x1l]=new void* [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_v2dim(void ***p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ****v3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  void ****p;
  p=new void*** [nx1] - x1l;
  p[x1l]=new void** [nx1*nx2] - x2l;
  p[x1l][x2l]=new void* [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_v3dim(void ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ****v3dim(int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  int nx1=x1h-x1l+1;
  int *nx2=new int [nx1]-x1l,nx12=0;
  int **nx3=i2dim(x1l,x1h,x2l,x2h),nx123=0;
  for(int x1=x1l;x1<=x1h;++x1){
    nx12+=(nx2[x1]=x2h[x1]-x2l+1);
    for(int x2=x2l;x2<=x2h[x1];++x2) nx123+=(nx3[x1][x2]=x3h[x1][x2]-x3l+1);
  }
//
  void ****p;
  p=new void*** [nx1]-x1l;
  p[x1l]=new void** [nx12]-x2l;
  p[x1l][x2l]=new void* [nx123]-x3l;
  for(int x2=x2l+1;x2<=x2h[x1l];++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l][x2-1];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2[x1-1];
    p[x1][x2l]=p[x1-1][x2h[x1-1]]+nx3[x1-1][x2h[x1-1]];
    for(int x2=x2l+1;x2<=x2h[x1];++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1][x2-1];
  }
//
  del_i2dim(nx3,x1l,x1h,x2l,x2h);
  delete[] (nx2+x1l);
//
  return p;
}

void del_v3dim(void ****p,int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ***v2dim(int x1l,int x1h,int x2l,int *x2h)
{
  int nx1=x1h-x1l+1,nxt=0;
  void ***p;
  p=new void** [nx1] - x1l;
  for(int x1=x1l;x1<=x1h;++x1) nxt+=x2h[x1]-x2l+1;
  p[x1l]=new void* [nxt] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+(x2h[x1-1]-x2l+1);
  return p;
}

void del_v2dim(void ***p,int x1l, int x1h, int x2l, int *x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ****v3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nxt=0;
  void ****p;
  p=new void*** [nx1] - x1l;
  p[x1l]=new void** [nx1*nx2] - x2l;
  int *nx3=new int [nx1] -x1l;
  for(int x1=x1l;x1<=x1h;++x1) nxt+=(nx3[x1]=x3h[x1]-x3l+1);
  p[x1l][x2l]=new void* [nx2*nxt] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3[x1-1];
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1];
  }
  delete[] (nx3+x1l);
  return p;
}

void del_v3dim(void ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t **ft2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  fp_t **p;
  p=new fp_t* [nx1] - x1l;
  p[x1l]=new fp_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_ft2dim(fp_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ***ft3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  fp_t ***p;
  p=new fp_t** [nx1] - x1l;
  p[x1l]=new fp_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new fp_t [nx1*nx2*nx3] - x3l;
  
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
  
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_ft3dim(fp_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t **ft2dim(int x1l,int x1h,int x2l,int *x2h)
{
  int nx1=x1h-x1l+1,nx12=0;
  fp_t **p;
  p=new fp_t* [nx1] - x1l;
  for(int x1=x1l;x1<=x1h;++x1) nx12+=(x2h[x1]-x2l+1);
  p[x1l]=new fp_t [nx12] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+(x2h[x1-1]-x2l+1);
  return p;
}

void del_ft2dim(fp_t **p,int x1l,int x1h,int x2l,int *x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ***ft3dim(int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  int nx1=x1h-x1l+1;
  int *nx2=new int [nx1]-x1l,nx12=0;
  int **nx3=i2dim(x1l,x1h,x2l,x2h),nx123=0;
  for(int x1=x1l;x1<=x1h;++x1){
    nx12+=(nx2[x1]=x2h[x1]-x2l+1);
    for(int x2=x2l;x2<=x2h[x1];++x2) nx123+=(nx3[x1][x2]=x3h[x1][x2]-x3l+1);
  }
//
  fp_t ***p;
  p=new fp_t** [nx1]-x1l;
  p[x1l]=new fp_t* [nx12]-x2l;
  p[x1l][x2l]=new fp_t [nx123]-x3l;
  for(int x2=x2l+1;x2<=x2h[x1l];++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l][x2-1];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2[x1-1];
    p[x1][x2l]=p[x1-1][x2h[x1-1]]+nx3[x1-1][x2h[x1-1]];
    for(int x2=x2l+1;x2<=x2h[x1];++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1][x2-1];
  }
//
  del_i2dim(nx3,x1l,x1h,x2l,x2h);
  delete[] (nx2+x1l);
//
  return p;
}

void del_ft3dim(fp_t ***p,int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ***ft3dim(int x1l,int x1h,int x2l,int *x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx3=x3h-x3l+1;
  int nx12=0,*nx2=new int [nx1]-x1l;
  for(int x1=x1l;x1<=x1h;++x1) nx12+=(nx2[x1]=x2h[x1]-x2l+1);
  fp_t ***p;
  p=new fp_t** [nx1] - x1l;
  p[x1l]=new fp_t* [nx12] - x2l;
  p[x1l][x2l]=new fp_t [nx12*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h[x1l];++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2[x1-1];
    p[x1][x2l]=p[x1-1][x2l]+nx2[x1-1]*nx3;
    for(int x2=x2l+1;x2<=x2h[x1];++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  delete[] (nx2+x1l);
  return p;
}

void del_ft3dim(fp_t ***p,int x1l,int x1h,int x2l,int *x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ***ft3dim(int x1l,int x1h,int x2l,int *x2h,int x3l,int *x3h)
{
  int nx1=x1h-x1l+1;
  int nx12=0,*nx2=new int [nx1]-x1l;
  int nx123=0,*nx3=new int [nx1]-x1l;
  for(int x1=x1l;x1<=x1h;++x1){
    nx12+=(nx2[x1]=x2h[x1]-x2l+1);
    nx123+=(nx3[x1]=x3h[x1]-x3l+1)*nx2[x1];
  }
  fp_t ***p;
  p=new fp_t** [nx1] - x1l;
  p[x1l]=new fp_t* [nx12] - x2l;
  p[x1l][x2l]=new fp_t [nx123] - x3l;
  for(int x2=x2l+1;x2<=x2h[x1l];++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2[x1-1];
    p[x1][x2l]=p[x1-1][x2l]+nx2[x1-1]*nx3[x1-1];
    for(int x2=x2l+1;x2<=x2h[x1];++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1];
  }
  delete[] (nx3+x1l);
  delete[] (nx2+x1l);
  return p;
}

void del_ft3dim(fp_t ***p,int x1l,int x1h,int x2l,int *x2h,int x3l,int *x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ***ft3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nxt=0;
  fp_t ***p;
  p=new fp_t** [nx1] - x1l;
  p[x1l]=new fp_t* [nx1*nx2] - x2l;
  for(int x1=x1l;x1<=x1h;++x1) nxt+=x3h[x1]-x3l+1;
  p[x1l][x2l]=new fp_t [nx2*nxt] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+(x3h[x1l]-x3l+1);
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*(x3h[x1-1]-x3l+1);
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+(x3h[x1]-x3l+1);
  }
  return p;
}

void del_ft3dim(fp_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ****ft4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1;
  fp_t ****p;
  p=new fp_t*** [nx1] - x1l;
  p[x1l]=new fp_t** [nx1*nx2] - x2l;
  p[x1l][x2l]=new fp_t* [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new fp_t [nx1*nx2*nx3*nx4] - x4l;
  for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  return p;
}

void del_ft4dim(fp_t ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ****ft4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx23t=0,nx4=x4h-x4l+1;
  fp_t ****p;
  p=new fp_t*** [nx1] - x1l;
  p[x1l]=new fp_t** [nx1*nx2] - x2l;
  int *nx3=new int [nx2]-x2l;
  for(int x2=x2l;x2<=x2h;++x2) nx23t+=(nx3[x2]=x3h[x2]-x3l+1); // equivalent to nx2*nx3
  p[x1l][x2l]=new fp_t* [nx1*nx23t] - x3l;
  p[x1l][x2l][x3l]=new fp_t [nx1*nx23t*nx4] - x4l;
  for(int x3=x3l+1;x3<=x3h[x2l];++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3[x2-1];
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3[x2-1]*nx4;
    for(int x3=x3l+1;x3<=x3h[x2];++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx23t;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx23t*nx4;
    for(int x3=x3l+1;x3<=x3h[x2l];++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3[x2-1];
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3[x2-1]*nx4;
      for(int x3=x3l+1;x3<=x3h[x2];++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  delete[] (nx3+x2l);
  return p;
}

void del_ft4dim(fp_t ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h,int x4l,int x4h)
{
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

// ------------------------------------------------------------------------------------------

// New function, for 5-dimensional array

/*fp_t ****ft4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1;
  fp_t ****p;
  p=new fp_t*** [nx1] - x1l;
  p[x1l]=new fp_t** [nx1*nx2] - x2l;
  p[x1l][x2l]=new fp_t* [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new fp_t [nx1*nx2*nx3*nx4] - x4l;

  for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  return p;
}*/

fp_t ***** ft5dim(int x1l, int x1h, int x2l, int x2h, int x3l, int x3h, int x4l, int x4h, int x5l, int x5h){
  
  int nx1 = x1h - x1l+1, nx2 = x2h - x2l+1, nx3 = x3h - x3l+1, nx4 = x4h - x4l+1, nx5 = x5h - x5l+1; 

 // printf("%d %d %d %d %d \n", nx1, nx2, nx3, nx4, nx5);

  fp_t ***** p;

  // Milic: Let's make sure that I understand this :) 

  p = new fp_t **** [nx1] - x1l; // Now we have nx1 long array of fp_t ****
  p[x1l] = new fp_t *** [nx1*nx2] - x2l; // Now we have nx1 x nx2 array of fp_t *** ...
  p[x1l][x2l] = new fp_t ** [nx1*nx2*nx3] - x3l; // fp_t **...
  p[x1l][x2l][x3l] = new fp_t * [nx1*nx2*nx3*nx4] - x4l; // fp_t *...
  p[x1l][x2l][x3l][x4l] = new fp_t [nx1 * nx2 * nx3 * nx4 * nx5] - x5l; //  fp_t ...

  // All memory is now allocated now what needs to be done tho is to point everyone where it is needed. :)
  // Let's do it the way MvN is doing that, but later convice ourselves that it MUST be done that way. 

  // We start with pointers to  *** 

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1] = p[x1-1] + nx2;

  // Then the pointer to  **

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l] = p[x1-1][x2l] + nx2 * nx3;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2] = p[x1][x2-1] + nx3;

  // We then continue further witn the pointers to *

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l][x3l] = p[x1-1][x2l][x3l] + nx2 * nx3 * nx4;

  
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2][x3l] = p[x1][x2-1][x3l] + nx3 * nx4;

  
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l+1; x3<=x3h; ++x3)
        p[x1][x2][x3] = p[x1][x2][x3-1] + nx4;

  // And finally, the pointers to double:

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l][x3l][x4l] = p[x1-1][x2l][x3l][x4l] + nx2 * nx3 * nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2][x3l][x4l] = p[x1][x2-1][x3l][x4l] + nx3 * nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l+1; x3<=x3h; ++x3)
        p[x1][x2][x3][x4l] = p[x1][x2][x3-1][x4l] + nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l; x3<=x3h; ++x3)
        for (int x4 = x4l+1; x4<=x4h; ++x4)
          p[x1][x2][x3][x4] = p[x1][x2][x3][x4-1] + nx5;

  //printf ("Hush-hush. \n What what? \n");

  return p;
}

void del_ft5dim(fp_t *****p, int x1l, int x1h, int x2l, int x2h, int x3l, int x3h, int x4l, int x4h, int x5l, int x5h){

  // Analogous to the first four:
  delete[] (p[x1l][x2l][x3l][x4l] + x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ******ft6dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1,nx6=x6h-x6l+1;
  fp_t ******p;
  p=new fp_t***** [nx1] - x1l;
  p[x1l]=new fp_t**** [nx1*nx2] - x2l;
  p[x1l][x2l]=new fp_t*** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new fp_t** [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new fp_t* [nx1*nx2*nx3*nx4*nx5] - x5l;
  p[x1l][x2l][x3l][x4l][x5l]=new fp_t [nx1*nx2*nx3*nx4*nx5*nx6] - x6l;
  for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4l][x5]=p[x1l][x2l][x3l][x4l][x5-1]+nx6;
  for(int x4=x4l+1;x4<=x4h;++x4){
    p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;    
    p[x1l][x2l][x3l][x4][x5l]=p[x1l][x2l][x3l][x4-1][x5l]+nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4][x5]=p[x1l][x2l][x3l][x4][x5-1]+nx6;
  }
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    p[x1l][x2l][x3][x4l][x5l]=p[x1l][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4l][x5]=p[x1l][x2l][x3][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5;    
      p[x1l][x2l][x3][x4][x5l]=p[x1l][x2l][x3][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4][x5]=p[x1l][x2l][x3][x4][x5-1]+nx6;
    }
  }
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    p[x1l][x2][x3l][x4l][x5l]=p[x1l][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4l][x5]=p[x1l][x2][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;    
      p[x1l][x2][x3l][x4][x5l]=p[x1l][x2][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4][x5]=p[x1l][x2][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      p[x1l][x2][x3][x4l][x5l]=p[x1l][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4l][x5]=p[x1l][x2][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;    
        p[x1l][x2][x3][x4][x5l]=p[x1l][x2][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4][x5]=p[x1l][x2][x3][x4][x5-1]+nx6;
      }
    }
  }
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    p[x1][x2l][x3l][x4l][x5l]=p[x1-1][x2l][x3l][x4l][x5l]+nx2*nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4l][x5]=p[x1][x2l][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;    
      p[x1][x2l][x3l][x4][x5l]=p[x1][x2l][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4][x5]=p[x1][x2l][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      p[x1][x2l][x3][x4l][x5l]=p[x1][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4l][x5]=p[x1][x2l][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;    
        p[x1][x2l][x3][x4][x5l]=p[x1][x2l][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4][x5]=p[x1][x2l][x3][x4][x5-1]+nx6;
      }
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      p[x1][x2][x3l][x4l][x5l]=p[x1][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4l][x5]=p[x1][x2][x3l][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;    
        p[x1][x2][x3l][x4][x5l]=p[x1][x2][x3l][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4][x5]=p[x1][x2][x3l][x4][x5-1]+nx6;
      }
      for(int x3=x3l+1;x3<=x3h;++x3){
        p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
        p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
        p[x1][x2][x3][x4l][x5l]=p[x1][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4l][x5]=p[x1][x2][x3][x4l][x5-1]+nx6;
        for(int x4=x4l+1;x4<=x4h;++x4){
          p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;    
          p[x1][x2][x3][x4][x5l]=p[x1][x2][x3][x4-1][x5l]+nx5*nx6;
          for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4][x5]=p[x1][x2][x3][x4][x5-1]+nx6;
        }
      }
    }
  }
  return p;
}

void del_ft6dim(fp_t ******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  delete[] (p[x1l][x2l][x3l][x4l][x5l]+x6l);
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

fp_t ******* ft7dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h, int x7l, int x7h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1,nx6=x6h-x6l+1, nx7=x7h-x7l+1;
  
  // Actually allocate all the memory. 
  fp_t *******p;
  p=new fp_t****** [nx1] - x1l;
  p[x1l]=new fp_t***** [nx1*nx2] - x2l;
  p[x1l][x2l]=new fp_t**** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new fp_t*** [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new fp_t** [nx1*nx2*nx3*nx4*nx5] - x5l;
  p[x1l][x2l][x3l][x4l][x5l]=new fp_t* [nx1*nx2*nx3*nx4*nx5*nx6] - x6l;
  p[x1l][x2l][x3l][x4l][x5l][x6l]=new fp_t [nx1*nx2*nx3*nx4*nx5*nx6*nx7] - x7l;

  // We start with pointers to  ***** 

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1] = p[x1-1] + nx2;

  // Then the pointer to  ****

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l] = p[x1-1][x2l] + nx2 * nx3;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2] = p[x1][x2-1] + nx3;

  // We then continue further witn the pointers to ***

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l][x3l] = p[x1-1][x2l][x3l] + nx2 * nx3 * nx4;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2][x3l] = p[x1][x2-1][x3l] + nx3 * nx4;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l+1; x3<=x3h; ++x3)
        p[x1][x2][x3] = p[x1][x2][x3-1] + nx4;

  // Then the pointers to **

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l][x3l][x4l] = p[x1-1][x2l][x3l][x4l] + nx2 * nx3 * nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2][x3l][x4l] = p[x1][x2-1][x3l][x4l] + nx3 * nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l+1; x3<=x3h; ++x3)
        p[x1][x2][x3][x4l] = p[x1][x2][x3-1][x4l] + nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l; x3<=x3h; ++x3)
        for (int x4 = x4l+1; x4<=x4h; ++x4)
          p[x1][x2][x3][x4] = p[x1][x2][x3][x4-1] + nx5;

  // The next ones are the pointers to * 
  for (int x1 = x1l; x1<=x1h; ++x1){
    if (x1>x1l) p[x1][x2l][x3l][x4l][x5l] = p[x1-1][x2l][x3l][x4l][x5l] + nx2*nx3*nx4*nx5*nx6;
    for (int x2=x2l;x2<=x2h;++x2){
      if (x2>x2l) p[x1][x2][x3l][x4l][x5l] = p[x1][x2-1][x3l][x4l][x5l] + nx3*nx4*nx5*nx6;
      for (int x3=x3l;x3<=x3h;++x3){
        if (x3>x3l) p[x1][x2][x3][x4l][x5l] = p[x1][x2][x3-1][x4l][x5l] + nx4*nx5*nx6;
        for (int x4=x4l;x4<=x4h;++x4){
          if(x4>x4l) p[x1][x2][x3][x4][x5l] = p[x1][x2][x3][x4-1][x5l] + nx5*nx6;
          for (int x5=x5l+1;x5<=x5h;++x5)
            p[x1][x2][x3][x4][x5] = p[x1][x2][x3][x4][x5-1] + nx6;
        }
      }
    }
  }

  // And finally, the pointers to fp_t itself

  for (int x1=x1l;x1<=x1h;++x1){
    if (x1>x1l) p[x1][x2l][x3l][x4l][x5l][x6l] = p[x1-1][x2l][x3l][x4l][x5l][x6l] + nx2*nx3*nx4*nx5*nx6*nx7;
    for (int x2=x2l;x2<=x2h;++x2){
      if (x2>x2l) p[x1][x2][x3l][x4l][x5l][x6l] = p[x1][x2-1][x3l][x4l][x5l][x6l] + nx3*nx4*nx5*nx6*nx7;
      for (int x3=x3l;x3<=x3h;++x3){
        if (x3>x3l) p[x1][x2][x3][x4l][x5l][x6l] = p[x1][x2][x3-1][x4l][x5l][x6l] + nx4*nx5*nx6*nx7;
        for (int x4=x4l;x4<=x4h;++x4){
          if (x4>x4l) p[x1][x2][x3][x4][x5l][x6l] = p[x1][x2][x3][x4-1][x5l][x6l] + nx5*nx6*nx7;
          for (int x5=x5l;x5<=x5h;++x5){
            if (x5>x5l) p[x1][x2][x3][x4][x5][x6l] = p[x1][x2][x3][x4][x5-1][x6l] + nx6*nx7;
            for (int x6=x6l+1;x6<=x6h;++x6)
              p[x1][x2][x3][x4][x5][x6] = p[x1][x2][x3][x4][x5][x6-1] + nx7;
          }
        }
      }
    }
  }
  return p;
} 

void del_ft7dim(fp_t *******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h,int x7l, int x7h)
{
  delete[] (p[x1l][x2l][x3l][x4l][x5l][x6l]+x7l);
  delete[] (p[x1l][x2l][x3l][x4l][x5l]+x6l);
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
} 

fp_t ******** ft8dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,
  int x6l,int x6h, int x7l, int x7h, int x8l, int x8h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1,nx6=x6h-x6l+1, nx7=x7h-x7l+1, nx8=x8h-x8l+1;
  
  // Actually allocate all the memory. 
  fp_t ********p;
  p=new fp_t******* [nx1] - x1l;
  p[x1l]=new fp_t****** [nx1*nx2] - x2l;
  p[x1l][x2l]=new fp_t***** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new fp_t**** [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new fp_t*** [nx1*nx2*nx3*nx4*nx5] - x5l;
  p[x1l][x2l][x3l][x4l][x5l]=new fp_t** [nx1*nx2*nx3*nx4*nx5*nx6] - x6l;
  p[x1l][x2l][x3l][x4l][x5l][x6l]=new fp_t * [nx1*nx2*nx3*nx4*nx5*nx6*nx7] - x7l;
  p[x1l][x2l][x3l][x4l][x5l][x6l][x7l]=new fp_t [nx1*nx2*nx3*nx4*nx5*nx6*nx7*nx8] - x8l;

  // We start with pointers to  ****** 

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1] = p[x1-1] + nx2;

  // Then the pointer to  *****

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l] = p[x1-1][x2l] + nx2 * nx3;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2] = p[x1][x2-1] + nx3;

  // We then continue further witn the pointers to ****

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l][x3l] = p[x1-1][x2l][x3l] + nx2 * nx3 * nx4;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2][x3l] = p[x1][x2-1][x3l] + nx3 * nx4;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l+1; x3<=x3h; ++x3)
        p[x1][x2][x3] = p[x1][x2][x3-1] + nx4;

  // Then the pointers to ***

  for (int x1 = x1l+1; x1<=x1h; ++x1)
    p[x1][x2l][x3l][x4l] = p[x1-1][x2l][x3l][x4l] + nx2 * nx3 * nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l+1; x2<=x2h; ++x2)
      p[x1][x2][x3l][x4l] = p[x1][x2-1][x3l][x4l] + nx3 * nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l+1; x3<=x3h; ++x3)
        p[x1][x2][x3][x4l] = p[x1][x2][x3-1][x4l] + nx4 * nx5;
  for (int x1 = x1l; x1<=x1h; ++x1)
    for (int x2 = x2l; x2<=x2h; ++x2)
      for (int x3 = x3l; x3<=x3h; ++x3)
        for (int x4 = x4l+1; x4<=x4h; ++x4)
          p[x1][x2][x3][x4] = p[x1][x2][x3][x4-1] + nx5;

  // The next ones are the pointers to ** 
  for (int x1 = x1l; x1<=x1h; ++x1){
    if (x1>x1l) p[x1][x2l][x3l][x4l][x5l] = p[x1-1][x2l][x3l][x4l][x5l] + nx2*nx3*nx4*nx5*nx6;
    for (int x2=x2l;x2<=x2h;++x2){
      if (x2>x2l) p[x1][x2][x3l][x4l][x5l] = p[x1][x2-1][x3l][x4l][x5l] + nx3*nx4*nx5*nx6;
      for (int x3=x3l;x3<=x3h;++x3){
        if (x3>x3l) p[x1][x2][x3][x4l][x5l] = p[x1][x2][x3-1][x4l][x5l] + nx4*nx5*nx6;
        for (int x4=x4l;x4<=x4h;++x4){
          if(x4>x4l) p[x1][x2][x3][x4][x5l] = p[x1][x2][x3][x4-1][x5l] + nx5*nx6;
          for (int x5=x5l+1;x5<=x5h;++x5)
            p[x1][x2][x3][x4][x5] = p[x1][x2][x3][x4][x5-1] + nx6;
        }
      }
    }
  }

  // The pointers to fp_t *

  for (int x1=x1l;x1<=x1h;++x1){
    if (x1>x1l) p[x1][x2l][x3l][x4l][x5l][x6l] = p[x1-1][x2l][x3l][x4l][x5l][x6l] + nx2*nx3*nx4*nx5*nx6*nx7;
    for (int x2=x2l;x2<=x2h;++x2){
      if (x2>x2l) p[x1][x2][x3l][x4l][x5l][x6l] = p[x1][x2-1][x3l][x4l][x5l][x6l] + nx3*nx4*nx5*nx6*nx7;
      for (int x3=x3l;x3<=x3h;++x3){
        if (x3>x3l) p[x1][x2][x3][x4l][x5l][x6l] = p[x1][x2][x3-1][x4l][x5l][x6l] + nx4*nx5*nx6*nx7;
        for (int x4=x4l;x4<=x4h;++x4){
          if (x4>x4l) p[x1][x2][x3][x4][x5l][x6l] = p[x1][x2][x3][x4-1][x5l][x6l] + nx5*nx6*nx7;
          for (int x5=x5l;x5<=x5h;++x5){
            if (x5>x5l) p[x1][x2][x3][x4][x5][x6l] = p[x1][x2][x3][x4][x5-1][x6l] + nx6*nx7;
            for (int x6=x6l+1;x6<=x6h;++x6)
              p[x1][x2][x3][x4][x5][x6] = p[x1][x2][x3][x4][x5][x6-1] + nx7;
          }
        }
      }
    }
  }

  // And finally the pointers to fp_t
  for (int x1=x1l;x1<=x1h;++x1){
    if (x1>x1l) p[x1][x2l][x3l][x4l][x5l][x6l][x7l] = p[x1-1][x2l][x3l][x4l][x5l][x6l][x7l] + nx2*nx3*nx4*nx5*nx6*nx7*nx8;
    for (int x2=x2l;x2<=x2h;++x2){
      if (x2>x2l) p[x1][x2][x3l][x4l][x5l][x6l][x7l] = p[x1][x2-1][x3l][x4l][x5l][x6l][x7l] + nx3*nx4*nx5*nx6*nx7*nx8;
      for (int x3=x3l;x3<=x3h;++x3){
        if (x3>x3l) p[x1][x2][x3][x4l][x5l][x6l][x7l] = p[x1][x2][x3-1][x4l][x5l][x6l][x7l] + nx4*nx5*nx6*nx7*nx8;
        for (int x4=x4l;x4<=x4h;++x4){
          if (x4>x4l) p[x1][x2][x3][x4][x5l][x6l][x7l] = p[x1][x2][x3][x4-1][x5l][x6l][x7l] + nx5*nx6*nx7*nx8;
          for (int x5=x5l;x5<=x5h;++x5){
            if (x5>x5l) p[x1][x2][x3][x4][x5][x6l][x7l] = p[x1][x2][x3][x4][x5-1][x6l][x7l] + nx6*nx7*nx8;
            for (int x6=x6l;x6<=x6h;++x6){
              if (x6>x6l) p[x1][x2][x3][x4][x5][x6][x7l] = p[x1][x2][x3][x4][x5][x6-1][x7l] + nx7*nx8;
                for (int x7=x7l+1;x7<=x7h;++x7)
                  p[x1][x2][x3][x4][x5][x6][x7] = p[x1][x2][x3][x4][x5][x6][x7-1] + nx8;
            }
          }
        }
      }
    }
  }
  return p;
}

void del_ft8dim(fp_t ********p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,
  int x6l,int x6h,int x7l, int x7h, int x8l, int x8h)
{
  delete[] (p[x1l][x2l][x3l][x4l][x5l][x6l][x7l]+x8l);
  delete[] (p[x1l][x2l][x3l][x4l][x5l][x6l]+x7l);
  delete[] (p[x1l][x2l][x3l][x4l][x5l]+x6l);
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}  


complex_t **ct2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  complex_t **p;
  p=new complex_t* [nx1] - x1l;
  p[x1l]=new complex_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_ct2dim(complex_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

complex_t ***ct3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  complex_t ***p;
  p=new complex_t** [nx1] - x1l;
  p[x1l]=new complex_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new complex_t [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_ct3dim(complex_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int08_t **i08t2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int08_t **p;
  p=new int08_t* [nx1] - x1l;
  p[x1l]=new int08_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_i08t2dim(int08_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int16_t **i16t2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int16_t **p;
  p=new int16_t* [nx1] - x1l;
  p[x1l]=new int16_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

int16_t **i16t2dim(int16_t *data,int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int16_t **p;
  p=new int16_t* [nx1] - x1l;
  p[x1l]=data - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_i16t2dim(int16_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int16_t ***i16t3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int16_t ***p;
  p=new int16_t** [nx1] - x1l;
  p[x1l]=new int16_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new int16_t [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

int16_t ***i16t3dim(int16_t *data,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int16_t ***p;
  p=new int16_t** [nx1] - x1l;
  p[x1l]=new int16_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=data - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_i16t3dim(int16_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

//

int16_t **i16t2dim(int x1l,int x1h,int x2l,int *x2h)
{
  int nx1=x1h-x1l+1,nx12=0;
  int16_t **p;
  p=new int16_t* [nx1] - x1l;
  for(int x1=x1l;x1<=x1h;++x1) nx12+=(x2h[x1]-x2l+1);
  p[x1l]=new int16_t [nx12] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+(x2h[x1-1]-x2l+1);
  return p;
}

void del_i16t2dim(int16_t **p,int x1l,int x1h,int x2l,int *x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

//

int32_t **i32t2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int32_t **p;
  p=new int32_t* [nx1] - x1l;
  p[x1l]=new int32_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

int32_t **i32t2dim(int32_t *data,int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int32_t **p;
  p=new int32_t* [nx1] - x1l;
  p[x1l]=data - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_i32t2dim(int32_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}


int32_t ***i32t3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int32_t ***p;
  p=new int32_t** [nx1] - x1l;
  p[x1l]=new int32_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new int32_t [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

int32_t ***i32t3dim(int32_t *data,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int32_t ***p;
  p=new int32_t** [nx1] - x1l;
  p[x1l]=new int32_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=data - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_i32t3dim(int32_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

//

uint32_t **ui32t2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  uint32_t **p;
  p=new uint32_t* [nx1] - x1l;
  p[x1l]=new uint32_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

uint32_t **ui32t2dim(uint32_t *data,int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  uint32_t **p;
  p=new uint32_t* [nx1] - x1l;
  p[x1l]=data - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_ui32t2dim(uint32_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}


uint32_t ***ui32t3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  uint32_t ***p;
  p=new uint32_t** [nx1] - x1l;
  p[x1l]=new uint32_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new uint32_t [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

uint32_t ***ui32t3dim(uint32_t *data,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  uint32_t ***p;
  p=new uint32_t** [nx1] - x1l;
  p[x1l]=new uint32_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=data - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_ui32t3dim(uint32_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

//


int64_t **i64t2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int64_t **p;
  p=new int64_t* [nx1] - x1l;
  p[x1l]=new int64_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

int64_t **i64t2dim(int64_t *data,int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int64_t **p;
  p=new int64_t* [nx1] - x1l;
  p[x1l]=data - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_i64t2dim(int64_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}


int64_t ***i64t3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int64_t ***p;
  p=new int64_t** [nx1] - x1l;
  p[x1l]=new int64_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new int64_t [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

int64_t ***i64t3dim(int64_t *data,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int64_t ***p;
  p=new int64_t** [nx1] - x1l;
  p[x1l]=new int64_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=data - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_i64t3dim(int64_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

//

float32_t **f32t2dim(int x1l,int x1h,int x2l,int *x2h)
{
  int nx1=x1h-x1l+1,nx12=0;
  float32_t **p;
  p=new float32_t* [nx1] - x1l;
  for(int x1=x1l;x1<=x1h;++x1) nx12+=(x2h[x1]-x2l+1);
  p[x1l]=new float32_t [nx12] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+(x2h[x1-1]-x2l+1);
  return p;
}

void del_f32t2dim(float32_t **p,int x1l,int x1h,int x2l,int *x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

float32_t **f32t2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  float32_t **p;
  p=new float32_t* [nx1] - x1l;
  p[x1l]=new float32_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

float32_t **f32t2dim(float32_t *data,int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  float32_t **p;
  p=new float32_t* [nx1] - x1l;
  p[x1l]=data - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_f32t2dim(float32_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

float32_t ***f32t3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  float32_t ***p;
  p=new float32_t** [nx1] - x1l;
  p[x1l]=new float32_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new float32_t [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

float32_t ***f32t3dim(float32_t *data,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  float32_t ***p;
  p=new float32_t** [nx1] - x1l;
  p[x1l]=new float32_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=data - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_f32t3dim(float32_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

float32_t ****f32t4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1;
  float32_t ****p;
  p=new float32_t*** [nx1] - x1l;
  p[x1l]=new float32_t** [nx1*nx2] - x2l;
  p[x1l][x2l]=new float32_t* [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new float32_t [nx1*nx2*nx3*nx4] - x4l;
  for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  return p;
}

void del_f32t4dim(float32_t ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

//

float64_t **f64t2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  float64_t **p;
  p=new float64_t* [nx1] - x1l;
  p[x1l]=new float64_t [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

float64_t **f64t2dim(float64_t *data,int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  float64_t **p;
  p=new float64_t* [nx1] - x1l;
  p[x1l]=data - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_f64t2dim(float64_t **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

float64_t ***f64t3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  float64_t ***p;
  p=new float64_t** [nx1] - x1l;
  p[x1l]=new float64_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=new float64_t [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

float64_t ***f64t3dim(float64_t *data,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  float64_t ***p;
  p=new float64_t** [nx1] - x1l;
  p[x1l]=new float64_t* [nx1*nx2] - x2l;
  p[x1l][x2l]=data - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_f64t3dim(float64_t ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}


//

fp_t **ft2dim(fp_t *t,int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  fp_t **p=new fp_t* [nx1]-x1l;
  p[x1l]=t-x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

fp_t ***ft3dim(fp_t *t,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  fp_t ***p=new fp_t** [nx1]-x1l;
  p[x1l]=new fp_t* [nx1*nx2]-x2l;
  p[x1l][x2l]=t-x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

fp_t ****ft4dim(fp_t *t,int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx23t=0,nx4=x4h-x4l+1;
  fp_t ****p;
  p=new fp_t*** [nx1] - x1l;
  p[x1l]=new fp_t** [nx1*nx2] - x2l;
  int *nx3=new int [nx2]-x2l;
  for(int x2=x2l;x2<=x2h;++x2) nx23t+=(nx3[x2]=x3h[x2]-x3l+1); // equivalent to nx2*nx3
  p[x1l][x2l]=new fp_t* [nx1*nx23t] - x3l;
  p[x1l][x2l][x3l]=t - x4l;
  for(int x3=x3l+1;x3<=x3h[x2l];++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3[x2-1];
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3[x2-1]*nx4;
    for(int x3=x3l+1;x3<=x3h[x2];++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx23t;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx23t*nx4;
    for(int x3=x3l+1;x3<=x3h[x2l];++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3[x2-1];
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3[x2-1]*nx4;
      for(int x3=x3l+1;x3<=x3h[x2];++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  delete[] (nx3+x2l);
  return p;
}

byte **b2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  byte **p;
  p=new byte* [nx1] - x1l;
  p[x1l]=new byte [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_b2dim(byte **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

byte ***b3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  byte ***p;
  p=new byte** [nx1] - x1l;
  p[x1l]=new byte* [nx1*nx2] - x2l;
  p[x1l][x2l]=new byte [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_b3dim(byte ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

byte **b2dim(int x1l,int x1h,int x2l,int *x2h)
{
  int nx1=x1h-x1l+1,nx12=0;
  byte **p;
  p=new byte* [nx1] - x1l;
  for(int x1=x1l;x1<=x1h;++x1) nx12+=(x2h[x1]-x2l+1);
  p[x1l]=new byte [nx12] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+(x2h[x1-1]-x2l+1);
  return p;
}

void del_b2dim(byte **p,int x1l,int x1h,int x2l,int *x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

