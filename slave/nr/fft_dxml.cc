#include <dxmldef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "types.h"
#include "mem.h"
#include "fft.h"

void fft_init(int,int)
{
}

void fft_2d(complex_t **&data,int n1,int n2,int isign)
{
  int lda=n1,sx=1,sy=1;
  complex_t **res=ct2dim(1,n1,1,n2);
  char dir[3]="f";
  if(isign<0) dir[0]='b';
  switch(sizeof(fp_t)){
    case(8): zfft_2d("c","c",dir,data[1]+1,res[1]+1,&n1,&n2,&lda,&sx,&sy); break;
    case(4): cfft_2d("c","c",dir,data[1]+1,res[1]+1,&n1,&n2,&lda,&sx,&sy); break;
    default: exit(fprintf(stderr,"error: unknown floating point type (size=%d)\n",sizeof(fp_t)));
  }
  del_ct2dim(data,1,n1,1,n2);
  data=res;
}

void fft_done(int,int)
{
}

void fft_reorder(complex_t **f,int np)
{
  int nh=np/2;
  complex_t *buf=new complex_t [nh];
  for(int x=1;x<=nh;++x){
    memcpy(buf,f[x]+1,nh*sizeof(complex_t));
    memcpy(f[x]+1,f[x+nh]+nh+1,nh*sizeof(complex_t));
    memcpy(f[x+nh]+nh+1,buf,nh*sizeof(complex_t));
//
    memcpy(buf,f[x+nh]+1,nh*sizeof(complex_t));
    memcpy(f[x+nh]+1,f[x]+nh+1,nh*sizeof(complex_t));
    memcpy(f[x]+nh+1,buf,nh*sizeof(complex_t));
  }
  delete buf;
}


