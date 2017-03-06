#include <sys/types.h>
#include <sys/select.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "types.h"
#include "uts.h"

void delay(fp_t t) // Specify the delay time in s!
{
  struct timeval tv;
  tv.tv_sec=(time_t)t;
  tv.tv_usec=(int)((t-(fp_t)tv.tv_sec)*1.0E+06); 
  select(0,NULL,NULL,NULL,&tv);
}

void swap(char &a,char &b)
{
  char c=a;
  a=b;
  b=c;
}

void swap(int64_t *x,int n)
{
  int sz=sizeof(int64_t);
  int m=sz/2;
  sz-=1;
  for(int i=0;i<=n-1;++i){
    char *p=(char*)(x+i);
    for(int j=0;j<=m-1;++j) swap(p[j],p[sz-j]);
  }
}

void swap(int32_t *x,int n)
{
  int sz=sizeof(int32_t);
  int m=sz/2;
  sz-=1;
  for(int i=0;i<=n-1;++i){
    char *p=(char*)(x+i);
    for(int j=0;j<=m-1;++j) swap(p[j],p[sz-j]);
  }
}

void swap(int16_t *x,int n)
{
  int sz=sizeof(int16_t);
  int m=sz/2;
  sz-=1;
  for(int i=0;i<=n-1;++i){
    char *p=(char*)(x+i);
    for(int j=0;j<=m-1;++j) swap(p[j],p[sz-j]);
  }
}

void swap(void *x,int sz,int n)
{
  int m=sz/2-1;
  int szm=sz-1;
  for(int i=0;i<=sz*(n-1);i+=sz){
    char *p=(char*)x+i;
    for(int j=0;j<=m;++j) swap(p[j],p[szm-j]);
  }
}

void sort(int *x,int nx) // inefficient algorithm, but quick and easy and nx is small anyway
{
  for(int i=1;i<=nx-1;++i)
    for(int j=1;j<=nx-i;++j)
      if(x[j]>x[j+1]) swap(x[j],x[j+1]);
}

void sort(fp_t *x,int nx) // inefficient algorithm, but quick and easy and nx is small anyway
{
  for(int i=1;i<=nx-1;++i)
    for(int j=1;j<=nx-i;++j)
      if(x[j]>x[j+1]) swap(x[j],x[j+1]);
}

void sort(int *x,int *y,int nx) // inefficient algorithm, but quick and easy and nx is small anyway
{
  for(int i=1;i<=nx-1;++i)
    for(int j=1;j<=nx-i;++j)
      if(x[j]>x[j+1]){
        swap(x[j],x[j+1]);
        swap(y[j],y[j+1]);
      }
}

void sort(fp_t *x,int *y,int nx) // inefficient algorithm, but quick and easy and nx is small anyway
{
  for(int i=1;i<=nx-1;++i)
    for(int j=1;j<=nx-i;++j)
      if(x[j]>x[j+1]){
        swap(x[j],x[j+1]);
        swap(y[j],y[j+1]);
      }
}

int get_numbers(const char *numbers,int *&img_num,int &ni)
{
  char *nums=strcpy(new char [strlen(numbers)+1],numbers);
  char **elems=new char* [2];
  elems[0]=nums;
  elems[1]=0;
  int n=0;
  while(char *p=strchr(elems[n],',')){
    char **tmp=new char* [n+3];
    memcpy(tmp,elems,(n+1)*sizeof(char*));
    delete[] elems;
    elems=tmp;
    elems[++n]=p+1;
    p[0]=0;
    elems[n+1]=0;
  }
  ni=0;
  img_num=new int [1] -1;
  for(int m=0;elems[m];++m)
    if(char *p=strchr(elems[m],'-')){ // element is range
      p[0]=0;
      int nmin,nmax;
      if(!sscanf(elems[m],"%d",&nmin)){
        fprintf(stderr,"error: could not convert \"%s\" to integer\n",elems[m]);
        return -1;
      }
      if(!sscanf(p+1,"%d",&nmax)){
        fprintf(stderr,"could not convert \"%s\" to integer\n",p+1);
        return -1;
      }
      if(nmin>nmax) swap(nmin,nmax);
      int nn=nmax-nmin+1;
      int *tmp=new int [ni+nn]-1;
      memcpy(tmp+1,img_num+1,ni*sizeof(int));
      delete[] (img_num+1);
      img_num=tmp;
      for(int i=nmin;i<=nmax;++i) img_num[i-nmin+ni+1]=i;
      ni+=nn;
    }else{                            // element is number (I hope)
      int nn;
      if(!sscanf(elems[m],"%d",&nn)){
        fprintf(stderr,"error: could not convert \"%s\" to integer\n",elems[m]);
        return -1;
      }
      int *tmp=new int [ni+1]-1;
      memcpy(tmp+1,img_num+1,ni*sizeof(int));
      delete[] (img_num+1);
      img_num=tmp;
      img_num[++ni]=nn;
    }
  delete[] elems;
  delete[] nums;
//  if(dosort) sort(img_num,ni);
  return ni;
}


int get_numbers(const char *numbers,fp_t *&nout,int &ni)
{
  char *nums=strcpy(new char [strlen(numbers)+1],numbers);
  char **elems=new char* [2];
  elems[0]=nums;
  elems[1]=0;
  int n=0;
  while(char *p=strchr(elems[n],',')){
    char **tmp=new char* [n+3];
    memcpy(tmp,elems,(n+1)*sizeof(char*));
    delete[] elems;
    elems=tmp;
    elems[++n]=p+1;
    p[0]=0;
    elems[n+1]=0;
  }
  ni=0;
  nout=new fp_t [1] -1;
  for(int m=0;elems[m];++m){
    fp_t nn;
    if(!sscanf(elems[m],"%lE",&nn)){
      fprintf(stderr,"error: could not convert \"%s\" to floating point number\n",elems[m]);
      return -1;
    }
    fp_t *tmp=new fp_t [ni+1]-1;
    if(ni) memcpy(tmp+1,nout+1,ni*sizeof(fp_t));
    delete[] (nout+1);
    nout=tmp;
    nout[++ni]=nn;
  }
  delete[] elems;
  delete[] nums;
  return ni;
}

int get_numbers(const char *numbers,fp_t *&nout,int &ni,fp_t dn)
{
  char *nums=strcpy(new char [strlen(numbers)+1],numbers);
  char **elems=new char* [2];
  elems[0]=nums;
  elems[1]=0;
  int n=0;
  while(char *p=strchr(elems[n],',')){
    char **tmp=new char* [n+3];
    memcpy(tmp,elems,(n+1)*sizeof(char*));
    delete[] elems;
    elems=tmp;
    elems[++n]=p+1;
    p[0]=0;
    elems[n+1]=0;
  }
  ni=0;
  nout=new fp_t [1] -1;
  for(int m=0;elems[m];++m){
    if(char *p=strchr(elems[m],'-')){ // element is range
      p[0]=0;
      fp_t nmin,nmax;
      if(!sscanf(elems[m],"%lE",&nmin)){
        fprintf(stderr,"error: could not convert \"%s\" to floating point value\n",elems[m]);
        return -1;
      }
      if(!sscanf(p+1,"%lE",&nmax)){
        fprintf(stderr,"could not convert \"%s\" to floating point value\n",p+1);
        return -1;
      }
      int nn=0.5+(nmax-nmin)/dn;
      fp_t ndn=(nmax-nmin)/(fp_t)nn;
      fp_t *tmp=new fp_t [ni+nn+1]-1;
      memcpy(tmp+1,nout+1,ni*sizeof(fp_t));
      delete[] (nout+1);
      nout=tmp;
      nout[ni+1]=nmin;
      for(int i=ni+2;i<ni+nn;++i) nout[i]=nout[i-1]+ndn;
      nout[ni+=nn]=nmax;
    }else{                            // element is number (I hope)
      fp_t nn;
      if(!sscanf(elems[m],"%lE",&nn)){
        fprintf(stderr,"error: could not convert \"%s\" to floating point value\n",elems[m]);
        return -1;
      }
      fp_t *tmp=new fp_t [ni+1]-1;
      if(ni) memcpy(tmp+1,nout+1,ni*sizeof(fp_t));
      delete[] (nout+1);
      nout=tmp;
      nout[++ni]=nn;
    }
  }
  delete[] elems;
  delete[] nums;
  return ni;
}


int get_number(const char *number,int &num)
{
  if(!sscanf(number,"%d",&num)){
    fprintf(stderr,"error: could not convert \"%s\" to integer\n",number);
    return -1;
  }
  return 0;
}

int get_number(const char *number,fp_t &num)
{
  if(!sscanf(number,"%lE",&num)){
    fprintf(stderr,"error: could not convert \"%s\" to floating point number\n",number);
    return -1;
  }
  return 0;
}

int contains(int *in,int nn,int num)
{
  for(int n=1;n<=nn;++n) if(in[n]==num) return 1;
  return 0;
}

int where(int *in,int nn,int num)
{
  for(int n=1;n<=nn;++n) if(in[n]==num) return n;
  return -1;
}

static __inline int pred(int k,int L,int M)
{
  return (k%M)*L+k/M;
}

void in_place_transpose(float32_t *data,int L,int M)
{
  int LxM=L*M;
  for(int i=0,stillToMove=LxM;stillToMove>0;++i){
    int k,j;
    for(j=pred(i,L,M);j>i;j=pred(j,L,M));
    if(j<i) continue;
    for(k=i,j=pred(i,L,M);j!=i;k=j,j=pred(j,L,M)){
      swap(data[k],data[j]);
      --stillToMove;
    }
    --stillToMove;
  }
}

void in_place_transpose(int32_t *data,int L,int M)
{
  int LxM=L*M;
  for(int i=0,stillToMove=LxM;stillToMove>0;++i){
    int k,j;
    for(j=pred(i,L,M);j>i;j=pred(j,L,M));
    if(j<i) continue;
    for(k=i,j=pred(i,L,M);j!=i;k=j,j=pred(j,L,M)){
      swap(data[k],data[j]);
      --stillToMove;
    }
    --stillToMove;
  }
} 

void in_place_transpose(int16_t *data,int L,int M)
{
  int LxM=L*M;
  for(int i=0,stillToMove=LxM;stillToMove>0;++i){
    int k,j;
    for(j=pred(i,L,M);j>i;j=pred(j,L,M));
    if(j<i) continue;
    for(k=i,j=pred(i,L,M);j!=i;k=j,j=pred(j,L,M)){
      swap(data[k],data[j]);
      --stillToMove;
    }
    --stillToMove;
  }
} 
