#include <math.h>
#include <string.h>

#include "types.h"
#include "fft.h"

void swap(fp_t &a,fp_t &b)
{
  fp_t c=a;
  a=b;
  b=c;
}

void fft(complex_t *data_in,int nn,int isign)
{
  fp_t *data=(fp_t*)data_in;
  int n=nn<<1,j=1;
  for(int i=1;i<n;i+=2){
    if(j>i){
      swap(data[j],data[i]);
      swap(data[j+1],data[i+1]);
    }
    int m=n>>1;
    while((m>=2)&&(j>m)){
      j-=m;
      m>>=1;
    }
    j+=m;
  }
  int mmax=2;
  while(n>mmax){
    int istep=mmax<<1;
    fp_t theta=isign*(6.28318530717959/mmax);
    fp_t wtemp=sin(0.5*theta);
    fp_t wpr=-2.0*wtemp*wtemp;
    fp_t wpi=sin(theta);
    fp_t wr=1.0;
    fp_t wi=0.0;
    for(int m=1;m<mmax;m+=2){
      for(int i=m;i<=n;i+=istep){
        j=i+mmax;
        fp_t tempr=wr*data[j]-wi*data[j+1];
        fp_t tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i]+=tempr;
        data[i+1]+=tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

#include <stdio.h>

void fft_n(complex_t *data_in,int *nn,int nd,int isign)
{ 
  fp_t *data=(fp_t*)&(data_in[1])-1;
  int ip1,ip2,ip3,nt=1;
  for(int id=1;id<=nd;++id) nt*=nn[id];
  int nprev=1;
  for(int id=nd;id>=1;--id){
    int n=nn[id];
    int nrem=nt/(n*nprev);
    ip1=nprev<<1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    int i2rev=1;
    for(int i2=1;i2<=ip2;i2+=ip1){
      if(i2<i2rev)
	for(int i1=i2;i1<=i2+ip1-2;i1+=2)
	  for(int i3=i1;i3<=ip3;i3+=ip2){
            int i3rev=i2rev+i3-i2;
	    swap(data[i3],data[i3rev]);
	    swap(data[i3+1],data[i3rev+1]);
	  }
      int ibit=ip2>>1;
      while(ibit>=ip1&&i2rev>ibit){
	i2rev-=ibit;
	ibit>>=1;
      }
      i2rev+=ibit;
    }
    int ifp1=ip1;
    while(ifp1<ip2){
      int ifp2=ifp1<<1;
      fp_t theta=isign*6.28318530717959/(ifp2/ip1);
      fp_t wtemp=sin(0.5*theta);
      fp_t wpr= -2.0*wtemp*wtemp;
      fp_t wpi=sin(theta);
      fp_t wr=1.0;
      fp_t wi=0.0;
      for(int i3=1;i3<=ifp1;i3+=ip1){
	for(int i1=i3;i1<=i3+ip1-2;i1+=2)
	  for(int i2=i1;i2<=ip3;i2+=ifp2){
	    int k1=i2;
	    int k2=k1+ifp1;
	    fp_t tempr=(fp_t)wr*data[k2]-(fp_t)wi*data[k2+1];
	    fp_t tempi=(fp_t)wr*data[k2+1]+(fp_t)wi*data[k2];
	    data[k2]=data[k1]-tempr;
	    data[k2+1]=data[k1+1]-tempi;
	    data[k1]+=tempr;
	    data[k1+1]+=tempi;
	  }
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev*=n;
  }
}

void fft_init(int,int)
{
}

void fft_2d(complex_t **&data,int n1,int n2,int isign)
{
  int nn[3]={0,n1,n2};
  fft_n(data[1],nn,2,isign);
  if(isign<0){
    fp_t n=(fp_t)(n1*n2);
    for(int x1=1;x1<=n1;++x1)
      for(int x2=1;x2<=n2;++x2){
        data[x1][x2].re/=n;
        data[x1][x2].im/=n;
      }
  }
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

void fft_reorder_2d(complex_t **f,int nx,int ny)
{
  int nhx=nx/2,nhy=ny/2;
  complex_t *buf=new complex_t [nhy];
  for(int x=1;x<=nhx;++x){
    memcpy(buf,f[x]+1,nhy*sizeof(complex_t));
    memcpy(f[x]+1,f[x+nhx]+nhy+1,nhy*sizeof(complex_t));
    memcpy(f[x+nhx]+nhy+1,buf,nhy*sizeof(complex_t));
//
    memcpy(buf,f[x+nhx]+1,nhy*sizeof(complex_t));
    memcpy(f[x+nhx]+1,f[x]+nhy+1,nhy*sizeof(complex_t));
    memcpy(f[x]+nhy+1,buf,nhy*sizeof(complex_t));
  }
  delete[] buf;
}


