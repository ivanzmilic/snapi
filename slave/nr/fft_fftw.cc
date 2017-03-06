#include <fftw3.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "types.h"
#include "mem.h"
#include "fft.h"

//
// WARNING!!! fftw returns the FT^* and not the FT (i.e. the complex 
// part has different sign). 
//

struct fft_plan_str{
  fftw_plan fftf_plan,fftb_plan;
  int n1,n2;
};

struct fft_plan_str **plans=0;

void fft_init(int n1,int n2)
{
  if(!plans){ // first call to fft_init
    plans=new fft_plan_str* [1];
    plans[0]=0;
  }
  for(int p=0;plans[p];++p)
    if((plans[p]->n1==n1)&&(plans[p]->n2==n2)) return;
  int q;
  for(q=0;plans[q];++q);
  fft_plan_str **tmp=new fft_plan_str* [q+2];
  memcpy(tmp,plans,q*sizeof(fft_plan_str*));
  delete[] plans;
  plans=tmp;
  plans[q]=new fft_plan_str;
  plans[q]->n1=n1;
  plans[q]->n2=n2;
  complex_t **in=ct2dim(1,n1,1,n2);
  complex_t **out=ct2dim(1,n1,1,n2);
  plans[q]->fftf_plan=fftw_plan_dft_2d(n1,n2,(fftw_complex*)(in[1]+1),(fftw_complex*)(out[1]+1),FFTW_FORWARD,FFTW_MEASURE);
  plans[q]->fftb_plan=fftw_plan_dft_2d(n1,n2,(fftw_complex*)(in[1]+1),(fftw_complex*)(out[1]+1),FFTW_BACKWARD,FFTW_MEASURE);  
  del_ct2dim(out,1,n1,1,n2);
  del_ct2dim(in,1,n1,1,n2);
  plans[q+1]=0;
}

void fft_2d(complex_t **&data,int n1,int n2,int isign)
{
  for(int p=0;plans[p];++p)
    if((plans[p]->n1==n1)&&(plans[p]->n2==n2)){
      complex_t **res=ct2dim(1,n1,1,n2);
      if(isign>0){
  	fftw_execute_dft(plans[p]->fftf_plan,(fftw_complex*)(data[1]+1),(fftw_complex*)(res[1]+1));
  	for(int x=1;x<=n1;++x)
  	  for(int y=1;y<=n2;++y) res[x][y].im=-res[x][y].im;
      }else{ // inverse
  	for(int x=1;x<=n1;++x)
  	  for(int y=1;y<=n2;++y) data[x][y].im=-data[x][y].im;
  	fftw_execute_dft(plans[p]->fftb_plan,(fftw_complex*)(data[1]+1),(fftw_complex*)(res[1]+1));
  	fp_t N=(fp_t)(n1*n2);
  	for(int x=1;x<=n1;++x)
  	  for(int y=1;y<=n2;++y){
  	    res[x][y].re/=N;
  	    res[x][y].im/=N;
  	  }
      }
      del_ct2dim(data,1,n1,1,n2);
      data=res;
      return;
    }
  fprintf(stderr,"Error: dimension mismatch (%dx%d) in fft (forgot to call fft_init)?\n",n1,n2);
}

void fft_done(int n1,int n2)
{
  if(plans){
    for(int p=0;plans[p];++p)
      if((plans[p]->n1==n1)&&(plans[p]->n2==n2)){
        fftw_destroy_plan(plans[p]->fftf_plan);
        fftw_destroy_plan(plans[p]->fftb_plan);
        int q;
        for(q=0;plans[q];++q);
        fft_plan_str **tmp=new fft_plan_str* [q];
        memcpy(tmp,plans,p*sizeof(fft_plan_str*));
        memcpy(tmp+p,plans+p+1,(q-p)*sizeof(fft_plan_str*));
        delete plans[p];
        delete[] plans;
        plans=tmp;
        p=q-2; // break the loop
      }
    if(!plans[0]){ // last call to fft_done
      delete[] plans;
      plans=0;
    }
  }
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
  delete[] buf;
}


