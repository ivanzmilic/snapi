#include <string.h>
#include <errno.h>
#include "types.h"
#include "io.h"
#include "fileio.h"
#include "ana_io.h"
#include "obs.h"
#include "mem.h"
#include "pack.h"

observable::observable(int ns_in):ns(ns_in)
{
  S=0;
  nlambda=0;
  nx=0;ny=0;
}

observable::observable(int nx_in,int ny_in,int ns_in)
{
  S=0;
  nlambda=0;
  ns=ns_in;
  nx=nx_in;ny=ny_in;
}
observable::observable(int nx_in,int ny_in,int ns_in, int nlambda_in)
{
  
  nlambda=nlambda_in;
  ns=ns_in;
  nx=nx_in;ny=ny_in;
  S=ft4dim(1,nx,1,ny,1,ns,1,nlambda);
  memset(S[1][1][1]+1,0,nx*ny*ns*nlambda*sizeof(fp_t));
  lambda = new fp_t[nlambda]-1;
  memset(lambda+1,0,nlambda*sizeof(fp_t));
}


observable::~observable(void){
  if(nlambda){
    del_ft4dim(S,1,nx,1,ny,1,ns,1,nlambda);
    delete[] (lambda+1);
  }
}

int32_t observable::size(io_class &io_in){

  int32_t sz = 4*sizeof(int); // ns, nlambda,nx,ny
  sz += nlambda*sizeof(fp_t); // lambda
  sz += nx*ny*nlambda*ns*sizeof(fp_t); // actual observation
  return sz;
}

int32_t observable::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in){

  int offs=::pack(buf+offs,nx,do_swap);
  offs+=::pack(buf+offs,ny,do_swap);
  offs+=::pack(buf+offs,ns,do_swap);
  offs+=::pack(buf+offs,nlambda,do_swap);
  
  offs+=::pack(buf+offs,lambda,1,nlambda,do_swap);
  offs+=::pack(buf+offs,S,1,nx,1,ny,1,ns,1,nlambda,do_swap);
}

void observable::add(fp_t *S_in,fp_t lambda_in)
{
  /*fp_t **tmp=ft2dim(1,ns,1,nlambda+1);
  fp_t *tmpl=new fp_t[nlambda+1]-1;
  if(nlambda){
    for(int s=1;s<=ns;++s) memcpy(tmp[s]+1,S[s]+1,nlambda*sizeof(fp_t));
    memcpy(tmpl+1,lambda+1,nlambda*sizeof(fp_t));
    del_ft2dim(S,1,ns,1,nlambda);
    delete[] (lambda+1);
  }
  S=tmp;
  lambda=tmpl;
  ++nlambda;
  for(int s=1;s<=ns;++s) S[s][nlambda]=S_in[s];
  lambda[nlambda]=lambda_in;*/ // Not sure I like this done this way.
}

void observable::set(fp_t * S_in, fp_t lambda_in, int i, int j, int l){

  lambda[l] = lambda_in;
  for (int s=1;s<=ns;++s)
    S[i][j][s][l] = S_in[s];
}

void observable::set(fp_t **** S_in){
  memcpy(S[1][1][1]+1,S_in[1][1][1]+1,nx*ny*ns*nlambda*sizeof(fp_t));
  /*for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      for (int s=1;s<=ns;++s)
        for (int l=1;l<=nlambda;++l)
          S[i][j][s][l] = S_in[l][s][j][i];*/
}

void observable::setlambda(fp_t * lambda_in){
  memcpy(lambda+1,lambda_in+1,nlambda*sizeof(fp_t));
}

fp_t **** observable::get_S(){

  fp_t **** S_copy;
  S_copy = ft4dim(1,nx,1,ny,1,ns,1,nlambda);
  memcpy(S_copy[1][1][1]+1,S[1][1][1]+1,ns*nlambda*nx*ny*sizeof(fp_t));

  return S_copy;
}

fp_t ** observable::get_S(int i, int j){

  fp_t ** S_copy;
  S_copy = ft2dim(1,ns,1,nlambda);
  memcpy(S_copy[1]+1,S[i][j][1]+1,nlambda*ns*sizeof(fp_t)); 

  return S_copy;
}

fp_t * observable::get_lambda(){
  fp_t * lambda_copy;
  lambda_copy = new fp_t [nlambda]-1;
  //printf("I am trying to copy an array of length %d \n", nlambda);
  memcpy(lambda_copy+1,lambda+1,nlambda*sizeof(fp_t));
  //printf("Does it work?\n");
  return lambda_copy;
}

int observable::get_n_lambda(){
  return nlambda;
}

void observable::write(const char *name,io_class &io,int i, int j)
{
  if(FILE *f=fopen(name,"w")){
    for(int l=1;l<=nlambda;++l){
      fprintf(f, "%.15e ", lambda[l]);
      for(int s=1;s<=ns;++s) fprintf(f, "%.15e ",S[i][j][s][l]);
      fprintf(f, "\n");
    }
    fclose(f);
  }else io.msg(IOL_WARN,"failed to open file \"%s\":\"%s\"\n",name,strerror(errno));
}

observable * observable::extract(int xl,int xh, int yl, int yh, int ll, int lh){

  int nxn = xh-xl+1,nyn=yh-yl+1,nln=lh-ll+1;
  observable * obs_small = new observable(nxn,nyn,ns,nln);
  fp_t **** S_small = ft4dim(1,nxn,1,nyn,1,ns,1,nln);
  for (int i=xl;i<=xh;++i)
    for (int j=yl;j<=yh;++j)
      for (int s=1;s<=ns;++s)
        for (int l=ll;l<=lh;++l)
          S_small[i-xl+1][j-yl+1][s][l-ll+1] = S[i][j][s][l];
  obs_small->set(S_small);

  del_ft4dim(S_small,1,nxn,1,nyn,1,ns,1,nln);

  fp_t * lambda_small = new fp_t [nln]-1;
  for (int l=ll;l<=lh;++l)
    lambda_small[l-ll+1] = lambda[l];
  obs_small->setlambda(lambda_small);
  delete[](lambda_small+1);
  return obs_small;

}

void observable::normalize(){
  // Normalizes already arranged observable to physical units.

  fp_t qs = 3.076E14; // quiet sun reference continuum

  // Average continuum in all field, assume here for simplicity that point #30 is continuum:
  int l_ref = 30;
  fp_t mean = 0.0;
  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      mean += S[i][j][1][30] / nx / ny;

  fp_t scale = qs/mean;

  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      for (int l=1;l<=nlambda;++l){
        S[i][j][1][l] *= scale;
        for (int s=2;s<=4;++s)
          S[i][j][s][l] *= S[i][j][1][l];
      }
}

void observable::read(char * name, io_class &io){

  // First, delete the arrays if they exist:
  //if (S)
  //  del_ft4dim(S,1,nx,1,ny,1,ns,1,nlambda);

  // WHY DOES THIS NOT WORK?
  //int n1,n2,n3,n4;
  //fp_t **** test = read_file(name,n1,n2,n3,n4,io);
  //S = read_file(name,nx,ny,ns,nlambda,io);
  //printf("I read an array with dimensions : %d %d %d %d \n", nx,ny,ns,nlambda);
}
