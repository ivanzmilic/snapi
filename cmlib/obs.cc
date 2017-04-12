#include <string.h>
#include <errno.h>
#include "types.h"
#include "io.h"
#include "fileio.h"
#include "ana_io.h"
#include "obs.h"
#include "mem.h"

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
  for (int s=1;s<=4;++s)
    S[i][j][s][l] = S_in[s];
}

void observable::set(fp_t **** S_in){
  memcpy(S[1][1][1]+1,S_in[1][1][1]+1,nx*ny*ns*nlambda*sizeof(fp_t));
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
  memcpy(S_copy[1]+1,S[i][j][i]+1,nlambda*ns*sizeof(fp_t)); 

  return S_copy;
}

fp_t * observable::get_lambda(){
  fp_t * lambda_copy;
  lambda = new fp_t [nlambda]-1;
  memcpy(lambda_copy+1,lambda+1,nlambda*sizeof(fp_t));
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
