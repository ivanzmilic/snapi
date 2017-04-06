#include <string.h>
#include <errno.h>
#include "types.h"
#include "io.h"
#include "obs.h"
#include "mem.h"

observable::observable(int ns_in):ns(ns_in)
{
  S=0;
  nlambda=0;
}

void observable::readsingle(const char* input_file){
  S = 0;
  nlambda =0;
  float l_in;
  float I_in, Q_in, U_in, V_in;
  fp_t * I_tmp = new fp_t [4]-1;
  FILE * input;
  input = fopen(input_file,"r");
  while(!feof(input)){
    //fscanf(input,"%e %e %e %e %e", &l_in, &I_in, &Q_in, &U_in, &V_in);
    I_tmp[1] = I_in; I_tmp[2] = Q_in; I_tmp[3] = U_in; I_tmp[4] = V_in;
    add(I_tmp, fp_t(l_in));
  }
  printf("I read spectra %d lines long.\n", nlambda);
}

observable::~observable(void){
  if(nlambda){
    del_ft2dim(S,1,ns,1,nlambda);
    delete[] (lambda+1);
  }
}

void observable::add(fp_t *S_in,fp_t lambda_in)
{
  fp_t **tmp=ft2dim(1,ns,1,nlambda+1);
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
  lambda[nlambda]=lambda_in;
}

fp_t ** observable::get_S(){

  fp_t ** S_copy;
  S_copy = ft2dim(1,ns,1,nlambda);
  memcpy(S_copy[1]+1,S[1]+1,ns*nlambda*sizeof(fp_t));

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

void observable::write(const char *name,io_class &io)
{
  if(FILE *f=fopen(name,"w")){
    for(int l=1;l<=nlambda;++l){
      fprintf(f, "%.15e ", lambda[l]);
      for(int s=1;s<=ns;++s) fprintf(f, "%.15e ",S[s][l]);
      fprintf(f, "\n");
    }
    fclose(f);
  }else io.msg(IOL_WARN,"failed to open file \"%s\":\"%s\"\n",name,strerror(errno));
}
