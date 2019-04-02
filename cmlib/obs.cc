#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include "types.h"
#include "io.h"
#include "fileio.h"
#include "ana_io.h"
#include "obs.h"
#include "mem.h"
#include "pack.h"
#include "mathtools.h"

observable::observable(int ns_in):ns(ns_in)
{
  S=0;
  nlambda=0;
  nx=0;ny=0;
  w_stokes = new fp_t[4];
}

observable::observable(int nx_in,int ny_in,int ns_in)
{
  S=0;
  nlambda=0;
  ns=ns_in;
  nx=nx_in;ny=ny_in;
  w_stokes = new fp_t[4];
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
  mask = new fp_t[nlambda]-1;
  memset(mask+1,0,nlambda*sizeof(fp_t));
  scattered_light=0.0;
  spectral_broadening=0.0;
  synth_qs = 1.0; obs_qs = 1.0;
  el=0.0;az=0.0;
  to_invert=0;no_iterations=0;start_lambda=1.0;
  w_stokes = new fp_t [4];
  memset(w_stokes,0,4*sizeof(fp_t));
}

observable::observable(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in){
  offs+=unpack(buf+offs,do_swap,io_in);
}


observable::~observable(void){
  if(nlambda){
    del_ft4dim(S,1,nx,1,ny,1,ns,1,nlambda);
    delete[] (lambda+1);
    delete[] (mask+1);
  }
  delete []w_stokes;
}

int32_t observable::size(io_class &io_in){

  int32_t sz = 4*sizeof(int); // ns, nlambda,nx,ny
  sz += 7*sizeof(fp_t); // scattered light,broadening,continuum,el,az, starting lambda for lm
  sz += 2*sizeof(int); // whether to invert or no, max number of iterations
  sz += 2*nlambda*sizeof(fp_t); // lambda,mask
  sz += nx*ny*nlambda*ns*sizeof(fp_t); // actual observation
  sz += 4*sizeof(fp_t); // weights for stokes;
  return sz;
}

int32_t observable::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in){

  int32_t offs=::pack(buf,nx,do_swap);
  offs+=::pack(buf+offs,ny,do_swap);
  offs+=::pack(buf+offs,ns,do_swap);
  offs+=::pack(buf+offs,nlambda,do_swap);
  offs+=::pack(buf+offs,scattered_light,do_swap);
  offs+=::pack(buf+offs,spectral_broadening,do_swap);
  offs+=::pack(buf+offs,obs_qs,do_swap);
  offs+=::pack(buf+offs,synth_qs,do_swap);
  offs+=::pack(buf+offs,el,do_swap);
  offs+=::pack(buf+offs,az,do_swap);
  offs+=::pack(buf+offs,to_invert,do_swap);
  offs+=::pack(buf+offs,no_iterations,do_swap);
  offs+=::pack(buf+offs,start_lambda,do_swap);
  offs+=::pack(buf+offs,lambda,1,nlambda,do_swap);
  offs+=::pack(buf+offs,mask,1,nlambda,do_swap);
  offs+=::pack(buf+offs,w_stokes,0,3,do_swap);
  offs+=::pack(buf+offs,S,1,nx,1,ny,1,ns,1,nlambda,do_swap);

  return offs;
}

int32_t observable::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in){

  int32_t offs=::unpack(buf,nx,do_swap);
  offs+=::unpack(buf+offs,ny,do_swap);
  offs+=::unpack(buf+offs,ns,do_swap);
  offs+=::unpack(buf+offs,nlambda,do_swap);
  offs+=::unpack(buf+offs,scattered_light,do_swap);
  offs+=::unpack(buf+offs,spectral_broadening,do_swap);
  offs+=::unpack(buf+offs,obs_qs,do_swap);
  offs+=::unpack(buf+offs,synth_qs,do_swap);
  offs+=::unpack(buf+offs,el,do_swap);
  offs+=::unpack(buf+offs,az,do_swap);
  offs+=::unpack(buf+offs,to_invert,do_swap);
  offs+=::unpack(buf+offs,no_iterations,do_swap);
  offs+=::unpack(buf+offs,start_lambda,do_swap);

  lambda = new fp_t [nlambda]-1;
  mask = new fp_t [nlambda]-1;
  S=ft4dim(1,nx,1,ny,1,ns,1,nlambda);
  w_stokes = new fp_t [4];

  offs+=::unpack(buf+offs,lambda,1,nlambda,do_swap);
  offs+=::unpack(buf+offs,mask,1,nlambda,do_swap);
  offs+=::unpack(buf+offs,w_stokes,0,3,do_swap);
  offs+=::unpack(buf+offs,S,1,nx,1,ny,1,ns,1,nlambda,do_swap);

  return offs;

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
}

void observable::setlambda(fp_t * lambda_in){
  memcpy(lambda+1,lambda_in+1,nlambda*sizeof(fp_t));
}

void observable::setmask(fp_t * mask_in){
  memcpy(mask+1,mask_in+1,nlambda*sizeof(fp_t));
}

void observable::set_inv_parameters(fp_t sl_in, fp_t sb_in, fp_t obs_qs_in, fp_t synth_qs_in){
  scattered_light = sl_in;
  spectral_broadening = sb_in;
  obs_qs = obs_qs_in;
  synth_qs = synth_qs_in;
}

void observable::set_viewing_angle(fp_t el_in, fp_t az_in){
  el=el_in;
  az=az_in;
}

void observable::set_to_invert(int to_invert_in){
  to_invert=to_invert_in;
}

void observable::set_no_iterations(int input){
  no_iterations = input;
}

void observable::set_start_lambda(fp_t input){
  start_lambda = input;
}

void observable::set_w_stokes(fp_t * w_stokes_input){
  memcpy(w_stokes,w_stokes_input,4*sizeof(fp_t));
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

fp_t ** observable::get_S_to_fit(int i, int j){

  fp_t ** S_copy;
  int nl_to_fit = get_n_lambda_to_fit();
  S_copy = ft2dim(1,ns,1,nl_to_fit);
  int lf=1;
  for (int l=1;l<=nlambda;++l){
    if (mask[l]){
      for (int s=1;s<=4;++s) S_copy[s][lf] = S[i][j][s][l];
      ++lf;
    }
  }

  return S_copy;
}


fp_t * observable::get_lambda(){
  fp_t * lambda_copy;
  lambda_copy = new fp_t [nlambda]-1;
  memcpy(lambda_copy+1,lambda+1,nlambda*sizeof(fp_t));
  return lambda_copy;
}

fp_t * observable::get_lambda_to_fit(){
  fp_t * lambda_copy;
  int nl_to_fit = get_n_lambda_to_fit();
  lambda_copy = new fp_t [nl_to_fit]-1;
  int lf=1;
  for (int l=1;l<=nlambda;++l)
    if (mask[l]){
      lambda_copy[lf] = lambda[l];
      ++lf;
    }
  return lambda_copy;
}

fp_t * observable::get_mask(){
  fp_t * mask_copy;
  mask_copy = new fp_t [nlambda]-1;
  memcpy(mask_copy+1,mask+1,nlambda*sizeof(fp_t));
  return mask_copy;
}

fp_t * observable::get_w_stokes(){
  fp_t * w_stokes_copy;
  w_stokes_copy = new fp_t[4];
  memcpy(w_stokes_copy,w_stokes,4*sizeof(fp_t));
  return w_stokes_copy;
}

int observable::get_n_lambda(){
  return nlambda;
}
int observable::get_n_lambda_to_fit(){
  int nl=0;
  for (int l=1;l<=nlambda;++l)
    if (mask[l]) nl+=1;
  return nl;
}

fp_t observable::get_scattered_light(){
  return scattered_light;
}

fp_t observable::get_spectral_broadening(){
  return spectral_broadening;
}

fp_t observable::get_synth_qs(){
  return synth_qs;
}

fp_t observable::get_el(){
  return el;
}

fp_t observable::get_az(){
  return az;
}
int observable::get_to_invert(){
  return to_invert;
}
int observable::get_no_iterations(){
  return no_iterations;
}
fp_t observable::get_start_lambda(){
  return start_lambda;
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
  obs_small->set_inv_parameters(scattered_light,spectral_broadening,obs_qs,synth_qs);

  fp_t * mask_small = new fp_t [nln]-1;
  for (int l=ll;l<=lh;++l)
    mask_small[l-ll+1] = mask[l];
  obs_small->setmask(mask_small);
  delete[](mask_small+1);
  return obs_small;

}

// ================================================================================================

void observable::normalize(){
  // Normalizes already arranged observable to physical units. 
  
  fp_t scale = synth_qs/obs_qs;

  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      for (int l=1;l<=nlambda;++l)
        for (int s=1;s<=4;++s)
        S[i][j][s][l] *= scale;
}

// ================================================================================================

void observable::add_scattered_light(fp_t fraction, fp_t continuum_level){

  // Here we basically subtract the fixed ratio of scattered light from all the data. 
  // We DON't want to do this externally since maybe we will want to fit for
  // scattered light at some point. 
  // Also quite ad hoc

  // This is performed AFTER the normalization:

  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      for (int l=1;l<=nlambda;++l){
        // THis looks wierd to me:
        //S[i][j][1][l] = S[i][j][1][l]/(1.0+fraction) + fraction/(1.0+fraction) * continuum_level;
        S[i][j][1][l] = S[i][j][1][l]+ fraction * continuum_level;
  } // wvl, spatial 
}

// ================================================================================================

void observable::spectral_convolve(fp_t width, int i, int j){

  convolve_spectra_with_gauss(S[i][j],lambda,nlambda,width);
}


// ================================================================================================

void observable::read(char * name, io_class &io){
}

observable * obs_new(int nx,int ny,int ns,int nlambda){
  return new observable(nx,ny,ns,nlambda);
}

observable * obs_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in){
  return new observable(buf,offs,do_swap,io_in);
}