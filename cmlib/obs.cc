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
  mask = new fp_t[nlambda]-1;
  memset(mask+1,0,nlambda*sizeof(fp_t));
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
}

int32_t observable::size(io_class &io_in){

  int32_t sz = 4*sizeof(int); // ns, nlambda,nx,ny
  sz += 2*nlambda*sizeof(fp_t); // lambda,mask
  sz += nx*ny*nlambda*ns*sizeof(fp_t); // actual observation
  return sz;
}

int32_t observable::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in){

  //printf("Packing : %d %d %d %d \n",nx,ny,ns,nlambda);
  int32_t offs=::pack(buf,nx,do_swap);
  offs+=::pack(buf+offs,ny,do_swap);
  offs+=::pack(buf+offs,ns,do_swap);
  offs+=::pack(buf+offs,nlambda,do_swap);
  
  offs+=::pack(buf+offs,lambda,1,nlambda,do_swap);
  offs+=::pack(buf+offs,mask,1,nlambda,do_swap);
  offs+=::pack(buf+offs,S,1,nx,1,ny,1,ns,1,nlambda,do_swap);

  return offs;
}

int32_t observable::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in){

  int32_t offs=::unpack(buf,nx,do_swap);
  offs+=::unpack(buf+offs,ny,do_swap);
  offs+=::unpack(buf+offs,ns,do_swap);
  offs+=::unpack(buf+offs,nlambda,do_swap);

  //printf("Unpacking : %d %d %d %d \n",nx,ny,ns,nlambda);

  lambda = new fp_t [nlambda]-1;
  mask = new fp_t [nlambda]-1;
  S=ft4dim(1,nx,1,ny,1,ns,1,nlambda);

  offs+=::unpack(buf+offs,lambda,1,nlambda,do_swap);
  offs+=::unpack(buf+offs,mask,1,nlambda,do_swap);
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
  /*for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      for (int s=1;s<=ns;++s)
        for (int l=1;l<=nlambda;++l)
          S[i][j][s][l] = S_in[l][s][j][i];*/
}

void observable::setlambda(fp_t * lambda_in){
  memcpy(lambda+1,lambda_in+1,nlambda*sizeof(fp_t));
}

void observable::setmask(fp_t * mask_in){
  memcpy(mask+1,mask_in+1,nlambda*sizeof(fp_t));
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
    //printf("%d %f %e \n",l,mask[l],S[i][j][1][l]);
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

int observable::get_n_lambda(){
  return nlambda;
}
int observable::get_n_lambda_to_fit(){
  int nl=0;
  for (int l=1;l<=nlambda;++l)
    if (mask[l]) nl+=1;
  return nl;
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
  // This is quite ad-hoc at the moment and we will need this as an input.

  fp_t qs = 3.275E14; // quiet sun reference continuum

  // Average continuum in all field, assume here for simplicity that point #30 is continuum:
  int l_ref = 30;
  fp_t mean = 0.0;
  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      mean += S[i][j][1][30] / nx / ny;

  mean = 28.69; // We did the above externally. Probably we neeed a way to do this better.
  fp_t scale = qs/mean;

  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      for (int l=1;l<=nlambda;++l){
        S[i][j][1][l] *= scale;
        for (int s=2;s<=4;++s)
          S[i][j][s][l] *= S[i][j][1][l];
      }
}

// ================================================================================================

void observable::correct_for_scattered_light(fp_t ratio){

  // Here we basically subtract the fixed ratio of scattered light from all the data. 
  // We DON't want to do this externally since maybe we will want to fit for
  // scattered light at some point. 
  // Also quite ad hoc

  fp_t qs = 3.275E14; // quiet sun reference continuum

  // This is performed AFTER the normalization:

  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j)
      for (int l=1;l<=nlambda;++l){
        S[i][j][1][l] -= qs*ratio;
        if (S[i][j][1][l] < 0) S[i][j][1][l] = 0.0; // Can't be negative!
  } // wvl, spatial 

}

// ================================================================================================

void observable::spectral_convolution(fp_t width){

  // Here we convolve the observation with a gaussian of given width. 
  // There has to be a function in all those libraries we are including which does this
  // For the moment it is implemented manually. Since lambda mash can be uneven, we first
  // interpolate (monothonically) the spectrum on a finer, evenly spaced grid. 
  // Convolve, in the wavelength space (i.e. no FFT or anything)
  // Interpolate back to the original grid.

  fp_t * lambda_fine; 
  int N_lambda_fine;

  // Guesstimate N_lambda_fine;
  fp_t l_min = lambda[1] - 5.0*width;
  fp_t l_max = lambda[nlambda] + 5.0*width;
  N_lambda_fine = int((l_max-l_min)/(0.2*width))+1;
  fp_t step = (l_min-l_max) / (N_lambda_fine-1);
  lambda_fine = new fp_t [N_lambda_fine]-1;
  for (int l=1;l<=N_lambda_fine;++l) 
    lambda_fine[l] = lambda_min + (l-1)*step;

  // Set-up a array for  keeping the kernel:
  fp_t * kernel = new fp_t [N_lambda_fine]-1;

  // Interpolate S to the existing grid:
  fp_t ** S_intepolated;
  S_intepolated = ft2dim(1,4,1,N_lambda_fine);
  for (int s=1;s<=4;++s)
    for (int l=1;l<=N_lambda_fine;++l)
      // This interpolation also extrapolates using constant values:
      S_intepolated[s][l] = interpol_1d(S[s]+1,lambda+1,nlambda,lambda_fine[l]);

  fp_t ** S_convolved;
  S_convolved = ft2dim(1,4,1,N_lambda_fine);

  fp_t * w_lambda = new fp_t [N_lambda_fine];
  w_lambda[1] = w_lambda[N_lambda_fine] = step*0.5; 
  for (int l=2;l<N_lambda_fine;++l)
    w_lambda[l] = step;

  // Then finally, convolve:
  for (int l=1;l<=N_lambda_fine;++l){ // for each lambda, convolve

    fp_t temp_to_integrate[N_lambda_fine]; 
    // Compute kernel:
    for (int ll=1;ll<=N_lambda_fine;++ll)
      kernel[ll] = 0.56418958354/width * exp(-(lambda[ll]-lambda[l])*(lambda[ll]-lambda[l])/width/width);
    // Multiply:
    for (int s=1;s<=4;++s){
      S_convolved[s][l] = 0.0;
      for (int ll=1;ll<=N_lambda_fine;++ll){
        temp_to_integrate[ll] = S_intepolated[s][ll] * kernel[ll];
        S_convolved += temp_to_integrate[ll] * w_lambda[ll];
      } // each lambda'
    } // each stokes component
  } // each lambda

  // Interpolate back:
  for (int s=1;s<=4;++s)
    for (int l=1;l<=nlambda;++l)
      S[s][l] = interpol_1d(S_convolved[s],lambda_fine,N_lambda_fine,lambda[l]);

  //cleanup:
  delete[](w_lambda+1);
  delete[](lambda_fine+1);
  delete[](kernel+1);
  del_ft2dim(S_convolved,1,4,1,N_lambda_fine);
  del_ft2dim(S_intepolated,1,4,1,N_lambda_fine);
}


// ================================================================================================

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

observable * obs_new(int nx,int ny,int ns,int nlambda){
  return new observable(nx,ny,ns,nlambda);
}

observable * obs_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in){
  return new observable(buf,offs,do_swap,io_in);
}