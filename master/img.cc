#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include "fftw3.h"

#include "types.h"
#include "fft.h"
#include "fileio.h"
#include "mem.h"
#include "io.h"
#include "uts.h"
#include "conf.h"
#include "const.h"
#include "ana_io.h"
#include "imgtools.h"
#include "fillpix.h"

#include "img_i16t.h"
#include "img_i32t.h"
#include "img_f32t.h"
#include "img_c64t.h"

#include "img.h"

char *imagepath(const char *dir,const char *fntmpl,int num,io_class &io)
{
  char fname[1000],*path=0;
  if(fntmpl[0]=='!') sprintf(fname,fntmpl+1,num); else sprintf(fname,fntmpl,num);
  if(dir){
    struct stat stat_buf;
    if(stat(dir,&stat_buf)<0) io.msg(IOL_ERROR,"data directory \"%s\" not found.\n",dir);  
    path=new char [strlen(dir)+strlen(fname)+2];
    if(dir[strlen(dir)-1]=='/')
      sprintf(path,"%s%s",dir,fname);
    else
      sprintf(path,"%s/%s",dir,fname);
  }else path=strcpy(new char [strlen(fname)+1],fname);
  return path;
}

image_t *readimage(const char *dir,const char *fntmpl,int num,io_class &io)
{
  char *path=imagepath(dir,fntmpl,num,io);
  struct stat stat_buf;
  if(stat(path,&stat_buf)<0){
    io.msg(IOL_ERROR,"input file \"%s\" not found.\n",path);
    return 0;
  }
  const char *names[]=MFBD_TYPE_NAMES;
  switch(int t=data_type(path,io)){
    case(MFBD_I16T): return new image_i16t(path,io);
    case(MFBD_I32T): return new image_i32t(path,io);
    case(MFBD_F32T): return new image_f32t(path,io);
    case(MFBD_I08T):
    case(MFBD_I64T):
    case(MFBD_F64T):{
      io.msg(IOL_ERROR,"readimage: input file \"%s\" has unsupported type %s.\n",path,names[t]);
      break;
    }
    default: io.msg(IOL_ERROR,"readimage: input file \"%s\" has unknown type.\n",path);
  }
  return 0;
}

image_t::image_t(void){
  nx=0;
  ny=0;
  header=0;
  path=0;
}

image_t::~image_t(void){ 
  if(header) delete[] header;
  if(path) delete[] path;
}
  
int64_t extract_time(char *header)
{
  char *p=strchr(header,':');
  if(!p) return -1;
  int hh,m,s,ms;
  sscanf(p-2,"%d:%d:%d.%d",&hh,&m,&s,&ms);
  return (int64_t)(3600000*hh+60000*m+1000*s+ms);
}

void image_preproc(int &nx,int &ny,image_t ****image,fp_t ***noise_sigma,fp_t **nf,int no,int *nc,int **nim,char *&time_obs,int border_clip,io_class &io)
// Normalize image ratio with linear fit to compensate for accelerating shutter.
{
  char *header=0;
  int64_t atime=0;
  int nt=0;
  image[1][1][1]->info(nx,ny,header);
  if(header)
    if(strcmp(header,"")&&(atime>=0)){
      atime+=extract_time(header);
      ++nt;
    }
  if(atime<=0) atime=-1;
  for(int o=1;o<=no;++o)
    for(int l=1;l<=nc[o];++l)
      for(int i=1;i<=nim[o][l];++i){
        int tnx,tny;
        image[o][l][i]->info(tnx,tny,header);
        if(header){
          if(strcmp(header,"")&&(atime>=0)){
            atime+=extract_time(header);
            ++nt;
          }else atime=-1;
        }else atime=-1;
        if((tnx!=nx)||(tny!=ny))
          io.msg(IOL_WARN,"img size mismatch: image (1,1,1) = %dx%d but image (%d,%d,%d) = %dx%d\n",nx,ny,o,k,i,tnx,tny);
      }
  io.msg(IOL_INFO,"raw image size = %dx%d\n",nx,ny);
  if(!strlen(time_obs)){
    if(atime>=0){
      atime/=(int64_t)nt;
      int hh=(int)(atime/(int64_t)3600000);
      int mm=(int)(atime-(int64_t)(3600000*hh))/(int64_t)60000;
      int ss=(int)(atime-(int64_t)(3600000*hh)-(int64_t)(60000*mm))/(int64_t)1000;
      int ms=(int)(atime-(int64_t)(3600000*hh)-(int64_t)(60000*mm)-(int64_t)(1000*ss));
      delete[] time_obs;
      time_obs=new char [13];
      sprintf(time_obs,"%02d:%02d:%02d.%03d",hh,mm,ss,ms);
      io.msg(IOL_INFO,"average time=%s\n",time_obs);
    }else io.msg(IOL_WARN,"average time not available from input data.\n");
  }
// normalisation (warning! allowed within each object only!)
  for(int o=1;o<=no;++o){
    fp_t **means=ft2dim(1,nc[o],1,nim[o]);
    fp_t maxmean=image[o][1][1]->mean(border_clip);
    for(int l=1;l<=nc[o];++l)
      for(int i=1;i<=nim[o][l];++i)
        maxmean=max(means[l][i]=image[o][l][i]->mean(border_clip),maxmean);
//
    for(int l=1;l<=nc[o];++l)
      for(int i=1;i<=nim[o][l];++i)
        if(means[l][i])
          image[o][l][i]->mult(maxmean/means[l][i]);
        else
          io.msg(IOL_WARN,"cannot normalize: mean[%d][%d][%d] is 0. (is border_clip too large)?\n",o,l,i);
    del_ft2dim(means,1,nc[o],1,nim[o]);
  }
// noise estimate
  fft_init(nx-2*border_clip,ny-2*border_clip);
  for(int o=1;o<=no;++o)
    for(int l=1;l<=nc[o];++l)
      if(nf[o][l])
        for(int i=1;i<=nim[o][l];++i) noise_sigma[o][l][i]=nf[o][l]*image[o][l][i]->noise_sigma(border_clip);
      else
        memset(noise_sigma[o][l]+1,0,nim[o][l]*sizeof(fp_t));
  fft_done(nx-2*border_clip,ny-2*border_clip);
//
  int nhx=nx/2;
  int nhy=ny/2;
  int size=512;
  int sh=size/2;
}

void image_t::descatter(fp_t **data,image_t *psf,image_t *bgain,io_class &io)
{
  image_c64t *p=(image_c64t*)psf;
  fp_t **img=ft2dim(1,2*nx,1,2*ny);
  memset(img[1]+1,0,4*nx*ny*sizeof(fp_t));
  fp_t **im=new fp_t* [nx]-1;
  for(int x=1;x<=nx;++x){
    im[x]=img[x+nx/2]+ny/2;
    for(int y=1;y<=ny;++y) im[x][y]=1.0;
  }
  im=bgain->mult(im); // g
  im=bgain->mult(im); // g*g
  for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) im[x][y]=data[x][y]/(1.0+im[x][y]);
//
  fp_t **cimg=ft2dim(1,2*nx,1,2*ny),cs;
  fp_t **cim=new fp_t* [nx]-1;
  for(int x=1;x<=nx;++x) cim[x]=cimg[x+nx/2]+ny/2;
//
  int i=0,iter_max=50;
  do{
    memcpy(cimg[1]+1,img[1]+1,4*nx*ny*sizeof(fp_t));
    cim=bgain->mult(cim);
    cimg=p->convolve(cimg);
    cim=bgain->mult(cim);
    cs=0.0;
    for(int x=1;x<=nx;++x)
      for(int y=1;y<=ny;++y){
        fp_t newimg=data[x][y]-cim[x][y];
        cs+=sqr(im[x][y]-newimg);
        im[x][y]=newimg;
      }
    cs/=(fp_t)(nx*ny);
    io.msg(IOL_XNFO,"image_t::descatter iteration %d: cs=%E\n",i,cs);
  }while((cs>1E-8)&&(++i<iter_max));
  delete[] (cim+1);
  del_ft2dim(cimg,1,2*nx,1,2*ny);
  for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) data[x][y]=im[x][y];
  delete[] (im+1);
  del_ft2dim(img,1,2*nx,1,2*ny);
}

void image_t::fft_reorder(fp_t **f,int Nx,int Ny)
{
  int nh=Nx/2;
  double *buf=new double [nh];
  for(int x=1;x<=nh;++x){
    memcpy(buf,f[x]+1,nh*sizeof(double));
    memcpy(f[x]+1,f[x+nh]+nh+1,nh*sizeof(double));
    memcpy(f[x+nh]+nh+1,buf,nh*sizeof(double));
//
    memcpy(buf,f[x+nh]+1,nh*sizeof(double));
    memcpy(f[x+nh]+1,f[x]+nh+1,nh*sizeof(double));
    memcpy(f[x]+nh+1,buf,nh*sizeof(double));
  }
  delete[] buf;
}

image_t *image_t::FT(io_class &io)
{
  int Nx=2*nx,Ny=2*ny;
  fp_t **img=ft2dim(1,Nx,1,Ny);
  memset(img[1]+1,0,Nx*Ny*sizeof(fp_t));
  fp_t **im=new fp_t* [nx]-1;
  for(int x=1;x<=nx;++x) im[x]=img[x+nx/2]+ny/2;
  this->add(im);
  delete[] (im+1);
// normalize
  fp_t integral=0.0;
  for(int x=1;x<=Nx;++x) for(int y=1;y<=Ny;++y) integral+=img[x][y];
  for(int x=1;x<=Nx;++x) for(int y=1;y<=Ny;++y) img[x][y]/=integral;
//
  fft_reorder(img,Nx,Ny);
//
  complex_t **ft=ct2dim(1,Nx,1,Ny/2+1);  
  fftw_plan fplan=fftw_plan_dft_r2c_2d(Nx,Ny,img[1]+1,(fftw_complex*)(ft[1]+1),FFTW_FORWARD|FFTW_ESTIMATE);
  fftw_plan bplan=fftw_plan_dft_c2r_2d(Nx,Ny,(fftw_complex*)(ft[1]+1),img[1]+1,FFTW_BACKWARD|FFTW_ESTIMATE);
  fftw_execute_dft_r2c(fplan,img[1]+1,(fftw_complex*)(ft[1]+1));
  del_ft2dim(img,1,Nx,1,Ny);
  image_t *i=new image_c64t(ft,Nx,Ny,fplan,bplan,io);
  del_ct2dim(ft,1,Nx,1,Ny/2+1);
  return i;
}

void image_t::flatfield(image_t *dark,image_t *resp,image_t *gain,image_t *psf,image_t *bgain,byte method,io_class &io)
{
  char *hh;
  int nxd,nyd,nzd;
  dark->info(nxd,nyd,hh);
  if((nx!=nxd)||(ny!=nyd)){
    io.msg(IOL_WARN,"image_t::flatfield: not flatfielding %s: dark field dimensions (%dx%d) are different from image dimensions (%dx%d)\n",path,nxd,nyd,nx,ny);
    return;
  }
  if(resp){
    resp->info(nxd,nyd,hh);
    if((nx!=nxd)||(ny!=nyd))
      io.msg(IOL_WARN,"image_t::flatfield: not correcting detector response for %s: response file dimensions (%dx%d) are different from image dimensions (%dx%d)\n",path,nxd,nyd,nx,ny);
  }
  gain->info(nxd,nyd,hh);
  if((nx!=nxd)||(ny!=nyd)){
    io.msg(IOL_WARN,"image_t::flatfield: not flatfielding %s: gain file dimensions (%dx%d) are different from image dimensions (%dx%d)\n",path,nxd,nyd,nx,ny);
    return;
  }
  io.msg(IOL_INFO,"image_t::flatfield: darkfielding %s\n",path);
// subtract dark current
  fp_t **data=get();
  data=dark->sub(data);
//
  io.msg(IOL_INFO,"image_t::flatfield: response correction for %s\n",path);
  if(!(data=resp->mult(data))){ // correct for the detector response (this should not contain the gain correction and must be done before descattering)
    io.msg(IOL_ERROR,"image_t::flatfield: response correction returned no data (bad CCD response file?)!\n");
    return;
  }
  if(bgain&&psf){          // apply backscatter correction
    io.msg(IOL_INFO,"image_t::flatfield: CCD transparency correction %s\n",path);
    descatter(data,psf,bgain,io);
  }
  io.msg(IOL_INFO,"image_t::flatfield: flatfielding %s\n",path);
  if(data=gain->mult(data)){                       // apply gain correction
    if(!(data=fillpix(data,1,nx,1,ny,method,io))){ // fill bad pixels
      del_ft2dim(data,1,nx,1,ny);
      return;
    }
  }else{
    io.msg(IOL_ERROR,"flatfielding returned no data (bad gain file?)!\n");
    return;
  }
//
  put(data);
  del_ft2dim(data,1,nx,1,ny);
}

