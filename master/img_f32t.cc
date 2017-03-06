#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <errno.h>

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
#include "img.h"
//#include "img_c64t.h"
#include "img_f32t.h"

image_f32t::image_f32t(int nx_in,int ny_in)
{
  pic=0;
  nx=nx_in;
  ny=ny_in;
  pic=f32t2dim(1,nx,1,ny);
  memset(pic[1]+1,0,nx*ny*sizeof(float32_t));
}

image_f32t::image_f32t(image_t *p,io_class &io)
{
  char *hh;
  p->info(nx,ny,hh);
  fp_t **tmp=ft2dim(1,nx,1,ny);
  memset(tmp[1]+1,0,nx*ny*sizeof(fp_t));
  p->add(tmp);
  pic=f32t2dim(1,nx,1,ny);
  for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) pic[x][y]=tmp[x][y];
  del_ft2dim(tmp,1,nx,1,ny);
}

image_f32t::image_f32t(const char *path_in,io_class &io)
{
  path=strcpy(new char [strlen(path_in)+1],path_in);
  int type;
  byte *data=read_image(path,nx,ny,type,header,io);
  if(data)
    if(type==MFBD_F32T){
      in_place_transpose((float32_t*)data,ny,nx);
      swap(nx,ny);
      pic=f32t2dim((float32_t*)data,1,nx,1,ny);
    }else{ // wrong type
      const char *names[]=MFBD_TYPE_NAMES;
      io.msg(IOL_ERROR,"image_f32t::image_f32t: input file \"%s\" has the wrong type (%s, expected %s).\n",path,names[type],names[MFBD_F32T]);
      delete[] data;
      nx=ny=0;
      if(header) delete[] header;
      header=0;
    }
}

image_f32t::image_f32t(const char *dir,const char *fnm,io_class &io)
{
  pic=0;
  struct stat stat_buf;
  char fname[1000];
  if(fnm[0]=='!') strcpy(fname,fnm+1); else strcpy(fname,fnm);
  if(dir){
    if(stat(dir,&stat_buf)<0) io.msg(IOL_ERROR,"data directory \"%s\" not found.\n",dir);  
    path=new char [strlen(dir)+strlen(fname)+2];
    if(dir[strlen(dir)-1]=='/')
      sprintf(path,"%s%s",dir,fname);
    else
      sprintf(path,"%s/%s",dir,fname);
  }else path=strcpy(new char [strlen(fname)+1],fname);
  if(stat(path,&stat_buf)<0){
    io.msg(IOL_ERROR,"input file \"%s\" not found.\n",path);
    pic=0;
    header=0;
    delete[] path;
    path=0;
    nx=ny=0;
  }else{
    int type;
    byte *data=read_image(path,nx,ny,type,header,io);
    if(data)
      if(type==MFBD_F32T){
        in_place_transpose((float32_t*)data,ny,nx);
        swap(nx,ny);
        pic=f32t2dim((float32_t*)data,1,nx,1,ny);
      }else{ // wrong type
        const char *names[]=MFBD_TYPE_NAMES;
        io.msg(IOL_ERROR,"input file \"%s\" has the wrong type (%s, expected %s).\n",path,names[type],names[MFBD_F32T]);
        delete[] data;
        nx=ny=0;
        if(header) delete[] header;
        header=0;
        if(path) delete[] path;
        path=0;
      }
  }
}


image_f32t::image_f32t(const char *dir,const char *fntmpl,int num,io_class &io)
{
  pic=0;
  struct stat stat_buf;
  char fname[1000];
  if(fntmpl[0]=='!') sprintf(fname,fntmpl+1,num); else sprintf(fname,fntmpl,num);
  if(dir){
    if(stat(dir,&stat_buf)<0) io.msg(IOL_ERROR,"data directory \"%s\" not found.\n",dir);  
    path=new char [strlen(dir)+strlen(fname)+2];
    if(dir[strlen(dir)-1]=='/')
      sprintf(path,"%s%s",dir,fname);
    else
      sprintf(path,"%s/%s",dir,fname);
  }else path=strcpy(new char [strlen(fname)+1],fname);
  if(stat(path,&stat_buf)<0){
    io.msg(IOL_ERROR,"input file \"%s\" not found.\n",path);  
    pic=0;
    header=0;
    delete[] path;
    path=0;
    nx=ny=0;
  }else{
    int type;
    byte *data=read_image(path,nx,ny,type,header,io);
    if(data)
      if(type==MFBD_F32T){
        in_place_transpose((float32_t*)data,ny,nx);
        swap(nx,ny);
        pic=f32t2dim((float32_t*)data,1,nx,1,ny);
      }else{ // wrong type
        const char *names[]=MFBD_TYPE_NAMES;
        io.msg(IOL_ERROR,"input file \"%s\" has the wrong type (%s, expected %s).\n",path,names[type],names[MFBD_F32T]);
        delete[] data;
        nx=ny=0;
        if(header) delete[] header;
        header=0;
        if(path) delete[] path;
        path=0;
      }
  }
}

int sum_extract(char *hh,const char *key)
{
  if(char *p=strstr(hh,key)){ // image is a sum
    --p;
    while((p>hh)&&(p[0]==' ')) --p; // remove whitespace
    while((p>hh)&&(p[0]!=' ')) --p; // seek beginning of number
    int n;
    if(sscanf(p,"%d",&n)) return n; else return 1;
  }
  return 1;
}

image_f32t::image_f32t(image_t **data,int n,io_class &io)
{
  pic=0;
  if(n){
    char *hh;
    data[1]->info(nx,ny,hh);
    if(hh) header=strcpy(new char [strlen(hh)+1],hh);
    fp_t **tmp=ft2dim(1,nx,1,ny);
    memset(tmp[1]+1,0,nx*ny*sizeof(fp_t));
    int N=0; // total # of images
    for(int i=1;i<=n;++i){
      int n_x,n_y;
      data[i]->info(n_x,n_y,hh);
      if((nx==n_x)&&(ny==n_y)){                // same dimensions?
        data[i]->add(tmp);                    // add image data to pic
        int nn=sum_extract(hh," ... SUM ... "); // image itself a sum?
        if(nn==1) io.msg(IOL_WARN,"dark field header did not contain SUM keyword! Assuming sum of 1\n");
        N+=nn;
      }else{
        io.msg(IOL_ERROR,"dimension mismatch: unable to add %dx%d image to %dx%d image\n",n_x,n_y,nx,ny);
        del_ft2dim(tmp,1,nx,1,ny);
        pic=0;
        nx=ny=0;
        if(header) delete[] header;
        header=0;
        return;
      }
    } // if we made it this far, we should have sum of N images: calculate average
    io.msg(IOL_XNFO,"calculating the average of %d images\n",N);
    pic=f32t2dim(1,nx,1,ny);
    for(int x=1;x<=nx;++x)
      for(int y=1;y<=ny;++y) pic[x][y]=tmp[x][y]/(fp_t)N;
    del_ft2dim(tmp,1,nx,1,ny);
  }
}

image_f32t::~image_f32t(void)
{
  if(pic) del_f32t2dim(pic,1,nx,1,ny);
}

void image_f32t::clip(io_class *io)
// (destructively) clip the image by shifting data and updating pointers
{
  float32_t *xs=new float32_t [nx] -1,*ys=new float32_t [ny] -1; 
  memset(xs+1,0,nx*sizeof(float32_t));
  memset(ys+1,0,ny*sizeof(float32_t));
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y){
      xs[x]+=pic[x][y];
      ys[y]+=pic[x][y];
    }
  int xl,xh,yl,yh;
  for(xl=1;xs[xl]==0;++xl);
  for(xh=nx;xs[xh]==0;--xh);
  for(yl=1;ys[yl]==0;++yl);
  for(yh=ny;ys[yh]==0;--yh);
  delete[] (xs+1);
  delete[] (ys+1);
  nx=xh-xl+1;
  ny=yh-yl+1;
  for(int x=xl;x<=xh;++x) memmove(pic[1]+(x-xl)*ny+1,pic[x]+yl,ny*sizeof(float32_t));
  for(int x=2;x<=nx;++x) pic[x]=pic[x-1]+ny; // update pointers
}

fp_t **image_f32t::add(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]+=pic[x][y];
  return data;
}

fp_t **image_f32t::sub(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]-=pic[x][y];
  return data;
}

fp_t **image_f32t::safediv(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) if(pic[x][y]!=0.0) data[x][y]/=pic[x][y]; else data[x][y]=0.0;
  return data;
}

fp_t **image_f32t::mult(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]*=pic[x][y];
  return data;
}

static __inline fp_t sf(fp_t x) // x=0..1
{
  return 0.5*(1.0-cos(pi*x));
}

void image_f32t::paste(fp_t **sub_img,int x,int y,int *xl,int *xh,int *yl,int *yh,int mx,int my,fp_t angle,io_class &io)
{
  int snx=xh[x]-xl[x]+1,sny=yh[y]-yl[y]+1;
  fp_t **sim=sub_img;
  if(angle) sim=rotate(sub_img,1,snx,1,sny,-angle);
  int margin=8;
  int offx=(xh[x]-xl[x]+1)/margin,offy=(yh[y]-yl[y]+1)/margin;
  int exl=xl[x]+offx,exh=xh[x]-offx,eyl=yl[y]+offy,eyh=yh[y]-offy;
  if((exl<1)||(eyl<1)||(exh>nx)||(eyh>ny)){
    io.msg(IOL_WARN,"warning: cannot paste subimage (%d,%d,%d,%d) into image (1,%d,1,%d). clipping...\n",xl[x],xh[x],yl[y],yh[y],nx,ny);
    if(exl<1) exl=1;
    if(eyl<1) eyl=1;
    if(exh>nx) exh=nx;
    if(eyh>ny) eyh=ny;
  }
  int ixl=(x> 1)?xh[x-1]-(xh[x-1]-xl[x-1]+1)/margin+1:exl;
  int ixh=(x<mx)?xl[x+1]+(xh[x+1]-xl[x+1]+1)/margin-1:exh;
  int iyl=(y> 1)?yh[y-1]-(yh[y-1]-yl[y-1]+1)/margin+1:eyl;
  int iyh=(y<my)?yl[y+1]+(yh[y+1]-yl[y+1]+1)/margin-1:eyh;
//
  fp_t **window=ft2dim(exl,exh,eyl,eyh);
  for(int ix=exl;ix<=exh;++ix)
    for(int iy=eyl;iy<=eyh;++iy) window[ix][iy]=1.0;
  for(int ix=exl;ix<=ixl-1;++ix){ // lower x margin
    fp_t xx=(fp_t)(ix-exl+1)/(fp_t)(ixl-exl+1);
    for(int iy=eyl;iy<=eyh;++iy) window[ix][iy]*=sf(xx)/(sf(xx)+sf(1.0-xx));
  }
  for(int ix=ixh+1;ix<=exh;++ix){ // upper x margin
    fp_t xx=1.0-(fp_t)(ix-ixh)/(fp_t)(exh-ixh+1);
    for(int iy=eyl;iy<=eyh;++iy) window[ix][iy]*=sf(xx)/(sf(xx)+sf(1.0-xx));
  }
  for(int iy=eyl;iy<=iyl-1;++iy){ // lower y margin
    fp_t yy=(fp_t)(iy-eyl+1)/(fp_t)(iyl-eyl+1);
    for(int ix=exl;ix<=exh;++ix) window[ix][iy]*=sf(yy)/(sf(yy)+sf(1.0-yy));
  }
  for(int iy=iyh+1;iy<=eyh;++iy){ // upper y margin
    fp_t yy=1.0-(fp_t)(iy-iyh)/(fp_t)(eyh-iyh+1);
    for(int ix=exl;ix<=exh;++ix) window[ix][iy]*=sf(yy)/(sf(yy)+sf(1.0-yy));
  }
  for(int ix=exl;ix<=exh;++ix)
    for(int iy=eyl;iy<=eyh;++iy)
      pic[ix][iy]+=(float32_t)(window[ix][iy]*sim[ix-xl[x]+1][iy-yl[y]+1]);
  del_ft2dim(window,exl,exh,eyl,eyh);
  if(angle) del_ft2dim(sim,1,snx,1,sny);
}

#include <fitsio.h>
#include <unistd.h>
#include "jnfo.h"

void image_f32t::write(char *fname,struct jnfo &ji,int o,io_class *io)
{
}

fp_t image_f32t::mean(int offs)
{
  fp_t sum=0.0;
  for(int x=offs;x<=nx-offs;++x)
    for(int y=offs;y<=ny-offs;++y) sum+=pic[x][y];
  return sum/(fp_t)((nx-2*offs)*(ny-2*offs));
}

void image_f32t::mult(fp_t f)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) pic[x][y]=f*pic[x][y];
}

fp_t **image_f32t::subimage(int xl,int xh,int yl,int yh)
{
  int snx=xh-xl+1,sny=yh-yl+1;
  fp_t **data=ft2dim(1,snx,1,sny);
  memset(data[1]+1,0,snx*sny*sizeof(fp_t)); // pad with 0
  int ixl=max(xl,1),ixh=min(xh,nx),iyl=max(yl,1),iyh=min(yh,ny);
  for(int x=ixl;x<=ixh;++x) for(int y=iyl;y<=iyh;++y) data[x-xl+1][y-yl+1]=pic[x][y];
  return data;
}

void image_f32t::clip(int16_t &xl,int16_t &xh,int16_t &yl,int16_t &yh,io_class *io)
// (destructively) clip the image by shifting data and updating pointers
{
  if((!xl)&&(!xh)&&(!yl)&&(!yh)) return; // don't clip 
  if((xl==1)&&(xh==nx)&&(yl==1)&&(yh==ny)) return; // don't clip
  int flip_x=0,flip_y=0;
  if(xl>xh){
    swap(xl,xh);
    flip_x=1;
  }
  if(yl>yh){
    swap(yl,yh);
    flip_y=1;
  }
  io->msg(IOL_XNFO,"clipping %s to %d,%d,%d,%d, x-flip=%s, y-flip=%s\n",path,xl,xh,yl,yh,(flip_x)?"yes":"no",(flip_y)?"yes":"no");
  if((xl<1)||(yl<1)||(xh>nx)||(yh>ny)){
    io->msg(IOL_ERROR,"cannot extend image (1,1,%d,%d) to (%d,%d,%d,%d)\n",nx,ny,xl,xh,yl,yh);
    xl=1;
    yl=1;
    xh=nx;
    yh=ny;
    return;
  }
  nx=xh-xl+1;
  ny=yh-yl+1;
  for(int x=xl;x<=xh;++x) memmove(pic[1]+(x-xl)*ny+1,pic[x]+yl,ny*sizeof(float32_t));
  for(int x=2;x<=nx;++x) pic[x]=pic[x-1]+ny; // update pointers
  if(flip_x){
    for(int x=1;x<=nx/2;++x)
      for(int y=1;y<=ny;++y) swap(pic[x][y],pic[nx-x+1][y]);
  }
  if(flip_y){
    for(int x=1;x<=nx;++x)
      for(int y=1;y<=ny/2;++y) swap(pic[x][y],pic[x][ny-y+1]);  
  }
}

fp_t **image_f32t::get(void)
{
  fp_t **data=ft2dim(1,nx,1,ny);
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]=(fp_t)pic[x][y];
  return data;
}

void image_f32t::put(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) pic[x][y]=(float32_t)data[x][y];
}

fp_t image_f32t::noise_sigma(int offs)
{
  int Nx=nx-2*offs,Ny=ny-2*offs;
  for(int n=0;n<=30;++n){
    int mask=~(1<<n);
    if(Nx&mask) Nx=Nx&mask;
    if(Ny&mask) Ny=Ny&mask;
  }
  int xl=((nx-Nx)/2)+1,yl=((ny-Ny)/2)+1;
  int xh=xl+Nx-1,yh=yl+Ny-1;
    
  complex_t **ft=ct2dim(1,Nx,1,Ny);
  memset(ft[1]+1,0,Nx*Ny*sizeof(complex_t));
  fp_t sum=0.0;
  for(int x=xl;x<=xh;++x)
    for(int y=yl;y<=yh;++y) sum+=ft[x-xl+1][y-yl+1].re=pic[x][y];
  sum/=(fp_t)(Nx*Ny);
  fp_t **window=make_window(Nx,Ny,8);
  for(int x=1;x<=Nx;++x)
    for(int y=1;y<=Ny;++y) ft[x][y].re=(ft[x][y].re-sum)*window[x][y];
  del_ft2dim(window,1,Nx,1,Ny);
//
  fft_2d(ft,Nx,Ny,1);
  fft_reorder_2d(ft,Nx,Ny);
//
  fp_t noise_power=0;
  int N=0,xx,yy;
  fp_t fxy=(fp_t)Nx/(fp_t)Ny;
  for(int x=1;x<=Nx;++x)
    if(abs(xx=x-Nx/2-1)>(Nx/6))
      for(int y=1;y<=Ny;++y)
        if(abs(yy=y-Ny/2-1)>(Ny/6))
          if(((fp_t)sqr(xx)+fxy*(fp_t)sqr(yy))>((fp_t)sqr(Nx)/4.0)){
            noise_power+=ft[x][y].re*ft[x][y].re+ft[x][y].im*ft[x][y].im;
            ++N;
          }
//
  del_ct2dim(ft,1,Nx,1,Ny);
  return sqrt(noise_power/((fp_t)N*(fp_t)Nx*(fp_t)Ny));
}

image_t *image_f32t::clone(io_class &io)
{
  return new image_f32t(this,io);
}
