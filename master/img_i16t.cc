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

#include "img_i16t.h"

image_i16t::image_i16t(int nx_in,int ny_in)
{
  pic=0;
  nx=nx_in;
  ny=ny_in;
  pic=i16t2dim(1,nx,1,ny);
  memset(pic[1]+1,0,nx*ny*sizeof(int16_t));
}

image_i16t::image_i16t(image_t *p,io_class &io)
{
  char *hh;
  p->info(nx,ny,hh);
  fp_t **tmp=ft2dim(1,nx,1,ny);
  memset(tmp[1]+1,0,nx*ny*sizeof(fp_t));
  p->add(tmp);
  pic=i16t2dim(1,nx,1,ny);
  for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) pic[x][y]=(int16_t)tmp[x][y];
  del_ft2dim(tmp,1,nx,1,ny);
}

image_i16t::image_i16t(const char *path_in,io_class &io)
{
  path=strcpy(new char [strlen(path_in)+1],path_in);
  int type;
  byte *data=read_image(path,nx,ny,type,header,io);
  if(data)
    if(type==MFBD_I16T){
      in_place_transpose((int16_t*)data,ny,nx);
      swap(nx,ny);
      pic=i16t2dim((int16_t*)data,1,nx,1,ny);
    }else{ // wrong type
      const char *names[]=MFBD_TYPE_NAMES;
      io.msg(IOL_ERROR,"input file \"%s\" has the wrong type (%s, expected %s).\n",path,names[type],names[MFBD_I16T]);
      delete[] data;
      nx=ny=0;
      if(header) delete[] header;
      header=0;
    }
}

/*
image_i16t::image_i16t(const char *dir,const char *fntmpl,io_class &io)
{
  pic=0;
  struct stat stat_buf;
  char fname[1000];
  if(fntmpl[0]=='!') sprintf(fname,fntmpl+1); else sprintf(fname,fntmpl);
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
      if(type==IOL_I16T){
        in_place_transpose((int16_t*)data,ny,nx);
        swap(nx,ny);
        pic=i16t2dim((int16_t*)data,1,nx,1,ny);
      }else{ // wrong type
        const char *names[]=IOL_TYPE_NAMES;
        io.msg(IOL_ERROR,"input file \"%s\" has the wrong type (%s, expected %s).\n",path,names[type],names[IOL_I16T]);
        delete[] data;
        nx=ny=0;
        if(header) delete[] header;
        header=0;
        if(path) delete[] path;
        path=0;
      }
  }
}
*/

image_i16t::image_i16t(const char *dir,const char *fntmpl,int num,io_class &io)
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
      if(type==MFBD_I16T){
        in_place_transpose((int16_t*)data,ny,nx);
        swap(nx,ny);
        pic=i16t2dim((int16_t*)data,1,nx,1,ny);
      }else{ // wrong type
        const char *names[]=MFBD_TYPE_NAMES;
        io.msg(IOL_ERROR,"input file \"%s\" has the wrong type (%s, expected %s).\n",path,names[type],names[MFBD_I16T]);
        delete[] data;
        nx=ny=0;
        if(header) delete[] header;
        header=0;
        if(path) delete[] path;
        path=0;
      }
  }
}

image_i16t::~image_i16t(void)
{
  if(pic) del_i16t2dim(pic,1,nx,1,ny);
}

fp_t **image_i16t::subimage(int xl,int xh,int yl,int yh)
{
  int snx=xh-xl+1,sny=yh-yl+1;
  fp_t **data=ft2dim(1,snx,1,sny);
  memset(data[1]+1,0,snx*sny*sizeof(fp_t)); // pad with 0
  int ixl=max(xl,1),ixh=min(xh,nx),iyl=max(yl,1),iyh=min(yh,ny);
  for(int x=ixl;x<=ixh;++x) for(int y=iyl;y<=iyh;++y) data[x-xl+1][y-yl+1]=pic[x][y];
  return data;
}

void image_i16t::clip(int xl,int xh,int yl,int yh)
{
  if((xl<1)||(yl<1)||(xh>nx)||(yh>ny)) exit(fprintf(stderr,"error: cannot extend image (1,1,%d,%d) to (%d,%d,%d,%d)\n",nx,ny,xl,xh,yl,yh));
  int dx=xl-1,dy=yl-1;
  int nnx=xh-xl+1,nny=yh-yl+1;
  if(dx||dy||(nx>nnx)||(ny>nny)){
    int16_t **tmp=i16t2dim(pic[1]+1,1,nnx,1,nny);
    for(int x=1;x<=nnx;++x) memmove(tmp[x]+1,pic[x+dx]+dy,nny*sizeof(int16_t));
    delete (pic+1);
    pic=tmp;
    nx=nnx;
    ny=nny;
  }
}

void image_i16t::clip(int16_t &xl,int16_t &xh,int16_t &yl,int16_t &yh,io_class *io)
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
  for(int x=xl;x<=xh;++x) memmove(pic[1]+(x-xl)*ny+1,pic[x]+yl,ny*sizeof(int16_t));
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

void image_i16t::clip(io_class *io)
// (destructively) clip the image by shifting data and updating pointers
{
  int32_t *xs=new int32_t [nx] -1,*ys=new int32_t [ny] -1; 
  memset(xs+1,0,nx*sizeof(int32_t));
  memset(ys+1,0,ny*sizeof(int32_t));
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
  for(int x=xl;x<=xh;++x) memmove(pic[1]+(x-xl)*ny+1,pic[x]+yl,ny*sizeof(int16_t));
  for(int x=2;x<=nx;++x) pic[x]=pic[x-1]+ny; // update pointers
}

static __inline fp_t sf(fp_t x) // x=0..1
{
  return 0.5*(1.0-cos(pi*x));
}

void image_i16t::paste(fp_t **sub_img,int x,int y,int *xl,int *xh,int *yl,int *yh,int mx,int my,fp_t angle,io_class &io)
{
  int snx=xh[x]-xl[x]+1,sny=yh[y]-yl[y]+1;
  fp_t **sim=sub_img;
  if(angle){
    io.msg(IOL_XNFO,"rotating (%d,%d), angle=%E\n",x,y,angle);
    sim=rotate(sub_img,1,snx,1,sny,-angle);
  }
  int margin=8;
  int offx=(xh[x]-xl[x]+1)/margin,offy=(yh[y]-yl[y]+1)/margin;
  int exl=xl[x]+offx,exh=xh[x]-offx,eyl=yl[y]+offy,eyh=yh[y]-offy;
  if((exl<1)||(eyl<1)||(exh>nx)||(eyh>ny)){
    io.msg(IOL_WARN,"cannot paste subimage (%d,%d,%d,%d) into image (1,%d,1,%d). clipping...\n",xl[x],xh[x],yl[y],yh[y],nx,ny);
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
      pic[ix][iy]+=(int16_t)(window[ix][iy]*sim[ix-xl[x]+1][iy-yl[y]+1]);
  del_ft2dim(window,exl,exh,eyl,eyh);
  if(angle) del_ft2dim(sim,1,snx,1,sny);
}

#include <fitsio.h>
#include <unistd.h>
#include "jnfo.h"

void image_i16t::write(char *fname,struct jnfo &ji,int o,io_class *io)
{
}

/*
void image_i16t::write(struct gcfg &gc,struct ocfg &cfg,char *ext)
{
  int status=0;	                         // initialize status before calling fitsio routines
  fitsfile *fptr;	                 // pointer to the FITS file; defined in fitsio.h
  char *fname=new char [strlen(cfg.output_file_name)+strlen("fits")+2];
  char extention[100];
  if(ext) sprintf(extention,"%s.fits",ext); else sprintf(extention,".fits");
  if(cfg.force_write)
    sprintf(fname,"!%s%s",cfg.output_file_name,extention);
  else
    sprintf(fname,"%s%s",cfg.output_file_name,extention);
  fits_create_file(&fptr,fname,&status); // create new file
  delete fname;
  long naxes[2]={nx,ny};                 // image width,height
  fits_create_img(fptr,SHORT_IMG,2,naxes,&status); // Create the primary array image (16-bit short integer pixels
//  long exposure = 1500.;
  char *str="Swedish 1-meter Solar Telescope";
  fits_update_key(fptr,TSTRING,"TELESCOP",str,"Telescope Description",&status); // Write a keyword; must pass the ADDRESS of the value
  fits_update_key(fptr,TSTRING,"TIME-OBS",gc.time_obs,"Average Observation Time",&status);
  fits_update_key(fptr,TSTRING,"DATE-OBS",gc.date_obs,"Observation Date",&status);
  time_t now=time(0);
  struct tm *tim=localtime(&now);
  str=new char [1000];
  sprintf(str,"%04d-%02d-%02d %02d:%02d:%02d",1900+tim->tm_year,1+tim->tm_mon,tim->tm_mday,tim->tm_hour,tim->tm_min,tim->tm_sec);
  fits_update_key(fptr,TSTRING,"DATE",str,"File Creation Date and Time",&status);
//
  str[0]=0;
  for(int k=1;k<=cfg.nc-1;++k){
    strcat(str,cfg.ch[k]->filename_template);
    strcat(str,",");
  }
  strcat(str,cfg.ch[cfg.nc]->filename_template);
  fits_update_key(fptr,TSTRING,"MFBD-TPL",str,"MFBD-config: image filename template",&status);
  sprintf(str,"Image #: ");
  for(int i=1;i<=gc.nim;++i){
    char *tmp=new char[1000];
    sprintf(tmp,"%s %d",str,cfg.image_num[i]);
    delete str;
    str=tmp;
  }
  fits_update_key(fptr,TSTRING,"MFBD-SEQ",str,"MFBD-config: image sequence numbers",&status);
  sprintf(str,"%E",cfg.lambda);
  fits_update_key(fptr,TSTRING,"MFBD-LAM",str,"MFBD-config: wavelength",&status);
  sprintf(str,"%E",gc.telescope_d);
  fits_update_key(fptr,TSTRING,"MFBD-D",str,"MFBD-config: telescope diameter",&status);
  sprintf(str,"%E",gc.arcsecperpix);
  fits_update_key(fptr,TSTRING,"MFBD-APP",str,"MFBD-config: arcseconds per pixel",&status);
  const char *basis_names[3]={"","Zernike","Karhunen-Loeve"};
  sprintf(str,"%s",basis_names[gc.basis]);
  fits_update_key(fptr,TSTRING,"MFBD-BAS",str,"MFBD-config: atmospheric model basis",&status);
  sprintf(str,"%d-%d",cfg.mode_num[1],cfg.mode_num[cfg.nm]);
  fits_update_key(fptr,TSTRING,"MFBD-MOD",str,"MFBD-config: atmospheric modes modelled",&status);
  sprintf(str,"%d,%d,0x%08X",gc.min_iter,gc.max_iter,gc.n_done_mask);
  fits_update_key(fptr,TSTRING,"MFBD-CNV",str,"MFBD-config: various iteration settings",&status);
  sprintf(str,"%d",gc.border_clip);
  fits_update_key(fptr,TSTRING,"MFBD-CLP",str,"MFBD-config: width of border not used in image normalisation",&status);
//
  long fpixel=1;
  in_place_transpose(pic[1]+1,ny,nx);
  fits_write_img(fptr,TSHORT,fpixel,nx*ny,pic[1]+1,&status); // Write the array of integers to the image
  fits_close_file(fptr,&status);        // close the file
  fits_report_error(stderr,status);     // print out any error messages
}
*/

fp_t image_i16t::noise_sigma(int offs)
{
  int Nx=nx-2*offs,Ny=ny-2*offs;
  for(int n=0;n<=30;++n){
    int mask=~(1<<n);
    if(Nx&mask) Nx=Nx&mask;
    if(Ny&mask) Ny=Ny&mask;
  }
  int xl=((nx-Nx)/2)+1,yl=((ny-Ny)/2)+1;
  int xh=xl+Nx-1,yh=yl+Ny-1;
  
//  fprintf(stderr,"nx=%d,ny=%d (%d-%d),(%d-%d)\n",Nx,Ny,xl,xh,yl,yh);
  
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

fp_t image_i16t::mean(int offs)
{
  fp_t sum=0.0;
  for(int x=offs;x<=nx-offs;++x)
    for(int y=offs;y<=ny-offs;++y) sum+=pic[x][y];
  return sum/(fp_t)((nx-2*offs)*(ny-2*offs));
}

void image_i16t::mult(fp_t f)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) pic[x][y]=(int16_t)(f*(fp_t)pic[x][y]);
}

fp_t **image_i16t::get(void)
{
  fp_t **data=ft2dim(1,nx,1,ny);
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]=(fp_t)pic[x][y];
  return data;
}

void image_i16t::put(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) pic[x][y]=(int16_t)data[x][y];
}
/*
void image_i16t::flatfield(image_t *dark,image_t *gain,image_t *psf,image_t *bgain,byte method,io_class &io)
{
  char *h;
  int nxd,nyd;
  dark->info(nxd,nyd,h);
  if((nx!=nxd)||(ny!=nyd)){
    io.msg(IOL_WARN,"not flatfielding %s: dark field dimensions (%dx%d) are different from image dimensions (%dx%d)\n",path,nxd,nyd,nx,ny);
    return;
  }
  gain->info(nxd,nyd,h);
  if((nx!=nxd)||(ny!=nyd)){
    io.msg(IOL_WARN,"not flatfielding %s: gain file dimensions (%dx%d) are different from image dimensions (%dx%d)\n",path,nxd,nyd,nx,ny);
    return;
  }
  io.msg(IOL_INFO,"darkfielding %s\n",path);
  fp_t **data=ft2dim(1,nx,1,ny);
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]=(fp_t)pic[x][y];
  data=dark->sub(data);
// 
  if(bgain&&psf){
    io.msg(IOL_INFO,"transparency correction %s\n",path);
    descatter(data,psf,bgain);
  }
//
  io.msg(IOL_INFO,"flatfielding %s\n",path);
// compute gain?
//
  data=gain->mult(data);
// fill bad pixels
  if(!(data=fillpix(data,1,nx,1,ny,method,io))){
    del_ft2dim(data,1,nx,1,ny);
    return;
  }
//
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) pic[x][y]=(int16_t)data[x][y];
  del_ft2dim(data,1,nx,1,ny);
}
*/
fp_t **image_i16t::mult(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]*=pic[x][y];
  return data;
}

fp_t **image_i16t::sub(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]-=pic[x][y];
  return data;
}

fp_t **image_i16t::add(fp_t **data)
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y) data[x][y]+=pic[x][y];
  return data;
}

image_t *image_i16t::clone(io_class &io)
{
  return new image_i16t(this,io);
}
