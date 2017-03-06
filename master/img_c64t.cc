#include <string.h>
#include "fftw3.h"

#include "types.h"
#include "io.h"
#include "img.h"
#include "mem.h"
#include "img_c64t.h"

image_c64t::image_c64t(complex_t **data_in,int nx_in,int ny_in,fftw_plan fplan_in,fftw_plan bplan_in,io_class &io):image_t()
{
  nx=nx_in;
  ny=ny_in;
  fplan=fplan_in;
  bplan=bplan_in;
  pic=ct2dim(1,nx,1,ny/2+1);
  memcpy(pic[1]+1,data_in[1]+1,nx*(ny/2+1)*sizeof(complex_t));
}

image_c64t::~image_c64t(void)
{
  if(pic) del_ct2dim(pic,1,nx,1,ny/2+1);
  fftw_destroy_plan(fplan);
  fftw_destroy_plan(bplan);
}

fp_t **image_c64t::convolve(fp_t **cimg)
{
  complex_t **ft=ct2dim(1,nx,1,ny/2+1);
  fftw_execute_dft_r2c(fplan,cimg[1]+1,(fftw_complex*)(ft[1]+1));
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny/2+1;++y){
      fp_t tmp=ft[x][y].re;
      ft[x][y].re=pic[x][y].re*tmp-pic[x][y].im*ft[x][y].im;
      ft[x][y].im=pic[x][y].im*tmp+pic[x][y].re*ft[x][y].im;
    }
  fftw_execute_dft_c2r(bplan,(fftw_complex*)(ft[1]+1),cimg[1]+1);
  del_ct2dim(ft,1,nx,1,ny/2+1);
// renormalize
  fp_t np=nx*ny;
  for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) cimg[x][y]/=np;
//
  return cimg;
}

image_t *image_c64t::clone(io_class &io)
{
  int nx_o=nx,ny_o=ny;
  complex_t **pic_o=ct2dim(1,nx,1,ny/2+1);
  fp_t **img=ft2dim(1,nx,1,ny);
  fftw_plan fplan_o=fftw_plan_dft_r2c_2d(nx,ny,img[1]+1,(fftw_complex*)(pic_o[1]+1),FFTW_FORWARD|FFTW_ESTIMATE);
  fftw_plan bplan_o=fftw_plan_dft_c2r_2d(nx,ny,(fftw_complex*)(pic_o[1]+1),img[1]+1,FFTW_BACKWARD|FFTW_ESTIMATE);
  del_ft2dim(img,1,nx,1,ny);
  image_t *i=new image_c64t(pic_o,nx_o,ny_o,fplan_o,bplan_o,io);
  del_ct2dim(pic_o,1,nx,1,ny/2+1);
  return i;
}
