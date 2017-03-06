#include <string.h>
#include <sys/stat.h>

#include "types.h"
#include "io.h"
#include "img.h"
#include "mem.h"
#include "uts.h"
#include "fileio.h"
#include "img_f32t.h"
#include "img_mat.h"

image_mat::image_mat(const char *dir,const char *fnm,io_class &io):image_t()
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
    nx=ny=nz=0;
  }else{
    int type;
    if(byte *data=read_file(path,nx,ny,nz,type,header,io))
      if(type==MFBD_F32T){
        float32_t ***p=f32t3dim((float32_t*)data,1,nz,1,ny,1,nx);
        pic=f32t3dim(1,nx,1,ny,1,nz);
        for(int x=1;x<=nx;++x)
          for(int y=1;y<=ny;++y)
            for(int z=1;z<=nz;++z) pic[x][y][z]=p[z][y][x];
        del_f32t3dim(p,1,nz,1,ny,1,nx);
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

image_mat::~image_mat(void)
{
  if(pic) del_f32t3dim(pic,1,nx,1,ny,1,nz);
}

image_t *image_mat::modgain(image_t *gain,fp_t *S,int row,int width)
{
  fp_t **data=ft2dim(1,nx,1,ny);
  memset(data[1]+1,0,nx*ny*sizeof(fp_t));
  data=gain->add(data);
//
  fp_t N=0.0,sum=0.0;
  int zoffs=(row-1)*width;
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y)
      if(data[x][y]){ // non-zero?
        fp_t polI=0.0;
        for(int z=1;z<=width;++z) polI+=S[z]*pic[x][y][zoffs+z];
        data[x][y]/=polI;
        sum+=data[x][y];
        N+=1.0;
      }
// renormalize (do we really want this?)
  for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) data[x][y]*=N/sum;
//
  image_t *res=new image_f32t(nx,ny);
  res->put(data);
//
  del_ft2dim(data,1,nx,1,ny);
//
  return res;
}

fp_t **image_mat::mult(fp_t **data) // apply O(nz) polynomial response function 
{
  for(int x=1;x<=nx;++x)
    for(int y=1;y<=ny;++y){
      fp_t sum=pic[x][y][nz];
      for(int z=nz-1;z>=1;--z) sum=pic[x][y][z]+data[x][y]*sum;
      data[x][y]=sum;
    }
  return data;
}
