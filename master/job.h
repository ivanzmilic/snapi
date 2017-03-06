#ifndef __JOB_H__ // __JOB_H__
#define __JOB_H__

#include <pthread.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include "types.h"
#include "net.h"
#include "io.h"
#include "mem.h"
#include "img.h"
#include "img_i16t.h"
#include "img_mat.h"
#include "slvcfg.h"
#include "jnfo.h"

#define PRI_MAX 100

struct mcs{
  int32_t nc,nm;
  int32_t *nim;
  int32_t xl,xh,yl,yh,*dx,*dy,nph;
  int32_t npsf,nobj,nres,nalpha,ndiv;
  fp_t **img,***psf,***obj,***res,**alpha,***diversity;
//
  mcs(int xl_in,int xh_in,int yl_in,int yh_in,int *dx_in,int *dy_in,int nc_in,int *nim_in){
    xl=xl_in;
    xh=xh_in;
    yl=yl_in;
    yh=yh_in;
    int nx=xh-xl+1,ny=yh-yl+1;
    nc=nc_in;
    dx=new int32_t [nc]-1;
    memcpy(dx+1,dx_in+1,nc*sizeof(int32_t));
    dy=new int32_t [nc]-1;
    memcpy(dy+1,dy_in+1,nc*sizeof(int32_t));
    nim=new int32_t [nc]-1;
    memcpy(nim+1,nim_in+1,nc*sizeof(int32_t));
    img=0;
    alpha=0;
    diversity=psf=obj=res=0;
    nalpha=npsf=nobj=nres=ndiv=0;
  }
  ~mcs(void){
    int nx=xh-xl+1,ny=yh-yl+1;
    if(dx) delete[] (dx+1);
    if(dy) delete[] (dy+1);
    if(nim) delete[] (nim+1);
    if(img) del_ft2dim(img,1,nx,1,ny);
    if(psf){
      for(int i=1;i<=npsf;++i) del_ft2dim(psf[i],1,nx,1,ny);
      delete[] (psf+1);
    }
    if(obj){
      for(int i=1;i<=nobj;++i) del_ft2dim(obj[i],1,nx,1,ny);
      delete[] (obj+1);
    }
    if(res){
      for(int i=1;i<=nres;++i) del_ft2dim(res[i],1,nx,1,ny);
      delete[] (res+1);
    }
    if(diversity){
      for(int i=1;i<=ndiv;++i) del_ft2dim(diversity[i],1,nx,1,ny);
      delete[] (diversity+1);
    }
    if(alpha){
      for(int i=1;i<=nalpha;++i) delete[] (alpha[i]+1);
      delete[] (alpha+1);
    }
  }
  void add_img(fp_t **img_in,int nx_in,int ny_in){
    int nx=xh-xl+1,ny=yh-yl+1;
    img=ft2dim(1,nx_in,1,ny_in);
    memcpy(img[1]+1,img_in[1]+1,nx*ny*sizeof(fp_t));
  }
  fp_t ***add_qty(fp_t ***p,int &np,fp_t **qty_in,int nx_in,int ny_in){
    fp_t ***tmp=new fp_t** [np+1]-1;
    if(np){
      memcpy(tmp+1,p+1,np*sizeof(fp_t**));
      delete[] (p+1);
    }
    ++np;
    tmp[np]=ft2dim(1,nx_in,1,ny_in);
    memcpy(tmp[np][1]+1,qty_in[1]+1,nx_in*ny_in*sizeof(fp_t));
    return tmp;
  }
  fp_t ***add_qty_avg(fp_t ***p,int &np,fp_t **qty_in,int nx_in,int ny_in,int n_add){
    fp_t ***tmp=p;
    if(!np){
      tmp=new fp_t** [np=1]-1;
      tmp[np]=ft2dim(1,nx_in,1,ny_in);
      memset(tmp[np][1]+1,0,nx_in*ny_in*sizeof(float32_t));
    }
    fp_t f=1.0/(fp_t)n_add;
    for(int x=1;x<=nx_in;++x) for(int y=1;y<=ny_in;++y) tmp[np][x][y]=(1.0-f)*tmp[np][x][y]+f*qty_in[x][y]; // type insensitive copy...
    return tmp;
  }
  void add_psf(fp_t **psf_in,int nx_in,int ny_in){
    psf=add_qty(psf,npsf,psf_in,nx_in,ny_in);
  }
  void add_psf_avg(fp_t **psf_in,int nx_in,int ny_in,int n_add){
    psf=add_qty_avg(psf,npsf,psf_in,nx_in,ny_in,n_add);
  }
  void add_obj(fp_t **obj_in,int nx_in,int ny_in){
    obj=add_qty(obj,nobj,obj_in,nx_in,ny_in);
  }
  void add_res(fp_t **res_in,int nx_in,int ny_in){
    res=add_qty(res,nres,res_in,nx_in,ny_in);
  }
  void add_diversity(fp_t **div_in,int nx_in,int ny_in){
    diversity=add_qty(diversity,ndiv,div_in,nx_in,ny_in);
  }
  void add_alpha(fp_t *alpha_in,int nm_in){
    fp_t **tmp=new fp_t* [nalpha+1]-1;
    if(nalpha){
      memcpy(tmp+1,alpha+1,nalpha*sizeof(fp_t*));
      delete[] (alpha+1);
    }
    alpha=tmp;
    ++nalpha;
    alpha[nalpha]=new fp_t [nm_in]-1;
    memcpy(alpha[nalpha]+1,alpha_in+1,nm_in*sizeof(fp_t));
    nm=nm_in;
  }
  off_t write(int fd,void *buf,int sz,io_class *io){
    if(::write(fd,buf,sz)!=sz) io->msg(IOL_ERROR,"mcs::write: failed to write: %s\n",strerror(errno));
    return sz;
  }
  off_t write(int fd,io_class *io){
    off_t sz=0;
    int nx=xh-xl+1,ny=yh-yl+1;
    float32_t **tmp=f32t2dim(1,nx,1,ny);
//
    sz+=write(fd,&xl,sizeof(int32_t),io);
    sz+=write(fd,&xh,sizeof(int32_t),io);
    sz+=write(fd,&yl,sizeof(int32_t),io);
    sz+=write(fd,&yh,sizeof(int32_t),io);
    sz+=write(fd,&nc,sizeof(int32_t),io);
    sz+=write(fd,nim+1,nc*sizeof(int32_t),io);
    sz+=write(fd,dx+1,nc*sizeof(int32_t),io);
    sz+=write(fd,dy+1,nc*sizeof(int32_t),io);
    uint08_t has_img=(img!=0);
//    fprintf(stderr,"image: %X %d\n",img,has_img);
    sz+=write(fd,&has_img,1,io);
    if(has_img){
      for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) tmp[x][y]=img[x][y];
      sz+=write(fd,tmp[1]+1,nx*ny*sizeof(float32_t),io);
    }
    sz+=write(fd,&npsf,sizeof(int32_t),io);
    if(npsf){
      for(int n=1;n<=npsf;++n){
        for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) tmp[x][y]=psf[n][x][y];
        sz+=write(fd,tmp[1]+1,nx*ny*sizeof(float32_t),io);
      }
    }
    sz+=write(fd,&nobj,sizeof(int32_t),io);
    if(nobj){
      for(int n=1;n<=nobj;++n){
        for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) tmp[x][y]=obj[n][x][y];
        sz+=write(fd,tmp[1]+1,nx*ny*sizeof(float32_t),io);
      }
    }
    sz+=write(fd,&nres,sizeof(int32_t),io);
    if(nres){
      for(int n=1;n<=nres;++n){
        for(int x=1;x<=nx;++x) for(int y=1;y<=ny;++y) tmp[x][y]=res[n][x][y];
        sz+=write(fd,tmp[1]+1,nx*ny*sizeof(float32_t),io);
      }
    }
    sz+=write(fd,&nalpha,sizeof(int32_t),io);
    if(nalpha){
      sz+=write(fd,&nm,sizeof(int32_t),io);
      float32_t *a=new float32_t [nm]-1;
      for(int n=1;n<=nalpha;++n){
        for(int m=1;m<=nm;++m) a[m]=alpha[n][m];
        sz+=write(fd,a+1,nm*sizeof(float32_t),io);
      }
      delete[] (a+1);
    }
//
    sz+=write(fd,&ndiv,sizeof(int32_t),io);
    if(ndiv){
      float32_t **tmph=f32t2dim(1,nx/2,1,ny/2);
      for(int n=1;n<=ndiv;++n){ // FIXME: add limits (lowest,highest realization)?
        for(int x=1;x<=nx/2;++x) for(int y=1;y<=ny/2;++y) tmph[x][y]=diversity[n][x][y];
        sz+=write(fd,tmph[1]+1,nx*ny*sizeof(float32_t)/4,io);
      }
      del_f32t2dim(tmph,1,nx/2,1,ny/2);
    }
//
    del_f32t2dim(tmp,1,nx,1,ny);
//
    return sz;
  }
};

struct chunk{
  int x,y,xl,xh,yl,yh,**dx,**dy;
  time_t tstart,tend;
  void *sid;
  struct slvcfg *cfg;
// input/output
  byte *buf;
  int bsz;
  int fswap;
  off_t offset;
  pthread_mutex_t *mutex;
  char *slv_id;
  fp_t s_time,u_time;
//
  chunk(int,int,int,int,int,int,struct slvcfg*);
  ~chunk(void);
  int pack(image_t ****image,image_i16t ***offx,image_i16t ***offy,fp_t ***noise_sigma,int no,int *nc,int **nim,int mls,int swapfile,off_t &swapfile_offset,pthread_mutex_t *swapfile_lock,int clvl,io_class &io);
  void put(byte *buf_in,int bsz_in){
    buf=buf_in;
    bsz=bsz_in;
    fswap=-1;
  }
  int put(byte *buf_in,int bsz_in,int fswap_in,off_t &offset_in,pthread_mutex_t *mutex_in){
    fswap=fswap_in;
    mutex=mutex_in;
    pthread_mutex_lock(mutex);
    offset=offset_in;
    lseek(fswap,offset,SEEK_SET);
    if(write(fswap,buf_in,bsz_in)<0){ // failed to write: keep in RAM
      pthread_mutex_unlock(mutex);
      buf=buf_in;
      bsz=bsz_in;
      fswap=-1;
      return -1;
    }
    offset_in+=bsz_in;
    pthread_mutex_unlock(mutex);
    bsz=bsz_in;
    delete[] buf_in;
    buf=0;
    return 0;
  }
  int get(byte *&buf_o,int &bsz_o){
    if((!buf)&&(fswap>=0)){
      buf=new byte [bsz];
      pthread_mutex_lock(mutex);
      lseek(fswap,offset,SEEK_SET);
      if(read(fswap,buf,bsz)<bsz){
        pthread_mutex_unlock(mutex);
        return -1;
      }
      pthread_mutex_unlock(mutex);
    }
//    fprintf(stderr,"chunk_get: %X %X\n",this,buf);
    bsz_o=bsz;
    buf_o=buf;
    return 0;
  }
  void result(byte *buf_i,int bsz_i){
    if(buf){
//      fprintf(stderr,"chunk_free: %X %X\n",this,buf);
      delete[] buf;
    }
    buf=buf_i;
    bsz=bsz_i;
  }
  void free(void){
    if((buf)&&(fswap>=0)){ // free buffer when reading from fswap...
//      fprintf(stderr,"chunk_free: %X %X\n",this,buf);
      delete[] buf;
      buf=0;
    }
  }
};

struct stats{
  struct stat *data;
  int tstart,tcpu; // starting time, cpu time on all slaves
};

class job_class{
  int id,pri;
  int verb_level,active,time_stamping;
  char *wd,*lgfname,*jname,*swapfilename;
//
  struct jnfo ji;                 // job setup info
  struct chunk **raw,**act,**fin; // chunk info 
//
  struct slvcfg *cfg;             // compressed slave config data needed for this job
  struct tm t_sub; 
  time_t t_start;
//
  io_class *io;
  int swapfile;
//
  int nx,ny;
  fp_t ppfrac;
//
  pthread_mutex_t active_lock;
  pthread_mutex_t swapfile_lock;
  int spair[2];
//  int abort(int,int,image_i16t***,image_i16t***,image_t***,uint08_t**,image_t***,uint08_t**,image_t***,image_mat***,uint08_t**,image_t***,uint08_t**,image_t***,uint08_t**,image_t****,fp_t***,jobthrd**);
public:
  job_class(sock_class*,int,struct slvcfg**&,io_class&);
  ~job_class(void);
  struct chunk *get_chunk(void*); // get work
  int unget_chunk(struct chunk*); // work was not done (slave failed?)
  int put_chunk(struct chunk*);   // give back processed data
  int state(void){ return active; } // completion info
  int status(void);               // status info
  int priority(void){
    return (this)?pri:PRI_MAX+1;  // identify end of queue by less than maximum priority number
  }
  void msg(int level,const char *mesg){
    io->msg(level,mesg);
  }
  int stat(byte*);                // status info
  int activate(void);             // activate task
  int start(void);                // start task
  int finalize(void);             // finalize task
  int stop(void);                 // finish task
  int kill(void);                 // kill task
  int jid(void){
    return (this)?id:-1;
  }
  void task_return(uint08_t cmd){
    ::send(spair[1],&cmd,1,0);
  }
  int08_t task_recv(void){
    int08_t cmd;
    if(::recv(spair[0],&cmd,1,0)==1) return cmd;
    return 0;
  }
};

#endif          // __JOB_H__
