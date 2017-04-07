#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "types.h"
#include "version.h"
#include "fileio.h"
#include "conf.h"
#include "net.h"
#include "jnfo.h"
#include "uts.h"
#include "pack.h"
#include "img_i16t.h"
#include "img_i32t.h"
#include "img_f32t.h"
#include "img_mat.h"
#include "img.h"
#include "mem.h"
#include "array.h"
#include "slvcfg.h"
#include "compress.h"
#include "ana_io.h"
#include "job.h"
#include "mathtools.h"
#include "obs.h"
#include "atmos.h"

// -------------------------------------


chunk::chunk(int x_in,int y_in,int xl_in,int xh_in,int yl_in,int yh_in,struct slvcfg *cfg_in)
{
  memset(this,0,sizeof(chunk));
  x=x_in;
  y=y_in;
  xl=xl_in;
  xh=xh_in;
  yl=yl_in;
  yh=yh_in;
  dx=dy=0;
  offset=0;
  fswap=-1;
  cfg=cfg_in->ref(this);
}

chunk::~chunk(void)
{
  if(cfg) cfg->unref(this);
  if(buf) delete[] buf;
  if(slv_id) delete[] slv_id;
  if(dx) del_i2dim(dx,1,1,1,1);
  if(dy) del_i2dim(dy,1,1,1,1);
}

int chunk::pack(image_t ****image,image_i16t ***offx,image_i16t ***offy,fp_t ***noise_sigma,int no,int *nc,int **nim,int mls,int swapfile,off_t &swapfile_offset,pthread_mutex_t *swapfile_lock,int clvl,io_class &io)
{
  int dxl=-mls,dxh=mls,dyl=-mls,dyh=mls,offs=0,sz=0;
  for(int o=1;o<=no;++o)
    for(int k=1;k<=nc[o];++k) sz+=nim[o][k]*(2*sizeof(int)+3*sizeof(fp_t)+(xh+dxh-xl-dxl+1)*(yh+dyh-yl-dyl+1)*sizeof(float32_t));
  byte *data=new byte [sz];
  dx=i2dim(1,no,1,nc);
  dy=i2dim(1,no,1,nc);
  for(int o=1;o<=no;++o)
    for(int k=1;k<=nc[o];++k){
      fp_t fdx=0.0,fdy=0.0;
      if(offx[o][k]){
        fp_t **ox=offx[o][k]->subimage(xl,xh,yl,yh);
        int snx=xh-xl+1,sny=yh-yl+1;
        for(int ix=1;ix<=snx;++ix)
          for(int iy=1;iy<=sny;++iy) fdx+=ox[ix][iy];
        fdx/=(fp_t)(snx*sny)*100.0;
        del_ft2dim(ox,1,snx,1,sny);
      }
      if(offy[o][k]){
        fp_t **oy=offy[o][k]->subimage(xl,xh,yl,yh);
        int snx=xh-xl+1,sny=yh-yl+1;
        for(int ix=1;ix<=snx;++ix)
          for(int iy=1;iy<=sny;++iy) fdy+=oy[ix][iy];
        fdy/=(fp_t)(snx*sny)*100.0;
        del_ft2dim(oy,1,snx,1,sny);
      }
      dx[o][k]=i_round(fdx);
      dy[o][k]=i_round(fdy);
      int ixl=xl+dxl+dx[o][k]+1,ixh=xh+dxh+dx[o][k]+1;
      int iyl=yl+dyl+dy[o][k]+1,iyh=yh+dyh+dy[o][k]+1;
      fp_t ixo=rest(fdx)-(fp_t)dxl; // initial tilt will cause placement at dxl-rest(fdx)?
      fp_t iyo=rest(fdy)-(fp_t)dyl;
      for(int i=1;i<=nim[o][k];++i){
        fp_t **tmp=image[o][k][i]->subimage(ixl,ixh,iyl,iyh);
        int snx=ixh-ixl+1,sny=iyh-iyl+1;
        float32_t **img=f32t2dim(1,snx,1,sny);
        for(int ix=1;ix<=snx;++ix) for(int iy=1;iy<=sny;++iy) img[ix][iy]=tmp[ix][iy];
        del_ft2dim(tmp,1,snx,1,sny);
        offs+=::pack(data+offs,snx,0);
        offs+=::pack(data+offs,sny,0);
        offs+=::pack(data+offs,ixo,0);
        offs+=::pack(data+offs,iyo,0);
        offs+=::pack(data+offs,noise_sigma[o][k][i],0);
        offs+=::pack(data+offs,img,1,snx,1,sny,0);
        del_f32t2dim(img,1,snx,1,sny);
      }
    }
  if(offs!=sz) io.msg(IOL_ERROR,"chunk::pack: inaccurate buffer size estimate! (actual: %d > estimate: %d)\n",offs,sz);
  byte *cbuf=z_compress(data,sz,clvl,0,io);
  delete[] data;
  if(swapfile>=0)
    if(put(cbuf,sz,swapfile,swapfile_offset,swapfile_lock)<0) io.msg(IOL_WARN,"chunk::pack: failed to write to swap file (keep your fingers crossed)\n",offs,sz); 
  else 
    put(cbuf,sz);
  return 0;
}

//
// job
//

job_class::job_class(sock_class *sock,int id_in,struct slvcfg **&cfg_in,io_class &io_in)
{
  memset(this,0,sizeof(job_class));
  struct id sid(sock),mid;
  if(sid.status()<0){
    io_in.msg(IOL_WARN,"in job_class::job_class: invalid handshake (not jsub)?!\n");
    return;
  }
  if(sock->send_id(mid)<0){
    io_in.msg(IOL_WARN,"in job_class::job_class: interrupted handshake (not jsub)?!\n");
    return;
  }
  id=id_in;
  int isz;
  byte *idta=sock->recv(isz);  // job info data
  if(!idta){
    id=-1;
    io_in.msg(IOL_WARN,"in job_class::job_class: no job config data, connection broken during transfer?\n");
    return;
  }
  int asz;
  byte *adta=sock->recv(asz);  // job info data
  if(!adta){
    delete[] idta;
    id=-1;
    io_in.msg(IOL_WARN,"job_class::job_class: no ancilliary data, connection broken during transfer?\n");
    return;
  }
// auxiliary info
  int offs=0;
  offs+=unpack(adta+offs,pri,0);           // priority level
  offs+=unpack(adta+offs,verb_level,0);    // verbosity
  offs+=unpack(adta+offs,time_stamping,0); // time_stamping
  offs+=unpack(adta+offs,wd);              // working directory
  offs+=unpack(adta+offs,lgfname);         // log file name
  offs+=unpack(adta+offs,jname);           // job name
  delete[] adta;
  pri=min(PRI_MAX,pri);
//
  struct jnfo tjnfo(idta,0,io_in);
  delete[] idta;
  memcpy(&ji,&tjnfo,sizeof(struct jnfo));
  memset(&tjnfo,0,sizeof(struct jnfo));
//
  io=0;
//
  raw=new chunk* [1];
  act=new chunk* [1];
  fin=new chunk* [1];
  raw[0]=act[0]=fin[0]=0;
//
  active=0;
  pthread_mutex_init(&active_lock,0);
  pthread_mutex_init(&swapfile_lock,0);
//
  time_t now=time(0);
  localtime_r(&now,&t_sub);
  swapfilename=0;
  swapfile=-1;
//
  ppfrac=0.0;
//
  spair[0]=spair[1]=-1;
//  if(socketpair(AF_UNIX,SOCK_STREAM,0,spair)<0)
//    io->msg(IOL_ERROR,"job_class::job_class: failed to create socket pair: %s\n",strerror(errno));
}

job_class::~job_class(void)
{
  if(spair[0]>=0) close(spair[0]);
  if(spair[1]>=0) close(spair[1]);
  pthread_mutex_destroy(&active_lock);
  if(raw){
    for(int n=0;raw[n];++n) delete raw[n];
    delete[] raw;
  }
  if(act){
    for(int n=0;act[n];++n) delete act[n];
    delete[] act;
  }
  if(fin){
    for(int n=0;fin[n];++n) delete fin[n];
    delete[] fin;
  }
  if(wd) delete[] wd;
  if(lgfname) delete[] lgfname;
  if(jname) delete[] jname;
  pthread_mutex_lock(&swapfile_lock);
  if(swapfile>=0){
    close(swapfile);
    swapfile=-1;
    io->msg(IOL_INFO,"removing swapfile \"%s\" ... ",swapfilename);
    unlink(swapfilename);
    io->msg(IOL_INFO|IOL_NOID,"done\n");
    delete[] swapfilename;
  }
  pthread_mutex_unlock(&swapfile_lock);
  pthread_mutex_destroy(&swapfile_lock);
  if(io){
    if(active<5) io->msg(IOL_ERROR,"Killed\n");
    delete io;
  }
  if(cfg) cfg->unref(this);
}

struct chunk *job_class::get_chunk(void *sid_in)
{
  if(((active!=2)&&(active!=3))||(!raw[0])) return 0; // not ready/no data
  active=3; // 'A' for active
  if(struct chunk *chnk=raw[0]){
    chnk->sid=sid_in;
    array_mv(chnk,raw,act); // move chunk to the active data list
    return chnk;
  }
  return 0;
}

int job_class::unget_chunk(struct chunk *chnk)
{
  chnk->sid=0;
  array_mv(chnk,act,raw); // move chunk back to the raw data list
  if(!act[0]&&(active==6)) active=5; // no more incoming stuff: job can be killed safely
  return 0;
}

int job_class::put_chunk(struct chunk *chnk)
{
  array_mv(chnk,act,fin); // move chunk to the finished data list
  io->msg(IOL_INFO|IOL_NOID,(chnk->buf)?"[%2d,%2d]":"R[%2d,%2d]",chnk->x,chnk->y);
  if(!act[0])
    switch(active){
      case(3):{
        if(!raw[0]) return this->finalize(); // image is done...
        break;
      }
      case(6):{         // killed
        active=5;       // no more incoming stuff: job can be killed safely
        break;
      }
    }
  return 0;
}

void *activate(void *p)
{
  job_class *job=(job_class*)p;
  int policy;
  struct sched_param parm;
  pthread_getschedparam(pthread_self(),&policy,&parm);
  if(parm.sched_priority>10)
    parm.sched_priority-=10;
  else
    parm.sched_priority=0;
  pthread_setschedparam(pthread_self(),policy,&parm);
  job->start();
  return 0;
}

int job_class::activate(void)
{
  if(!this) return -1;
  if(active) return 1;
  pthread_t thread;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
  ++active;
//
  if(socketpair(AF_UNIX,SOCK_STREAM,0,spair)<0)
    io->msg(IOL_ERROR,"job_class::job_class: failed to create socket pair: %s\n",strerror(errno));
//
  if(pthread_create(&thread,&attr,::activate,this)) return -2;
  return 0;
}

#include <sys/times.h>
int job_class::start(void)
{

  t_start=time(0);
//
  char id_str[10];
  sprintf(id_str,"%06d",id);
  io=new io_class(wd,lgfname,id_str,ji.uid,ji.gid,verb_level,time_stamping);
//
  io->msg(IOL_INFO,"job ID = %s\n",id_str);

//
  for(int a=0;a<ji.na;++a){
    ji.atmos[a]->init(wd,io); // setup structure

    //ji.atmos[a]->build_from_nodes(current_model);  
//
    // ---------------------------------------------------------------------------------------------------------------------------------
    /// Time computation
    ji.atmos[a]->set_grid(0); // Sets grid to tau if positive, otherwise to geometrical scale
    io->msg(IOL_INFO,"Calculating observables...\n");
    for(int o=0;o<ji.no;++o){
      int tickspersec=sysconf(_SC_CLK_TCK);
      struct tms t_strt;
      clock_t t1=times(&t_strt);
      io->msg(IOL_INFO,"nlambda=%d, lambda[0]=%E, lambda[%d]=%E\n",ji.nlambda[o],ji.lambda[o][0],ji.nlambda[o]-1,ji.lambda[o][ji.nlambda[o]-1]);
      
      class observable *obs = ji.atmos[a]->obs_stokes(ji.el[o],ji.az[o],ji.lambda[o],ji.nlambda[o]);

      //ji.atmos[a]->obs_stokes_num_responses(ji.el[o],ji.az[o],ji.lambda[o],ji.nlambda[o],0);      
      //class observable * obs;
      //obs = new observable(4);
      //obs->readsingle("spectrum_to_fit_short.dat");
      obs->write(ji.name[o],*io);
      
      // Here we execute the fitting procedure
      //ji.atmos[a]->stokes_lm_fit(obs,ji.el[o],ji.az[o],ji.lambda[o], ji.nlambda[o]);

      //delete obs;
      struct tms t_end;
      clock_t t2=times(&t_end);
      int user=t_end.tms_utime-t_strt.tms_utime;
      int sys=t_end.tms_stime-t_strt.tms_stime;
      double utime=((double)user)/(double)tickspersec;
      double stime=((double)sys)/(double)tickspersec;
      double ttime=((double)user+(double)sys)/(double)tickspersec;
      delete obs;
      fprintf(stderr,"job time = user: %7.5f  system: %7.5f  total: %7.5f\n",utime,stime,ttime);
    }
    exit(1);
  }
 
//
  io->msg(IOL_INFO,"init done.\n");
//
  pthread_mutex_lock(&active_lock);
  ++active;
  pthread_mutex_unlock(&active_lock);
  return 0;
}

void *finalize(void *p)
{
  job_class *job=(job_class*)p;
  job->stop();
  return 0;
}

int job_class::finalize(void)
{
  pthread_t thread;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
  active=4;
  io->msg(IOL_INFO|IOL_NOID,"\n");
  if(pthread_create(&thread,&attr,::finalize,this)) io->msg(IOL_ERROR,"could not create thread\n");
  return 1;
}

struct slv_data{
  char *id;
  fp_t ut,st;
  int pc;
};

int job_class::stop(void)
{
/******************************
 * statistics                 *
 ******************************/
  time_t now=time(0);
  int64_t day=86400,hr=3600,mn=60;
  int64_t tot=(int64_t)(now-t_start);
  char *tt=new char [1000];
  sprintf(tt,"%"I64FMT"-%02"I64FMT":%02"I64FMT":%02"I64FMT,tot/day,(tot%day)/hr,(tot%hr)/mn,tot%mn);
  io->msg(IOL_INFO|IOL_NOID,"total reduction time: %s\n",tt);
  delete[] tt;
  io->msg(IOL_INFO|IOL_NOID,"breakdown by slave:\n");
  int ids=0;
//
  pthread_mutex_lock(&active_lock);
  active=5;
  pthread_mutex_unlock(&active_lock);
  pthread_mutex_lock(&swapfile_lock);
  if(swapfile>=0){
    close(swapfile);
    swapfile=-1;
    io->msg(IOL_INFO,"removing swapfile \"%s\" ... ",swapfilename);
    unlink(swapfilename);
    io->msg(IOL_INFO|IOL_NOID,"done\n");
    delete[] swapfilename;
    swapfilename=0;
  }
  pthread_mutex_unlock(&swapfile_lock);
  return 0;
}

int job_class::status(void)
{
  if(id<0) return -1;
  return 0;
}

int job_class::stat(byte *buf)
{
  int offs=0;
  pthread_mutex_lock(&active_lock);
  offs+=pack(buf+offs,active,0);
  fp_t g=ppfrac;
  pthread_mutex_unlock(&active_lock);
  offs+=pack(buf+offs,id,0);
  offs+=pack(buf+offs,pri,0);
  int n=-1,m=-1;
  while(fin[++n]);
  while(raw[++m]);
  fp_t f=(n)?1.0:0.0;
  offs+=pack(buf+offs,f,0);
  offs+=pack(buf+offs,g,0);
  offs+=pack(buf+offs,ji.uname);
  offs+=pack(buf+offs,jname);
  offs+=pack(buf+offs,t_sub.tm_sec,0);
  offs+=pack(buf+offs,t_sub.tm_min,0);
  offs+=pack(buf+offs,t_sub.tm_hour,0);
  offs+=pack(buf+offs,t_sub.tm_year,0);
  offs+=pack(buf+offs,t_sub.tm_mon,0);
  offs+=pack(buf+offs,t_sub.tm_mday,0);
  return offs;
}

int job_class::kill(void)
{
  pthread_mutex_lock(&active_lock);
  switch(active){
    case(0): active=5; break;  
    case(1): active+=5; break; // schedule for kill
    case(2): active=5; break;  
    case(3): active=(act[0])?6:5; break;  // if active patches -> wait for them, else: kill
  }
  pthread_mutex_unlock(&active_lock);
  return 0;
}

