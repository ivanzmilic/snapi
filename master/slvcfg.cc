#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include "array.h"
//#include "config.h"
//#include "modes.h"
#include "pack.h"
#include "slvcfg.h"

slvcfg::slvcfg(byte *data_in,int size_in,slvcfg **&p_in,io_class &io):p(p_in)
{
  size=size_in;
  memcpy(data=new byte [size],data_in,size);
//
/*
  struct config *conf=new config(data,0,io); // no byte swapping (we are master)...
  struct klmc *kl_cfs=0;
  if(conf->kl_min_mode||conf->kl_max_mode){
    io.msg(IOL_INFO,"calculating Karhunen-Loeve coefficients\n");
    kl_cfs=kl_cfg(conf->kl_min_mode,conf->kl_max_mode);
  }
  struct modes **mode=new modes* [conf->no]-1;
  for(int o=1;o<=conf->no;++o)
    mode[o]=new modes(kl_cfs,conf->lambda[o],conf->lim_freq[o]/2.0,conf->nph[o],conf->basis,conf->nm,conf->mode_num,conf->nc[o],conf->ndo[o],conf->dorder[o],conf->dtype[o],conf->kl_min_mode,conf->kl_max_mode,conf->svd_reg,conf->angle[o],conf->pupil[o],io);
//
  msize=0;
  mdata=0;
  for(int o=1;o<=conf->no;++o){
    int32_t sz;
    byte *q=mode[o]->compress(sz,1,0,io); // compression level 1, no byte swapping (we are master)
    delete mode[o];
    byte *tmp=new byte [msize+sz+sizeof(int32_t)];
    if(msize){
      memcpy(tmp,mdata,msize);
      delete[] mdata;
    }
    mdata=tmp;
    int32_t offs=pack(mdata+msize,sz,0);
    memcpy(mdata+msize+offs,q,sz);
    delete[] q;
    msize+=sz+offs;
  }
  if(kl_cfs) kl_uncfg(kl_cfs,conf->kl_min_mode,conf->kl_max_mode);
  delete[] (mode+1);
  delete conf;
*/
//
  refs=new void* [1];
  refs[0]=0;
  pthread_mutex_init(&active_lock,0);
}

slvcfg::slvcfg(byte *data_in,int size_in,slvcfg **&p_in):p(p_in)
{
  size=size_in;
  memcpy(data=new byte [size],data_in,size);
  msize=0;
  mdata=0;
//
  refs=new void* [1];
  refs[0]=0;
  pthread_mutex_init(&active_lock,0);
}
  
slvcfg::~slvcfg(void)
{
  pthread_mutex_destroy(&active_lock);
  if(refs) delete[] refs;
  if(data) delete[] data;
  if(mdata) delete[] mdata;
}
  
slvcfg *slvcfg::ref(void *r)
{
  pthread_mutex_lock(&active_lock); // can be from any thread...
  array_add(r,refs);
  pthread_mutex_unlock(&active_lock);
  return this;
}
  
void slvcfg::unref(void *r)
{
  pthread_mutex_lock(&active_lock); // can be from any thread...
  array_del(r,refs);                // called from main but different job thread may be running too
  pthread_mutex_unlock(&active_lock);
  if(!refs[0]){        // no more references left
    array_del(this,p); // remove from pool, no lock needed (we are in main)
    delete this;       // delete oneself
  }
}

struct slvcfg *slvcfg::same(byte *data_in,int size_in)
{
  if(size_in==size)
    if(!memcmp(data,data_in,size)) return this;
  return 0;
}

int slvcfg::send(sock_class *sock,int thr_sock)
{
  if(sock->send(data,size,thr_sock)<0) return -1;
  return sock->send(mdata,msize,thr_sock);
}

/*
struct modes *slvcfg::mde(int o,io_class &io)
{
  int offs=0;
  for(int i=1;i<o;++i){
    int32_t sz;
    offs+=unpack(mdata+offs,sz,0);
    if((offs+sz)>msize){
      io.msg(IOL_ERROR,"slvcfg::modes: no data to read mode %d\n",o);
      return 0;
    }
    offs+=sz;
  }
  int32_t sz;
  offs+=unpack(mdata+offs,sz,0); // the size of the one we want
  if((offs+sz)>msize) return 0;
  return new modes(mdata+offs,0,io);
}

void slvcfg::modeinfo(io_class &io)
{
  int32_t offs=0;
  for(int o=1;offs<msize;++o){
    int32_t sz;
    offs+=unpack(mdata+offs,sz,0);
    offs+=sz;
    io.msg(IOL_XNFO,"modes[%d]=%d bytes\n",o,sz);
  }
}

struct config *slvcfg::cfg(io_class &io)
{
  return new config(data,0,io);
}
*/
