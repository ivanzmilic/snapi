#include <sys/types.h>
#include <sys/socket.h>
#include <string.h>
#include <pthread.h>
#include <time.h>

#include "job.h"
#include "net.h"
#include "pack.h"
#include "mem.h"
#include "compress.h"
#include "slv.h"

void *slv_thread_start(void *p)
{
  slv_class *slv=(slv_class*)p;
  slv->helper();
  return 0;
}

slv_class::slv_class(int id_in,void *(*func)(void*),io_class &io_in):io(io_in)
{
  t_start=time(0);
  t_last=t_start;
//
  sock=0;
  job=0;
  cfg=0;
  chnk=0;
  flop=0;
  flops=-1;
  utime=0.0;
  stime=0.0;
  myid=id_in;
//
  int sck[2];
  socketpair(AF_UNIX,SOCK_STREAM,0,sck);
  mgr_sock=sck[0];
  thr_sock=sck[1];
//
  pthread_t thread;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
  if(pthread_create(&thread,&attr,(func)?func:slv_thread_start,this))
    io.msg(IOL_ERROR,"could not create thread\n");
//
  pthread_mutex_init(&time_lock,0);
}

slv_class::slv_class(sock_class *sock_in,int id_in,int clvl_in,io_class &io_in):sid(sock_in),io(io_in)
{
  t_start=time(0);
  t_last=t_start;
//
  sock=sock_in;
  struct id mid;
  sock->send_id(mid);
  job=0;
  cfg=0;
  chnk=0;
  clvl=clvl_in;
  flop=0;
  flops=-1;
  utime=0.0;
  stime=0.0;
  myid=id_in;
//
  int sck[2];
  socketpair(AF_UNIX,SOCK_STREAM,0,sck);
  mgr_sock=sck[0];
  thr_sock=sck[1];
//
  pthread_t thread;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
  if(pthread_create(&thread,&attr,slv_thread_start,this))
    io.msg(IOL_ERROR,"could not create thread\n");
//
  pthread_mutex_init(&time_lock,0);
}

void slv_class::send_cmd(int s,byte cmd)
{
  ::send(s,&cmd,1,0);
}

byte slv_class::recv_cmd(int s)
{
  byte cmd=0;
  ::recv(s,&cmd,1,0);
  return cmd;
}

slv_class::~slv_class(void)
{
  if(chnk&&job) job->unget_chunk(chnk); // were we doing work?: put back chunk for someone else
  chnk=0;
  job=0;
//
  pthread_mutex_destroy(&time_lock);
  close(mgr_sock);
  close(thr_sock);
  delete sock;
  if(cfg) cfg->unref(this);
  time_t now=time(0);
  struct tm *tstr=localtime(&t_start);
  char *ts=new char [1000];
  sprintf(ts,"%04d.%02d.%02d %02d:%02d:%02d",1900+tstr->tm_year,1+tstr->tm_mon,tstr->tm_mday,tstr->tm_hour,tstr->tm_min,tstr->tm_sec);
  struct tm *tend=localtime(&now);
  char *te=new char [1000];
  sprintf(te,"%04d.%02d.%02d %02d:%02d:%02d",1900+tend->tm_year,1+tend->tm_mon,tend->tm_mday,tend->tm_hour,tend->tm_min,tend->tm_sec);
  int64_t day=86400,hr=3600,mn=60;
  int64_t tot=(int64_t)(now-t_start);
  char *tt=new char [1000];
  sprintf(tt,"%" I64FMT "-%02" I64FMT ":%02" I64FMT ":%02" I64FMT,tot/day,(tot%day)/hr,(tot%hr)/mn,tot%mn);
  int64_t istime=(int64_t)stime;
  char *tsys=new char [1000];
  sprintf(tsys,"%" I64FMT "-%02" I64FMT ":%02" I64FMT ":%02" I64FMT,istime/day,(istime%day)/hr,(istime%hr)/mn,istime%mn);
  int64_t iutime=(int64_t)utime;
  char *tusr=new char [1000];
  sprintf(tusr,"%" I64FMT "-%02" I64FMT ":%02" I64FMT ":%02" I64FMT,iutime/day,(iutime%day)/hr,(iutime%hr)/mn,iutime%mn);
  io.msg(IOL_INFO,"%s:  %s pid=%d (%s) - (%s) [%s] user: %s sys: %s\n",te,sid.name,sid.pid,ts,te,tt,tusr,tsys);
  delete[] ts;
  delete[] te;
  delete[] tt;
  delete[] tsys;
  delete[] tusr;
}

void slv_class::reset_time(void)
{
  pthread_mutex_lock(&time_lock);
  t_last=time(0);
  pthread_mutex_unlock(&time_lock);
}

int slv_class::idle_time(void)
{
  time_t now=time(0);
  pthread_mutex_lock(&time_lock);
  time_t t_idle=now-t_last;
  pthread_mutex_unlock(&time_lock);
  return t_idle;
}

void slv_class::kill(void)
{
  send_cmd(mgr_sock,CMD_DEL_SLV);
}

void slv_class::time_out(void)
{
  send_cmd(mgr_sock,CMD_QUIT);
}

int slv_class::get_task(job_class *job_in)
{
  if(chnk=job_in->get_chunk(this)){
    job=job_in;
    if(cfg!=chnk->cfg){ // not same config: reconfigure slave
      if(cfg) cfg->unref(this);
      cfg=chnk->cfg->ref(this);
      send_cmd(mgr_sock,CMD_THREAD_CFG);
    }else send_cmd(mgr_sock,CMD_THREAD_TSK);
    return 1;
  }
  return 0;
}

void slv_class::put_task(void)
{
  char tmp_str[1000];
  sprintf(tmp_str,"%s[%d]",sid.name,sid.pid);
  chnk->slv_id=strcpy(new char [strlen(tmp_str)+1],tmp_str);
  if(chnk->bsz){
    int32_t csize;
    unpack(chnk->buf,csize,0);
    if(chnk->bsz>(csize+2*sizeof(int32_t))){ // extract some extra info
      fp_t u_time,s_time;
      int64_t lflops;
      int offs=csize+2*sizeof(int32_t);
      offs+=unpack(chnk->buf+offs,u_time,0);
      offs+=unpack(chnk->buf+offs,s_time,0);
      offs+=unpack(chnk->buf+offs,lflops,0);
      stime+=s_time;
      utime+=u_time;
      chnk->s_time=s_time;
      chnk->u_time=u_time;
      flop+=(int64_t)((fp_t)lflops*(s_time+u_time));
      flops=lflops;
    }
  }
//
  job->put_chunk(chnk); // give processed data back to job
  chnk=0;
  job=0;
}

#include <stdio.h>

byte slv_class::io_cmd(void)
{
  return recv_cmd(mgr_sock);
}

int slv_class::stat(byte *buf,const char rs)
{
  int offs=0;
  int s_id=id();
  offs+=pack(buf+offs,s_id,0);
  offs+=pack(buf+offs,sid.name);
  offs+=pack(buf+offs,sid.pid,0);
  offs+=pack(buf+offs,sid.version);
  offs+=pack(buf+offs,sid.os);
  offs+=pack(buf+offs,sid.machine);
  offs+=pack(buf+offs,flops,0);
  buf[offs++]=rs;
  return offs;
}

void slv_class::helper(void) // thread main loop
{
  fdset fds;
  fds.add(thr_sock,&thr_sock);
  sock->add(fds,sock);
  int quit=0;
  do{
    fds.act();
    for(void *token=fds.next();token;token=fds.next()){
      if(token==&thr_sock){ // message from manager
	switch(int cmd=recv_cmd(thr_sock)){
	  case(CMD_QUIT):{
            quit=1;
            while(token=fds.next());   // empty the queue
	    break;
	  }
	  case(CMD_DEL_SLV):{
            sock->send_cmd(CMD_QUIT);  // kill slave if it is listening
            quit=1;
            while(token=fds.next());   // empty the queue
	    break;
	  }
	  case(CMD_THREAD_CFG):{
            reset_time();
            sock->send_cmd(CMD_SLV_CFG);
            io.msg(IOL_INFO,"sending task\n");
//            cfg->modeinfo(io);
            if(cfg->send(sock,thr_sock)<0){
              quit=1;
              while(token=fds.next());   // empty the queue
	    }
            break;
	  }
	  case(CMD_THREAD_TSK):{
            reset_time();
            sock->send_cmd(CMD_SLV_TSK);
            byte *buf;
            int bsz;
            chnk->get(buf,bsz,io);
            if(sock->send(buf,bsz,thr_sock)<0){
              quit=1;
              while(token=fds.next());   // empty the queue
	    }
            chnk->free();
            break;
	  }
	  default: io.msg(IOL_WARN,"ignoring unknown command from manager %d\n",cmd);
	}
      }
      if(token==sock){ // message from slave
	switch(sock->recv_cmd()){
          case(CMD_QUIT):
          case(CMD_CRASH):{  // slave died!
            quit=1;
            while(token=fds.next());   // empty the queue
	    break;
          }
	  case(CMD_SLV_IO):{
            reset_time();
            int size,offs=0,type;
            char *mesg;
            byte *data=sock->recv(size);
            offs+=unpack(data+offs,type,0);
            offs+=unpack(data+offs,mesg);
            char *xmesg=new char [strlen(mesg)+strlen(sid.name)+8];
            sprintf(xmesg,"%s(%d):",sid.name,sid.pid,mesg);
            delete[] data;
            if(job) job->msg(type,xmesg); else io.msg(type,xmesg);
            delete[] xmesg;
            delete[] mesg;
            break;
          }
          case(CMD_SLV_RES):{ // receive data
            reset_time();
	    int size;
	    byte *buf=sock->recv(size,thr_sock);
            if(!buf){         // active slave (connection) died! image restoration was not done...
              quit=1;
            }else{ // result is here
              chnk->result(buf,size);
              send_cmd(thr_sock,CMD_THREAD_RES);
	    }
	    break;
	  }
          case(CMD_SLV_REJ):{ // slave rejected the data, don't try again
	    reset_time();
            chnk->result(0,0);
            send_cmd(thr_sock,CMD_THREAD_RES);
	    break;
	  }
	  case(CMD_SLV_CFG):{
            reset_time();
            sock->send_cmd(CMD_SLV_TSK);
            byte *buf;
            int bsz;
            chnk->get(buf,bsz,io);
            if(sock->send(buf,bsz,thr_sock)<0){
              quit=1;
              while(token=fds.next());   // empty the queue
	    }
            chnk->free();
            break;
	  }
          default: job->msg(IOL_WARN,"ignoring unknown command\n");
        }
      }
    }
  }while(!quit);
  sock->del(fds);
  fds.rem(thr_sock);
  send_cmd(thr_sock,CMD_QUIT);
}

