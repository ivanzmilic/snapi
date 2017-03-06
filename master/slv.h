#ifndef __SLV_H__ // __SLV_H__
#define __SLV_H__

#include <pthread.h>
#include <time.h>
#include "job.h"
#include "net.h"
#include "slvcfg.h"
#include "fdset.h"

#define CMD_THREAD_CFG  1
#define CMD_THREAD_TSK  2
#define CMD_THREAD_RES  3

class slv_class{
  job_class *job;
  struct chunk *chnk;
  struct slvcfg *cfg;        // compressed slave config data needed for this chunk
  sock_class *sock;
protected:
  io_class &io;
  struct id sid;
//
  int64_t flops,flop;
  fp_t utime,stime;
  int myid;
//
  int mgr_sock,thr_sock;
//
  int clvl;                  // compression level
//
  pthread_mutex_t time_lock;
  time_t t_start,t_last;
  void send_cmd(int,byte);
  byte recv_cmd(int);
public:
  slv_class(int id_in,void *(*func)(void*),io_class &io_in);
  slv_class(sock_class *sock_in,int sid_in,int clvl_in,io_class &io_in);
  virtual ~slv_class(void);
  void add(fdset &fds){
    fds.add(mgr_sock,this);
  }
  void del(fdset &fds){
    fds.rem(mgr_sock);
  }
  int64_t id(void){
    return myid;
  }
  void reset_time(void);
  virtual int idle_time(void);
  virtual int get_task(job_class *job_in);
  virtual void put_task(void);
  virtual int stat(byte*,const char);
  virtual byte is_busy(void){ return 1; }
  void kill(void);
  void time_out(void);
  byte io_cmd(void);
  void helper(void);
};

#endif             // __SLV_H__
