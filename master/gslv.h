#ifndef __GSLV_H__ // __GSLV_H__
#define __GSLV_H__

#include <pthread.h>
#include <time.h>
#include "job.h"
#include "net.h"
#include "slvcfg.h"
#include "fdset.h"
#include "pack.h"

#include "slv.h"

struct tsk{
  job_class *job;      // job pointer
  struct chunk *chnk;  // job chunk data
  struct slvcfg *cfg;  // compressed slave config data for this chunk
//
  tsk(job_class *job_in,struct chunk *chnk_in){
    job=job_in;
    chnk=chnk_in;
    cfg=(chnk)?chnk->cfg->ref(this):0;
  }
  ~tsk(void){
    if(chnk&&job) job->unget_chunk(chnk); // were we doing work?: put back chunk for someone else
    if(cfg) cfg->unref(this);
  }
  struct slvcfg *ref(void){ return (chnk)?chnk->cfg->ref(this):0; }
  void unref(void){ chnk->cfg->unref(this); }
  void put(struct id &sid,fp_t &stime,fp_t &utime,int64_t &flop,int64_t &flops){
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
// return result
    job->put_chunk(chnk); // give processed data back to job
    chnk=0;
    job=0;
  }
};

class gslv_class:public slv_class{
  struct tsk **task;
public:
  gslv_class(int sid_in,io_class &io_in);
  virtual ~gslv_class(void);
  virtual int idle_time(void){ return 0; }; // don't time out (ever)
  virtual int get_task(job_class *job_in);
  virtual void put_task(void);
  virtual int stat(byte *data,const char rs){ return slv_class::stat(data,'R'); };
  virtual byte is_busy(void){ return 0; }; // never busy...
  void helper(void); // thread main loop
};

#endif             // __GSLV_H__
