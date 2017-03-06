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
#include "gslv.h"

void *gslv_thread_start(void *p)
{
  gslv_class *slv=(gslv_class*)p;
  slv->helper();
  return 0;
}

gslv_class::gslv_class(int id_in,io_class &io_in):slv_class(id_in,gslv_thread_start,io)
{
  if(sid.machine) delete[] sid.machine;
  sid.machine=strcpy(new char [strlen("GRID")+1],"GRID");
  if(sid.name) delete[] sid.name;
  sid.name=strcpy(new char [strlen("GRID client")+1],"GRID client");
//
  task=new tsk *[1];
  task[0]=0;
}

gslv_class::~gslv_class(void)
{
  for(int t=0;task[t];++t){ // should be empty....
    // cancel GRID job
    delete task[t]; // put back and delete task
  }
  delete[] task;
}

int gslv_class::get_task(job_class *job_in)
{
  if(struct chunk *c=job_in->get_chunk(this)){
    struct tsk *t=new tsk(job_in,c);
    send_cmd(mgr_sock,CMD_THREAD_TSK);
    ::send(mgr_sock,&t,sizeof(struct tsk*),0);
    return 1;
  }
  return 0;
}

void gslv_class::put_task(void)
{ 
  struct tsk *t;
  ::recv(mgr_sock,&t,sizeof(struct tsk*),0);
  t->put(sid,stime,utime,flop,flops);
  delete t;
}

void gslv_class::helper(void) // thread main loop
{
  fdset fds;
  fds.add(thr_sock,&thr_sock);
/*
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
	  case(CMD_THREAD_TSK):{
            reset_time();
            sock->send_cmd(CMD_SLV_CFG);
            if(sock->send(cfg->data,cfg->size,thr_sock)<0){
              quit=1;
              while(token=fds.next());   // empty the queue
	    }
            sock->send_cmd(CMD_SLV_TSK);
            byte *buf;
            int bsz;
            chnk->get(buf,bsz);
            if(sock->send(buf,bsz,thr_sock)<0){
              quit=1;
              while(token=fds.next());   // empty the queue
	    }
            chnk->free();
            break;
	  }
	  default: io.msg(MFBD_WARN,"ignoring unknown command from manager %d\n",cmd);
	}
      }
      if(token==sock){ // message from GRID API
	switch(sock->recv_cmd()){
          case(CMD_QUIT):
          case(CMD_CRASH):{  // slave died!
            quit=1;
            while(token=fds.next());   // empty the queue
	    break;
          }
	  case(CMD_SLV_IO):{
            int size,offs=0,type;
            char *mesg;
            byte *data=sock->recv(size);
            offs+=unpack(data+offs,type,0);
            offs+=unpack(data+offs,mesg);
            delete[] data;
            if(job) job->msg(type,mesg); else io.msg(type,mesg);
            delete[] mesg;
            break;
          }
          case(CMD_SLV_RES):{ // receive data
	    int size;
	    byte *buf=sock->recv(size,thr_sock);
            if(!buf){         // active slave (connection) died! image restoration was not done...
              quit=1;
            }else{ // result is here: give it to the corresponding task
              chnk->result(buf,size);
              send_cmd(thr_sock,CMD_THREAD_RES);
              ::send(thr_sock,&t,sizeof(struct tsk*),0);
	    }
	    break;
	  }
          case(CMD_SLV_REJ):{ // slave rejected the data, don't try again
            chnk->result(0,0);
            send_cmd(thr_sock,CMD_THREAD_RES);
	    break;
	  }
	  case(CMD_SLV_CFG):{
            sock->send_cmd(CMD_SLV_TSK);
            byte *buf;
            int bsz;
            chnk->get(buf,bsz);
            if(sock->send(buf,bsz,thr_sock)<0){
              quit=1;
              while(token=fds.next());   // empty the queue
	    }
            chnk->free();
            break;
	  }
          default: job->msg(MFBD_WARN,"ignoring unknown command\n");
        }
      }
    }
  }while(!quit);
  sock->del(fds);
  fds.rem(thr_sock);
  send_cmd(thr_sock,CMD_QUIT);
*/
}

