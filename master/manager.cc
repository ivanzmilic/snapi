#include <stdio.h>
#include <string.h>

#include "array.h"
#include "net.h"
#include "job.h"
#include "slv.h"
#include "gslv.h"
#include "slvcfg.h"
#include "io.h"
#include "pack.h"
#include "fdset.h"
#include "cmdln.h"
#include "cmdcfg.h"
#include "version.h"

int main(int argc,char *argv[])
{
// setup recognised options (define defaults, then redefine with command line options)
  const char *options[]=CMDOPTS,*maps[]=CMDMAPS;    // recognised options
  char *hlp;
  cmdln cmd(argc,argv,options,maps);    // parse command line
  if(cmd.get("--help",hlp)){
    fprintf(stdout,"%s\n",(cmd.get(hlp))?cmd.help(hlp):cmd.help("--help"));
    return 0;
  }
  int version;
  cmd.get("--version",version);
  if(version){
    fprintf(stdout,"MOMFBD version: %s\n",VERSION_STR);
    return 0;
  }
  int port;
  cmd.get("-p",port);  
  int time_out;
  cmd.get("-to",time_out);  
  int clvl=3,quit=0;                       // compression level
  io_class io(cmd);
  int jid=0,sid=0;                         // job ID,slave ID
  job_class **queue=new job_class* [1];
  queue[0]=0;                              // queue is empty
  slv_class **islave=new slv_class* [1];
  islave[0]=0;                             // no inactive slaves
  slvcfg **cfg=new slvcfg* [1];
  cfg[0]=0;                                // no slave config data
  slv_class **aslave=new slv_class* [1];
  aslave[0]=0;                             // no active slaves
  sock_class ls(port);                     // open listening socket
  fdset fds;
  ls.add(fds);
//
  int gs;
  if(cmd.get("-g",gs)){
    fprintf(stdout,"Starting GRID slave\n");
    slv_class *slave=new gslv_class(++sid,io);
    slave->add(fds);
    array_add(slave,islave);
  }

  while(!quit){
    fds.act(1.00);
    for(void *nxt=fds.next();nxt;nxt=fds.next()){ // handle all events (if any)
      if(nxt==&ls){                               // new connection
        sock_class *ns=new sock_class(ls);        // accept connection: may be adding slaves
        switch(ns->recv_cmd()){
          case(CMD_NEW_SLV):{
            slv_class *slave=new slv_class(ns,++sid,clvl,io);
            slave->add(fds);
            array_add(slave,islave);
            break;
          }
          case(CMD_NEW_JOB):{
            job_class *job=new job_class(ns,++jid,cfg,io);
            if(job->status()>=0){                       // initialisation successful
              int pri=job->priority(),pos=0;
              while(queue[pos]->priority()<=pri) ++pos; // find the first job with higher priority number in the queue
              array_ins(job,queue,pos);                 // insert job in queue
            }else{                                      // job initialization failed!
              delete job;
              --jid;
            }
            delete ns;
            break;
          }
          case(CMD_DEL_SLV):{
            struct id slvid(ns);
            struct id mid;
            ns->send_id(mid);
//
            int size,offs=0,*ksid,nks;
            byte *data=ns->recv(size);
            offs+=unpack(data+offs,nks,0);
            offs+=unpack(data+offs,ksid=new int [nks]-1,1,nks,0);
//
            for(int s=0;aslave[s];++s)
              for(int k=1;k<=nks;++k)
                if(ksid[k])
                  if(aslave[s]->id()==ksid[k]){
                    aslave[s]->kill();
                    ksid[k]=0;
                  }
            for(int s=0;islave[s];++s)
              for(int k=1;k<=nks;++k)
                if(ksid[k])
                  if(islave[s]->id()==ksid[k]){
                    islave[s]->kill();
                    ksid[k]=0;
                  }
//
            offs=0;
            offs+=pack(data+offs,nks,0);
            offs+=pack(data+offs,ksid,1,nks,0);
            delete[] (ksid+1);
            ns->send(data,size);
            delete[] data;
            delete ns;
            break;
          }
          case(CMD_DEL_JOB):{
            struct id slvid(ns);
            struct id mid;
            ns->send_id(mid);
//
            int size,offs=0,*kjid,nkj;
            byte *data=ns->recv(size);
            offs+=unpack(data+offs,nkj,0);
            offs+=unpack(data+offs,kjid=new int [nkj]-1,1,nkj,0);
 //
            for(int q=0;queue[q];++q)
              for(int k=1;k<=nkj;++k)
                if(kjid[k])
                  if(queue[q]->jid()==kjid[k]){
                    queue[q]->kill();
                    kjid[k]=0;
                  }
//
            offs=0;
            offs+=pack(data+offs,nkj,0);
            offs+=pack(data+offs,kjid,1,nkj,0);
            delete[] (kjid+1);
            ns->send(data,size);
            delete[] data;
            delete ns;
            break;
          }
          case(CMD_MOD_JOB):{
            delete ns;
            break;
          }
          case(CMD_STAT):{
//            send_info
            switch(ns->recv_cmd()){
              case(CMD_STAT_JOBQ):{
                struct id slvid(ns);
                struct id mid;
                ns->send_id(mid);
//
                int size=0;
                byte *data=new byte [1];
                byte *tmp=new byte [10000]; // hopefully sufficient
                for(int q=0;queue[q];++q){
                  int sz=queue[q]->stat(tmp);
                  byte *t=new byte [size+sz];
                  if(size) memcpy(t,data,size);
                  delete[] data;
                  data=t;
                  if(sz) memcpy(t+size,tmp,sz);
                  size+=sz;
                }
                if(!ns->test()) ns->send(data,size);
                delete[] data;
                delete[] tmp;
                break;
              }
              case(CMD_STAT_SLVQ):{
                struct id slvid(ns);
                struct id mid;
                ns->send_id(mid);
//
                int size=0;
                byte *data=new byte [1];
                byte *tmp=new byte [10000]; // hopefully sufficient
                for(int q=0;aslave[q];++q){
                  int sz=aslave[q]->stat(tmp,'A');
                  byte *t=new byte [size+sz+1];
                  if(size) memcpy(t,data,size);
                  delete[] data;
                  data=t;
                  if(sz) memcpy(t+size,tmp,sz);
                  size+=sz;
                }
                for(int q=0;islave[q];++q){
                  int sz=islave[q]->stat(tmp,'W');
                  byte *t=new byte [size+sz+1];
                  if(size) memcpy(t,data,size);
                  delete[] data;
                  data=t;
                  if(sz) memcpy(t+size,tmp,sz);
                  size+=sz;
                }
                if(!ns->test()) ns->send(data,size);
                delete[] data;
                delete[] tmp;
                break;
              }
            }
//
            delete ns;
            break;
          }
          default: delete ns;                   // default: ignore silently
        }
      }else{                                    // established connection activity: results?
        slv_class *slave=(slv_class*)nxt;
        switch(slave->io_cmd()){
          case(CMD_QUIT):{
            array_del(slave,islave);            // remove slave from inactive list
            array_del(slave,aslave);            // remove slave from active list
            slave->del(fds);
            delete slave;                       // this will put back task, if any
	    break;
	  }
          case(CMD_THREAD_RES):{
	    slave->put_task();
            if(slave->is_busy()){
	      array_add(slave,islave);
              array_del(slave,aslave);          // remove slave from active list
            }
	    break;
	  }
          default: io.msg(IOL_ERROR,"ignoring unknown command from slave\n");
        }
      }
    } // for(fds.next)
    for(int s=0;aslave[s];++s) if(aslave[s]->idle_time()>time_out) aslave[s]->time_out();
    if(queue[0]->activate()==-2) io.msg(IOL_ERROR,"Could not create thread for job!\n"); // start things off
    for(int q=0;queue[q];++q){  // clean up the queue 
      switch(queue[q]->state()){
        case(5):
        case(7):{ // dead job
          delete queue[q];
          array_del(queue[q],queue);
          --q;
          break;
        }
      }
    }
    for(int q=0;queue[q];++q)
      if(queue[q]->state()==0){ // first inactive job in queue
        if(q) if(queue[q-1]->state()>=3) if(queue[q]->activate()==-2) io.msg(IOL_ERROR,"Could not create thread for job!\n"); // if previous job>=active : activate job 
        while(queue[q+1]) ++q;  // seek end of queue
      }
//
    for(int q=0;queue[q]&&islave[0];++q){   // hand out more work 
      while(islave[0]->get_task(queue[q]))
        if(islave[0]->is_busy()){
          array_add(islave[0],aslave);      // inactive -> active
          array_del(islave[0],islave);      // delete from inactive list
          if(!islave[0]) break;             // all slaves active
        } // while
    } // for(q)
  } // while(1)
  return 0;
}
