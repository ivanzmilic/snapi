#include <sys/times.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "types.h"
#include "nr/fft.h"
//#include "config.h"
#include "cmdln.h"
#include "mem.h"
#include "pack.h"
#include "uts.h"
#include "net.h"
#include "conf.h"
#include "version.h"

#include "cmdcfg.h"
#include "atmos.h"

int main(int argc,char *argv[])
{
// setup recognised options (define defaults, then redefine with command line options)
  const char *options[]=CMDOPTS,*maps[]=CMDMAPS;    // recognised options
  char *hlp;
  cmdln cmd(argc,argv,options,maps);    // parse command line
  if(cmd.get("--help",hlp)){
    fprintf(stdout,"%s\n",(cmd.get(hlp))?cmd.help(hlp):cmd.help("--help"));
    exit(0);
  }
  int version;
  cmd.get("--version",version);
  if(version){
    fprintf(stdout,"INVERT version: %s\n",VERSION_STR);
    exit(0);
  }
  int tickspersec=sysconf(_SC_CLK_TCK);

  int clvl=3,quit=0;
  char *host;
  cmd.get("--master",host);
  int port;
  cmd.get("--port",port);  
  int nthread;
  cmd.get("--nthread",nthread);
  fp_t timeout;
  cmd.get("--timeout",timeout);
//  fprintf(stderr,"computing on %d threads\n",nthread);
  class atmosphere *atmos=0;
//
  while(!quit){
    sock_class sock(host,port,timeout);         // open connection
    sock.send_cmd(CMD_NEW_SLV);                 // slave_connect
    struct id sid;
    sock.send_id(sid);                          // identify
    struct id mid(&sock);                       // master ID

    byte swap_endian=(sid.endian!=mid.endian);  // do we need to byteswap?
    sock.set_swap(swap_endian);

    io_class io(sock,swap_endian,CMD_SLV_IO,2,0); // I/O to master, default verbosity level = 2

    int alive=1;

/*
    slvthrd **thread=0;
    if(nthread>1){
      thread=new slvthrd* [nthread+1];
      for(int s=0;s<nthread;++s) thread[s]=new slvthrd(0,io);
      thread[nthread]=0;
    }
*/

    while(alive){
      switch(uint08_t msg=sock.recv_cmd()){
        case(CMD_CRASH):{    // master/connection died accidentally
          alive=0;
          break;
        }
        case(CMD_QUIT):{     // master sent quit
          alive=0;
          quit=1;
          break;
        }
        case(CMD_SLV_CFG):{  // new config
          int32_t size,offs=0;
          uint08_t *buf=sock.recv(size);
          if(atmos) delete atmos;         // cleanup old structure
          atmos=new atmosphere(buf,offs,swap_endian,io);  // create atmospheric structure
          delete[] buf;
          break;
        }
        case(CMD_SLV_TSK):{  // new patch
          break;
        }
        default:{
          io.msg(IOL_ERROR,"slave: unknown command 0x%X\n",msg);
        }
      }
    }
/*
    if(thread){
      for(int s=0;s<nthread;++s) delete thread[s];
      delete[] thread;
    }
*/
    quit=1;
  }
  return 0;
}
