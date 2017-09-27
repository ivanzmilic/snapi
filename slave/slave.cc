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
#include "compress.h"
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
  class model *mod=0;
  class observable *obs=0;
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
          int64_t flops=0;
          struct tms t_start;
          clock_t t0=times(&t_start);

          int32_t csz;
          byte *buf=sock.recv(csz);

          int32_t size;
          uint08_t *ubuf=z_uncompress(buf,size,0,io); // decompress results
          delete[] buf;
          int32_t offs=0;
          atmos=atmos_new(ubuf,offs,swap_endian,io);  // create atmospheric structure
          obs=obs_new(ubuf,offs,swap_endian,io);  // create observbable
          int to_invert = obs->get_to_invert();
          if (to_invert)
            mod=model_new(ubuf,offs,swap_endian,io);  // create model
          delete[] ubuf;
          if(offs!=size) io.msg(IOL_ERROR,"inaccurate buffer size estimate! (actual: %d > estimate: %d)\n",offs,size);

          fp_t * lambda = obs->get_lambda();
          int nlambda = obs->get_n_lambda();
          fp_t el = obs->get_el();
          fp_t az = obs->get_az();

          uint08_t *data;
          int32_t rsz;

          if (to_invert){

            class observable *fit=atmos->stokes_lm_fit(obs,el,az,mod);
            rsz=atmos->size(io);
            rsz+=mod->size(io);
            rsz+=fit->size(io);
            rsz+=3*sizeof(int32_t);
            data=new uint08_t [rsz];
            offs=atmos->pack(data,0,io);
            offs+=fit->pack(data+offs,0,io);
            offs+=mod->pack(data+offs,0,io);
          
            // I prefer deleting immediately after. 
            delete fit;
            delete mod;
            delete atmos;
            delete obs;
            delete[](lambda+1);
          }
          else {

            atmos->set_grid(0);
            class observable *synth=atmos->obs_stokes(el,az,lambda,nlambda);
            rsz=atmos->size(io);
            rsz+=synth->size(io);
            rsz+=3*sizeof(int32_t);
            data=new uint08_t [rsz];
            offs=atmos->pack(data,0,io);
            offs+=synth->pack(data+offs,0,io);
          
            // I prefer deleting immediately after. 
            delete synth;
            delete atmos;
            delete obs;
            delete[](lambda+1);
          }
          
//
          struct tms t_end;
          clock_t t1=times(&t_end);
//
          int32_t user=t_end.tms_utime-t_start.tms_utime;
          int32_t sys=t_end.tms_stime-t_start.tms_stime;
          int32_t clock=tickspersec;
          offs+=pack(data+offs,user,swap_endian);
          offs+=pack(data+offs,sys,swap_endian);
          offs+=pack(data+offs,clock,swap_endian);
          if(offs!=rsz) io.msg(IOL_WARN,"in main: packed %d bytes, but buffer is %d bytes!\n",offs,size); // sanity check...
          size=rsz;
          buf=z_compress(data,size,clvl,swap_endian,io);       // +(2)
          delete[] data;                                       // -(1)
//
          data=new byte [size+2*sizeof(fp_t)+sizeof(int64_t)]; // +(3)
          memcpy(data,buf,size);
          delete[] buf;                                        // -(2)
          offs=size;
          fp_t utime=(fp_t)user/(fp_t)clock;
          fp_t stime=(fp_t)sys/(fp_t)clock;
          offs+=pack(data+offs,utime,swap_endian);
          offs+=pack(data+offs,stime,swap_endian);
          offs+=pack(data+offs,flops,swap_endian);
          size=offs;
// send data to master
          int conn_stat=sock.test();
          if(!conn_stat){
            sock.send_cmd(CMD_SLV_RES);
            sock.send(data,size);
          }else{
            alive=0;
            quit=(conn_stat<0);
          }
          delete[] data;                                       // -(3)
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
