#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "version.h"
#include "types.h"
#include "cmdcfg.h"
#include "cmdln.h"
#include "net.h"
#include "io.h"
#include "pack.h"

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
    fprintf(stdout,"MOMFBD version: %s\n",VERSION_STR);
    exit(0);
  }
  io_class io(cmd);
//
  int *jid,nj=0;
  if(argc<2) io.msg(IOL_ERROR|IOL_FATAL,"job ID required\n");
  for(int a=1;a<=argc-1;++a){
    int id;
    if(!sscanf(argv[a],"%d",&id)) io.msg(IOL_ERROR|IOL_FATAL,"could not convert %s to job ID\n",argv[1]);
    int *tmp=new int [nj+1]-1;
    if(nj){
      memcpy(tmp+1,jid+1,nj*sizeof(int));
      delete[] (jid+1);
    }
    jid=tmp;
    jid[++nj]=id;
  }
//
  char *hostname;
  cmd.get("--master",hostname);
  int port;
  cmd.get("--port",port);
  sock_class sock(hostname,port);
// 
  sock.send_cmd(CMD_DEL_JOB); // get statistics
//
  struct id sid;
  sock.send_id(sid);
  struct id mid(&sock);
  byte swap_endian=(mid.endian!=sid.endian);
  sock.set_swap(swap_endian);
//
  int size=(1+nj)*sizeof(int),offs=0;
  byte *buf=new byte [size];
  offs+=pack(buf+offs,nj,swap_endian);
  offs+=pack(buf+offs,jid,1,nj,swap_endian);
  sock.send(buf,size);        // send jid
  delete[] (jid+1);
  delete[] buf;
//
  buf=sock.recv(size);
  offs=0;
  offs+=unpack(buf+offs,nj,swap_endian);
  offs+=unpack(buf+offs,jid=new int [nj]-1,1,nj,swap_endian);
  for(int j=1;j<=nj;++j)
    if(jid[j]) io.msg(IOL_ERROR,"no job with ID %d found\n",jid[j]);
  delete[] buf;
  return 0;
}
