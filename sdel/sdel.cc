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
  int *slv_id,ns=0;
  if(argc<2) io.msg(MFBD_ERROR|MFBD_FATAL,"slave ID required\n");
  for(int a=1;a<=argc-1;++a){
    int id;
    if(!sscanf(argv[a],"%d",&id)) io.msg(MFBD_ERROR|MFBD_FATAL,"could not convert %s to slave ID\n",argv[1]);
    int *tmp=new int [ns+1]-1;
    if(ns){
      memcpy(tmp+1,slv_id+1,ns*sizeof(int));
      delete[] (slv_id+1);
    }
    slv_id=tmp;
    slv_id[++ns]=id;
  }
//
  char *hostname;
  cmd.get("--master",hostname);
  int port;
  cmd.get("--port",port);
  sock_class sock(hostname,port);
// 
  sock.send_cmd(CMD_DEL_SLV); // tell manager what we want
//
  struct id sid;
  sock.send_id(sid);
  struct id mid(&sock);
  byte swap_endian=(mid.endian!=sid.endian);
  sock.set_swap(swap_endian);
//
  int size=(1+ns)*sizeof(int),offs=0;
  byte *buf=new byte [size];
  offs+=pack(buf+offs,ns,swap_endian);
  offs+=pack(buf+offs,slv_id,1,ns,swap_endian);
  sock.send(buf,size);        // send jid
  delete[] (slv_id+1);
  delete[] buf;
//
  buf=sock.recv(size);
  offs=0;
  offs+=unpack(buf+offs,ns,swap_endian);
  offs+=unpack(buf+offs,slv_id=new int [ns]-1,1,ns,swap_endian);
  for(int s=1;s<=ns;++s)
    if(slv_id[s]) io.msg(MFBD_ERROR,"no slave with ID %d found\n",slv_id[s]);
  delete[] buf;
  return 0;
}
