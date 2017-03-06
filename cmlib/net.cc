#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/param.h>
#include <sys/utsname.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "types.h"
#include "uts.h"
#include "fdset.h"
#include "pack.h"
#include "version.h"
#include "net.h"

id::id(void)
{
  varsize[VT_PTR]=sizeof(void*);
  varsize[VT_PID]=sizeof(pid_t);
  char hostname[MAXHOSTNAMELEN+1];
  gethostname(hostname,MAXHOSTNAMELEN+1);
  struct hostent *h=gethostbyname(hostname);
  name=strcpy(new char [strlen(h->h_name)+1],h->h_name);
  version=strcpy(new char [strlen(VERSION_STR)+1],VERSION_STR);
  struct utsname nm;
  if(uname(&nm)<0){
    sprintf(os=new char [2],"-");
    sprintf(machine=new char [2],"-");
  }else{
    sprintf(os=new char [strlen(nm.sysname)+strlen(nm.release)+2],"%s %s",nm.sysname,nm.release);
    sprintf(machine=new char [strlen(nm.machine)+1],"%s",nm.machine);
  }
  pid=getpid();
  int one=1;
  endian=(*(char*)&one==0);
}

id::id(sock_class *sock)
{
  if(sock->recv_id(*this)<0) memset(this,0,sizeof(struct id));
}

id::~id(void)
{
  delete[] name;
  delete[] version;
  delete[] os;
  delete[] machine;
}

struct sockaddr_in sock_data(const char *name,int port)
{
  struct sockaddr_in sd;
  if(strcmp(name,"any_host")){
    int IP=0;
    struct hostent *h=gethostbyname(name);
    if(!h) if((IP=inet_addr(name))==-1){
      fprintf(stderr,"Couldn't resolve host name\n");
      memset(&sd,0,sizeof(sd));
      return sd;
    }
    memset(&sd,'\0',sizeof(sd));
    if(IP) bcopy((const void *)&IP,(void *)&(sd.sin_addr),sizeof(IP));
      else bcopy(h->h_addr,(void *)&(sd.sin_addr),h->h_length);
  }else sd.sin_addr.s_addr=INADDR_ANY;
  sd.sin_family=AF_INET;
  sd.sin_port=htons(port);
  return sd;
}

sock_class::sock_class(int port_in)
{
  name=0;
  port=port_in;
  sd=sock_data((char*)"any_host",port);
  sds=sizeof(sd);
  sock=socket(AF_INET,SOCK_STREAM,0);
  int one=1;
  setsockopt(sock,SOL_SOCKET,SO_REUSEADDR,&one,sizeof(one));
  setsockopt(sock,SOL_SOCKET,SO_KEEPALIVE,&one,sizeof(one));
  if(bind(sock,(const struct sockaddr*)&(sd),sds)<0) exit(fprintf(stderr,"error: port %d already in use?\n",port));
  listen(sock,MAX_QUEUE);
  type=0;
  swap_endian=0;
}

void sock_class::outconn(const char *name_in,int port_in,fp_t timeout)
{
  port=port_in;
  name=strcpy(new char [strlen(name_in)+1],name_in);
  sd=sock_data(name,port);
  sds=sizeof(sd);
  int rv;
  do{
    sock=socket(AF_INET,SOCK_STREAM,0);
    if((rv=connect(sock,(const struct sockaddr *)&sd,sds))<0){
      fprintf(stderr,".");
      delay(timeout);
      close(sock);
    }
  }while(rv<0);
  type=1;
  swap_endian=0;
  fcntl(sock,F_SETFL,O_NONBLOCK);
}

sock_class::sock_class(const char *name_in,int port_in)
{
  outconn(name_in,port_in,30.0);
}

sock_class::sock_class(const char *name_in,int port_in,fp_t timeout_in)
{
  outconn(name_in,port_in,timeout_in);
}

sock_class::sock_class(sock_class &ls)
{
  sock=ls.accept(sd,sds);
  port=-1;
  name=0;
  type=1;
  swap_endian=0;
  fcntl(sock,F_SETFL,O_NONBLOCK);
}

sock_class::~sock_class(void)
{
  close(sock);
  if(name) delete[] name;
}

#ifndef HAVE_SOCKLEN_T
typedef int socklen_t;
#endif

int sock_class::accept(struct sockaddr_in &sd_in,size_t &sds_in)
{
  return ::accept(sock,(struct sockaddr *)&(sd_in),(socklen_t*)&(sds_in=sizeof(sd_in)));
}

int sock_class::test(void)
{
  fd_set fds;
  FD_ZERO(&fds);
  FD_SET(sock,&fds);
  struct timeval tv;
  tv.tv_sec=0;
  tv.tv_usec=0; 
  select(sock+1,&fds,NULL,NULL,&tv);
  if(FD_ISSET(sock,&fds)){                        // ready for read...
    char cmd;
    if(::recv(sock,&cmd,1,MSG_PEEK)<=0) return 1; // empty read or error!
    if(cmd==CMD_QUIT) return -1;                  // master wants us dead!
  }
  return 0;
}

byte sock_class::recv_cmd(void)
{
  byte cmd;
  fd_set rfds;
  int rv=0;
  while(!rv){
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    select(sock+1,&rfds,NULL,NULL,NULL);
    rv=::recv(sock,&cmd,1,0);
    if(!rv) if(errno!=EWOULDBLOCK) return CMD_CRASH; // crash
  }
  return cmd;
}

void sock_class::send_cmd(byte cmd)
{
  ::send(sock,&cmd,1,0);
}

byte *sock_class::recv(int &size)
{
  fd_set rfds;
  int offs=0,doffs;
  while(sizeof(size)>offs){
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    select(sock+1,&rfds,NULL,NULL,NULL);
    offs+=((doffs=::recv(sock,(byte*)(&size)+offs,sizeof(size)-offs,0))>0)?doffs:0;
    if(!doffs)
      if(errno!=EWOULDBLOCK) return 0; // slave crash
  }
  if(swap_endian) swap(&size,sizeof(size),1);
  if(!size) return 0;
  byte *buf=new byte [size];
  offs=0;
  while(size>offs){
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    select(sock+1,&rfds,NULL,NULL,NULL);
    offs+=((doffs=::recv(sock,buf+offs,size-offs,0))>0)?doffs:0;
    if(!doffs)
      if(errno!=EWOULDBLOCK){          // slave crash
        size=0;
        delete[] buf;
        return 0;
      }
  }
  return buf;
}

byte *sock_class::recv(int &size,int isock)
{
  fd_set rfds;
  int offs=0,doffs,maxsock=max(sock,isock);
  while(sizeof(size)>offs){
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    FD_SET(isock,&rfds);
    select(maxsock+1,&rfds,NULL,NULL,NULL);
    if(FD_ISSET(isock,&rfds)) return 0; // interruption
    offs+=((doffs=::recv(sock,(byte*)(&size)+offs,sizeof(size)-offs,0))>0)?doffs:0;
    if(!doffs)
      if(errno!=EWOULDBLOCK) return 0;    // slave crash
  }
  if(swap_endian) swap(&size,sizeof(size),1);
  if(!size) return 0;
  byte *buf=new byte [size];
  offs=0;
  while(size>offs){
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    FD_SET(isock,&rfds);
    select(maxsock+1,&rfds,NULL,NULL,NULL);
    if(FD_ISSET(isock,&rfds)){ // interruption
      size=0;
      delete[] buf;
      return 0;
    }
    offs+=((doffs=::recv(sock,buf+offs,size-offs,0))>0)?doffs:0;
    if(!doffs)
      if(errno!=EWOULDBLOCK){    // slave crash
        size=0;
        delete[] buf;
        return 0;
      }
  }
  return buf;
}

int sock_class::send(byte *data,int size)
{
  int sz=size;
  if(swap_endian) swap(&sz,sizeof(sz),1);
  ::send(sock,&sz,sizeof(sz),0);
  fd_set sfds,rfds;
  int offs=0,doffs;
  while(size>offs){
    FD_ZERO(&sfds);
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    FD_SET(sock,&sfds);
    select(sock+1,&rfds,&sfds,NULL,NULL);
    if(FD_ISSET(sock,&rfds)){
      byte cmd;
      int s=::recv(sock,&cmd,1,MSG_PEEK);
      if(!s) return -1;
    }else{
      offs+=((doffs=::send(sock,data+offs,size-offs,0))>0)?doffs:0;
      if(doffs<0) fprintf(stderr,"send error!\n");
    }
  }
  return 0;
}

int sock_class::send(byte *data,int size,int isock)
{
  int sz=size;
  if(swap_endian) swap(&sz,sizeof(sz),1);
  ::send(sock,&sz,sizeof(sz),0);
  fd_set sfds,rfds;
  int offs=0,doffs,maxsock=max(sock,isock);
  while(size>offs){
    FD_ZERO(&sfds);
    FD_ZERO(&rfds);
    FD_SET(sock,&sfds);
    FD_SET(sock,&rfds);
    FD_SET(isock,&rfds);
    select(maxsock+1,&rfds,&sfds,NULL,NULL);
    if(FD_ISSET(isock,&rfds)) return -1; // interruption
    if(FD_ISSET(sock,&rfds)){            // slave crash?
      byte cmd;
      int s=::recv(sock,&cmd,1,MSG_PEEK);
      if(!s) return -1;                  // only return on crash
    }else{                               // send!
      offs+=((doffs=::send(sock,data+offs,size-offs,0))>0)?doffs:0;
      if(doffs<0) fprintf(stderr,"send error!?\n");
    }
  }
  return 0;
}

int sock_class::send_id(struct id &id_in)
{
  int size=2+sizeof(pid_t)+strlen(id_in.name)+strlen(id_in.version)+strlen(id_in.os)+strlen(id_in.machine)+4;
  byte *data=new byte [size];
  int offs=0,doffs;
  offs+=pack(data+offs,id_in.varsize[VT_PTR]);
  offs+=pack(data+offs,id_in.varsize[VT_PID]);
  offs+=pack(data+offs,id_in.name);
  offs+=pack(data+offs,id_in.pid,0);
  offs+=pack(data+offs,id_in.version);
  offs+=pack(data+offs,id_in.os);
  offs+=pack(data+offs,id_in.machine);
//
  ::send(sock,&(id_in.endian),1,0);
  ::send(sock,&size,sizeof(size),0);
  fd_set sfds,rfds;
  offs=0;
  while(size>offs){
    FD_ZERO(&sfds);
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    FD_SET(sock,&sfds);
    select(sock+1,&rfds,&sfds,NULL,NULL);
    if(FD_ISSET(sock,&rfds)){ // slave died! (already?)
      byte cmd;
      int s=::recv(sock,&cmd,1,MSG_PEEK);
      if(!s){ // slave died! (already?)
        delete[] data;
        return -1;
      }
    }else{
      offs+=((doffs=::send(sock,data+offs,size-offs,0))>0)?doffs:0;
      if(doffs<0){
        fprintf(stderr,"send error in sock_class::send_id!\n");
        return -1;
      }
    }
  }
  delete[] data;
  return 0;
}

int sock_class::recv_id(struct id &id_in)
{
  fd_set rfds;
  int size=0,one=1;
  byte endian=(*(char*)&one==0); // big endian=1, little endian=0
  FD_ZERO(&rfds);
  FD_SET(sock,&rfds);
  select(sock+1,&rfds,NULL,NULL,NULL);
  if(!::recv(sock,&(id_in.endian),1,0)) return -1;
  int offs=0,doffs;
  while(sizeof(size)>offs){
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    select(sock+1,&rfds,NULL,NULL,NULL);
    offs+=((doffs=::recv(sock,(byte*)(&size)+offs,sizeof(size)-offs,0))>0)?doffs:0;
    if(!doffs) return -1;
  }
  byte lswap_endian=(id_in.endian!=endian);
  if(lswap_endian) swap(&size,sizeof(size),1);
  byte *data=new byte [size];
  offs=0;
  while(size>offs){
    FD_ZERO(&rfds);
    FD_SET(sock,&rfds);
    select(sock+1,&rfds,NULL,NULL,NULL);
    offs+=((doffs=::recv(sock,data+offs,size-offs,0))>0)?doffs:0;
    if(!doffs){                   // aiiiii!!!
      fprintf(stderr,"recv error in sock_class::recv_id!\n");
      delete[] data;
      return -1;
    }
  }
//
  offs=0;
  offs+=unpack(data+offs,id_in.varsize[VT_PTR]);
  offs+=unpack(data+offs,id_in.varsize[VT_PID]);
  offs+=unpack(data+offs,id_in.name);
  offs+=unpack(data+offs,id_in.pid,lswap_endian); // only here...
  offs+=unpack(data+offs,id_in.version);
  if(offs<size){
    offs+=unpack(data+offs,id_in.os);
    offs+=unpack(data+offs,id_in.machine);
  }else{
    strcpy(id_in.os=new char [2],"-");
    strcpy(id_in.machine=new char [2],"-");
  }
  delete[] data;
  return 0;
}
