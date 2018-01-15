#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/param.h>
#include <sys/utsname.h>
#include <sys/time.h>
#include <poll.h>
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
  if(name) delete[] name;
  if(version) delete[] version;
  if(os) delete[] os;
  if(machine) delete[] machine;
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
  struct pollfd fds;
  fds.fd=sock;
  fds.events=POLLIN;
  fds.revents=0;
//
  poll(&fds,1,0);
  if(fds.revents&POLLIN){                        // ready for read...
    char cmd;
    if(::recv(sock,&cmd,1,MSG_PEEK)<=0) return 1; // empty read or error!
    if(cmd==CMD_QUIT) return -1;                  // master wants us dead!
  }
  return 0;
}

byte sock_class::recv_cmd(void)
{
  struct pollfd fds;
  fds.fd=sock;
  fds.events=POLLIN;
//
  byte cmd;
  int rv=0;
  while(!rv){
    fds.revents=0;
    poll(&fds,1,-1);
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
  struct pollfd rfds;
  rfds.fd=sock;
  rfds.events=POLLIN;
//
  int offs=0,doffs;
  while(sizeof(size)>(uint32_t)offs){
    rfds.revents=0;
    poll(&rfds,1,-1);
    offs+=((doffs=::recv(sock,(byte*)(&size)+offs,sizeof(size)-offs,0))>0)?doffs:0;
    if(!doffs)
      if(errno!=EWOULDBLOCK) return 0; // slave crash
  }
  if(swap_endian) swap(&size,sizeof(size),1);
  if(!size) return 0;
  byte *buf=new byte [size];
  offs=0;
  while(size>offs){
    rfds.revents=0;
    poll(&rfds,1,-1);
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
  struct pollfd fds[2];
  fds[0].fd=sock;
  fds[1].fd=isock;
  fds[0].events=fds[1].events=POLLIN;
//
  int offs=0,doffs;
  while(sizeof(size)>(uint32_t)offs){
    fds[0].revents=fds[1].revents=0;
    poll(fds,2,-1);
    if(fds[1].revents&POLLIN) return 0; // interruption
    offs+=((doffs=::recv(sock,(byte*)(&size)+offs,sizeof(size)-offs,0))>0)?doffs:0;
    if(!doffs)
      if(errno!=EWOULDBLOCK) return 0;    // slave crash
  }
  if(swap_endian) swap(&size,sizeof(size),1);
  if(!size) return 0;
  byte *buf=new byte [size];
  offs=0;
  while(size>offs){
    fds[0].revents=fds[1].revents=0;
    poll(fds,2,-1);
    if(fds[1].revents&POLLIN){ // interruption
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
  ::send(sock,&sz,sizeof(sz),MSG_NOSIGNAL);
  if(errno==EPIPE) fprintf(stderr,"sock_class::send: remote connection shutdown\n");
//
  struct pollfd fds;
  fds.fd=sock;
  fds.events=POLLIN|POLLOUT;
//
  int offs=0,doffs;
  while(size>offs){
    fds.revents=0;
    poll(&fds,1,-1);
    if(fds.revents&POLLIN){ // something to read...
      byte cmd;
      int s=::recv(sock,&cmd,1,MSG_PEEK);
      if(!s) return -1;
    }else{
      offs+=((doffs=::send(sock,data+offs,size-offs,MSG_NOSIGNAL))>0)?doffs:0;
      if((doffs<0)||(errno==EPIPE)) fprintf(stderr,"sock_class::send: send error!?\n");
    }
  }
  return 0;
}

int sock_class::send(byte *data,int size,int isock)
{
  int sz=size;
  if(swap_endian) swap(&sz,sizeof(sz),1);
  ::send(sock,&sz,sizeof(sz),MSG_NOSIGNAL);
  if(errno==EPIPE) fprintf(stderr,"sock_class::send: remote connection shutdown\n");
//
  struct pollfd fds[2];
  fds[0].fd=sock;
  fds[0].events=POLLIN|POLLOUT;
  fds[1].fd=isock;
  fds[1].events=POLLIN;
//
  int offs=0,doffs;
  while(size>offs){
    fds[0].revents=fds[1].revents=0;
    poll(fds,2,-1);
    if(fds[1].revents&POLLIN) return -1; // interruption
    if(fds[0].revents&POLLIN){           // slave crash?
      byte cmd;
      int s=::recv(sock,&cmd,1,MSG_PEEK);
      if(!s) return -1;                  // only return on crash
    }else{          // send!
      offs+=((doffs=::send(sock,data+offs,size-offs,MSG_NOSIGNAL))>0)?doffs:0;
      if((doffs<0)||(errno==EPIPE)) fprintf(stderr,"sock_class::send: send error!?\n");
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
  ::send(sock,&(id_in.endian),1,MSG_NOSIGNAL);
  ::send(sock,&size,sizeof(size),MSG_NOSIGNAL);
  if(errno==EPIPE){
    fprintf(stderr,"send error in sock_class::send_id!\n");
    delete[] data;
    return -1;
  }
//
  struct pollfd fds;
  fds.fd=sock;
  fds.events=POLLIN|POLLOUT;
//
  offs=0;
  while(size>offs){
    fds.revents=0;
    poll(&fds,1,-1);
    if(fds.revents&POLLIN){ // slave data...
      byte cmd;
      int s=::recv(sock,&cmd,1,MSG_PEEK);
      if(!s){ // slave died! (already?)
        delete[] data;
        return -1;
      }
    }else{
      offs+=((doffs=::send(sock,data+offs,size-offs,MSG_NOSIGNAL))>0)?doffs:0;
      if((doffs<0)||(errno==EPIPE)){
        fprintf(stderr,"send error in sock_class::send_id!\n");
        delete[] data;
        return -1;
      }
    }
  }
  delete[] data;
  return 0;
}

int sock_class::recv_id(struct id &id_in)
{
  struct pollfd rfds;
  rfds.fd=sock;
  rfds.events=POLLIN;

  int size=0,one=1;
  byte endian=(*(char*)&one==0); // big endian=1, little endian=0
  rfds.revents=0;
  poll(&rfds,1,-1);
  if(!::recv(sock,&(id_in.endian),1,0)) return -1;
  int offs=0,doffs;
  while(sizeof(size)>(uint32_t)offs){
    rfds.revents=0;
    poll(&rfds,1,-1);
    offs+=((doffs=::recv(sock,(byte*)(&size)+offs,sizeof(size)-offs,0))>0)?doffs:0;
    if(!doffs) return -1;
  }
  byte lswap_endian=(id_in.endian!=endian);
  if(lswap_endian) swap(&size,sizeof(size),1);
  byte *data=new byte [size];
  offs=0;
  while(size>offs){
    rfds.revents=0;
    poll(&rfds,1,-1);
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
