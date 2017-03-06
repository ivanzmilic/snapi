#ifndef __NET_H__ // __NET_H__
#define __NET_H__

#define VT_PTR 0
#define VT_PID 1
#define NUM_VT 2

#define CMD_QUIT       0
#define CMD_NEW_SLV    1
#define CMD_NEW_JOB    2
#define CMD_MOD_JOB    3
#define CMD_DEL_JOB    4
#define CMD_STAT       5
#define CMD_SLV_CFG    6
#define CMD_SLV_TSK    7
#define CMD_SLV_IO     8
#define CMD_SLV_RES    9
#define CMD_SLV_REJ   10
#define CMD_DEL_SLV   11
#define CMD_CRASH    255

#define CMD_STAT_JOBQ    1
#define CMD_STAT_SLVQ    2

#define MAX_QUEUE                  255

#include <netinet/in.h>
#include "types.h"
#include "fdset.h"

class sock_class;

struct id{
  byte varsize[NUM_VT];
  byte endian;
  pid_t pid;
  char *name,*version,*os,*machine;
//
  id(void);
  id(sock_class*);
  ~id(void);
  int status(void){
    if(name||version||os||machine||pid) return 0;
    return -1;
  }
};

class sock_class{
  struct sockaddr_in sd;
  size_t sds;
  int sock,port;
  char *name;
  byte type,swap_endian;
  int accept(struct sockaddr_in&,size_t&);
  void outconn(const char*,int,fp_t);
public:
  sock_class(int);             // listening connection
  sock_class(const char*,int); // outgoing connection
  sock_class(const char*,int,fp_t); // outgoing connection with connection timeout
  sock_class(sock_class&);     // accept connection
  ~sock_class(void);
  void add(fdset &fds){
    fds.add(sock,(void*)this);
  }
  void add(fdset &fds,void *token){
    fds.add(sock,token);
  }
  void del(fdset &fds){
    fds.rem(sock);
  }
  int test(void);
  void send_cmd(byte);
  byte recv_cmd(void);
  int send(byte*,int);
  int send(byte*,int,int);
  byte *recv(int&);
  byte *recv(int&,int);
  int send_id(struct id&);
  int recv_id(struct id&);
  void set_swap(byte swap_endian_in){ swap_endian=swap_endian_in; }
};

#endif            // __NET_H__
