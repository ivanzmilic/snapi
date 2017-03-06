#ifndef __IO_H__
#define __IO_H__

#include <stdio.h>
#include <unistd.h>
#include "cmdln.h"
#include "net.h"
#include "types.h"

#define IOL_ERROR        0x00000001
#define IOL_WARN         0x00000002
#define IOL_INFO         0x00000004
#define IOL_XNFO         0x00000008
#define IOL_DEB1         0x00000010
#define IOL_DEB2         0x00000020
#define IOL_DEB3         0x00000040
#define IOL_LEVEL_MASK   0x000000FF

#define IOL_NOID         0x00000100
#define IOL_NOTM         0x00000200
#define IOL_FATAL        0x00000400

#define IOL_SCOPE_START  0x00010000
#define IOL_SCOPE_END    0x00020000

class io_class{
  FILE *f;
  sock_class *ssock;
  int sock;
  byte iocmd,swap_endian;
  int verb,level_mask,ts;
public:
  io_class(void);
  io_class(int,byte,byte,cmdln&);
  io_class(int,byte,byte,int,int);
  io_class(sock_class&,byte,byte,int,int);
  io_class(cmdln&);
  io_class(char*,char*&,const char*,int,int);
  io_class(char*,char*&,const char*,uid_t,gid_t,int,int);
  ~io_class(void);
  void reconf(int,byte,int,int);
  void reconf(int);
  int msg(int,const char*,...);
  int vlevel(void){ return verb; }
};

#endif
