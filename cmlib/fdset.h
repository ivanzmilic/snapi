#ifndef __FDSET_H__ // __FDSET_H__
#define __FDSET_H__

class fdset;

#include <sys/types.h>
#include "types.h"

class fdset{
  fd_set fds,rd_fds;
  int mfds,idx;
  int *fd;
  void **token;
  int n;
public:
  fdset(void);
  ~fdset(void);
  void rem(int);
  void add(int,void*);
  void act(void);
  void act(fp_t);
  void *next(void);
};

#endif              // __FDSET_H__
