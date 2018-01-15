#ifndef __FDSET_H__ // __FDSET_H__
#define __FDSET_H__

class fdset;

#include <sys/types.h>
#include <poll.h>
#include "types.h"

class fdset{
  struct pollfd *fds;
  void **token;
  int idx;
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
