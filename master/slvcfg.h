#ifndef __SLVCFG_H__ // __SLVCFG_H__
#define __SLVCFG_H__ 

#include <pthread.h>
#include "types.h"
#include "io.h"
#include "modes.h"
//#include "config.h"

struct slvcfg{
protected:
  int size;
  byte *data;
  int msize;
  byte *mdata;
  void **refs;
  slvcfg **&p;
  pthread_mutex_t active_lock;
public:
  slvcfg(byte *data_in,int size_in,slvcfg**&,io_class&);
  slvcfg(byte *data_in,int size_in,slvcfg**&);
  ~slvcfg(void);
  struct slvcfg *ref(void *r);
  void unref(void *r);
  struct slvcfg *same(byte *data_in,int size_in);
  int send(sock_class*,int);
//  struct modes *mde(int o,io_class &io);
//  struct config *cfg(io_class &io);
//  void modeinfo(io_class &io);
};

#endif               // __SLVCFG_H__
