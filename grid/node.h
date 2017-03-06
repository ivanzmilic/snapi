#ifndef __NODE_H__  // __NODE_H__
#define __NODE_H__

#include <pthread.h>

#include "types.h"
#include "net.h"
#include "io.h"

//
// ********************************************
// * A node is a member of a network          *
// * It's task is to maintain the connections *
// * with it's neighbours and maintain the    *
// * basic communication mechanisms           *
// ********************************************
//

class node:public sock_class{
protected:
  io_class &io;
//
  int comm[2];
  pthread_t thread;
  pthread_mutex_t mutex;
//
  uint08_t ndim;  // grid dimensionality
  
  sock_class ****nghb_in; // incoming connections
  sock_class ****nghb_out; // outgoing connections
//
  uint08_t nw_init(void);   // network init
  uint08_t nw_clup(void);   // network cleanup
public:
//  node(const char**);
  node(uint08_t,io_class&);
  node(uint08_t*,int32_t&,uint08_t,io_class&);
  virtual ~node(void);
//
  virtual int32_t size(void);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);
//
  void *network_handler(void);
};

#endif              // __NODE_H__
