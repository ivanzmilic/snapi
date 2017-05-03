#include <pthread.h>
#include <string.h>
#include <errno.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "pack.h"

#include "node.h"

void *start_network_handler(void *arg)
{
  return ((node*)arg)->network_handler();
}

node::node(uint08_t ndim_in,io_class &io_in):sock_class(0),io(io_in)
{
  if((ndim=ndim_in)) if(nw_init()) io.msg(IOL_ERROR,"node::node: failed to initialize network handler\n");
}

node::node(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):sock_class(0),io(io_in)
{
// extract ndim?
  offs+=unpack(buf+offs,do_swap,io); 
  if(ndim) if(nw_init()) io.msg(IOL_ERROR,"node::node: failed to initialize network handler\n");
}

node::~node(void)
{
  if(ndim) if(nw_clup()) io.msg(IOL_ERROR,"node::~node: failed to shutdown network handler\n");
}

int32_t node::size(void)
{
  return sizeof(uint08_t);
}

int32_t node::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{

  int offs=::pack(buf,ndim);
  io_in.msg(IOL_DEB1,"node::pack: ndim=%d\n",ndim);
  return offs;
}

int32_t node::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
  int offs=::unpack(buf,ndim);
  io_in.msg(IOL_DEB1,"node::unpack: ndim=%d\n",ndim);
  return offs;
}

uint08_t node::nw_init(void) // init
{
  nghb_in=(sock_class****)v3dim(-1,1,-1,1,-1,1);
  memset(nghb_in[-1][-1]-1,0,9*sizeof(void*));
  nghb_in[0][0][0]=this;
//
  nghb_out=(sock_class****)v3dim(-1,1,-1,1,-1,1);
  memset(nghb_out[-1][-1]-1,0,9*sizeof(void*));
  nghb_out[0][0][0]=this;
//
  pthread_mutex_init(&mutex,0);
  socketpair(AF_UNIX,SOCK_STREAM,0,comm);
  if(pthread_create(&thread,0,start_network_handler,this)<0)
    io.msg(IOL_ERROR,"node::nw_init: failed to create communication main thread: %s!\n",strerror(errno));
  else
    io.msg(IOL_INFO,"node::nw_init: succesfully initialized communication interface main thread\n");
// sending is taken care of by the send_thread to avoid blocking conditions
  
  return 0;
}

uint08_t node::nw_clup(void) // cleanup...
{
  void *rv;
  pthread_join(thread,&rv);
  pthread_join(thread,&rv);
// threads are gone...
  pthread_mutex_destroy(&mutex);
  if(nghb_in){
    for(int x=-1;x<=1;++x)
      for(int y=-1;y<=1;++y)
        for(int z=-1;z<=1;++z)
          if(nghb_in[x][y][z]!=this) delete nghb_in[x][y][z];
    del_v3dim((void****)nghb_in,-1,1,-1,1,-1,1);
  }
//
  if(nghb_out){
    for(int x=-1;x<=1;++x)
      for(int y=-1;y<=1;++y)
        for(int z=-1;z<=1;++z)
          if(nghb_out[x][y][z]!=this) delete nghb_out[x][y][z];
    del_v3dim((void****)nghb_out,-1,1,-1,1,-1,1);
  }
  return 0;
}

void *node::network_handler(void)
{
  fdset fds;
  fds.add(comm[1],this);
//
  for(int x=-1;x<=1;++x)
    for(int y=-1;y<=1;++y)
      for(int z=-1;z<=1;++z)
        if(nghb_out[x][y][z]) nghb_out[x][y][z]->add(fds);
//
  while(1){
    fds.act();
    for(void *nxt=fds.next();nxt;nxt=fds.next()){
      if(nxt==this){  // from main thread
        io.msg(IOL_INFO,"node::network_handler: got data from main...\n");
//        ::recv(comm[1],&cif,sizeof(struct com*),0);
      
      }else{          // from neighbour
        io.msg(IOL_INFO,"node::network_handler: got data from neighbour...\n");
//        socket->recv(cif->data,cif->sz);

      }
    }
  }
//
  return 0;
}

