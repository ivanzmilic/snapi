#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <poll.h>

#include "types.h"
#include "fdset.h"

fdset::fdset(void)
{
  idx=-1;
  fds=new struct pollfd [1]-1;
  token=new void*[1]-1;
  n=0;
}

fdset::~fdset(void)
{
  delete[] (fds+1);
  delete[] (token+1);
}

void fdset::add(int fd_in,void *token_in)
{
  for(int i=1;i<=n;++i) if(fds[i].fd==fd_in) return; // already registered
//
  struct pollfd *tfds=new struct pollfd [n+1]-1;
  void **tto=new void* [n+1]-1;
  memcpy(tfds+1,fds+1,n*sizeof(struct pollfd));
  memcpy(tto+1,token+1,n*sizeof(void*));
  delete[] (fds+1);
  delete[] (token+1);
  fds=tfds;
  token=tto;
  ++n;
  fds[n].fd=fd_in;
  fds[n].events=POLLIN;
  fds[n].revents=0; // very important for on the fly adding
  token[n]=token_in;
}

void fdset::act(void)
{
  for(int i=1;i<=n;++i) fds[i].revents=0;
  poll(fds+1,n,-1); // indefinite poll
  idx=0;
}

void fdset::act(fp_t t)
{
  for(int i=1;i<=n;++i) fds[i].revents=0;
  poll(fds+1,n,1000.0*t);
  idx=0;
}

void *fdset::next(void)
{
  if(idx<0) return 0;
  while(idx<n){
    ++idx;
    if(fds[idx].revents&POLLIN) return token[idx];
  }
  idx=-1;
  return 0;
}

void fdset::rem(int fd_in)
{
  for(int i=1;i<=n;++i)
    if(fds[i].fd==fd_in){
      if(n>1){
        struct pollfd *tfds=new struct pollfd [n-1]-1;
        void **tto=new void* [n-1]-1;
        if(i>1){
          memcpy(tfds+1,fds+1,(i-1)*sizeof(struct pollfd));
          memcpy(tto+1,token+1,(i-1)*sizeof(void*));
        }
        if(n>i){
          memcpy(tfds+i,fds+i+1,(n-i)*sizeof(struct pollfd));
          memcpy(tto+i,token+i+1,(n-i)*sizeof(void*));
        }
        delete[] (fds+1);
        delete[] (token+1);
        fds=tfds;
        token=tto;
        --n;
        return;
      }else --n;
    }
}

