#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "types.h"
#include "fdset.h"

fdset::fdset(void)
{
  mfds=-1;
  idx=-1;
  FD_ZERO(&fds);
  fd=new int[1]-1;
  token=new void*[1]-1;
  n=0;
}

fdset::~fdset(void)
{
  delete[] (fd+1);
  delete[] (token+1);
}

void fdset::add(int fd_in,void *token_in)
{
  FD_SET(fd_in,&fds);
  if(fd_in>mfds) mfds=fd_in;
  for(int i=1;i<=n;++i) if(fd[i]==fd_in) return; // already registered
  int *tfd=new int [n+1]-1;
  void **tto=new void* [n+1]-1;
  memcpy(tfd+1,fd+1,n*sizeof(int));
  memcpy(tto+1,token+1,n*sizeof(void*));
  delete[] (fd+1);
  delete[] (token+1);
  fd=tfd;
  token=tto;
  ++n;
  fd[n]=fd_in;
  token[n]=token_in;
}

void fdset::act(void)
{
  rd_fds=fds;
  select(mfds+1,&rd_fds,NULL,NULL,NULL);
  idx=0;
}

void fdset::act(fp_t t)
{
  struct timeval tv;
  tv.tv_sec=(time_t)t;
  tv.tv_usec=(int)((t-(fp_t)tv.tv_sec)*1.0E+06); 
  rd_fds=fds;
  select(mfds+1,&rd_fds,NULL,NULL,&tv);
  idx=0;
}

void *fdset::next(void)
{
  if(idx<0) return 0;
  while(idx<n){
    ++idx;
    if(FD_ISSET(fd[idx],&rd_fds)) return token[idx];
  }
  idx=-1;
  return 0;
}

void fdset::rem(int fd_in)
{
  FD_CLR(fd_in,&fds);
  for(int i=1;i<=n;++i)
    if(fd[i]==fd_in)
      if(n>1){
        int *tfd=new int [n-1]-1;
        void **tto=new void* [n-1]-1;
        if(i>1){
          memcpy(tfd+1,fd+1,(i-1)*sizeof(int));
          memcpy(tto+1,token+1,(i-1)*sizeof(void*));
        }
        if(n>i){
          memcpy(tfd+i,fd+i+1,(n-i)*sizeof(int));
          memcpy(tto+i,token+i+1,(n-i)*sizeof(void*));
        }
        delete[] (fd+1);
        delete[] (token+1);
        fd=tfd;
        token=tto;
        --n;
        return;
      }else --n;
  if(fd_in>=mfds)
    while((!FD_ISSET(mfds,&fds))&&(mfds>0)) --mfds;
}

