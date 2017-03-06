#include <stdio.h>
#include <string.h>

#include "types.h"
#include "pack.h"
#include "mem.h"
#include "compress.h"
#include "io.h"
#include "jnfo.h"

jnfo::jnfo(byte *buf,byte swap_endian,io_class &io)
{
  memset(this,0,sizeof(jnfo));
// unpack the buffer
  int size;
  byte *data=z_uncompress(buf,size,swap_endian,io);
//
  int offs=0;
 //
  offs+=unpack(data+offs,na,swap_endian);
  if(na){
    atmos=new atmosphere* [na];
    for(int a=0;a<na;++a) atmos[a]=atmos_new(data,offs,swap_endian,io);
  }else atmos=0;
//
  offs+=unpack(data+offs,no,swap_endian);  
  if(no){
    offs+=unpack(data+offs,az=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,el=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,nlambda=new int [no],0,no-1,swap_endian);
    lambda=new fp_t* [no];
    name=new char* [no];
    for(int o=0;o<no;++o){
      offs+=unpack(data+offs,lambda[o]=new fp_t [nlambda[o]],0,nlambda[o]-1,swap_endian);
      offs+=unpack(data+offs,name[o]);
    }
 }
//
  offs+=unpack(data+offs,uid,swap_endian);
  offs+=unpack(data+offs,gid,swap_endian);
  offs+=unpack(data+offs,uname);
//
  offs+=unpack(data+offs,cdcl);
  offs+=unpack(data+offs,nmth);
  offs+=unpack(data+offs,nsth);
  offs+=unpack(data+offs,nsln);
//
  if(offs!=size) io.msg(IOL_WARN,"jnfo::jnfo: number of unpacked bytes was %d but %d bytes received (possible version conflict)?!\n",offs,size);
  delete[] data;
}

jnfo::jnfo(void)
{
  memset(this,0,sizeof(jnfo));
}

jnfo::~jnfo(void)
{
  if(na){
    for(int a=0;a<na;++a) delete atmos[a];
    delete[] atmos;
  }
  if(no){
    if(az) delete[] az;
    if(el) delete[] el;
    if(nlambda) delete[] nlambda;
    if(lambda) for(int o=0;o<no;++o) if(lambda[o]) delete[] lambda[o];
    if(lambda) delete[] lambda;
    if(name) for(int o=0;o<no;++o) delete[] name[o];
    if(name) delete[] name;
  }
  if(uname) delete[] uname;
  
}

byte *jnfo::compress(int &size,int level,byte swap_endian,io_class &io)
{
  int32_t usize=sizeof(int);                     // na
  for(int a=0;a<na;++a) usize+=atmos[a]->size(io);
//  io.msg(IOL_INFO,"jnfo::compress: usize %d\n",usize);
//
  usize+=sizeof(int);                            // no
  if(no){
    usize+=2*no*sizeof(fp_t);                    // az,el
    usize+=no*sizeof(int);                       // nlambda
    for(int o=0;o<no;++o){
      usize+=nlambda[o]*sizeof(fp_t);  // lambda
      usize+=strlen(name[o])+1;        // file name
    }
  }
  usize+=sizeof(uid_t)+sizeof(gid_t);            // uid,gid
  usize+=strlen(uname)+1;                        // uname
  usize+=4*sizeof(uint08_t);                     // cdcl,nmth,nsth,nsln
//
  byte *data=new byte [usize];
//
  int offs=pack(data,na,swap_endian);
  for(int a=0;a<na;++a) offs+=atmos[a]->pack(data+offs,swap_endian,io);
//  io.msg(IOL_INFO,"jnfo::compress: %d %d\n",offs,usize);
//
  offs+=pack(data+offs,no,swap_endian);
  if(no){
    offs+=pack(data+offs,az,0,no-1,swap_endian);
    offs+=pack(data+offs,el,0,no-1,swap_endian);
    offs+=pack(data+offs,nlambda,0,no-1,swap_endian);
    for(int o=0;o<no;++o){
      offs+=pack(data+offs,lambda[o],0,nlambda[o]-1,swap_endian);
      offs+=pack(data+offs,name[o]);
    }
  }
//
  offs+=pack(data+offs,uid,swap_endian);
  offs+=pack(data+offs,gid,swap_endian);
  offs+=pack(data+offs,uname);
//
  offs+=pack(data+offs,cdcl); // transfer compression level
  offs+=pack(data+offs,nmth); // # master threads
  offs+=pack(data+offs,nsth); // # slave threads
  offs+=pack(data+offs,nsln); // # slave nodes
//  io.msg(IOL_INFO,"jnfo::compress: %d %d\n",offs,usize);
//
  if(offs!=usize) io.msg(IOL_WARN,"jnfo::compress: number of packed bytes was %d but %d bytes allocated\n",offs,usize);
//
  size=usize;
  byte *buf=z_compress(data,size,level,swap_endian,io);
  delete[] data;
  return buf;
}
