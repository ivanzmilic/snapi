#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "types.h"
#include "pack.h"
#include "io.h"
#include "uts.h"
#include "compress.h"

byte *z_compress(byte *data,int &size,int clvl,byte swap_endian,io_class &io)
{
  int32_t usize=size,csize;
  uLongf dl=(uLongf)(((int64_t)101*(int64_t)usize)/(int64_t)100+(int64_t)12),sl=usize;
  Bytef *dest=new Bytef [dl],*source=data;
  int succes=0; // assume failure
  switch(compress2(dest,&dl,source,sl,clvl)){
    case(Z_OK): succes=1; break;
    case(Z_MEM_ERROR): io.msg(IOL_ERROR,"compressing data: out of memory\n"); break;
    case(Z_BUF_ERROR): io.msg(IOL_ERROR,"compressing data: buffer not large enough\n"); break;
    default: io.msg(IOL_ERROR,"compressing data: unknown reason\n"); break;
  }
  if(succes){
    byte *buf=new byte [(csize=dl)+2*sizeof(int32_t)];
    int offs=0;
    offs+=pack(buf+offs,csize,swap_endian);
    offs+=pack(buf+offs,usize,swap_endian);
    memcpy(buf+offs,dest,csize);
    delete[] dest;
    size=csize+2*sizeof(int32_t);
    return buf;
  }
  delete[] dest;
  size=0;
  return 0;
}

byte *z_uncompress(byte *buf,int &size,byte swap_endian,io_class &io)
{
  int offs=0;
  int32_t csize,usize;
  offs+=unpack(buf+offs,csize,swap_endian);
  offs+=unpack(buf+offs,usize,swap_endian);
  byte *data=new byte [usize];             // config data block
  Bytef *source=buf+offs,*dest=data;
  uLong dl=usize,sl=csize;
  int succes=0;                           // assume failure
  switch(uncompress(dest,&dl,source,sl)){
    case(Z_OK): succes=1; break;
    case(Z_DATA_ERROR): io.msg(IOL_ERROR,"decompressing data: corrupt buffer content\n"); break;
    case(Z_MEM_ERROR): io.msg(IOL_ERROR,"decompressing data: out of memory\n"); break;
    case(Z_BUF_ERROR): io.msg(IOL_ERROR,"decompressing data: buffer not large enough\n"); break;
    default: io.msg(IOL_ERROR,"decompressing data: unknown reason\n"); break;
  }
  if(succes){
    size=dl;
    return data;
  }
  delete[] data;
  size=0;
  return 0;
}
