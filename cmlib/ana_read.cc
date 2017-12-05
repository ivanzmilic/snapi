#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include "types.h"
#include "io.h"
#include "anadecompress.h"

#include "uts.h"
#include "fileio.h"
#include "ana_io.h"

int ck_synch_hd(int fin,struct fzhead *&fho,int t_endian,io_class &io)
{
  int wwflag=0;
  struct fzhead fh;
  if(read(fin,&fh,sizeof(struct fzhead))!=sizeof(struct fzhead)) io.msg(IOL_ERROR,"error in fzread while reading header\n");
//
  int syncpat=(fh.synch_pattern==0x5555aaaa);
  int revsyncpat=(fh.synch_pattern==0xaaaa5555);
  if(!(syncpat||revsyncpat)){
    close(fin);
    fho=0;
    io.msg(IOL_ERROR,"error: file does not have the F0 synch pattern (found 0x%x instead)\n",fh.synch_pattern);
    return -1;
  }
//
  fho=new fzhead [fh.nhb];
  fho[0]=fh;
//
  if(syncpat==t_endian){
    io.msg(IOL_WARN,"warning: reversed F0 synch pattern\n");
    wwflag=1; 
  }
  if(fh.nhb>1){                     // if the header is long, read in the rest now
    if(fh.nhb>15){
      close(fin);
      io.msg(IOL_ERROR,"error: annot handle header more than 16 blocks!\n");
      return -1;
    }
    int size=(fh.nhb-1)*sizeof(struct fzhead);
    if(read(fin,fho+1,size)!=size){
      close(fin);
      io.msg(IOL_ERROR,"error reading file: %s\n",strerror(errno));
      return -1;
    }
  }
  if(t_endian) swap(fho->dim,fho->ndim); // for big endian machines
  return wwflag; 
}

int ana_data_type(const char *file_name,io_class &io)
{
  if(int fin=open(file_name,O_RDONLY)){
    int AT2MFBDT[]={MFBD_I08T,MFBD_I16T,MFBD_I32T,MFBD_F32T,MFBD_F64T,MFBD_I64T};
    int one=1;
    int t_endian=(*(char*)&one==0);                      // an endian detector, taken from SL's tiff library 
    struct fzhead *fh;
    int sef;
    if((sef=ck_synch_hd(fin,fh,t_endian,io))<0) return -1; // error: type unknown
    close(fin);
    int rv=AT2MFBDT[fh->datyp];
    delete[] fh;
    return rv;
  }
  io.msg(IOL_WARN,"file \"%s\" not found!\n",file_name);
  return -1; // error: type unknown
}

byte *ana_fzread(char *file_name,int *&ds,int &nd,char *&header,int &type,io_class &io) // fzread subroutine	
{
  int AT2MFBDT[]={MFBD_I08T,MFBD_I16T,MFBD_I32T,MFBD_F32T,MFBD_F64T,MFBD_I64T};
  struct stat stat_buf;
  if(stat(file_name,&stat_buf)<0) io.msg(IOL_ERROR,"file \"%s\" not found.\n",file_name);  
  io.msg(IOL_INFO,"reading ANA file \"%s\"\n",file_name);
  int type_sizes[]=ANA_VAR_SZ;
  int one=1;
  int t_endian=(*(char*)&one==0);                      // an endian detector, taken from SL's tiff library 
//
  int fin=open(file_name,O_RDONLY);
  if(!fin){
    io.msg(IOL_WARN,"file \"%s\" not found!\n",file_name);
    return 0;
  }
  struct fzhead	*fh;
  int sef;
  if((sef=ck_synch_hd(fin,fh,t_endian,io))<0) return 0;
//
  (header=new char [1])[0]=0;
  for(int i=0;i<fh->nhb;++i){
    char *tmp=new char [strlen(header)+strlen(fh[i].txt)+1];
    sprintf(tmp,"%s%s",header,fh[i].txt);
    delete[] header;
    header=tmp;
  }
//
  ds=new int [nd=fh->ndim]-1;
  for(int d=1;d<=nd;++d) ds[d]=fh->dim[d-1];

  int n_elem=1;
  for(int d=0;d<=fh->ndim-1;++d) n_elem*=fh->dim[d]; // compute size of array
  type=fh->datyp;
  int f_endian=(fh->subf>=128);                     // the top bit of the byte fh->subf denotes endian type, 0 for little and 1 for big
  int swap_endian=(f_endian!=t_endian);             // file has different endianness
  if(sef) swap_endian=(!swap_endian);               // file contains strange data
  int compressed=(fh->subf&1);                       // file is in compressed format
//
  delete[] fh;
//
  if(compressed){                                  // compressed format
    struct compresshead ch;
    if(read(fin,&ch,14)<14) io.msg(IOL_ERROR,"error reading in compression header\n");
// header by default little-endian?
    if(t_endian){ // big endian platform
      swap(&ch.tsize,1);
      swap(&ch.nblocks,1);
      swap(&ch.bsize,1);
    }
// read data
    int size=ch.tsize-14;
    byte *buf=new byte [size+4]; // +4 bytes to avoid decrunch read beyond data boundary segfault
    if(read(fin,buf,size)<size) io.msg(IOL_ERROR,"error reading in compressed data\n");
//    if(posix_fadvise(fin,0,stat_buf.st_size,POSIX_FADV_DONTNEED)<0) io.msg(IOL_ERROR,"ana_fzread: releasing page cache: %s\n",strerror(errno));
    close(fin);
//
    if(ch.bsize*ch.nblocks>n_elem){ // fix a possible problem with ch.nblocks
      io.msg(IOL_WARN,"warning, bad ch.nblocks = %d\ncorrecting to %d, hope this is right!\n",ch.nblocks,n_elem/ch.bsize);
      ch.nblocks=n_elem/ch.bsize;
    }
    if(ch.type%2==type) io.msg(IOL_ERROR,"inconsistent compression type\n"); // consistency check
    int rv;
    byte *out=new byte [n_elem*type_sizes[type]];
    switch(ch.type){
      case(0): rv=anadecrunch(buf,(int16_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==MFBD_LITTLE_ENDIAN); break;
      case(1): rv=anadecrunch8(buf,(int08_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==MFBD_LITTLE_ENDIAN); break;
      case(2): rv=anadecrunchrun(buf,(int16_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==MFBD_LITTLE_ENDIAN); break;
      case(3): rv=anadecrunchrun8(buf,(int08_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==MFBD_LITTLE_ENDIAN); break;
      case(4): rv=anadecrunch32(buf,(int32_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==MFBD_LITTLE_ENDIAN); break;
      default: io.msg(IOL_ERROR,"error in data type for compressed data, type=%d\n",type);
    }
    delete[] buf;
    type=AT2MFBDT[type];
    return out;
  }else{                            // uncompressed
    size_t size=(size_t)n_elem*(size_t)type_sizes[type];
    byte *out=new byte [size];
//
    size_t rdsz=0;
    while(rdsz<size){
      size_t rv=read(fin,out+rdsz,size-rdsz);
      if(rv<0){
        close(fin);
        fprintf(stderr,"error: unexpected end of file\n");
      }else rdsz+=rv;
    }
//    size_t size_read =read(fin,out,size);
//    printf("Size = %ld \n",size_read);   
//    if(read(fin,out,size)<size){
//      close(fin);
//      io.msg(IOL_ERROR,"error: unexpected end of file\n");
//    }
//    if(posix_fadvise(fin,0,stat_buf.st_size,POSIX_FADV_DONTNEED)<0) io.msg(IOL_ERROR,"ana_fzread: releasing page cache: %s\n",strerror(errno));
    close(fin);
    if(swap_endian) // endianness is wrong
      switch(type){
        case(INT16): swap((int16_t*)out,n_elem); break;
        case(INT32):
        case(FLOAT32): swap((int32_t*)out,n_elem); break;
        case(FLOAT64): swap((int64_t*)out,n_elem); break;
      }
    type=AT2MFBDT[type];
    return out; 
  } // end if(compressed)
}
