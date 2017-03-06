#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "types.h"
#include "uts.h"
#include "io.h"
#include "anacompress.h"
#include "ana_io.h"

void ana_fcwrite(byte *data,char *file_name,int *ds,int nd,char *header,int type,io_class &io)	/* fcwrite subroutine */
{ // write standard f0 files, compressed format
  FILE *f=fopen(file_name,"w");
//
  int nhead=(header)?min((strlen(header)+254)/255,15):1;
  struct fzhead *fh=new fzhead [nhead];
  memset(fh,0,nhead*sizeof(struct fzhead));
//
  int one=1;
  int t_endian=(*(char*)&one==0);    // an endian detector, taken from SL's tiff library 
  if(t_endian){ // BIG_ENDIAN
    fh->synch_pattern=0xaaaa5555;
    fh->subf=129;
  }else{        // LITTLE_ENDIAN
    fh->synch_pattern=0x5555aaaa;
    fh->subf=1;
  }
  fh->source=0;
  fh->nhb=nhead;					// may be changed later
  fh->datyp=type;
  fh->ndim=nd;
//
  int n_elem=1;
  for(int n=1;n<=nd;++n) n_elem*=(fh->dim[n-1]=ds[n]);
  int nx=fh->dim[0];
  int ny=n_elem/nx;
  int type_sizes[]=ANA_VAR_SZ;
  int size=n_elem*type_sizes[type];
  if(t_endian){ // big endian platform
    switch(type){
      case(INT16): swap((int16_t*)data,n_elem); break;
      case(INT32):
      case(FLOAT32): swap((int32_t*)data,n_elem); break;
      case(FLOAT64): swap((int64_t*)data,n_elem); break;
    }
    swap(fh->dim,nd);
  }
  if(header){
    char *p=header;
    for(int i=0;i<nhead;++i){
      int len=min(strlen(p),255);
      strncpy(fh[i].txt,p,len);
      fh[i].txt[len]=0;
      p+=len;
    }
  }
/* now compress the array, must be a byte or short */
/* extended to 32 bits 2/4/96 */
  int res,crunch_slice=5,runlengthflag=0,limit=size+size/2; // reserve a bit extra just in case
  byte *q=new byte [limit];
  switch(type){
    case(0):{
      if(runlengthflag)
        res=anacrunchrun8(q,data,crunch_slice,nx,ny,limit,t_endian);
      else
        res=anacrunch8(q,data,crunch_slice,nx,ny,limit,t_endian);
      break;
    }
    case(1):{
      if(runlengthflag)
        res=anacrunchrun(q,(int16_t*)data,crunch_slice,nx,ny,limit,t_endian);
      else
        res=anacrunch(q,(int16_t*)data,crunch_slice,nx,ny,limit,t_endian);
      break;
    }
    case(2):{
      if(runlengthflag){
        io.msg(IOL_WARN,"FCRUNWRITE not supported for I*4 yet\n");
        fclose(f);
        delete q;
        delete[] fh;
        return;
      }else
        res=anacrunch32(q,(int32_t*)data,crunch_slice,nx,ny,limit,t_endian);
      break;
    }
    default:{
      io.msg(IOL_WARN,"FCWRITE: unsupported type\n");
      fclose(f);
      delete[] q;
      delete[] fh;
      return;
    }
  }
  if(res<0){
    io.msg(IOL_WARN,"not enough space allocated (%d bytes) for compressed array, trying uncompressed!\n",limit);
    delete[] q;
    delete[] fh;
    fclose(f);
    ana_fzwrite(data,file_name,ds,nd,header,type,io);
    return;
  }
  if(res>size){
    io.msg(IOL_WARN,"compressed data (%d bytes) larger than raw data (%d bytes), writing uncompressed!\n",limit,size);
    delete[] q;
    delete[] fh;
    fclose(f);
    ana_fzwrite(data,file_name,ds,nd,header,type,io);
    return;
  }
  size=res;
  if(t_endian) swap(&res,1);
  for(int i=0;i<=3;++i) fh->cbytes[i]=((byte*)&res)[i];
  fwrite(fh,nhead*sizeof(struct fzhead),1,f);  // write header
  fwrite(q,1,size,f);
  delete[] q;
  delete[] fh;
  fclose(f);
}

void ana_fzwrite(byte *data,char *file_name,int *ds,int nd,char *header,int type,io_class &io)	/* fcwrite subroutine */
{ // write standard f0 files, compressed format
  FILE *f=fopen(file_name,"w");
//
  int nhead=(header)?min((strlen(header)+254)/255,15):1;
  struct fzhead *fh=new fzhead [nhead];
  memset(fh,0,nhead*sizeof(struct fzhead));
//
  int one=1;
  int t_endian=(*(char*)&one==0);    // an endian detector, taken from SL's tiff library 
  if(t_endian){ // BIG_ENDIAN
    fh->synch_pattern=0xaaaa5555;
  }else{        // LITTLE_ENDIAN
    fh->synch_pattern=0x5555aaaa;
  }
  fh->source=0;
  fh->nhb=nhead;					// may be changed later
  fh->datyp=type;
  fh->ndim=nd;
//
  int n_elem=1;
  for(int n=1;n<=nd;++n) n_elem*=(fh->dim[n-1]=ds[n]);
  int nx=fh->dim[0];
  int ny=n_elem/nx;
  int type_sizes[]=ANA_VAR_SZ;
  int size=n_elem*type_sizes[type];
  if(t_endian){ // big endian platform
    switch(type){
      case(INT16): swap((int16_t*)data,n_elem); break;
      case(INT32):
      case(FLOAT32): swap((int32_t*)data,n_elem); break;
      case(FLOAT64): swap((int64_t*)data,n_elem); break;
    }
    swap(fh->dim,nd);
  }
  if(header){
    char *p=header;
    for(int i=0;i<nhead;++i){
      int len=min(strlen(p),255);
      strncpy(fh[i].txt,p,len);
      fh[i].txt[len]=0;
      p+=len;
    }
  }
  int res=size;
  if(t_endian) swap(&res,1);
  for(int i=0;i<=3;++i) fh->cbytes[i]=((byte*)&res)[i];
  fwrite(fh,nhead*sizeof(struct fzhead),1,f);  // write header
  fwrite(data,1,size,f);
  delete[] fh;
  fclose(f);
}

ana_file::ana_file(const char *file_name,int *ds,int nd,const char *header,int type_in,io_class &io_in):io(io_in)
{ // write standard f0 files, compressed format
  if(f=fopen(file_name,"w")){
    struct fzhead fh;
    memset(&fh,0,sizeof(struct fzhead));
    int one=1;
    t_endian=(*(char*)&one==0);    // an endian detector, taken from SL's tiff library 
    if(t_endian){ // BIG_ENDIAN
      fh.synch_pattern=0xaaaa5555;
    }else{        // LITTLE_ENDIAN
      fh.synch_pattern=0x5555aaaa;
    }
    fh.source=0;
    fh.nhb=1;					// may be changed later
    fh.datyp=type_in;
    fh.ndim=nd;
//
    elem=1;
    for(int n=1;n<=nd;++n) elem*=(fh.dim[n-1]=ds[n]);
    int type_sizes[]=ANA_VAR_SZ;
//
    type=type_in;
    tsize=type_sizes[type];
    size=elem*tsize;
    felem=0;
    fsize=0;
//
    if(t_endian) swap(fh.dim,nd);
    if(header){
      int len=min(strlen(header),255);
      strncpy(fh.txt,header,len);
      fh.txt[len]=0;
    }
    int res=size;
    if(t_endian) swap(&res,1);
    for(int i=0;i<=3;++i) fh.cbytes[i]=((byte*)&res)[i];
    fwrite(&fh,sizeof(struct fzhead),1,f);  // write header
  }
}

ana_file::~ana_file(void)
{
  if((fsize!=size)||(felem!=elem)) io.msg(IOL_ERROR,"ana_file::~ana_file: wrote %d elements (%d bytes) but wrote %d (%d bytes)\n",elem,size,felem,fsize);
  if(f) fclose(f);
}

int ana_file::append(void *data,int n,int type_in)
{
  if(f)
    if(type_in==type){
      int sz=n*tsize;
      if(t_endian){ // big endian platform
        switch(type){
          case(INT16): swap((int16_t*)data,n); break;
          case(INT32):
          case(FLOAT32): swap((int32_t*)data,n); break;
          case(FLOAT64): swap((int64_t*)data,n); break;
        }
      }
      if(fwrite(data,sz,1,f)==1){
        fsize+=sz;
        felem+=n;
        return sz;
      }
    }
  return -1;
}
