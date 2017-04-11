#include "types.h"
#include "mem.h"
#include "io.h"
#include "ana_io.h"
#include "idl_read.h"
#include "fits_read.h"
#include "fileio.h"

fp_t *float32to64conv(byte *data,int n_elem)
{
  fp_t *fdata=new fp_t [n_elem];
  for(int i=0;i<=n_elem-1;++i) fdata[i]=(fp_t)((float32_t*)data)[i];
  delete[] data;
  return fdata;
}

fp_t *float64to32conv(byte *data,int n_elem)
{
  fp_t *fdata=new fp_t [n_elem];
  for(int i=0;i<=n_elem-1;++i) fdata[i]=(fp_t)((float64_t*)data)[i]; // let's hope the data fit in 32 bit floats...
  delete[] data;
  return fdata;
}

fp_t ***read_file(char *file_name,int &nx1,int &nx2,int &nx3,io_class &io)
{
  char *header;
  int type,nd,*dn;
  byte *data=ana_fzread(file_name,dn,nd,header,type,io);
  if(!data) return 0;
  if(nd!=3) io.msg(IOL_ERROR,"Number of dimensions of data in \"%s\" is not 3!\n",file_name);
  delete[] header;
  nx1=dn[1];
  nx2=dn[2];  
  nx3=dn[3];  
  delete[] (dn+1);
  fp_t *fdata;
  int t_sz[]=MFBD_TYPE_SIZES;
  switch(type){
    case(MFBD_F32T): fdata=(t_sz[type]!=sizeof(fp_t))?float32to64conv(data,nx1*nx2*nx3):(fp_t*)data; break;
    case(MFBD_F64T): fdata=(t_sz[type]!=sizeof(fp_t))?float64to32conv(data,nx1*nx2*nx3):(fp_t*)data; break;
    default: io.msg(IOL_ERROR,"data type of \"%s\" is not floating point\n",file_name);
  }
  fp_t ***p=ft3dim(fdata,1,nx3,1,nx2,1,nx1);
  fp_t ***q=ft3dim(1,nx1,1,nx2,1,nx3);
  for(int x1=1;x1<=nx1;++x1)
    for(int x2=1;x2<=nx2;++x2)
      for(int x3=1;x3<=nx3;++x3) q[x1][x2][x3]=p[x3][x2][x1];
  del_ft3dim(p,1,nx3,1,nx2,1,nx1);
  return q;
}
fp_t ****read_file(char *file_name,int &nx1,int &nx2,int &nx3,int &nx4,io_class &io)
{
  char *header;
  int type,nd,*dn;
  byte *data=ana_fzread(file_name,dn,nd,header,type,io);

  if(!data) return 0;
  if(nd!=4){ 
    io.msg(IOL_ERROR,"Number of dimensions of data in \"%s\" is not 4!\n",file_name);
    return 0;
  }
  delete[] header;
  nx1=dn[1];
  nx2=dn[2];  
  nx3=dn[3];
  nx4=dn[4];  
  delete[] (dn+1);
  fp_t *fdata;
  int t_sz[]=MFBD_TYPE_SIZES;
  switch(type){
    case(MFBD_F32T): fdata=(t_sz[type]!=sizeof(fp_t))?float32to64conv(data,nx1*nx2*nx3*nx4):(fp_t*)data; break;
    case(MFBD_F64T): fdata=(t_sz[type]!=sizeof(fp_t))?float64to32conv(data,nx1*nx2*nx3*nx4):(fp_t*)data; break;
    default: io.msg(IOL_ERROR,"data type of \"%s\" is not floating point\n",file_name);
  }
  fp_t ****p=ft4dim(fdata,1,nx4,1,nx3,1,nx2,1,nx1);
  fp_t ****q=ft4dim(1,nx1,1,nx2,1,nx3,1,nx4);
  for(int x1=1;x1<=nx1;++x1)
    for(int x2=1;x2<=nx2;++x2)
      for(int x3=1;x3<=nx3;++x3)
        for (int x4=1;x4<=nx4;++x4) q[x1][x2][x3][x4]=p[x4][x3][x2][x1];
  del_ft4dim(p,1,nx4,1,nx3,1,nx2,1,nx1);
  return q;
}

fp_t **read_file(char *file_name,int &nx1,int &nx2,io_class &io)
{
  char *header;
  int type,nd,*dn;
  byte *data=ana_fzread(file_name,dn,nd,header,type,io);
  if(!data) return 0;
  if(nd!=2) io.msg(IOL_ERROR,"Number of dimensions of data in \"%s\" is not 2!\n",file_name);
  delete[] header;
  nx1=dn[2];
  nx2=dn[1];  
  delete[] (dn+1);
  fp_t *fdata;
  int t_sz[]=MFBD_TYPE_SIZES;
  switch(type){
    case(MFBD_F32T): fdata=(t_sz[type]!=sizeof(fp_t))?float32to64conv(data,nx1*nx2):(fp_t*)data; break;
    case(MFBD_F64T): fdata=(t_sz[type]!=sizeof(fp_t))?float64to32conv(data,nx1*nx2):(fp_t*)data; break;
    default: io.msg(IOL_ERROR,"data type of \"%s\" is not floating point\n",file_name);
  }
  return ft2dim(fdata,1,nx1,1,nx2);
/*
  fp_t **p=ft2dim(fdata,1,nx2,1,nx1);
  fp_t **q=ft2dim(1,nx1,1,nx2);
  for(int x1=1;x1<=nx1;++x1)
    for(int x2=1;x2<=nx2;++x2) q[x1][x2]=p[x2][x1];
  del_ft2dim(p,1,nx2,1,nx1);
  return q;
*/
}

fp_t *read_file(char *file_name,int &nx1,io_class &io)
{
  char *header;
  int type,nd,*dn;
  byte *data=ana_fzread(file_name,dn,nd,header,type,io);
  if(!data) return 0;
  if(nd!=1) io.msg(IOL_ERROR,"Number of dimensions of data in \"%s\" is not 1!\n",file_name);
  delete[] header;
  nx1=dn[1];
  delete[] (dn+1);
  fp_t *fdata;
  int t_sz[]=MFBD_TYPE_SIZES;
  switch(type){
    case(MFBD_F32T): fdata=(t_sz[type]!=sizeof(fp_t))?float32to64conv(data,nx1):(fp_t*)data; break;
    case(MFBD_F64T): fdata=(t_sz[type]!=sizeof(fp_t))?float64to32conv(data,nx1):(fp_t*)data; break;
    default: io.msg(IOL_ERROR,"data type of \"%s\" is not floating point\n",file_name);
  }
  return fdata-1;
}

#define FT_ERROR  -1
#define FT_ANA  1
#define FT_IDL  2
#define FT_FITS 3

int file_type(const char *name)
{
  byte h[4];
  FILE *f=fopen(name,"r");
  if(f){
    if(fread(h,1,4,f)<4){
      fclose(f);
      return FT_ERROR;
    }
    fclose(f);
    if((h[0]==170)&&(h[1]==170)&&(h[2]==85)&&(h[3]==85)) return FT_ANA;
    if((h[0]==83)&&(h[1]==73)&&(h[2]==77)&&(h[3]==80)) return FT_FITS;
    return FT_IDL;
  }
  return FT_ERROR;
}

byte *read_image(char *file_name,int &nx1,int &nx2,int &type,char *&header,io_class &io)
{
  int nd,*dn;
  byte *data;  
  switch(file_type(file_name)){
    case(FT_ANA): data=ana_fzread(file_name,dn,nd,header,type,io); break;
    case(FT_IDL): data=idl_read(file_name,dn,nd,header,type,io); break;
    case(FT_FITS): data=fits_read(file_name,dn,nd,header,type,io); break;
    case(FT_ERROR):{
      io.msg(IOL_ERROR,"file \"%s\" not found.\n",file_name);
      return 0;
    }
    default:{
      io.msg(IOL_ERROR,"file \"%s\" has unknown format.\n",file_name);
      return 0;
    }
  }
  if(!data) return 0;
  if(nd!=2) io.msg(IOL_ERROR,"Number of dimensions of data in \"%s\" is not 2!\n",file_name);
  nx2=dn[1];
  nx1=dn[2];
  delete[] (dn+1);
  return data;
}

byte *read_file(char *file_name,int &nx1,int &nx2,int &nx3,int &type,char *&header,io_class &io)
{
  int nd,*dn;
  byte *data;  
  switch(file_type(file_name)){
    case(FT_ANA): data=ana_fzread(file_name,dn,nd,header,type,io); break;
    case(FT_IDL): data=idl_read(file_name,dn,nd,header,type,io); break;
    case(FT_FITS): data=fits_read(file_name,dn,nd,header,type,io); break;
    case(FT_ERROR):{
      io.msg(IOL_ERROR,"file \"%s\" not found.\n",file_name);
      return 0;
    }
    default:{
      io.msg(IOL_ERROR,"file \"%s\" has unknown format.\n",file_name);
      return 0;
    }
  }
  if(!data) return 0;
  switch(nd){
    case(3):{
      nx3=dn[1];
      nx2=dn[2];
      nx1=dn[3];
      break;
    }
    case(2):{
      nx3=dn[1];
      nx2=dn[2];
      nx1=1;
      break;
    }
    default: io.msg(IOL_ERROR,"Incompatible number of dimensions of data in \"%s\" (%d)!\n",file_name,nd);
  }
  delete[] (dn+1);
  return data;
}

int data_type(const char *file_name,io_class &io)
{
  int nd,*dn;
  byte *data;  
  switch(file_type(file_name)){
    case(FT_ANA): return ana_data_type(file_name,io);
    case(FT_IDL): return MFBD_I16T; // the only supported type...
    case(FT_FITS): return fits_data_type(file_name,io);
    case(FT_ERROR):{
      io.msg(IOL_ERROR,"file \"%s\" not found.\n",file_name);
      return -1;
    }
    default:{
      io.msg(IOL_ERROR,"file \"%s\" has unknown format.\n",file_name);
      return -1;
    }
  }
  return -1;
}

int16_t **read_file(char *file_name,int &nx1,int &nx2,char *&header,io_class &io)
{
  int type,nd,*dn;
  int16_t *data;  
  switch(file_type(file_name)){
    case(FT_ANA): data=(int16_t*)ana_fzread(file_name,dn,nd,header,type,io); break;
    case(FT_IDL): data=(int16_t*)idl_read(file_name,dn,nd,header,type,io); break;
    case(FT_FITS): data=(int16_t*)fits_read(file_name,dn,nd,header,type,io); break;
    case(FT_ERROR):{
      io.msg(IOL_ERROR,"file \"%s\" not found.\n",file_name);
      return 0;
    }
    default:{
      io.msg(IOL_ERROR,"file \"%s\" has unknown format.\n",file_name);
      return 0;
    }
  }
  if(!data) return 0;
  if(nd!=2) io.msg(IOL_ERROR,"Number of dimensions of data in \"%s\" is not 2!\n",file_name);
  nx1=dn[1];
  nx2=dn[2];
  delete[] (dn+1);
  int16_t **p=i16t2dim(data,1,nx2,1,nx1);
  int16_t **q=i16t2dim(1,nx1,1,nx2);
  for(int x1=1;x1<=nx1;++x1)
    for(int x2=1;x2<=nx2;++x2) q[x1][x2]=p[x2][x1];
  del_i16t2dim(p,1,nx2,1,nx1);
  return q;
}

#undef FT_ANA
#undef FT_IDL

void write_file(char *name,fp_t ***data,int nx1,int nx2,int nx3,io_class &io)
{
  int *ds=new int [3]-1;
  ds[1]=nx1;
  ds[2]=nx2;
  ds[3]=nx3;
  float ***q=f3dim(1,nx3,1,nx2,1,nx1);
  for(int x1=1;x1<=nx1;++x1)
    for(int x2=1;x2<=nx2;++x2)
      for(int x3=1;x3<=nx3;++x3) q[x3][x2][x1]=data[x1][x2][x3];
  ana_fzwrite((byte*)(q[1][1]+1),name,ds,3,0,FLOAT32,io);
  del_f3dim(q,1,nx3,1,nx2,1,nx1);
  delete (ds+1);
}

void write_file(char *name,fp_t ****data,int nx1,int nx2,int *nx3,int nx4,io_class &io)
{
  int nxt=0;
  for(int x2=1;x2<=nx2;++x2) nxt+=nx3[x2];
  nxt*=nx4;
//
  int *ds=new int [2]-1;
  ds[1]=nxt;
  ds[2]=nx1;
  float **data_out=f2dim(1,nx1,1,nxt);
  for(int k=1;k<=nx1;++k)
    for(int x2=1,x=0;x2<=nx2;++x2)
      for(int x3=1;x3<=nx3[x2];++x3)
        for(int x4=1;x4<=nx4;++x4)
          data_out[k][++x]=data[k][x2][x3][x4];
  ana_fzwrite((byte*)(data_out[1]+1),name,ds,2,0,FLOAT32,io);
  del_f2dim(data_out,1,nx1,1,nxt);
  delete[] (ds+1);
}

void write_file(char *name,float32_t ****data,int nx1,int nx2,int nx3,int nx4,io_class &io)
{
  int *ds=new int [4]-1;
  ds[1]=nx4;
  ds[2]=nx3;
  ds[3]=nx2;
  ds[4]=nx1;
  ana_fzwrite((byte*)(data[1][1][1]+1),name,ds,4,0,FLOAT32,io);
  delete[] (ds+1);
}
