#include <stdio.h>
#include <string.h>
#include <fitsio.h>

#include "types.h"
#include "io.h"
#include "fits_read.h"
#include "fileio.h"
#include "uts.h"

int fits_data_type(const char *file_name,io_class &io)
{
  int FT2MFBDT[129];
  FT2MFBDT[BYTE_IMG+64]=MFBD_I08T;
  FT2MFBDT[SHORT_IMG+64]=MFBD_I16T;
  FT2MFBDT[LONG_IMG+64]=MFBD_I32T;
  FT2MFBDT[LONGLONG_IMG+64]=MFBD_I64T;
  FT2MFBDT[FLOAT_IMG+64]=MFBD_F32T;
  FT2MFBDT[DOUBLE_IMG+64]=MFBD_F64T;
  int tsz[]=MFBD_TYPE_SIZES;
//
  int status=0;  /* MUST initialize status */
  fitsfile *fptr;
  fits_open_file(&fptr,file_name,READONLY,&status);
  int bitpix,nd;
  long naxes[10];
  fits_get_img_param(fptr,10,&bitpix,&nd,naxes,&status);
  int type=FT2MFBDT[bitpix+64];
  fits_close_file(fptr,&status);
  return type;
}

byte *fits_read(char *file_name,int *&ds,int &nd,char *&header,int &type,io_class &io) // fzread subroutine	
{
  int FT2MFBDT[129];
  FT2MFBDT[BYTE_IMG+64]=MFBD_I08T;
  FT2MFBDT[SHORT_IMG+64]=MFBD_I16T;
  FT2MFBDT[LONG_IMG+64]=MFBD_I32T;
  FT2MFBDT[LONGLONG_IMG+64]=MFBD_I64T;
  FT2MFBDT[FLOAT_IMG+64]=MFBD_F32T;
  FT2MFBDT[DOUBLE_IMG+64]=MFBD_F64T;
  int tsz[]=MFBD_TYPE_SIZES;
//
  io.msg(IOL_INFO,"reading FITS file \"%s\"\n",file_name);
  int status=0;  /* MUST initialize status */
  fitsfile *fptr;
  fits_open_file(&fptr,file_name,READONLY,&status);
  int bitpix;
  long naxes[10];
  fits_get_img_param(fptr,10,&bitpix,&nd,naxes,&status);
  type=FT2MFBDT[bitpix+64];
  
  io.msg(IOL_XNFO,"%d bits per pixel\n",bitpix);
  io.msg(IOL_XNFO,"%d axes\n",nd);
  if(nd!=2) io.msg(IOL_WARN,"file data has %d dimensions, expected 2\n",nd);
  ds=new int [nd]-1;
  ds[1]=naxes[0];
  ds[2]=naxes[1];
  for(int d=1;d<=nd;++d) io.msg(IOL_XNFO,"axis %d: %d\n",d,ds[d]);
  long len=ds[1]*ds[2]*tsz[type],fpixel[10]={1,1};
  int anynul;
  byte *data=new byte [len];
  int ts[]=MFBD_TYPE_SIZES;
  int fits_types[]={TBYTE,TSHORT,TINT,TLONGLONG,TFLOAT,TDOUBLE};
  fits_read_pix(fptr,fits_types[type],fpixel,len/ts[type],0,data,&anynul,&status);
  fits_close_file(fptr,&status);
  fits_report_error(stderr,status);     // print out any error messages
  header=0;
  return data; 
}
