#include <stdio.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "idl_read.h"
#include "fileio.h"
#include "uts.h"

byte *idl_read(char *file_name,int *&ds,int &nd,char *&header,int &type,io_class &io) // fzread subroutine	
{
  int IDLT2MFBDT[]={0,MFBD_I08T,MFBD_I16T,MFBD_I32T,MFBD_F32T,MFBD_F64T};
  io.msg(IOL_INFO,"reading IDL file \"%s\"\n",file_name);
  int type_sizes[]=IDL_VAR_SZ;
  int one=1;
  int t_endian=(*(char*)&one==0); // an endian detector, taken from SL's tiff library: 0=LE, 1=BE
//
  FILE *fin=fopen(file_name,"r");
  if(!fin){
    io.msg(IOL_WARN,"file \"%s\" not found!\n",file_name);
    return 0;
  }
  header=new char [513];
  if(fread(header,1,512,fin)!=512) io.msg(IOL_ERROR,"error in idl_read while reading header\n");
  sscanf(strstr(header,"dims=")+5,"%d",&nd);
  ds=new int [nd]-1;
  sscanf(strstr(header,"type=")+5,"%d",&type);
  sscanf(strstr(header,"nx=")+3,"%d",&(ds[1]));
  sscanf(strstr(header,"ny=")+3,"%d",&(ds[2]));
  char end;
  char *p=strstr(header,"endian=");
  if(p) sscanf(p+3,"%c",&end); else end='l';
  int swap_endian=((end=='l')==t_endian);
//
  int n_elem=ds[1]*ds[2];
  int size=n_elem*type_sizes[type];
  byte *out=new byte [size];
  if(fread(out,1,size,fin)<size){
    fclose(fin);
    io.msg(IOL_ERROR,"error: unexpected end of file\n");
  }
  fclose(fin);
  if(swap_endian) // endianness is wrong
    switch(type){
      case(IDL_INT16): swap((int16_t*)out,n_elem); break;
      case(IDL_INT32):
      case(IDL_FLOAT32): swap((int32_t*)out,n_elem); break;
      case(IDL_FLOAT64): swap((int64_t*)out,n_elem); break;
    }
  type=IDLT2MFBDT[type];
  return out; 
}
