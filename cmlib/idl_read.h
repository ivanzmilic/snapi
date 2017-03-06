#ifndef __IDL_READ_H__      // __IDL_READ_H__
#define __IDL_READ_H__      // __IDL_READ_H__

#define IDL_INT8	        1
#define IDL_INT16	        2
#define IDL_INT32	        4
#define IDL_FLOAT32	        5
#define IDL_FLOAT64	        6

#define IDL_VAR_SZ {0,1,2,4,4,8}

byte *idl_read(char *file_name,int *&ds,int &nd,char *&header,int &type,io_class &io);

#endif                      // __ANA_READ_H__
