#ifndef __FITS_READ_H__      // __FITS_READ_H__
#define __FITS_READ_H__      // __FITS_READ_H__

byte *fits_read(char *file_name,int *&ds,int &nd,char *&header,int &type,io_class &io);
int fits_data_type(const char *file_name,io_class &io);

#endif                      // __FITS_READ_H__
