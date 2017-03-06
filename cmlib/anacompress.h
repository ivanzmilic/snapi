#ifndef __ANACOMPRESS_H__      // __ANACOMPRESS_H__
#define __ANACOMPRESS_H__      // __ANACOMPRESS_H__

int anacrunchrun8(byte *x,byte *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunch8(byte *x,byte *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunchrun(byte *x,int16_t *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunch(byte *x,int16_t *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunch32(byte *x,int32_t *array,int slice,int nx,int ny,int limit,int t_endian);

#endif                         // __ANACOMPRESS_H__
