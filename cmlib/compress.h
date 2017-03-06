#ifndef __COMPRESS_H__    // __COMPRESS_H__
#define __COMPRESS_H__

#include "types.h"
#include "io.h"

byte *z_compress(byte*,int&,int,byte,io_class&);
byte *z_uncompress(byte*,int&,byte,io_class&);

#endif                    // __COMPRESS_H__
