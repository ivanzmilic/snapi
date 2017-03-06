#ifndef __OBS_H__ // __OBS_H__
#define __OBS_H__

#include "types.h"
#include "io.h"

class observable{
  fp_t **S,*lambda; // Stokes vector (N_stokes components x N_wvl)
  int ns,nlambda;
public:
  observable(int ns_in);
  virtual ~observable(void);
  void add(fp_t*,fp_t);
  void readsingle(const char*); // This reads the spectrum from the 
  void write(const char*,io_class&);
  fp_t ** get_S();
  fp_t * get_lambda();
};

#endif          // __OBS_H__
