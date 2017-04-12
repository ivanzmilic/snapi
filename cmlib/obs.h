#ifndef __OBS_H__ // __OBS_H__
#define __OBS_H__

#include "types.h"
#include "io.h"
//#include "fileio.h"

class observable{
  fp_t ****S,*lambda; // Stokes vector (N_stokes components x N_wvl)
  int ns,nlambda,nx,ny;
public:
  observable(int ns_in);
  observable(int ns_in, int nx_in, int ny_in);
  observable(int ns_in,int nlambda_in,int nx_in,int ny_in);
  virtual ~observable(void);
  void add(fp_t*,fp_t); // Maybe we will not use this, maybe something like 'set' is better?
  void set(fp_t, int,int,int,int);
  void set(fp_t *,fp_t,int,int,int);
  void set(fp_t **,int,int);
  void write(const char*,io_class&,int,int);
  void write(const char*,io_class&){}; // This is the "full" one
  void read(char*,io_class&);
  fp_t **** get_S();
  fp_t ** get_S(int,int);
  fp_t * get_lambda();
  int get_n_lambda();
  int get_nx();
  int get_ny();
};

#endif          // __OBS_H__
