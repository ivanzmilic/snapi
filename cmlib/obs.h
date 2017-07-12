#ifndef __OBS_H__ // __OBS_H__
#define __OBS_H__

#include "types.h"
#include "io.h"
#include "fileio.h"
#include "ana_io.h"

class observable{
  fp_t ****S,*lambda,*mask; // Stokes vector (N_stokes components x N_wvl)
  int ns,nlambda,nx,ny;
public:
  observable(int ns_in);
  observable(int ns_in, int nx_in, int ny_in);
  observable(int ns_in,int nlambda_in,int nx_in,int ny_in);
  virtual ~observable(void);

  observable(uint08_t*,int32_t&,uint08_t,io_class&); // create from buffer

  virtual int32_t size(io_class&);
  virtual int32_t pack(uint08_t*,uint08_t,io_class&);
  virtual int32_t unpack(uint08_t*,uint08_t,io_class&);

  void add(fp_t*,fp_t); // Maybe we will not use this, maybe something like 'set' is better?
  void set(fp_t, int,int,int,int);
  void set(fp_t *,fp_t,int,int,int);
  void set(fp_t **,int,int);
  void set(fp_t ****);
  void setlambda(fp_t *);
  void setmask(fp_t*);
  void write(const char*,io_class&,int,int);
  void write(const char*,io_class&){}; // This is the "full" one
  void read(char*,io_class&);
  observable* extract(int,int,int,int,int,int);
  fp_t **** get_S();
  fp_t ** get_S(int,int);
  fp_t ** get_S_to_fit(int,int);
  fp_t * get_lambda();
  fp_t * get_lambda_to_fit();
  fp_t * get_mask();
  int get_n_lambda();
  int get_n_lambda_to_fit();
  int get_nx();
  int get_ny();

  void normalize();
  void correct_for_scattered_light(fp_t);
  void spectral_convolve(fp_t,int,int);
};

observable * obs_new(int nx,int ny,int ns,int nlambda);
observable * obs_new(uint08_t*,int32_t&,uint08_t,io_class&);

#endif          // __OBS_H__
