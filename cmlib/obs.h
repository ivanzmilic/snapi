#ifndef __OBS_H__ // __OBS_H__
#define __OBS_H__

#include "types.h"
#include "io.h"
#include "fileio.h"
#include "ana_io.h"

class observable{
  fp_t ****S,*lambda,*mask; // Stokes vector (N_stokes components x N_wvl)
  int ns,nlambda,nx,ny;
  fp_t scattered_light,spectral_broadening,obs_qs,synth_qs;
  
  int n_spsf; // Number of psf point is psf is externally prescribed
  fp_t * spsf; // actual psf file
  
  fp_t el,az;
  int to_invert;
  int no_iterations; // number of L-M iterations to use
  fp_t start_lambda; // starting value of lambda for LM
  fp_t stopping_chisq; // value of chisq (reduced) to stop LM iterations at
  fp_t * w_stokes;
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
  
  // Sets:
  void set(fp_t, int,int,int,int);
  void set(fp_t *,fp_t,int,int,int);
  void set(fp_t **,int,int);
  void set(fp_t ****);
  void set_lambda(fp_t *);
  void set_mask(fp_t*);
  void set_inv_parameters(fp_t,fp_t,fp_t,fp_t);
  void set_viewing_angle(fp_t, fp_t);
  void set_to_invert(int);
  void set_no_iterations(int);
  void set_start_lambda(fp_t);
  void set_stopping_chisq(fp_t);
  void set_w_stokes(fp_t *);
  void set_n_spsf(int);
  void set_spsf(fp_t *);
  
  int get_to_invert();
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
  fp_t get_scattered_light();
  fp_t get_spectral_broadening();
  fp_t get_synth_qs();
  fp_t get_el();
  fp_t get_az();
  int get_no_iterations();
  fp_t get_start_lambda();
  fp_t get_stopping_chisq();
  fp_t * get_w_stokes();
  int get_n_spsf();
  fp_t * get_spsf();

  void normalize();
  void add_scattered_light(fp_t, fp_t);
  void spectral_convolve(fp_t,int,int);
  void psf_convolve(int,fp_t*,int,int);
};

observable * obs_new(int nx,int ny,int ns,int nlambda);
observable * obs_new(uint08_t*,int32_t&,uint08_t,io_class&);

#endif          // __OBS_H__
