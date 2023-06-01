#include <stdio.h>
#include <string.h>

#include "types.h"
#include "pack.h"
#include "mem.h"
#include "compress.h"
#include "io.h"
#include "jnfo.h"

jnfo::jnfo(byte *buf,byte swap_endian,io_class &io)
{
  memset(this,0,sizeof(jnfo));
// unpack the buffer
  int size;
  byte *data=z_uncompress(buf,size,swap_endian,io);
//
  int offs=0;
 //
  offs+=unpack(data+offs,na,swap_endian);
  if(na){
    atmos=new atmosphere* [na];
    for(int a=0;a<na;++a) atmos[a]=atmos_new(data,offs,swap_endian,io);
  }else atmos=0;
  
  offs+=unpack(data+offs,nm,swap_endian);
  if(nm){
    models=new model* [nm];
    for(int m=0;m<nm;++m) models[m]=model_new(data,offs,swap_endian,io);

  }else models=0;
  
//
  offs+=unpack(data+offs,no,swap_endian);  
  if(no){
    offs+=unpack(data+offs,az=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,el=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,nlambda=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,to_invert=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,return_model=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,return_atmos=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,extra_settings=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,xl=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,xh=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,yl=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,yh=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,ll=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,lh=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,scattered_light=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,spectral_broadening=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,obs_qs=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,synth_qs=new fp_t [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,no_iterations=new int [no],0,no-1,swap_endian);
    offs+=unpack(data+offs,starting_lambda=new fp_t [no],0,no-1,swap_endian); 
    offs+=unpack(data+offs,stopping_chisq=new fp_t [no],0,no-1,swap_endian); 
    lambda=new fp_t* [no];
    weights = new fp_t*[no];
    w_stokes = new fp_t*[no];
    name=new char* [no];
    for(int o=0;o<no;++o){
      offs+=unpack(data+offs,lambda[o]=new fp_t [nlambda[o]],0,nlambda[o]-1,swap_endian);
      offs+=unpack(data+offs,weights[o]=new fp_t [nlambda[o]],0,nlambda[o]-1,swap_endian);
      offs+=unpack(data+offs,w_stokes[o]=new fp_t [4],0,3,swap_endian);
      offs+=unpack(data+offs,name[o]);
    }
  }
  
  if (nm){
    offs+=unpack(data+offs,read_model_from_file=new int [nm],0,nm-1,swap_endian);
    input_models = new char*[nm];
    for (int m=0;m<nm;++m)
      if (read_model_from_file[m])
        offs+=unpack(data+offs,input_models[m]);
      else 
        input_models[m]=0;
  }
//
  offs+=unpack(data+offs,uid,swap_endian);
  offs+=unpack(data+offs,gid,swap_endian);
  offs+=unpack(data+offs,uname);
//
  offs+=unpack(data+offs,cdcl);
  offs+=unpack(data+offs,nmth);
  offs+=unpack(data+offs,nsth);
  offs+=unpack(data+offs,nsln);
//
  if(offs!=size) io.msg(IOL_WARN,"jnfo::jnfo: number of unpacked bytes was %d but %d bytes received (possible version conflict)?!\n",offs,size);
  delete[] data;
}

jnfo::jnfo(void)
{
  memset(this,0,sizeof(jnfo));
}

jnfo::~jnfo(void)
{
  if(na){
    for(int a=0;a<na;++a) delete atmos[a];
    delete[] atmos;
  }
  if (nm){
    for (int m=0;m<nm;++m) delete models[m];
    delete[]models;
    for (int m=0;m<nm;++m)
      if (input_models[m])
        delete[]input_models[m];
    delete[]input_models;
    delete[]read_model_from_file;

  }
  if(no){
    if(az) delete[] az;
    if(el) delete[] el;
    if(nlambda) delete[] nlambda;
    if(to_invert) delete[] to_invert;
    if(return_model) delete[] return_model;
    if(return_atmos) delete[] return_atmos;
    if(extra_settings) delete[] extra_settings;
    if (xl) delete[] xl;
    if (xh) delete[] xh;
    if (yl) delete[] yl;
    if (yh) delete[] yh;
    if (ll) delete[] ll;
    if (lh) delete[] lh;
    if(lambda) for(int o=0;o<no;++o) if(lambda[o]) delete[] lambda[o];
    if(lambda) delete[] lambda;
    if(weights) for(int o=0;o<no;++o) if(weights[o]) delete[] weights[o];
    if(weights) delete[] weights;
    if (w_stokes) for (int o=0;o<no;++o) delete []w_stokes[o];
    if (w_stokes) delete []w_stokes;
    if(name) for(int o=0;o<no;++o) delete[] name[o];
    if(name) delete[] name;
    if (scattered_light) delete[]scattered_light;
    if (spectral_broadening) delete []spectral_broadening; 
    if (synth_qs) delete []synth_qs;
    if (obs_qs) delete []obs_qs;
    if (no_iterations) delete []no_iterations;
    if (starting_lambda) delete []starting_lambda;
    if (stopping_chisq) delete []stopping_chisq;
  }
  if(uname) delete[] uname;
  
}

byte *jnfo::compress(int &size,int level,byte swap_endian,io_class &io)
{
  int32_t usize=sizeof(int);                     // na
  for(int a=0;a<na;++a) usize+=atmos[a]->size(io);
  usize+=sizeof(int);                     // nm
  for(int m=0;m<nm;++m) usize+=models[m]->size(io);
    
//  io.msg(IOL_INFO,"jnfo::compress: usize %d\n",usize);
//
  usize+=sizeof(int);                            // no
  if(no){
    usize+=8*no*sizeof(fp_t);   // az,el,scat_l,broadedning,qs obs and synth,starting_lambda, stopping chisq
    usize+=2*no*sizeof(int);    // nlambda,no_iterations
    usize+=no*10*sizeof(int);   // to_invert,return_model,return_atmos
                                // extra_settings, xl,xh,yl,yh,ll,lh
    for(int o=0;o<no;++o){
      usize+=2*nlambda[o]*sizeof(fp_t);  // lambda,weights
      usize+=4*sizeof(fp_t); // w_stokes
      usize+=strlen(name[o])+1;        // file name
    }
  }
  if (nm){
    usize +=nm*sizeof(int); // read_model_from_file
    for (int m=0;m<nm;++m)
      if (input_models[m])
        usize+=strlen(input_models[m])+1;
  }
  usize+=sizeof(uid_t)+sizeof(gid_t);            // uid,gid
  usize+=strlen(uname)+1;                        // uname
  usize+=4*sizeof(uint08_t);                     // cdcl,nmth,nsth,nsln
//
  byte *data=new byte [usize];
//
  int offs=pack(data,na,swap_endian);
  for(int a=0;a<na;++a) offs+=atmos[a]->pack(data+offs,swap_endian,io);
  
  offs+=pack(data+offs,nm,swap_endian);
  for(int m=0;m<nm;++m) offs+=models[m]->pack(data+offs,swap_endian,io);
//  io.msg(IOL_INFO,"jnfo::compress: %d %d\n",offs,usize);
//
  offs+=pack(data+offs,no,swap_endian);
  if(no){
    offs+=pack(data+offs,az,0,no-1,swap_endian);
    offs+=pack(data+offs,el,0,no-1,swap_endian);
    offs+=pack(data+offs,nlambda,0,no-1,swap_endian);
    offs+=pack(data+offs,to_invert,0,no-1,swap_endian);
    offs+=pack(data+offs,return_model,0,no-1,swap_endian);
    offs+=pack(data+offs,return_atmos,0,no-1,swap_endian);
    offs+=pack(data+offs,extra_settings,0,no-1,swap_endian);
    offs+=pack(data+offs,xl,0,no-1,swap_endian);
    offs+=pack(data+offs,xh,0,no-1,swap_endian);
    offs+=pack(data+offs,yl,0,no-1,swap_endian);
    offs+=pack(data+offs,yh,0,no-1,swap_endian);
    offs+=pack(data+offs,ll,0,no-1,swap_endian);
    offs+=pack(data+offs,lh,0,no-1,swap_endian);
    offs+=pack(data+offs,scattered_light,0,no-1,swap_endian);
    offs+=pack(data+offs,spectral_broadening,0,no-1,swap_endian);
    offs+=pack(data+offs,obs_qs,0,no-1,swap_endian);
    offs+=pack(data+offs,synth_qs,0,no-1,swap_endian);
    offs+=pack(data+offs,no_iterations,0,no-1,swap_endian);
    offs+=pack(data+offs,starting_lambda,0,no-1,swap_endian);
    offs+=pack(data+offs,stopping_chisq,0,no-1,swap_endian);
    for(int o=0;o<no;++o){
      offs+=pack(data+offs,lambda[o],0,nlambda[o]-1,swap_endian);
      offs+=pack(data+offs,weights[o],0,nlambda[o]-1,swap_endian);
      offs+=pack(data+offs,w_stokes[o],0,3,swap_endian);
      offs+=pack(data+offs,name[o]);
    }
  }
  if (nm){
    offs+=pack(data+offs,read_model_from_file,0,nm-1,swap_endian);
    for (int m=0;m<nm;++m)
      if (read_model_from_file[m])
        offs+=pack(data+offs,input_models[m]);
  }

//
  offs+=pack(data+offs,uid,swap_endian);
  offs+=pack(data+offs,gid,swap_endian);
  offs+=pack(data+offs,uname);
//
  offs+=pack(data+offs,cdcl); // transfer compression level
  offs+=pack(data+offs,nmth); // # master threads
  offs+=pack(data+offs,nsth); // # slave threads
  offs+=pack(data+offs,nsln); // # slave nodes
//
  if(offs!=usize) io.msg(IOL_WARN,"jnfo::compress: number of packed bytes was %d but %d bytes allocated\n",offs,usize);
//
  size=usize;
  byte *buf=z_compress(data,size,level,swap_endian,io);
  delete[] data;
  return buf;
}
