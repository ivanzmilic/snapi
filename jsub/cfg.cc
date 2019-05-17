#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pwd.h>
#include <unistd.h>
#include "types.h"
//#include "const.h"
//#include "conf.h"
#include "io.h"
#include "uts.h"
#include "mem.h"
//#include "ana_io.h"
//#include "fileio.h"
//#include "imgtools.h"
//#include "conf.h"
#include "cfguts.h"
#include "struts.h"
#include "mathtools.h"

#include "cfg.h"

gcfg::gcfg(cmdln &cmd,io_class &io)
{
  memset(this,0,sizeof(gcfg));
  char *filename;
  cmd.get("-cfg",filename);
  struct stat stat_buf;
  if(stat(filename,&stat_buf)<0) io.msg(IOL_ERROR|IOL_FATAL,"config file \"%s\" not found!\n",filename);
  FILE *f=fopen(filename,"r");
  char *fc=file2buf(f);
  fclose(f);
  char **obss=scope_sep(fc,"obs",io);
  for(no=0;obss[no];++no);
  while(strlen(fc)&&(fc[strlen(fc)-1]=='\n')) fc[strlen(fc)-1]=0;
  char **atms=scope_sep(fc,"atmos",io);
  for(na=0;atms[na];++na);
  while(strlen(fc)&&(fc[strlen(fc)-1]=='\n')) fc[strlen(fc)-1]=0;
  char **mods=scope_sep(fc,"model",io);
  for(nm=0;mods[nm];++nm);
  while(strlen(fc)&&(fc[strlen(fc)-1]=='\n')) fc[strlen(fc)-1]=0;
//
  if(char *s=arg_test(fc)) io.msg(IOL_WARN,"global config: the following lines were not processed:%s\n",s);
  delete[] fc;
//
  atm=new acfg* [na];
  for(int a=0;a<na;++a) atm[a]=new acfg(atms[a],io);
  del_str(atms);
//
  obs=new ocfg* [no];
  for(int o=0;o<no;++o) obs[o]=new ocfg(obss[o],*this,io);
  del_str(obss);

  if (mods[0]){
    mod=new mcfg* [nm];
    for(int m=0;m<nm;++m) mod[m]=new mcfg(mods[m],*this,io);
    del_str(mods);
  }else io.msg(IOL_WARN,"global config: no model specified.\n");
}

gcfg::~gcfg(void)
{
  if(atm){
    for(int a=0;a<na;++a) if(atm[a]) delete atm[a];
    delete[] atm;
  }
  if(obs){
    for(int o=0;o<no;++o) if(obs[o]) delete obs[o];
    delete[] obs;
  }
}

ocfg::ocfg(char *odata,struct gcfg &gdata,io_class &io)
{
  memset(this,0,sizeof(ocfg));
  while(strlen(odata)&&(odata[strlen(odata)-1]=='\n')) odata[strlen(odata)-1]=0;
//
  if(!(id=get_arg(odata,"ID",0))) io.msg(IOL_ERROR|IOL_FATAL,"obs config: observable has no ID\n");

  fp_t tempfp;
  to_invert = 0; //default
  if(char *tmp_str=get_arg(odata,"INVERT",0)){
    get_number(tmp_str,tempfp);
    delete[] tmp_str;
    if (tempfp) to_invert = 1;
  }
  return_model = 0;
  if(char *tmp_str=get_arg(odata,"RETURN_MODEL",0)){
    get_number(tmp_str,tempfp);
    delete[] tmp_str;
    if (tempfp) return_model = 1;
  }
  return_atmos = 0;
  if(char *tmp_str=get_arg(odata,"RETURN_ATMOS",0)){
    get_number(tmp_str,tempfp);
    delete[] tmp_str;
    if (tempfp) return_atmos = 1;
  }
  //
  if (to_invert == 1){ // read inversion parameters only if we are in the invert mode
    if(char *tmp_str=get_arg(odata,"SCATTERED_LIGHT",0)){
      get_number(tmp_str,scattered_light);
      delete[] tmp_str;
    }else scattered_light = 0.0;//default
    if(char *tmp_str=get_arg(odata,"SPECTRAL_BROADENING",0)){
      get_number(tmp_str,spectral_broadening);
      delete[] tmp_str;
      spectral_broadening *= 1E-11;//convert to cm from mA
      spectral_broadening /= 2.35; // convert from FWHM to sigma
    }else spectral_broadening = 0;
    if(char *tmp_str=get_arg(odata,"OBSERVED_CONTINUUM",0)){
      get_number(tmp_str,obs_qs);
      delete[] tmp_str;
    } else io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: error extracting observed continuum level\n",id);
    if(char *tmp_str=get_arg(odata,"CGS_CONTINUUM",0)){
      get_number(tmp_str,synth_qs);
      delete[] tmp_str;
    } else io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: error extracting calculated continuum level\n",id);
    if(char *tmp_str=get_arg(odata,"ITERATIONS",0)){
      fp_t n_iter; // temp
      get_number(tmp_str,n_iter);
      no_iterations = int(n_iter);
      delete[] tmp_str;
    } else no_iterations = 10;
    if(char *tmp_str=get_arg(odata,"STARTING_LM",0)){
      get_number(tmp_str,starting_lambda);
      delete[] tmp_str;
    } else starting_lambda=1E3;
    int m = 0;
    if(char *val_str=get_arg(odata,"STOKES_WEIGHTS",0)){
      if(get_numbers(val_str,w_stokes,m)<0){
        io.msg(IOL_ERROR|IOL_FATAL,"ocfg::ocfg: failed to convert STOKES_WEIGTHS argument \"%s\" to floating point values\n",val_str);
        m=0;
        delete[]w_stokes;
      }
      if(m!=4){ 
        io.msg(IOL_ERROR|IOL_FATAL,"ocfg::ocfg: number of stokes weights differs from 4. Setting them to 1,0,0,1");
        m=0;
        delete[]w_stokes;
      }
    }
    else
      m =0;
    if (m==0){
      w_stokes = new fp_t [4]-1;
      w_stokes[1] = w_stokes[4] = 1.0;
      w_stokes[2] = w_stokes[3] = 0.0;
    }
    w_stokes+=1;
  }
  else {
    w_stokes = new fp_t [4]-1;
    memset(w_stokes+1,0,4*sizeof(fp_t));
    w_stokes+=1;
  }
  //
  if(char *tmp_str=get_arg(odata,"AZ",0)){
    get_number(tmp_str,az);
    delete[] tmp_str;
  }else az=0.0; // default
  if(char *tmp_str=get_arg(odata,"EL",0)){
    get_number(tmp_str,el);
    delete[] tmp_str;
  }else el=0.0; // default
  if(char *tmp_str=get_arg(odata,"LAMBDA",0)){
    if(char *dl_str=get_arg(odata,"DLAMBDA",0)){
      fp_t dl;
      if(get_number(dl_str,dl)<0)
        io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: error extracting wavelength step from \"%s\"!\n",id,dl_str);
      delete[] dl_str;
      if(get_numbers(tmp_str,lambda,nlambda,dl)<0)
        io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: error extracting wavelength points from \"\" using step size %f!\n",id,tmp_str,dl);
    }else{
      if(get_numbers(tmp_str,lambda,nlambda)<0)
        io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: error extracting wavelength points (did you specify a range but not a step size?)!\n",id);
    }
    lambda+=1; // get_numbers returns
    for(int l=0;l<nlambda;++l){
      lambda[l]*=1E-8; // A->cm
    }
    
    io.msg(IOL_XNFO,"obs: %d %E...%E\n",nlambda,lambda[0],lambda[nlambda-1]);
    delete[] tmp_str;
  }
  else if(char *tmp_str2=get_arg(odata,"LGRID",0)){
    FILE * tmpinpt = fopen(tmp_str2,"r");
    nlambda = 0;
    if (fscanf(tmpinpt,"%d", &nlambda)!=EOF){
      lambda = new fp_t [nlambda];
    }
    else lambda = 0; // SHOULD BE ERROR
    fp_t tmp;
    for (int i=0;i<nlambda;++i)
      if (fscanf(tmpinpt,"%lf",&tmp)!=EOF){
        lambda[i] = tmp*1E-8; // A -> cm
    }
    // should be else - error  
  fclose(tmpinpt);
  delete[] tmp_str2;
  //for (int i=0;i<nlambda;++i)
  //    fprintf(stderr,"%d %f \n",i,lambda[i]*1E8);
  }else io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: observation must have wavelength specified to be observable!\n",id);
  if(!(name=get_arg(odata,"NAME",0)))
      io.msg(IOL_ERROR|IOL_FATAL,"obs: no output file name specified!\n");
  
  if(char *tmp_str=get_arg(odata,"XRANGE",0)){
    fp_t * xrange;
    int temp=0;
    if(get_numbers(tmp_str,xrange,temp)<0) io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: failed to convert VALUE argument \"%s\" to floating point values\n",id,tmp_str);
    xrange +=1;
    xl = int(xrange[0]);xh=int(xrange[1]);
    delete[] tmp_str;
    delete[] xrange;
  }else xl=xh=0; // default
  if(char *tmp_str=get_arg(odata,"YRANGE",0)){
    fp_t * yrange;
    int temp=0;
    if(get_numbers(tmp_str,yrange,temp)<0) io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: failed to convert VALUE argument \"%s\" to floating point values\n",id,tmp_str);
    yrange +=1;
    yl = int(yrange[0]);yh=int(yrange[1]);
    delete[] tmp_str;
    delete[] yrange;
  }else yl=yh=0; // default
  if(char *tmp_str=get_arg(odata,"LRANGE",0)){
    fp_t * lrange;
    int temp=0;
    if(get_numbers(tmp_str,lrange,temp)<0) io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: failed to convert VALUE argument \"%s\" to floating point values\n",id,tmp_str);
    lrange +=1;
    ll = int(lrange[0]);lh=int(lrange[1]);
    delete[] tmp_str;
    delete[] lrange;
  }else ll=lh=0; // default
  if(char *tmp_str=get_arg(odata,"MASK",0)){
    FILE * tmpinpt = fopen(tmp_str,"r");
    int N_points = 0;
    if (fscanf(tmpinpt,"%d", &N_points)!=EOF)
      weight = new fp_t [N_points];
    else weight = 0;
    if (weight){
      fp_t tmp;
      for (int i=0;i<N_points;++i)
        if (fscanf(tmpinpt,"%lf",&tmp)!=EOF){
          weight[i] = tmp;
        }
        else {
          io.msg(IOL_WARN,"obs \"%s\" config: could not read wavelength mask propely. Using 1 everywhere\n",id);
          delete[](weight);
          weight = 0;
          break;
        }
    }
    fclose(tmpinpt);
    delete[] tmp_str;
  }else weight = 0;

//
  if(char *s=arg_test(odata)) io.msg(IOL_WARN,"obs \"%s\" config: the following lines were not processed:%s\n",id,s);
}

ocfg::~ocfg(void)
{
  if(id) delete[] id;
  if(lambda) delete[] lambda;
  if(name) delete[] name;
  if(weight) delete[]weight;
  delete []w_stokes;
}

mcfg::mcfg(char *mdata,struct gcfg &gdata,io_class &io)
{
  //memset(this,0,sizeof(mcfg)); // Why is this necessary?
  char **pars=scope_sep(mdata,"parameter",io);
  while(strlen(mdata)&&(mdata[strlen(mdata)-1]=='\n')) mdata[strlen(mdata)-1]=0;
//
  if(!(id=get_arg(mdata,"ID",0))) io.msg(IOL_ERROR|IOL_FATAL,"mod config: model has no ID\n");
//
  if(char *tmp_str=get_arg(mdata,"READ_FROM_FILE",0)){
    get_number(tmp_str,read_from_file);
    delete[] tmp_str;
  }else read_from_file=0; // default
  if(char *tmp_str=get_arg(mdata,"TAUGRID",0)){
    fp_t * taurange;
    int temp=0;
    if(get_numbers(tmp_str,taurange,temp)<0) io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: failed to convert VALUE argument \"%s\" to floating point values\n",id,tmp_str);
    taurange +=1;
    tau_min = taurange[0];
    tau_max = taurange[1];
    delete[] tmp_str;
    delete[] taurange;
  }else{
    tau_min=-5.0;
    tau_max=1.0;
  }
  if (read_from_file){
    if(!(filename=get_arg(mdata,"FILENAME",0))){ 
      io.msg(IOL_ERROR|IOL_FATAL,"mod config: filename needed but not specified\n");
      filename = 0;
    }
  }
  if (!read_from_file)
    filename=0;

  np=0;
  par=0;
  if(pars){
    for(np=0;pars[np];++np);
    par=new parcfg* [np];
    for(int p=0;p<np;++p) par[p]=new parcfg(pars[p],io);
    del_str(pars);
  }
 
  if(char *s=arg_test(mdata)) io.msg(IOL_WARN,"mod \"%s\" config: the following lines were not processed:%s\n",id,s);
}

mcfg::~mcfg(void)
{
  if(id) delete[] id;
  if (par){
    for (int p=0;p<np;++p)
      if (par[p]) delete par[p];
    delete []par;
  }
  if (filename) delete[] filename;

}

parcfg::parcfg(char* pardata, io_class &io){

  if(!(id=get_arg(pardata,"ID",0))) io.msg(IOL_ERROR|IOL_FATAL,"parcfg::parcfg: no parameter type specified!\n");

  if(char *tau_str=get_arg(pardata,"TAU",0)){
    //printf("%s\n", tau_str);
    if(get_numbers(tau_str,tau,n)<0) io.msg(IOL_ERROR|IOL_FATAL,"parcfg::parcfg: failed to convert TAU argument \"%s\" to floating point values\n",tau_str);
    tau+=1;
    delete[] tau_str;
  }else io.msg(IOL_ERROR|IOL_FATAL,"parcfg::parcfg: no tau specified for parameter nodes!\n");
  int m =0;
  if(char *val_str=get_arg(pardata,"VALUE",0)){
    if(get_numbers(val_str,value,m)<0) io.msg(IOL_ERROR|IOL_FATAL,"parcfg::parcfg: failed to convert VALUE argument \"%s\" to floating point values\n",val_str);
    value+=1;
    delete[] val_str;  
  }else io.msg(IOL_ERROR|IOL_FATAL,"parcfg::parcfg: no value specified for parameter nodes!\n");
  if(char *tmp_str=get_arg(pardata,"REG_TYPE",0)){
    get_number(tmp_str,reg_type);
    delete[] tmp_str;
  }else reg_type=0; // default
  if(char *tmp_str=get_arg(pardata,"REG_ALPHA",0)){
    get_number(tmp_str,reg_alpha);
    delete[] tmp_str;
  }else reg_alpha=0.0; // default

  if (strcmp(id,"THETA") == 0 || strcmp(id,"PHI") == 0)
    for (int nn=0;nn<n;++nn)
      value[nn] *= 3.141592653589793238462/180.0; // convert to radians

  //printf("%s \n", id);
  //for (int nn=0;nn<n;++nn)
    //printf("%e %e \n", tau[nn],value[nn]);
//
}

parcfg::~parcfg(void){
  if (id) delete[]id;
  if (tau) delete[]tau;
  if (value) delete[]value;
}

