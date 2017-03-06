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
//
// initialise gcfg first
//

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
    lambda+=1; // get_numbers returns range from 0...N-1
    for(int l=0;l<nlambda;++l){
      lambda[l]*=1E-8; // A->cm
      lambda[l] = airtovac(lambda[l]);
    }
    io.msg(IOL_XNFO,"obs: %d %E...%E\n",nlambda,lambda[0],lambda[nlambda-1]);
    if(!(name=get_arg(odata,"NAME",0)))
      io.msg(IOL_ERROR|IOL_FATAL,"obs: no output file name specified!\n");
  
    delete[] tmp_str;
  }else io.msg(IOL_ERROR|IOL_FATAL,"obs \"%s\" config: observation must have wavelength specified to be observable!\n",id);
//
  if(char *s=arg_test(odata)) io.msg(IOL_WARN,"obs \"%s\" config: the following lines were not processed:%s\n",id,s);
}

ocfg::~ocfg(void)
{
  if(id) delete[] id;
  if(lambda) delete[] lambda;
  if(name) delete[] name;
}

