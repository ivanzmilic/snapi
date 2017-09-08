#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <errno.h>
#include "version.h"
#include "types.h"
#include "cmdcfg.h"
#include "cmdln.h"
#include "cfg.h"
//#include "config.h"
#include "net.h"
#include "uts.h"
#include "pack.h"
#include "jnfo.h"
#include "mem.h"

int main(int argc,char *argv[])
{
// setup recognised options (define defaults, then redefine with command line options)
  const char *options[]=CMDOPTS,*maps[]=CMDMAPS;    // recognised options
  char *hlp;
  cmdln cmd(argc,argv,options,maps);    // parse command line
  if(cmd.get("--help",hlp)){
    fprintf(stdout,"%s\n",(cmd.get(hlp))?cmd.help(hlp):cmd.help("--help"));
    exit(0);
  }
  int version;
  cmd.get("--version",version);
  if(version){
    fprintf(stdout,"INVERT version: %s\n",VERSION_STR);
    exit(0);
  }
  io_class io(cmd);
  struct gcfg cfg(cmd,io);

//  io.msg(IOL_INFO,"%d wavelengths\n",cfg.obs[0]->nlambda);
//  for(int i=0;i<cfg.obs[0]->nlambda;++i) io.msg(IOL_INFO,"%f\n",cfg.obs[0]->lambda[i]);

//
  int pri;
  cmd.get("-pri",pri);
  int v,q;
  cmd.get("-v",v);
  cmd.get("-q",q);
  int verb_level=max(v-q+2,0);
  int time_stamping;
  cmd.get("-t",time_stamping);
  
  char wd[10000];
  if(getcwd(wd,10000)==NULL)
    io.msg(IOL_WARN,"failed to get current working directory: %s\n",strerror(errno));
  
  char *lgfname;
  if(!cmd.get("-lg",lgfname)) lgfname=strcpy(new char [1],"");
  char *jname;
  if(!cmd.get("-name",jname)) jname=strcpy(new char [1],"");
//
  struct jnfo ji;
// configure atmospheric structures
  ji.na=cfg.na;
  ji.atmos=new atmosphere* [ji.na];
  for(int a=0;a<ji.na;++a) ji.atmos[a]=atmos_new(cfg.atm[a],io);
//
  ji.no=cfg.no;
  ji.az=new fp_t [ji.no];
  ji.el=new fp_t [ji.no];
  ji.nlambda=new int [ji.no];
  ji.to_invert = new int [ji.no];
  ji.return_model = new int [ji.no];
  ji.return_atmos = new int [ji.no];
  ji.xl = new int [ji.no]; ji.xh = new int [ji.no];
  ji.yl = new int [ji.no]; ji.yh = new int [ji.no];
  ji.ll = new int [ji.no]; ji.lh = new int [ji.no];
  ji.scattered_light = new fp_t [ji.no];
  ji.spectral_broadening = new fp_t [ji.no];
  ji.obs_qs = new fp_t [ji.no];
  ji.synth_qs = new fp_t [ji.no];
  ji.lambda=new fp_t* [ji.no];
  ji.weights=new fp_t* [ji.no];
  ji.name=new char* [ji.no];
  for(int o=0;o<ji.no;++o){
    ji.az[o]=cfg.obs[o]->az;
    ji.el[o]=cfg.obs[o]->el;
    ji.nlambda[o]=cfg.obs[o]->nlambda;
    ji.to_invert[o]=cfg.obs[o]->to_invert;
    ji.return_model[o]=cfg.obs[o]->return_model;
    ji.return_atmos[o]=cfg.obs[o]->return_atmos;
    ji.xl[o] = cfg.obs[o]->xl;ji.xh[o] = cfg.obs[o]->xh;
    ji.yl[o] = cfg.obs[o]->yl;ji.yh[o] = cfg.obs[o]->yh;
    ji.ll[o] = cfg.obs[o]->ll;ji.lh[o] = cfg.obs[o]->lh;
    if (ji.to_invert[o]){
      ji.scattered_light[o] = cfg.obs[o]->scattered_light;
      ji.spectral_broadening[o] = cfg.obs[o]->spectral_broadening;
      ji.obs_qs[o] = cfg.obs[o]->obs_qs;
      ji.synth_qs[o] = cfg.obs[o]->synth_qs;
    }
    ji.lambda[o]=new fp_t [ji.nlambda[o]];
    memcpy(ji.lambda[o],cfg.obs[o]->lambda,ji.nlambda[o]*sizeof(fp_t));
    ji.weights[o] = new fp_t [ji.nlambda[o]];
    if (cfg.obs[o]->weight)
      memcpy(ji.weights[o],cfg.obs[o]->weight,ji.nlambda[o]*sizeof(fp_t));
    else 
      for (int i=0;i<ji.nlambda[o];++i) ji.weights[o][i] = 1.0;

    ji.name[o]=strcpy(new char [strlen(cfg.obs[o]->name)+1],cfg.obs[o]->name);
  }
//
  // Now do the same for the model as well 
  ji.nm = cfg.nm;
  ji.models = new model* [ji.nm];
  for (int m=0;m<ji.nm;++m) ji.models[m]=model_new(cfg.mod[m],io);

  ji.read_model_from_file = new int [ji.nm];
  ji.input_models = new char* [ji.nm];
  for (int m=0;m<ji.nm;++m){ 
    ji.read_model_from_file[m] = cfg.mod[m]->read_from_file;
    if (ji.read_model_from_file[m])
      ji.input_models[m]=strcpy(new char [strlen(cfg.mod[m]->filename)+1],cfg.mod[m]->filename);
    else 
      ji.input_models[m]=0;

  }
//
  ji.uid=getuid();
  ji.gid=getgid();
//
  int bsz=sysconf(_SC_GETPW_R_SIZE_MAX);
  char *buf=new char [bsz];
  struct passwd  p,*pwd;
  if(getpwuid_r(ji.uid,&p,buf,bsz,&pwd)==0)
    ji.uname=strcpy(new char [strlen(pwd->pw_name)+1],pwd->pw_name);
  else
    ji.uname=strcpy(new char [8],"unknown");
  delete[] buf;
//
  int tmp;
  cmd.get("--chunk_compress",tmp);
  ji.cdcl=tmp;
  cmd.get("--master_threads",tmp);
  ji.nmth=tmp;
  cmd.get("--slave_threads",tmp);
  ji.nsth=tmp;
  cmd.get("--num_slaves",tmp);
  ji.nsln=tmp;
// connect
  char *hostname;
  cmd.get("--master",hostname);
  int port;
  cmd.get("--port",port);
  sock_class sock(hostname,port);
  sock.send_cmd(CMD_NEW_JOB); // new job
// handshake
  struct id sid;
  sock.send_id(sid);
  struct id mid(&sock);
  byte swap_endian=(mid.endian!=sid.endian);
  sock.set_swap(swap_endian);
// compress data
  int level;
  cmd.get("-cl",level);
//
  int isz;
  byte *idta=ji.compress(isz,level,swap_endian,io);
  sock.send(idta,isz);        // job info
  delete[] idta;
//
  int asz=3*sizeof(int)+strlen(wd)+strlen(lgfname)+strlen(jname)+3;
  byte *adta=new byte [asz];
  int offs=0;
  offs+=pack(adta+offs,pri,swap_endian);
  offs+=pack(adta+offs,verb_level,swap_endian);
  offs+=pack(adta+offs,time_stamping,swap_endian);
  offs+=pack(adta+offs,wd);
  offs+=pack(adta+offs,lgfname);
  offs+=pack(adta+offs,jname);
  sock.send(adta,asz);
  delete[] adta;
//
  delete[] lgfname;
  delete[] jname;
//
  return 0;
}
