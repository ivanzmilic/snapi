#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include "version.h"
#include "types.h"
#include "cmdcfg.h"
#include "cmdln.h"
#include "net.h"
#include "io.h"
#include "uts.h"
#include "pack.h"

struct jinfo{
  int active,id,pri;
  fp_t f,g;
  char *uname,*jname;
  struct tm t_sub;
  jinfo(void);
  ~jinfo(void);
};

jinfo::jinfo(void)
{
  memset(this,0,sizeof(struct jinfo));
}

jinfo::~jinfo(void)
{
  if(uname) delete[] uname;
  if(jname) delete[] jname;
}

struct jinfo *j_unpack(byte *buf,int size,int &n,byte swap_endian)
{
  int offs=0;
  n=0;
  struct jinfo *job=new jinfo [1]-1;
  while(offs<size){
    struct jinfo *t=new jinfo [n+1]-1;
    if(n){
      memcpy(t+1,job+1,n*sizeof(struct jinfo));
      memset(job+1,0,n*sizeof(struct jinfo));
    }
    delete[] (job+1);
    job=t;
    ++n;
//
    offs+=unpack(buf+offs,job[n].active,swap_endian);
    offs+=unpack(buf+offs,job[n].id,swap_endian);
    offs+=unpack(buf+offs,job[n].pri,swap_endian);
    offs+=unpack(buf+offs,job[n].f,swap_endian);
    offs+=unpack(buf+offs,job[n].g,swap_endian);
    offs+=unpack(buf+offs,job[n].uname);
    offs+=unpack(buf+offs,job[n].jname);
    if(!strlen(job[n].jname)){
      delete[] job[n].jname;
      job[n].jname=strcpy(new char [3],"-");
    }
    offs+=unpack(buf+offs,job[n].t_sub.tm_sec,swap_endian);
    offs+=unpack(buf+offs,job[n].t_sub.tm_min,swap_endian);
    offs+=unpack(buf+offs,job[n].t_sub.tm_hour,swap_endian);
    offs+=unpack(buf+offs,job[n].t_sub.tm_year,swap_endian);
    offs+=unpack(buf+offs,job[n].t_sub.tm_mon,swap_endian);
    offs+=unpack(buf+offs,job[n].t_sub.tm_mday,swap_endian);
  }
  return job;
}

struct sinfo{
  int64_t flops;
  int s_id;
  char *name,*version,*os,*machine;
  byte state;
  pid_t pid;
  sinfo(void);
  ~sinfo(void);
};

sinfo::sinfo(void)
{
  memset(this,0,sizeof(struct sinfo));
}

sinfo::~sinfo(void)
{
  if(name) delete[] name;
  if(version) delete[] version;
  if(os) delete[] os;
  if(machine) delete[] machine;
}

struct sinfo *s_unpack(byte *buf,int size,int &n,byte swap_endian)
{
  int offs=0;
  n=0;
  struct sinfo *slv=new sinfo [1]-1;
  while(offs<size){
    struct sinfo *t=new sinfo [n+1]-1;
    if(n){
      memcpy(t+1,slv+1,n*sizeof(struct sinfo));
      memset(slv+1,0,n*sizeof(struct sinfo));
    }
    delete[] (slv+1);
    slv=t;
    ++n;
//
    offs+=unpack(buf+offs,slv[n].s_id,swap_endian);
    offs+=unpack(buf+offs,slv[n].name);
    offs+=unpack(buf+offs,slv[n].pid,swap_endian);
    offs+=unpack(buf+offs,slv[n].version);
    offs+=unpack(buf+offs,slv[n].os);
    offs+=unpack(buf+offs,slv[n].machine);
    offs+=unpack(buf+offs,slv[n].flops,swap_endian);
    offs+=unpack(buf+offs,slv[n].state);  
  }
  return slv;
}

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
    fprintf(stdout,"MOMFBD version: %s\n",VERSION_STR);
    exit(0);
  }
  io_class io(cmd);
//
  char *hostname;
  cmd.get("--master",hostname);
  int port;
  cmd.get("--port",port);
  sock_class sock(hostname,port);
// 
  int jobs;
  if(cmd.get("-j",jobs)){
    sock.send_cmd(CMD_STAT);      // get statistics
    sock.send_cmd(CMD_STAT_JOBQ); // get job statistics
//
    struct id sid;
    sock.send_id(sid);
    struct id mid(&sock);
    byte swap_endian=(sid.endian!=mid.endian);
    sock.set_swap(swap_endian);
//
    int sz,nj;
    byte *data=sock.recv(sz);
    struct jinfo *job=j_unpack(data,sz,nj,swap_endian);
    delete[] data;
    char stt[8]={'Q','P','R','A','W','D','K','K'};
    int i_w=0,id_w=2,u_w=4,jn_w=4,pri_w=3;
    char tmp[1000]; // this better be long enough...
    for(int i=1;i<=nj;++i){
      sprintf(tmp,"%d",i);
      i_w=max(strlen(tmp),i_w);
      sprintf(tmp,"%d",job[i].id);
      id_w=max(strlen(tmp),id_w);
      u_w=max(strlen(job[i].uname),u_w);
      jn_w=max(strlen(job[i].jname),jn_w);
      sprintf(tmp,"%d",job[i].pri);
      pri_w=max(strlen(tmp),pri_w);
    }
// header
    char fill[100];
    sprintf(tmp,"");
    for(int j=0;j<=(i_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
    fill[i_w-strlen(tmp)]=0;
    fprintf(stdout,"%s%s  ",fill,tmp);
//
    sprintf(tmp,"ID");
    for(int j=0;j<=(id_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
    fill[id_w-strlen(tmp)]=0;
    fprintf(stdout,"%s%s  ",fill,tmp);
//
    fprintf(stdout,"    DATE     TIME   ");
//
    sprintf(tmp,"NAME");
    for(int j=0;j<=(jn_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
    fill[jn_w-strlen(tmp)]=0;
    fprintf(stdout,"%s%s  ",fill,tmp);
//
    sprintf(tmp,"USER");
    for(int j=0;j<=(u_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
    fill[u_w-strlen(tmp)]=0;
    fprintf(stdout,"%s%s  ",fill,tmp);
//
    sprintf(tmp,"PRI");
    for(int j=0;j<=(pri_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
    fill[pri_w-strlen(tmp)]=0;
    fprintf(stdout,"%s%s  ",fill,tmp);
//
    fprintf(stdout,"STATUS\n");
//
    for(int i=1;i<=nj;++i){
      sprintf(tmp,"%d",i);
      for(int j=0;j<=(i_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[i_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s: ",fill,tmp);
//
      sprintf(tmp,"%d",job[i].id);
      for(int j=0;j<=(id_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[id_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      fprintf(stdout,"%04d-%02d-%02d %02d:%02d:%02d ",1900+job[i].t_sub.tm_year,1+job[i].t_sub.tm_mon,job[i].t_sub.tm_mday,job[i].t_sub.tm_hour,job[i].t_sub.tm_min,job[i].t_sub.tm_sec);
//
      sprintf(tmp,"%s",job[i].jname);
      for(int j=0;j<=(jn_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[jn_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      sprintf(tmp,"%s",job[i].uname);
      for(int j=0;j<=(u_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[u_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      sprintf(tmp,"%d",job[i].pri);
      for(int j=0;j<=(pri_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[pri_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      fprintf(stdout,"%c",stt[job[i].active]);
//
      if(stt[job[i].active]=='P'){
        const char pps[]={'R','F','C','?','?'};
        int st=(int)(job[i].g-1E-10);
        if((job[i].g>=0)&&(job[i].g<=4))
          fprintf(stdout,"[%c(%4.1f%%)]",pps[st],100.0*(job[i].g-(fp_t)st));
        else
          fprintf(stdout,"[err]");
        
      }else if(stt[job[i].active]=='A')
        fprintf(stdout,"[%4.1f%%]",100.0*job[i].f);
      fprintf(stdout,"\n");
    }
    delete[] (job+1);
  }
  char *slvs;
  if(cmd.get("-s",slvs)){
    sock.send_cmd(CMD_STAT);      // get statistics
    sock.send_cmd(CMD_STAT_SLVQ); // get slave statistics
//
    struct id sid;
    sock.send_id(sid);
    struct id mid(&sock);
    byte swap_endian=(sid.endian!=mid.endian);
    sock.set_swap(swap_endian);
//
    int sz,ns;
    byte *data=sock.recv(sz);
    struct sinfo *slv=s_unpack(data,sz,ns,swap_endian);
    delete[] data;
    int i_w=0,id_w=2,n_w=4,pid_w=3,ver_w=7,os_w=2,m_w=4,f_w=6;
    char tmp[1000]; // this better be long enough...
    int64_t flops=0,GF=1000000000,Gf=100000;
    for(int i=1;i<=ns;++i){
      sprintf(tmp,"%d",i);
      i_w=max(strlen(tmp),i_w);
      sprintf(tmp,"%d",slv[i].s_id);
      id_w=max(strlen(tmp),id_w);
      n_w=max(strlen(slv[i].name),n_w);
      ver_w=max(strlen(slv[i].version),ver_w);
      sprintf(tmp,"%d",slv[i].pid);
      pid_w=max(strlen(tmp),pid_w);
      sprintf(tmp,"%s",slv[i].os);
      os_w=max(strlen(tmp),os_w);
      sprintf(tmp,"%s",slv[i].machine);
      m_w=max(strlen(tmp),m_w);
      if(slv[i].flops>=0){
        flops+=slv[i].flops;
        sprintf(tmp,"%"I64FMT".%04"I64FMT,slv[i].flops/GF,(slv[i].flops%GF)/Gf);
      }else
        sprintf(tmp,"N/A");
      f_w=max(strlen(tmp),f_w);
    }
    fprintf(stdout,"%d slaves, %"I64FMT".%04"I64FMT" GFlops total\n",ns,flops/GF,(flops%GF)/Gf);
    char fill[100];
//
    memset(fill,' ',i_w);
    fill[i_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    sprintf(tmp,"ID");
    memset(fill,' ',id_w);
    memcpy(fill+(id_w-strlen(tmp))/2,tmp,strlen(tmp));
    fill[id_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    sprintf(tmp,"NAME");
    memset(fill,' ',n_w);
    memcpy(fill+(n_w-strlen(tmp))/2,tmp,strlen(tmp));
    fill[n_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    sprintf(tmp,"PID");
    memset(fill,' ',pid_w);
    memcpy(fill+(pid_w-strlen(tmp))/2,tmp,strlen(tmp));
    fill[pid_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    fprintf(stdout,"STATE  ");
//
    sprintf(tmp,"VERSION");
    memset(fill,' ',ver_w);
    memcpy(fill+(ver_w-strlen(tmp))/2,tmp,strlen(tmp));
    fill[ver_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    sprintf(tmp,"ARCH");
    memset(fill,' ',m_w);
    memcpy(fill+(m_w-strlen(tmp))/2,tmp,strlen(tmp));
    fill[m_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    sprintf(tmp,"OS");
    memset(fill,' ',os_w);
    memcpy(fill+(os_w-strlen(tmp))/2,tmp,strlen(tmp));
    fill[os_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    sprintf(tmp,"GFLOPS");
    memset(fill,' ',f_w);
    memcpy(fill+(f_w-strlen(tmp))/2,tmp,strlen(tmp));
    fill[f_w]=0;
    fprintf(stdout,"%s  ",fill);
//
    fprintf(stdout,"\n");
    for(int i=1;i<=ns;++i){
      sprintf(tmp,"%d",i);
      for(int j=0;j<=(i_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[i_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s: ",fill,tmp);
//
      sprintf(tmp,"%d",slv[i].s_id);
      for(int j=0;j<=(id_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[id_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",tmp,fill);
//
      sprintf(tmp,"%s",slv[i].name);
      for(int j=0;j<=(n_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[n_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",tmp,fill);
//
      sprintf(tmp,"%d",slv[i].pid);
      for(int j=0;j<=(pid_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[pid_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      fprintf(stdout,"  %c    ",(char)(slv[i].state));
//
      sprintf(tmp,"%s",slv[i].version);
      for(int j=0;j<=(ver_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[ver_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      sprintf(tmp,"%s",slv[i].machine);
      for(int j=0;j<=(m_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[m_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      sprintf(tmp,"%s",slv[i].os);
      for(int j=0;j<=(os_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
      fill[os_w-strlen(tmp)]=0;
      fprintf(stdout,"%s%s  ",fill,tmp);
//
      if(slv[i].flops>=0){
        sprintf(tmp,"%"I64FMT".%04"I64FMT,slv[i].flops/GF,(slv[i].flops%GF)/Gf);
        for(int j=0;j<=(f_w-(int)strlen(tmp)-1);++j) fill[j]=' ';
        fill[f_w-strlen(tmp)]=0;
        fprintf(stdout,"%s%s",fill,tmp);
      }else{
        sprintf(tmp,"N/A");
        memset(fill,' ',f_w);
        memcpy(fill+(f_w-strlen(tmp))/2,tmp,strlen(tmp));
        fill[f_w]=0;
        fprintf(stdout,"%s",fill);
      }
//
      fprintf(stdout,"\n");
    }
    delete[] (slv+1);
  }
  return 0;
}
