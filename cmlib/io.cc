#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/socket.h>
#include <errno.h>

#include "cmdln.h"
#include "uts.h"
#include "pack.h"
#include "net.h"
#include "io.h"

io_class::io_class(void)
{
  level_mask=~(0xFFFFFFFF<<2);
  f=stderr;
  sock=0;
  ssock=0;
  swap_endian=0;
  ts=0;
}

io_class::io_class(int sock_in,byte swap_endian_in,byte iocmd_in,cmdln &cmd)
{
  int v,q;
  cmd.get("-v",v);
  cmd.get("-q",q);
  level_mask=~(0xFFFFFFFF<<(verb=max(v-q+2,0)));
  sock=sock_in;
  ssock=0;
  f=0;
  iocmd=iocmd_in;
  swap_endian=swap_endian_in;
  cmd.get("-t",ts);
}

io_class::io_class(sock_class &sock_in,byte swap_endian_in,byte iocmd_in,int level,int ts_in)
{
  level_mask=~(0xFFFFFFFF<<(verb=max(level,0)));
  ssock=&sock_in;
  sock=0;
  f=0;
  iocmd=iocmd_in;
  swap_endian=swap_endian_in;
  ts=ts_in;
}

io_class::io_class(int sock_in,byte swap_endian_in,byte iocmd_in,int level,int ts_in)
{
  level_mask=~(0xFFFFFFFF<<(verb=max(level,0)));
  sock=sock_in;
  ssock=0;
  f=0;
  iocmd=iocmd_in;
  swap_endian=swap_endian_in;
  ts=ts_in;
}

io_class::io_class(cmdln &cmd)
{
  int v,q;
  cmd.get("-v",v);
  cmd.get("-q",q);
  level_mask=~(0xFFFFFFFF<<(verb=max(v-q+2,0)));
  f=stderr;
  sock=0;
  ssock=0;
  swap_endian=0;
  cmd.get("-t",ts);
}

io_class::~io_class(void)
{
  if((f)&&(f!=stdout)) fclose(f);
}

io_class::io_class(char *wd,char *&fname,const char *jid,uid_t uid,gid_t gid,int level,int ts_in)
{
  level_mask=~(0xFFFFFFFF<<(verb=max(level,0)));
  if(strlen(fname)){
    char *path=new char [strlen(wd)+strlen(fname)+2];
    if(wd[strlen(wd)-1]=='/')
      sprintf(path,"%s%s",wd,fname);
    else
      sprintf(path,"%s/%s",wd,fname);
    if(!(f=fopen(path,"w")))
      fprintf(stderr,"error opening file %s: %s\n",path,strerror(errno));
    if(chown(path,uid,gid)<0)
      fprintf(stderr,"error changing file permissions for %s: %s\n",path,strerror(errno));
    delete[] path;
  }else{
    delete[] fname;
    fname=new char [strlen("invlog.")+strlen(jid)+1];
    sprintf(fname,"invlog.%s",jid);
    char *path=new char [strlen(wd)+strlen(fname)+2];
    if(wd[strlen(wd)-1]=='/')
      sprintf(path,"%s%s",wd,fname);
    else
      sprintf(path,"%s/%s",wd,fname);
    if(f=fopen(path,"w")){
      if(chown(path,uid,gid)<0)
        fprintf(stderr,"error changing file permissions for %s: %s\n",path,strerror(errno));
    }else
      fprintf(stderr,"error opening file %s: %s\n",path,strerror(errno));
    delete[] path;
  }
  ts=ts_in;
}

io_class::io_class(char *wd,char *&fname,const char *jid,int level,int ts_in)
{
  level_mask=~(0xFFFFFFFF<<(verb=max(level,0)));
  if(strlen(fname)){
    char *path=new char [strlen(wd)+strlen(fname)+2];
    if(wd[strlen(wd)-1]=='/')
      sprintf(path,"%s%s",wd,fname);
    else
      sprintf(path,"%s/%s",wd,fname);
    f=fopen(path,"w");
    delete[] path;
  }else{
    delete[] fname;
    fname=new char [strlen("invlog.")+strlen(jid)+1];
    sprintf(fname,"invlog.%s",jid);
    char *path=new char [strlen(wd)+strlen(fname)+2];
    if(wd[strlen(wd)-1]=='/')
      sprintf(path,"%s%s",wd,fname);
    else
      sprintf(path,"%s/%s",wd,fname);
    f=fopen(path,"w");
    delete[] path;
  }
  ts=ts_in;
}

void io_class::reconf(int sock_in,byte iocmd_in,int level,int ts_in)
{
  level_mask=~(0xFFFFFFFF<<(verb=max(level,0)));
  sock=sock_in;
  f=0;
  iocmd=iocmd_in;
  ts=ts_in;
}

void io_class::reconf(int level)
{
  level_mask=~(0xFFFFFFFF<<(verb=max(level,0)));
}

int io_class::msg(int type,const char *formatstring,...)
{
  const char *message[]={"","error","warning","info","extra info","debug level 1","debug level 2","",""};
  int msg_type=type&level_mask;
  int level=0;
  while(msg_type>>level) ++level;
//
  if(msg_type){
    va_list ap;
    va_start(ap,formatstring);
    char *newformatstring;
    if(!(type&IOL_NOID)){
      if(ts){
        struct tm t;
        time_t now=time(0);
        localtime_r(&now,&t);
        newformatstring=new char [strlen(formatstring)+strlen(message[level])+23];
        sprintf(newformatstring,"%04d.%02d.%02d %02d:%02d:%02d %s: %s",t.tm_year+1900,t.tm_mon+1,t.tm_mday,t.tm_hour,t.tm_min,t.tm_sec,message[level],formatstring);
      }else{
        newformatstring=new char [strlen(formatstring)+strlen(message[level])+3];
        sprintf(newformatstring,"%s: %s",message[level],formatstring);
      }
    }else newformatstring=strcpy(new char [strlen(formatstring)+1],formatstring);
    if(f){
      vfprintf(f,newformatstring,ap);
      fflush(f);
    }else{
      char *str=new char [10000];    // fixed length: not so nice...
      vsprintf(str,newformatstring,ap);
      int size=strlen(str)+1+sizeof(int),offs=0;
      byte *buf=new byte [size];
      offs+=pack(buf+offs,type,swap_endian);
      offs+=pack(buf+offs,str);
      delete[] str;
      if(ssock){
        ssock->send_cmd(iocmd);
        ssock->send(buf,size);
      }else{
        ::send(sock,&iocmd,1,0);
        ::send(sock,&size,sizeof(size),0);
        ::send(sock,buf,size,0);
      }
      delete[] buf;
    }
    delete[] newformatstring;
    va_end(ap);
    if(type&IOL_FATAL) exit(1);
    if(msg_type&IOL_ERROR) return -1;
  }
  return 0;
}
