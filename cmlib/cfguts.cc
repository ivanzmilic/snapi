#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pwd.h>
#include "types.h"
#include "io.h"

char *clean_path(char *wd)
// Interpret ~name/ and remove ../ combinations...
{
  char *lwd=new char [strlen(wd)+1];
  char *s;

  strcpy(lwd,wd);
  if(lwd[0]=='~'){
    char *hd=0;
    if(lwd[1]=='/'){
      char *q;
      if(q=getenv("HOME")){ // $HOME exists...
        strcpy(hd=new char [strlen(q)+1],q);
      }else{ // $HOME does not exist... Look up user in passwd file
        struct passwd *pwd=getpwuid(geteuid());
        strcpy(hd=new char [strlen((*pwd).pw_dir)+1],(*pwd).pw_dir);
      }
    }else{ // Some other user's home-dir...
      int n=(int)(strchr(lwd,'/')-lwd-1);
      char *user=new char [n+1];
      strncpy(user,lwd+1,n)[n]=0;
      struct passwd *pwd=getpwnam(user);
      strcpy(hd=new char [strlen(pwd->pw_dir)+1],pwd->pw_dir);
      delete[] user;
    }
    if(hd){
      char *ns=new char [strlen(hd)+strlen(strchr(lwd,'/'))+1];
      sprintf(ns,"%s%s",hd,strchr(lwd,'/'));
      delete[] lwd;
      lwd=ns;
    }
  }
  s=lwd;
  // remove double '/'s
  while(s=strstr(lwd,"//")){
    char *q;
    s[0]=0;
    ++s;
    q=new char [strlen(lwd)+strlen(s)+1];
    sprintf(q,"%s%s",lwd,s);
    delete[] lwd;
    lwd=q;
  }
  // remove 'dir/../' occurences
  while(s=strstr(lwd,"/../")){
    char *v=s-1;
    while((v[0]!='/')&&(v>lwd)) --v; // find start or previous directory start
    if(v[0]=='/') ++v; // start of directory name
    if(strncmp(v,"..",s-v)&&strncmp(v,".",s-v)){
      v[0]=0;
      v=new char [strlen(lwd)+strlen(s)+1];
      sprintf(v,"%s%s",lwd,s+strlen("/../"));
      delete[] lwd;
      lwd=v;
    }
  }
  while(s=strstr(lwd,"/./")){
    char *v=new char [strlen(lwd)-1];
    s[1]=0;
    s+=3;
    sprintf(v,"%s%s",lwd,s);
    delete[] lwd;
    lwd=v;
  }
  if(lwd[strlen(lwd)-1]!='/'){
    char *v=new char [strlen(lwd)+2];
    sprintf(v,"%s/",lwd);
    delete[] lwd;
    lwd=v;
  }
  return lwd;
}

char *file2buf(FILE *f)
{
  char *fc=new char [1],s[1000];
  fc[0]=0;
  int theend;
  do theend=(fgets(s,1000,f)==NULL); while((s[0]=='#')&&(!feof(f))&&(!theend));
  do{ // read file in buffer
    char *tfc=new char [strlen(fc)+strlen(s)+1];
    strcpy(tfc,fc);
    strcat(tfc,s);
    delete[] fc;
    fc=tfc;
    do theend=(fgets(s,1000,f)==NULL); while((s[0]=='#')&&(!feof(f))&&(!theend));
  }while(!feof(f));
  return fc;
}

char **scope_sep(char *fc,const char *scope_type,io_class &io)
{
  int ns=0;
  char **scopes=new char* [ns+1];
  scopes[ns]=0;
  while(char *p=strstr(fc,scope_type)){
    char *q=p+strlen(scope_type);
    while(q[0]==' ') ++q; // skip spaces
    if(q[0]!='{') io.msg(IOL_ERROR|IOL_FATAL,"unexpected character in config file (expected \"{\" but found \"%c\")\n",q[0]);
    char *r=q+1;
    int level=1;
    level+=(r[0]=='{');
    level-=(r[0]=='}');
    while(level&&(r[0])){
      ++r;
      level+=(r[0]=='{');
      level-=(r[0]=='}');
    }
    if(!(r[0])) io.msg(IOL_ERROR|IOL_FATAL,"unexpected end of config file (forgot closing braces?)!");
    char **ts=new char* [ns+2];
    memcpy(ts,scopes,ns*sizeof(char*));
    delete[] scopes;
    scopes=ts;
    scopes[ns]=new char [r-q];
    memcpy(scopes[ns],q+1,(r-q-1)*sizeof(char));
    scopes[ns][r-q-1]=0;
    scopes[ns+1]=0;
    ++ns;
    memmove(p,r+1,strlen(r)*sizeof(char));
  }
  return scopes;
}

char **str_sep(const char *str,char sep)
{
  int nf=1;
  char **field=new char* [nf+1];
  field[nf]=0;
  field[nf-1]=strcpy(new char [strlen(str)+1],str);
  while(char *p=strchr(field[nf-1],sep)){
    char **tmp=new char* [nf+2];
    memcpy(tmp,field,nf*sizeof(char*));
    delete[] field;
    field=tmp;
    p[0]=0;
    field[nf]=p+1;
    field[++nf]=0;   
  }
  return field;
}

char *get_arg(char *args,const char *key,int isflag)
{
  if(char *p=strstr(args,key)){          // find keyword
    char *q=p+strlen(key);
    char *r=strchr(q,'\n');              // find end of line
    if(!r) r=q+strlen(q);                // last argument
// copy key+arg
    char *arg=new char [r-q+1];
    memcpy(arg,q,(r-q)*sizeof(char));
    arg[r-q]=0;
//
    if(!(q=strchr(arg,'='))){
      if(isflag){
        arg[0]=0;
        p[0]='*';
	return arg;
      }
      exit(fprintf(stderr,"error: keyword \"%s\" needs argument but is not followed by \"=\" symbol\n",key)); // find assignment operator  
    }
// assignment found
    q+=1;                                // start argument
    while(q[0]==' ') ++q;                // remove starting spaces
    r=q+strlen(q)-1;
    while(r[0]==' ') --r;                // remove trailing spaces
    r[1]=0;
    memmove(arg,q,strlen(q)+1);
    p[0]='*';
    return arg;
  }else return 0; // argument not found
}

void append(int *&in,int nl,int &nh,int num)
{
  int *tmp=new int [nh+1]-nl;
  if(nh-nl+1) memcpy(tmp+nl,in+nl,(nh-nl+1)*sizeof(int));
  delete[] (in+nl);
  in=tmp;
  in[++nh]=num;
}

char *arg_test(char *data)
{
  char *rs=new char [1];
  rs[0]=0;
  char *q=data;
  while(char *p=strchr(q,'\n')){
    if(p-q>1){
      int i=0;
      while(q[i]==' ') ++i;
      if(p-q>i)
        if(q[i]!='*'){
          char *s=new char [strlen(rs)+(p-q)+4];
          memcpy(s,rs,strlen(rs));
          s[strlen(rs)]='\n';
          s[strlen(rs)+1]='"';
          memcpy(s+strlen(rs)+2,q,p-q);
          s[strlen(rs)+(p-q)+2]='"';
          s[strlen(rs)+(p-q)+3]=0;
          delete[] rs;
          rs=s;
        }
    }
    q=p+1;
  }
  if(strlen(q)){
    int i=0;
    while(q[i]==' ') ++i;
    if(strlen(q)>i)
      if(q[i]!='*'){
        char *s=new char [strlen(rs)+strlen(q)+4];
        memcpy(s,rs,strlen(rs));
        s[strlen(rs)]='\n';
        s[strlen(rs)+1]='"';
        strcpy(s+strlen(rs)+2,q);
        s[strlen(rs)+strlen(q)+2]='"';
        s[strlen(rs)+strlen(q)+3]=0;
        delete[] rs;
        rs=s;
      }
  }
  if(!strlen(rs)){
    delete[] rs;
    return 0;
  }else return rs;
}

