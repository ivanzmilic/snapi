#include <string.h>

#include "types.h"
#include "io.h"
#include "uts.h"
#include "cfguts.h"
#include "struts.h"

#include "molcfg.h"

eqcfg::eqcfg(char *eqdata,const char *name,const char *id,io_class &io)
{
  memset(this,0,sizeof(eqcfg));
  if(type=get_arg(eqdata,"TYPE",0)){
    if(!strcmp(type,"CONST")){
      if(char *kstr=get_arg(eqdata,"K",0)){
        if(get_number(kstr,K)<0) io.msg(IOL_ERROR|IOL_FATAL,"eqcfg::eqcfg: failed to convert constant argument \"%s\" for molecule %s[%s] to floating point value\n",kstr,name,id);
        delete[] kstr;
      }else io.msg(IOL_ERROR|IOL_FATAL,"eqcfg::eqcfg: no equilibrium constant specified for molecule %s[%s]\n",name,id);
    }else if(!strcmp(type,"POLY")){
      if(char *cfstr=get_arg(eqdata,"COEFS",0)){
        int nc;
        if(get_numbers(cfstr,cfs,nc)<0) io.msg(IOL_ERROR|IOL_FATAL,"eqcfg::eqcfg: failed to convert coefficient argument \"%s\" for molecule %s[%s] to floating point values\n",cfstr,name,id);
        cfs+=1;
        ncfs=nc;
        delete[] cfstr;
      }else io.msg(IOL_ERROR|IOL_FATAL,"eqcfg::eqcfg: no equilibrium constant specified for molecule %s[%s]\n",name,id);
    
//    }else if(!strcmp(tmpstr,"SAHA")){
    
    }else if(!strcmp(type,"HMIN")){
    }else io.msg(IOL_ERROR|IOL_FATAL,"eqcfg::eqcfg: equilibrium constant type \"%s \" for molecule %s[%s] not recognised\n",type,name,id);
  }else io.msg(IOL_ERROR|IOL_FATAL,"eqcfg::eqcfg: no equilibrium constant type specified for molecule %s[%s]\n",name,id);
//
  if(char *s=arg_test(eqdata)) io.msg(IOL_WARN,"eqcfg::eqcfg: the following lines were not processed:%s\n",s);
}

eqcfg::~eqcfg(void)
{
  if(type) delete[] type;
  if(cfs) delete[] cfs;
}

molcfg::molcfg(char *mdata,io_class &io)
{
  char **ks=scope_sep(mdata,"eqc",io);
  while(strlen(mdata)&&(mdata[strlen(mdata)-1]=='\n')) mdata[strlen(mdata)-1]=0;
//
  if(!(name=get_arg(mdata,"NAME",0))) io.msg(IOL_ERROR|IOL_FATAL,"no name specified for atom!\n");
  if(!(id=get_arg(mdata,"ID",0))) io.msg(IOL_ERROR|IOL_FATAL,"no id specified for atom %s\n",name);
//
  if(char *tmpstr=get_arg(mdata,"COMPONENTS",0)){
    zid=cid=0;
    char **s=string2strings(tmpstr,io);
    for(int i=0;s[i];++i){
      if(char *p=strchr(s[i],'[')){
        if(char *q=strchr(p+1,']')){
          q[0]=0;
          zid=append(zid,p+1);
        }else  io.msg(IOL_ERROR|IOL_FATAL,"molcfg::molcfg: parse error molecule %s, component %d: %s\n",name,i,s[i]);
        p[0]=0;
      }else if((!strcmp(s[i],"H+"))||(!strcmp(s[i],"p"))){
        zid=append(zid,"II");
        s[i][0]='H';
        s[i][1]=0;
      }else if(!strcmp(s[i],"e")){
        zid=append(zid,"-");
      }else{
        zid=append(zid,"I");
      }
      cid=append(cid,s[i]);
    }
    del_str(s);
    delete[] tmpstr;
  }else  io.msg(IOL_ERROR|IOL_FATAL,"molcfg::molcfg: molecule %s: no components specified!\n",name);
//
  for(nk=0;ks[nk];++nk);
  K=new eqcfg* [nk];
  if(nk>1) io.msg(IOL_WARN,"molcfg::molcfg: %s[%s]: multiple (%d) constant descriptors found!\n",name,id,nk);
  for(int i=0;i<nk;++i) K[i]=new eqcfg(ks[i],name,id,io);
  del_str(ks);
//
  if(char *s=arg_test(mdata)) io.msg(IOL_WARN,"molcfg: %s[%s]: the following lines were not processed:%s\n",name,id,s);
}

molcfg::~molcfg(void)
{
  if(name) delete[] name;
  if(id) delete[] id;
  if(cid) del_str(cid);
  if(zid) del_str(zid);
  if(K){
    for(int i=0;i<nk;++i) delete K[i];
    delete[] K;
  }
}
