#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "uts.h"
#include "cfguts.h"
#include "struts.h"

#include "atmol/atom/atomcfg.h"
#include "atmol/mol/molcfg.h"

#include "acfg.h"

acfg::acfg(char *adata,io_class &io)
{
  memset(this,0,sizeof(acfg));
//
  if(!(id=get_arg(adata,"ID",0))) io.msg(IOL_ERROR|IOL_FATAL,"atmos config: atmos has no ID\n");
  if(!(geom=get_arg(adata,"GEOMETRY",0))){
    io.msg(IOL_WARN,"acfg::acfg: no geometry specified, assuming Cartesian 3D\n");
    geom=strcpy(new char [strlen("CARTESIAN3D")+1],"CARTESIAN3D");
  }
  if(!(rts=get_arg(adata,"RTSOLVER",0))){
    io.msg(IOL_WARN,"acfg::acfg: no radiative solver specified, assuming quadratic SC\n");
    rts=strcpy(new char [strlen("QuadSC")+1],"QuadSC");
  }
  filetype=get_arg(adata,"TYPE",0);
  filename=get_arg(adata,"FILE",0);
// cube dimensions
  if(char *tmp_str=get_arg(adata,"NX",0)){
    get_number(tmp_str,nx);
    delete[] tmp_str;
  }else nx=0;
  if(char *tmp_str=get_arg(adata,"NY",0)){
    get_number(tmp_str,ny);
    delete[] tmp_str;
  }else ny=0;
  if(char *tmp_str=get_arg(adata,"NZ",0)){
    get_number(tmp_str,nz);
    delete[] tmp_str;
  }else nz=0;
// grid parameters
  if(char *tmp_str=get_arg(adata,"DX",0)){
    get_number(tmp_str,dx);
    delete[] tmp_str;
  }else dx=0.0;
  if(char *tmp_str=get_arg(adata,"DY",0)){
    get_number(tmp_str,dy);
    delete[] tmp_str;
  }else dy=0.0;
  if(char *tmp_str=get_arg(adata,"DZ",0)){
    get_number(tmp_str,dz);
    delete[] tmp_str;
  }else dz=0.0;
// atoms
  char **atms=0;
  while(char *tmpstr=get_arg(adata,"ATOM",0)){
    FILE *f=fopen(tmpstr,"r");
    char *fc=file2buf(f);
    fclose(f);
//
    while(char *fname=get_arg(fc,"ATOM",0)){
      io.msg(IOL_INFO,"atmos config: ID=\"%s\": atomic data file \"%s\": including \"%s\"\n",id,tmpstr,fname);
      FILE *ff=fopen(fname,"r");
      char *ffc=file2buf(ff);
      fclose(ff);
      delete[] fname;
      char *tmp=new char [strlen(fc)+strlen(ffc)+2];
      sprintf(tmp,"%s\n%s",fc,ffc);
      delete[] fc;
      delete[] ffc;
      fc=tmp;
    }
//
    char **tmps=scope_sep(fc,"atom",io);
    for(int a=0;tmps[a];++a) atms=append(atms,tmps[a]);
    del_str(tmps);
//
    if(char *s=arg_test(fc)) io.msg(IOL_WARN,"atmos config: ID=\"%s\": included file %s (and below): the following lines were not processed: %s\n",id,tmpstr,s);
    
    delete[] tmpstr;
    delete[] fc;
  }
  if(atms){
    for(natm=0;atms[natm];++natm);
    atm=new atmcfg* [natm];
    for(int a=0;a<natm;++a) atm[a]=new atmcfg(atms[a],io);
    del_str(atms);
  }
// scale abundance to total=1
  fp_t sum=0.0;
  for(int a=0;a<natm;++a) sum+=atm[a]->abund;
  for(int a=0;a<natm;++a) atm[a]->abund/=sum; 
// molecules
  char **mols=0;
  while(char *tmpstr=get_arg(adata,"MOL",0)){
    FILE *f=fopen(tmpstr,"r");
    char *fc=file2buf(f);
    fclose(f);
//
    char **tmps=scope_sep(fc,"mol",io);
    for(int m=0;tmps[m];++m) mols=append(mols,tmps[m]);
    del_str(tmps);
//
    delete[] tmpstr;
    delete[] fc; 
  }
  if(mols){
    for(nmol=0;mols[nmol];++nmol);
    mol=new molcfg* [nmol];
    for(int m=0;m<nmol;++m) mol[m]=new molcfg(mols[m],io);
    del_str(mols);
  }
//
  // opacity fudge related stuff:
  if(char *tmp_str=get_arg(adata,"OF",0)){
    get_number(tmp_str,of);
    delete[] tmp_str;
  }else of=0;
  if (of)
    of_filename=get_arg(adata,"OF_FILE",0);
  else 
    of_filename=0;
//
  if(char *s=arg_test(adata)) io.msg(IOL_WARN,"atmos config: ID=%s: the following lines were not processed:%s\n",id,s);


}

acfg::~acfg(void)
{
  if(id) delete[] id;
  if(rts) delete[] rts;
  if(geom) delete[] geom;
  if(filetype) delete[] filetype;
  if(filename) delete[] filename;
  if(of_filename) delete[] of_filename;
  if(atm){
    for(int a=0;a<natm;++a) delete atm[a];
    delete[] atm;
  }
  if(mol){
    for(int m=0;m<nmol;++m) delete mol[m];
    delete[] mol;
  }
}

