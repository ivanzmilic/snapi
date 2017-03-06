
#include "types.h"
#include "io.h"
#include "uts.h"
#include "cfguts.h"
#include "struts.h"
#include "const.h"

#include "atomcfg.h"

tpfcfg::tpfcfg(char *pdata,io_class &io)
{
  if(char *tmpstr=get_arg(pdata,"XION",0)){
    if(get_number(tmpstr,ee)<0) io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: failed to convert XION argument \"%s\" to floating point value\n",tmpstr);
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: no ionization energy specified\n");
//
  if(char *tmpstr=get_arg(pdata,"G2",0)){
    if(get_number(tmpstr,g2)<0) io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: failed to convert G2 argument \"%s\" to floating point value\n",tmpstr);
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: no ionization energy specified\n");
//
  if(char *tmpstr=get_arg(pdata,"XL",0)){
    if(get_number(tmpstr,ll)<0) io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: failed to convert XL argument \"%s\" to floating point value\n",tmpstr);
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: no ionization energy specified\n");
//
  if(char *tmpstr=get_arg(pdata,"GAMMA",0)){
    int nvals;
    fp_t *vals;
    if(get_numbers(tmpstr,vals,nvals)<0) io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: failed to convert GAMMA argument %s to valid floating point values\n",tmpstr);
    g=vals+1;
    n=nvals;
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: no GAMMA values specified for Traving style partition function\n");
//
  if(char *tmpstr=get_arg(pdata,"ALPHA",0)){
    int nvals;
    fp_t *vals;
    if(get_numbers(tmpstr,vals,nvals)<0) io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg failed to convert ALPHA argument %s to valid floating point values\n",tmpstr);
    if(nvals!=n){
      io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: %d alpha values specified for %d gamma values\n",nvals,n);
    }else a=vals+1;
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: no ALPHA values specified for Traving style partition function\n");
//
  if(char *tmpstr=get_arg(pdata,"ASYMP",0)){
    asym=0;
    if(strstr(tmpstr,"Baschek")){
      asym=1;
    }else if(strstr(tmpstr,"Fischel")){
      asym=2;
    }else io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: invalid ASYMP argument \"%s\"\n",tmpstr);
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"tpfcfg::tpfcfg: no ionization energy specified\n");
//
  if(char *s=arg_test(pdata)) io.msg(IOL_WARN,"tpfcfg::tpfcfg: the following lines were not processed:%s\n",s);
}

tpfcfg::~tpfcfg(void)
{
  if(a) delete[] a;
  if(g) delete[] g;
}

pcfg::pcfg(char *pdata,io_class &io)
{
  if(!(pftype=get_arg(pdata,"TYPE",0))) io.msg(IOL_ERROR|IOL_FATAL,"pcfg::pcfg: no partition type specified!\n");
//
  if(!strcmp(pftype,"TRAV")){
    if(char *val_str=get_arg(pdata,"G0",0)){
      if(get_number(val_str,g0)<0) io.msg(IOL_ERROR|IOL_FATAL,"pcfg::pcfg: failed to convert G0 \"%s\" to floating point value\n",val_str);
      delete[] val_str;
    }else io.msg(IOL_ERROR|IOL_FATAL,"pcfg::pcfg: no value specified for constant partition function!\n");
//
    char **pfcs=scope_sep(pdata,"trav",io);
    while(strlen(pdata)&&(pdata[strlen(pdata)-1]=='\n')) pdata[strlen(pdata)-1]=0;
    if(pfcs){
      for(npfc=0;pfcs[npfc];++npfc);
      pfc=new tpfcfg* [npfc];
      for(int i=0;i<npfc;++i) pfc[i]=new tpfcfg(pfcs[i],io);
      del_str(pfcs);
    }else io.msg(IOL_ERROR|IOL_FATAL,"pcfg::pcfg: no parameters specified for TRAV type partition function!\n");
  }else{
    pfc=0;
    npfc=0;
  }
//
  if(!strcmp(pftype,"CONST")){
    if(char *val_str=get_arg(pdata,"VALUE",0)){
      if(get_number(val_str,value)<0) io.msg(IOL_ERROR|IOL_FATAL,"pcfg::pcfg: failed to convert VALUE \"%s\" to floating point value\n",val_str);
      delete[] val_str;
    }else io.msg(IOL_ERROR|IOL_FATAL,"pcfg::pcfg: no value specified for constant partition function!\n");
    
  }else value=0.0;
//
  if(char *s=arg_test(pdata)) io.msg(IOL_WARN,"pcfg::pcfg: the following lines were not processed:%s\n",s);
}

pcfg::~pcfg(void)
{
  if(pfc){
    for(int i=0;i<npfc;++i) delete pfc[i];
    delete[] pfc;
  }
  if(pftype) delete[] pftype;
}

bfcfg::bfcfg(char *bfdata,io_class &io)
{
  if(!(bftype=get_arg(bfdata,"TYPE",0))) io.msg(IOL_ERROR|IOL_FATAL,"bfcfg::bfcfg: no crossection type specified!\n");
//
  if((!strcmp(bftype,"HYDROGENIC"))||(!strcmp(bftype,"HYD"))||(!strcmp(bftype,"TABLE"))||(!strcmp(bftype,"TAB"))||(!strcmp(bftype,"PHY"))){
    if(char *lam_str=get_arg(bfdata,"LAMBDA",0)){
      if(get_numbers(lam_str,l,n)<0) io.msg(IOL_ERROR|IOL_FATAL,"bfcfg::bfcfg: failed to convert LAMBDA argument \"%s\" to floating point values\n",lam_str);
      l+=1;
      delete[] lam_str;
      for(int i=0;i<n;++i) l[i]/=1E8; // from AA to cm...
    }else io.msg(IOL_ERROR|IOL_FATAL,"bfcfg::bfcfg: no wavelength specified for crossection!\n");
//
    if(char *val_str=get_arg(bfdata,"VALUE",0)){
      int m;
      if(get_numbers(val_str,v,m)<0) io.msg(IOL_ERROR|IOL_FATAL,"bfcfg::bfcfg: failed to convert VALUE argument \"%s\" to floating point values\n",val_str);
      if(m!=n) io.msg(IOL_ERROR|IOL_FATAL,"bfcfg::bfcfg: number of wavelengths (%d) not equal to number of values (%d)!\n",n,m);
      v+=1;
      delete[] val_str;
    }else io.msg(IOL_ERROR|IOL_FATAL,"bfcfg::bfcfg: no value specified for crossection!\n");

  }else{
    io.msg(IOL_WARN,"bfcfg::bfcfg: unknown crossection type \"%s\" specified!\n",bftype);
    l=v=0;
    n=0;
  }
//
  if(char *s=arg_test(bfdata)) io.msg(IOL_WARN,"bfcfg::bfcfg: the following lines were not processed:%s\n",s);
}

bfcfg::~bfcfg(void)
{
  if(l) delete[] l;
  if(v) delete[] v;
  delete[] bftype;
}

crcfg::crcfg(char *crdata,io_class &io)
{
  if(!(crtype=get_arg(crdata,"TYPE",0))) io.msg(IOL_ERROR|IOL_FATAL,"crcfg::crcfg: no collisional rate type specified!\n");
//
  if((!strcmp(crtype,"HYDROGENIC"))||(!strcmp(crtype,"HYD"))||(!strcmp(crtype,"TABLE"))||(!strcmp(crtype,"TAB"))){
    if(char *tmp_str=get_arg(crdata,"TEMP",0)){
      if(get_numbers(tmp_str,t,n)<0) io.msg(IOL_ERROR|IOL_FATAL,"crcfg::crcfg: failed to convert TEMP argument \"%s\" to floating point values\n",tmp_str);
      t+=1;
      delete[] tmp_str;
    }else io.msg(IOL_ERROR|IOL_FATAL,"crcfg::crcfg: no temperatures specified for crossection!\n");
//
    if(char *val_str=get_arg(crdata,"VALUE",0)){
      int m;
      if(get_numbers(val_str,v,m)<0) io.msg(IOL_ERROR|IOL_FATAL,"crcfg::crcfg: failed to convert VALUE argument \"%s\" to floating point values\n",val_str);
      if(m!=n) io.msg(IOL_ERROR|IOL_FATAL,"crcfg::crcfg: number of wavelengths (%d) not equal to number of values (%d)!\n",n,m);
      v+=1;
      delete[] val_str;
    }else io.msg(IOL_ERROR|IOL_FATAL,"crcfg::crcfg: no value specified for crossection!\n");

  }else{
    io.msg(IOL_WARN,"crcfg::crcfg: unknown crossection type \"%s\" specified!\n",crtype);
    t=v=0;
    n=0;
  }
//
  if(char *s=arg_test(crdata)) io.msg(IOL_WARN,"crcfg::crcfg: the following lines were not processed:%s\n",s);
}

crcfg::~crcfg(void)
{
  if(t) delete[] t;
  if(v) delete[] v;
  delete[] crtype;
}

lcfg::lcfg(char *ldata,io_class &io)
{
  char **bfcss=scope_sep(ldata,"bf",io);
  char **crcss=scope_sep(ldata,"cr",io);
  while(strlen(ldata)&&(ldata[strlen(ldata)-1]=='\n')) ldata[strlen(ldata)-1]=0;
//
  if(char *tmpstr=get_arg(ldata,"EE",0)){
    if(sscanf(tmpstr,"%lf",&ee)!=1) io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: failed to convert energy argument \"%s\" to a valid floating point value\n",tmpstr);
    delete[] tmpstr;
    ee*=ev2eh*h; // assume eV: convert to energy
  }else io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: no energy specified for level \n");
//
  if(char *tmpstr=get_arg(ldata,"G",0)){
    int gi;
    if(sscanf(tmpstr,"%d",&gi)!=1) io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: failed to convert weight argument \"%s\" to a valid integer value\n",tmpstr);
    delete[] tmpstr;
    g=gi;
  }else io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: no statistical weight specified for level \n");
//
  if(char *tmpstr=get_arg(ldata,"J",0)){
    int ji;
    if(sscanf(tmpstr,"%d",&ji)!=1) io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: failed to convert angular quantum number argument \"%s\" to a valid integer value\n",tmpstr);
    delete[] tmpstr;
    j_qn=ji;
  }
  else{
    io.msg(IOL_WARN,"lcfg::lcfg: no angular quantum number specified for level, setting to zero \n");
    j_qn = 0;
  }
//
  if(char *tmpstr=get_arg(ldata,"L",0)){
    int li;
    if(sscanf(tmpstr,"%d",&li)!=1) io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: failed to convert orbital quantum number argument \"%s\" to a valid integer value\n",tmpstr);
    delete[] tmpstr;
    l_qn=li;
  }else{ 
    io.msg(IOL_WARN,"lcfg::lcfg: no orbital quantum number specified for level, setting to zero \n");
    l_qn = 0;
  }
//
    if(char *tmpstr=get_arg(ldata,"GL",0)){
    fp_t li;
    if(sscanf(tmpstr,"%lf",&li)!=1) io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: failed to convert lande factor argument \"%s\" to a valid integer value\n",tmpstr);
    delete[] tmpstr;
    g_lande=li;
  }else{ 
    io.msg(IOL_WARN,"lcfg::lcfg: no lande factor specified for level, setting to zero \n");
    g_lande = 0.0;
  }
//

  if(char *tmpstr=get_arg(ldata,"N",0)){
    if(sscanf(tmpstr,"%d",&n)!=1) io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: failed to convert principal quantum number argument \"%s\" to a valid integer\n",tmpstr);
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: no principal quantum number specified for level \n");
//
  if(char *trans_str=get_arg(ldata,"A",0)){
    if(get_numbers(trans_str,A,nt)<0) io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: failed to convert transition probability argument \"%s\" to valid numbers\n",trans_str);
    A+=1;
    delete[] trans_str;
  }else A=0;
//  else io.msg(IOL_ERROR|IOL_FATAL,"lcfg::lcfg: no principal quantum number specified for level \n");
//  
  bf=0;
  if(bfcss){
    if(bfcss[0]){
      if(bfcss[1]) io.msg(IOL_ERROR,"lcfg::lcfg: multiple bound-free crossection descriptors found for level energy %E: using the first one\n",ee);
      bf=new bfcfg(bfcss[0],io);
    }
    del_str(bfcss);
  }
//
  ncr=0;
  cr=0;
  if(crcss){
    for(ncr=0;crcss[ncr];++ncr);
    cr=new crcfg* [ncr];
    for(int l=0;l<ncr;++l) cr[l]=new crcfg(crcss[l],io);
    del_str(crcss);
  }
//
  if(char *s=arg_test(ldata)) io.msg(IOL_WARN,"lcfg::lcfg: the following lines were not processed:%s\n",s);
}

lcfg::~lcfg(void)
{
  if(bf) delete bf;
  if(A) delete[] A;
}

icfg::icfg(int z,char *idata,io_class &io)
{
  char **levels=scope_sep(idata,"level",io);
  while(strlen(idata)&&(idata[strlen(idata)-1]=='\n')) idata[strlen(idata)-1]=0;
  char **partitions=scope_sep(idata,"partition",io);
  while(strlen(idata)&&(idata[strlen(idata)-1]=='\n')) idata[strlen(idata)-1]=0;
//
  if(char *tmpstr=get_arg(idata,"CHARGE",0)){
    int tmpi;
    if(sscanf(tmpstr,"%d",&tmpi)!=1) io.msg(IOL_ERROR|IOL_FATAL,"icfg::icfg: failed to convert ion CHARGE argument \"%s\" to a integer\n",tmpstr);
    charge=tmpi;
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"icfg::icfg: no charge specified for ion\n");
//
  if(char *tmpstr=get_arg(idata,"EION",0)){
    if(sscanf(tmpstr,"%lf",&eion)!=1) io.msg(IOL_ERROR|IOL_FATAL,"icfg::icfg: failed to convert ion EION argument \"%s\" to a valid floating point value\n",tmpstr);
    if(charge==z){ 
      if(isfinite(eion))
        io.msg(IOL_ERROR|IOL_FATAL,"icfg::icfg: charge(%d)==atomic number(%d): ion cannot be ionized but has finite ionization potential (%E)!\n",charge,z,eion);
     }else eion*=ev2eh*h; // assume eV
     delete[] tmpstr;
  }else if(charge==z){
    eion=NAN;
  }else io.msg(IOL_ERROR|IOL_FATAL,"icfg::icfg: no ionization energy specified for ion\n");
//  
  pf=0;
  if(partitions){
    if(partitions[0]){
      if(partitions[1]) io.msg(IOL_ERROR,"icfg::icfg: multiple partition function descriptors found for ion charge %d: using the first one\n",charge);
      pf=new pcfg(partitions[0],io);
    }
    del_str(partitions);
  }
//  if(!pf) io.msg(IOL_ERROR|IOL_FATAL,"icfg::icfg: no partition function descriptor found for ion charge %d\n",charge);
  nl=0;
  level=0;
  if(levels){
    for(nl=0;levels[nl];++nl);
    level=new lcfg* [nl];
    for(int l=0;l<nl;++l) level[l]=new lcfg(levels[l],io);
    del_str(levels);
  }
//
  if(char *s=arg_test(idata)) io.msg(IOL_WARN,"icfg::icfg: the following lines were not processed:%s\n",s);
}

icfg::~icfg(void)
{
  if(pf) delete pf;
  for(int l=0;l<nl;++l) delete level[l];
  delete[] level;
}

atmcfg::atmcfg(char *adata,io_class &io)
{
  char **ions=scope_sep(adata,"ion",io);
  while(strlen(adata)&&(adata[strlen(adata)-1]=='\n')) adata[strlen(adata)-1]=0;
//
  if(!(name=get_arg(adata,"NAME",0))) io.msg(IOL_ERROR|IOL_FATAL,"no name specified for atom!\n");
  if(!(id=get_arg(adata,"ID",0))) io.msg(IOL_ERROR|IOL_FATAL,"no id specified for atom %s\n",name);
  if(char *tmpstr=get_arg(adata,"Z",0)){
    int zz;
    if(sscanf(tmpstr,"%d",&zz)!=1) io.msg(IOL_ERROR|IOL_FATAL,"failed to convert atomic number argument %s to integer for atom %s[%s]\n",tmpstr,name,id);
    z=(int08_t)zz;
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"no atomic number specified for atom %s[%s]\n",name,id);
//
  if(char *tmpstr=get_arg(adata,"MASS",0)){
    if(sscanf(tmpstr,"%lf",&mass)!=1) io.msg(IOL_ERROR|IOL_FATAL,"failed to convert atomic mass argument %s to a valid floating point value for atom %s[%s]\n",tmpstr,name,id);
    delete[] tmpstr;
  }else io.msg(IOL_ERROR|IOL_FATAL,"no mass specified for atom %s[%s]\n",name,id);
//
  if(char *tmpstr=get_arg(adata,"ABUND",0)){
    if(sscanf(tmpstr,"%lf",&abund)!=1) io.msg(IOL_ERROR|IOL_FATAL,"failed to convert abundance argument %s to a valid floating point value for atom %s[%s]\n",tmpstr,name,id);
    delete[] tmpstr;
// number fraction relative to a number of choice (will be renormalized later)
  }else io.msg(IOL_ERROR|IOL_FATAL,"no abundance specified for atom %s[%s]\n",name,id);

  if(char *tmpstr=get_arg(adata,"NLTE",0)){
    if(sscanf(tmpstr,"%d",&nlte)!=1) io.msg(IOL_ERROR|IOL_FATAL,"failed to convert nlte argument %s to a valid int value for atom %s[%s]\n",tmpstr,name,id);
    delete[] tmpstr;
  }else nlte = 0;
//
//
  ion=new icfg* [z+1];
  memset(ion,0,(z+1)*sizeof(struct icfg*));
  if(ions){
    for(ni=0;ions[ni];++ni);
    for(int i=0;i<ni;++i){
      icfg *tmp=new icfg(z,ions[i],io);
      int08_t zz=tmp->charge;
      if(ion[zz]) io.msg(IOL_ERROR|IOL_FATAL,"atmcfg::atmcfg: ion with charge %d already specified for atom %s[%s]\n",zz,name,id);
      ion[zz]=tmp;
    }
    del_str(ions);
  }

  if(char *s=arg_test(adata)) io.msg(IOL_WARN,"atmcfg::atmcfg: %s[%s]: the following lines were not processed:%s\n",name,id,s);
}

atmcfg::~atmcfg(void)
{
  if(ion){
    for(int i=0;i<=z;++i) if(ion[i]) delete ion[i];
    delete[] ion;
  }
  if(name) delete[] name;
  if(id) delete[] id;
}

