#include <math.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "pack.h"
#include "uts.h"
#include "const.h"

#include "atomcfg.h"

#include "bfcs.h"

// *******************************************
// * bound-free cross-sections               *
// *******************************************

class bfcs *bf_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  int08_t type;
  ::unpack(buf+offs,type);
  switch(type){
    case(BFCS_TYPE_HYD): return new hybf(buf,offs,do_swap,io);
    case(BFCS_TYPE_TAB): return new tabbf(buf,offs,do_swap,io);
    case(BFCS_TYPE_PHY): return new phbf(buf,offs,do_swap,io);   
  }
  return new bfcs(buf,offs,do_swap,io);
}

class bfcs *bf_new(fp_t ee,fp_t ip,fp_t n,fp_t Z,struct bfcfg *cfg,io_class &io_in)
{
  if(!cfg){
    return new bfcs(BFCS_TYPE_NONE);
  } 
  if(!strcmp(cfg->bftype,"HYD")){ 
    //printf("HYD \n");
    return new hybf(ee,ip,n,Z,cfg->l,cfg->v,cfg->n,io_in);
  }
  if(!strcmp(cfg->bftype,"TAB")){ 
    return new tabbf(ee,ip,n,Z,cfg->l,cfg->v,cfg->n,io_in);
    //printf("TAB! \n");
  }
  if(!strcmp(cfg->bftype,"PHY")){ 
    return new phbf(ee,ip,n,Z,cfg->l,cfg->v,cfg->n,io_in);
    //printf("PHY! \n");
  }
  return 0;
}

bfcs::bfcs(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

bfcs::bfcs(int08_t type_in)
{
  type=type_in;
}

int32_t bfcs::size(io_class &io)
{
  return sizeof(int08_t);
}

int32_t bfcs::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::pack(buf,type);
  return offs;
}

int32_t bfcs::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::unpack(buf,type);
  return offs;
}

// *************************************************************************
// * Hydrogenic: SI/CGS invariant version adapted from Mihalas, 2ed. p99,  *
// * with charge dependence and Rydberg constants retained.                *
// *************************************************************************

hybf::hybf(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):bfcs(buf,offs,do_swap,io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

hybf::hybf(fp_t ee,fp_t ip,fp_t np,fp_t Z,fp_t *l_in,fp_t *g_in,int n_in,io_class &io):bfcs(BFCS_TYPE_HYD)
{
  n=n_in;
  l=new fp_t [n];
  g=new fp_t [n];
  memcpy(l,l_in,n*sizeof(fp_t));
  memcpy(g,g_in,n*sizeof(fp_t));
// prepare some numbers
  fp_t fct=sqr(2.0*e*Ryd)/(sqrt(27.0)*pi*e0*me*c*c);
//  fp_t np=max(j[z][l],1);            // primary quantum number
//
  l_b=(h*c)/(ip-ee);        // ionization edge wavelength
  a_b=fct*Z*Z*Z*Z/(sqr(sqr(np))*np);
  //printf("%e %e %e \n", l_b, a_b, a_b*l_b*l_b*l_b);
  //for(int i=0;i<n;++i) fprintf(stderr,"hybf: %E %E %E\n",l_b,l[i],g[i]);
}

hybf::~hybf(void)
{
  if(l) delete[] l;
  if(g) delete[] g;
}

int32_t hybf::size(io_class &io)
{
  int sz=bfcs::size(io);
  sz+=sizeof(int32_t);  // n
  sz+=2*n*sizeof(fp_t); // l,g
  sz+=2*sizeof(fp_t);   // l_b,a_b
  return sz;
}

int32_t hybf::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=bfcs::pack(buf,do_swap,io);
// local stuff
  offs+=::pack(buf+offs,n,do_swap);
  offs+=::pack(buf+offs,l,0,n-1,do_swap);
  offs+=::pack(buf+offs,g,0,n-1,do_swap);
  offs+=::pack(buf+offs,l_b,do_swap);
  offs+=::pack(buf+offs,a_b,do_swap);
  return offs;
}

int32_t hybf::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only unpack local stuff
  int32_t offs=::unpack(buf,n,do_swap);
  offs+=::unpack(buf+offs,l=new fp_t [n],0,n-1,do_swap);
  offs+=::unpack(buf+offs,g=new fp_t [n],0,n-1,do_swap);
  offs+=::unpack(buf+offs,l_b,do_swap);
  offs+=::unpack(buf+offs,a_b,do_swap);
//
  return offs;
}

fp_t interpol(fp_t lam,fp_t *l,fp_t *g,int n)
{ // simple linear interpolation
  for(int i=1;i<n;++i) if(l[i]>lam) return (g[i]*(lam-l[i-1])+g[i-1]*(l[i]-lam))/(l[i]-l[i-1]);
  return g[n-1];
}

fp_t *hybf::U(fp_t *lam,int32_t nlam)
{
  fp_t *s=new fp_t [nlam]; // gII is in tabulated form taken from MULTI config files
  for(int i=0;i<nlam;++i)
    if(lam[i]<=l_b){
      fp_t gII=interpol(lam[i],l,g,n);
      s[i]=gII*a_b*lam[i]*lam[i]*lam[i];
    }else s[i]=0.0;
  return s;
}

fp_t hybf::U(fp_t lam)
{
  if(lam<=l_b){
    fp_t gII=interpol(lam,l,g,n);
    return gII*a_b*lam*lam*lam;
  }
  return 0.0;
}

fp_t *hybf::getlambda(int &np)
{
  if(l_b){
    fp_t *ll=new fp_t [np=n+1];
    memcpy(ll,l,n*sizeof(fp_t));
    ll[n]=l_b;
    return ll;
  }
  return 0;
}

fp_t hybf::getlambda_crit()
{
  return l_b;
}

/*
fp_t hybf::dU(fp_t,fp_t,io_class&)
{
  return 0.0;
}
*/
// ----------------------------------------------------------------------------------------
// SAme here but for tabulated cross-sections, let's see if we can keep up! 
// In this formalism, basically the crosssection is simply interpolated and that's it
// ----------------------------------------------------------------------------------------

// This constructor is simply unpacking 
tabbf::tabbf(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):bfcs(buf,offs,do_swap,io) 
{
  offs+=unpack(buf+offs,do_swap,io);
}

tabbf::tabbf(fp_t ee,fp_t ip,fp_t np,fp_t Z,fp_t *l_in,fp_t *g_in,int n_in,io_class &io):bfcs(BFCS_TYPE_TAB)
{
  n = n_in;
  lambda_tabulated=new fp_t [n];
  bf_tabulated=new fp_t [n];
  memcpy(lambda_tabulated,l_in,n*sizeof(fp_t));
  memcpy(bf_tabulated,g_in,n*sizeof(fp_t));

  fp_t lambda_crit = h * c / (ip - ee);
  lambda_tabulated[n-1] = lambda_crit;

}

tabbf::~tabbf(void)
{
  if(lambda_tabulated) delete[] lambda_tabulated;
  if(bf_tabulated) delete[] bf_tabulated;
}

int32_t tabbf::size(io_class &io)
{
  int sz=bfcs::size(io);
  sz+=sizeof(int32_t);  // n
  sz+=2*n*sizeof(fp_t); // l,g
  return sz;
}

int32_t tabbf::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=bfcs::pack(buf,do_swap,io);
// local stuff
  offs+=::pack(buf+offs,n,do_swap);
  offs+=::pack(buf+offs,lambda_tabulated,0,n-1,do_swap);
  offs+=::pack(buf+offs,bf_tabulated,0,n-1,do_swap);
  return offs;
}

int32_t tabbf::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only unpack local stuff
  int32_t offs=::unpack(buf,n,do_swap);
  offs+=::unpack(buf+offs,lambda_tabulated=new fp_t [n],0,n-1,do_swap);
  offs+=::unpack(buf+offs,bf_tabulated=new fp_t [n],0,n-1,do_swap);

//
  return offs;
}

fp_t tabbf::U(fp_t lam)
{
 
  fp_t cross_section = 0.0;
  if(lam<=lambda_tabulated[n-1] && lam >= lambda_tabulated[0])
    cross_section=interpol(lam,lambda_tabulated,bf_tabulated,n);  

  return cross_section;
}

fp_t *tabbf::getlambda(int &np)
{
  if(n){
    np = n;
    fp_t *ll=new fp_t [n];
    memcpy(ll,lambda_tabulated,n*sizeof(fp_t));
    return ll;
  }
  return 0;
}

fp_t tabbf::getlambda_crit(){
  return lambda_tabulated[n-1];
}

// ----------------------------------------------------------------------------------------
// SAme here but for pseudo-hydrogenic cross-sections, let's see if we can keep up! 
// In this formalism, ionization cross-section is given and then everything is scaled with lambda ** 3.0
// ----------------------------------------------------------------------------------------

// This constructor is simply unpacking 
phbf::phbf(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):bfcs(buf,offs,do_swap,io) 
{
  offs+=unpack(buf+offs,do_swap,io);
}

phbf::phbf(fp_t ee,fp_t ip,fp_t np,fp_t Z,fp_t *l_in,fp_t *g_in,int n_in,io_class &io):bfcs(BFCS_TYPE_PHY)
{
  lambda_crit = h * c / (ip - ee);  
  bf_crit=g_in[0];
}

phbf::~phbf(void)
{
}

int32_t phbf::size(io_class &io)
{
  int sz=bfcs::size(io);
  //sz+=sizeof(int32_t);  // n
  sz+=2*sizeof(fp_t); // sigma, lambda
  return sz;
}

int32_t phbf::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=bfcs::pack(buf,do_swap,io);
// local stuff
  offs+=::pack(buf+offs,lambda_crit,do_swap);
  offs+=::pack(buf+offs,bf_crit,do_swap);
  //offs+=::pack(buf+offs,lambda_tabulated,0,n-1,do_swap);
  //offs+=::pack(buf+offs,bf_tabulated,0,n-1,do_swap);
  return offs;
}

int32_t phbf::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only unpack local stuff
  
  int32_t offs=::unpack(buf,lambda_crit,do_swap);
  offs+=::unpack(buf+offs,bf_crit,do_swap);
  
//
  return offs;
}

fp_t phbf::U(fp_t lam)
{
 
  fp_t cross_section = 0.0;
  if(lam<=lambda_crit)
    cross_section = bf_crit * pow(lam/lambda_crit, 3.0);
  
  //printf("lam = %e lambda_tabulated[n-1] = %e cs = %e \n", lam, lambda_crit, cross_section);
  return cross_section;
}

fp_t *phbf::getlambda(int &np)
{
  np = 10;
  fp_t lambda_min = lambda_crit / 20.0;
  fp_t * l = new fp_t [np];
  fp_t lambda_step = exp(log(20.0) / 9.0);
  for (int i=0;i<10;++i){
    l[i] = lambda_crit / pow(lambda_step,i);
    //printf("%d %e \n", i, l[i]);
  }
  return l;
}

fp_t phbf::getlambda_crit(){
  return lambda_crit;
}
