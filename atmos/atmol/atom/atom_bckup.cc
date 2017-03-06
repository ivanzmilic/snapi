#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "uts.h"
#include "io.h"
#include "pack.h"
#include "struts.h"
#include "mem.h"
#include "profiles.h"
#include "const.h"
#include "atomcfg.h"
#include "mathtools.h"

#include "H/hydrogen.h"
#include "He/helium.h"
//#include "Li/lithium.h"
//#include "Be/beryllium.h"
//#include "B/boron.h"
//#include "C/carbon.h"
//#include "N/nitrogen.h"
//#include "O/oxygen.h"
//#include "F/fluorine.h"
//#include "Ne/neon.h"
//#include "Na/sodium.h"
//#include "Mg/magnesium.h"
//#include "Al/aluminium.h"
//#include "Si/silicon.h"
//#include "P/phosphorus.h"
//#include "S/sulfur.h"
//#include "Cl/chlorine.h"
//#include "Ar/argon.h"
//#include "K/potassium.h"
//#include "Ca/calcium.h"
//#include "Sc/scandium.h"
//#include "Ti/titanium.h"
//#include "V/vanadium.h"
//#include "Cr/chromium.h"
//#include "Mn/manganese.h"
//#include "Fe/iron.h"
//#include "Co/cobalt.h"
//#include "Ni/nickel.h"
//#include "Cu/copper.h"
//#include "Zn/zinc.h"

#include "atomcfg.h"
#include "partf.h"  // partition functions
#include "bfcs.h"   // bound-free crossections
#include "colr.h"   // collisional rates
#include "atom.h"

atmol *atom_new(atmcfg *cfg,io_class &io_in)
{
  switch(cfg->z){
    case( 1): return new hydrogen(cfg,io_in);
    case( 2): return new helium(cfg,io_in);
  }
  return new atom(cfg,io_in); // something generic?
}

atmol *atom_new(uint64_t numid,uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in)
{
  switch(numid){
    case( 1): return new hydrogen(buf,offs,do_swap,io_in);
    case( 2): return new helium(buf,offs,do_swap,io_in);
//    case( 3): return new lithium(buf,offs,do_swap,io_in);
//    case( 4): return new berilium(buf,offs,do_swap,io_in);
//    case( 5): return new boron(buf,offs,do_swap,io_in);
//    case( 6): return new carbon(buf,offs,do_swap,io_in);
//    case( 7): return new nitrogen(buf,offs,do_swap,io_in);
//    case( 8): return new oxygen(buf,offs,do_swap,io_in);
//    case( 9): return new fluorine(buf,offs,do_swap,io_in);
//    case(10): return new neon(buf,offs,do_swap,io_in);
//    case(11): return new sodium(buf,offs,do_swap,io_in);
//    case(12): return new magnesium(buf,offs,do_swap,io_in);
//    case(13): return new aluminium(buf,offs,do_swap,io_in);
//    case(14): return new silicon(buf,offs,do_swap,io_in);
//    case(15): return new phosphorus(buf,offs,do_swap,io_in);
//    case(16): return new sulfur(buf,offs,do_swap,io_in);
//    case(17): return new chlorine(buf,offs,do_swap,io_in);
//    case(18): return new argon(buf,offs,do_swap,io_in);
//    case(19): return new potassium(buf,offs,do_swap,io_in);
//    case(20): return new calcium(buf,offs,do_swap,io_in);
//    case(21): return new scandium(buf,offs,do_swap,io_in);
//    case(22): return new titanium(buf,offs,do_swap,io_in);
//    case(23): return new vanadium(buf,offs,do_swap,io_in);
//    case(24): return new chromium(buf,offs,do_swap,io_in);
//    case(25): return new manganese(buf,offs,do_swap,io_in);
//    case(26): return new iron(buf,offs,do_swap,io_in);
//    case(27): return new cobalt(buf,offs,do_swap,io_in);
//    case(28): return new nickel(buf,offs,do_swap,io_in);
//    case(29): return new copper(buf,offs,do_swap,io_in);
//    case(30): return new zinc(buf,offs,do_swap,io_in);
  }
  return new atom(buf,offs,do_swap,io_in); // something generic?
}


fp_t ***newtrans(uint08_t Z,uint16_t *nl)
{
  int nn=0,nt=0;
  for(uint08_t z=0;z<=Z;++z){
    nn+=nl[z];
    nt+=(nl[z]*nl[z]);
  }
//
  fp_t ***p=new fp_t** [Z+1];
  p[0]=new fp_t* [nn];
  p[0][0]=new fp_t [nt];
  int m=0,n=0;
  for(uint08_t z=1;z<=Z;++z){
    p[z]=p[0]+(m+=nl[z-1]);
    p[z][0]=p[0][0]+(n+=nl[z-1]*nl[z-1]);
  }
//
  for(uint08_t z=0;z<=Z;++z) for(uint16_t l=1;l<nl[z];++l) p[z][l]=p[z][l-1]+nl[z];
//
  return p;
}

void deltrans(fp_t ***p)
{
  delete[] (p[0][0]);
  delete[] (p[0]);
  delete[] (p);
}

atom::atom(atmcfg *cfg,io_class &io_in):atmol(cfg->name,cfg->id,io_in)
{
//  io_in.msg(IOL_INFO,"atom::atom: %s[%s]\n",cfg->name,cfg->id);
  Z=cfg->z;
  mass=mu*cfg->mass;
  abund=cfg->abund;
  ip=new fp_t [Z];
// setup each ionization stage
  nl=new uint16_t [Z+1];
  ee=new fp_t* [Z+1];
  g=new uint08_t* [Z+1];
  j=new uint08_t* [Z+1];
  flags=new uint08_t* [Z+1];
  // New ones:
  j_qn = new uint08_t* [Z+1]; // Angular quantum number
  l_qn = new uint08_t* [Z+1]; // Orbital quantum number

  partf=new pf* [Z+1];
  bf=new bfcs** [Z+1];
// unpack cfg structure
  nmap=0;
  for(int i=0;i<=Z;++i){
    if(cfg->ion[i]){
      if(i<Z) ip[i]=cfg->ion[i]->eion;
      partf[i]=pf_new(cfg->ion[i]->charge,cfg->ion[i]->pf);
// initialize the levels...
      nl[i]=max(cfg->ion[i]->nl,1); // at least one level...
      bf[i]=new bfcs* [nl[i]];
      ee[i]=new fp_t [nl[i]];
      g[i]=new uint08_t [nl[i]];
      j[i]=new uint08_t [nl[i]];
      flags[i]=new uint08_t [nl[i]];
      j_qn[i]=new uint08_t [nl[i]];
      l_qn[i]=new uint08_t [nl[i]];
      
      io_in.msg(IOL_INFO,"atom::atom: %s: ion=%d nl=%d\n",cfg->name,i,nl[i]);
      for(int l=0;l<nl[i];++l){
        if(l<cfg->ion[i]->nl){
          ee[i][l]=cfg->ion[i]->level[l]->ee;
          g[i][l]=cfg->ion[i]->level[l]->g;
          j[i][l]=cfg->ion[i]->level[l]->n;

          j_qn[i][l]=cfg->ion[i]->level[l]->j_qn;
          l_qn[i][l]=cfg->ion[i]->level[l]->l_qn;
          //printf("Z = %d j = %d l = %d \n",Z, j_qn[i][l], l_qn[i][l]);

          flags[i][l]=0;
          bf[i][l]=bf_new(ee[i][l],ip[i],j[i][l],i+1,cfg->ion[i]->level[l]->bf,io_in);
        }else{ // no levels specified
          ee[i][l]=0.0;
          g[i][l]=1;
          j[i][l]=0;
          j_qn[i][l]=0;
          l_qn[i][l]=0;
          flags[i][l]=FL_NONE; // not a level
          bf[i][l]=new bfcs(BFCS_TYPE_NONE);
        }
        io_in.msg(IOL_INFO,"atom::atom: %s: ion=%d nl=%d ee=%E j=%d g=%d\n",cfg->name,i,nl[i],ee[i][l],j[i][l],g[i][l]);
      }
    }else{
      if(i<Z) ip[i]=NAN;
      partf[i]=new pf(PF_TYPE_NONE);
      nl[i]=1;
      (bf[i]=new bfcs* [nl[i]])[0]=new bfcs(BFCS_TYPE_NONE); // unknown cross-section
      (ee[i]=new fp_t [nl[i]])[0]=0.0;   // ground state only
      (g[i]=new uint08_t [nl[i]])[0]=0;  // unknown...
      (j[i]=new uint08_t [nl[i]])[0]=0;  // unknown...
      (j_qn[i]=new uint08_t [nl[i]])[0]=0;  // unknown...
      (l_qn[i]=new uint08_t [nl[i]])[0]=0;  // unknown..
      (flags[i]=new uint08_t [nl[i]])[0]=FL_NONE;  // not a level...
    }
    nmap+=nl[i];
  }
// level maps...
  lmap=new uint16_t [nmap];
  zmap=new uint08_t [nmap];
  for(int z=0,m=0;z<=Z;++z)
    for(int l=0;l<nl[z];++l,++m){
      lmap[m]=l;
      zmap[m]=z;
    }
//
  A=newtrans(Z,nl);
  B=newtrans(Z,nl);
  col_dam_cross_section = newtrans(Z, nl); // Also the same with cross-section
  alpha_col_dam = newtrans(Z, nl);
  osc_str = newtrans(Z, nl); // And with the oscillator strength
//
  for(int z=0;z<=Z;++z) for(int i=0;i<nl[z];++i) for(int ii=0;ii<nl[z];++ii) A[z][i][ii]=B[z][i][ii]=0.0;
//
  for(int z=0;z<=Z;++z)
    for(int lu=1;lu<nl[z];++lu) // upper level
      for(int ll=0;ll<lu;++ll){ // lower level
        fp_t lam=h*c/(ee[z][lu]-ee[z][ll]);
        A[z][lu][ll]=cfg->ion[z]->level[lu]->A[ll]; // Aul
        B[z][lu][ll]=((lam*lam*lam)/(2.0*h*c))*A[z][lu][ll]; // Bul
        B[z][ll][lu]=((fp_t)(g[z][lu])/(fp_t)(g[z][ll]))*B[z][lu][ll]; // Blu

        // Modifications in order to have B in proper units.
        B[z][lu][ll] *= lam * lam / c;
        B[z][ll][lu] *= lam * lam / c;
        //
        osc_str[z][ll][lu] = osc_str[z][lu][ll] = 0.0; // This will be computed later, in get lambda function

        //col_dam_cross_section = compute_col_dam(z, lu, ll); // Function compute_col_dam not yet written
        col_dam_cross_section[z][lu][ll] = col_dam_cross_section[z][ll][lu] = 0.0;
        alpha_col_dam[z][lu][ll] = alpha_col_dam[z][ll][lu] = 0.0;

        //printf("Z = %d lu = %d ll = %d \n",Z, l_qn[z][lu], l_qn[z][ll]);

        compute_damp_col(z, lu, ll);

        //printf("Z = %d lu = %d ll = %d \n",Z, l_qn[z][lu], l_qn[z][ll]);
      }
//
// setup ionization populations
//  N=new fp_t [Z+1];
//  memset(N,0,(Z+1)*sizeof(fp_t));
//  n=new fp_t* [Z+1];
//  for(int i=0;i<=Z;++i){
//    n[i]=new fp_t [nl[i]];
//    for(int l=0;l<nl[i];++l) n[i][l]=0.0;
//  }
//
  numid=Z;
  io_in.msg(IOL_INFO,"atom::atom: %s = 0x%016lX\n",name,numid);
  pop=0;
}

atom::atom(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):atmol(buf,offs,do_swap,io_in)
{
  offs+=unpack(buf+offs,do_swap,io_in);
  pop=0;
}

atom::~atom(void)
{
//  if(n){
//    for(int z=0;z<=Z;++z) delete[] n[z];
//    delete[] n;
//  }
//  if(N) delete[] N;
//
  if(pop) popclean(x1l,x1h,x2l,x2h,x3l,x3h);
  if(ip) delete[] ip;
  if(partf){
    for(int16_t i=0;i<=Z;++i) if(partf[i]) delete partf[i];
    delete[] partf;
  }
  if(flags){
    for(int z=0;z<=Z;++z) delete[] flags[z];
    delete[] flags;
  }
  if(j){
    for(int z=0;z<=Z;++z) delete[] j[z];
    delete[] j;
  }
  if(g){
    for(int z=0;z<=Z;++z) delete[] g[z];
    delete[] g;
  }
   if(j_qn){
    for(int z=0;z<=Z;++z) delete[] j_qn[z];
    delete[] j_qn;
  }
   if(j){
    for(int z=0;z<=Z;++z) delete[] l_qn[z];
    delete[] l_qn;
  }
  if(A) deltrans(A);
  if(B) deltrans(B);

  if(osc_str) deltrans(osc_str);
  if(col_dam_cross_section) deltrans(col_dam_cross_section);
  if(alpha_col_dam) deltrans(alpha_col_dam);
  
  if(ee){
    for(int z=0;z<=Z;++z) delete[] ee[z];
    delete[] ee;
  }
  if(bf){
    for(int z=0;z<=Z;++z){
      for(int l=0;l<nl[z];++l) delete bf[z][l];
      delete[] bf[z];
    }
    delete[] bf;
  }
  if(zmap) delete[] zmap;
  if(lmap) delete[] lmap;
  if(nl) delete[] nl;
}

int32_t atom::size(io_class &io_in)
{
  int32_t sz=atmol::size(io_in);
  sz+=sizeof(uint08_t); // z
  sz+=(Z+1)*sizeof(fp_t); // abund,ip
  sz+=(Z+1)*sizeof(uint16_t); // nl
  sz+=(nmap+1)*sizeof(uint16_t); // nmap,lmap
  sz+=nmap*sizeof(uint08_t); // zmap
  for(uint08_t i=0;i<=Z;++i){
    sz+=(sizeof(fp_t)+5*sizeof(int08_t))*(int32_t)(nl[i]);  // ee,j,g,flags
    sz+=partf[i]->size(io_in); // partition function data
    for(uint16_t l=0;l<nl[i];++l) sz+=bf[i][l]->size(io_in);  // bound-free crossection data
    for(uint16_t l=0;l<nl[i];++l) sz+=5*nl[i]*sizeof(fp_t); // A,B, osc_str, dam_col_cross_section, alpha
  }
//
//  sz+=(Z+1)*sizeof(fp_t); // N
//  for(uint08_t i=0;i<=Z;++i) sz+=sizeof(fp_t)*(int32_t)(nl[i]);  // n
//
//  io.msg(IOL_INFO,"atom::size %d\n",sz);
//
  return sz;
}

int32_t atom::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
  int32_t offs=atmol::pack(buf,do_swap,io);
  offs+=::pack(buf+offs,Z);
  offs+=::pack(buf+offs,abund,do_swap);
  offs+=::pack(buf+offs,ip,0,Z-1,do_swap);
  offs+=::pack(buf+offs,nl,0,Z,do_swap);
  offs+=::pack(buf+offs,nmap,do_swap);
  offs+=::pack(buf+offs,lmap,0,nmap-1,do_swap);
  offs+=::pack(buf+offs,zmap,0,nmap-1);
  for(uint08_t i=0;i<=Z;++i){
    offs+=::pack(buf+offs,ee[i],0,nl[i]-1,do_swap);
    offs+=::pack(buf+offs,g[i],0,nl[i]-1);
    offs+=::pack(buf+offs,j[i],0,nl[i]-1);
    offs+=::pack(buf+offs,j_qn[i],0,nl[i]-1);
    offs+=::pack(buf+offs,l_qn[i],0,nl[i]-1);
    offs+=::pack(buf+offs,flags[i],0,nl[i]-1);
    offs+=partf[i]->pack(buf+offs,do_swap,io);
    for(uint16_t l=0;l<nl[i];++l) offs+=bf[i][l]->pack(buf+offs,do_swap,io);
    for(uint16_t l=0;l<nl[i];++l){
      offs+=::pack(buf+offs,A[i][l],0,nl[i]-1,do_swap);
      offs+=::pack(buf+offs,B[i][l],0,nl[i]-1,do_swap);
      offs+=::pack(buf+offs,osc_str[i][l],0,nl[i]-1,do_swap);
      offs+=::pack(buf+offs,col_dam_cross_section[i][l],0,nl[i]-1,do_swap);
      offs+=::pack(buf+offs,alpha_col_dam[i][l],0,nl[i]-1,do_swap);

    }
  }
//
//  offs+=::pack(buf+offs,N,0,Z,do_swap);
//  for(uint08_t i=0;i<=Z;++i) offs+=::pack(buf+offs,n[i],0,nl[i]-1,do_swap);
//
  return offs;
}

int32_t atom::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
// only unpack local stuff
  int32_t offs=::unpack(buf,Z);
  offs+=::unpack(buf+offs,abund,do_swap);
  offs+=::unpack(buf+offs,ip=new fp_t [Z],0,Z-1,do_swap);
  offs+=::unpack(buf+offs,nl=new uint16_t [Z+1],0,Z,do_swap);
  offs+=::unpack(buf+offs,nmap,do_swap);
  offs+=::unpack(buf+offs,lmap=new uint16_t [nmap],0,nmap-1,do_swap);
  offs+=::unpack(buf+offs,zmap=new uint08_t [nmap],0,nmap-1);
//
  ee=new fp_t* [Z+1];
  g=new uint08_t* [Z+1];
  j=new uint08_t* [Z+1];
  j_qn=new uint08_t* [Z+1];
  l_qn=new uint08_t* [Z+1];
  flags=new uint08_t* [Z+1];
  partf=new pf* [Z+1];
  bf=new bfcs** [Z+1];
  A=newtrans(Z,nl);
  B=newtrans(Z,nl);
  osc_str = newtrans(Z, nl);
  col_dam_cross_section = newtrans(Z, nl);
  alpha_col_dam = newtrans(Z, nl);
  for(uint08_t i=0;i<=Z;++i){
    offs+=::unpack(buf+offs,ee[i]=new fp_t [nl[i]],0,nl[i]-1,do_swap);
    offs+=::unpack(buf+offs,g[i]=new uint08_t [nl[i]],0,nl[i]-1);
    offs+=::unpack(buf+offs,j[i]=new uint08_t [nl[i]],0,nl[i]-1);
    offs+=::unpack(buf+offs,j_qn[i]=new uint08_t [nl[i]],0,nl[i]-1);
    offs+=::unpack(buf+offs,l_qn[i]=new uint08_t [nl[i]],0,nl[i]-1);
    offs+=::unpack(buf+offs,flags[i]=new uint08_t [nl[i]],0,nl[i]-1);
    partf[i]=pf_new(buf,offs,do_swap,io);
    bf[i]=new bfcs* [nl[i]];
    for(uint16_t l=0;l<nl[i];++l) bf[i][l]=bf_new(buf,offs,do_swap,io);
    for(uint16_t l=0;l<nl[i];++l){
      offs+=::unpack(buf+offs,A[i][l],0,nl[i]-1,do_swap);
      offs+=::unpack(buf+offs,B[i][l],0,nl[i]-1,do_swap);
      offs+=::unpack(buf+offs,osc_str[i][l],0,nl[i]-1,do_swap);
      offs+=::unpack(buf+offs,col_dam_cross_section[i][l],0,nl[i]-1,do_swap);
      offs+=::unpack(buf+offs,alpha_col_dam[i][l],0,nl[i]-1,do_swap);
      
    }
  }
//
//  n=new fp_t* [Z+1];
//  offs+=::unpack(buf+offs,N=new fp_t [Z+1],0,Z,do_swap);
//  for(uint08_t i=0;i<=Z;++i) offs+=::unpack(buf+offs,n[i]=new fp_t [nl[i]],0,nl[i]-1,do_swap);
//
  return offs;
}

fp_t * atom::add(fp_t *src,fp_t *dst,int32_t nn){
  if(src){
    for(int32_t i=0;i<nn;++i) dst[i]+=src[i];
    delete[] src;
  }
  return dst;
}
  
fp_t *** atom::add(fp_t ***src,fp_t ***dst,int32_t ll1,int32_t ul1,int32_t ll2,int32_t ul2,int32_t ll3,int32_t ul3){
  if(src){
    int32_t nn=(ul1-ll1+1)*(ul2-ll2+1)*(ul3-ll3+1);
    for(int32_t i=0;i<nn;++i) dst[ll1][ll2][ll3+i]+=src[ll1][ll2][ll3+i];
    del_ft3dim(src,ll1,ul1,ll2,ul2,ll3,ul3);
  }
  return dst;
}

fp_t **** atom::add(fp_t **** src, fp_t **** dst, int32_t ll1, int32_t ul1, int32_t ll2, int32_t ul2, int32_t ll3, int32_t ul3, int32_t ll4, int32_t ul4){

  if(src){
    int32_t nn = (ul1 - ll1 + 1) * (ul2 - ll2 + 1) * (ul3 - ll3 + 1) * (ul4 - ll4 + 1);
    for (int32_t i=0;i<nn;++i)
      dst[ll1][ll2][ll3][ll4+i] += src[ll1][ll2][ll3][ll4+i];
    del_ft4dim(src, ll1, ul1, ll2, ul2, ll3, ul3, ll4, ul4);
  }
  return dst;
}

fp_t ***** atom::add(fp_t *****src, fp_t *****dst, int32_t ll1, int32_t ul1, int32_t ll2, int32_t ul2, int32_t ll3, int32_t ul3, int32_t ll4, int32_t ul4, int32_t ll5, int32_t ul5){

  if(src){
    int32_t nn = (ul1 - ll1 + 1) * (ul2 - ll2 + 1) * (ul3 - ll3 + 1) * (ul4 - ll4 + 1) * (ul5 - ll5 + 1);
    for (int32_t i =0 ; i<nn; ++i)
      dst[ll1][ll2][ll3][ll4][ll5+i] += src[ll1][ll2][ll3][ll4][ll5+i];
    del_ft5dim(src, ll1, ul1, ll2, ul2, ll3, ul3, ll4, ul4, ll5, ul5);
  }

  return dst;
}

// 
fp_t *atom::rayleigh_em(fp_t *lambda,int32_t nlambda)
// *************************************************************************
// Rayleigh scattering emission=redistributed(J)/absorption?
// *************************************************************************
{
  fp_t *em=new fp_t [nlambda]; // Milic changed from op to em just for the clarity, does not change much
  memset(em,0,nlambda*sizeof(fp_t));
  return em;
}

fp_t *atom::rayleigh_op(fp_t *lambda,int32_t nlambda)
// *************************************************************************
// Rayleigh scattering absorption
// *************************************************************************
{
  fp_t *op=new fp_t [nlambda];
  memset(op,0,nlambda*sizeof(fp_t));
  return op;
}

fp_t *atom::freefree_em(fp_t T,fp_t ne,fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i)
// *************************************************************************
// bremsstrahlung: B[Te/Ti]/alpha?
// *************************************************************************
{
  fp_t *em=new fp_t [nlambda],gff=1.0;
  memset(em,0,nlambda*sizeof(fp_t));
  fp_t a=3.7E8*gff*ne*c*c*c/sqrt(T);
  for(int z=0;z<=Z;++z){
    for(int i=0;i<nlambda;++i) em[i]+=a*z*z*pop[x1i][x2i][x3i].N[z];
  }
  fp_t b=-(h*c)/(k*T);
  for(int i=0;i<nlambda;++i){
    em[i]*=lambda[i]*lambda[i]*lambda[i]*(1.0-exp(b/lambda[i]));
//    em[i]=planck_nu(T,c/lambda[i])/em[i];
  }
  return em;
}

fp_t *atom::freefree_op(fp_t T,fp_t ne,fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i)
// *************************************************************************
// * Free-Free opacities:
// * default: Hydrogenic: SI/CGS invariant?
// * Hydrogenic expression from Mihalas, 2ed. p101, with charge dependence
// * and Rydberg constants retained.
// *************************************************************************
{
  fp_t *op=new fp_t [nlambda],gff=1.1; // Berger, ApJ 1956
  memset(op,0,nlambda*sizeof(fp_t));
  fp_t fct=(e*e*c*h*h*Ryd)/(6*pi*pi*pi*e0*me*c)*sqrt((2.0*pi)/(3.0*k*me*me*me));
  fp_t a=fct*gff*ne/(sqrt(T)*c*c*c);
  for(int z=1;z<=Z;++z) for(int l=0;l<nlambda;++l) op[l]+=a*sqr((fp_t)z)*pop[x1i][x2i][x3i].N[z];
  fp_t b=-(h*c)/(k*T); // stimulated emission correction
  for(int l=0;l<nlambda;++l) op[l]*=lambda[l]*sqr(lambda[l])*(1.0-exp(b/lambda[l]));
  return op;
}

fp_t *atom::boundfree_em(fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i) // emissivity
// *************************************************************************
// * Bound-Free recombination...                                           *
// * Only recombination to particles in the ground state is considered.    *
// * This is considered reasonable since recombination with an excited     *
// * particle is usually followed by immediate autoionization (not always) *
// *************************************************************************
{
  fp_t *em=new fp_t [nlambda],gbf=1.0;
  memset(em,0,nlambda*sizeof(fp_t));
//
  for(int z=0;z<Z;++z){ // final stage cannot be ionized
    for(int l=0;l<nl[z];++l){
      if(fp_t *sigma=bf[z][l]->U(lambda,nlambda)){ // probability
        for(int i=0;i<nlambda;++i) em[i]+=pop[x1i][x2i][x3i].n[z+1][0]*sigma[i]; // emissivity=probability*number-density
        delete[] sigma; // this is a waste: we may want to keep it around...
      }
    }
  }
  return em;
}

fp_t *atom::boundfree_op(fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i) // absorption
// *************************************************************************
// * Bound-Free opacities as specified for each level...                   *
// *************************************************************************
{
  fp_t *op=new fp_t [nlambda];
  memset(op,0,nlambda*sizeof(fp_t));
//
  for(int z=0;z<Z;++z){ // final stage cannot be ionized
    for(int l=0;l<nl[z];++l){
      if(fp_t *sigma=bf[z][l]->U(lambda,nlambda)){ // cross-section
        for(int i=0;i<nlambda;++i) op[i]+=pop[x1i][x2i][x3i].n[z][l]*sigma[i]; // opacity=cross-section*number-density
        delete[] sigma; // this is a waste: we may want to keep it around...
      }
    }
  }
//
  return op;
}

fp_t *atom::boundbound_em(fp_t Temp_in,fp_t Ne_in,fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i) // absorption
{
  fp_t *em=new fp_t [nlambda];
  memset(em,0,nlambda*sizeof(fp_t));
  for(int z=0;z<Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
//
        fp_t a=damp_rad(A[z],ii,i);
        a+=damp_col(lam,Temp_in,Ne_in,ee[z][ii],ee[z][i],z);
//
        fp_t dld=broad_dop(lam,Temp_in,0.0); // last arg=vturb..
        a/=4.0*pi*c*dld/(lam*lam);           // a=a/(4 pi dnu_D)
//
        fp_t Aul=A[z][i][ii];
        for(int l=0;l<nlambda;++l){
          fp_t ndw=(lambda[l]-lam)/dld;
          em[l]+=h*c*fvoigt(ndw,a)*pop[x1i][x2i][x3i].n[z][i]*Aul/(4.0*pi*lambda[l]);
        }
      }
  return em;
}

fp_t *atom::boundbound_op(fp_t Temp_in,fp_t Ne_in,fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i) // absorption
{
  fp_t *op=new fp_t [nlambda];
  memset(op,0,nlambda*sizeof(fp_t));
  for(int z=0;z<Z;++z)
    for(int i=1;i<nl[z];++i) // upper level
      for(int ii=0;ii<i;++ii){ // lower level
//        fp_t lam=(1.0+Vr/c)*h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
        fp_t lam=h*c/(ee[z][i]-ee[z][ii]);   // transition wavelength: may be shifted by the local velocity
//
        fp_t a=damp_rad(A[z],ii,i);
        a+=damp_col(lam,Temp_in,Ne_in,ee[z][ii],ee[z][i],z);
        fp_t dld=broad_dop(lam,Temp_in,0.0); // last arg=vturb..
        a/=4.0*pi*c*dld/(lam*lam);           // a=a/(4 pi dnu_D)
//
        fp_t gr=((fp_t)(g[z][i])/(fp_t)(g[z][ii]));
        fp_t Blu=B[z][ii][i],Bul=B[z][i][ii];
        for(int l=0;l<nlambda;++l){
          fp_t ndw=(lambda[l]-lam)/dld;
          op[l]+=dld*fvoigt(ndw,a)*(Blu*pop[x1i][x2i][x3i].n[z][ii]-Bul*pop[x1i][x2i][x3i].n[z][i])*(h*c/(4.0*pi*lambda[l]));
        }
      }
  return op;
}

fp_t *atom::opacity(fp_t T,fp_t ne,fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i)
{
  fp_t *op=freefree_op(T,ne,lambda,nlambda,x1i,x2i,x3i);             // free-free opacity
  op=add(op,boundfree_op(lambda,nlambda,x1i,x2i,x3i),nlambda);       // bound-free ionization
  op=add(op,boundbound_op(T,ne,lambda,nlambda,x1i,x2i,x3i),nlambda); // bound-bound transitions
  return op;
}

fp_t *atom::emissivity(fp_t T,fp_t ne,fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i)
{
  fp_t *em=freefree_em(T,ne,lambda,nlambda,x1i,x2i,x3i);
  em=add(em,boundfree_em(lambda,nlambda,x1i,x2i,x3i),nlambda);
  em=add(em,boundbound_em(T,ne,lambda,nlambda,x1i,x2i,x3i),nlambda); // bound-bound transitions
  return em;
}


// ----------------------------------------------------------------------------------------------------------------------------

fp_t atom::boundbound_op(uint08_t z,uint16_t ll,uint16_t lu,fp_t T,fp_t Ne,fp_t Vr,fp_t Vt,fp_t lambda,struct pps &p)
{
  fp_t lam=h*c/(ee[z][lu]-ee[z][ll]); // transition wavelength: may be shifted by the local velocity
  fp_t Blu=B[z][ll][lu],Bul=B[z][lu][ll];

  fp_t dld=broad_dop(lam,T,Vt);                     // Doppler (thermal) broadening (last arg=vturb..)
  fp_t sf=(lam*lam)/(4.0*pi*c*dld);                 // sf=1.0/(4 pi dnu_D)= lam^2/(4 pi c dlam_D)

  fp_t ar=damp_rad(A[z],ll,lu);                     // radiative (natural) broadening
  fp_t ac=damp_col(lam,T,Ne,ee[z][ll],ee[z][lu],z); // collisional broadening
  fp_t a=sf*(ar+ac);                                // a=sf*gamma

  fp_t x=(lambda-(1.0+Vr/c)*lam)/dld;               // Doppler shifted wavelength in Doppler widths...
  fp_t prof=(lam*lam/c)*fvoigt(x,a)/dld;                        // normalized to the Doppler width...
  return prof*(Blu*p.n[z][ll]-Bul*p.n[z][lu])*(h*c/(4.0*pi*lam));
//  return prof;
//  fprintf(stderr,"boundbound_op: %d %d %d %E %E %E %E %E\n",z,ll,lu,x,a,fvoigt(x,a),lambda,lam);
//  fp_t fct=me*c/(8.0*pi*pi*e*e);
//  fp_t f=fct*lam*lam*A[z][lu][ll];
//  return fvoigt(x,a)*(Blu*p.n[z][ll]-Bul*p.n[z][lu])*(sqrt(pi)*e*e*f/(me*c*dld));
}


// **********************************************************************
// * routines dealing with the level populations                        *
// **********************************************************************

void atom::lte(fp_t ***T,fp_t ***Ne)
{ // Saha/Boltzmann distribution
 
  for(int z=0;z<=Z;++z)
    for(int x1i=x1l;x1i<=x1h;++x1i)
      for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i){
          fp_t U=partf[z]->U(T[x1i][x2i][x3i],Ne[x1i][x2i][x3i],io);
          for(int l=0;l<nl[z];++l){ 
            pop[x1i][x2i][x3i].n[z][l]=pop[x1i][x2i][x3i].N[z]*((g[z][l])?(fp_t)(g[z][l])
            *exp(-ee[z][l]/(k*T[x1i][x2i][x3i]))/U:1.0);

          }
    }        
}

void atom::lte(fp_t T, fp_t Ne, int x1i, int x2i, int x3i){

  for(int z=0;z<=Z;++z){
    fp_t U=partf[z]->U(T,Ne,io);
    for(int l=0;l<nl[z];++l){ 
      pop[x1i][x2i][x3i].n[z][l]=pop[x1i][x2i][x3i].N[z]*((g[z][l])?(fp_t)(g[z][l]) * exp(-ee[z][l]/(k*T))/U:1.0);
    }
  }
}

void atom::compute_active_population(fp_t *** T, fp_t *** Ne){

for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        pop[x1i][x2i][x3i].Na = 0.0;
      }

  for(int z=0;z<=Z;++z)
    for(int x1i=x1l;x1i<=x1h;++x1i)
      for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i){
          fp_t U=partf[z]->U(T[x1i][x2i][x3i],Ne[x1i][x2i][x3i],io);
          for(int l=0;l<nl[z];++l){ 
            pop[x1i][x2i][x3i].Na += pop[x1i][x2i][x3i].N[z]*((g[z][l])?(fp_t)(g[z][l])
            *exp(-ee[z][l]/(k*T[x1i][x2i][x3i]))/U:1.0);
          }
        }
}

fp_t atom::derivative_active_population(int x1i, int x2i, int x3i){

  // We only want to find the derivative of the active population in this point.

  fp_t delta_T = 1.0;
  fp_t derivative = 0.0;
  
  // Ok first we perturb the temperature.
  fp_t local_T = fetch_temperature(x1i, x2i, x3i);
  
  parent_atm->set_Temp(x1i, x2i, x3i, local_T + 0.5 * delta_T);
  parent_atm->execute_chemeq_for_point(x1i, x2i, x3i);
  fp_t local_ne = fetch_Ne(x1i, x2i, x3i);
  for(int z=0;z<=Z;++z){
    fp_t U=partf[z]->U(local_T+0.5*delta_T,local_ne,io);
    for(int l=0;l<nl[z];++l) 
      derivative += pop[x1i][x2i][x3i].N[z]*((g[z][l])?(fp_t)(g[z][l])
      *exp(-ee[z][l]/(k*(local_T+0.5*delta_T)))/U:1.0);
  }
  parent_atm->set_Temp(x1i, x2i, x3i, local_T - 0.5 * delta_T);
  parent_atm->execute_chemeq_for_point(x1i, x2i, x3i);
  local_ne = fetch_Ne(x1i, x2i, x3i);
  for(int z=0;z<=Z;++z){
    fp_t U=partf[z]->U(local_T-0.5*delta_T,local_ne,io);
    for(int l=0;l<nl[z];++l)
      derivative -= pop[x1i][x2i][x3i].N[z]*((g[z][l])?(fp_t)(g[z][l])
      *exp(-ee[z][l]/(k*(local_T-0.5*delta_T)))/U:1.0);
  }
  parent_atm->set_Temp(x1i, x2i, x3i, local_T);
  parent_atm->execute_chemeq_for_point(x1i, x2i, x3i); 
  
  derivative /= delta_T;

  return derivative;
}

void atom::ionfrc(fp_t Temp,fp_t Ne_in,fp_t *&fk,fp_t *&df,int &nk)
// **********************************************************************
// * return the ionization fractions and the derivative to the electron *
// * number density. This is easy in LTE but we need a NLTE formalism   *
// * to compute fk and df...                                            *
// **********************************************************************
{
  nk=Z+1;
  if(!fk) fk=new fp_t [nk];                        // allocate appropriate length array
  if(!df) df=new fp_t [nk];
//
  fp_t fct=2.0*pow((2.0*pi*me*k*Temp)/(h*h),1.5);  // constant factor
  fp_t Up=partf[0]->U(Temp,Ne_in,io);
  fp_t dUp=partf[0]->dU(Temp,Ne_in,io)/Up;
//
  fp_t fsum=(fk[0]=1.0);                           // population fraction N[i]/N[0]
  fp_t dfsum=(df[0]=0.0);                          // partial derivative d fk[i]/d Ne
//
  for(int z=0;z<Z;++z){
    fp_t Un=partf[z+1]->U(Temp,Ne_in,io);          // partition function for the next stage
    fp_t dUn=partf[z+1]->dU(Temp,Ne_in,io)/Un;     // partial derivative d fk[i]/d Ne
    fp_t f=fct*(Un/Up)*exp(-ip[z]/(k*Temp))/Ne_in; // ratio between successive stages
    fp_t fd=f*((dUn)-(dUp*Ne_in+1.0)/Ne_in);       // derivative of ratio between successive stages
    fsum+=(fk[z+1]=fk[z]*f);                       // factor w.r.t. the ground state
    dfsum+=(df[z+1]=(df[z]*f+fk[z]*fd));           // use the chain rule: (f*g)'=f'*g+f*g' and neglect the [weak?] dependency of Un/Up and ip[z] on Ne
    Up=Un;                                         // save the partition function for the next iteration
    dUp=dUn;                                       // save the partition function derivative for the next iteration
  }
  for(int z=0;z<=Z;++z){
    df[z]=(df[z]-dfsum*fk[z]/fsum)/fsum;           // chain rule
    fk[z]/=fsum;                                   // sum_k f_k = 1
  }
}

fp_t atom::damp_rad(fp_t **Az,uint16_t ll,uint16_t ul)
{ 
  fp_t gu=0.0,gl=0.0;
  for(uint16_t l=0;l<ll;++l){
    gu+=Az[ul][l];
    gl+=Az[ll][l];
  }
  for(uint16_t l=ll;l<ul;++l) gu+=Az[ul][l];
  return (gu+gl); 
}

// Fresh, brand new, collisional damping method, accounting for Van der Waals and resonance broadening

fp_t atom::damp_col(int ix1, int ix2, int ix3, int z, int i_from, int i_to, fp_t Temp, fp_t Ne, fp_t lambda_0){

  // When you are bored, do things right, because mistakes come easy:
  // int ix1, ix2, ix3 = coordinates of the point, we will need it in order to fetch relevant populations
  // int z             = ionization state
  // int i_from, i_to  = indices of levels involved into the transition.
  // Temperature? Electron density? Let's keep them to be consistent with other functions but in general we do not need it 
  // fp_t Temp         = local temperature
  // fp_t Ne           = local electron density, might be needed
  // fp_t lambda_0     = line center of the transition

  // First we will account for Van der Waals broadening. 

  fp_t alpha = alpha_col_dam[z][i_from][i_to];

  fp_t red_mass = 1.66E-24 / (1.0/1.008 + 1.66E-24/mass); // Reduced mass which goes into equation for mean relative velocity
  fp_t vmean = sqrt(8.0*k*Temp/pi/red_mass); // Mean relative velocity
  fp_t w = fetch_population(ix1, ix2, ix3, 0, 0, 0) * vmean * pow(vmean/1E6, -alpha) * col_dam_cross_section[z][i_from][i_to];
  
  //if (Z == 2)
    //printf("Z = %d d = %d from = %d to = %d  w = %e colliders = %e \n",Z, ix3, i_from, i_to,  w, fetch_population(ix1, ix2, ix3, 0, 0, 0));

  // And then we add resonance broadening
  // Ground level of the same ion. Should be, if I am not mistaken? 
  w += fetch_population(ix1, ix2, ix3, Z-1, z, 0) * 4.0 * e * e / mass * lambda_0 / c;
  
  return w;
}

int atom::compute_damp_col(int z, int lu, int ll){

  // This one computes cross-section for collisions with electrons

  // z  - ionization state
  // lu - index of the upper level
  // ll - index of the lower level 
  // In come parts of the code we use to and from in some we use upper and lower. Keep this in mind

  // These should be, of course, the result of the interpolation.

  if (Z == 1 && z == 0) { // Hydrogen
    // For hydrogen we can (more or less safely) assume that the dominant transition is the one satisfying:
    // l_qn_lower = ll - 1 and l_qn_upper = lu-1
    l_qn[z][ll] = ll-1;
    l_qn[z][lu] = lu-1;

  }

  fp_t n_eff_l = sqrt((z+1) * (z+1) * 2.18E-11 / (ip[z] - ee[z][ll]));
  fp_t n_eff_u = sqrt((z+1) * (z+1) * 2.18E-11 / (ip[z] - ee[z][lu]));

  //printf("%f %f \n", n_eff_l, n_eff_u);

  // Now these are upper and lower but we also have to find which one is s which one is p and which one is s
  if (l_qn[z][lu] + l_qn[z][ll] == 1) // sp transition
  	interpolate_col_damp_sp(n_eff_u, n_eff_l, z, lu, ll);
  else if (l_qn[z][lu] + l_qn[z][ll] == 3) // pd transition
  	interpolate_col_damp_pd(n_eff_u, n_eff_l, z, lu, ll);
  else if (l_qn[z][lu] + l_qn[z][ll] == 5) // df transition
  	interpolate_col_damp_sp(n_eff_u, n_eff_l, z, lu, ll);
  else{
  	// We will use Van der Waals approximation as given in Stell Atm. 3rd edition.
    // page 239, eq. 8.55 and table 8.1
  	alpha_col_dam[z][lu][ll] = alpha_col_dam[z][ll][lu] = 0.0;
  	col_dam_cross_section[z][lu][ll] = col_dam_cross_section[z][ll][lu] = 0.0;
  }

  return 0;
}

int atom::interpolate_col_damp_sp(fp_t n_eff_u, fp_t n_eff_l, int z, int lu, int ll){

	// Cross-section table:
	fp_t sigma_sp_table[21][18] = {{126,   140,   165,  202,  247,  299,  346,  383,  435,  491,  553,  617,  685,  769,  838,  925, 1011, 1082},
                                 {140,   150,   162,  183,  218,  273,  327,  385,  440,  501,  557,  620,  701,  764,  838,  923, 1025, 1085},
                                 {154,   167,   175,  192,  216,  251,  299,  357,  423,  487,  549,  617,  684,  759,  834,  910, 1014, 1064},
                                 {166,   180,   192,  206,  226,  253,  291,  339,  397,  459,  532,  600,  676,  755,  832,  896, 1002, 1055},
                                 {208,   194,   207,  223,  242,  265,  296,  335,  384,  445,  511,  583,  656,  726,  817,  889,  988, 1044},
                                 {262,   254,   220,  239,  261,  283,  310,  344,  388,  442,  496,  568,  635,  725,  791,  890,  970, 1036},
                                 {311,   306,   299,  251,  280,  304,  330,  361,  396,  443,  500,  563,  630,  704,  796,  880,  951, 1033},
                                 {358,   359,   350,  338,  293,  323,  352,  381,  416,  455,  511,  566,  635,  706,  780,  859,  946, 1039},
                                 {411,   409,   405,  392,  370,  340,  375,  406,  439,  478,  525,  580,  644,  714,  790,  873,  961, 1050},
                                 {462,   463,   459,  450,  443,  400,  394,  432,  467,  501,  546,  595,  650,  711,  786,  873,  963, 1050},
                                 {522,   525,   529,  524,  516,  518,  438,  454,  495,  532,  565,  621,  671,  741,  813,  874,  951, 1034},
                                 {589,   593,   590,  583,  579,  568,  565,  483,  517,  560,  600,  644,  691,  752,  821,  904,  978, 1048},
                                 {658,   655,   666,  657,  649,  653,  649,  587,  549,  592,  674,  674,  728,  782,  833,  902,  992, 1084},
                                 {738,   742,   747,  725,  721,  729,  699,  730,  626,  622,  668,  721,  765,  809,  887,  938, 1001, 1109},
                                 {838,   838,   810,  809,  790,  800,  769,  815,  757,  679,  704,  755,  806,  854,  901,  974, 1034, 1105},
                                 {942,   946,   925,  901,  918,  895,  919,  897,  933,  890,  785,  797,  859,  908,  976, 1020, 1115, 1173},
                                 {1059,  1061,  1056, 1061, 1074, 1031, 1036, 1036,  993, 1038,  932,  852,  878,  943, 1003, 1074, 1131, 1200},
                                 {1069,  1076,  1083, 1095, 1102, 1091, 1126, 1156, 1103, 1149, 1157, 1036,  972, 1007, 1064, 1124, 1209, 1283},
                                 {1338,  1350,  1356, 1354, 1324, 1301, 1312, 1318, 1257, 1239, 1297, 1233, 1089, 1059, 1106, 1180, 1218, 1317},
                                 {1409,  1398,  1367, 1336, 1313, 1313, 1409, 1354, 1317, 1287, 1353, 1386, 1279, 1158, 1141, 1188, 1260, 1335},
                                 {1328,  1332,  1342, 1369, 1405, 1451, 1502, 1524, 1506, 1477, 1522, 1594, 1572, 1436, 1328, 1325, 1382, 1446}};
  // Alpha table:
  fp_t alpha_sp_table[21][18] = {{.268, .269, .335, .377, .327, .286, .273, .270, .271, .268, .267, .264, .264, .264, .261, .256, .248, .245},
                                {.261, .256, .254, .282, .327, .355, .321, .293, .287, .271, .267, .273, .270, .270, .268, .268, .264, .263},
                                {.266, .264, .257, .252, .267, .289, .325, .339, .319, .301, .292, .284, .281, .281, .277, .282, .276, .274},
                                {.262, .274, .258, .251, .247, .254, .273, .291, .316, .322, .320, .302, .294, .290, .287, .292, .283, .277}, 
                                {.322, .275, .264, .259, .250, .245, .273, .255, .271, .284, .294, .308, .296, .299, .288, .289, .282, .278},
                                {.267, .300, .260, .268, .254, .242, .243, .242, .239, .246, .267, .277, .280, .290, .282, .281, .274, .271},
                                {.259, .274, .275, .252, .265, .248, .249, .237, .283, .236, .247, .254, .254, .271, .268, .267, .258, .262},
                                {.260, .255, .268, .268, .268, .264, .248, .239, .229, .240, .236, .234, .238, .244, .252, .251, .244, .255},
                                {.255, .255, .244, .247, .317, .246, .255, .244, .237, .231, .227, .231, .235, .232, .235, .241, .237, .245},
                                {.256, .254, .254, .249, .227, .319, .253, .253, .240, .237, .238, .233, .231, .230, .228, .234, .227, .241},
                                {.257, .254, .252, .235, .253, .240, .284, .251, .246, .241, .235, .228, .222, .225, .225, .219, .228, .233},
                                {.244, .240, .245, .238, .248, .230, .283, .252, .244, .244, .238, .235, .234, .236, .228, .224, .225, .231},
                                {.244, .241, .244, .237, .237, .249, .219, .324, .239, .245, .242, .242, .232, .233, .221, .227, .231, .218},
                                {.241, .245, .249, .239, .243, .250, .217, .254, .308, .237, .247, .244, .234, .228, .233, .224, .227, .226},
                                {.243, .243, .232, .227, .235, .253, .227, .220, .320, .270, .243, .252, .248, .238, .234, .241, .225, .227},
                                {.225, .226, .234, .230, .226, .233, .249, .225, .216, .300, .286, .237, .240, .247, .243, .234, .231, .238},
                                {.268, .260, .247, .238, .233, .241, .254, .248, .207, .227, .315, .260, .226, .237, .240, .239, .239, .240},
                                {.248, .246, .238, .226, .213, .221, .226, .226, .204, .194, .248, .316, .234, .216, .236, .233, .221, .230},
                                {.200, .202, .198, .194, .206, .207, .227, .224, .207, .185, .198, .275, .315, .233, .229, .231, .233, .236},
                                {.202, .209, .221, .226, .230, .245, .202, .257, .246, .225, .215, .246, .320, .321, .244, .239, .251, .253}, 
                                {.246, .248, .255, .265, .274, .285, .292, .284, .273, .250, .225, .239, .295, .352, .320, .258, .260, .269}};

    // Convert them to fp_t **
    fp_t ** cs_table = ft2dim(0, 20, 0, 17);
    fp_t ** alpha_table = ft2dim(0, 20, 0, 17);
    for (int i = 0; i<21; ++i)
      for (int ii = 0; ii<18; ++ii){
        cs_table[i][ii] = sigma_sp_table[i][ii];
        alpha_table[i][ii] = alpha_sp_table[i][ii];
      }
    fp_t s_state[21];
  	fp_t p_state[18];
  	for (int i = 0; i<21; ++i)
  		s_state[i] = 1.0 + i*0.1;
  	for (int i = 0; i<18; ++i)
  		p_state[i] = 1.3 + i*0.1;              
    // First we determine which one is s which one is p
    fp_t p_for_interpolation, s_for_interpolation;

    if (l_qn[z][lu] == 1){
    	p_for_interpolation = n_eff_u;
    	s_for_interpolation = n_eff_l;
    }
    else {
    	s_for_interpolation = n_eff_u;
    	p_for_interpolation = n_eff_l;
    }

    //printf("s_for interpolation = %f p_for interpolation = %f \n", s_for_interpolation, p_for_interpolation);

    
    fp_t cross_section = interpol_2d(cs_table, s_state, p_state, 21, 18, s_for_interpolation, p_for_interpolation);
    fp_t alpha = interpol_2d(alpha_table, s_state, p_state, 21, 18, s_for_interpolation, p_for_interpolation);

    //printf("Z = %d Cross section = %e Alpha = %e \n",Z, cross_section, alpha);

    cross_section *= 5.62E-17 * pow((4.0 / pi), alpha/2.0) * tgamma(2.0 - alpha / 2.0);

    alpha_col_dam[z][lu][ll] = alpha_col_dam[z][ll][lu] = alpha;
    col_dam_cross_section[z][lu][ll] = col_dam_cross_section[z][ll][lu] = cross_section;
    del_ft2dim(cs_table, 0, 20, 0, 17);
    del_ft2dim(alpha_table, 0, 20, 0, 17);

    //printf("sp!\n");
    
	return 0;
}

int atom::interpolate_col_damp_pd(fp_t n_eff_u, fp_t n_eff_l, int z, int lu, int ll){

  // Cross-section table:
  fp_t sigma_pd_table[18][18] = {{425,  461,  507,  566,  630,  706,  799,  889,  995, 1083, 1191, 1334, 1478, 1608, 1790, 1870, 1936, 2140},
                                  {429,  460,  505,  565,  633,  704,  795,  896,  985, 1082, 1199, 1340, 1487, 1611, 1795, 1872, 1937, 2136},
                                  {419,  451,  501,  556,  627,  700,  785,  891,  977, 1088, 1212, 1346, 1493, 1604, 1793, 1863, 1930, 2144},
                                  {402,  437,  489,  544,  614,  695,  779,  875,  975, 1102, 1221, 1350, 1488, 1591, 1774, 1844, 1919, 2126},
                                  {384,  418,  467,  529,  595,  674,  769,  856,  976, 1108, 1224, 1338, 1467, 1570, 1743, 1817, 1900, 2118},
                                  {366,  397,  443,  505,  576,  651,  755,  841,  973, 1095, 1210, 1308, 1435, 1545, 1702, 1786, 1878, 2081},
                                  {356,  387,  432,  489,  562,  635,  722,  841,  961, 1078, 1175, 1273, 1397, 1517, 1672, 1763, 1863, 2034},
                                  {359,  388,  431,  479,  545,  624,  707,  834,  943, 1059, 1158, 1256, 1368, 1490, 1647, 1747, 1849, 1998},
                                  {361,  394,  436,  483,  547,  615,  704,  817,  920, 1027, 1124, 1238, 1358, 1465, 1624, 1736, 1838, 1978},
                                  {400,  382,  440,  489,  546,  610,  690,  817,  897,  998, 1115, 1201, 1351, 1453, 1599, 1728, 1829, 1953},
                                  {474,  461,  416,  491,  549,  612,  701,  806,  883,  974, 1078, 1194, 1310, 1456, 1569, 1716, 1818, 1925},
                                  {531,  518,  507,  463,  547,  615,  694,  784,  881,  958, 1047, 1153, 1297, 1432, 1547, 1688, 1809, 1901},
                                  {594,  585,  577,  564,  513,  615,  695,  779,  879,  949, 1041, 1145, 1264, 1388, 1544, 1644, 1804, 1879},
                                  {675,  659,  651,  639,  632,  576,  695,  782,  879,  957, 1046, 1141, 1254, 1391, 1524, 1614, 1793, 1871},
                                  {739,  734,  726,  719,  715,  708,  663,  776,  901,  971, 1022, 1117, 1232, 1355, 1478, 1616, 1766, 1887},
                                  {819,  821,  805,  784,  773,  761,  736,  761,  888,  958, 1044, 1145, 1237, 1346, 1487, 1614, 1721, 1891},
                                  {899,  895,  871,  852,  856,  861,  854,  759,  883,  984, 1027, 1113, 1226, 1355, 1467, 1568, 1703, 1885},
                                  {973,  946,  955,  925,  939,  927,  902,  920,  870,  987, 1061, 1145, 1234, 1319, 1439, 1552, 1722, 1859}};
  // Alpha table:
  fp_t alpha_pd_table[18][18] = {{.281, .288, .283, .282, .278, .281, .272, .274, .268, .257, .251, .243, .246, .251, .254, .268, .304, .308},  
                                  {.290, .297, .291, .290, .286, .282, .277, .275, .267, .254, .252, .244, .250, .257, .260, .274, .308, .312}, 
                                  {.294, .299, .293, .294, .288, .289, .281, .276, .265, .256, .251, .247, .258, .264, .268, .283, .318, .317},  
                                  {.297, .298, .302, .300, .289, .295, .290, .276, .264, .256, .260, .258, .268, .277, .281, .292, .330, .327},  
                                  {.305, .311, .313, .315, .305, .304, .299, .279, .271, .272, .273, .276, .285, .290, .293, .302, .340, .340},  
                                  {.292, .294, .303, .305, .301, .307, .290, .277, .274, .278, .287, .288, .295, .302, .306, .312, .343, .346},  
                                  {.268, .277, .279, .285, .285, .290, .279, .278, .280, .283, .295, .296, .305, .310, .313, .315, .342, .346},  
                                  {.288, .285, .280, .278, .278, .277, .272, .271, .279, .288, .297, .305, .310, .313, .311, .310, .335, .338},  
                                  {.314, .304, .292, .282, .275, .275, .262, .272, .290, .293, .299, .307, .308, .310, .303, .302, .325, .328},  
                                  {.346, .329, .313, .295, .283, .275, .264, .274, .288, .302, .307, .310, .306, .307, .292, .296, .315, .320},  
                                  {.320, .295, .326, .318, .294, .277, .275, .271, .293, .303, .305, .309, .309, .303, .294, .294, .310, .313},  
                                  {.304, .310, .297, .320, .317, .297, .283, .274, .298, .305, .308, .311, .313, .300, .290, .293, .305, .306},  
                                  {.314, .313, .308, .297, .325, .314, .293, .276, .292, .309, .314, .308, .303, .296, .286, .291, .301, .302},  
                                  {.308, .311, .307, .312, .288, .340, .305, .285, .294, .310, .315, .309, .296, .285, .281, .288, .298, .295},  
                                  {.313, .310, .315, .303, .313, .294, .331, .286, .294, .307, .320, .316, .303, .281, .278, .285, .290, .292},  
                                  {.315, .306, .308, .297, .295, .283, .334, .297, .280, .294, .314, .321, .313, .291, .280, .279, .287, .290},  
                                  {.308, .304, .305, .297, .279, .285, .251, .278, .278, .284, .297, .314, .307, .289, .274, .274, .274, .291}, 
                                  {.301, .299, .298, .285, .265, .279, .241, .285, .260, .286, .302, .306, .302, .288, .277, .263, .271, .293}};

    // Convert them to fp_t **
    fp_t ** cs_table = ft2dim(0, 17, 0, 17);
    fp_t ** alpha_table = ft2dim(0, 17, 0, 17);
    for (int i = 0; i<18; ++i)
      for (int ii = 0; ii<18; ++ii){
        cs_table[i][ii] = sigma_pd_table[i][ii];
        alpha_table[i][ii] = alpha_pd_table[i][ii];
      }
    fp_t p_state[18];
    fp_t d_state[18];
    for (int i = 0; i<18; ++i)
      p_state[i] = 1.3 + i*0.1;
    for (int i = 0; i<18; ++i)
      d_state[i] = 2.3 + i*0.1;              

    // First we determine which one is s which one is p
    fp_t p_for_interpolation, d_for_interpolation;

    if (l_qn[z][lu] == 2){
      d_for_interpolation = n_eff_u;
      p_for_interpolation = n_eff_l;
    }
    else {
      p_for_interpolation = n_eff_u;
      d_for_interpolation = n_eff_l;
    }

   
    fp_t cross_section = interpol_2d(cs_table, p_state, d_state, 18, 18, p_for_interpolation, d_for_interpolation);
    fp_t alpha = interpol_2d(alpha_table, p_state, d_state, 18, 18, p_for_interpolation, d_for_interpolation);

    cross_section *= 5.62E-17 * pow((4.0 / pi), alpha/2.0) * tgamma(2.0 - alpha / 2.0);

    alpha_col_dam[z][lu][ll] = alpha_col_dam[z][ll][lu] = alpha;
    col_dam_cross_section[z][lu][ll] = col_dam_cross_section[z][ll][lu] = cross_section;
   
    // Now time to write a 2D interpolation procedure:
    del_ft2dim(cs_table, 0, 17, 0, 17);
    del_ft2dim(alpha_table, 0, 17, 0, 17);
    
  return 0;
}

int atom::interpolate_col_damp_df(fp_t n_eff_u, fp_t n_eff_l, int z, int lu, int ll){

  // Cross-section table:
  fp_t sigma_df_table[18][18] = {{808,  873,  958, 1059, 1175, 1306, 1453, 1615, 1793, 1979, 2121, 2203, 2461, 2604, 2764, 2757, 2784, 3156},
                                  {798,  866,  953, 1052, 1172, 1299, 1450, 1606, 1776, 1967, 2114, 2196, 2451, 2601, 2763, 2767, 2783, 3142},
                                  {781,  848,  934, 1030, 1149, 1276, 1416, 1596, 1751, 1944, 2100, 2188, 2436, 2594, 2767, 2777, 2795, 3123},
                                  {766,  831,  915, 1010, 1124, 1239, 1398, 1564, 1729, 1912, 2083, 2180, 2426, 2585, 2776, 2790, 2808, 3106},
                                  {750,  814,  897,  987, 1097, 1201, 1355, 1530, 1718, 1875, 2060, 2171, 2414, 2575, 2779, 2809, 2820, 3103},
                                  {733,  797,  872,  950, 1049, 1166, 1326, 1502, 1670, 1851, 2026, 2165, 2396, 2562, 2779, 2827, 2832, 3099},
                                  {726,  786,  853,  936, 1011, 1128, 1303, 1472, 1649, 1844, 1979, 2159, 2371, 2548, 2778, 2840, 2848, 3103},
                                  {709,  783,  847,  912, 1002, 1093, 1270, 1419, 1606, 1787, 1951, 2139, 2335, 2533, 2775, 2847, 2863, 3104},
                                  {758,  721,  838,  907, 1010, 1066, 1211, 1401, 1600, 1774, 1972, 2098, 2313, 2528, 2781, 2857, 2892, 3121},
                                  {869,  882,  820,  870, 1003, 1098, 1165, 1368, 1527, 1735, 1896, 2030, 2288, 2534, 2776, 2844, 2902, 3123},
                                  {970,  967,  934,  938,  918, 1130, 1194, 1287, 1507, 1679, 1821, 2021, 2271, 2525, 2732, 2786, 2882, 3085},
                                 {1079, 1043, 1056, 1007, 1014, 1021, 1200, 1326, 1424, 1668, 1818, 1988, 2242, 2493, 2672, 2719, 2853, 3035},
                                 {1174, 1173, 1127, 1154, 1104, 1099, 1169, 1288, 1442, 1580, 1704, 1882, 2136, 2400, 2561, 2648, 2832, 2994},
                                 {1285, 1278, 1269, 1225, 1252, 1229, 1116, 1343, 1380, 1594, 1710, 1874, 2054, 2309, 2484, 2607, 2813, 2932},
                                 {1440, 1408, 1422, 1380, 1383, 1341, 1361, 1192, 1448, 1454, 1675, 1873, 2069, 2246, 2432, 2610, 2811, 2878},
                                 {1572, 1545, 1553, 1517, 1481, 1502, 1469, 1349, 1373, 1561, 1586, 1781, 2072, 2301, 2490, 2626, 2754, 2832},
                                 {1698, 1701, 1694, 1641, 1617, 1651, 1566, 1600, 1374, 1547, 1698, 1749, 1989, 2289, 2511, 2594, 2689, 2774},
                                 {1870, 1841, 1786, 1752, 1777, 1757, 1666, 1732, 1522, 1533, 1707, 1817, 1928, 2194, 2435, 2574, 2665, 2742}};
  // Alpha table:
  fp_t alpha_df_table[18][18] = {{.295, .286, .299, .300, .307, .310, .311, .311, .316, .319, .325, .351, .364, .369, .372, .379, .373, .351},
                                  {.295, .295, .301, .302, .311, .316, .314, .314, .320, .321, .324, .349, .361, .365, .368, .374, .368, .349},  
                                  {.286, .298, .302, .304, .311, .323, .321, .319, .324, .323, .323, .345, .355, .358, .362, .367, .361, .343}, 
                                  {.290, .295, .307, .316, .322, .329, .326, .325, .329, .324, .321, .343, .350, .351, .354, .360, .358, .337},  
                                  {.292, .299, .307, .321, .327, .336, .333, .330, .330, .320, .321, .338, .344, .344, .345, .352, .352, .332},  
                                  {.291, .299, .309, .323, .335, .339, .335, .333, .327, .323, .319, .333, .336, .336, .336, .344, .345, .329}, 
                                  {.297, .302, .312, .321, .340, .338, .333, .327, .325, .319, .318, .324, .329, .330, .330, .336, .337, .325}, 
                                  {.319, .314, .317, .327, .334, .344, .339, .327, .323, .318, .312, .318, .319, .322, .322, .326, .327, .316},  
                                  {.333, .328, .339, .325, .359, .351, .332, .325, .322, .311, .309, .310, .311, .316, .314, .317, .321, .313}, 
                                  {.274, .273, .323, .412, .318, .339, .359, .328, .324, .311, .309, .325, .322, .315, .318, .319, .325, .314},  
                                  {.297, .296, .273, .302, .436, .325, .354, .335, .326, .311, .314, .330, .323, .324, .325, .323, .330, .314}, 
                                  {.284, .295, .296, .280, .300, .438, .322, .348, .332, .318, .320, .332, .335, .334, .335, .331, .333, .309}, 
                                  {.280, .278, .285, .297, .279, .320, .445, .319, .320, .324, .328, .338, .348, .346, .345, .336, .328, .300}, 
                                  {.280, .273, .267, .273, .284, .268, .343, .390, .323, .308, .318, .325, .343, .348, .346, .337, .311, .286},  
                                  {.277, .270, .260, .266, .276, .263, .294, .408, .337, .324, .299, .308, .331, .334, .345, .327, .315, .280}, 
                                  {.270, .262, .258, .260, .273, .273, .262, .375, .410, .298, .312, .294, .313, .331, .328, .322, .307, .270},  
                                  {.271, .267, .262, .264, .274, .269, .261, .323, .351, .359, .294, .325, .310, .318, .321, .315, .291, .268}, 
                                  {.275, .276, .272, .276, .279, .270, .264, .295, .393, .340, .319, .287, .320, .330, .316, .302, .280, .261}};

    // Convert them to fp_t **
    fp_t ** cs_table = ft2dim(0, 17, 0, 17);
    fp_t ** alpha_table = ft2dim(0, 17, 0, 17);
    for (int i = 0; i<18; ++i)
      for (int ii = 0; ii<18; ++ii){
        cs_table[i][ii] = sigma_df_table[i][ii];
        alpha_table[i][ii] = alpha_df_table[i][ii];
      }
    fp_t d_state[18];
    fp_t f_state[18];
    for (int i = 0; i<18; ++i)
      d_state[i] = 1.3 + i*0.1;
    for (int i = 0; i<18; ++i)
      f_state[i] = 2.3 + i*0.1;              

    // First we determine which one is s which one is p
    fp_t d_for_interpolation, f_for_interpolation;

    if (l_qn[z][lu] == 3){
      f_for_interpolation = n_eff_u;
      d_for_interpolation = n_eff_l;
    }
    else {
      d_for_interpolation = n_eff_u;
      f_for_interpolation = n_eff_l;
    }
    
    fp_t cross_section = interpol_2d(cs_table, d_state, f_state, 18, 18, d_for_interpolation, f_for_interpolation);
    fp_t alpha = interpol_2d(alpha_table, d_state, f_state, 18, 18, d_for_interpolation, f_for_interpolation);

    cross_section *= 5.62E-17 * pow((4.0 / pi), alpha/2.0) * tgamma(2.0 - alpha / 2.0);

    alpha_col_dam[z][lu][ll] = alpha_col_dam[z][ll][lu] = alpha;
    col_dam_cross_section[z][lu][ll] = col_dam_cross_section[z][ll][lu] = cross_section;
   
    // Now time to write a 2D interpolation procedure:
    del_ft2dim(cs_table, 0, 17, 0, 17);
    del_ft2dim(alpha_table, 0, 17, 0, 17);
    
  return 0;
}

fp_t atom::broad_dop(fp_t lambda_0,fp_t Temp_in,fp_t vturb_in)
{
  return sqrt(2.0 * k*Temp_in/mass+ vturb_in*vturb_in)*lambda_0/c;
}

fp_t *gen_points(fp_t ll,fp_t ul,int np) // generate wavelength points
{
  fp_t *p=new fp_t [np],dl=(ul-ll)/(fp_t)(np-1);
  for(int l=0;l<np;++l) p[l]=ll+dl*(fp_t)l;
  return p;
}

fp_t *add_points(fp_t *a,int32_t &n,fp_t *b,int32_t m) // add points to array
{
//
// must stay sorted and observe similar step size to before
//
  fp_t *cc=new fp_t [n+m];       // allocate new array
  if(n){
    memcpy(cc,a,n*sizeof(fp_t)); // copy existing points
    delete[] a;
  }
  memcpy(cc+n,b,m*sizeof(fp_t));
  n+=m;
  delete[] b;
  return cc;
}

fp_t *atom::getlambda(fp_t *lambda,int32_t &nlambda,fp_t Temp_in,fp_t Nt_in,fp_t Ne_in)
// *******************************************************
// * Check what wavelengths are required by the          *
// * atomic model and add those points to the wavelength *
// * array if necessary (based on existing points)       *
// * it is important not to overdo this, only the        *
// * accuracy of the rates should be considered, not the *
// * spectrum (this is an important distinction for      *
// * strong lines with extended wings)                   *
// *******************************************************
{
  if(1){                                       // check if NLTE levels present <-- This does not actually check!
    int32_t nnlambda=0;                        // calculate the number of required wavelength points
    fp_t *tmplamb=new fp_t [nnlambda+1];         // compute sorted wavelength array
    fp_t fct=me*c/(8.0*pi*pi*e*e);             // common factor for line transitions

// add points, avoiding overlap and duplicates
    for(int z=0;z<Z;++z)
      for(int ii=0;ii<nl[z]-1;++ii){ // lower level
// b-f frequency points?
        if(bf[z][ii]){
          int np=0;
          if(fp_t *tmp=bf[z][ii]->getlambda(np)) tmplamb=add_points(tmplamb,nnlambda,tmp,np);
        }
// b-b frequency points
        for(int i=ii+1;i<nl[z];++i){                           // for each upper level
          
          fp_t lam=h*c/(ee[z][i]-ee[z][ii]);                   // transition wavelength: may be shifted by the local velocity, which we did not account for here!
          
          fp_t gf=((fp_t)(g[z][i]) / g[z][ii])*fct*lam*lam*A[z][i][ii]; // Oscillator strength which we need why ? 
          osc_str[z][i][ii] = osc_str[z][ii][i] = gf;
          
          fp_t a=damp_rad(A[z],ii,i);                          // Take radiative damping.
          a+=damp_col(x1h, x2h, x3h, z,i, ii, Temp_in, Ne_in, lam); // Add to it collisional damping
          fp_t dld_max=broad_dop(lam,10.0E4,5E5);  // thermal broadening for some maximum parameters
          fp_t dld_min=broad_dop(lam,4.0E4,1E5);   // thermal broadening for some minimum parameters  
          a/=4.0*pi*dld_min*c/(lam*lam);                           // a=a/(4 pi dnu_D)=a lam^2 /(4 pi c dl_D)
//       
          a = 1E-3; // 'Cannonical value'
          // We compute how many doppler widths away from the center we want to be. Maximum is 250. We stop when the profile reachers E-15.
          // Better criterion should be relative opacity of continuum and the line, but ok, this is not so bad. 
          
          fp_t ndw,eps=-10;
          for(ndw=0.0;ndw<=250.0;ndw+=1.0) if(log10(fvoigt(ndw,a))<eps) break;
//
          // We then compute, profile values in the center and two outermost points:
          fp_t *vv=new fp_t [3],*ll=new fp_t [3];
          ll[0]=-ndw;
          ll[1]=0.0;
          ll[2]=ndw;
          vv[0]=log10(fvoigt(ll[0],a));
          vv[1]=log10(fvoigt(ll[1],a));
          vv[2]=log10(fvoigt(ll[2],a));
//
          // So far so good, what happens now, is that we try to somehow assemble wavelengths inside.
          int np=3,inc;
          fp_t epss=0.75,epsx=0.1*exp(vv[1]/0.43429448);
          
          do{
            inc=0;
            for(int jj=1;jj<np;++jj) if((fabs(vv[jj-1]-vv[jj])>epss)||(fabs(exp(vv[jj-1]/0.43429448)-exp(vv[jj]/0.43429448))>epsx)) ++inc;
            fp_t *tvv=new fp_t [np+inc],*tll=new fp_t [np+inc];
            for(int idx=0,jj=1;jj<np;++jj,++idx){
              if((fabs(vv[jj-1]-vv[jj])>epss)||(fabs(exp(vv[jj-1]/0.43429448)-exp(vv[jj]/0.43429448))>epsx)){
                tll[idx]=ll[jj-1];
                tvv[idx]=vv[jj-1];
                ++idx;
                tll[idx]=0.5*(ll[jj-1]+ll[jj]);
                tvv[idx]=0.43429448*log(fvoigt(tll[idx],a));
              }else{
                tll[idx]=ll[jj-1];
                tvv[idx]=vv[jj-1];
              }
              tll[idx+1]=ll[jj];
              tvv[idx+1]=vv[jj];
            }
            np+=inc;
            delete[] vv;
            delete[] ll;
            vv=tvv;
            ll=tll;            
          }while(inc);
          delete[] vv;

          // Debug mode
          delete []ll;

          // Until we understand the above, we make our own grid. With blackjack.
          // There will be a conversion here:
          ndw = 70; // Seventy times dld_max;
          fp_t res = 0.2; // Resolution of 0.2 dld_min;
          np = 2 * int(ndw * dld_max / dld_min / 0.2) + 1.0;
          ll = new fp_t [np];
          for (int jj = 0; jj<np;++jj)
            ll[jj] = -ndw + 2.0*jj * ndw / fp_t(np-1);


         // printf("Ok, number of points is: %d \n", np);
          //for (int jj=0;jj<np;++jj){
           // printf("%d %e \n",jj,ll[jj]);
          //}

          for(int jj=0;jj<np;++jj){
            //printf("%d %f \n", jj, ll[jj]);
            ll[jj]=lam+ll[jj]*dld;           
          } 

          io.msg(IOL_INFO,"%s %d %d %9.2f A [%9.2f - %9.2f] A ndw=%E a=%E A=%E log(gf)=%6.3f\n",id,ii,i,1E8*lam,1E8*ll[0],1E8*ll[np-1],ndw,a,A[z][i][ii],log10(gf));
          io.msg(IOL_DEB1,"%s %d %d %9.2f A (%E %E %E %E) %E\n",id,ii,i,1E8*lam,6.669E+15*gf*0.25/(1E16*lam*lam),damp_rad(A[z],ii,i),damp_col(x1h, x2h, x3h, z,i, ii, Temp_in, Ne_in, lam),dld,gf);
          
          tmplamb=add_points(tmplamb,nnlambda,ll,np);
        }
      }
// merge with the new wavelength array
    lambda=add_points(lambda,nlambda,tmplamb,nnlambda);
  }
  sort(lambda-1,nlambda);
  return lambda;
}

void atom::popsetup(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
{
  if(!pop){
    x1l=x1l_in;
    x1h=x1h_in;
    x2l=x2l_in;
    x2h=x2h_in;
    x3l=x3l_in;
    x3h=x3h_in;
//
    int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
//
    pop=new pps** [nx1] - x1l;
    pop[x1l]=new pps* [nx1*nx2] - x2l;
    pop[x1l][x2l]=new pps [nx1*nx2*nx3] - x3l;
    for(int x2=x2l+1;x2<=x2h;++x2) pop[x1l][x2]=pop[x1l][x2-1]+nx3;
    for(int x1=x1l+1;x1<=x1h;++x1){
      pop[x1]=pop[x1-1]+nx2;
      pop[x1][x2l]=pop[x1-1][x2l]+nx2*nx3;
      for(int x2=x2l+1;x2<=x2h;++x2) pop[x1][x2]=pop[x1][x2-1]+nx3;
    }
//
    int32_t nn=0;
    for(int i=0;i<=Z;++i) nn+=nl[i];
//
    for(int x1i=x1l;x1i<=x1h;++x1i)
      for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i){
          pop[x1i][x2i][x3i].N=new fp_t [Z+1];
          pop[x1i][x2i][x3i].n=new fp_t* [Z+1];
          pop[x1i][x2i][x3i].n[0]=new fp_t [nn];
          for(int i=1;i<=Z;++i) pop[x1i][x2i][x3i].n[i]=pop[x1i][x2i][x3i].n[i-1]+nl[i-1];
        }
  }
}

void atom::popclean(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
{
  delete[] (pop[x1l][x2l]+x3l);
  delete[] (pop[x1l]+x2l);
  delete[] (pop+x1l);
//
  pop=0;
  x1l=x1h=x2l=x2h=x3l=x3h=0;
}

void atom::pupdate(fp_t N_in,fp_t *fk,int,int32_t x1i,int32_t x2i,int32_t x3i)
{
  pop[x1i][x2i][x3i].Nt=N_in; // total free number density
  for(int z=0;z<=Z;++z) pop[x1i][x2i][x3i].N[z]=pop[x1i][x2i][x3i].Nt*fk[z];
}

uint08_t atom::rtsetup(fp_t *theta,fp_t *phi,int ntp,fp_t *lambda,int32_t nlambda,
                       int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
// *********************************************************
// * Setup contributions to rate equations and approximate *
// * lambda operator                                       *
// *********************************************************
{
  if((x1l!=x1l_in)||(x1h!=x1h_in)||(x2l!=x2l_in)||(x2h!=x2h_in)||(x3l!=x3l_in)||(x3h!=x3h_in)) return EC_DIM_MISMATCH;
// allocate
  tmap=new uint32_t** [Z+1];
  ntr=0;
  for(int z=0;z<=Z;++z){
    tmap[z]=ui32t2dim(0,nl[z],0,nl[z]);
    memset(tmap[z][0],0,sqr(nl[z]+1)*sizeof(uint32_t));
    for(int i=0;i<nl[z];++i){
      if((~flags[z][i])&FL_NONE){ // only explicitly specified levels...
        for(int ii=i+1;ii<nl[z];++ii) if(A[z][ii][i]) tmap[z][i][ii]=tmap[z][ii][i]=++ntr; // radiative transition
        tmap[z][i][nl[z]]=tmap[z][nl[z]][i]=++ntr; // bound-free "transition"
      }
    }
  }

  // While we are here let us set up rmap, just to have it.
  rmap = new uint16_t * [Z+1];
  for (int z=0;z<=Z;++z)
    rmap[z] = new uint16_t[nl[z]];

  for (int i=0;i<nmap;++i){
    int z = zmap[i];
    int l = lmap[i];
    rmap[z][l] = i;
  }


  // What are we going to need is also inverse mapping, to be completely clear we have to do it only once, and here
  if (ntr){

    inverse_tmap = ui32t2dim(1, ntr, 1, 4);
    int ntr_temp = 0;
   
    // Search the matrix:
    for (int z = 0; z<=Z; ++z)
      for (int i = 0; i<nl[z]; ++i){ // Lower level

        if((~flags[z][i])&FL_NONE){ // only explicitly specified levels...
          for (int ii = i+1; ii<nl[z]; ++ii){ // Upper level

            if (A[z][ii][i]) ++ntr_temp;
            //printf("%d \n", ntr_temp);
            inverse_tmap[ntr_temp][1] = ntr_temp;
            inverse_tmap[ntr_temp][2] = z;
            inverse_tmap[ntr_temp][3] = i;
            inverse_tmap[ntr_temp][4] = ii;
          }
          // Bound-free transition:
          ++ntr_temp;
          //printf("%d \n", ntr_temp);
          inverse_tmap[ntr_temp][1] = ntr_temp;
          inverse_tmap[ntr_temp][2] = z;
          inverse_tmap[ntr_temp][3] = i; // Lower level
          inverse_tmap[ntr_temp][4] = nl[z]; // Uper level
        }

      }
  }

  if(ntr){ // at least one level explicitly defined
    Jb=ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
    Ju=ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
    Ls=ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
    norm = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);

  // And additional thingies:
    current_profile = ft4dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, ntr);
    //F_old = ft4dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, ntr-1);
    
    // Broyden:

    //Broyden = ft5dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, ntr-1, 1, ntr-1);
    //memset(Broyden[x1l][x2l][x3l][1] + 1, 0, (x3h - x3l + 1) * (ntr-1) * (ntr-1) *sizeof(fp_t));

  }else{
   Jb=Ju=Ls=current_profile=norm=0;
   //Broyden=0;
   //F_old = 0;
   inverse_tmap = 0;
 }
  return EC_OK;
}

uint08_t atom::rtclean(int ntp,int32_t nlambda,
                       int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
// *********************************************************
// * Setup contributions to rate equations and approximate *
// * lambda operator                                       *
// *********************************************************
{
// cleanup 
  if(Jb) del_ft4dim(Jb,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  if(Ju) del_ft4dim(Ju,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  if(Ls) del_ft4dim(Ls,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  if(norm) del_ft4dim(norm,x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);
  if(current_profile) del_ft4dim(current_profile, x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);

  //if(Broyden) del_ft5dim(Broyden, x1l, x1h, x2l, x2h, x3l, x3h, 1, ntr-1, 1, ntr-1);
  
  //if (F_old) del_ft4dim(F_old, x1l, x1h, x2l, x2h, x3l, x3h, 1, ntr-1);

  if(inverse_tmap) del_ui32t2dim(inverse_tmap, 1, ntr, 1, 4);

  for(int z=0;z<=Z;++z) del_ui32t2dim(tmap[z],0,nl[z],0,nl[z]);
  delete[] tmap;
  return EC_OK;
}

void atom::rtinit(void)
{
// set to 0: start of new transition integration cycle
  if(Jb) memset(Jb[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
  if(Ju) memset(Ju[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
  if(Ls) memset(Ls[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
  if(norm) memset(norm[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));
}

void atom::prof_init(void){
  if(current_profile) memset (current_profile[x1l][x2l][x3l]+1,0, (x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*ntr*sizeof(fp_t));

}

void atom::prof_setup(void){

tmap=new uint32_t** [Z+1];
ntr=0;
  for(int z=0;z<=Z;++z){
    tmap[z]=ui32t2dim(0,nl[z],0,nl[z]);
    memset(tmap[z][0],0,sqr(nl[z]+1)*sizeof(uint32_t));
    for(int i=0;i<nl[z];++i){
      if((~flags[z][i])&FL_NONE){ // only explicitly specified levels...
        for(int ii=i+1;ii<nl[z];++ii) if(A[z][ii][i]) tmap[z][i][ii]=tmap[z][ii][i]=++ntr; // radiative transition
        tmap[z][i][nl[z]]=tmap[z][nl[z]][i]=++ntr; // bound-free "transition"
      }
    }
}

if (ntr)
  current_profile = ft4dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, ntr);
else 
  current_profile = 0;

}

void atom::prof_clear(void){

  if (current_profile)
    del_ft4dim(current_profile, x1l,x1h,x2l,x2h,x3l,x3h,1,ntr);


}

void atom::popnorm(){

  // Renormalize all the populations with respect to N_t 
  // THIS IS NOT USED ANY MORE. (At least as much as I know)

  for (int x1i = x1l; x1i <= x1h; ++x1i)
    for (int x2i = x2l; x2i <= x2h; ++x2i)
      for (int x3i = x3l; x3i <= x3h; ++x3i){

        fp_t actual_sum = 0.0;
        for (int i = 0; i<nmap; ++i)
          actual_sum += pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]];

        fp_t normalize = actual_sum / pop[x1i][x2i][x3i].Nt;

        for (int i=0;i<nmap;++i)
          pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] /= normalize;
      } 
}

fp_t atom::pops(atmol **atm,uint16_t natm,fp_t Temp,fp_t ne,int32_t x1i,int32_t x2i,int32_t x3i)
// *********************************************************
// * Solve the rate equations for given temperature and    *
// * electron density. The ionization fractions and        *
// * derivative to the electron density must be returned   *
// * by ionfrc                                             *
// *********************************************************

// In this particular function we do it by means of MALI as given in Rybicki & Hummer. Equations 
// of SE are linear with respect to level populations.

{

  int x3i_control = x3h+1;

  //if (x3i == x3i_control)
  //for (int i = 1; i<=ntr; ++i)
    //io.msg(IOL_INFO, "Tansition = %d from z =  %d, i = %d to j = %d, J_lu = %e J_ul = %e \n", i, inverse_tmap[i][2], inverse_tmap[i][3], inverse_tmap[i][4], Jb[x1i][x2i][x3i][i], Ju[x1i][x2i][x3i][i]);

  if(Jb){
    fp_t *J_lu=Jb[x1i][x2i][x3i]; // integrated intensity
    fp_t *J_ul=Ju[x1i][x2i][x3i];
    fp_t *L=Ls[x1i][x2i][x3i]; // approximate lambda operator
    fp_t *nrm = norm[x1i][x2i][x3i]; // norm of the profile: <---- Who thought that THIS would be thing you need to save? Frequency-by-frequency approach is a bitch

    fp_t **M=ft2dim(1,nmap,1,nmap); // Allocate the rate matrix
    memset(M[1]+1,0,nmap*nmap*sizeof(fp_t)); // Set all elements to zero

    fp_t* b = new fp_t[nmap];
    memset(b,0,nmap*sizeof(fp_t));
  
    // Setup the rate equation:
    for(int i=0;i<nmap-1;++i){ // For all the 'levels' but the last one 
      uint16_t l=lmap[i]; // "current lvl
      uint08_t z=zmap[i]; // apropriate ionization stage

      // Here we plot these rates/intensities:

      // First account for b-b transitions within the particular ionization stage:
      for(int ll=0; ll<nl[z]; ++ll){ // For all the levels:

        fp_t JJ = (tmap[z][l][ll]) ? J_lu[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] :-1.0;   // angular and frequency integrated intensity
        fp_t LL = (tmap[z][l][ll]) ? L[tmap[z][l][ll]] / nrm[tmap[z][l][ll]] : 0.0;
   
        // Transitions from this level
        fp_t Radiative_rates = 1.0 * R_ij_local_ALO(z, l, ll, JJ, LL, pop[x1i][x2i][x3i].n[z][l], pop[x1i][x2i][x3i].n[z][ll]);
        fp_t Collisional_rates = C_ij(z, l, ll, Temp, ne);
        //if(l>ll)
        //  printf("%e\n", Collisional_rates);
        //Collisional_rates += C_ij_H(z, l, ll, Temp, fetch_population(x1i, x2i, x3i, 0, 0, 0)); // Modify for resonant broadening:

        M[i+1][i+1] -= (Radiative_rates + 1.0 * Collisional_rates); 
        
        // Transitions to this level:
        Radiative_rates = 1.0 * R_ij_local_ALO(z, ll, l, JJ, LL, pop[x1i][x2i][x3i].n[z][ll], pop[x1i][x2i][x3i].n[z][l]);
        Collisional_rates = C_ij(z, ll, l, Temp, ne);
        //Collisional_rates += C_ij_H(z, ll, l, Temp, fetch_population(x1i, x2i, x3i, 0, 0, 0));
        int dl = ll - l;

        M[i+1][i+1+dl] += (Radiative_rates + Collisional_rates);  
      }
      // Then account for b-f and f-b transitions:
      if(int dl=nl[z]-l){              // ground level of next ionization stage (if not this level)
        
        // ionization out of level z,l
        fp_t JJ=(tmap[z][l][nl[z]])?J_lu[tmap[z][l][nl[z]]]:-1.0;     // angular and frequency integrated intensity
        
        // Transitions from this level:
        fp_t Radiative_rates = 0.0 * R_i_cont(z, l, JJ, Temp);
        fp_t Collisional_rates = C_i_cont(z, l, Temp, ne);

        M[i+1][i+1] -= (Radiative_rates + 1.0 * Collisional_rates);
        // But now the same these rates should populate the continuum:
        M[i+1+dl][i+1] += (Radiative_rates + 1.0 * Collisional_rates);
        
        JJ=(tmap[z][l][nl[z]])?J_ul[tmap[z][l][nl[z]]]:-1.0;     // angular and frequency integrated intensity
        
        Radiative_rates = 0.0 * R_cont_i(z, l, JJ, Temp, ne);
        Collisional_rates = C_cont_i(z, l, Temp, ne);
        
        M[i+1][i+1+dl] += (Radiative_rates + 1.0 * Collisional_rates);
        // And these same rates should depopulate the continuum:
        M[i+1+dl][i+1+dl] -= (Radiative_rates + 1.0 * Collisional_rates);
      }
  }
  for(int ii=1;ii<=nmap;++ii) M[nmap][ii]=1.0;  // final equation: sum(n)=cst (number of particles remains the same)
  for (int i =0;i<nmap;++i) b[nmap-1] += pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]];
  
  // Print the rate matrix:
  
  if (x3i == x3i_control){
    printf("Rate matrix:\n");
    //io.msg(IOL_INFO,"atom::pops: %s\n",name);
    for(int i=1;i<=nmap;++i){
      for(int ii=1;ii<=nmap;++ii){
        fprintf(stderr,"%6.6E ",M[i][ii]);
      
      }
      fprintf(stderr,"%6.6E \n", b[i-1]);
    }
  }

  // Let us invert the matrix: 
  fp_t * M_to_solve = M[1] +1;
  fp_t * M_LU = new fp_t [nmap * nmap];
  fp_t * solution = new fp_t [nmap];  
  Crout(nmap,M_to_solve, M_LU);
  solveCrout(nmap,M_LU,b,solution);
  // Correct:
  fp_t delta = 0.0;
  for (int i = 0; i<nmap; ++i){
    fp_t rel_delta = fabs(solution[i] - pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]]) / pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]];
    delta = (rel_delta > delta) ? rel_delta : delta;
  }
  for (int i = 0; i<nmap; ++i)
    pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] = solution[i];

  // Now here you can immediately get the response, if you want:
  // This agrees perfectly with the scheme from atom::responses. You can keep this here, but basically, it is not necessary. 
  //memset(b,0,nmap*sizeof(fp_t));
  //for (int i=0;i<nmap-1;++i)// For each i
    //for (int ii=0;ii<nmap;++ii) // For each j
      //b[i] += pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]] * derivative_collisions_full_temp(x1i,x2i,x3i,i,ii) - 
        //pop[x1i][x2i][x3i].n[zmap[ii]][lmap[ii]] * derivative_collisions_full_temp(x1i,x2i,x3i,ii,i);

  // Last one, of course needs to be the derivative of total population:
  //b[nmap-1] = derivative_active_population(x1i,x2i,x3i);
  //solveCrout(nmap,M_LU,b,solution);
  //for (int i=0;i<nmap;++i)
  //  printf("%d %d %e \n", x3i, i, solution[i]/pop[x1i][x2i][x3i].n[zmap[i]][lmap[i]]);
  
  
  delete []M_LU;
  delete []b;
  delete []solution;
  return delta;
  }
return 0.0;
}


// *********************************************************
// * modify angular/frequency integration                  *
// *********************************************************
fp_t **atom::weight(fp_t op,fp_t T,fp_t Ne,fp_t Vr,fp_t Vt,struct pps &pp,fp_t lambda,fp_t t,fp_t p) // purely geometric stuff here
{
  fp_t **w=ft2dim(0,0,1,ntr);
  memset(w[0]+1,0,ntr*sizeof(fp_t));
  for(int i=0;i<=0;++i)
    for(uint08_t z=0;z<Z;++z)
      for(uint16_t ll=0;ll<nl[z];++ll){ // lower level
        for(uint16_t ul=ll+1;ul<nl[z];++ul) w[i][tmap[z][ll][ul]]=boundbound_op(z,ll,ul,T,Ne,Vr,Vt,lambda,pp)/op; // bound-bound
        w[i][tmap[z][ll][nl[z]]]=boundfree_op(z,ll,Vr,lambda,pp)/op; // bound-free
//        if(tmap[z][ll][nl[z]]==10) io.msg(IOL_INFO|IOL_FATAL,"%d,%E %E %E\n",ll,w[i][tmap[z][ll][nl[z]]],lambda,pp.n[z][ll]);
      }
  return w;
}

void atom::add(fp_t ***I,fp_t ***L,fp_t ***op,fp_t ***T,fp_t ***Ne,fp_t ***Vr,fp_t bin,fp_t lambda,fp_t t,fp_t p) // I-only version
// *********************************************************
// * Add intensity and approximate lambda operator         *
// * contribution for this atom to the relevant transition *
// * weighted intensity                                    *
// *********************************************************
// *********************************************************
// * For atomic polarization, the changes are confined     *
// * to this class (unless strong coupling between species *
// * exists)                                               *
// *********************************************************
{
//
// should the radiative quantities be reordered for this? <-- Milic: I do not understand what does this mean.
//
  int ncmp=1; // Total number of Stokes components
// there are multiple contributions here for the angle dependent quantities?
  if(Jb&&Ls){
    for(int x1i=x1l;x1i<=x1h;++x1i)
      for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i){
          
          // These are obviously integrating weigths for the mean intensity
          fp_t **w=weight(op[x1i][x2i][x3i],T[x1i][x2i][x3i],Ne[x1i][x2i][x3i],Vr[x1i][x2i][x3i],0.0,pop[x1i][x2i][x3i],lambda,t,p);

          fp_t *Jm=Jb[x1i][x2i][x3i];                  // integrated/redistributed intensity array
          fp_t *Lm=Ls[x1i][x2i][x3i];                  // approximate lambda array
          for(int i=0;i<ncmp;++i)                      // Stokes index
            for(int tr=1;tr<=ntr;++tr){                // transition index
              //Jm[tr]+=w[i][tr]*bin;
              Jm[tr]+=w[i][tr] * I[x1i][x2i][x3i];
              Lm[tr]+=w[i][tr]*bin*L[x1i][x2i][x3i];
            }
          del_ft2dim(w,0,0,1,ntr);
        }
  }
}

void atom::add(fp_t *** I, fp_t *** L, fp_t *** opp, fp_t lambda, fp_t lambda_w, fp_t angular_weight){

  // Ok Ivan, you wrote this more or less from the scratch. So see if you really want it to look like this.

  int ncmp = 1; // total number of Stokes components

  // If this is the proper transition, i.e. if it has mean intensiy, approximate operator and `norm' 
  if (Jb && Ls && norm){
    for(int x1i=x1l;x1i<=x1h;++x1i)
      for(int x2i=x2l;x2i<=x2h;++x2i)
        for(int x3i=x3l;x3i<=x3h;++x3i){

          // Fetch them in the more convenient arrays:
          fp_t *Jm = Jb[x1i][x2i][x3i];
          fp_t *Jn = Ju[x1i][x2i][x3i];
          fp_t *Lm = Ls[x1i][x2i][x3i];
          fp_t *nrm = norm[x1i][x2i][x3i];

          // For each Stokes component and for each transition. 
          for (int i=0; i<ncmp; ++i)
            for (int tr=1; tr<=ntr; ++tr){

              // Add appropriate contribution to norm, in case the integral of the absorption profile is not one.
              nrm[tr] += lambda_w * angular_weight * current_profile[x1i][x2i][x3i][tr] * 0.25 / pi;

              // Now, for the general case which involves the presence of the continuum and the overlapping lines, the approximate lambda operator for the transition 'tr'
              // is actually first multiplied by the ratio of the opacity due to the transition in question to the total opacity. We need to extract the transition.

              int z_state = inverse_tmap[tr][2];
              int lower_level = inverse_tmap[tr][3];
              int upper_level = inverse_tmap[tr][4]; 

              if (upper_level < nl[z_state]){

                fp_t line_energy = ee[z_state][upper_level] - ee[z_state][lower_level];

                // This is a bit slow, no?
                fp_t line_opacity = (pop[x1i][x2i][x3i].n[z_state][lower_level] * B[z_state][lower_level][upper_level] - pop[x1i][x2i][x3i].n[z_state][upper_level] * B[z_state][upper_level][lower_level]) 
                  * current_profile[x1i][x2i][x3i][tr] * line_energy * 0.25 / pi;
                fp_t line_factor = line_opacity / opp[x1i][x2i][x3i];

                fp_t elementary_contribution = lambda_w * angular_weight * I[x1i][x2i][x3i] * current_profile[x1i][x2i][x3i][tr] * line_factor * 0.25 / pi;

                Jm[tr] += elementary_contribution;
                Jn[tr] += elementary_contribution;
                Lm[tr] += lambda_w * angular_weight * L[x1i][x2i][x3i] * current_profile[x1i][x2i][x3i][tr] * line_factor * 0.25 / pi ;
              }

              // If it is a bound-free transition. There is no profile. Instead, compute rates directly, so:

              else if (upper_level == nl[z_state]){
                fp_t sigma = (bf[z_state][lower_level]) ? bf[z_state][lower_level]->U(lambda) : 0.0;
                Jm[tr] += lambda_w * angular_weight * sigma * lambda / h /c * I[x1i][x2i][x3i];
                // This might need double-checking
                Jn[tr] += lambda_w * angular_weight * sigma * lambda / h /c * (I[x1i][x2i][x3i] + 2.0 * h * c * c / pow(lambda, 5.0)) * exp(-h * c / lambda / k / fetch_temperature(x1i, x2i, x3i));

              }
            }
          }
  }
}

void atom::info(void)
{
  for(int08_t z=0;z<=Z;++z){
    char blank[]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',0};
    char *r=roman(z+1);
    blank[10-strlen(r)-strlen(id)]=0;
//    io.msg(IOL_INFO,"atom::lte: %s %s%s N=%E\n",id,r,blank,N[z]);
//    for(int l=0;l<nl[z];++l) io.msg(IOL_INFO,"atom::lte: %s %s%s   n[%d]=%E\n",id,r,blank,l,n[z][l]);
    delete[] r;
  }
}

fp_t atom::get_pop(int x1i, int x2i, int x3i, int index_ion){
  return (index_ion > Z || index_ion < 0) ? 0 : pop[x1i][x2i][x3i].N[index_ion];
  //return pop[x1i][x2i][x3i].Nt;
}

fp_t atom::get_pop(int x1i, int x2i, int x3i, int index_ion, int index_lvl){
  return (index_ion > Z || index_ion < 0 || index_lvl < 0 || index_lvl >= nl[index_ion]) ? 0 : pop[x1i][x2i][x3i].n[index_ion][index_lvl];
  //return pop[x1i][x2i][x3i].Nt;
}

fp_t atom::get_J(int x1i, int x2i, int x3i, int transition){

  return Jb[x1i][x2i][x3i][transition];
}

fp_t atom::get_L(int x1i, int x2i, int x3i, int transition){

  return Ls[x1i][x2i][x3i][transition];
}

fp_t atom::get_norm(int x1i, int x2i, int x3i, int transition){

  return norm[x1i][x2i][x3i][transition];
}

fp_t atom::get_active_pop(int x1i, int x2i, int x3i){
  return pop[x1i][x2i][x3i].Na;
}


// =========================================================================================================================================================================================

// Old stuff:

fp_t atom::damp_col(fp_t lambda_0,fp_t Temp_in,fp_t Ne_in,fp_t El,fp_t Eu, fp_t N_6, int08_t z)
{

  printf("Warning: obsolete function fp_t atom::damp_col(fp_t lambda_0,fp_t Temp_in,fp_t Ne_in,fp_t El,fp_t Eu, fp_t N_6, int08_t z)!\n");
  return 0;
}
fp_t atom::damp_col(fp_t lambda_0,fp_t Temp_in,fp_t Ne_in,fp_t El,fp_t Eu, int08_t z)
{
  // Here we might have some approximate treatment.

  printf("Warning: obsolete function fp_t atom::damp_col(fp_t lambda_0,fp_t Temp_in,fp_t Ne_in,fp_t El,fp_t Eu, int08_t z) !\n");

  return 0;
}