#include <string.h>
#include <math.h>

#include "types.h"
#include "uts.h"
#include "const.h"
#include "io.h"
#include "struts.h"
#include "pack.h"
#include "mem.h"
#include "mathtools.h"

#include "molcfg.h"

#include "moleqc.h"

#include "mol.h"

atmol *mol_new(molcfg *cfg,atmol **atm,int na,io_class &io_in)
{
  if (!strcmp("H-", cfg->id))
    return new h_minus_mol(cfg, atm, na, io_in);
  
  return new mol(cfg,atm,na,io_in);
}

atmol *mol_new(uint32_t numid,uint08_t *buf,int32_t &offs,uint08_t do_swap,atmol **atm,int na,io_class &io_in)
{
  if (numid == 256)
    return new h_minus_mol(buf,offs,do_swap,atm,na,io_in);
  
  return new mol(buf,offs,do_swap,atm,na,io_in);
}

atmol *h_minus_mol_new(molcfg *cfg,atmol **atm,int na,io_class &io_in)
{
  return new h_minus_mol(cfg,atm,na,io_in);
}

atmol *h_minus_mol_new(uint32_t numid,uint08_t *buf,int32_t &offs,uint08_t do_swap,atmol **atm,int na,io_class &io_in)
{
  return new h_minus_mol(buf,offs,do_swap,atm,na,io_in);
}

mol::mol(molcfg *cfg,atmol **atm,int na,io_class &io_in):atmol(cfg->name,cfg->id,io_in)
{
  int len=0;
  nc=0;
  for(int i=0;cfg->cid[i];++i,++nc) len+=strlen(cfg->cid[i])+strlen(cfg->zid[i])+3;
  comp=new ion [nc+1];
//
  char *lhs=new char [len];
  lhs[0]=0;
  numid=0;
  for(int i=0;cfg->cid[i];++i){
    if(i) strcat(lhs," + ");
    strcat(lhs,cfg->cid[i]);
    strcat(lhs,cfg->zid[i]);

    comp[i].atom=0;                   // link to element...
    for(int a=0;a<na;++a)
      if(atm[a]->has_id(cfg->cid[i])){
        comp[i].atom=atm[a];
        break;
      }
//
    if(!comp[i].atom)
      if(strcmp(cfg->cid[i],"e")){
        const char* postfix[]={"st","nd","rd","th"};
        io_in.msg(IOL_ERROR|IOL_FATAL,"mol::mol: %s: %d%s component \"%s\" not found!\n",name,i+1,postfix[min(i,3)],cfg->cid[i]);
      }
//
    if((comp[i].z=unroman(cfg->zid[i]))>0) comp[i].z-=1;
    uint64_t cid=(comp[i].atom)?comp[i].atom->get_id():0;
    numid|=cid<<(8*(i+1));
  }
//
  comp[nc].atom=0;
  io_in.msg(IOL_INFO,"mol::mol: %s <=> %s\n",lhs,cfg->id);
//  numid=((((uint64_t)this)>>1)<<1)+1;
  io_in.msg(IOL_INFO,"mol::mol: %s numid => 0x%016lX\n",name,numid);
//
  if(cfg->nk){
    K=new eqc* [nk=cfg->nk];
    for(int i=0;i<cfg->nk;++i) K[i]=eqc_new(cfg->K[i],io);
  }else io_in.msg(IOL_ERROR|IOL_FATAL,"mol::mol: %s [%s]: no equilibrium constant specified!\n",name,id);
//
  delete[] lhs;
}

mol::mol(uint08_t *buf,int32_t &offs,uint08_t do_swap,atmol **atm,int na,io_class &io_in):atmol(buf,offs,do_swap,io_in)
{
  offs+=unpack(buf+offs,do_swap,atm,na,io_in);
}

mol::~mol(void)
{
  if(comp) delete[] comp;
  if(K){
    for(int i=0;i<nk;++i) delete K[i];
    delete[] K;
  }
}

int32_t mol::size(io_class &io_in)
{
  int32_t sz=atmol::size(io);
  sz+=sizeof(uint08_t);      // nc
  sz+=nc*sizeof(uint64_t);   // cid
  sz+=nc*sizeof(uint08_t);   // zid
// dissociation constant stuff
  sz+=sizeof(uint08_t);      // nk
  for(int08_t i=0;i<nk;++i) sz+=K[i]->size(io); // dissociation coefficient data
//
//  io.msg(IOL_INFO,"mol::size %d\n",sz);
//
  return sz;
}

int32_t mol::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
  int32_t offs=atmol::pack(buf,do_swap,io);
  offs+=::pack(buf+offs,nc);
  for(int i=0;i<nc;++i){
    uint64_t cid=(comp[i].atom)?comp[i].atom->get_id():0; // electrons have null pointer
    offs+=::pack(buf+offs,cid,do_swap);
    offs+=::pack(buf+offs,comp[i].z);
  }
// dissociation constant stuff
  offs+=::pack(buf+offs,nk);
  for(int08_t i=0;i<nk;++i) offs+=K[i]->pack(buf+offs,do_swap,io); // dissociation coefficient data
//
  return offs;
}

int32_t mol::unpack(uint08_t *buf,uint08_t do_swap,atmol **atm,int na,io_class &io_in)
{
// only unpack local stuff
  int32_t offs=::unpack(buf,nc);
  comp=new ion [nc];
  for(int i=0;i<nc;++i){
    uint64_t cid;
    offs+=::unpack(buf+offs,cid,do_swap);
    offs+=::unpack(buf+offs,comp[i].z);
    comp[i].atom=0; // if unchanged electron is assumed
    for(int a=0;a<na;++a)
      if(atm[a]->has_id(cid)){
        comp[i].atom=atm[a];
        break;
      }
  }
// dissociation constant stuff
  offs+=::unpack(buf+offs,nk);
  K=new eqc* [nk];
  for(int08_t i=0;i<nk;++i) K[i]=eqc_new(buf,offs,do_swap,io); // dissociation coefficient data
//
  return offs;
}

fp_t mol::dissoc(fp_t T,fp_t ne)
{
  return K[0]->K(T,ne);
}

void mol::pupdate(fp_t N_in,fp_t*,int)
{
  //printf("You are now running the pupdate from molecule class. \n");
  //N=N_in;
}

/*
H- BF               DOUG. FR
   ILOGL=0   KVADL=1   MINEX=0   MAXEX=0  NLATB=19
     00000      1000      1500      2000      3000      4000
      5000      6000      7000      8000      9000     10000
     11000     12000     13000     14000     15000     16000
     16421
   ILOGT=0   KVADT=0   MINET=0   MAXET=0  NTETB=1   ITETA=0
     0.000     0.390     0.629     0.912     1.567     2.219
     2.810     3.306     3.676     3.887     3.913     3.741
     3.377     2.842     2.176     1.437     0.708     0.124
     0.000      .         .         .         .         .
H- FF               S&C. J
   ILOGL=0   KVADL=1   MINEX=0   MAXEX=0  NLATB=22
     00000      3038      4556      5063      5695      6509
      7594      9113     10125     11391     13018     15188
     18225     22782     30376     45563     91127    200000
    400000    800000   1600000   3200000
   ILOGT=0   KVADT=1   MINET=0   MAXET=0  NTETB=16   ITETA=1
       0.5       0.6       0.7       0.8       0.9       1.0
       1.1       1.2       1.3       1.4       1.5       1.6
       1.7       1.8       1.9       2.0
  0.00 E00  3.44 E-2  7.70 E-2  9.59 E-2  1.21 E-1  1.56 E-1
  2.10 E-1  2.98 E-1  3.65 E-1  4.58 E-1  5.92 E-1  7.98 E-1
  1.14 E00  1.77 E00  3.10 E00  6.92 E00  2.72 E01  1.31 E02
  5.23 E02  2.09 E03  8.36 E03  3.35 E04
  0.00 E00  4.18 E-2  9.41 E-2  1.16 E-1  1.45 E-1  1.88 E-1
  2.53 E-1  3.59 E-1  4.39 E-1  5.50 E-1  7.11 E-1  9.58 E-1
  1.36 E00  2.11 E00  3.71 E00  8.27 E00  3.25 E01  1.56 E02
  6.26 E02  2.50 E03  1.00 E04  4.00 E04
  0.00 E00  4.91 E-2  1.10 E-1  1.35 E-1  1.69 E-1  2.18 E-1
  2.93 E-1  4.16 E-1  5.09 E-1  6.37 E-1  8.24 E-1  1.11 E00
  1.58 E00  2.44 E00  4.29 E00  9.56 E00  3.77 E01  1.81 E02
  7.24 E02  2.90 E03  1.16 E04  4.64 E04
  0.00 E00  5.65 E-2  1.25 E-1  1.53 E-1  1.92 E-1  2.47 E-1
  3.32 E-1  4.70 E-1  5.75 E-1  7.21 E-1  9.31 E-1  1.25 E00
  1.78 E00  2.75 E00  4.84 E00  1.08 E01  4.26 E01  2.05 E02
  8.19 E02  3.27 E03  1.31 E04  5.24 E04
  0.00 E00  6.39 E-2  1.40 E-1  1.72 E-1  2.14 E-1  2.76 E-1
  3.69 E-1  5.22 E-1  6.39 E-1  8.00 E-1  1.03 E00  1.39 E00
  1.98 E00  3.05 E00  5.37 E00  1.19 E01  4.73 E01  2.28 E02
  9.10 E02  3.64 E03  1.46 E04  5.82 E04
  0.00 E00  7.13 E-2  1.56 E-1  1.90 E-1  2.36 E-1  3.03 E-1
  4.06 E-1  5.73 E-1  7.00 E-1  8.76 E-1  1.13 E00  1.52 E00
  2.17 E00  3.34 E00  5.87 E00  1.31 E01  5.19 E01  2.49 E02
  9.97 E02  3.99 E03  1.60 E04  6.38 E04
  0.00 E00  7.87 E-2  1.71 E-1  2.08 E-1  2.58 E-1  3.31 E-1
  4.41 E-1  6.21 E-1  7.58 E-1  9.49 E-1  1.23 E00  1.65 E00
  2.34 E00  3.62 E00  6.36 E00  1.42 E01  5.63 E01  2.71 E02
  1.08 E03  4.33 E03  1.73 E04  6.92 E04
  0.00 E00  8.62 E-2  1.86 E-1  2.25 E-1  2.80 E-1  3.57 E-1
  4.75 E-1  6.68 E-1  8.15 E-1  1.02 E00  1.32 E00  1.77 E00
  2.52 E00  3.89 E00  6.83 E00  1.52 E01  6.06 E01  2.91 E02
  1.16 E03  4.66 E03  1.86 E04  7.45 E04
  0.00 E00  9.36 E-2  2.01 E-1  2.43 E-1  3.01 E-1  3.84 E-1
  5.09 E-1  7.15 E-1  8.71 E-1  1.09 E00  1.40 E00  1.89 E00
  2.68 E00  4.14 E00  7.28 E00  1.62 E01  6.47 E01  3.11 E02
  1.24 E03  4.97 E03  1.99 E04  7.96 E04
  0.00 E00  1.01 E-1  2.16 E-1  2.61 E-1  3.22 E-1  4.10 E-1
  5.43 E-1  7.60 E-1  9.25 E-1  1.15 E00  1.49 E00  2.00 E00
  2.84 E00  4.39 E00  7.72 E00  1.72 E01  6.87 E01  3.30 E02
  1.32 E03  5.28 E03  2.11 E04  8.45 E04
  0.00 E00  1.08 E-1  2.31 E-1  2.78 E-1  3.43 E-1  4.36 E-1
  5.76 E-1  8.04 E-1  9.77 E-1  1.22 E00  1.57 E00  2.11 E00
  3.00 E00  4.63 E00  8.14 E00  1.82 E01  7.26 E01  3.49 E02
  1.40 E03  5.58 E03  2.23 E04  8.93 E04
  0.00 E00  1.16 E-1  2.45 E-1  2.96 E-1  3.64 E-1  4.62 E-1
  6.08 E-1  8.47 E-1  1.03 E00  1.28 E00  1.65 E00  2.21 E00
  3.15 E00  4.86 E00  8.55 E00  1.91 E01  7.64 E01  3.67 E02
  1.47 E03  5.87 E03  2.35 E04  9.40 E04
  0.00 E00  1.23 E-1  2.60 E-1  3.13 E-1  3.85 E-1  4.87 E-1
  6.40 E-1  8.90 E-1  1.08 E00  1.34 E00  1.73 E00  2.32 E00
  3.29 E00  5.08 E00  8.95 E00  2.00 E01  8.01 E01  3.85 E02
  1.54 E03  6.16 E03  2.46 E04  9.85 E04
  0.00 E00  1.30 E-1  2.75 E-1  3.30 E-1  4.06 E-1  5.12 E-1
  6.72 E-1  9.32 E-1  1.13 E00  1.40 E00  1.80 E00  2.42 E00
  3.43 E00  5.30 E00  9.33 E00  2.09 E01  8.37 E01  4.02 E02
  1.61 E03  6.44 E03  2.57 E04  1.03 E05
  0.00 E00  1.38 E-1  2.89 E-1  3.47 E-1  4.26 E-1  5.37 E-1
  7.03 E-1  9.73 E-1  1.18 E00  1.46 E00  1.88 E00  2.51 E00
  3.57 E00  5.51 E00  9.71 E00  2.17 E01  8.73 E01  4.19 E02
  1.68 E03  6.71 E03  2.68 E04  1.07 E05
  0.00 E00  1.45 E-1  3.03 E-1  3.64 E-1  4.46 E-1  5.62 E-1
  7.34 E-1  1.01 E00  1.23 E00  1.52 E00  1.95 E00  2.61 E00
  3.70 E00  5.71 E00  1.01 E01  2.25 E01  9.07 E01  4.36 E02
  1.74 E03  6.97 E03  2.79 E04  1.11 E05
*/


// -----------------------------------------------------------------------------------------------------------------------------------------

// Below all this, some functions for H minus molecule. At least we try to do it:


h_minus_mol::h_minus_mol(molcfg *cfg,atmol **atm,int na,io_class &io_in):mol(cfg, atm, na, io_in)
{
  N = 0;
  }
//
 
h_minus_mol::h_minus_mol(uint08_t *buf,int32_t &offs,uint08_t do_swap,atmol **atm,int na,io_class &io_in):mol(buf, offs, do_swap, atm, na, io_in)
{
  N = 0;
}

// Here now popsetup:

void h_minus_mol::popsetup(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in){

  if(!N){
    x1l=x1l_in;
    x1h=x1h_in;
    x2l=x2l_in;
    x2h=x2h_in;
    x3l=x3l_in;
    x3h=x3h_in;
//
    int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
//
    N=new fp_t** [nx1] - x1l;
    N[x1l]=new fp_t* [nx1*nx2] - x2l;
    N[x1l][x2l]=new fp_t [nx1*nx2*nx3] - x3l;
    for(int x2=x2l+1;x2<=x2h;++x2) N[x1l][x2]=N[x1l][x2-1]+nx3;
    for(int x1=x1l+1;x1<=x1h;++x1){
      N[x1]=N[x1-1]+nx2;
      N[x1][x2l]=N[x1-1][x2l]+nx2*nx3;
      for(int x2=x2l+1;x2<=x2h;++x2) N[x1][x2]=N[x1][x2-1]+nx3;
    }
  }
}

//

void h_minus_mol::popclean(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in){
  delete[] (N[x1l_in][x2l_in]+x3l_in);
  delete[] (N[x1l_in]+x2l_in);
  delete[] (N+x1l_in);
//
  N=0;
}

//
void h_minus_mol::pupdate(fp_t N_in,fp_t*,int,int32_t x1i,int32_t x2i,int32_t x3i){
  N[x1i][x2i][x3i] =  N_in;
}

fp_t ***h_minus_mol::add(fp_t ***src,fp_t ***dst,int32_t ll1,int32_t ul1,int32_t ll2,int32_t ul2,int32_t ll3,int32_t ul3)
{
  if(src){
    int32_t nn=(ul1-ll1+1)*(ul2-ll2+1)*(ul3-ll3+1);
    for(int32_t i=0;i<nn;++i) dst[ll1][ll2][ll3+i]+=src[ll1][ll2][ll3+i];
    del_ft3dim(src,ll1,ul1,ll2,ul2,ll3,ul3);
  }
  return dst;
}

fp_t ***** h_minus_mol::add(fp_t *****src, fp_t *****dst, int32_t ll1, int32_t ul1, int32_t ll2, int32_t ul2, int32_t ll3, int32_t ul3, int32_t ll4, int32_t ul4, int32_t ll5, int32_t ul5){
  if(src){
    int32_t nn = (ul1 - ll1 + 1) * (ul2 - ll2 + 1) * (ul3 - ll3 + 1) * (ul4 - ll4 + 1) * (ul5 - ll5 + 1);
    for (int32_t i=0;i<nn;++i)
      dst[ll1][ll2][ll3][ll4][ll5+i] += src[ll1][ll2][ll3][ll4][ll5+i];
    del_ft5dim(src, ll1, ul1, ll2, ul2, ll3, ul3, ll4, ul4, ll5, ul5);
  }
  return dst;
}

// And now the ones for the opacity and emissivity:

fp_t *** h_minus_mol::emissivity(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){
  //return 0;
  fp_t *** em = opacity(T, Ne, Vlos, Vt, B, theta, phi, lambda);
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
       em[x1i][x2i][x3i] *= Planck_f(lambda, T[x1i][x2i][x3i]);
      }
  return em;
}

fp_t *** h_minus_mol::opacity(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){
  fp_t *** op = boundfree_op(Vlos, lambda);
  op = add(freefree_op(T, Ne, Vlos, lambda), op, x1l, x1h, x2l, x2h, x3l, x3h);
  return op;
}

fp_t ***** h_minus_mol::emissivity_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){
  fp_t ***** em_pert = opacity_pert(T, Ne, Vlos, Vt, B, theta, phi, lambda);
  fp_t *** op = opacity(T, Ne, Vlos, Vt, B, theta, phi, lambda);
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        fp_t Planck_local = Planck_f(lambda, T[x1i][x2i][x3i]);
        em_pert[1][x3i][x1i][x2i][x3i] *= Planck_local;
        em_pert[1][x3i][x1i][x2i][x3i] += op[x1i][x2i][x3i] * Planck_f_derivative(lambda, T[x1i][x2i][x3i]);
          
        em_pert[2][x3i][x1i][x2i][x3i] *= Planck_local;
      }

  del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
  
  return em_pert;
}

fp_t ***** h_minus_mol::opacity_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){
  fp_t ***** op_pert = boundfree_op_pert(Vlos, lambda);
  op_pert = add(freefree_op_pert(T, Ne, Vlos, lambda),op_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  
  return op_pert;
}

// -----------------------------------------------------------------------------------------------------------------------------------

// Very simple functions to compute population responses. So far, only to Temperature:

void h_minus_mol::compute_lte_population_responses(){
  // This one will be made in a really simple way:
  dN = ft4dim(1,7,x1l,x1h,x2l,x2h,x3l,x3h); 
  memset(dN[1][x1l][x2l]+x3l,0,7*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)     
      for (int x3i=x3l;x3i<=x3h;++x3i){

        // Temperature
        fp_t local_T = fetch_temperature(x1i,x2i,x3i);
        // Perturb Temperature:
        parent_atm->set_Temp(x1i,x2i,x3i,local_T+delta_T * 0.5);
        // Solve chemeq for the local point only:
        parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
        dN[1][x1i][x2i][x3i] = N[x1i][x2i][x3i];
        // Perturb temperature again:
        parent_atm->set_Temp(x1i,x2i,x3i,local_T-delta_T * 0.5);
        parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
        dN[1][x1i][x2i][x3i] -= N[x1i][x2i][x3i];
        dN[1][x1i][x2i][x3i] /= delta_T;
        parent_atm->set_Temp(x1i,x2i,x3i,local_T);
        
        // Density:
        fp_t Nt = fetch_Nt(x1i,x2i,x3i);
        fp_t Nt_step = delta_Nt_frac * Nt;
        parent_atm->set_Nt(x1i,x2i,x3i,Nt+Nt_step*0.5);
        // Solve chemeq for the local point only:
        parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
        dN[2][x1i][x2i][x3i] = N[x1i][x2i][x3i];
        // Perturb temperature again:
        parent_atm->set_Nt(x1i,x2i,x3i,Nt-Nt_step*0.5);
        parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
        dN[2][x1i][x2i][x3i] -= N[x1i][x2i][x3i];
        dN[2][x1i][x2i][x3i] /= Nt_step;
        parent_atm->set_Nt(x1i,x2i,x3i,Nt);
        parent_atm->execute_chemeq_for_point(x1i,x2i,x3i);
      }
}

fp_t h_minus_mol::get_population_response(int parameter, int x3k, int x1i, int x2i, int x3i, int z){
    if (x3k==x3i)
      return dN[parameter][x1l][x2l][x3i];
    return 0;
}

void h_minus_mol::compute_nlte_population_responses(){
  // So far nothing better then lte!
  compute_lte_population_responses();
}

int h_minus_mol::responses_clear(){
  del_ft4dim(dN,1,7,x1l,x1h,x2l,x2h,x3l,x3h);
  return 0;
}


// -----------------------------------------------------------------------------------------------------------------------------------

fp_t **** h_minus_mol::emissivity_vector(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){

  fp_t **** em;
  fp_t *** em_scalar = emissivity(T, Ne, Vlos, Vt, B, theta, phi, lambda);
  em = ft4dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, 4);
  memset(em[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  for (int x1i = x1l; x1i <= x1h; ++x1i)
    for (int x2i = x2l; x2i <= x2h; ++x2i)
      for (int x3i = x3l; x3i <= x3h; ++x3i)
          em[x1i][x2i][x3i][1] = em_scalar[x1i][x2i][x3i];

  del_ft3dim(em_scalar,x1l,x1h,x2l,x2h,x3l,x3h);
  return em;
}

fp_t *** h_minus_mol::emissivity_vector_synth(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,
  fp_t* lambda, fp_t nlambda){

  fp_t *** em_vector = ft3dim(1,nlambda,x3l,x3h,1,4);
  fp_t lambda_m = (lambda[1] + lambda[nlambda]) * 0.5;

  fp_t **** em_vector_mean = emissivity_vector(T,Ne,Vlos,Vt,B,theta,phi,lambda_m);

  

}

fp_t ***** h_minus_mol::opacity_vector(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){

  fp_t ***** op;
  fp_t *** op_scalar = opacity(T, Ne, Vlos, Vt, B, theta, phi, lambda);
  op = ft5dim(x1l, x1h, x2l, x2h, x3l, x3h, 1, 4, 1, 4);
  memset(op[x1l][x2l][x3l][1]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));
  for (int x1i = x1l; x1i <= x1h; ++x1i)
    for (int x2i = x2l; x2i <= x2h; ++x2i)
      for (int x3i = x3l; x3i <= x3h; ++x3i)
        for (int s = 1; s<=4; ++s)
          op[x1i][x2i][x3i][s][s] = op_scalar[x1i][x2i][x3i];

  del_ft3dim(op_scalar,x1l,x1h,x2l,x2h,x3l,x3h);

  return op;
}

fp_t ****** h_minus_mol::emissivity_vector_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){

  fp_t ****** em_pert = ft6dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  memset(em_pert[1][x3l][x1l][x2l][x3l]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));

  fp_t ***** em_scalar_pert = emissivity_pert(T, Ne, Vlos, Vt, B, theta, phi, lambda);

  for (int p=1;p<=7;++p)
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i){
            em_pert[p][x3k][x1i][x2i][x3i][1] = em_scalar_pert[p][x3k][x1i][x2i][x3i];
  }

  del_ft5dim(em_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  return em_pert;
}

fp_t ******* h_minus_mol::opacity_vector_pert(fp_t ***T,fp_t ***Ne,fp_t ***Vlos,fp_t ***Vt, fp_t **** B, fp_t theta,fp_t phi,fp_t lambda){

  fp_t ******* op_pert = ft7dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  memset(op_pert[1][x3l][x1l][x2l][x3l][1]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  fp_t ***** op_scalar_pert = opacity_pert(T, Ne, Vlos, Vt, B, theta, phi, lambda);

  for (int p=1;p<=7;++p)
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i){
            op_pert[p][x3k][x1i][x2i][x3i][1][1] = op_pert[p][x3k][x1i][x2i][x3i][2][2] = op_pert[p][x3k][x1i][x2i][x3i][3][3] = 
              op_pert[p][x3k][x1i][x2i][x3i][4][4] = op_scalar_pert[p][x3k][x1i][x2i][x3i];
  }

  del_ft5dim(op_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  
  return op_pert;
}




// ----------------------------------------------------------------------------------------------------------------------------------

fp_t *** h_minus_mol::freefree_op(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t lambda){

  // 05/11/2014: Ok here we will now try to implement the tables and exponents from John(1988). 

  // First we write down polynomial coefficients for wavelengts greated then 0.3645 micrometer.

  fp_t C_L[6][6]; // The first index goes from A to F, second one from 1 to 6

  C_L[0][0] = 0.0; C_L[0][1] = 2483.3460; C_L[0][2] = -3449.889; C_L[0][3] = 2200.04; C_L[0][4] = -696.271; C_L[0][5] = 88.283;
  C_L[1][0] = 0.0; C_L[1][1] = 285.827; C_L[1][2] = -1158.382; C_L[1][3] = 2427.719; C_L[1][4] = -1841.4; C_L[1][5] = 444.517;
  C_L[2][0] = 0.0; C_L[2][1] = -2054.291; C_L[2][2] = 8746.523; C_L[2][3] = -13651.105; C_L[2][4] = 8624.97; C_L[2][5] = -1863.864;
  C_L[3][0] = 0.0; C_L[3][1] = 2827.776; C_L[3][2] = -11485.632; C_L[3][3] = 16755.524; C_L[3][4] = -10051.53; C_L[3][5] = 2095.288;
  C_L[4][0] = 0.0; C_L[4][1] = -1341.537; C_L[4][2] = 5303.609; C_L[4][3] = -7510.494; C_L[4][4] = 4400.067; C_L[4][5] = -901.788;
  C_L[5][0] = 0.0; C_L[5][1] = 208.952; C_L[5][2] = -812.939; C_L[5][3] = 1132.738; C_L[5][4] = -655.02; C_L[5][5] = 132.985;

  fp_t ***op = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  fp_t T_P = 0.0; // 5040./T 
  fp_t lam = lambda * 1E4;

  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        op[x1i][x2i][x3i] = 0.0;
        T_P = 5040. / T[x1i][x2i][x3i];
        for (int n =0; n<6; ++n){
          op[x1i][x2i][x3i] += pow(T_P, (n+2) * 0.5) * (C_L[0][n] * lam * lam + C_L[1][n] + C_L[2][n] / lam +
            C_L[3][n] / lam / lam + C_L[4][n] / lam / lam / lam + C_L[5][n] / lam / lam / lam / lam);
        }
        op[x1i][x2i][x3i] *= N[x1i][x2i][x3i] * 1.33E-29 * pow(T[x1i][x2i][x3i], 2.5) * exp(-8750.0 / T[x1i][x2i][x3i]);
      }
  return op;
}

fp_t ***** h_minus_mol::freefree_op_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t lambda){

  // 05/11/2014: Ok here we will now try to implement the tables and exponents from John(1988). 

  // First we write down polynomial coefficients for wavelengts greated then 0.3645 micrometer.

  fp_t C_L[6][6]; // The first index goes from A to F, second one from 1 to 6

  C_L[0][0] = 0.0; C_L[0][1] = 2483.3460; C_L[0][2] = -3449.889; C_L[0][3] = 2200.04; C_L[0][4] = -696.271; C_L[0][5] = 88.283;
  C_L[1][0] = 0.0; C_L[1][1] = 285.827; C_L[1][2] = -1158.382; C_L[1][3] = 2427.719; C_L[1][4] = -1841.4; C_L[1][5] = 444.517;
  C_L[2][0] = 0.0; C_L[2][1] = -2054.291; C_L[2][2] = 8746.523; C_L[2][3] = -13651.105; C_L[2][4] = 8624.97; C_L[2][5] = -1863.864;
  C_L[3][0] = 0.0; C_L[3][1] = 2827.776; C_L[3][2] = -11485.632; C_L[3][3] = 16755.524; C_L[3][4] = -10051.53; C_L[3][5] = 2095.288;
  C_L[4][0] = 0.0; C_L[4][1] = -1341.537; C_L[4][2] = 5303.609; C_L[4][3] = -7510.494; C_L[4][4] = 4400.067; C_L[4][5] = -901.788;
  C_L[5][0] = 0.0; C_L[5][1] = 208.952; C_L[5][2] = -812.939; C_L[5][3] = 1132.738; C_L[5][4] = -655.02; C_L[5][5] = 132.985;

  fp_t *****op_pert = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  fp_t T_P = 0.0; // 5040./T 
  fp_t lam = lambda * 1E4;

  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        
        fp_t sigma = 0.0;
        fp_t d_sigma_dT = 0.0;

        T_P = 5040. / T[x1i][x2i][x3i];
        for (int n =0; n<6; ++n)
          d_sigma_dT += ((n+2.0) * 0.5)* (-T_P) / T[x1i][x2i][x3i] * pow(T_P, n * 0.5) * (C_L[0][n] * lam * lam + C_L[1][n] + C_L[2][n] / lam +
            C_L[3][n] / lam / lam + C_L[4][n] / lam / lam / lam + C_L[5][n] / lam / lam / lam / lam);

        
        for (int n =0; n<6; ++n)
          sigma += pow(T_P, (n+2.0) * 0.5) * (C_L[0][n] * lam * lam + C_L[1][n] + C_L[2][n] / lam +
            C_L[3][n] / lam / lam + C_L[4][n] / lam / lam / lam + C_L[5][n] / lam / lam / lam / lam);
       
        op_pert[1][x3i][x1i][x2i][x3i] = sigma * dN[1][x1i][x2i][x3i] * 1.33E-29 * pow(T[x1i][x2i][x3i], 2.5) * exp(-8750.0 / T[x1i][x2i][x3i]);
        op_pert[1][x3i][x1i][x2i][x3i] += N[x1i][x2i][x3i] * 1.33E-29 * pow(T[x1i][x2i][x3i], 2.5) * exp(-8750.0 / T[x1i][x2i][x3i]) * d_sigma_dT;
        op_pert[1][x3i][x1i][x2i][x3i] += sigma * N[x1i][x2i][x3i] * 1.33E-29 * 2.5 * pow(T[x1i][x2i][x3i], 1.5) * exp(-8750.0 / T[x1i][x2i][x3i]);
        op_pert[1][x3i][x1i][x2i][x3i] += sigma * N[x1i][x2i][x3i] * 1.33E-29 * pow(T[x1i][x2i][x3i], 2.5) * exp(-8750.0 / T[x1i][x2i][x3i]) 
          * 8750.0 / T[x1i][x2i][x3i] / T[x1i][x2i][x3i];

        op_pert[2][x3i][x1i][x2i][x3i] = sigma * dN[2][x1i][x2i][x3i] * 1.33E-29 * pow(T[x1i][x2i][x3i], 2.5) * exp(-8750.0 / T[x1i][x2i][x3i]);
      }

  return op_pert;
}

fp_t *** h_minus_mol::freefree_em(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t lambda){

  // FIX!
  return 0;
}

// 31/10/2014, Halloween: We are implementing b-f opacity and emissivity from a H- ion:

fp_t *** h_minus_mol::boundfree_op(fp_t*** Vlos, fp_t lambda){
  
  // Allocate space for the opacity:
  fp_t ***op=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  // For each point in the medium compute the opacity:  
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
          //op[x1i][x2i][x3i] = 3.95E-17 * exp(- (lambda - 8.5E-5) * (lambda - 8.5E-5) /6.02 / 6.02 * 1E6) * N[x1i][x2i][x3i];
        fp_t h_nu = 19.8581e-17 / lambda;
        fp_t arg_temp = ((h_nu - 1.20823E-12) / h_nu / h_nu);
         op[x1i][x2i][x3i] = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * N[x1i][x2i][x3i];
      }
//
  return op;
}

fp_t ***** h_minus_mol::boundfree_op_pert(fp_t*** Vlos, fp_t lambda){
  
  // Allocate space for the opacity:
  fp_t *****op_pert=ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  // For each point in the medium compute the opacity:  
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
          //op[x1i][x2i][x3i] = 3.95E-17 * exp(- (lambda - 8.5E-5) * (lambda - 8.5E-5) /6.02 / 6.02 * 1E6) * N[x1i][x2i][x3i];
        fp_t h_nu = 19.8581e-17 / lambda;
        fp_t arg_temp = ((h_nu - 1.20823E-12) / h_nu / h_nu);
        op_pert[1][x3i][x1i][x2i][x3i] = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * dN[1][x1i][x2i][x3i];
        op_pert[2][x3i][x1i][x2i][x3i] = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * dN[2][x1i][x2i][x3i];
      }
//
  return op_pert;
}


fp_t h_minus_mol::opacity_continuum(fp_t T, fp_t Ne, fp_t lambda, int x1i, int x2i, int x3i){

  // Basically we just add the contributions for ff and b-f opacity together in the body of the function. 
  // There has to be more elegant way to work on this. But ok, let it be.
  fp_t op_bf=0;
  // Add b-f
  fp_t h_nu = 19.8581e-17 / lambda;
  fp_t arg_temp = ((h_nu - 1.20823E-12) / h_nu / h_nu);
  op_bf = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * N[x1i][x2i][x3i];
  // And then f-f

  fp_t C_L[6][6]; // The first index goes from A to F, second one from 1 to 6
  C_L[0][0] = 0.0; C_L[0][1] = 2483.3460; C_L[0][2] = -3449.889; C_L[0][3] = 2200.04; C_L[0][4] = -696.271; C_L[0][5] = 88.283;
  C_L[1][0] = 0.0; C_L[1][1] = 285.827; C_L[1][2] = -1158.382; C_L[1][3] = 2427.719; C_L[1][4] = -1841.4; C_L[1][5] = 444.517;
  C_L[2][0] = 0.0; C_L[2][1] = -2054.291; C_L[2][2] = 8746.523; C_L[2][3] = -13651.105; C_L[2][4] = 8624.97; C_L[2][5] = -1863.864;
  C_L[3][0] = 0.0; C_L[3][1] = 2827.776; C_L[3][2] = -11485.632; C_L[3][3] = 16755.524; C_L[3][4] = -10051.53; C_L[3][5] = 2095.288;
  C_L[4][0] = 0.0; C_L[4][1] = -1341.537; C_L[4][2] = 5303.609; C_L[4][3] = -7510.494; C_L[4][4] = 4400.067; C_L[4][5] = -901.788;
  C_L[5][0] = 0.0; C_L[5][1] = 208.952; C_L[5][2] = -812.939; C_L[5][3] = 1132.738; C_L[5][4] = -655.02; C_L[5][5] = 132.985;

  fp_t lam = lambda * 1E4;
  fp_t T_P = 5040. / T;
  fp_t op_ff = 0.0;
  for (int n =0; n<6; ++n)
    op_ff += pow(T_P, (n+2) * 0.5) * (C_L[0][n] * lam * lam + C_L[1][n] + C_L[2][n] / lam +
      C_L[3][n] / lam / lam + C_L[4][n] / lam / lam / lam + C_L[5][n] / lam / lam / lam / lam);
  
  op_ff *= N[x1i][x2i][x3i] * 1.33E-29 * pow(T, 2.5) * exp(-8750.0 / T);

  return op_bf + op_ff;

}

fp_t** h_minus_mol::opacity_continuum_pert(fp_t T, fp_t Ne, fp_t lambda, int x1i, int x2i, int x3i){

  fp_t ** op_continuum_pert = ft2dim(1,7,x3l,x3h);
  memset(op_continuum_pert[1]+x3l,0,7*(x3h-x3l+1)*sizeof(fp_t));
  fp_t op_cont = opacity_continuum(T, Ne, lambda, x1i, x2i, x3i);
  op_continuum_pert[1][x3i] = op_cont * dN[1][x1i][x2i][x3i] / N[x1i][x2i][x3i];
  op_continuum_pert[2][x3i] = op_cont * dN[2][x1i][x2i][x3i] / N[x1i][x2i][x3i];
  return op_continuum_pert;
}

fp_t *** h_minus_mol::boundfree_em(fp_t*** Vlos, fp_t lambda){
  return 0;
}
