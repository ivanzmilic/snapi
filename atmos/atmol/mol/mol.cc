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
}


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

void h_minus_mol::compute_lte_population_responses_analytical(fp_t ***, fp_t ***){
  compute_lte_population_responses();
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

int h_minus_mol::op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B, fp_t theta,fp_t phi,
  fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  // This is function which computes opacity and emissivity at wavelengths simultaneously assuming H- opacity is flat. 
  // (There is no real jump as we assume that F-F opacity is sort of off-seting that)

  fp_t lambda_m = (lambda[nlambda] + lambda[1]) * 0.5;

  fp_t *** op_scalar = opacity(T,Ne,Vlos,Vt,B,theta,phi,lambda_m);
  // For emissivity we assume it is in LTE with opacity
  fp_t *** em_scalar = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  for (int x1i = x1l; x1i <= x1h; ++x1i)
    for (int x2i = x2l; x2i <= x2h; ++x2i)
      for (int x3i = x3l; x3i <= x3h; ++x3i)
        em_scalar[x1i][x2i][x3i] = op_scalar[x1i][x2i][x3i] * Planck_f(lambda_m,T[x1i][x2i][x3i]);
  

  // Now we just copy scalar, monochromatic opacity to our vector:
  for (int l=1;l<=nlambda;++l)
    for (int x1i = x1l; x1i <= x1h; ++x1i)
      for (int x2i = x2l; x2i <= x2h; ++x2i)
        for (int x3i = x3l; x3i <= x3h; ++x3i){
          op_vector[l][x1i][x2i][x3i][1][1] += op_scalar[x1i][x2i][x3i];
          em_vector[l][x1i][x2i][x3i][1]    += em_scalar[x1i][x2i][x3i];
  }
  del_ft3dim(op_scalar,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(em_scalar,x1l,x1h,x2l,x2h,x3l,x3h);

  return 0;
}

int h_minus_mol::op_em_vector_plus_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B, fp_t theta,fp_t phi,
  fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector,fp_t ********op_vector_pert, fp_t *******em_vector_pert){

  // This is function which computes opacity and emissivity at wavelengths simultaneously assuming H- opacity is flat. 
  // (There is no real jump as we assume that F-F opacity is sort of off-seting that)

  fp_t lambda_m = (lambda[nlambda] + lambda[1]) * 0.5;

  fp_t *** op_scalar = opacity(T,Ne,Vlos,Vt,B,theta,phi,lambda_m);
  fp_t ***** op_scalar_pert = opacity_pert(T,Ne,Vlos,Vt,B,theta,phi,lambda_m);
  // For emissivity we assume it is in LTE with opacity
  fp_t *** em_scalar = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  fp_t ***** em_scalar_pert = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(em_scalar_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1));
  
  for (int x1i = x1l; x1i <= x1h; ++x1i)
    for (int x2i = x2l; x2i <= x2h; ++x2i)
      for (int x3i = x3l; x3i <= x3h; ++x3i)
        em_scalar[x1i][x2i][x3i] = op_scalar[x1i][x2i][x3i] * Planck_f(lambda_m,T[x1i][x2i][x3i]);

  for (int p=1;p<=2;++p)
    for (int x1i = x1l; x1i <= x1h; ++x1i)
      for (int x2i = x2l; x2i <= x2h; ++x2i)
        for (int x3i = x3l; x3i <= x3h; ++x3i){
          em_scalar_pert[p][x3i][x1i][x2i][x3i] = em_scalar[x1i][x2i][x3i]/op_scalar[x1i][x2i][x3i]
            * op_scalar_pert[p][x3i][x1i][x2i][x3i];
          if (p==1)
            em_scalar_pert[p][x3i][x1i][x2i][x3i] += op_scalar[x1i][x2i][x3i] * Planck_f_derivative(lambda_m,T[x1i][x2i][x3i]);
        }
    
  // Now we just copy scalar, monochromatic opacity to our vector:
  for (int l=1;l<=nlambda;++l)
    for (int x1i = x1l; x1i <= x1h; ++x1i)
      for (int x2i = x2l; x2i <= x2h; ++x2i)
        for (int x3i = x3l; x3i <= x3h; ++x3i){
          op_vector[l][x1i][x2i][x3i][1][1] += op_scalar[x1i][x2i][x3i];
          em_vector[l][x1i][x2i][x3i][1]    += em_scalar[x1i][x2i][x3i];
  }

  // And the same for responses:
  for (int p=1;p<=2;++p)
    for (int l=1;l<=nlambda;++l)
      for (int x1i = x1l; x1i <= x1h; ++x1i)
        for (int x2i = x2l; x2i <= x2h; ++x2i)
          for (int x3i = x3l; x3i <= x3h; ++x3i){
            op_vector_pert[l][p][x3i][x1i][x2i][x3i][1][1] += op_scalar_pert[p][x3i][x1i][x2i][x3i];
            em_vector_pert[l][p][x3i][x1i][x2i][x3i][1]    += em_scalar_pert[p][x3i][x1i][x2i][x3i];
  }
  del_ft3dim(op_scalar,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft3dim(em_scalar,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft5dim(op_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  del_ft5dim(em_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);

  return 0;
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
  // We gave this up on 25/06/2017 to implement the method from Gray:

  
  fp_t ***op = ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  
  // Now approach from Grey p. 155-156:
  
  fp_t l = log10(lambda*1E8); // lambda in Angstroms:
  fp_t f_0 = -2.2763 - 1.6850*l + 0.76661*l*l - 0.053346*l*l*l;
  fp_t f_1 = 15.2827 - 9.2846*l + 1.99381*l*l - 0.142631*l*l*l;
  fp_t f_2 = -197.789 + 190.266*l - 67.9775*l*l + 10.6913*l*l*l - 0.625151*l*l*l*l;

  fp_t T_crit = h*c/1.6421E-4/k;
  
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        fp_t T_P = 5040. / T[x1i][x2i][x3i];
        op[x1i][x2i][x3i] = 1E-26 * Ne[x1i][x2i][x3i]*k*T[x1i][x2i][x3i] * 
          pow(10.0,f_0+f_1*log10(T_P)+f_2*log10(T_P)*log10(T_P))*fetch_population(x1i,x2i,x3i,0,0);
  } // points in the atmosphere

  return op;
}

// ======================================================================================================

fp_t ***** h_minus_mol::freefree_op_pert(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t lambda){

  // Will be the same approach as above (Gray book):
  fp_t *****op_pert = ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  // These do not change with the temperature, this is only lambda dependence:
  fp_t l = log10(lambda*1E8); // lambda in Angstrom
  fp_t f_0 = -2.2763 - 1.6850*l + 0.76661*l*l - 0.053346*l*l*l;
  fp_t f_1 = 15.2827 - 9.2846*l + 1.99381*l*l - 0.142631*l*l*l;
  fp_t f_2 = -197.789 + 190.266*l - 67.9775*l*l + 10.6913*l*l*l - 0.625151*l*l*l*l;

  // Now the temperature/density dependence comes:
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        // Do some step-by-step analytical population of derivatives:
        // Temperature: 
        fp_t T_P = 5040. / T[x1i][x2i][x3i];
        fp_t d_T_P_dT = -T_P/T[x1i][x2i][x3i]; // derivative w.r.t T
        fp_t d_logTP_dT = 1.0/2.302585/T_P * d_T_P_dT;
        fp_t exp_factor = pow(10.0,f_0+f_1*log10(T_P)+f_2*log10(T_P)*log10(T_P));
        fp_t d_exp_factor_dT = exp_factor*2.302585*(f_1 * d_logTP_dT + 2.0*f_2*log10(T_P)*d_logTP_dT);

        fp_t N_H = fetch_population(x1i,x2i,x3i,0,0); // Get population of neutral hydrogen
        fp_t N_e = Ne[x1i][x2i][x3i]; // Electrons
        fp_t const_factor = 1E-26*k; // This is the part which does not depend on the temperature.

        // From multiplication by T:
        op_pert[1][x3i][x1l][x2l][x3i] = const_factor * N_e * exp_factor * N_H;
        // From the exp_factor derivative:
        op_pert[1][x3i][x1l][x2l][x3i] += const_factor * T[x1i][x2i][x3i] * N_e * N_H * d_exp_factor_dT;
        // From the derivatives of N_e and N_H
        fp_t dN_e_dT = parent_atm->get_ne_lte_derivative(1,x1i,x2i,x3i);
        fp_t dN_H_dT = parent_atm->get_neutral_H_derivative_lte(1,x1i,x2i,x3i);
        op_pert[1][x3i][x1l][x2l][x3i] += const_factor * T[x1i][x2i][x3i] * exp_factor * 
          (N_e * dN_H_dT + N_H * dN_e_dT);

        // Density:
        // Only derivatives of N_e and N_H
        fp_t dN_e_dN = parent_atm->get_ne_lte_derivative(2,x1i,x2i,x3i);
        fp_t dN_H_dN = parent_atm->get_neutral_H_derivative_lte(2,x1i,x2i,x3i);
        op_pert[2][x3i][x1l][x2l][x3i] = const_factor * T[x1i][x2i][x3i] * exp_factor * 
          (N_e * dN_H_dT + N_H * dN_e_dT);
  }
  return op_pert;
}

// ======================================================================================================

fp_t *** h_minus_mol::freefree_em(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t lambda){

  fprintf(stderr,"h_minus_mol::freefree_em : You are calling an unfinished function. Please check.\n");
  return 0;
}

// ======================================================================================================

fp_t *** h_minus_mol::boundfree_op(fp_t*** Vlos, fp_t lambda){
  
  // Allocate space for the opacity:
  fp_t ***op=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  //
  // We are computing the opacity following the book of Gray:
  fp_t a[7]={1.99654,-1.18267E-5,2.64243E-6,-4.40524E-10,3.23992E-14,-1.39568E-18,2.78701E-23};
  fp_t alpha = 0.0;
  for (int i=0;i<7;++i)
    alpha += a[i] * pow(lambda*1E8,i);
  
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        fp_t T = fetch_temperature(x1i,x2i,x3i);
        fp_t T_P = 5040./T;
        op[x1i][x2i][x3i] = 4.158E-28 * alpha * fetch_Ne(x1i,x2i,x3i)*k*T *
          pow(T_P,2.5) * pow(10.0,0.754*T_P)*fetch_population(x1i,x2i,x3i,0,0)* (1.0 - exp(-h*c/lambda/k/T));

        //fp_t h_nu = 19.8581e-17 / lambda;
        //fp_t arg_temp = ((h_nu - 1.20823E-12) / h_nu / h_nu);
          //op[x1i][x2i][x3i] = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * N[x1i][x2i][x3i];
  } // points in the atmosphere
  return op;
}

fp_t ***** h_minus_mol::boundfree_op_pert(fp_t*** Vlos, fp_t lambda){
  
  // Allocate space for the opacity:
  fp_t *****op_pert=ft5dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

  fp_t a[7]={1.99654,-1.18267E-5,2.64243E-6,-4.40524E-10,3.23992E-14,-1.39568E-18,2.78701E-23};
  fp_t alpha = 0.0;
  for (int i=0;i<7;++i)
    alpha += a[i] * pow(lambda*1E8,i);
  
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        
        fp_t T = fetch_temperature(x1i,x2i,x3i);
        fp_t T_P = 5040.0/T;

        fp_t op = 4.158E-28 * alpha * fetch_Ne(x1i,x2i,x3i)*k*T *
          pow(T_P,2.5) * pow(10.0,0.754*T_P)*fetch_population(x1i,x2i,x3i,0,0)* (1.0 - exp(-h*c/lambda/k/T));
        
        op_pert[1][x3i][x1i][x2i][x3i]  = op * (1.0/T + 2.5 * 1.0/T_P * (-5040.0/T/T)); // direct T dependence
        op_pert[1][x3i][x1i][x2i][x3i] += op * parent_atm->get_ne_lte_derivative(1,x1i,x2i,x3i)/fetch_Ne(x1i,x2i,x3i);
        op_pert[1][x3i][x1i][x2i][x3i] += op * parent_atm->get_neutral_H_derivative_lte(1,x1i,x2i,x3i)/fetch_population(x1i,x2i,x3i,0,0);
        op_pert[1][x3i][x1i][x2i][x3i] += op * (-2.302585*0.754*5040.0/T/T);
        op_pert[1][x3i][x1i][x2i][x3i] += op * (- exp(-h*c/lambda/k/T) * h*c/lambda/k/T/T)/ (1.0 - exp(-h*c/lambda/k/T));

        op_pert[2][x3i][x1i][x2i][x3i]  = op * parent_atm->get_ne_lte_derivative(2,x1i,x2i,x3i)/fetch_Ne(x1i,x2i,x3i);
        op_pert[2][x3i][x1i][x2i][x3i] += op * parent_atm->get_neutral_H_derivative_lte(2,x1i,x2i,x3i)/fetch_population(x1i,x2i,x3i,0,0);

        //fp_t h_nu = 19.8581e-17 / lambda;
        //fp_t arg_temp = ((h_nu - 1.20823E-12) / h_nu / h_nu);
        //op_pert[1][x3i][x1i][x2i][x3i] = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * dN[1][x1i][x2i][x3i];
        //op_pert[2][x3i][x1i][x2i][x3i] = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * dN[2][x1i][x2i][x3i];
      }
//
  return op_pert;
}


fp_t h_minus_mol::opacity_continuum(fp_t T, fp_t Ne, fp_t lambda, int x1i, int x2i, int x3i){

  // Basically we just add the contributions for ff and b-f opacity together in the body of the function. 
  // There has to be more elegant way to work on this. But ok, let it be.
  
  // Grey:
  fp_t a[7]={1.99654,-1.18267E-5,2.64243E-6,-4.40524E-10,3.23992E-14,-1.39568E-18,2.78701E-23};
  fp_t alpha = 0.0;
  for (int i=0;i<7;++i)
    alpha += a[i] * pow(lambda*1E8,i);  
  fp_t T_P = 5040./T;
  fp_t op_bf = 4.158E-28 * alpha * Ne*k*T * pow(T_P,2.5) * pow(10.0,0.754*T_P) * fetch_population(x1i,x2i,x3i,0,0)
    *(1.0 - exp(-8764.2/T));

  fp_t h_nu = 19.8581e-17 / lambda;
  fp_t arg_temp = ((h_nu - 1.20823E-12) / h_nu / h_nu);
  op_bf = 4.19648E-34 * sqrt(arg_temp * arg_temp * arg_temp) * N[x1i][x2i][x3i];


  // Grey:
  fp_t l = log10(lambda*1E8);
  fp_t f_0 = -2.2763 - 1.6850*l + 0.76661*l*l - 0.053346*l*l*l;
  fp_t f_1 = 15.2827 - 9.2846*l + 1.99381*l*l - 0.142631*l*l*l;
  fp_t f_2 = -197.789 + 190.266*l - 67.9775*l*l + 10.6913*l*l*l - 0.625151*l*l*l*l;

  T_P = 5040. / T;
  fp_t op_ff = 1E-26 * Ne*k*T * pow(10.0,f_0+f_1*log10(T_P)+f_2*log10(T_P)*log10(T_P))*
    fetch_population(x1i,x2i,x3i,0,0);

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
