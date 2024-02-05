#include <string.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "pack.h"
#include "struts.h"
#include "const.h"
#include "math.h"

//#include "obs/obs.h"
#include "atmol/atmol.h"

#include "acfg.h"

#include "plpar/atmos_pp.h"
#include "atmos.h"

atmosphere *atmos_new(acfg *cfg,io_class &io_in)
{
  if(!strcmp(cfg->geom,"PLANEPARALLEL")) return atmos_pp_new(cfg,io_in);
  if(!strcmp(cfg->geom,"CARTESIAN3D")) return new atmosphere(cfg,io_in);
  return new atmosphere(cfg,io_in);
}

atmosphere *atmos_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in)
{
  uint08_t gtype;
  offs+=::unpack(buf+offs,gtype);
  int08_t bla;
  ::unpack(buf+offs,bla);
  switch(gtype){
    case(ATMOS_GEOM_PP): return atmos_pp_new(buf,offs,do_swap,io_in);
//    case(ATMOS_GEOM_2D): return new atmos_2d(buf,offs,do_swap,io_in);
    case(ATMOS_GEOM_3D):{ // this should call atmos_3d_new
      int08_t rtstype;
      offs+=::unpack(buf,rtstype);
      return new atmosphere(buf,offs,do_swap,io_in); // this is also the default
    }
  }
  int08_t rtstype;
  offs+=::unpack(buf,rtstype);
  return new atmosphere(buf,offs,do_swap,io_in);
}

atmosphere::atmosphere(acfg *cfg,io_class &io_in):grid(io_in),flags(ATMOS_FLAG_MASK)
{
// init
  flags.clear();
// some output
  io.msg(IOL_INFO,"atmosphere::atmosphere: ID=%s\n",cfg->id);
  id=dup_str(cfg->id);
  io.msg(IOL_INFO,"atmosphere::atmosphere: FILE=%s\n",cfg->filename);
  fname=dup_str(cfg->filename);
//
  io.msg(IOL_INFO,"atmosphere::atmosphere: TYPE=%s\n",cfg->filetype);
  struct typemap typemap[]=FILETYPES;
  for(int i=0;typemap[i].type_name;++i) if((ftype=typemap[i].equals(cfg->filetype))) break;
  if(!ftype) io_in.msg(IOL_ERROR,"atmosphere::atmosphere: unknown file type \"%s\"\n",cfg->filetype);
// check type and dimensions, read the whole thing?
  grid::resize(1,cfg->nx,1,cfg->ny,1,cfg->nz,
               0.0,(fp_t)(cfg->nx-1)*cfg->dx,
               0.0,(fp_t)(cfg->ny-1)*cfg->dy,
               0.0,(fp_t)(cfg->nz-1)*cfg->dz);
  io.msg(IOL_INFO,"atmosphere::atmosphere: dims=%d..%d, %d..%d, %d..%d\n",x1l,x1h,x2l,x2h,x3l,x3h);
  io.msg(IOL_INFO,"atmosphere::atmosphere: grid=%E x %E x %E\n",cfg->dx,cfg->dy,cfg->dz);
// initialize variables
  fp_t ****p[]={&T,&rho,&Nt,&Ne,&Bx,&By,&Bz,&Vx,&Vy,&Vz,&Vt,&tau_referent,&op_referent,0};
  if (x1h-x1l+1 && x2h-x2l+1 && x3h-x3l+1)
    for(int i=0;p[i];++i){
      (*(p[i]))=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
      memset((*(p[i]))[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  }
// setup atomic and molecular data
  natm=cfg->natm+cfg->nmol;
  boundary_condition_for_rt = -1; // Semi-infinite atmosphere
  io.msg(IOL_INFO,"atmosphere::atmosphere: NATOM=%d\n",cfg->natm);
  atml=new atmol* [natm+1];
  for(int a=0;a<cfg->natm;++a){ 
    atml[a]=atmol_new(cfg->atm[a],io_in);
  }
  io.msg(IOL_INFO,"atmosphere::atmosphere: NMOL=%d\n",cfg->nmol);
  for(int m=0;m<cfg->nmol;++m){
    atml[m+cfg->natm]=atmol_new(cfg->mol[m],atml,cfg->natm,io_in); // This makes a new object of type atmol which is either atom or mol (or H-)
  }
  if (cfg->of){ //if fudge is on
    if(FILE *f=fopen(cfg->of_filename,"r")){
      char line[1000];
      if(!fgets(line,1000,f)) io_in.msg(IOL_WARN,"atmosphere::atmosphere: when reading line 1\n");
      sscanf(line,"%d",&N_of);
      io_in.msg(IOL_INFO,"atmosphere::atmosphere: %d opacity fudge points\n",N_of);

      if (N_of){
        lambda_of = new fp_t [N_of]-1;
        value_of = new fp_t [N_of]-1;
      }

      for(int i=1;i<=N_of;++i){
        if(!fgets(line,1000,f)) io_in.msg(IOL_WARN,"atmosphere::atmosphere: when reading line %d %s \n",i+1,cfg->of_filename);
        if(feof(f)){
          io_in.msg(IOL_ERROR,"atmosphere::atmosphere: unexpected end of file reading opacity fudge file line %d\n",i+1);
          fclose(f);
        }
        float32_t lambda_in, value_in;
        sscanf(line,"%E %E",&lambda_in,&value_in);
        lambda_of[i] = lambda_in;
        value_of[i] = value_in;
      }
      fclose(f);
    }
    else 
      io_in.msg(IOL_ERROR,"atmosphere::atmosphere: error opening file: %s\n",cfg->of_filename);
  }
  else 
    N_of=0;
  conserve_charge = cfg->conserve_charge;
  tau_grid = cfg->tau_grid;
  use_atm_lvls = cfg->use_atm_lvls;
  n_lvls = 0; // We are going to calculate the amount of atomic levels either way
  for (int a=0; a<natm; ++a)
      n_lvls += atml[a]->get_total_lvls();
  n_lvls += 1; // This is for electrons 

  atm_lvl_pops = 0; // Still not allocated
//
  flags.set(ATMOS_FLAG_DEF);
//
  gtype=ATMOS_GEOM_3D;
  rtstype=ATMOS_RTS_QSC;
}

atmosphere::atmosphere(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):grid(buf,offs,do_swap,io_in),flags(ATMOS_FLAG_MASK)
{
  gtype=ATMOS_GEOM_3D;
  rtstype=ATMOS_RTS_QSC;
  offs+=unpack(buf+offs,do_swap,io_in);
  
}

atmosphere::~atmosphere(void)
{
  if(atml){
    for(int a=0;a<natm;++a) delete atml[a];
    delete[] atml;
  }
//
  fp_t ****p[]={&T,&rho,&Nt,&Ne,&Bx,&By,&Bz,&Vx,&Vy,&Vz,&Vt,&tau_referent,&op_referent,0};
  if((x1l<=x1h)&&(x2l<=x2h)&&(x3l<=x3h))
    for(int i=0;p[i];++i) if(*(p[i])) del_ft3dim(*(p[i]),x1l,x1h,x2l,x2h,x3l,x3h);
//
  if(id) delete[] id;
  if(fname) delete[] fname;

  if(N_of){
    delete[](lambda_of+1);
    delete[](value_of+1);
  }
  if (use_atm_lvls == 2 && n_lvls)
    del_ft4dim(atm_lvl_pops,x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls);
}

void atmosphere::set_grid(int input){
  // If you gave positive number, set it to tau. Otherwise, leave it intact
  if (input){
    tau_grid = 1;
    rt_grid = tau_referent[x1l][x2l];
  }
  else {
    tau_grid = 0;
    rt_grid = x3;
  }
}

int08_t atmosphere::resize(int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
{
  fp_t ****p[]={&T,&rho,&Nt,&Ne,&Bx,&By,&Bz,&Vx,&Vy,&Vz,&Vt,&tau_referent,&op_referent,0};
  for(int i=0;p[i];++i) if(*(p[i])) del_ft3dim(*(p[i]),x1l,x1h,x2l,x2h,x3l,x3h); // ANOTHER PROBLEM HERE, INVESTIGATE THIS 
  //
  grid::resize(x1l_in,x1h_in,x2l_in,x2h_in,x3l_in,x3h_in);
//
  for(int i=0;p[i];++i){
    (*(p[i]))=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h);
    memset((*(p[i]))[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  }
  return 0;
}

int32_t atmosphere::size(io_class &io_in)
{
  int32_t sz=sizeof(gtype);
  sz+=sizeof(rtstype);
//
  sz+=grid::size();
// add local stuff
  sz+=flags.size();
//
  sz+=strlen(id)+1;
  sz+=strlen(fname)+1;
  sz+=sizeof(ftype);
//
  sz+=sizeof(natm);
  sz+=sizeof(boundary_condition_for_rt);
  for(int a=0;a<natm;++a) sz+=atml[a]->size(io_in);
//
  fp_t ****p[]={&T,&rho,&Nt,&Ne,&Bx,&By,&Bz,&Vx,&Vy,&Vz,&Vt,&tau_referent,&op_referent,0};
  for(int i=0;p[i];++i) sz+=(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t);

  sz+=sizeof(int);// N_of (opacity_fudge)
  sz+=2*N_of*sizeof(fp_t);
  sz+=sizeof(int);// whether to conserve charge
  sz+=sizeof(int);// whether to use tau or h as the grid
  sz+=2*sizeof(int);// whether to use atomic level populations, and how many levels there are
  
  if (use_atm_lvls == 2)
    sz+=(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*n_lvls*sizeof(fp_t);
  
  return sz;
}

int32_t atmosphere::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
// selection variables
  int32_t offs=::pack(buf,gtype);
  offs+=::pack(buf+offs,rtstype);
// grid...
  offs+=grid::pack(buf+offs,do_swap,io_in);
// local stuff
  offs+=flags.pack(buf+offs,io_in);
//
  offs+=::pack(buf+offs,id);
  offs+=::pack(buf+offs,fname);
  offs+=::pack(buf+offs,ftype);
//
  offs+=::pack(buf+offs,natm,do_swap);
  offs+=::pack(buf+offs,boundary_condition_for_rt,do_swap);
  for(int a=0;a<natm;++a) offs+=atml[a]->pack(buf+offs,do_swap,io_in);
//
  fp_t ****p[]={&T,&rho,&Nt,&Ne,&Bx,&By,&Bz,&Vx,&Vy,&Vz,&Vt,&tau_referent,&op_referent,0};
  if((x1l<=x1h)&&(x2l<=x2h)&&(x3l<=x3h))
    for(int i=0;p[i];++i) offs+=::pack(buf+offs,*(p[i]),x1l,x1h,x2l,x2h,x3l,x3h,do_swap);

  offs+=::pack(buf+offs,N_of,do_swap);
  if (N_of){
    offs+=::pack(buf+offs,lambda_of,1,N_of,do_swap);
    offs+=::pack(buf+offs,value_of,1,N_of,do_swap);
  } 
  offs+=::pack(buf+offs,conserve_charge,do_swap);
  offs+=::pack(buf+offs,tau_grid,do_swap);
  offs+=::pack(buf+offs,use_atm_lvls,do_swap);
  offs+=::pack(buf+offs,n_lvls,do_swap);

  if (use_atm_lvls == 2 && n_lvls) // Only try to pack if it's actually allocated, hence == 2
    offs+=::pack(buf+offs,atm_lvl_pops,x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls, do_swap);

  //
  return offs;
}

atmosphere * atmosphere::extract(int i, int j,io_class &io_in){

  // Extracts the given column, packs it into a buffer and then extracts it a fresh 1d atmosphere
  // This seemed genius a while back, is it still?
  
  atmosphere * column;
  uint08_t do_swap = 0;
  grid * column_grid = grid::extract_grid(i,j,io_in);
  
  int32_t sz=sizeof(gtype);
  sz+=sizeof(rtstype);
//
  sz+=column_grid->size();
// add local stuff
  sz+=flags.size();
//
  sz+=strlen(id)+1;
  sz+=strlen(fname)+1;
  sz+=sizeof(ftype);
//
  sz+=sizeof(natm);
  sz+=sizeof(boundary_condition_for_rt);
  for(int a=0;a<natm;++a) sz+=atml[a]->size(io_in);

  sz+=sizeof(N_of);
  sz+=2.0*N_of*sizeof(fp_t);
  sz+=sizeof(conserve_charge);
  sz+=sizeof(tau_grid);
  sz+=sizeof(use_atm_lvls);
  sz+=sizeof(n_lvls);

  if (use_atm_lvls && n_lvls){
    sz+=(x3h-x3l+1)*n_lvls*sizeof(fp_t);
  }
//
  fp_t ****p[]={&T,&rho,&Nt,&Ne,&Bx,&By,&Bz,&Vx,&Vy,&Vz,&Vt,&tau_referent,&op_referent,0};
  for(int ii=0;p[ii];++ii) sz+=(x3h-x3l+1)*sizeof(fp_t);
  
  uint08_t *buf=new uint08_t [sz];

  // We start by packing everything that we want to be in the small atmosphere:
  int32_t offs=::pack(buf,gtype); // gtype is the same
  offs+=::pack(buf+offs,rtstype); // rtype is the same 
// grid...
  offs+=column_grid->pack(buf+offs,do_swap,io_in); 
// local stuff
  offs+=flags.pack(buf+offs,io_in); // flags are the same 
//
  offs+=::pack(buf+offs,id); // id is the same
  offs+=::pack(buf+offs,fname); // fname can be the same, we will not init it any more
  offs+=::pack(buf+offs,ftype); // ftype can also be the same
//
  offs+=::pack(buf+offs,natm,do_swap); // same
  offs+=::pack(buf+offs,boundary_condition_for_rt,do_swap); // same
  for(int a=0;a<natm;++a) offs+=atml[a]->pack(buf+offs,do_swap,io_in); //same

  // make 1x1xNZ arrays for all the parameters and pack them up:
  fp_t *** T_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(T_small[1][1]+x3l,T[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** rho_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(rho_small[1][1]+x3l,rho[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Nt_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Nt_small[1][1]+x3l,Nt[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Ne_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Ne_small[1][1]+x3l,Ne[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Bx_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Bx_small[1][1]+x3l,Bx[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** By_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(By_small[1][1]+x3l,By[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Bz_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Bz_small[1][1]+x3l,Bz[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Vx_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Vx_small[1][1]+x3l,Vx[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Vy_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Vy_small[1][1]+x3l,Vy[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Vz_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Vz_small[1][1]+x3l,Vz[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** Vt_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(Vt_small[1][1]+x3l,Vt[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** tau_referent_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(tau_referent_small[1][1]+x3l,tau_referent[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));
  fp_t *** op_referent_small = ft3dim(1,1,1,1,x3l,x3h);
  memcpy(op_referent_small[1][1]+x3l,op_referent[i][j]+x3l,(x3h-x3l+1)*sizeof(fp_t));

  
  if (use_atm_lvls == 2 && n_lvls){ // Again only extract if it's actually allocated
    fp_t **** atm_lvl_pops_small = ft4dim(1,1,1,1,x3l,x3h,1,n_lvls);
    memcpy(atm_lvl_pops_small[1][1][x3l]+1, atm_lvl_pops[i][j][x3l]+1,(x3h-x3l+1)*n_lvls*sizeof(fp_t));
    del_ft4dim(atm_lvl_pops_small,1,1,1,1,x3l,x3h,1,n_lvls);
  }
//
  fp_t **** pp[]={&T_small,&rho_small,&Nt_small,&Ne_small,&Bx_small,&By_small,&Bz_small,&Vx_small,
    &Vy_small,&Vz_small,&Vt_small,&tau_referent_small,&op_referent_small,0};
  if((x1l<=x1h)&&(x2l<=x2h)&&(x3l<=x3h))
    for(int ii=0;p[ii];++ii) offs+=::pack(buf+offs,*(pp[ii]),1,1,1,1,x3l,x3h,do_swap);

  offs+=::pack(buf+offs,N_of,do_swap);
  if (N_of){
    offs+=::pack(buf+offs,lambda_of,1,N_of,do_swap);
    offs+=::pack(buf+offs,value_of,1,N_of,do_swap);
  }
  offs+=::pack(buf+offs,conserve_charge,do_swap); 
  offs+=::pack(buf+offs,tau_grid,do_swap);
  offs+=::pack(buf+offs,use_atm_lvls,do_swap);
  offs+=::pack(buf+offs,n_lvls,do_swap);

  if (use_atm_lvls == 2 && n_lvls) // Pack this if it exists
    offs+=::pack(buf+offs,atm_lvl_pops,x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls, do_swap);

  del_ft3dim(T_small,1,1,1,1,x3l,x3h);
  del_ft3dim(rho_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Nt_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Ne_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Bx_small,1,1,1,1,x3l,x3h);
  del_ft3dim(By_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Bz_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Vx_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Vy_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Vz_small,1,1,1,1,x3l,x3h);
  del_ft3dim(Vt_small,1,1,1,1,x3l,x3h);
  del_ft3dim(tau_referent_small,1,1,1,1,x3l,x3h);
  del_ft3dim(op_referent_small,1,1,1,1,x3l,x3h);

  delete column_grid;

  offs = 0;
  column = atmos_new(buf,offs,0,io_in);
  
  return column;
}

int32_t atmosphere::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
  int32_t offs=flags.unpack(buf,io_in);
//
  offs+=::unpack(buf+offs,id);
  offs+=::unpack(buf+offs,fname);
  offs+=::unpack(buf+offs,ftype);
//
  offs+=::unpack(buf+offs,natm,do_swap);
  offs+=::unpack(buf+offs,boundary_condition_for_rt,do_swap);
  io_in.msg(IOL_DEB1,"atmosphere::atmosphere: natmol=%d\n",natm);
  atml=new atmol* [natm];
  for(int a=0;a<natm;++a) atml[a]=atmol_new(buf,offs,do_swap,atml,a,io_in);
//
  fp_t ****p[]={&T,&rho,&Nt,&Ne,&Bx,&By,&Bz,&Vx,&Vy,&Vz,&Vt,&tau_referent,&op_referent,0};
  if((x1l<=x1h)&&(x2l<=x2h)&&(x3l<=x3h))
    for(int i=0;p[i];++i) offs+=::unpack(buf+offs,(*(p[i]))=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h),x1l,x1h,x2l,x2h,x3l,x3h,do_swap);
  else
    for(int i=0;p[i];++i) (*(p[i]))=0;

  offs+=::unpack(buf+offs,N_of,do_swap);
  if (N_of){
    lambda_of = new fp_t[N_of]-1;
    value_of = new fp_t[N_of]-1;
    offs+=::unpack(buf+offs,lambda_of,1,N_of,do_swap);
    offs+=::unpack(buf+offs,value_of,1,N_of,do_swap);
  }
  offs+=::unpack(buf+offs,conserve_charge,do_swap);
  offs+=::unpack(buf+offs,tau_grid,do_swap);
  offs+=::unpack(buf+offs,use_atm_lvls,do_swap);
  offs+=::unpack(buf+offs,n_lvls,do_swap);

  if (use_atm_lvls == 2 && n_lvls) // Only unpack if it it alocated
    offs+=::unpack(buf+offs, atm_lvl_pops=ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls),x1l,x1h,x2l,x2h,x3l,x3h,1,n_lvls, do_swap);
  
  for (int a = 0; a<natm; ++a)
    atml[a]->set_parent_atmosphere(this);
  return offs;
}

fp_t * atmosphere::add(fp_t *src,fp_t *dst,int32_t nn){
  if(src){
    for(int32_t i=0;i<nn;++i) dst[i]+=src[i];
    delete[] src;
  }
  return dst;
}

  
fp_t ** atmosphere::add(fp_t **src,fp_t **dst,int32_t ll1,int32_t ul1,int32_t ll2,int32_t ul2){
  if(src){
    int32_t nn=(ul1-ll1+1)*(ul2-ll2+1);
    for(int32_t i=0;i<nn;++i) dst[ll1][ll2+i]+=src[ll1][ll2+i];
    del_ft2dim(src,ll1,ul1,ll2,ul2);
  }
  return dst;
}
  
fp_t *** atmosphere::add(fp_t ***src,fp_t ***dst,int32_t ll1,int32_t ul1,int32_t ll2,int32_t ul2,int32_t ll3,int32_t ul3){
  if(src){
    int32_t nn=(ul1-ll1+1)*(ul2-ll2+1)*(ul3-ll3+1);
    for(int32_t i=0;i<nn;++i) dst[ll1][ll2][ll3+i]+=src[ll1][ll2][ll3+i];
    del_ft3dim(src,ll1,ul1,ll2,ul2,ll3,ul3);
  }
  return dst;
}

fp_t **** atmosphere::add(fp_t **** src, fp_t **** dst, int32_t ll1, int32_t ul1, int32_t ll2, int32_t ul2, int32_t ll3, int32_t ul3, int32_t ll4, int32_t ul4){

  if(src){
    int32_t nn = (ul1 - ll1 + 1) * (ul2 - ll2 + 1) * (ul3 - ll3 + 1) * (ul4 - ll4 + 1);
    for (int32_t i=0;i<nn;++i)
      dst[ll1][ll2][ll3][ll4+i] += src[ll1][ll2][ll3][ll4+i];
    del_ft4dim(src, ll1, ul1, ll2, ul2, ll3, ul3, ll4, ul4);
  }
  return dst;
}

fp_t ***** atmosphere::add(fp_t *****src, fp_t *****dst, int32_t ll1, int32_t ul1, int32_t ll2, int32_t ul2, int32_t ll3, int32_t ul3, int32_t ll4, int32_t ul4, int32_t ll5, int32_t ul5){

  if(src){
    int32_t nn = (ul1 - ll1 + 1) * (ul2 - ll2 + 1) * (ul3 - ll3 + 1) * (ul4 - ll4 + 1) * (ul5 - ll5 + 1);
    for (int32_t i =0 ; i<nn; ++i)
      dst[ll1][ll2][ll3][ll4][ll5+i] += src[ll1][ll2][ll3][ll4][ll5+i];
    del_ft5dim(src, ll1, ul1, ll2, ul2, ll3, ul3, ll4, ul4, ll5, ul5);
  }

  return dst;
}

fp_t ****** atmosphere::add(fp_t ******src, fp_t ******dst, int32_t ll1, int32_t ul1, int32_t ll2, int32_t ul2, int32_t ll3, int32_t ul3, int32_t ll4, int32_t ul4, int32_t ll5, int32_t ul5, int32_t ll6, int32_t ul6){
  if(src){
    int32_t nn = (ul1 - ll1 + 1) * (ul2 - ll2 + 1) * (ul3 - ll3 + 1) * (ul4 - ll4 + 1) * (ul5 - ll5 + 1) * (ul6-ll6+1);
    for (int32_t i=0;i<nn;++i)
      dst[ll1][ll2][ll3][ll4][ll5][ll6+i] += src[ll1][ll2][ll3][ll4][ll5][ll6+i];
    del_ft6dim(src, ll1, ul1, ll2, ul2, ll3, ul3, ll4, ul4, ll5, ul5,ll6,ul6);
  }
  return dst;
}

fp_t ******* atmosphere::add(fp_t *******src, fp_t *******dst, int32_t ll1, int32_t ul1, int32_t ll2, int32_t ul2, int32_t ll3, int32_t ul3, int32_t ll4, int32_t ul4, int32_t ll5, int32_t ul5, int32_t ll6, int32_t ul6, int32_t ll7, int32_t ul7){
  if(src){
    int32_t nn = (ul1 - ll1 + 1) * (ul2 - ll2 + 1) * (ul3 - ll3 + 1) * (ul4 - ll4 + 1) * (ul5 - ll5 + 1) * (ul6-ll6+1)*(ul7-ll7+1);
    for (int32_t i=0;i<nn;++i)
      dst[ll1][ll2][ll3][ll4][ll5][ll6][ll7+i] += src[ll1][ll2][ll3][ll4][ll5][ll6][ll7+i];
    del_ft7dim(src, ll1, ul1, ll2, ul2, ll3, ul3, ll4, ul4, ll5, ul5,ll6,ul6,ll7,ul7);
  }
  return dst;
}
   
int08_t atmosphere::init(const char *wd,io_class *io_in)
{
  if(flags.is_clr(ATMOS_FLAG_DEF)) return -1; // problem not defined!
  if(flags.is_set(ATMOS_FLAG_INIT)){
    io_in->msg(IOL_INFO,"atmosphere::atmosphere: %s: already initialized!\n",id);
    return 0; // problem already initialized
  }
  io_in->msg(IOL_INFO,"atmosphere::atmosphere: ID=%s\n",id);
  io_in->msg(IOL_INFO,"atmosphere::atmosphere: FILE=%s\n",fname);
//
  read_atmos(wd,fname,ftype,io_in);
//
  return 0;
}

void atmosphere::set_Temp(int x1i, int x2i, int x3i, fp_t Temp_in){

  T[x1i][x2i][x3i] = Temp_in;
}

fp_t atmosphere::get_Nt(int x1i, int x2i, int x3i){
  return Nt[x1i][x2i][x3i];
}

void atmosphere::set_Nt(int x1i, int x2i, int x3i, fp_t Nt_in){
  Nt[x1i][x2i][x3i] = Nt_in;
}

fp_t atmosphere::get_vt(int x1i, int x2i, int x3i){
  return Vt[x1i][x2i][x3i];
}

int atmosphere::get_N_depths(){
  return x3h-x3l+1;
}

fp_t atmosphere::get_partf(int who, int z, fp_t T_in, fp_t Ne_in){
  return atml[who]->get_partf(z,T_in,Ne_in);
}

void atmosphere::set_conserve_charge(int input){
  conserve_charge = input;
}
int atmosphere::get_conserve_charge(){
  return conserve_charge;
}

fp_t ** atmosphere::return_as_array(){
  fp_t ** atmos;
  int ND = x3h-x3l+1;
  int NP = 12;
  atmos = ft2dim(1,NP,1,ND);
  for (int x3i=x3l;x3i<=x3h;++x3i){
    int i=x3i-x3l+1;
    atmos[1][i] = log10(-tau_referent[x1l][x2l][x3i]);
    atmos[2][i] = x3[x3i];
    atmos[3][i] = T[x1l][x2l][x3i];
    atmos[4][i] = Nt[x1l][x2l][x3i]*k*T[x1l][x2l][x3i]; // We want this as total gaspressure
    atmos[5][i] = Ne[x1l][x2l][x3i]*k*T[x1l][x2l][x3i]; // We want this as electron pressure
    atmos[6][i] = 0.0;
    atmos[7][i] = op_referent[x1l][x2l][x3i];
    fp_t B = sqrt(Bx[x1l][x2l][x3i]*Bx[x1l][x2l][x3i] + By[x1l][x2l][x3i]*By[x1l][x2l][x3i] + Bz[x1l][x2l][x3i]*Bz[x1l][x2l][x3i]);
    atmos[8][i] = B;
    atmos[9][i] = Vt[x1l][x2l][x3i];
    atmos[10][i] = Vz[x1l][x2l][x3i];
    atmos[11][i] = acos(Bz[x1l][x2l][x3i]/B);
    atmos[12][i] = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3l]);
  }
  return atmos;
}

int atmosphere::copy_from_array(fp_t ** array){
  
  int ND = x3h-x3l+1;
  int NP = 12;
  for (int x3i=x3l;x3i<=x3h;++x3i){
    int i=x3i-x3l+1;
    tau_referent[x1l][x2l][x3i] = -pow(10.0,array[1][i]);
    x3[x3i] = array[2][i];
    T[x1l][x2l][x3i] = array[3][i];
    Nt[x1l][x2l][x3i] = array[4][i];
    Ne[x1l][x2l][x3i] = array[5][i];
    op_referent[x1l][x2l][x3i] = array[7][i];
    //array[8][i] = B;
    Vt[x1l][x2l][x3i] = array[9][i];
    Vz[x1l][x2l][x3i] = array[10][i];
    //array[11][i] = acos(Bz[x1l][x2l][x3i]/B);
    //array[12][i] = atan(By[x1l][x2l][x3i]/Bx[x1l][x2l][x3l]);
  }
  return 0;
}

int atmosphere::build_from_nodes(model *){
  return 0;
}
int atmosphere::interpolate_from_nodes(model *){
  return 0;
}

// ----------------------------------------------------------------------------------

// Some additional useful i/o:

void atmosphere::print_atmos(){

  // Prints a 1D atmosphere or a first column of a 3D atmosphere on a stderr

  for (int x3i=x3l; x3i<=x3h; ++x3i)

    fprintf(stderr, "%1.5e %1.5e %1.5e %1.5e %1.5e \n",x3[x3i] , T[x1l][x2l][x3i],
      Nt[x1l][x2l][x3i], Vz[x1l][x2l][x3i], Vt[x1l][x2l][x3i]);


}
