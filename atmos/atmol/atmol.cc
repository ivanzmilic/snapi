#include <string.h>

#include "pack.h"
#include "struts.h"

#include "atmol.h"
#include "mol/mol.h"
#include "mol/molcfg.h"
#include "atom/atom.h"
#include "atom/atomcfg.h"

atmol **append(atmol *q,atmol **p,int &len)
{
  //printf("atmol::append started, len = %d \n ",len);
  for(int i=0;i<len;++i) if(p[i]==q) return p;
  atmol **tmp=new atmol* [len+1];
  if(len){
    memcpy(tmp,p,len*sizeof(atmol*));
    delete[] p;
  }
  tmp[len++]=q;
  return tmp;
}

atmol *atmol_new(atmcfg *cfg,io_class &io_in)
{
  return atom_new(cfg,io_in); 
}

atmol *atmol_new(molcfg *cfg,atmol **atm,int na,io_class &io_in)
{
  return mol_new(cfg,atm,na,io_in);
}

atmol *atmol_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,atmol **atm,int na,io_class &io_in)
{
  uint64_t numid; // z atom 1, z atom 2 * 255, etc...
  ::unpack(buf+offs,numid,do_swap);
//  io_in.msg(IOL_INFO,"atmol_new: %lX %d %d\n",buf,offs,numid);
  if(numid>>8) return mol_new(numid,buf,offs,do_swap,atm,na,io_in); // molecules have high numid
  return atom_new(numid,buf,offs,do_swap,io_in); // atoms have even numid
}

atmol::atmol(const char *name_in,const char *id_in,io_class &io_in):io(io_in)
{
  name=dup_str(name_in);
  id=dup_str(id_in);
  io_in.msg(IOL_DEB1,"atmol::atmol: %s [%s]\n",name,id);
}

atmol::atmol(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):io(io_in)
{
  offs+=unpack(buf+offs,do_swap,io_in);
  io_in.msg(IOL_DEB1,"atmol::atmol: %s [%s]\n",name,id);
}

atmol::~atmol(void)
{
  if(id) delete[] id;
  if(name) delete[] name;
}

int32_t atmol::size(io_class &io_in)
{
  int32_t sz=sizeof(uint64_t); // numid
  sz+=(id)?strlen(id)+1:1;     // id
  sz+=(name)?strlen(name)+1:1; // name
  sz+=sizeof(fp_t); // mass
//  io.msg(IOL_INFO,"atmol::size: %s %s %d\n",name,id,sz);
  return sz;
}

int32_t atmol::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
  int32_t offs=::pack(buf,numid,do_swap);
  offs+=::pack(buf+offs,name);
  offs+=::pack(buf+offs,id);
  offs+=::pack(buf+offs,mass,do_swap);

//  io.msg(IOL_INFO,"atmol::pack: %d\n",offs);
  
  return offs;
}

int32_t atmol::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{
  int32_t offs=::unpack(buf,numid,do_swap);
  offs+=::unpack(buf+offs,name);
  offs+=::unpack(buf+offs,id);
  offs+=::unpack(buf+offs,mass,do_swap);
  
//  io.msg(IOL_INFO,"atmol::unpack: %d\n",offs);

  return offs;
}

// --------------------------------------------------------------------------------------------------------------------------------------------------

// Procedures dealing with parent atmosphere:
int atmol::set_parent_atmosphere(atmosphere * atm_in){

  parent_atm = atm_in;

  return 0;
}

int atmol::clear_parent_atmosphere(){

  parent_atm = 0;

  return 0;
}

fp_t atmol::fetch_population(int x1i, int x2i, int x3i, int species, int z, int i){

  return parent_atm->get_pop(x1i, x2i, x3i, species, z, i);
}

fp_t atmol::fetch_population(int x1i, int x2i, int x3i, int species, int z){
  return parent_atm->get_pop(x1i,x2i,x3i, species, z);
}

fp_t atmol::fetch_temperature(int x1i, int x2i, int x3i){
  return parent_atm->get_T(x1i, x2i, x3i);
}

fp_t atmol::fetch_Ne(int x1i, int x2i, int x3i){
  return parent_atm->get_Ne(x1i, x2i, x3i);
}

fp_t atmol::fetch_vt(int x1i, int x2i, int x3i){
  return parent_atm->get_vt(x1i, x2i, x3i);
}

fp_t * atmol::fetch_magnetic_field(int x1i, int x2i, int x3i){
  return parent_atm->get_magnetic_field(x1i, x2i, x3i);
}

fp_t  atmol::fetch_Nt(int x1i, int x2i, int x3i){
  return parent_atm->get_Nt(x1i,x2i,x3i);
}

//
//
//

fp_t *atmol::opacity(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t)
{
  return 0;
}

fp_t *atmol::emissivity(fp_t,fp_t,fp_t*,int32_t,int32_t,int32_t,int32_t)
{
  return 0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------

fp_t ***atmol::opacity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t)
{
//  io.msg(IOL_INFO,"atmol::opacity: %s\n",id);
  return 0;
}


fp_t ***atmol::emissivity(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t)
{
//  io.msg(IOL_INFO,"atmol::emissivity: %s\n",id);
  return 0;
}

fp_t ***** atmol::opacity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t){
  return 0;
}

fp_t ***** atmol::emissivity_pert(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t){
  return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------

fp_t *****atmol::opacity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t){
  return 0;
}

fp_t ****atmol::emissivity_vector(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t)
{
//  io.msg(IOL_INFO,"atmol::emissivity: %s\n",id);
  return 0;
}

uint08_t atmol::rtsetup(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t)
{
  return 0;
}

uint08_t atmol::rtclean(int,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t)
{
  return 0;
}

void atmol::rtinit(void)
{
}

void atmol::prof_init(void){
}

fp_t atmol::pops(atmol**,uint16_t,fp_t,fp_t,int32_t,int32_t,int32_t, int)
{
  return 0;
}
