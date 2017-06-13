#include <math.h>
#include <errno.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "const.h"

#include "atmos.h"

int08_t atmosphere::read_atmos(const char *wd_in,const char *fname_in,uint08_t ftype_in,io_class *io_in)
{
  switch(ftype_in){
    case(ATMOS_TYPE_SPINOR): return read_spinor(wd_in,fname_in,io_in);
    case(ATMOS_TYPE_MHD): return read_mhd(wd_in,fname_in,io_in);
    case(ATMOS_TYPE_MURAM): return read_muram(wd_in,fname_in,io_in);
  }
  return -1;
}

int08_t atmosphere::read_spinor(const char *wd,const char *filename,io_class *io_in)
{
  char *path=new char [strlen(wd)+strlen(filename)+2];
  sprintf(path,"%s/%s",wd,filename);
  io_in->msg(IOL_INFO,"atmosphere::read_spinor: %s\n",path);
  if(FILE *f=fopen(path,"r")){
    char line[1000];
    if(!fgets(line,1000,f)) io_in->msg(IOL_WARN,"atmosphere::read_spinor: when reading line 1: %s\n",strerror(errno));
    int32_t nz;
    sscanf(line,"%d",&nz);
    io_in->msg(IOL_INFO,"atmosphere::read_spinor: %d depth points\n",nz);
//
// resize the grid...
//
    this->resize(1,1,1,1,1,nz); // spinor input is strictly 1-D 
//
    for(int ix3=x3l;ix3<=x3h;++ix3){
      if(!fgets(line,1000,f)) io_in->msg(IOL_WARN,"atmosphere::read_spinor: when reading line %d: %s\n",ix3-x3l+2,strerror(errno));
      if(feof(f)){
        io_in->msg(IOL_ERROR,"atmosphere::read_spinor: unexpected end of file reading spinor file line %d\n",ix3-x3l+2);
        delete[] path;
        fclose(f);
        return -1;
      }
      float32_t lt,hh,tm,pg,pe,aa,rr,bb,vt,vd,el=0.0,az=0.0;
      int nn=sscanf(line,"%E %E %E %E %E %E %E %E %E %E %E %E",&lt,&hh,&tm,&pg,&pe,&aa,&rr,&bb,&vt,&vd,&el,&az);
      if(nn<10) io_in->msg(IOL_WARN,"atmosphere::read_spinor: when reading line %d: only %d successful conversions\n",ix3,nn);
//
      io_in->msg(IOL_XNFO,"atmosphere::read_spinor: %E %E %E %E %E %E %E %E %E %E %E %E\n",lt,hh,tm,pg,pe,aa,rr,bb,vt,vd,el,az);
//
      tau_referent[x1l][x2l][ix3] = -pow(10.0,lt);
      x3[ix3]=hh;
      T[x1l][x2l][ix3]=tm;
      rho[x1l][x2l][ix3]=rr;
      Nt[x1l][x2l][ix3]=pg/(k*tm);
      Ne[x1l][x2l][ix3]=pe/(k*tm);
      Bx[x1l][x2l][ix3]=bb*sin(el)*cos(az);
      By[x1l][x2l][ix3]=bb*sin(el)*sin(az);
      Bz[x1l][x2l][ix3]=bb*cos(el);
      Vx[x1l][x2l][ix3]=0.0;
      Vy[x1l][x2l][ix3]=0.0;
      Vz[x1l][x2l][ix3]=vd;
      Vt[x1l][x2l][ix3]=vt;
    }
    fclose(f);
    delete[] path;
    return 0;
  }
  io_in->msg(IOL_ERROR,"atmosphere::read_spinor: error opening file: %s\n",strerror(errno));
  delete[] path;
  return -1;
}

int08_t atmosphere::read_mhd(const char *wd,const char *filename,io_class *io_in)
{
  io_in->msg(IOL_INFO,"atmosphere::read_mhd: reading MHD cube [not implemented yet]\n");
  return 0;
}

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

int08_t atmosphere::read_muram(const char *wd,const char *number,io_class *io_in)
{
  const char *prfx[]={"eosP","eosT","result_0","result_1","result_2","result_3","result_4","result_5","result_6","result_7",0};
  io_in->msg(IOL_INFO,"atmosphere::read_muram: snapshot %s\n",number);
  io_in->msg(IOL_INFO,"atmosphere::read_muram: %d..%d  %d..%d  %d..%d\n",x1l,x1h,x2l,x2h,x3l,x3h);
  for(int i=0;prfx[i];++i){
    char *path=new char [strlen(wd)+strlen(prfx[i])+strlen(number)+3];
    sprintf(path,"%s/%s.%s",wd,prfx[i],number);
    io_in->msg(IOL_XNFO,"atmosphere::read_muram: \"%s\"\n",path);
    struct stat status;
    if(stat(path,&status)<0){
      io_in->msg(IOL_ERROR,"atmosphere::read_muram: failed to recover file status for \"%s\": %s\n",path,strerror(errno));
      delete[] path;
      return -1;
    }
    delete[] path;
  }
//
  fp_t ***p[]={Nt,T,rho,Bx,By,Bz,Vx,Vy,Vz,0};
  float32_t ***q=f32t3dim(x1l,x1h,x2l,x2h,x3l,x3h);
  for(int i=0;p[i];++i){
    char *path=new char [strlen(wd)+strlen(prfx[i])+strlen(number)+3];
    sprintf(path,"%s/%s.%s",wd,prfx[i],number);
    int fd=open(path,O_RDONLY);
    if(fd<0){
      io_in->msg(IOL_ERROR,"atmosphere::read_muram: failed to open \"%s\": %s\n",path,strerror(errno));
      delete[] path;
      return -1;
    }
    uint64_t nelem=(uint64_t)(x1h-x1l+1)*(uint64_t)(x2h-x2l+1)*(uint64_t)(x3h-x3l+1);
    off_t rdsz=read(fd,q[x1l][x2l]+x3l,nelem*sizeof(float32_t));
    if(rdsz<(nelem*sizeof(float32_t))) io_in->msg(IOL_WARN,"atmosphere::read_muram: short read of \"%s\" (%ld out of %ld bytes read): %s\n",path,rdsz,(nelem*sizeof(float32_t)),strerror(errno));
    if(close(fd)) io_in->msg(IOL_WARN,"atmosphere::read_muram: failed to close \"%s\": %s\n",path,strerror(errno));
    for(int ix1=x1l;ix1<=x1h;++ix1)
      for(int ix2=x2l;ix2<=x2h;++ix2)
        for(int ix3=x3l;ix3<=x3h;++ix3) p[i][ix1][ix2][ix3]=q[ix1][ix2][ix3]; // type conversion here
  }
  for(int ix1=x1l;ix1<=x1h;++ix1)
    for(int ix2=x2l;ix2<=x2h;++ix2)
      for(int ix3=x3l;ix3<=x3h;++ix3){
        Nt[ix1][ix2][ix3]/=(k*T[ix1][ix2][ix3]);
        Ne[ix1][ix2][ix3]=1.0;
        Vt[x1l][x2l][ix3]=0.0;
      }
  del_f32t3dim(q,x1l,x1h,x2l,x2h,x3l,x3h);
  io_in->msg(IOL_INFO,"atmosphere::read_muram: done.\n");
  return 0;
}
