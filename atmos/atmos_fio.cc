#include <math.h>
#include <errno.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "const.h"
#include "ana_io.h"
#include "fileio.h"
#include "mathtools.h"

#include "atmos.h"

int08_t atmosphere::read_atmos(const char *wd_in,const char *fname_in,uint08_t ftype_in,io_class *io_in)
{
  switch(ftype_in){
    case(ATMOS_TYPE_SPINOR): return read_spinor(wd_in,fname_in,io_in);
    case(ATMOS_TYPE_SPINOR3D): return read_spinor3d(wd_in,fname_in,io_in);
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

int08_t atmosphere::read_spinor3d(const char *wd,const char *filename,io_class *io_in)
{
  char *path=new char [strlen(wd)+strlen(filename)+2];
  sprintf(path,"%s/%s",wd,filename);
  io_in->msg(IOL_INFO,"atmosphere::read_spinor3d: %s\n",path);
  int32_t n1,n2,n3,n4;
  fp_t **** atmoscube = read_file(path,n1,n2,n3,n4,*io_in);
  atmoscube = transpose(atmoscube,n1,n2,n3,n4);
  int np = n4;
  int nx = n3;
  int ny = n2;
  int nz = n1;
  io_in->msg(IOL_INFO,"atmosphere::read_spinor3d: n_parameters=%d, nx=%d, ny=%d, nz=%d\n",np,nx,ny,nz);
  if (np!=12){
    io_in->msg(IOL_ERROR,"atmosphere::read_spinor3d: data cube does not have exactly 12 parameters.\n");
    return -1;  
  }
  
  this->resize(1,nx,1,ny,1,nz);

  for (int x3i=1;x3i<=nz;++x3i)
    x3[x3i] = atmoscube[2][1][1][x3i];

  for (int x1i=1;x1i<=nx;++x1i)
    for (int x2i=1;x2i<=ny;++x2i)
      for (int x3i=1;x3i<=nz;++x3i){
        tau_referent[x1i][x2i][x3i] = -pow(10.0,atmoscube[1][x1i][x2i][x3i]);
        T[x1i][x2i][x3i] = atmoscube[3][x1i][x2i][x3i];
        Nt[x1i][x2i][x3i] = atmoscube[4][x1i][x2i][x3i]/(k*T[x1i][x2i][x3i]);
        Ne[x1i][x2i][x3i] = atmoscube[5][x1i][x2i][x3i]/(k*T[x1i][x2i][x3i]);
        rho[x1i][x2i][x3i] = atmoscube[6][x1i][x2i][x3i];
        op_referent[x1i][x2i][x3i] = atmoscube[7][x1i][x2i][x3i];
        fp_t B = atmoscube[8][x1i][x2i][x3i];
        fp_t el = atmoscube[11][x1i][x2i][x3i];
        fp_t az = atmoscube[12][x1i][x2i][x3i];
        Bx[x1i][x2i][x3i]=B*sin(el)*cos(az);
        By[x1i][x2i][x3i]=B*sin(el)*sin(az);
        Bz[x1i][x2i][x3i]=B*cos(el);
        Vt[x1i][x2i][x3i] = atmoscube[9][x1i][x2i][x3i];;
        Vx[x1i][x2i][x3i] = 0.0;
        Vy[x1i][x2i][x3i] = 0.0;
        Vz[x1i][x2i][x3i] = atmoscube[10][x1i][x2i][x3i];
  }
  // Test:
  //for (int x3i=x3l;x3i<=x3h;++x3i)
  //   printf("%e %e %e %e %e \n",tau_referent[x1l][x2l][x3i],T[x1l][x2l][x3i],Nt[x1l][x2l][x3i],Bz[x1l][x2l][x3i],Vz[x1l][x2l][x3i]);
  del_ft4dim(atmoscube,1,np,1,nx,1,ny,1,nz);
  return 0;
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
