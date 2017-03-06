//
//  C++ translation of the original mfbd code written in ANA
//
//  Author: Michiel van Noort
//  Created: 2003.12.03
//  Original Files:     mfbd.ana recipes.ana mfbd_read_images.ana init.ana
//	   Created:    <2002-01-05 ~tE:~tm:~ts ~ui>
//	    Author:     mats@astro.su.se
//      Time-stamp:    <2003-06-19 13:52:25 mats>
//

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "uts.h"
#include "mem.h"
#include "fileio.h"
#include "io.h"
#include "cmdln.h"
#include "nr/gammln.h"
#include "nr/svdcmp.h"
#include "conf.h"
#include "consts.h"
#include "modes.h"
#include "pack.h"
#include "compress.h"

//
// wavefront modes
//

void zernike_mn(int j,int &m,int &n) // j=1...
{
  n=0;
  int len=1;
  for(int i=1;len<j;++i) len+=(n=i)+1;
  int dl=n+1-len+j;
  m=2*((dl+(n%2))/2)+!(n%2)-1;
}

fp_t zernike_covar(int j)
{
  if(j<2) return 0.0;
  int m,n;
  zernike_mn(j,m,n);
  int n1=n;
//  ; Now deal with the numerical terms: Dai
  fp_t tmp;
  fp_t k=pow(4.8*exp(gammln(6.0/5.0,tmp)),5.0/6.0)*exp(gammln(14.0/3.0,tmp)+2.0*gammln(11.0/6.0,tmp))/(pow(2.0,(8.0/3.0))*PI);
  k*=pow(-1.0,(fp_t)((n+n1-2*m)/2))*sqrt((fp_t)(n+1)*(n1+1));
  fp_t g1_sgn,g2_sgn,g3_sgn,g4_sgn;
  fp_t g1=gammln(((fp_t)(n+n1)- 5.0/3.0)/2.0,g1_sgn);
  fp_t g2=gammln(((fp_t)(n-n1)+17.0/3.0)/2.0,g2_sgn);
  fp_t g3=gammln(((fp_t)(n1-n)+17.0/3.0)/2.0,g3_sgn);
  fp_t g4=gammln(((fp_t)(n+n1)+23.0/3.0)/2.0,g4_sgn);
  return k*exp(g1-g2-g3-g4)*g1_sgn*g2_sgn*g3_sgn*g4_sgn;
}

fp_t zernike_covar(int i,int j)
{
  if((i<2)||(j<2)) return 0.0;
  int m,n,o,p;
  zernike_mn(i,m,n);
  zernike_mn(j,o,p);
  if(m!=o) return 0.0;
  if(m) if((i+j)%2) return 0.0;
//  ; Now deal with the numerical terms: Dai
  fp_t tmp;
  fp_t k=pow(4.8*exp(gammln(6.0/5.0,tmp)),5.0/6.0)*exp(gammln(14.0/3.0,tmp)+2.0*gammln(11.0/6.0,tmp))/(pow(2.0,(8.0/3.0))*PI);
  k*=pow(-1.0,(fp_t)((n+p-2*m)/2))*sqrt((fp_t)((n+1)*(p+1)));
  fp_t g1_sgn,g2_sgn,g3_sgn,g4_sgn;
  fp_t g1=gammln(((fp_t)(n+p)- 5.0/3.0)/2.0,g1_sgn);
  fp_t g2=gammln(((fp_t)(n-p)+17.0/3.0)/2.0,g2_sgn);
  fp_t g3=gammln(((fp_t)(p-n)+17.0/3.0)/2.0,g3_sgn);
  fp_t g4=gammln(((fp_t)(n+p)+23.0/3.0)/2.0,g4_sgn);
  return k*exp(g1-g2-g3-g4)*g1_sgn*g2_sgn*g3_sgn*g4_sgn;
}

fp_t kl_covar(int j,struct klmc *cfs)
{
  return cfs[j].v;
}

fp_t *RC(int m,int n)
{
  int nmm=(n-m)/2,npm=(n+m)/2;
  int nmax=max(npm,n);
  fp_t *f=new fp_t [nmax+1];
  f[0]=1.0;
  for(int i=1;i<=nmax;++i) f[i]=(fp_t)i*f[i-1];
  fp_t *res=new fp_t [nmm+1]; 
  for(int s=0,pm=-1;s<=nmm;++s)
    res[s]=(fp_t)((pm*=-1)*f[n-s])/(f[s]*f[npm-s]*f[nmm-s]);
  delete[] f;
  return res;
}

fp_t **zernike_function(fp_t lambda,fp_t r_c,int nph,int j,fp_t angle)
{
  int xo=1+nph/2,yo=1+nph/2;
  fp_t **z=ft2dim(1,nph,1,nph);
  if(j==1)
    for(int x=1;x<=nph;++x)
      for(int y=1;y<=nph;++y) z[x][y]=1.0;
  else{                                       // j>1
    int m,n;
    zernike_mn(j,m,n);
    fp_t *rc=RC(m,n);
    fp_t **r=ft2dim(1,nph,1,nph);
    fp_t **rs=ft2dim(1,nph,1,nph);
    for(int x=1;x<=nph;++x)                   // s=0
      for(int y=1;y<=nph;++y){
        fp_t rr=sqrt((fp_t)sqr(x-xo)+(fp_t)sqr(y-yo))/r_c; 
        rs[x][y]=rr*rr;
        r[x][y]=pow(rr,n);
        z[x][y]=r[x][y]*rc[0];
      }
    rs[xo][yo]=1.0;                           // avoid division by 0
    for(int s=1;s<=(n-m)/2;++s){
      for(int x=1;x<=nph;++x)
        for(int y=1;y<=nph;++y)
          z[x][y]+=(r[x][y]/=rs[x][y])*rc[s];
      if(!(n-2*s)) z[xo][yo]+=rc[s];          // dividing 0 by 1 will never give 1...
    }
    del_ft2dim(rs,1,nph,1,nph);
    del_ft2dim(r,1,nph,1,nph);
    if(m){                                    // m!=0
      fp_t sf=sqrt((fp_t)(2*(n+1)));
      if(j%2)                                 // odd
        for(int x=1;x<=nph;++x)
          for(int y=1;y<=nph;++y) z[x][y]*=sf*sin(((fp_t)m)*(atan2((fp_t)(y-yo),(fp_t)(x-xo))+angle));
      else                                    // even
        for(int x=1;x<=nph;++x)
          for(int y=1;y<=nph;++y) z[x][y]*=sf*cos(((fp_t)m)*(atan2((fp_t)(y-yo),(fp_t)(x-xo))+angle));
    }else{                                    // m==0
      fp_t sf=sqrt((fp_t)(n+1));
      for(int x=1;x<=nph;++x)
        for(int y=1;y<=nph;++y) z[x][y]*=sf;
    }
    delete[] rc;
  }
  fp_t sum=0.0,N=0.0,rcs=r_c*r_c,dx=0.5/r_c,dy=0.5/r_c;
  for(int x=1;x<=nph;++x){
    fp_t xl=fabs((fp_t)(x-xo))/r_c-dx,xh=fabs((fp_t)(x-xo))/r_c+dx;
    fp_t xhs=sqr(xh);
    for(int y=1;y<=nph;++y){
      fp_t yl=fabs((fp_t)(y-yo))/r_c-dy,yh=fabs((fp_t)(y-yo))/r_c+dy;
      fp_t yhs=sqr(yh);
      fp_t rsl=sqr(xl)+sqr(yl),rsh=xhs+yhs;
      int ti=(rsl<rcs)+(rsh<rcs);
      if(rsl<=1.0)       // good pixel
        if(rsh<1.0){     // full pixel
          sum+=sqr(z[x][y]);
          N+=1.0;
        }else{           // partial pixel
          fp_t x2=sqrt(max(1.0-yhs,(fp_t)0.0));
          fp_t y3=sqrt(max(1.0-xhs,(fp_t)0.0));
          fp_t f=(xh>yh)?(yh-yl)*(min(xh,max(xl,x2))-xl)/(4*dx*dy):
                         (xh-xl)*(min(yh,max(yl,y3))-yl)/(4*dx*dy);
          sum+=f*sqr(z[x][y]);
          N+=f;
        }
    }
  }
  sum/=N;
  for(int x=1;x<=nph;++x)
    for(int y=1;y<=nph;++y) z[x][y]/=lambda*sqrt(sum);
  return z;
}

fp_t **kl_function(fp_t lambda,fp_t r_c,int nph,int j,struct klmc *cfs,fp_t cutoff,fp_t ***z,fp_t angle)
{
  struct klmc cf=cfs[j];
  fp_t **kl=ft2dim(1,nph,1,nph),c;
  memset(kl[1]+1,0,nph*nph*sizeof(fp_t));
  for(int m=1;m<=cf.nm;++m)
    if(fabs(c=cf.c[m])>=cutoff){
      int mode=cf.m[m];
      if(!z[mode]) z[mode]=zernike_function(lambda,r_c,nph,mode,angle);
      for(int x=1;x<=nph;++x)
        for(int y=1;y<=nph;++y) kl[x][y]+=c*z[mode][x][y];
    }
  return kl;
}

int *mapping_vector(fp_t func(int,int),int nl,int nh)
{                               // Mapping to make matrix block diagonal
  int n=nh-nl+1;
  int *map=new int [n]-nl;
  for(int i=nl;i<=nh;++i) map[i]=i;
  for(int i=nl;i<=nh-1;++i){
    fp_t fp=func(map[i],map[nh]);
    for(int j=nh-1;j>=i;--j){
      fp_t f=func(map[i],map[j]);
      if((!f)&&(fp))
        swap(map[j],map[j+1]);
      else
        fp=f;
    }
  }
  return map;
}

int *reverse_mapping_vector(int *map,int nl,int nh)
{
  int n=nh-nl+1;
  int *rev_map=new int [n]-nl;
  for(int i=nl;i<=nh;++i) rev_map[map[i]]=i;
  return rev_map;
}

fp_t ***reordered_matrix(fp_t func(int,int),int *map,int nl,int nh,int &nsm,int *&nnl,int *&nnh)
{
  nsm=1;
  fp_t *off_diag=new fp_t [nh-nl+1]-nl-1;
  for(int i=nl+1;i<=nh;++i)
    if((off_diag[i]=func(map[i-1],map[i]))==0) ++nsm;
  off_diag[nh+1]=0;
  nnl=new int [nsm]-1;
  nnh=new int [nsm]-1;
  fp_t ***sm=new fp_t** [nsm]-1;
  for(int s=1,n=nl;s<=nsm;++s){
    nnl[s]=n++;
    while((off_diag[n])&&(n<=nh)) ++n;
    nnh[s]=n-1;
    sm[s]=ft2dim(nnl[s],nnh[s],nnl[s],nnh[s]);
    for(int i=nnl[s];i<=nnh[s];++i) sm[s][i][i]=func(map[i],map[i]);
    for(int i=nnl[s]+1;i<=nnh[s];++i) sm[s][i-1][i]=sm[s][i][i-1]=off_diag[i];
    for(int i=nnl[s];i<=nnh[s]-2;++i)
      for(int j=i+2;j<=nnh[s];++j) sm[s][i][j]=sm[s][j][i]=func(map[i],map[j]);
  }
  delete[] (off_diag+nl+1);
  return sm;
}

struct klmc *kl_coefs(fp_t ***sm,fp_t *wr,int *nnl,int *nnh,int *map,int *rmap,int nl,int nh)
{
  struct klmc *cfs=new klmc [nh-nl+1]-nl;
  for(int i=nl;i<=nh;++i){
    int im=rmap[i],s;
    for(s=1;(nnl[s]>im)||(nnh[s]<im);++s);
    int n=nnh[s]-nnl[s]+1;
    cfs[i].nm=n;
    cfs[i].v=wr[im];
    cfs[i].c=new fp_t [n]-1;
    cfs[i].m=new int [n]-1;
    for(int m=1;m<=n;++m){
      int j=m+nnl[s]-1;
      cfs[i].m[m]=map[j];
      cfs[i].c[m]=sm[s][j][im];
    }
  }
  return cfs;
}

void svdcmp_block(fp_t ***sm,fp_t *wr,int nsm,int *nnl,int *nnh)
{ // block-wise SVDCMP routine for use with block-diagonal square symmetrical matrix
  for(int s=1;s<=nsm;++s){
    int n=nnh[s]-nnl[s]+1;
    if(n>1){
      fp_t **u=ft2dim(1,n,1,n);
      fp_t **v=ft2dim(1,n,1,n);
      fp_t *w=new fp_t [n]-1;
      int *idx=new int [n]-1;
//
      memcpy(u[1]+1,sm[s][nnl[s]]+nnl[s],n*n*sizeof(fp_t));
      svdcmp(u,n,n,w,v);
      for(int i=1;i<=n;++i) idx[i]=i;
      sort(w,idx,n);
      for(int i=1;i<=n;++i) wr[i+nnl[s]-1]=w[n-i+1];
      for(int i=1;i<=n;++i) memcpy(sm[s][i+nnl[s]-1]+nnl[s],u[idx[n-i+1]]+1,n*sizeof(fp_t));
//
      delete[] (idx+1);
      delete[] (w+1);
      del_ft2dim(v,1,n,1,n);
      del_ft2dim(u,1,n,1,n);
    }else{
      wr[nnl[s]]=sm[s][nnl[s]][nnl[s]];
      sm[s][nnl[s]][nnl[s]]=1.0;
    }
  }
}

struct klmc *kl_cfg(int kl_min_mode,int kl_max_mode)
{
  int *map=mapping_vector(zernike_covar,kl_min_mode,kl_max_mode);
  int nsm,*nnl,*nnh;
  fp_t ***rm=reordered_matrix(zernike_covar,map,kl_min_mode,kl_max_mode,nsm,nnl,nnh);
  fp_t *wr=new fp_t [kl_max_mode-kl_min_mode+1]-kl_min_mode;
  svdcmp_block(rm,wr,nsm,nnl,nnh);
  int *rmap=reverse_mapping_vector(map,kl_min_mode,kl_max_mode);
  struct klmc *cfs=kl_coefs(rm,wr,nnl,nnh,map,rmap,kl_min_mode,kl_max_mode);
  delete[] (wr+kl_min_mode);
  for(int s=1;s<=nsm;++s) del_ft2dim(rm[s],nnl[s],nnh[s],nnl[s],nnh[s]);
  delete[] (rm+1);
  delete[] (nnl+1);
  delete[] (nnh+1);
  delete[] (map+kl_min_mode);
  delete[] (rmap+kl_min_mode);
  return cfs;
}

void kl_uncfg(struct klmc *cfs,int kl_min_mode,int kl_max_mode)
{
  for(int m=kl_min_mode;m<=kl_max_mode;++m){
    delete[] (cfs[m].m+1);
    delete[] (cfs[m].c+1);
  }
  delete[] (cfs+kl_min_mode);
}

modes::modes(struct klmc *kl_cfs,fp_t lambda,fp_t r_c,int nph_in,int basis,int nm,int *mode_num,int nch,int *ndo,int **dorder,int **dtype,int kl_min_mode,int kl_max_mode,fp_t svd_reg,fp_t angle,fp_t **pupil_in,io_class &io)
{
  nph=nph_in;
  zmin=0x7FFFFFFF;
  zmax=0x00000000;
  kmin=0x7FFFFFFF;
  kmax=0x00000000;
  mmin=0x7FFFFFFF;
  mmax=0x00000000;
  for(int c=1;c<=nch;++c){
    for(int m=1;m<=nm;++m){ // basis modes
      if(basis==CFG_ZERNIKE){
        zmin=min(zmin,mode_num[m]);
        zmax=max(zmax,mode_num[m]);
      }
      if(basis==CFG_KARHUNEN_LOEVE){
        kmin=min(kmin,mode_num[m]);
        kmax=max(kmax,mode_num[m]);
      }
      mmin=min(mmin,mode_num[m]);
      mmax=max(mmax,mode_num[m]);
    }
    for(int m=1;m<=ndo[c];++m){ // diversity modes
      if(dtype[c][m]==CFG_ZERNIKE){
        zmin=min(zmin,dorder[c][m]);
        zmax=max(zmax,dorder[c][m]);
      }
      if(dtype[c][m]==CFG_KARHUNEN_LOEVE){
        kmin=min(kmin,dorder[c][m]);
        kmax=max(kmax,dorder[c][m]);
      }
    }
  }
  mode=new fp_t*** [CFG_KARHUNEN_LOEVE+1];
  mode[0]=new fp_t** [mmax-mmin+1]-mmin;
  memset(mode[0]+mmin,0,(mmax-mmin+1)*sizeof(fp_t**));
  if(zmin<=zmax){
    mode[CFG_ZERNIKE]=new fp_t** [zmax-zmin+1]-zmin;
    memset(mode[CFG_ZERNIKE]+zmin,0,(zmax-zmin+1)*sizeof(fp_t**));
  }else mode[CFG_ZERNIKE]=0;
  if(kmin<=kmax){
    mode[CFG_KARHUNEN_LOEVE]=new fp_t** [kmax-kmin+1]-kmin;
    memset(mode[CFG_KARHUNEN_LOEVE]+kmin,0,(kmax-kmin+1)*sizeof(fp_t**));
  }else mode[CFG_KARHUNEN_LOEVE]=0;
  covar=new fp_t* [CFG_KARHUNEN_LOEVE+1];
  covar[0]=new fp_t [mmax-mmin+1]-mmin;
  if(zmin<=zmax) covar[CFG_ZERNIKE]=new fp_t [zmax-zmin+1]-zmin; else covar[CFG_ZERNIKE]=0;
  if(kmin<=kmax) covar[CFG_KARHUNEN_LOEVE]=new fp_t [kmax-kmin+1]-kmin; else covar[CFG_KARHUNEN_LOEVE]=0;
// use temporary storage z to avoid recomputing too many Zernikes
  fp_t ***z=new fp_t** [kl_max_mode-kl_min_mode+1]-kl_min_mode;
  memset(z+kl_min_mode,0,(kl_max_mode-kl_min_mode+1)*sizeof(fp_t**));
  for(int c=1;c<=nch;++c)
    for(int m=1,mn;m<=nm;++m)
     if(!mode[0][mn=mode_num[m]])
        if(basis==CFG_ZERNIKE){                     // Zernikes
          mode[0][mn]=(mode[CFG_ZERNIKE][mn]=zernike_function(lambda,r_c,nph,mn,angle));
          covar[0][mn]=(covar[CFG_ZERNIKE][mn]=sqrt(zernike_covar(mn)))*sqr(lambda);
	}else{                                          // Karhunen-Loeves for the rest
          if((mn==2)||(mn==3)){                         // use Zernikes for the tilts
            mode[0][mn]=(mode[CFG_KARHUNEN_LOEVE][mn]=zernike_function(lambda,r_c,nph,mn,angle));
            covar[0][mn]=(covar[CFG_KARHUNEN_LOEVE][mn]=sqrt(zernike_covar(mn)))*sqr(lambda);
          }else{
            mode[0][mn]=(mode[CFG_KARHUNEN_LOEVE][mn]=kl_function(lambda,r_c,nph,mn,kl_cfs,svd_reg,z,angle));
            covar[0][mn]=(covar[CFG_KARHUNEN_LOEVE][mn]=sqrt(kl_covar(mn,kl_cfs)))*sqr(lambda);
          }
        }
  for(int c=1;c<=nch;++c)                 // compute diversity orders
    for(int m=1,mn,type;m<=ndo[c];++m){
      if(!mode[type=dtype[c][m]][mn=dorder[c][m]])
        if(type==CFG_ZERNIKE){               // use Zernikes for the tilts
          mode[CFG_ZERNIKE][mn]=zernike_function(lambda,r_c,nph,mn,angle);
          covar[CFG_ZERNIKE][mn]=sqrt(zernike_covar(mn))*sqr(lambda);
	}else{                               // Karhunen-Loeves for the rest
          mode[CFG_KARHUNEN_LOEVE][mn]=kl_function(lambda,r_c,nph,mn,kl_cfs,svd_reg,z,angle);
          covar[CFG_KARHUNEN_LOEVE][mn]=sqrt(kl_covar(mn,kl_cfs))*sqr(lambda);
        }
    }
  for(int m=kl_min_mode;m<=kl_max_mode;++m) if(z[m]) del_ft2dim(z[m],1,nph,1,nph);
  delete[] (z+kl_min_mode);
// initialise pupil
  pupil=ft2dim(1,nph,1,nph);
  if(pupil_in){
    memcpy(pupil[1]+1,pupil_in[1]+1,nph*nph*sizeof(fp_t));
  }else{
    int xo=1+nph/2,yo=1+nph/2;
    fp_t rcs=r_c*r_c,dx=0.5/r_c,dy=0.5/r_c;
    for(int x=1;x<=nph;++x){
      fp_t xl=fabs((fp_t)(x-xo))/r_c-dx,xh=fabs((fp_t)(x-xo))/r_c+dx;
      fp_t xhs=sqr(xh);
      for(int y=1;y<=nph;++y){
        fp_t yl=fabs((fp_t)(y-yo))/r_c-dy,yh=fabs((fp_t)(y-yo))/r_c+dy;
        fp_t yhs=sqr(yh);
        fp_t rsl=sqr(xl)+sqr(yl),rsh=xhs+yhs;
        int ti=(rsl<rcs)+(rsh<rcs);
        if(rsl<=1.0){      // inside pixel
          if(rsh<1.0)      // full pixel
            pupil[x][y]=1.0;
          else{            // partial pixel
            fp_t x2=sqrt(max(1.0-yhs,(fp_t)0.0));
            fp_t y3=sqrt(max(1.0-xhs,(fp_t)0.0));
            fp_t f=(xh>yh)?(yh-yl)*(min(xh,max(xl,x2))-xl)/(4*dx*dy):
                           (xh-xl)*(min(yh,max(yl,y3))-yl)/(4*dx*dy);
            pupil[x][y]=f;
          }
        }else pupil[x][y]=0.0; // outside pixel
      }
    }
  }
//
  area=0.0;
  for(int x=1;x<=nph;++x)
    for(int y=1;y<=nph;++y) area+=pupil[x][y];
}

modes::modes(byte *buf,byte swap_endian,io_class &io)
{
  memset(this,0,sizeof(struct modes));
  io.msg(MFBD_INFO,"modes:modes: unpacking mode data\n");
  int size;
  byte *p=z_uncompress(buf,size,swap_endian,io);
// unpack p
  int offs=0;
  offs+=unpack(p+offs,area,swap_endian);
  offs+=unpack(p+offs,nph,swap_endian);
  mode=new float64_t*** [CFG_KARHUNEN_LOEVE+1];
  memset(mode,0,(CFG_KARHUNEN_LOEVE+1)*sizeof(float64_t***));
// Zernikes
  offs+=unpack(p+offs,zmin,swap_endian);
  offs+=unpack(p+offs,zmax,swap_endian);
  if(zmin<=zmax){
    mode[CFG_ZERNIKE]=new float64_t** [zmax-zmin+1]-zmin;
    for(int m=zmin;m<=zmax;++m){
      mode[CFG_ZERNIKE][m]=0;
      byte is_used;
      offs+=unpack(p+offs,is_used);
      if(is_used) offs+=unpack(p+offs,mode[CFG_ZERNIKE][m]=f64t2dim(1,nph,1,nph),1,nph,1,nph,swap_endian);
   }
  }
// Karhune-Loeves
  offs+=unpack(p+offs,kmin,swap_endian);
  offs+=unpack(p+offs,kmax,swap_endian);
  if(kmin<=kmax){
    mode[CFG_KARHUNEN_LOEVE]=new float64_t** [kmax-kmin+1]-kmin;
    for(int m=kmin;m<=kmax;++m){
      mode[CFG_KARHUNEN_LOEVE][m]=0;
      byte is_used;
      offs+=unpack(p+offs,is_used);
      if(is_used) offs+=unpack(p+offs,mode[CFG_KARHUNEN_LOEVE][m]=f64t2dim(1,nph,1,nph),1,nph,1,nph,swap_endian);
    }
  }
// link info
  offs+=unpack(p+offs,mmin,swap_endian);
  offs+=unpack(p+offs,mmax,swap_endian);
  if(mmin<=mmax){
    mode[0]=new float64_t** [mmax-mmin+1]-mmin;
    for(int m=mmin;m<=mmax;++m){
      mode[0][m]=0;
      byte is_used;
      offs+=unpack(p+offs,is_used);
      if(is_used){
        byte type=0;
        int index=-1;
        offs+=unpack(p+offs,type);
        offs+=unpack(p+offs,index,swap_endian);
        if(type) mode[0][m]=mode[type][index];
      }
    }
  }
//
  covar=new float64_t* [CFG_KARHUNEN_LOEVE+1];
  byte has_covar;
  offs+=unpack(p+offs,has_covar);
  if(has_covar) offs+=unpack(p+offs,covar[0]=new float64_t [mmax-mmin+1]-mmin,mmin,mmax,swap_endian); else covar[0]=0;
  offs+=unpack(p+offs,has_covar);
  if(has_covar) offs+=unpack(p+offs,covar[CFG_ZERNIKE]=new float64_t [zmax-zmin+1]-zmin,zmin,zmax,swap_endian); else covar[CFG_ZERNIKE]=0;
  offs+=unpack(p+offs,has_covar);
  if(has_covar) offs+=unpack(p+offs,covar[CFG_KARHUNEN_LOEVE]=new float64_t [kmax-kmin+1]-kmin,kmin,kmax,swap_endian); else covar[CFG_KARHUNEN_LOEVE]=0;
//
  byte has_pupil;
  offs+=unpack(p+offs,has_pupil);
  if(has_pupil) offs+=unpack(p+offs,pupil=f64t2dim(1,nph,1,nph),1,nph,1,nph,swap_endian); else pupil=0;
// end unpack
  if(offs!=size) io.msg(MFBD_WARN,"modes:modes: unpacked %d bytes but should be %d\n",offs,size);
  delete[] p;
}

modes::~modes(void)
{
  if(mode[0]) delete[] (mode[0]+mmin); // don't delete mode[0][m], they are co-pointers!
  for(int m=zmin;m<=zmax;++m) if(mode[CFG_ZERNIKE][m]) del_ft2dim(mode[CFG_ZERNIKE][m],1,nph,1,nph);
  if(mode[CFG_ZERNIKE]) delete[] (mode[CFG_ZERNIKE]+zmin);
  for(int m=kmin;m<=kmax;++m) if(mode[CFG_KARHUNEN_LOEVE][m]) del_ft2dim(mode[CFG_KARHUNEN_LOEVE][m],1,nph,1,nph);
  if(mode[CFG_KARHUNEN_LOEVE]) delete[] (mode[CFG_KARHUNEN_LOEVE]+kmin);
  delete[] mode;
//
  if(covar[0]) delete[] (covar[0]+mmin);
  if(covar[CFG_ZERNIKE]) delete[] (covar[CFG_ZERNIKE]+zmin);
  if(covar[CFG_KARHUNEN_LOEVE]) delete[] (covar[CFG_KARHUNEN_LOEVE]+kmin);
  delete[] covar;
  if(pupil) del_ft2dim(pupil,1,nph,1,nph);
}

byte *modes::compress(int &size,int level,byte swap_endian,io_class &io)
{
// calculate length
  int32_t usize=sizeof(float64_t)+sizeof(int32_t); // area,nph
// Zernikes
  usize+=2*sizeof(int32_t); // zmin zmax
  for(int m=zmin;m<=zmax;++m){
    usize+=sizeof(byte);
    if(mode[CFG_ZERNIKE][m]!=0) usize+=nph*nph*sizeof(float64_t); //mode[CFG_ZERNIKE][m] 
  }
// Karhune-Loeves
  usize+=2*sizeof(int32_t); // kmin kmax
  for(int m=kmin;m<=kmax;++m){
    usize+=sizeof(byte);
    if(mode[CFG_KARHUNEN_LOEVE][m]!=0) usize+=nph*nph*sizeof(float64_t); //mode[CFG_KARHUNEN_LOEVE][m] 
  }
// mode links
  usize+=2*sizeof(int32_t); // mmin mmax
  for(int m=mmin;m<=mmax;++m){
    usize+=sizeof(byte);    // is_used
    if(mode[0][m]!=0) usize+=sizeof(byte)+sizeof(int32_t); // type,index
  }
// covariance arrays
  usize+=sizeof(byte);
  if(covar[0]!=0) usize+=(mmax-mmin+1)*sizeof(float64_t); // covar[0]
  usize+=sizeof(byte);
  if(covar[CFG_ZERNIKE]!=0) usize+=(zmax-zmin+1)*sizeof(float64_t); // covar[CFG_ZERNIKE]
  usize+=sizeof(byte);
  if(covar[CFG_KARHUNEN_LOEVE]!=0) usize+=(kmax-kmin+1)*sizeof(float64_t); // covar[CFG_KARHUNEN_LOEVE]
// pupil
  usize+=sizeof(byte);
  if(pupil!=0) usize+=nph*nph*sizeof(float64_t);
// all done (I hope)
  size=usize;
  byte *p=new byte [size];
  int offs=0;
// pack p
  offs+=pack(p+offs,area,swap_endian);
  offs+=pack(p+offs,nph,swap_endian);
// Zernikes
  offs+=pack(p+offs,zmin,swap_endian);
  offs+=pack(p+offs,zmax,swap_endian);
  for(int m=zmin;m<=zmax;++m){
    byte is_used=(mode[CFG_ZERNIKE][m]!=0);
    offs+=pack(p+offs,is_used);
    if(is_used) offs+=pack(p+offs,mode[CFG_ZERNIKE][m],1,nph,1,nph,swap_endian);
  }
// Karhune-Loeves
  offs+=pack(p+offs,kmin,swap_endian);
  offs+=pack(p+offs,kmax,swap_endian);
  for(int m=kmin;m<=kmax;++m){
    byte is_used=(mode[CFG_KARHUNEN_LOEVE][m]!=0);
    offs+=pack(p+offs,is_used);
    if(is_used) offs+=pack(p+offs,mode[CFG_KARHUNEN_LOEVE][m],1,nph,1,nph,swap_endian);
  }
// link info
  offs+=pack(p+offs,mmin,swap_endian);
  offs+=pack(p+offs,mmax,swap_endian);
  for(int m=mmin;m<=mmax;++m){
    byte is_used=(mode[0][m]!=0);
    offs+=pack(p+offs,is_used);
    if(is_used){
      byte type=0;
      int index=-1;
//
      for(int n=zmin;n<=zmax;++n)
        if(mode[0][m]==mode[CFG_ZERNIKE][n]){
          index=n;
          type=CFG_ZERNIKE;
          break;
        }
      if(index<0)
        for(int n=kmin;n<=kmax;++n)
          if(mode[0][m]==mode[CFG_KARHUNEN_LOEVE][n]){
            index=n;
            type=CFG_KARHUNEN_LOEVE;
            break;
          }
//
      if(index>=0){ // found!
        offs+=pack(p+offs,type);
        offs+=pack(p+offs,index,swap_endian);
      }else{ // custom mode? Send the whole thing...?
        offs+=pack(p+offs,type);
        offs+=pack(p+offs,index,swap_endian);
      }
    } // not used...
  }
//
  byte has_covar=(covar[0]!=0);
  offs+=pack(p+offs,has_covar);
  if(has_covar) offs+=pack(p+offs,covar[0],mmin,mmax,swap_endian);
  has_covar=(covar[CFG_ZERNIKE]!=0);
  offs+=pack(p+offs,has_covar);
  if(has_covar) offs+=pack(p+offs,covar[CFG_ZERNIKE],zmin,zmax,swap_endian);
  has_covar=(covar[CFG_KARHUNEN_LOEVE]!=0);
  offs+=pack(p+offs,has_covar);
  if(has_covar) offs+=pack(p+offs,covar[CFG_KARHUNEN_LOEVE],kmin,kmax,swap_endian);
//
  byte has_pupil=(pupil!=0);
  offs+=pack(p+offs,has_pupil);
  if(has_pupil) offs+=pack(p+offs,pupil,1,nph,1,nph,swap_endian);
//
  if(offs!=usize) io.msg(MFBD_WARN,"modes::compress: data size mismatch; packed %d bytes but buffer was %d bytes!\n",offs,usize);
// end pack
  byte *buf=z_compress(p,size,level,swap_endian,io);
  delete[] p;
  return buf;
}

