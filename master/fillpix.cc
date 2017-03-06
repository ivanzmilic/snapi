#include <math.h>
#include <string.h>
#include "types.h"
#include "conf.h"
#include "mem.h"
#include "io.h"
#include "uts.h"
#include "fillpix.h"

#define TINY 1.0E-20;
int ludcmp(fp_t **a,int n,int *indx,fp_t &d,io_class &io)
{
  int i,imax,j,k;
  fp_t big,dum,sum,temp;
  fp_t *vv=new fp_t [n]-1;
  
  d=1.0;
  for(i=1;i<=n;i++){
    big=0.0;
    for(j=1;j<=n;j++)
      if((temp=fabs(a[i][j]))>big) big=temp;
    if(big==0.0){
      io.msg(IOL_ERROR,"Singular matrix in routine ludcmp\n");
      return -1;
    }
    vv[i]=1.0/big;
  }
  for(j=1;j<=n;j++){
    for(i=1;i<j;i++){
      sum=a[i][j];
      for(k=1;k<i;k++) sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<=n;i++){
      sum=a[i][j];
      for(k=1;k<j;k++)
        sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum))>=big){
        big=dum;
        imax=i;
      }
    }
    if(j!=imax){
      for(k=1;k<=n;k++){
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      d=-d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if(a[j][j]==0.0) a[j][j]=TINY;
    if(j!=n){
      dum=1.0/(a[j][j]);
      for(i=j+1;i<=n;i++) a[i][j]*=dum;
    }
  }
  delete[] (vv+1);
  return 0;
}
#undef TINY

void lubksb(fp_t **a,int n,int *indx,fp_t *b)
{
  int i,ii=0,ip,j;
  fp_t sum;

  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

void bad_region(fp_t **p,int &xl,int &xh,int &yl,int &yh,int x_s,int y_s)
{ // grow rectangle around bad area
  int xxl=x_s-(x_s>xl),xxh=x_s+(x_s<xh),yyl=y_s-(y_s>yl),yyh=y_s+(y_s<yh);
  int xl_bad=1,xh_bad=1,yl_bad=1,yh_bad=1;
  do{
    if(xl_bad){
      xl_bad=0;
      for(int y=yyl;y<=yyh;++y) xl_bad+=(p[xxl][y]==0);
      if(xl_bad) if(xxl<=xl) xl_bad=0; else --xxl;
    }      
    if(xh_bad){
      xh_bad=0;
      for(int y=yyl;y<=yyh;++y) xh_bad+=(p[xxh][y]==0);
      if(xh_bad) if(xxh>=xh) xh_bad=0; else ++xxh;
    }      
    if(yl_bad){
      yl_bad=0;
      for(int x=xxl;x<=xxh;++x) yl_bad+=(p[x][yyl]==0);
      if(yl_bad) if(yyl<=yl) yl_bad=0; else --yyl;
    }      
    if(yh_bad){
      yh_bad=0;
      for(int x=xxl;x<=xxh;++x) yh_bad+=(p[x][yyh]==0);
      if(yh_bad) if(yyh>=yh) yh_bad=0; else ++yyh;
    }      
  }while(xl_bad||xh_bad||yl_bad||yh_bad);
  xl=xxl;
  xh=xxh;
  yl=yyl;
  yh=yyh;
}

int poly_fill_region(fp_t **a,int dx,int dy,int xl,int xh,int yl,int yh,io_class &io)
{
  fp_t *d=new fp_t [xh-xl+1]-xl,*xx=new fp_t [2*dx+2];
  fp_t **x=ft2dim(xl,xh,0,dx);
  memset(xx,0,(2*dx+2)*sizeof(fp_t));
  for(int i=xl;i<=xh;++i) d[i]=(fp_t)(i-xl)/(fp_t)(xh-xl)-0.5; // -0.5 .. 0.5
  for(int i=xl;i<=xh;++i){
    xx[0]+=(x[i][0]=1.0);
    for(int j=1;j<=dx;++j) xx[j]+=(x[i][j]=x[i][j-1]*d[i]);
    for(int j=dx+1;j<=2*dx;++j) xx[j]+=x[i][j-dx]*x[i][dx];
  }
  delete[] (d+xl);
  d=new fp_t [yh-yl+1]-yl;
  fp_t *yy=new fp_t [2*dy+2];
  fp_t **y=ft2dim(yl,yh,0,dy);
  memset(yy,0,(2*dy+2)*sizeof(fp_t));
  for(int i=yl;i<=yh;++i) d[i]=(fp_t)(i-yl)/(fp_t)(yh-yl)-0.5; // -0.5 .. 0.5
  for(int i=yl;i<=yh;++i){
    yy[0]+=(y[i][0]=1.0);
    for(int j=1;j<=dy;++j) yy[j]+=(y[i][j]=y[i][j-1]*d[i]);
    for(int j=dy+1;j<=2*dy;++j) yy[j]+=y[i][j-dy]*y[i][dy];
  }
  delete[] (d+yl);
//
  int N=(dx+1)*(dy+1);
  fp_t **alpha=ft2dim(1,N,1,N);
  fp_t *beta=new fp_t [N]-1;
  fp_t *row=new fp_t [yh-yl+1]-yl;
//
  int inc1=1;
  for(int i=0;i<=dx;++i)
    for(int j=0;j<=dy;++j){
      if(!j)
        for(int iy=yl;iy<=yh;++iy){
          row[iy]=0.0;
          for(int ix=xl;ix<=xh;++ix) if(a[ix][iy]) row[iy]+=x[ix][i]*a[ix][iy];
        }
      beta[inc1]=0.0;
      for(int iy=yl;iy<=yh;++iy) beta[inc1]+=row[iy]*y[iy][j];
      int inc2=1;
      for(int k=0;k<=dx;++k)
        for(int l=0;l<=dy;++l)
 	  alpha[inc1][inc2++]=xx[i+k]*yy[j+l];	  // Matrix to be LU decomposed.
      ++inc1;
    }
  del_ft2dim(x,xl,xh,0,dx);
  del_ft2dim(y,yl,yh,0,dy);
  delete[] (row+yl);
  delete[] xx;
  delete[] yy;
  fp_t det;
  int *idx=new int [N]-1;
  if(ludcmp(alpha,N,idx,det,io)<0) return -1;      // Replaces alfa with its LU decompose
  lubksb(alpha,N,idx,beta);                        // Solves the set of linear equations.
  delete[] (idx+1);
  del_ft2dim(alpha,1,N,1,N);
  fp_t **cfs=ft2dim(beta+1,0,dx,0,dy);
// fill in the bad pixels
  for(int ix=xl;ix<=xh;++ix){
    fp_t lx=(fp_t)(ix-xl)/(fp_t)(xh-xl)-0.5;
    for(int iy=yl;iy<=yh;++iy){
      fp_t ly=(fp_t)(iy-yl)/(fp_t)(yh-yl)-0.5;
      if(!a[ix][iy]){
        fp_t xsum=0.0;
        for(int jx=dx;jx>=0;--jx){
          xsum*=lx;
          fp_t ysum=0.0;
          for(int jy=dy;jy>=0;--jy){
            ysum*=ly;
            ysum+=cfs[jx][jy];
          }
          xsum+=ysum;
        }
        a[ix][iy]=(float32_t)xsum;
      }
    }
  }
//
  del_ft2dim(cfs,0,dx,0,dy);
  return 0;
}

fp_t median(fp_t **a,int xl,int xh,int yl,int yh)
{ // only edge needs to be tested, all other pixels are 0
  float32_t *b=new float32_t [2*(xh-xl+yh-yl)]-1;
  int n=0;
  for(int y=yl;y<=yh;++y){
    if(a[xl][y]) b[++n]=a[xl][y];
    if(a[xh][y]) b[++n]=a[xh][y];
  }
  for(int x=xl+1;x<=xh-1;++x){
    if(a[x][yl]) b[++n]=a[x][yl];
    if(a[x][yh]) b[++n]=a[x][yh];
  }
  if(!n){
    delete[] (b+1);
    return 0.0;
  }
// sort
  for(int i=1;i<=n-1;++i)
    for(int j=1;j<=n-i;++j)
      if(b[j]>b[j+1]) swap(b[j],b[j+1]);
//
  fp_t rv=(n%2)?b[n/2+1]:0.5*(b[n/2]+b[n/2+1]);
  delete[] (b+1);
  return rv;
}

int median_fill_region(fp_t **a,int xl,int xh,int yl,int yh,io_class &io)
{
  fp_t **b=ft2dim(xl,xh,yl,yh);
  memset(b[xl]+yl,0,(xh-xl+1)*(yh-yl+1)*sizeof(float32_t));
  for(int x=xl;x<=xh;++x)
    for(int y=yl;y<=yh;++y)
      if(a[x][y])
        b[x][y]=a[x][y];
      else
        for(int r=1;!b[x][y];++r) b[x][y]=median(a,max(x-r,xl),min(x+r,xh),max(y-r,yl),min(y+r,yh));
  for(int x=xl;x<=xh;++x) memcpy(a[x]+yl,b[x]+yl,(yh-yl+1)*sizeof(fp_t));
  del_ft2dim(b,xl,xh,yl,yh);
  return 0;
}

#define beta 4.0
#define deltasqr 4.0

fp_t inv_dist_wght(fp_t **a,int xl,int xh,int yl,int yh,int xb,int yb)
{
  fp_t weight=0.0,res=0.0;
  for(int x=xl;x<=xh;++x)
    for(int y=yl;y<=yh;++y)
      if(a[x][y]){
        fp_t c=pow(sqrt((fp_t)sqr(x-xb)+(fp_t)sqr(y-yb)+deltasqr),-beta);
        res+=c*a[x][y];
        weight+=c;
      }
  return res/weight;
}

#undef beta
#undef deltasqr

fp_t **fillpix(fp_t **pic,int xl,int xh,int yl,int yh,byte method,io_class &io)
{
  int badpix=0;
  for(int x=xl;x<=xh;++x)
    for(int y=yl;y<=yh;++y) if(!pic[x][y]) ++badpix;
  if(!badpix) return pic;
  fp_t *bp=new fp_t [badpix]-1;
  memset(bp+1,0,badpix*sizeof(fp_t));
  for(int x=xl,n=0;x<=xh;++x)
    for(int y=yl;y<=yh;++y)
      if(!pic[x][y]){
        ++n;
        switch(method){
          case(CFG_FPM_MEDIAN):{
            for(int r=1;!bp[n];++r) bp[n]=median(pic,max(x-r,xl),min(x+r,xh),max(y-r,yl),min(y+r,yh));
            break;
          }
          case(CFG_FPM_INVDISTWEIGHT):{
            bp[n]=inv_dist_wght(pic,max(x-64,xl),min(x+64,xh),max(y-64,yl),min(y+64,yh),x,y);
            break;
          }
        }
      }
  for(int x=xl,n=0;x<=xh;++x)
    for(int y=yl;y<=yh;++y) if(!pic[x][y]) pic[x][y]=bp[++n];
  delete[] (bp+1);
  return pic;
}
