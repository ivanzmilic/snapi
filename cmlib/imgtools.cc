#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "types.h"
#include "mem.h"
#include "const.h"

#include "imgtools.h"

fp_t linipol(fp_t **sub_img,int xl,int xh,int yl,int yh,int ix,int iy,fp_t x,fp_t y)
{
  if((ix<xl)||(iy<yl)||(ix>xh-1)||(iy>yh-1)) return 0.0; // out of bounds
  return (x*sub_img[ix+1][iy+1]+(1.0-x)*sub_img[ix][iy+1])*y+(x*sub_img[ix+1][iy]+(1.0-x)*sub_img[ix][iy])*(1.0-y);
}

fp_t bcsipol(fp_t **sub_img,int xl,int xh,int yl,int yh,int ix,int iy,fp_t x,fp_t y)
{
  if((ix<xl)||(iy<yl)||(ix>xh-1)||(iy>yh-1)) return 0.0; // out of bounds
  fp_t im[3][3][16][16]=BCSPLMAT;
  int xb=1+(ix==xh-1)-(ix==xl),yb=1+(iy==yh-1)-(iy==yl);
  int xx=ix+1-xb,yy=iy+1-yb;
  fp_t gx[8]={sub_img[xx+1][yy-1]-sub_img[xx-1][yy-1],sub_img[xx+1][yy  ]-sub_img[xx-1][yy  ],
              sub_img[xx+1][yy+1]-sub_img[xx-1][yy+1],sub_img[xx+1][yy+2]-sub_img[xx-1][yy+2],
              sub_img[xx+2][yy-1]-sub_img[xx  ][yy-1],sub_img[xx+2][yy  ]-sub_img[xx  ][yy  ],
              sub_img[xx+2][yy+1]-sub_img[xx  ][yy+1],sub_img[xx+2][yy+2]-sub_img[xx  ][yy+2]};
  fp_t gy[4]={sub_img[xx  ][yy+1]-sub_img[xx  ][yy-1],sub_img[xx+1][yy+1]-sub_img[xx+1][yy-1],
              sub_img[xx  ][yy+2]-sub_img[xx  ][yy  ],sub_img[xx+1][yy+2]-sub_img[xx+1][yy  ]};
  fp_t v[16]={sub_img[xx  ][yy  ],sub_img[xx+1][yy  ],sub_img[xx  ][yy+1],sub_img[xx+1][yy+1],
              gx[1],gx[2],gx[5],gx[6],gy[0],gy[1],gy[2],gy[3],
              gx[2]-gx[0],gx[3]-gx[1],gx[6]-gx[4],gx[7]-gx[5]};
  fp_t cfs[16];
  for(int i=0;i<=15;++i){
    cfs[i]=0.0;
    for(int j=0;j<=15;++j) cfs[i]+=im[xb][yb][i][j]*v[j];
  }
  return (cfs[0]+x*(cfs[1]+x*(cfs[2]+x*cfs[3]))+
          y*(cfs[4]+x*(cfs[5]+x*(cfs[6]+x*cfs[7]))+
          y*(cfs[8]+x*(cfs[9]+x*(cfs[10]+x*cfs[11]))+
          y*(cfs[12]+x*(cfs[13]+x*(cfs[14]+x*cfs[15]))))))/16.0;
}

fp_t **rotate(fp_t **sub_img,int xl,int xh,int yl,int yh,fp_t angle)
{
  fp_t **sim=ft2dim(xl,xh,yl,yh);
  int xo=(xl+xh)/2,yo=(yl+yh)/2;
  fp_t sa=sin(-angle),ca=cos(-angle);
  for(int ix=xl;ix<=xh;++ix)
    for(int iy=yl;iy<=yh;++iy){
      fp_t xx=(fp_t)(ix-xo)*ca-(fp_t)(iy-yo)*sa+(fp_t)xo;
      fp_t yy=(fp_t)(ix-xo)*sa+(fp_t)(iy-yo)*ca+(fp_t)yo;
      int x=(int)xx-(xx<0);
      int y=(int)yy-(yy<0);
      fp_t dx=xx-(fp_t)x;
      fp_t dy=yy-(fp_t)y;
      sim[ix][iy]=bcsipol(sub_img,xl,xh,yl,yh,x,y,dx,dy);
    }
  return sim;
}

fp_t **resample(fp_t **img,int xli,int xhi,int yli,int yhi,fp_t range,int xlo,int xho,int ylo,int yho)
{
  fp_t **result=ft2dim(xlo,xho,ylo,yho);
  memset(result[xlo]+ylo,0,(xho-xlo+1)*(yho-ylo+1)*sizeof(fp_t));
  fp_t o_range=0.5*(fp_t)(xhi-xli+1);
  int xio=(xli+xhi)/2,yio=(yli+yhi)/2;
  int xoo=(xlo+xho)/2,yoo=(ylo+yho)/2;
  for(int xo=xlo,xi=0;xo<=xho;++xo){
    fp_t x=(fp_t)(xo-xoo)/range;
    while(x>(fp_t)(xi-xio+1)/o_range) ++xi;
    fp_t dx=x*o_range-(fp_t)(xi-xio);
    for(int yo=ylo,yi=0;yo<=yho;++yo){
      fp_t y=(fp_t)(yo-yoo)/range;
      while(y>(fp_t)(yi-yio+1)/o_range) ++yi;
      fp_t dy=y*o_range-(fp_t)(yi-yio);
      result[xo][yo]=bcsipol(img,xli,xhi,yli,yhi,xi,yi,dx,dy);
    }
  }
  return result;
}

fp_t **rebin(fp_t **img,int xli,int xhi,int yli,int yhi,fp_t range,int xlo,int xho,int ylo,int yho)
{
  fp_t **result=ft2dim(xlo,xho,ylo,yho);
  memset(result[xlo]+ylo,0,(xho-xlo+1)*(yho-ylo+1)*sizeof(fp_t));
  fp_t **weight=ft2dim(xlo,xho,ylo,yho);
  memset(weight[xlo]+ylo,0,(xho-xlo+1)*(yho-ylo+1)*sizeof(fp_t));
  fp_t o_range=0.5*(fp_t)(xhi-xli+1);
//
  int xio=1+(xli+xhi)/2,yio=1+(yli+yhi)/2;
  int xoo=1+(xlo+xho)/2,yoo=1+(ylo+yho)/2;
//
  for(int xi=xli;xi<=xhi;++xi){
    fp_t xil=range*((fp_t)(xi-xio)-0.5)/o_range+(fp_t)xoo; // location in output pixels
    fp_t xih=range*((fp_t)(xi-xio)+0.5)/o_range+(fp_t)xoo;
    if((int)xil<(int)xih){ // split contribution
      fp_t fx=(xih-(fp_t)((int)xih))/(xih-xil);
      for(int yi=yli;yi<=yhi;++yi){
        fp_t yil=range*((fp_t)(yi-yio)-0.5)/o_range+(fp_t)yoo; // location in output pixels
        fp_t yih=range*((fp_t)(yi-yio)+0.5)/o_range+(fp_t)yoo;
        if((int)yil<(int)yih){ // split contribution
          fp_t fy=(yih-(fp_t)((int)yih))/(yih-yil);
          result[(int)xih][(int)yih]+=fx*fy*img[xi][yi];
          weight[(int)xih][(int)yih]+=fx*fy;
          result[(int)xil][(int)yih]+=(1.0-fx)*fy*img[xi][yi];
          weight[(int)xil][(int)yih]+=(1.0-fx)*fy;
          result[(int)xih][(int)yil]+=fx*(1.0-fy)*img[xi][yi];
          weight[(int)xih][(int)yil]+=fx*(1.0-fy);
          result[(int)xil][(int)yil]+=(1.0-fx)*(1.0-fy)*img[xi][yi];
          weight[(int)xil][(int)yil]+=(1.0-fx)*(1.0-fy);
        }else{                 // all in one pixel
          result[(int)xih][(int)yih]+=fx*img[xi][yi];
          weight[(int)xih][(int)yih]+=fx;
          result[(int)xil][(int)yih]+=(1.0-fx)*img[xi][yi];
          weight[(int)xil][(int)yih]+=(1.0-fx);
        }
      }
    }else{                     // all in one pixel
      for(int yi=yli;yi<=yhi;++yi){
        fp_t yil=range*((fp_t)(yi-yio)-0.5)/o_range+(fp_t)yoo; // location in output pixels
        fp_t yih=range*((fp_t)(yi-yio)+0.5)/o_range+(fp_t)yoo;
        if((int)yil<(int)yih){ // split contribution
          fp_t fy=(yih-(fp_t)((int)yih))/(yih-yil);
          result[(int)xih][(int)yih]+=fy*img[xi][yi];
          weight[(int)xih][(int)yih]+=fy;
          result[(int)xih][(int)yil]+=(1.0-fy)*img[xi][yi];
          weight[(int)xih][(int)yil]+=(1.0-fy);
        }else{                 // all in one pixel
          result[(int)xih][(int)yih]+=img[xi][yi];
          weight[(int)xih][(int)yih]+=1.0;
        }
      }
    }
  }
  for(int xo=xlo;xo<=xho;++xo)
    for(int yo=ylo;yo<=yho;++yo) if(weight[xo][yo]) result[xo][yo]/=weight[xo][yo]; else result[xo][yo]=0.0;
  del_ft2dim(weight,xlo,xho,ylo,yho);
  return result;
}

fp_t **make_window(int nx,int ny,int sm)
{
  if(2*sm>nx) exit(fprintf(stderr,"in make_window: margins too large in x-direction\n"));
  if(2*sm>ny) exit(fprintf(stderr,"in make_window: margins too large in y-direction\n"));
  fp_t **window=ft2dim(1,nx,1,ny);
  for(int y=1;y<=sm;++y) window[1][y]=(window[1][ny-y+1]=0.5*(1.0-cos(2.0*pi*(fp_t)(y-1)/(fp_t)(2*sm-1))));
  for(int y=sm+1;y<=ny-sm;++y) window[1][y]=1.0;
  for(int x=nx/2+1;x<=nx;++x)
    for(int y=1;y<=ny;++y) window[x][y]=window[1][x]*window[1][y];
  for(int x=1;x<=nx/2;++x) memcpy(window[x]+1,window[nx-x+1]+1,ny*sizeof(fp_t));
  return window;
}

