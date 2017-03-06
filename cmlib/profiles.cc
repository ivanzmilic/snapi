#include <math.h>
#include <stdlib.h>

#include "types.h"
#include "io.h"

void fvoigt(fp_t x,fp_t a,fp_t &hh,fp_t &ff,fp_t &dh,fp_t &df)
{
  fp_t v=x;
  fp_t s=fabs(v)+a;
//
  complex_t t(a,-v);
  complex_t u(a*a-v*v,-2.0*a*v);
//
  complex_t w4;
//
  if(s>=15.0){                        // Region 1 
    w4=t*0.5641896/(0.5+u);
  }else if(s>5.5){                    // Region 2
    w4=t*(1.410474+u*0.5641896)/(0.75+u*(3.0+u));
  }else if(a>(0.195*fabs(v)-0.176)){  // Region 3
    w4=(16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*0.5642236))))/(16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))));
  }else{                              // Region 4: a v<(0.195*fabs(v)-0.176) v
    w4=t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419))))));
    complex_t v4=(32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))));
    fp_t su=sin(u.im);
    w4=exp(u.re)*complex_t(sqrt(1.0-su*su),su)-w4/v4; // cos=sqrt(1-sin^2) was used
//    w4=exp(u.re)*complex_t(cos(u.re)),sin(u.re)))-w4/v4; // cos=sqrt(1-sin^2) was used
  }
  hh=w4.re/1.77245386399759;
  ff=w4.im/1.77245386399759;
//    Armstrong & Nicholls (1972): dw/dz = -2*z*w(z)+ 2i/sqrt(pi)
  complex_t z(v,a);
  complex_t dw=-2.0*z*w4;
  dh=dw.re/1.77245386399759;
  df=(dw.im+1.1283791671)/1.77245386399759;
}

void fvoigtn(fp_t x,fp_t a,fp_t &hh,fp_t &ff,fp_t &dh,fp_t &df)
{
  fp_t v=x;
  fp_t s=fabs(v)+a;
  fp_t ds = fabs(v)/v + a;
//
  complex_t t(a,-v);
  complex_t dt(0,-1.0);
  complex_t u(a*a-v*v,-2.0*a*v);
  complex_t du(-2.0*v,-2.0*a);
//
  complex_t w4,v4;
  complex_t dw4,dv4;
//
  if(s>=15.0){                        // Region 1
    //printf("Region 1!\n"); 
    w4=t*0.5641896/(0.5+u);
    dw4 = dt*0.5641896/(0.5+u) - t*0.5641896/(0.5+u)/(0.5+u)*du;
  }else if(s>5.5){                    // Region 2
    //printf("Region 2!\n"); 
    w4=t*(1.410474+u*0.5641896)/(0.75+u*(3.0+u));
    dw4 = dt*(1.410474+u*0.5641896)/(0.75+u*(3.0+u)) + t*du*0.5641896/(0.75+u*(3.0+u));
    dw4 = dw4 - t*du*(1.410474+u*0.5641896)/(0.75+u*(3.0+u))/(0.75+u*(3.0+u)) * (3.0+2.0*u);
  }else if(a>(0.195*fabs(v)-0.176)){  // Region 3
    //printf("Region 3!\n"); 
    w4=(16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*0.5642236))))/(16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))));
    complex_t upper = (16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*0.5642236))));
    complex_t lower = (16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))));
    complex_t der_upper = 20.20933 + 11.96482*2.0*t + 3.778987*3.0*t*t + 0.5642236*4.0*t*t*t;
    complex_t der_lower = 38.82363 + 39.27121*2.0*t + 21.69274*3.0*t*t + 6.699398*4.0*t*t*t + 5.0*t*t*t*t;
    dw4 = dt*(der_upper*lower - der_lower*upper)/lower/lower;
  }else{                              // Region 4: a v<(0.195*fabs(v)-0.176) v
    //printf("Region 4!\n"); 
    
    // Numerical verification:
    /*
    v += 0.001;
    complex_t t_up(a,-v);
    complex_t u_up(a*a-v*v,-2.0*a*v);
    t = t_up; u=u_up;
    w4=t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419))))));
    complex_t dw4_num = w4;
    v4=(32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))));
    complex_t dv4_num = v4;
    v -= 0.002;
    complex_t t_down(a,-v);
    complex_t u_down(a*a-v*v,-2.0*a*v);
    t = t_down; u=u_down;
    w4=t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419))))));
    dw4_num = dw4_num - w4;
    v4=(32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))));
    dv4_num = dv4_num - v4;
    dw4_num = dw4_num / 0.002;
    dv4_num = dv4_num / 0.002;
    complex_t dtn = (t_up - t_down)/0.002;
    complex_t dun = (u_up - u_down)/0.002;
    printf("Numerical:\n");
    printf("dt = %e %e \n", dtn.re, dtn.im);
    printf("du = %e %e \n", dun.re, dun.im);
    printf("dw4 = %e %e \n", dw4_num.re, dw4_num.im);
    printf("dv4 = %e %e \n", dv4_num.re, dv4_num.im);
    v+=0.001;
    complex_t t_final(a,-v);
    complex_t u_final(a*a-v*v,-2.0*a*v);
    t = t_final; u=u_final;
    */

    w4=t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419))))));
    dw4 = dt*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419))))));
    dw4 = dw4 - t*du*(3321.9905-u*(1540.787*2.0-u*(219.0313*3.0-u*(35.76683*4.0-u*(1.320522*5.0-u*0.56419*6.0)))));
    /*printf("Analytical:\n");
    printf("dt = %e %e \n", dt.re, dt.im);
    printf("du = %e %e \n", du.re, du.im);
    printf("dw4 = %e %e \n", dw4.re, dw4.im);*/
    

    v4=(32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))));

    dv4 = -1.0*du*(24322.84-u*(9022.228*2.0-u*(2186.181*3.0-u*(364.2191*4.0-u*(61.57037*5.0-u*(1.841439*6.0-u*7.0))))));

    //printf("dv4 = %e %e \n", dv4.re, dv4.im); 
    
    fp_t su=sin(u.im);
    fp_t dsu = cos(u.im) * du.im;

    complex_t dw4_full = exp(u.re)*du.re * complex_t(sqrt(1.0-su*su),su) + exp(u.re)*complex_t(-su*du.im,dsu);
    //complex_t dw4_full = exp(u.re)*du.re * complex_t(sqrt(1.0-su*su),su) + exp(u.re) * complex_t(-su*dsu/sqrt(1.0-su*su),dsu);
    dw4_full = dw4_full - (dw4*v4 - dv4*w4)/v4/v4;
    dw4=dw4_full;
    w4=exp(u.re)*complex_t(sqrt(1.0-su*su),su)-w4/v4; // cos=sqrt(1-sin^2) was used
//    w4=exp(u.re)*complex_t(cos(u.re)),sin(u.re)))-w4/v4; // cos=sqrt(1-sin^2) was used
  }
  hh=w4.re/1.77245386399759;
  ff=w4.im/1.77245386399759;
  dh = dw4.re/1.77245386399759;
  df = dw4.im/1.77245386399759;
  //printf("%e %e %e %e %e \n", x, hh,ff,dh,df);
//    Armstrong & Nicholls (1972): dw/dz = -2*z*w(z)+ 2i/sqrt(pi)
  /*complex_t z(v,a);
  complex_t dw=-2.0*z*w4;
  dh=dw.re/1.77245386399759;
  df=(dw.im+1.1283791671)/1.77245386399759;
  printf("%e %e \n", dh,df);*/
}



fp_t fvoigt(fp_t v,fp_t a)
{
  fp_t s=fabs(v)+a;
//
  complex_t t(a,-v);
  complex_t u(a*a-v*v,-2.0*a*v);
//
  complex_t w4;
//
  if(s>=15.0){                        // Region 1 
    w4=t*0.5641896/(0.5+u);
  }else if(s>5.5){                    // Region 2
    w4=t*(1.410474+u*0.5641896)/(0.75+u*(3.0+u));
  }else if(a>(0.195*fabs(v)-0.176)){  // Region 3
    w4=(16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*0.5642236))))/(16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))));
  }else{                              // Region 4: a v<(0.195*fabs(v)-0.176) v
    w4=t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419))))));
    complex_t v4=(32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))));
    fp_t su=sin(u.im);
    w4=exp(u.re)*complex_t(sqrt(1.0-su*su),su)-w4/v4; // cos=sqrt(1-sin^2) was used
//    w4=exp(u.re)*complex_t(cos(u.re)),sin(u.re)))-w4/v4; // cos=sqrt(1-sin^2) was used
  }
  return w4.re/1.77245386399759;
  
}

void voigt_num_der(fp_t x, fp_t a,fp_t &hh,fp_t &ff,fp_t &dh,fp_t &df){

  fp_t x_up = x+0.001;
  fp_t H,F,DH,DF;
  fvoigt(x_up, a, H,F,DH,DF);
  dh = H;
  df = F;
  fp_t x_down = x-0.001;
  fvoigt(x_down, a, H,F,DH,DF);
  dh -= H;
  df -= F;
  dh /= 0.002;
  df /= 0.002;
}
