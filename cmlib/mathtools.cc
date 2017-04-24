#include <math.h>
#include <string.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <float.h>

#include "types.h"
#include "io.h"
#include "uts.h"

#include "mathtools.h"
#include "mem.h"

#define TINY 1.0E-20

int ludcmp(fp_t **a,int ll,int ul,int *indx,fp_t &d,io_class &io)
{
  int imax;
  fp_t *vv=new fp_t [ul-ll+1]-ll;
  
  d=1.0;
  for(int i=ll;i<=ul;++i){
    fp_t big=0.0,temp;
    for(int j=ll;j<=ul;++j) if((temp=fabs(a[i][j]))>big) big=temp;
    if(big==0.0){
      io.msg(IOL_ERROR,"Singular matrix in routine ludcmp\n");
      return -1;
    }
    vv[i]=1.0/big;
  }
  for(int j=ll;j<=ul;++j){
    for(int i=ll;i<j;++i){
      fp_t sum=a[i][j];
      for(int k=ll;k<i;++k) sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    fp_t big=0.0;
    for(int i=j;i<=ul;++i){
      fp_t sum=a[i][j],dum;
      for(int k=1;k<j;++k) sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum))>=big){
        big=dum;
        imax=i;
      }
    }
    if(j!=imax){
      for(int k=ll;k<=ul;++k){
        fp_t dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      d=-d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if(a[j][j]==0.0) a[j][j]=TINY;
    if(j!=ul){
      fp_t dum=1.0/(a[j][j]);
      for(int i=j+1;i<=ul;i++) a[i][j]*=dum;
    }
  }
  delete[] (vv+ll);
  return 0;
}

#undef TINY

void lubksb(fp_t **a,int ll,int ul,int *indx,fp_t *b)
{
  int ii=ll-1;

  for(int i=ll;i<=ul;++i){
    int ip=indx[i];
    fp_t sum=b[ip];
    b[ip]=b[i];
    if(ii>=ll)
      for(int j=ii;j<=i-1;++j) sum-=a[i][j]*b[j];
    else if(sum) ii=i;
    b[i]=sum;
  }
  for(int i=ul;i>=ll;--i){
    fp_t sum=b[i];
    for(int j=i+1;j<=ul;++j) sum-=a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

void Shipley_Inversion(fp_t **MI,int ll,int ul)
{
  for(int k=ll;k<=ul;k++){
    MI[k][k]=-1.0/MI[k][k];         // the pivot element 
    for(int i=ll;i<=ul;++i) if(i!=k) MI[i][k]*=MI[k][k];//the pivot column 
    for(int i=ll;i<=ul;++i)           //elements not in a pivot row or column
      if(i!=k)
        for(int j=ll;j<=ul;++j)
          if(j!=k)
            MI[i][j]+=MI[i][k]*MI[k][j];
    for(int i=ll;i<=ul;++i)           //elements in a pivot row
      if(i!=k)
        MI[k][i]*=MI[k][k];
  }
  for(int i=ll;i<=ul;++i)
    for(int j=ll;j<=ul;++j) MI[i][j]=-MI[i][j];
}

void invert(fp_t **M,fp_t **MI,int m,int n)
{
  if(m==n){ // square
    memcpy(MI[1]+1,M[1]+1,n*m*sizeof(fp_t));
/*
    memset(M[1]+1,0,n*m*sizeof(fp_t));
    for(int i=1;i<=m;++i) 
      for(int j=1;j<=m;++j)
        for(int k=1;k<=m;++k) M[i][k]+=MI[i][j]*MI[j][k];
    gaussj(MI,n,M,m);
*/
    Shipley_Inversion(MI,1,m);
  }else{ // non-square
    fprintf(stderr,"error: non square matrix inversion not implemented!\n");
  }
}

fp_t pythag(fp_t a, fp_t b)
{
  fp_t absa,absb;
  //return sqrt(a*a + b*b);
  absa=fabs(a);
  absb=fabs(b);
  if(absa>absb)
    return absa*sqrt(1.0+sqr(absb/absa));
  else
    return(absb==0.0)?0.0:absb*sqrt(1.0+sqr(absa/absb));
}

int svdcmp(fp_t **a, int m, int n, fp_t *w, fp_t **v)
{
  fp_t g,scale,anorm;
  fp_t *rv1=new fp_t [n]-1;
  g=scale=anorm=0.0;
  for(int i=1;i<=n;++i){
    int l=i+1;
    rv1[i]=scale*g;
    g=scale=0.0;
    if(i<=m){
      for(int k=i;k<=m;++k) scale+=fabs(a[k][i]);
      if(scale){
        fp_t s=0.0;
        for(int k=i;k<=m;++k){
          a[k][i]/=scale;
          s+=a[k][i]*a[k][i];
        }
        fp_t f=a[i][i];
        g=-sign(sqrt(s),f);
        fp_t h=f*g-s;
        a[i][i]=f-g;
        for(int j=l;j<=n;++j){
          fp_t sum=0.0;
          for(int k=i;k<=m;++k) sum+=a[k][i]*a[k][j];
          fp_t fct=sum/h;
          for(int k=i;k<=m;++k) a[k][j]+=fct*a[k][i];
        }
        for(int k=i;k<=m;++k) a[k][i]*=scale;
      }
    }
    w[i]=scale*g;
    g=scale=0.0;
    if((i<=m)&&(i!=n)){
      for(int k=l;k<=n;++k) scale+=fabs(a[i][k]);
      if(scale){
        fp_t s=0.0;
        for(int k=l;k<=n;++k){
          a[i][k]/=scale;
          s+=a[i][k]*a[i][k];
        }
        fp_t f=a[i][l];
        g=-sign(sqrt(s),f);
        fp_t h=f*g-s;
        a[i][l]=f-g;
        for(int k=l;k<=n;++k) rv1[k]=a[i][k]/h;
        for(int j=l;j<=m;++j){
          fp_t sum=0.0;
          for(int k=l;k<=n;++k) sum+=a[j][k]*a[i][k];
          for(int k=l;k<=n;++k) a[j][k]+=sum*rv1[k];
        }
        for(int k=l;k<=n;++k) a[i][k]*=scale;
      }
    }
    anorm=max(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  {
    fp_t f;
    for(int i=n,l;i>=1;--i){
      if(i<n){
        if(f){
          for(int j=l;j<=n;++j) v[j][i]=(a[i][j]/a[i][l])/f;
          for(int j=l;j<=n;++j){
            fp_t sum=0.0;
            for(int k=l;k<=n;++k) sum+=a[i][k]*v[k][j];
            for(int k=l;k<=n;++k) v[k][j]+=sum*v[k][i];
          }
        }
        for(int j=l;j<=n;++j) v[i][j]=v[j][i]=0.0;
      }
      v[i][i]=1.0;
      f=rv1[i];
      l=i;
    }
  }
  for(int i=min(m,n);i>=1;--i){
    int l=i+1;
    g=w[i];
    for(int j=l;j<=n;++j) a[i][j]=0.0;
    if(g){
      g=1.0/g;
      for(int j=l;j<=n;++j){
        fp_t sum=0.0;
        for(int k=l;k<=m;++k) sum+=a[k][i]*a[k][j];
        fp_t f=(sum/a[i][i])*g;
        for(int k=i;k<=m;++k) a[k][j]+=f*a[k][i];
      }
      for(int j=i;j<=m;++j) a[j][i]*=g;
    }else for(int j=i;j<=m;++j) a[j][i]=0.0;
    ++a[i][i];
  }
  for(int k=n;k>=1;--k){
    for(int its=1;its<=1000;++its){
      int flag=1,nm,l;
      for(l=k;l>=1;--l){
        nm=l-1;
        if((fp_t)(fabs(rv1[l])+anorm)==anorm){
          flag=0;
          break;
        }
        if((fp_t)(fabs(w[nm])+anorm)==anorm) break;
      }
      if(flag){
        fp_t c=0.0,s=1.0;
        for(int i=l;i<=k;++i){
          fp_t f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if((fp_t)(fabs(f)+anorm)==anorm) break;
          g=w[i];
          fp_t h=pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s=-f*h;
          for(int j=1;j<=m;++j){
            fp_t y=a[j][nm];
            fp_t z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      fp_t z=w[k];
      if(l==k){
        if(z<0.0){
          w[k]=-z;
          for(int j=1;j<=n;++j) v[j][k]=-v[j][k];
        }
        break;
      }
      if(its==00) return -1;
      fp_t x=w[l];
      fp_t y=w[nm=k-1];
      g=rv1[nm];
      fp_t h=rv1[k];
      fp_t f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
      fp_t c=1.0,s=1.0;
      for(int j=l;j<=nm;++j){
        int i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g*=c;
        z=pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y*=c;
        for(int jj=1;jj<=n;++jj){
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j]=z;
        if(z){
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for(int jj=1;jj<=m;++jj){
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  delete[] (rv1+1);
  return 0;
}

int svd_analyze(fp_t * M, int start, int N){

  // Analyze square matrix with svd to spot bad conditioning

  fp_t ** M_to_svd = ft2dim(1,N,1,N);
  memcpy(M_to_svd[1]+1, M+start, N*N*sizeof(fp_t));

  fp_t * w = new fp_t [N]-1; 
  fp_t ** V = ft2dim(1,N,1,N);

  if(svdcmp(M_to_svd,N,N,w,V)<0) printf("SVD: singular matrix?\n"); 

  for (int i=1;i<=N;++i)
    printf("%d %e \n", i, w[i]);

  del_ft2dim(V,1,N,1,N);
  del_ft2dim(M_to_svd,1,N,1,N);
  delete[](w+1);

  return 0;
}

fp_t ** solve_svd(fp_t ** Matrix, fp_t ** rhs, int start, int N, int N_equations){

  fp_t ** M_to_svd = ft2dim(1,N,1,N);
  memcpy(M_to_svd[1]+1, Matrix[start]+start, N*N*sizeof(fp_t));

  fp_t * w = new fp_t [N]-1; 
  fp_t ** V = ft2dim(1,N,1,N);

  if(svdcmp(M_to_svd,N,N,w,V)<0) printf("SVD: singular matrix?\n"); 

  // Now we have the decomposition now it is time to use it to solve:

  fp_t crit = 1E13;

  fp_t wmax=w[1];
  for(int i=1;i<=N;++i){ 
    //printf("%d %e \n", i, w[i]);
    wmax=max(wmax,w[i]);
  }
  for(int i=1;i<=N;++i){
    if((w[i])&&((wmax/w[i])<=crit)){
      for(int j=1;j<=N;++j) V[j][i]/=w[i];
    }else{
      for(int j=1;j<=N;++j) V[j][i]=0.0;
      //printf("SVD: w[%d]: limit exceeded (%E<=%E*%E)!\n",i,wmax,w[i],crit);
    }
  }

  fp_t ** result = ft2dim(1,N_equations,1,N);

  for (int ie = 1; ie<=N_equations; ++ie){
    fp_t *t=new fp_t [N]-1;
    for(int i=1;i<=N;++i){
      t[i]=0.0;
      for(int j=1;j<=N;++j) t[i]+=M_to_svd[j][i]*rhs[ie][j];
    }
    for(int i=1;i<=N;++i){
      result[ie][i]=0.0;
      for(int j=1;j<=N;++j) result[ie][i]+=V[i][j]*t[j];
    }
    delete[] (t+1);
  }
  return result;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------

// More functions:

int multiply_4x4(fp_t ** A, fp_t ** B, fp_t ** result){

  // This is a trivial, for-the-debug function which multiplies matrices, no debug lines, no check no nothing. 

  for (int i = 1; i<=4; ++i)
    for (int j = 1; j<=4; ++j)
    {
      result[i][j] = 0.0;
      for (int k = 1; k<=4; ++k)
        result[i][j] += A[i][k] * B[k][j];
    }

  return 0;
}

int multiply_4x4(fp_t ** A, fp_t * B, fp_t * result){

  // This is a trivial, for-the-debug function which multiplies matrices, no debug lines, no check no nothing. 

  //printf("You ran me biatch!\n");

  for (int i = 1; i<=4; ++i)
    {
      result[i] = 0.0;
      for (int k = 1; k<=4; ++k)
        result[i] += A[i][k] * B[k];
    }

  return 0;
}

fp_t max_1d(fp_t * array, int begin, int end){
  fp_t max = array[begin];
  for (int i = begin+1; i<=end; ++i)
    if (array[i] > max) max = array[i];

  return max;
}

fp_t min_1d(fp_t * array, int begin, int end){
  fp_t min = array[begin];
  for (int i = begin+1; i<=end; ++i)
    if (array[i] < min) min = array[i];

  return min;
}

int  max_1d_index(fp_t * array, int begin, int end){
  int index = begin;
  fp_t max = array[begin];
  for (int i = begin+1; i<=end; ++i)
    if (array[i] > max){
      index = i;
      max = array[i];
    }

  return index;
}

int  max_1d_index(int * array, int begin, int end){
  int index = begin;
  int max = array[begin];
  for (int i = begin+1; i<=end; ++i)
    if (array[i] > max){
      index = i;
      max = array[i];
    }

  return index;
}



// Otheer versions of matrix inverse:

void Crout(int d,fp_t*S,fp_t*D){
   for(int k=0;k<d;++k){
      for(int i=k;i<d;++i){
         fp_t sum=0.;
         for(int p=0;p<k;++p)sum+=D[i*d+p]*D[p*d+k];
         D[i*d+k]=S[i*d+k]-sum; // not dividing by diagonals
      }
      for(int j=k+1;j<d;++j){
         fp_t sum=0.;
         for(int p=0;p<k;++p)sum+=D[k*d+p]*D[p*d+j];
         D[k*d+j]=(S[k*d+j]-sum)/D[k*d+k];
      }
   }
}

void solveCrout(int d,fp_t*LU,fp_t*b,fp_t*x){
   fp_t y[d];
   for(int i=0;i<d;++i){
      fp_t sum=0.;
      for(int k=0;k<i;++k)sum+=LU[i*d+k]*y[k];
      y[i]=(b[i]-sum)/LU[i*d+i];
   }
   for(int i=d-1;i>=0;--i){
      fp_t sum=0.;
      for(int k=i+1;k<d;++k)sum+=LU[i*d+k]*x[k];
      x[i]=(y[i]-sum); // not dividing by diagonals
   }
}

// Same as the one above but long double, i.e. __float128
void Crout(int d,__float128 *S,__float128 *D){
   for(int k=0;k<d;++k){
      for(int i=k;i<d;++i){
         fp_t sum=0.0;
         for(int p=0;p<k;++p)sum+=D[i*d+p]*D[p*d+k];
         D[i*d+k]=S[i*d+k]-sum; // not dividing by diagonals
      }
      for(int j=k+1;j<d;++j){
         fp_t sum=0.;
         for(int p=0;p<k;++p)sum+=D[k*d+p]*D[p*d+j];
         D[k*d+j]=(S[k*d+j]-sum)/D[k*d+k];
      }
   }
}

void solveCrout(int d,__float128 * LU, __float128 *b, __float128 *x){
   fp_t y[d];
   for(int i=0;i<d;++i){
      fp_t sum=0.;
      for(int k=0;k<i;++k)sum+=LU[i*d+k]*y[k];
      y[i]=(b[i]-sum)/LU[i*d+i];
   }
   for(int i=d-1;i>=0;--i){
      fp_t sum=0.;
      for(int k=i+1;k<d;++k)sum+=LU[i*d+k]*x[k];
      x[i]=(y[i]-sum); // not dividing by diagonals
   }
}




fp_t ** multiply_square(fp_t ** A, fp_t ** B, int dim){

  fp_t ** result = ft2dim(1, dim, 1, dim);
  for (int i = 1; i<=dim; ++i)
    for (int ii = 1; ii<=dim; ++ii){
      result[i][ii] = 0.0;
      for (int k = 1; k<=dim; ++k)
        result[i][ii] += A[i][k] * B[k][ii];
    }

  return result;
}

fp_t * multiply_vector(fp_t **A, fp_t * B, int dim){

  fp_t * result = new fp_t [dim] - 1;
  for (int i = 1; i<=dim; ++i){
    result [i] = 0.0;
    for (int ii = 1; ii<=dim; ++ii)
      result[i] += A[i][ii] * B[ii];
  }
  return result;
}

fp_t * multiply_vector(fp_t ** A, fp_t * B, int N_rows, int N_columns){

  fp_t * result = new fp_t [N_rows] - 1;
  memset(result+1,0,N_rows*sizeof(fp_t));
  for (int i=1;i<=N_rows;++i)
    for (int j=1;j<=N_columns;++j)
      result[i] += A[i][j] * B[j];

  return result;
}

fp_t * multiply_with_vector(fp_t * A, fp_t ** B, int dim){

  fp_t * result = new fp_t [dim] - 1;
  for (int i = 1; i<=dim; ++i){
    result [i] = 0.0;
    for (int ii = 1; ii<=dim; ++ii)
      result[i] += A[ii] * B[ii][i];
  }
  return result;

}

fp_t ** transpose(fp_t ** A, int N_rows, int N_columns){

  fp_t ** B = ft2dim(1, N_columns, 1, N_rows);
  for (int i=1;i<=N_rows;++i) for (int j=1;j<=N_columns;++j) B[j][i] = A[i][j];
  return B;
}

fp_t **** transpose(fp_t **** A, int N1, int N2, int N3, int N4){

  fp_t **** B = ft4dim(1,N4,1,N3,1,N2,1,N1);
  for (int i1=1;i1<=N1;++i1)
    for (int i2=1;i2<=N2;++i2)
      for (int i3=1;i3<=N3;++i3)
        for (int i4=1;i4<=N4;++i4)
          B[i4][i3][i2][i1] = A[i1][i2][i3][i4];
  del_ft4dim(A,1,N1,1,N2,1,N3,1,N4);
  return B;
}

fp_t ** multiply_with_transpose(fp_t ** A, int N_rows, int N_columns){

  fp_t ** B = ft2dim(1,N_columns,1,N_columns);
  memset(B[1]+1,0,N_columns*N_columns*sizeof(fp_t));
  for (int i=1;i<=N_columns;++i)
    for (int j=1;j<=N_columns;++j)
      for (int k=1;k<=N_rows;++k)
        B[i][j] += A[k][i] * A[k][j];

  return B;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fp_t * solve(fp_t ** A, fp_t * rhs, int from, int to){

  int size = to-from+1;
  fp_t * M_to_solve = A[from]+from;
  fp_t * M_LU = new fp_t [size*size];
  fp_t * solution = new fp_t [size];
  fp_t * b = rhs+from;  
  Crout(size,M_to_solve, M_LU);
  solveCrout(size,M_LU,b,solution);
  Crout(size,M_to_solve, M_LU);
  solveCrout(size,M_LU,b,solution);

  return solution-from;
}




// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fp_t aps(fp_t x){
  return (x>0) ? x : -x;
}

int Broyden_update (fp_t ** Broyden, fp_t * defect, fp_t * correction, int size){} // OBSOLETE

// Stolen implementation of E1 exponential integral, taken from: http://www.mymathlib.com/c_source/functions/exponential_integrals/exponential_integral_Ei.c

//                         Internally Defined Routines                        //



//                         Internally Defined Constants                       //
static const long double epsilon = 10.0 * LDBL_EPSILON;

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Ei( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Exponential_Integral_Ei( x );                                      //
////////////////////////////////////////////////////////////////////////////////
double Exponential_Integral_Ei( double x )
{
   return (double) xExponential_Integral_Ei( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xExponential_Integral_Ei( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the exponential integral Ei().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xExponential_Integral_Ei( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xExponential_Integral_Ei( long double x )
{
   if ( x < -5.0L ) return Continued_Fraction_Ei(x);
   if ( x == 0.0L ) return -DBL_MAX;
   if ( x < 6.8L )  return Power_Series_Ei(x);
   if ( x < 50.0L ) return Argument_Addition_Series_Ei(x);
   return Continued_Fraction_Ei(x);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Continued_Fraction_Ei( long double x )                  //
//                                                                            //
//  Description:                                                              //
//     For x < -5 or x > 50, the continued fraction representation of Ei      //
//     converges fairly rapidly.                                              //
//                                                                            //
//     The continued fraction expansion of Ei(x) is:                          //
//        Ei(x) = -exp(x) { 1/(-x+1-) 1/(-x+3-) 4/(-x+5-) 9/(-x+7-) ... }.    //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Continued_Fraction_Ei( long double x )
{
   long double Am1 = 1.0L;
   long double A0 = 0.0L;
   long double Bm1 = 0.0L;
   long double B0 = 1.0L;
   long double a = expl(x);
   long double b = -x + 1.0L;
   long double Ap1 = b * A0 + a * Am1;
   long double Bp1 = b * B0 + a * Bm1;
   int j = 1;

   a = 1.0L;
   while ( fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1) ) {
      if ( fabsl(Bp1) > 1.0L) {
         Am1 = A0 / Bp1;
         A0 = Ap1 / Bp1;
         Bm1 = B0 / Bp1;
         B0 = 1.0L;
      } else {
         Am1 = A0;
         A0 = Ap1;
         Bm1 = B0;
         B0 = Bp1;
      }
      a = -j * j;
      b += 2.0L;
      Ap1 = b * A0 + a * Am1;
      Bp1 = b * B0 + a * Bm1;
      j += 1;
   }
   return (-Ap1 / Bp1);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Power_Series_Ei( long double x )                        //
//                                                                            //
//  Description:                                                              //
//     For -5 < x < 6.8, the power series representation for                  //
//     (Ei(x) - gamma - ln|x|)/exp(x) is used, where gamma is Euler's gamma   //
//     constant.                                                              //
//     Note that for x = 0.0, Ei is -inf.  In which case -DBL_MAX is          //
//     returned.                                                              //
//                                                                            //
//     The power series expansion of (Ei(x) - gamma - ln|x|) / exp(x) is      //
//        - Sum(1 + 1/2 + ... + 1/j) (-x)^j / j!, where the Sum extends       //
//        from j = 1 to inf.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Power_Series_Ei( long double x )
{ 
   long double xn = -x;
   long double Sn = -x;
   long double Sm1 = 0.0L;
   long double hsum = 1.0L;
   long double g = 0.5772156649015328606065121L;
   long double y = 1.0L;
   long double factorial = 1.0L;
  
   if ( x == 0.0L ) return (long double) -DBL_MAX;
 
   while ( fabsl(Sn - Sm1) > epsilon * fabsl(Sm1) ) {
      Sm1 = Sn;
      y += 1.0L;
      xn *= (-x);
      factorial *= y;
      hsum += (1.0 / y);
      Sn += hsum * xn / factorial;
   }
   return (g + logl(fabsl(x)) - expl(x) * Sn);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Argument_Addition_Series_Ei(long double x)              //
//                                                                            //
//  Description:                                                              //
//     For 6.8 < x < 50.0, the argument addition series is used to calculate  //
//     Ei.                                                                    //
//                                                                            //
//     The argument addition series for Ei(x) is:                             //
//     Ei(x+dx) = Ei(x) + exp(x) Sum j! [exp(j) expj(-dx) - 1] / x^(j+1),     //
//     where the Sum extends from j = 0 to inf, |x| > |dx| and expj(y) is     //
//     the exponential polynomial expj(y) = Sum y^k / k!, the Sum extending   //
//     from k = 0 to k = j.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////
static long double Argument_Addition_Series_Ei(long double x)
{
   static long double ei[] = {
      1.915047433355013959531e2L,  4.403798995348382689974e2L,
      1.037878290717089587658e3L,  2.492228976241877759138e3L,
      6.071406374098611507965e3L,  1.495953266639752885229e4L,
      3.719768849068903560439e4L,  9.319251363396537129882e4L,
      2.349558524907683035782e5L,  5.955609986708370018502e5L,
      1.516637894042516884433e6L,  3.877904330597443502996e6L,
      9.950907251046844760026e6L,  2.561565266405658882048e7L,
      6.612718635548492136250e7L,  1.711446713003636684975e8L,
      4.439663698302712208698e8L,  1.154115391849182948287e9L,
      3.005950906525548689841e9L,  7.842940991898186370453e9L,
      2.049649711988081236484e10L, 5.364511859231469415605e10L,
      1.405991957584069047340e11L, 3.689732094072741970640e11L,
      9.694555759683939661662e11L, 2.550043566357786926147e12L,
      6.714640184076497558707e12L, 1.769803724411626854310e13L,
      4.669055014466159544500e13L, 1.232852079912097685431e14L,
      3.257988998672263996790e14L, 8.616388199965786544948e14L,
      2.280446200301902595341e15L, 6.039718263611241578359e15L,
      1.600664914324504111070e16L, 4.244796092136850759368e16L,
      1.126348290166966760275e17L, 2.990444718632336675058e17L,
      7.943916035704453771510e17L, 2.111342388647824195000e18L,
      5.614329680810343111535e18L, 1.493630213112993142255e19L,
      3.975442747903744836007e19L, 1.058563689713169096306e20L
   };
   int  k = (int) (x + 0.5);
   int  j = 0;
   long double xx = (long double) k;
   long double dx = x - xx;
   long double xxj = xx;
   long double edx = expl(dx);
   long double Sm = 1.0L;
   long double Sn = (edx - 1.0L) / xxj;
   long double term = DBL_MAX;
   long double factorial = 1.0L;
   long double dxj = 1.0L;

   while (fabsl(term) > epsilon * fabsl(Sn) ) {
      j++;
      factorial *= (long double) j;
      xxj *= xx;
      dxj *= (-dx);
      Sm += (dxj / factorial);
      term = ( factorial * (edx * Sm - 1.0L) ) / xxj;
      Sn += term;
   }
   
   return ei[k-7] + Sn * expl(xx); 
}

fp_t interpol_2d(fp_t ** f, fp_t * x, fp_t * y, int Nx, int Ny, fp_t x_to_interpol, fp_t y_to_interpol){

  // Now this one is computed using interpol_1d, and we do not care a lot what is interpol 1_d like.

  // Let us first find where is the location of the point in x:
  int i;
  for (i = 0; i<=Nx-2; ++i)
    if (x_to_interpol < x[i+1])
      break;

  // now we got i

  fp_t y_interpolated[3];
  fp_t x_restricted[3];

  if (i > 0 && i < Nx-1) {
    y_interpolated[0] = interpol_1d(f[i-1], y, Ny, y_to_interpol);
    y_interpolated[1] = interpol_1d(f[i], y, Ny, y_to_interpol);
    y_interpolated[2] = interpol_1d(f[i+1], y, Ny, y_to_interpol);
    x_restricted[0] = x[i-1]; x_restricted[1] = x[i]; x_restricted[2] = x[i+2]; 
  }
  else if (i == 0){
    y_interpolated[0] = interpol_1d(f[0], y, Ny, y_to_interpol);
    y_interpolated[1] = interpol_1d(f[1], y, Ny, y_to_interpol);
    y_interpolated[2] = interpol_1d(f[2], y, Ny, y_to_interpol);
    x_restricted[0] = x[0]; x_restricted[1] = x[1]; x_restricted[2] = x[2];
  }
  else  {
    y_interpolated[0] = interpol_1d(f[i-2], y, Ny, y_to_interpol);
    y_interpolated[1] = interpol_1d(f[i-1], y, Ny, y_to_interpol);
    y_interpolated[2] = interpol_1d(f[i], y, Ny, y_to_interpol);
    x_restricted[0] = x[i-2]; x_restricted[1] = x[i-1]; x_restricted[2] = x[i];
  }

  return interpol_1d(y_interpolated, x_restricted, 3, x_to_interpol);
}

fp_t interpol_2d(fp_t * f, fp_t * x, fp_t * y, int Nx, int Ny, fp_t x_to_interpol, fp_t y_to_interpol){

  fp_t ** fp;
  fp = ft2dim(0, Nx-1, 0, Ny-1);
  for (int i = 0; i<Nx; ++i)
    for (int j = 0; j<Ny; ++j)
      fp[i][j] = f[i * Ny + j];
  fp_t result = interpol_2d(fp, x, y, Nx, Ny, x_to_interpol, y_to_interpol);
  del_ft2dim(fp, 0, Nx-1, 0, Ny-1);
  return result;
}



fp_t interpol_1d(fp_t * f, fp_t * x, int N, fp_t x_to_interpol){

  // Function for 1D Bezier style interpolation
  // Assumes that x is monotonously increasing

  if (x_to_interpol < x[N-1] && x_to_interpol > x[0]){

    // Here is where the quadratic interpolation happens:
    int i;
    for (i = 0; i<=N-2; ++i)
      if (x_to_interpol < x[i+1])
        break;

    // Now we are interpolating between i and i+1
    // Quadratic interpolation, with derivatives from Fritsch and Butland
    // Now it is question what derivatives do we have:
    fp_t h1, h2;
    fp_t d1, d2;
    fp_t f_der_1, f_der_2;
    f_der_1 = f_der_2 = 0;
    fp_t C0, C1;
    C0 = C1 = 0.0;

    // Now we need to compute the derivatives in i and i+1, the question is can we do it. 
    // In i we can do it if it is greater than zero (it is, by default smaller then N-1)
    if (i > 0){ // Fritsch & Butland (1984)
      h1 = x[i] - x[i-1];
      h2 = x[i+1] - x[i];
      d1 = (f[i] - f[i-1])/h1;
      d2 = (f[i+1] - f[i])/h2;
      fp_t alpha = 1.0/3.0 * (1.0 + h2/(h1+h2));
      f_der_1 = (d1 * d2 > 0) ? d1*d2/(alpha*d2+(1.0-alpha)*d1) : 0.0;
    }
    else // linear
      f_der_1 = (f[1] - f[0]) / (x[1] - x[0]);
    // And one can compute C0
    C0 = f[i] + h2 * 0.5 * f_der_1;
    // If i+1 we can do it if i+2 exists, that is if i+1 < N-1 (this is more readable than i < N-2)
    if (i+1 < N-1){
      h1 = x[i+1] - x[i];
      h2 = x[i+2] - x[i+1];
      d1 = (f[i+1] - f[i])/h1;
      d2 = (f[i+2] - f[i+1])/h2;
      fp_t alpha = 1.0/3.0 * (1.0 + h2/(h1+h2));
      f_der_2 = (d1 * d2 > 0) ? d1*d2/(alpha*d2+(1.0-alpha)*d1) : 0.0;
    }
    else 
      f_der_2 = (f[N-1] - f[N-2]) / (x[N-1] - x[N-2]);
    
    C1 = f[i+1] - h1 * 0.5 * f_der_2;

    fp_t C = (C0 + C1) * 0.5;
    
    fp_t u = (x_to_interpol - x[i]) / (x[i+1] - x[i]);
    return (1.0-u) * (1.0-u) * f[i] + u * u * f[i+1] + u * 2.0 * (1.0 - u) * C;
  }
  else if (x_to_interpol >= x[N-1]) // If it is higher then highhest
    return f[N-1];
  else
    return f[0]; // Else it is lower then the lowest, return it
}

fp_t interpol_1d_linear(fp_t * f, fp_t * x, int N, fp_t x_to_interpol){

  if (x_to_interpol < x[0])
    return (f[1]-f[0]) / (x[1]-x[0]) * (x_to_interpol-x[0]);
  else if (x_to_interpol > x[N-1])
    return (f[N-1] - f[N-2]) / (x[N-1] - x[N-2]) * (x_to_interpol - x[N-1]);
  else {

    int k;
    for (k=0; k<N-1;k++)
      if (x_to_interpol < x[k+1]) break;

    return (f[k+1] - f[k]) / (x[k+1] - x[k]) * (x_to_interpol - x[k]);

  }  
}

// A function which in itself solves a determined (hopefully determined) linear system by means of svd. 
// We are doing this to see if we can get well behaved convergence. 


fp_t Planck_f(fp_t lambda, fp_t T){
  return 2.0 * 6.62606957 * 8.9875517E-7 / lambda / lambda / lambda / lambda / lambda / (exp(1.4387769 / lambda / T) - 1.0);
}

fp_t Planck_f_derivative(fp_t lambda, fp_t T){
  fp_t exponent = (1.4387769/lambda/T);
  fp_t temp = (exponent > 5.0) ? 1.0 : exp(exponent)/(exp(exponent) - 1.0);
  return Planck_f(lambda, T) * temp * 1.4387769/lambda/T/T;
}

fp_t compute_E(fp_t y0, fp_t y1, fp_t y2, fp_t h1, fp_t h2){
  
  // First compute local derivative. 

  fp_t d1 = (y1-y0) / h1;
  fp_t d2 = (y2-y1) / h2;

  fp_t derivative;

  if (d1 * d2 > 0 ){
    fp_t alpha = 0.3333333 * (1.0 + h2 / (h1+h2));
    derivative = d1 * d2 / (alpha * d2 + (1.0-alpha) * d1);
  }
  else
    derivative = 0.0;

  return y1 + h2 * 0.3333333 * derivative;

}

fp_t compute_F(fp_t y0, fp_t y1, fp_t y2, fp_t h1, fp_t h2){

  fp_t derivative = (y2 - y1) / h2;
  return y2 - h2 * 0.3333333 * derivative;
}


fp_t compute_w0(fp_t tb){

  fp_t etb = exp(-tb);

  if (tb < 1E-3) // Expansion
    return (tb * (0.5 - tb * (1.0/3.0  - tb/8.0)));
  else
    return (1.0 / tb - etb * (1.0 + 1.0 / tb));
}

fp_t compute_w1(fp_t tb){

  fp_t etb = exp(-tb);

  if (tb < 1E-3) // Expansion
    return (0.5 * tb * (1.0 - tb/3.0 * (1.0 - tb/4.0)));
  else 
    return (1.0 - 1.0 / tb + etb / tb);
    
}

int atmospheric_interpolation(fp_t * node_tau, fp_t * node_value, int N_nodes, fp_t * tau_grid, fp_t * quantity, int from, int to){

  // This is the function which performs a specific kind of interpolation for atmospheric quantities. We will use third order, strictly monotoneous interpolation. 
  // Best choice for this are probably Bezier splines as in De La Cruz Rodriguez & Piskunov (2013)

  // First step is to check amount of nodes, if we have more then one, there will be some work:
  if (N_nodes > 1){

    // First thing is to compute the derivatives everywhere: 
    fp_t * derivatives = new fp_t [N_nodes] - 1;

    derivatives[1] = 0.0;
    derivatives[N_nodes] = 0.0; // Derivative zero on the edges

    // In between:
    for (int i=2; i<=N_nodes-1;++i){
      fp_t d1,d2,h1,h2;
      h1 = node_tau[i] - node_tau[i-1];
      h2 = node_tau[i+1] - node_tau[i];
      d1 = (node_value[i] - node_value[i-1]) / h1;
      d2 = (node_value[i+1] - node_value[i]) / h2;
      derivatives[i] = 0.0;
      if (d1 * d2 > 0){
        fp_t alpha = 0.3333333 * (1.0 + h2 / (h1+h2));
        derivatives[i] = d1 * d2 / (alpha * d2 + (1.0-alpha) * d1);
      }
    }

    // And then use the derivative to interpolate:
    for (int i=from;i<=to;++i){

      // If we are out of bounds, just set to const:
      if (tau_grid[i] < node_tau[1])
        quantity[i] = node_value[1];

      else if (tau_grid[i] > node_tau[N_nodes])
        quantity[i] = node_value[N_nodes];

      else{

        // You first need to find where you are:
        int ii = 1;
        for (ii=1;ii<N_nodes;++ii){
          if (tau_grid[i] <= node_tau[ii+1])
            break;
          
        }
    
        fp_t E = node_value[ii] + (node_tau[ii+1] - node_tau[ii]) / 3.0 * derivatives[ii];
        fp_t F = node_value[ii+1] - (node_tau[ii+1] - node_tau[ii]) / 3.0 * derivatives[ii+1];
        fp_t u = (tau_grid[i] -  node_tau[ii]) / (node_tau[ii+1] - node_tau[ii]);
        quantity[i] = (1.0-u)*(1.0-u)*(1.0-u) * node_value[ii] + u*u*u*node_value[ii+1] + 3.0*u*(1.0-u)*(1.0-u)*E + 3.0*u*u*(1.0-u)*F; 
      }
    }
    delete[](derivatives+1);
    return 0;
  }
  else if (N_nodes == 1){ // Everything is constant
    for (int i=from;i<=to;++i)
      quantity[i] = node_value[1];
    return 0;
  }
  else if (N_nodes == 0)
    return 0;

  printf("Error in atmospheric interpolation : you are not supposed to pass negative number of nodes.\n");
  return 1;
}

fp_t vactoair(fp_t lambda_vac){

  fp_t s = 1E4/(lambda_vac*1E8);
  fp_t n =  1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s*s) + 0.0001599740894897 / (38.92568793293 - s*s);
  return lambda_vac / n;
}


fp_t airtovac(fp_t lambda_air){

  fp_t s = 1E4/(lambda_air*1E8);
  fp_t n = 1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s*s) + 0.0001599740894897 / (38.92568793293 - s*s);
  return lambda_air * n;
}

fp_t * airtovac(fp_t * lambda_air, int N){
  fp_t * lambda_out = new fp_t [N];
  for (int l=0;l<N;++l)
    lambda_out[l] = airtovac(lambda_air[l]);
  delete[]lambda_air;
  return lambda_out;
}