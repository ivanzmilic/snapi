#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "types.h"
#include "io.h"
#include "mem.h"
#include "uts.h"
#include "mathtools.h"

#include "const.h"

#include "atmos.h"

uint08_t converged(fp_t eps,fp_t *N,fp_t *d,int natoms)
{
  for(int i=0;i<=natoms;++i) if(fabs(d[i]/N[i])>eps) return 0;
  return 1;
}

int08_t atomindex(class atmol *a,class atmol **b,int08_t n)
{
  if(!a) return 0; // electron
  for(int08_t i=1;i<=n;++i) if(a==b[i]) return i; // atom found
  return -1;       // not found!
}

uint08_t atmosphere::chemeq(class atmol **atml_in,int natm_in,fp_t T_in,fp_t Nt_in,fp_t &Ne_in,
                        int32_t x1i,int32_t x2i,int32_t x3i)
{
// ***************************************************************************************
// * chemical equilibruim: electrons are computed using ionization ratios and derivatives.
// * This cannot be correctly done in NLTE, so it must be converged iteratively
// * A more elegant solution is to use rate equations and solve the NLTE rate equations 
// * and chemical rate equations simultaneously, but this requires dissociation rates and 
// * we only have equilibrium constants
// ***************************************************************************************
  
  class atmol **atoms=0,**mols=0;
  int natoms=0,nmol=0,ncomp=0;
  for(int a=0,nc;a<natm_in;++a){
    if(struct ion *cmp=atml[a]->components(nc)){  // molecule
      ncomp+=nc;
      mols=append(atml[a],mols,nmol);     // atom
    }else atoms=append(atml[a],atoms,natoms);     // atom
  }
  if(nmol+natoms!=natm_in) io.msg(IOL_WARN,"atmosphere::chemeq: %d + %d != %d\n",nmol,natoms,natm_in);
  atoms-=1; // 1..natoms
//
// abundance number fractions (arbitrary scaling)
//
  fp_t *abund=new fp_t [natoms]-1;
  for(int08_t a=1;a<=natoms;++a) abund[a]=atoms[a]->abundance();   // return N_k/N_0
//
// molecular variables
//
  int *nc;
  int08_t **cmpidx,**cmpion,**molmul;
  fp_t *K,*q,*nm,**dm;
//
  if(nmol){  // we have some molecules
    mols-=1; // 1..nmol
// allocate 
    nc=new int [nmol]-1;
    cmpidx=new int08_t* [nmol]-1;
    cmpion=new int08_t* [nmol]-1;
    molmul=i08t2dim(1,nmol,0,natoms); // elemental multiplicity (needed for derivatives)
    K=new fp_t [nmol]-1;              // dissociation constant
    q=new fp_t [nmol]-1;              // molecular charge
    nm=new fp_t [nmol]-1;             // molecular number density
    dm=ft2dim(1,nmol,0,natoms);       // d nm[m]/d N[a]
// irregular matrices
    cmpidx[1]=new int08_t [ncomp];
    cmpion[1]=new int08_t [ncomp];
// initialize
    memset(cmpidx[1],0,ncomp*sizeof(int08_t));
    memset(cmpion[1],0,ncomp*sizeof(int08_t));
    memset(molmul[1],0,nmol*(natoms+1)*sizeof(int08_t));
    memset(q+1,0,nmol*sizeof(fp_t));
// formulate the equilibrium equations
    for(int08_t m=1;m<=nmol;++m)
      if(struct ion *cmp=mols[m]->components(nc[m])){ // this is a molecule
        if(m>1){
          cmpidx[m]=cmpidx[m-1]+nc[m-1];
          cmpion[m]=cmpion[m-1]+nc[m-1];
        }
        for(int i=0;i<nc[m];++i){
          int08_t idx=atomindex(cmp[i].atom,atoms,natoms);
          molmul[m][idx]+=1;                          // multiplicity
          cmpidx[m][i]=idx;                           // atomic index
          cmpion[m][i]=max(cmp[i].z,0);               // ionization index: limit to >=0
          q[m]+=(fp_t)(cmp[i].z);                     // total molecular charge
        }
        K[m]=mols[m]->dissoc(T_in,Ne_in);             // compute the equilibrium constant
      }
  }
// ***************************************
// * number densities starting values... *
// * This is important for stability of  *
// * the scheme. Can we stabilize this?  *
// ***************************************
  fp_t *N=new fp_t [natoms+1];     // Ne, N1..NN
  fp_t Np=Nt_in-(N[0]=Ne_in); // number of particles
//
  for(int08_t a=1;a<=natoms;++a) N[a]=Np*abund[a];
// variables for atoms and electrons
  fp_t *d=new fp_t [natoms+1],*v=new fp_t [natoms+1];
  fp_t **A=ft2dim(0,natoms,0,natoms),det;
  int *idx=new int [natoms+1];
  int *nk=new int [natoms]-1;
  fp_t **fk=new fp_t* [natoms+1];  // ionization fraction
  fp_t **df=new fp_t* [natoms+1];  // dfk/dne
  (fk[0]=new fp_t [1])[0]=1.0;     // all electrons in the "ground state"
  (df[0]=new fp_t [1])[0]=0.0;     // fraction is constant
  memset(fk+1,0,natoms*sizeof(fp_t*));
  memset(df+1,0,natoms*sizeof(fp_t*));
//
// compute suitable initialization values
// assuming only ionization and no molecules
//
//  Ne >= 1
  fp_t F1=0.0,S1=0.0;
  for(int a=1;a<=natoms;++a){
    atoms[a]->ionfrc(T_in,1.0,fk[a],df[a],nk[a]);   // f_k=N_k/N_1 and d f_k/d N[0]; nk=Z-1
    fp_t sum=0.0;
    for(int i=0;i<nk[a];++i) sum+=i*fk[a][i]; // count electrons
    F1+=abund[a]*sum;
    S1+=abund[a];
  }
  io.msg(IOL_DEB1,"atmosphere::chemeq: %E %E\n",1.0,F1/S1);
//  Ne <= Nt
  fp_t F3=0.0,S3=0.0;
  for(int a=1;a<=natoms;++a){
    atoms[a]->ionfrc(T_in,Nt_in,fk[a],df[a],nk[a]);   // f_k=N_k/N_1 and d f_k/d N[0]
    fp_t sum=0.0;
    for(int i=0;i<nk[a];++i) sum+=i*fk[a][i]; // count electrons
    F3+=abund[a]*sum;
    S3+=abund[a];
  }
  io.msg(IOL_DEB1,"atmosphere::chemeq: %E %E\n",Nt_in,F3/S3);
// solve for Ne assuming a linear interpolation for Ne
  fp_t aa=F3/(S3*Nt_in);
  fp_t bb=F1/S1-aa;
  fp_t cc=-bb*Nt_in;
  bb+=1.0-aa*Nt_in;
  N[0]=max((-bb+sqrt(bb*bb-4.0*aa*cc))/(2.0*aa),1.0);          // positive root is the first approximation for Ne
  io.msg(IOL_DEB1,"atmosphere::chemeq: Ne=%E (NR=%E)\n",N[0],(-bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa));
// "polish" the rough value
  fp_t esum;
  do{ // divide and conquer
    esum=0.0;
    for(int a=1;a<=natoms;++a){
      atoms[a]->ionfrc(T_in,N[0],fk[a],df[a],nk[a]);   // f_k=N_k/N_1 and d f_k/d N[0]
// why don't we solve this properly? (set up a linear system?)
      N[a]=(Nt_in-N[0])*abund[a]/S3;
      for(int i=0;i<nk[a];++i) esum+=fk[a][i]*N[a]*(fp_t)i;
//      io.msg(IOL_DEB1,"atmosphere::chemeq: %d %E\n",a,N[a]);
    }
    io.msg(IOL_DEB1,"atmosphere::chemeq: %d %E %E\n",0,N[0],esum);
    N[0]=0.5*(N[0]+esum); // the truth is somewhere in the middle
  }while((1E-2*(N[0]+esum))<fabs(N[0]-esum));
//
// we now have a fairly accurate ionization balance without molecules
//
  io.msg(IOL_XNFO,"atmosphere::chemeq: starting values: T=%E, Nt=%E:\n",T_in,Nt_in);
  for(int a=0;a<=natoms;++a) io.msg(IOL_XNFO,"atmosphere::chemeq: %d %E\n",a,N[a]);
// SVD variables
  fp_t **U=ft2dim(1,natoms+1,1,natoms+1);
  fp_t **V=ft2dim(1,natoms+1,1,natoms+1);
  fp_t *w=new fp_t [natoms+1]-1;
//
  int iter=0;
  do{ // the main Newton-Raphson loop
// ****************************************************
// * calculate ionization fractions and derivatives   *
// ****************************************************
    for(int a=1;a<=natoms;++a) atoms[a]->ionfrc(T_in,N[0],fk[a],df[a],nk[a]); // f_k=N_k/N_1 and d f_k/d N[0]
// ****************************************************
// * molecule number densities: this is a chemical    *
// * equilibrium which is a pressure equilibrium      *
// * (hence the factors kT)                           *
// ****************************************************
    for(int a=0;a<=natoms;++a) io.msg(IOL_XNFO,"atmosphere::chemeq: %d N[%d]=%E\n",iter,a,N[a]);
    for(int m=1;m<=nmol;++m){                  // molecule number densities
      nm[m]=1.0/(K[m]*k*T_in);                 // P=N*k*T was used
      memset(dm[m],0,(natoms+1)*sizeof(fp_t));
      for(int08_t i=0;i<nc[m];++i){            // loop over the components
        int08_t a=cmpidx[m][i];                // element index
        int08_t z=cmpion[m][i];                // charge
        nm[m]*=fk[a][z]*N[a]*k*T_in;           // P=N*k*T was used
        dm[m][a]+=1.0/N[a];                    // direct dependence on the N[a]
        dm[m][0]+=df[a][z]/fk[a][z];           // indirect dependence on the electron density [through the ionization balance]
      }
      for(int a=0;a<=natoms;++a) dm[m][a]*=nm[m];
      io.msg(IOL_DEB2,"atmosphere::chemeq: %d: moleq: nm[%d]=%E K=%E %E\n",iter,m,nm[m],K[m],1.0/(K[m]*k*T_in));
    }
// perhaps we should dampen the molecular populations to avoid destabilizing the N-R scheme?

// **************************
// * the defect(s)          *
// **************************
// charge conservation
    d[0]=-N[0];
    for(int a=1;a<=natoms;++a) for(int i=0;i<nk[a];++i) d[0]+=N[a]*fk[a][i]*(fp_t)i;    // electrons from ionization
    for(int m=1;m<=nmol;++m) d[0]+=q[m]*nm[m];                                          // charged molecules (HOWTO get q?)  
// specified total particle density
    d[1]=N[0]-Nt_in;
    for(int a=1;a<=natoms;++a) d[1]+=N[a];                                              // ions
    for(int m=1;m<=nmol;++m) d[1]+=nm[m];                                               // molecules  
// abundance ratios: reference element = 1
    fp_t Nr=N[1];                                                    // free reference atoms
    for(int m=1;m<=nmol;++m) Nr+=nm[m]*(fp_t)(molmul[m][1]);         // particles tied up in molecules
    for(int a=2;a<=natoms;++a){                                      // count all other elements
      fp_t Nc=N[a];                                                  // free atomic particles
      for(int m=1;m<=nmol;++m) Nc+=nm[m]*(fp_t)(molmul[m][a]);       // particles tied up in molecules
      d[a]=abund[a]*Nr-abund[1]*Nc;                                  // the defect
    }
    for(int08_t a=0;a<=natoms;++a) io.msg(IOL_DEB1,"atmosphere::chemeq: %d: defect: %d %E %E\n",iter,a,d[a],fabs(d[a]/N[a]));    
// **************************
// * calculate the Jacobian *
// **************************
// d [charge conservation]/d Ne
    A[0][0]=-1.0;
    for(int a=1;a<=natoms;++a) for(int i=1;i<nk[a];++i) A[0][0]+=N[a]*df[a][i]*(fp_t)i; // electrons from ionization
    for(int m=1;m<=nmol;++m) A[0][0]+=q[m]*dm[m][0];                                    // charged molecules (HOWTO get q?)
// d [total particle density]/d Ne
    A[1][0]=1.0;
//    for(int a=1;a<=natoms;++a) for(int i=0;i<nk[a];++i) A[1][0]+=N[a]*df[a][i]; // ions+electrons
    for(int m=1;m<=nmol;++m) A[1][0]+=dm[m][0];                                 // molecules  
// d [abundance ratios]/d Ne
    fp_t dNr=0.0;                                                    // count all reference atoms
    for(int m=1;m<=nmol;++m) dNr+=dm[m][0]*(fp_t)(molmul[m][1]);     // particles tied up in molecules
//
    for(int a=2;a<=natoms;++a){                                      // count all other atoms
      fp_t dNc=0.0;
      for(int m=1;m<=nmol;++m) dNc+=dm[m][0]*(fp_t)(molmul[m][a]);   // particles tied up in molecules
      A[a][0]=abund[a]*dNr-abund[1]*dNc;                             // the derivative
    }
    for(int b=1;b<=natoms;++b){                                      // derivative w.r.t. N[1..natoms]
// d [charge conservation]/d N[b]
      A[0][b]=0.0; // is it really 0?
      for(int i=0;i<nk[b];++i) A[0][b]+=fk[b][i]*(fp_t)i;            // electrons from ionization
      for(int m=1;m<=nmol;++m) A[0][b]+=q[m]*dm[m][b];               // electrons from charged molecules  
// d [total particle density]/d N[b]
      A[1][b]=1.0;
      for(int m=1;m<=nmol;++m) A[1][b]+=dm[m][b];                    // molecules  
// d [abundance ratios]/d N[b]
      dNr=0.0;                                                       // count all reference atoms
      for(int m=1;m<=nmol;++m) dNr+=dm[m][b]*(fp_t)(molmul[m][1]);   // particles tied up in molecules
      if(b==1) dNr+=1.0;                                             // atomic particles
      for(int a=2;a<=natoms;++a){                                    // count all other atoms
        fp_t dNc=0.0;
        for(int m=1;m<=nmol;++m) dNc+=dm[m][b]*(fp_t)(molmul[m][a]); // particles tied up in molecules
        if(b==a) dNc+=1.0;                                           // atomic particles
        A[a][b]=abund[a]*dNr-abund[1]*dNc;                           // the derivative
      }
    }
//
// SVD solution
//    
    memcpy(U[1]+1,A[0],(natoms+1)*(natoms+1)*sizeof(fp_t));
    if(svdcmp(U,natoms+1,natoms+1,w,V)<0) io.msg(IOL_WARN,"SVD: singular matrix?\n"); 
    for(int a=1;a<=natoms+1;++a) io.msg(IOL_DEB1,"SVD: w[%d]=%E\n",a,w[a]);
    fp_t wmax=w[1];
    for(int i=1;i<=natoms+1;++i) wmax=max(wmax,w[i]);
    for(int i=1;i<=natoms+1;++i){
      if((w[i])&&((wmax/w[i])<=1E30)){
        for(int j=1;j<=natoms+1;++j) V[j][i]/=w[i];
      }else{
        for(int j=1;j<=natoms+1;++j) V[j][i]=0.0;
        io.msg(IOL_WARN,"SVD: w[%d]: limit exceeded (%E<=%E*%E)!\n",i,wmax,w[i],1E30);
        exit(1);
      }
    }
    fp_t *t=new fp_t [natoms+1]-1;
    for(int i=1;i<=natoms+1;++i){
      t[i]=0.0;
      for(int j=1;j<=natoms+1;++j) t[i]+=U[j][i]*d[j-1];
    }
    for(int i=1;i<=natoms+1;++i){
      v[i-1]=0.0;
      for(int j=1;j<=natoms+1;++j) v[i-1]+=V[i][j]*t[j];
    }
    delete[] (t+1);      
//
// add the correction
//
    io.msg(IOL_DEB1,"atmosphere::chemeq: iteration %d:\n",++iter);
    fp_t fct=1.0;
    for(int a=0;a<=natoms;++a){ // using a "damping" constant significantly reduces the chance of overshoot
      io.msg(IOL_DEB1,"atmosphere::chemeq: %E f=%E v=%E f*v=%E %E %E\n",N[a],fct,v[a],fct*v[a],N[a]-fct*v[a],fabs(v[a]/(N[a]-fct*v[a])));
      N[a]-=fct*v[a];
      if(N[a]<1E-2) N[a]=1E-2; // avoid negative solutions: chose this factor wisely: not too low so that things can't recover
    }
// **********************************************************************************
// * convergence fraction: this gives problems at low temperatures where the electron
// * density is so low that it is of the order of the numerical accuracy of particle 
// * conservation perhaps it should be treated separately in that regime...?
// **********************************************************************************
    iter++;
    if (iter > 20) break;
  }while(!converged(1E-6,N,d,natoms));
//
  io.msg(IOL_XNFO,"atmosphere::chemeq: final values: T=%E, Nt=%E:\n",T_in,Nt_in);
  Ne_in=N[0];
  io.msg(IOL_XNFO,"atmosphere::chemeq: %d %E\n",0,N[0]);
  for(int a=1;a<=natoms;++a){
    atoms[a]->pupdate(N[a],fk[a],nk[a],x1i,x2i,x3i);
// pretty-print
    char *s=new char [nk[a]*13];
    s[0]=0;
    char tmp[20];
    for(int i=0;i<nk[a];++i){
      sprintf(tmp,"%s%10.8lf",(i)?",":"",fk[a][i]);
      strcat(s,tmp);
    }
    io.msg(IOL_XNFO,"atmosphere::chemeq: %d %E [%s]\n",a,N[a],s);
    delete[] s;
  }
  io.msg(IOL_XNFO,"atmosphere::chemeq: molecules:\n");
  for(int08_t m=1;m<=nmol;++m){ // this is a molecule
    nm[m]=1.0/(K[m]*k*T_in);
    for(int08_t i=0;i<nc[m];++i) nm[m]*=fk[cmpidx[m][i]][cmpion[m][i]]*N[cmpidx[m][i]]*k*T_in;
    mols[m]->pupdate(nm[m],fk[0],0,x1i,x2i,x3i);
// pretty-print
    io.msg(IOL_XNFO,"atmosphere::chemeq: %d %s: %E \n",m,mols[m]->get_frm(),nm[m]);
  }
  if(nmol){
    delete[] (K+1);
    delete[] (q+1);
    del_i08t2dim(molmul,1,nmol,0,natoms);
    del_i08t2dim(cmpidx,1,nmol,0,1);
    del_i08t2dim(cmpion,1,nmol,0,1);
    delete[] (nc+1);
    delete[] (nm+1);
    del_ft2dim(dm,1,nmol,0,natoms);
    delete[] (mols+1);
  }
//
  delete[] (w+1);
  del_ft2dim(U,1,natoms+1,1,natoms+1);
  del_ft2dim(V,1,natoms+1,1,natoms+1);
  delete[] v;
  delete[] d;
  del_ft2dim(A,0,natoms,0,natoms);
  delete[] idx;
  for(int i=0;i<=natoms;++i){ // change to atom internal arrays?
    delete[] fk[i];
    delete[] df[i];
  }
  delete[] (fk);
  delete[] (df);
  delete[] (nk+1);
  delete[] (atoms+1);
  delete[] (abund+1);
  delete[] N;
  return 0;     // return error code if necessary
}

int atmosphere::execute_chemeq_for_point(int x1i, int x2i, int x3i){

  chemeq(atml,natm,T[x1i][x2i][x3i],Nt[x1i][x2i][x3i],Ne[x1i][x2i][x3i],x1i,x2i,x3i);

  return 0;
}