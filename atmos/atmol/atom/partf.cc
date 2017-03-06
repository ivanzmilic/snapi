#include <math.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "pack.h"
#include "uts.h"
#include "const.h"

#include "atomcfg.h"

#include "partf.h"

class pf *pf_new(int08_t z,struct pcfg *cfg)
{
  if(!cfg) return new pf(PF_TYPE_NONE);
  if(!strcmp(cfg->pftype,"TRAV")) return new trav(z,cfg->g0,cfg->pfc,cfg->npfc);
  if(!strcmp(cfg->pftype,"CONST")) return new cpf(cfg->value);
  return new pf(PF_TYPE_NONE);
}

class pf *pf_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  int08_t type;
  ::unpack(buf+offs,type);
  switch(type){
    case(PF_TYPE_TRAV): return new trav(buf,offs,do_swap,io);
    case(PF_TYPE_CONST): return new cpf(buf,offs,do_swap,io);
  }
  return new pf(buf,offs,do_swap,io);
}

//
// partition function base class
//

pf::pf(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

pf::pf(int08_t type_in)
{
  type=type_in;
}

//pf::pf(void)
//{
//  type=PF_TYPE_NONE;
//}

int32_t pf::size(io_class &io)
{
  return sizeof(int08_t);
}

int32_t pf::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::pack(buf,type);
  return offs;
}

int32_t pf::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=::unpack(buf,type);
  return offs;
}

//
// constant partition function
//

cpf::cpf(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):pf(buf,offs,do_swap,io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

cpf::cpf(fp_t val_in):val(val_in),pf(PF_TYPE_CONST)
{
}

int32_t cpf::size(io_class &io)
{
  int sz=pf::size(io);
  sz+=sizeof(fp_t);
  return sz;
}

int32_t cpf::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=pf::pack(buf,do_swap,io);
  offs+=::pack(buf+offs,val,do_swap);
  return offs;
}

int32_t cpf::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only unpack local stuff
  int32_t offs=::unpack(buf,val,do_swap);
  return offs;
}

//
// Traving et. al. style partition function
//

trav::trav(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io):pf(buf,offs,do_swap,io)
{
  offs+=unpack(buf+offs,do_swap,io);
}

trav::trav(int08_t z_in,fp_t g0_in,tpfcfg **cfg,int nc):pf(PF_TYPE_TRAV) // Traving et al (1965ZA.....61...92S)
{
  z=(fp_t)z_in+1.0;
  g0=g0_in;
  n=nc;
  g2=new fp_t [n];
  ee=new fp_t [n];
  ll=new fp_t [n];
  a=new fp_t* [n];
  g=new fp_t* [n];
  m=new int [n];
  asym=new int08_t [n];
  for(int i=0;i<n;++i){
    m[i]=cfg[i]->n;
    a[i]=new fp_t [m[i]];
    g[i]=new fp_t [m[i]];
    for(int mi=0;mi<m[i];++mi){
      a[i][mi]=cfg[i]->a[mi];
      g[i][mi]=cfg[i]->g[mi];
    }
    g2[i]=cfg[i]->g2;
    ee[i]=cfg[i]->ee;
    ll[i]=cfg[i]->ll;
    asym[i]=cfg[i]->asym;
  }
}

trav::~trav(void)
{
  for(int i=0;i<n;++i){
    delete[] a[i];
    delete[] g[i];
  }
  delete[] m;
  delete[] a;
  delete[] g;
  delete[] g2;
  delete[] ee;
  delete[] ll;
  delete[] asym;
}

//

int32_t trav::size(io_class &io)
{
  int sz=pf::size(io);
  sz+=(n+1)*sizeof(int); // n,*m
  sz+=n*sizeof(int08_t); // *asym
  sz+=2*sizeof(fp_t);    // z,g0
  sz+=3*n*sizeof(fp_t);  // *g2,*ee,*ll;
  for(int i=0;i<n;++i) sz+=2*m[i]*sizeof(fp_t); // **a,**g
  return sz;
}

int32_t trav::pack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
  int32_t offs=pf::pack(buf,do_swap,io);
// local stuff
  offs+=::pack(buf+offs,z,do_swap);
  offs+=::pack(buf+offs,g0,do_swap);
  offs+=::pack(buf+offs,n,do_swap);
  for(int i=0;i<n;++i){
    offs+=::pack(buf+offs,g2[i],do_swap);
    offs+=::pack(buf+offs,ee[i],do_swap);
    offs+=::pack(buf+offs,ll[i],do_swap);
    offs+=::pack(buf+offs,m[i],do_swap);
    offs+=::pack(buf+offs,g[i],0,m[i]-1,do_swap);
    offs+=::pack(buf+offs,a[i],0,m[i]-1,do_swap);
    offs+=::pack(buf+offs,asym[i]);
  }
  return offs;
}

int32_t trav::unpack(uint08_t *buf,uint08_t do_swap,io_class &io)
{
// only unpack local stuff
  int32_t offs=::unpack(buf,z,do_swap);
  offs+=::unpack(buf+offs,g0,do_swap);
  offs+=::unpack(buf+offs,n,do_swap);
// allocate arrays
  g2=new fp_t [n];
  ee=new fp_t [n];
  ll=new fp_t [n];
  g=new fp_t* [n];
  a=new fp_t* [n];
  m=new int [n];
  asym=new int08_t [n];
//
  for(int i=0;i<n;++i){
    offs+=::unpack(buf+offs,g2[i],do_swap);
    offs+=::unpack(buf+offs,ee[i],do_swap);
    offs+=::unpack(buf+offs,ll[i],do_swap);
    offs+=::unpack(buf+offs,m[i],do_swap);
    offs+=::unpack(buf+offs,g[i]=new fp_t [m[i]],0,m[i]-1,do_swap);
    offs+=::unpack(buf+offs,a[i]=new fp_t [m[i]],0,m[i]-1,do_swap);
    offs+=::unpack(buf+offs,asym[i]);
  }
  return offs;
}

//

fp_t trav::U(fp_t T,fp_t ne,io_class &io)
{
  fp_t x=5040.0/T;
// p is the highest existing level accoring to Fischel and Sparks (1971ApJ...164..359F)
// this should be looked at again at a later time...
//  fp_t f=4.2E3; // MULTI CGS
  fp_t f=1.405169521381833474E+05*((4.0*pi*e0)/sqr(e*Ryd*Ryd));
  fp_t p=z*f/pow(ne,1.0/6.0);   // MULTI CGS/SI invariant (hopefully)
  fp_t alpha=31.321*z*z*x;      // CGS/SI invariant
  fp_t sum=g0;
  for(int i=0;i<n;++i){
    for(int l=0;l<m[i];++l) sum+=a[i][l]*exp(-2.302585*g[i][l]*x);
    sum+=g2[i]*exp(-2.302585*ee[i]*x)*qas_f(ll[i],p,alpha);
  }
/*
  fp_t nhi=sqrt((13.595*z)/(4.98E-4*x*sqrt(10.0*ne*k*T))); // MULTI=Beckers in SI
  fp_t nhi=sqrt((13.595*z)/(4.98E-4*x*sqrt(ne*k*T)));    // MULTI=Beckers in CGS
  fp_t nhi=sqrt((4E8*zp*sqrt(T/ne)));                // Mihalas 2ed p112
    sum+=g2[i]*exp(-2.302585*ee[i]*x)*qas_b(ll[i],nhi,alpha);
*/
  return sum;
}

fp_t trav::dU(fp_t T,fp_t ne,io_class &io)
{
  fp_t x=5040.0/T;
  fp_t f=1.405169521381833474E+05*((4.0*pi*e0)/sqr(e*Ryd*Ryd));
  fp_t p=z*f/pow(ne,1.0/6.0);   // MULTI CGS/SI invariant (hopefully)
  fp_t alpha=31.321*z*z*x;      // CGS/SI invariant
  fp_t dsum=0.0;
  for(int i=0;i<n;++i)
    dsum=-g2[i]*exp(-2.302585*ee[i]*x)*dqas_f(ll[i],p,alpha)*z*f/(6.0*pow(ne,7.0/6.0));
  return dsum;
}

fp_t trav::qas_b(fp_t l,fp_t lh,fp_t alpha)
{
// Baschek et al. Abh. Hamb. VIII, 26 (1966)
  return (lh*(lh+1.0)*(lh+0.5)-l*(l+1.0)*(l+0.5))/3.0+alpha*(lh-l)+0.5*alpha*alpha*(lh-l)/(lh*l);
}

fp_t trav::qas_f(fp_t l,fp_t p,fp_t Dz)
{
// Fischel and Sparks (1971ApJ...164..359F.pdf)
// Simplified Eq. (26)
   if(p>l){
     fp_t llm=l-1.0;
     return ((1.3333333*p+0.5)*p+0.16666667+1.33333333*Dz)*p-0.4*Dz*Dz/p-
            ((0.33333333*llm-0.5)*llm-0.16666667-Dz)*llm+0.5*Dz*Dz/l;
   }
// Simplified Eq. (27)
  fp_t axl2=Dz/(l*l);
  return (1.0+axl2*(0.33333333+0.1*axl2))*p*p*p*p/l;
}

fp_t trav::dqas_f(fp_t l,fp_t p,fp_t Dz)
{
// Fischel and Sparks (1971ApJ...164..359F.pdf)
// Simplified Eq. (26)
   if(p>l){
     fp_t llm=l-1.0;
     return ((3.0*1.3333333*p+1.0)*p+0.16666667+1.33333333*Dz)+0.4*Dz*Dz/(p*p);
   }
// Simplified Eq. (27)
  fp_t axl2=Dz/(l*l);
  return 4.0*(1.0+axl2*(0.33333333+0.1*axl2))*p*p*p/l;
}

/* 
C
C***********************************************************************
C
      FUNCTION QTRAV(TETA,HP,J,JA)
C
C        HERE THE PARTITION FUNCTIONS ACCORDING TO TRAVING ET AL., ABH. HAMB.
C        STERNW. VIII, 1 (1966) ARE COMPUTED. THE SYMBOLS ARE GIVEN
C        IN THE COMMENTS AT THE BEGINNING OF SUBROUTINE INJON.
C        FUNCTION QAS IS CALLED.
C
C        DIMENSIONS NECESSARY
C        A(5),ASDE(KMAX),H(5),QPRIM(KMAX)
C        KMAX IS THE TOTAL NUMBER OF ELECTRON CONFIGURATIONS.
C        DIMENSIONS OF ARRAYS IN COMMON /CI3/ ARE COMMENTED ON IN SUBROUTINE
C        INJON.
C
C: QTRAV  92-11-27  MODIFICATIONS: (MATS CARLSSON)
C:        SAVE STATEMENT ADDED
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAMO'
      DIMENSION ASDE(MKMAX),QPRIM(MKMAX)
C
      COMMON/CI3/ALFA(MLMAX),GAMMA(MLMAX),G0(MJMAX),G2(MKMAX),
     * XION(MKMAX),XL(MKMAX),JBBEG(MJMAX),JCBEG(MJMAX),NK(MJMAX),
     * NL(MKMAX),IFISH
C
      COMMON/CI7/A(5),PFISH,ITP
C
      SAVE
C
C
C        STATEMENT FUNCTION FOR 10.**
      EXP10(X)=EXP(2.302585*X)
C
      FLJ=J
      JB=JBBEG(JA)
      JC1=JCBEG(JA)
      NKP=NK(JA)
      QSUM=0.
C
C        WE START THE LOOP OVER DIFFERENT ELECTRON CONFIGURATIONS ('THE K-LOOP'
      DO5 K=1,NKP
      JC2=NL(JB)+JC1-1
C
C        IS TETA=PRECEDING TETA
      IF(ITP.GT.0)GO TO 4
      PRA=XION(JB)*TETA
      IF(PRA.LT.12.)GO TO 1
      ASDE(JB)=0.
      GO TO 2
    1 ASDE(JB)=G2(JB)*EXP10(-PRA)
C
    2 QPRIM(JB)=0.
      IF(NL(JB).LE.0)GO TO 4
      DO3 L=JC1,JC2
      PRE=GAMMA(L)*TETA
      IF(PRE.GT.12.)GO TO 3
      QPRIM(JB)=QPRIM(JB)+ALFA(L)*EXP10(-PRE)
    3 CONTINUE
    4 JC1=JC2+1
      QSUM=QPRIM(JB)+ASDE(JB)*QAS(HP,XL(JB),A(J),FLJ,PFISH,IFISH)
     *             +QSUM
    5 JB=JB+1
C        END OF 'THE K-LOOP'
      QTRAV=G0(JA)+QSUM
C
      RETURN
      END
      FUNCTION QAS(H,XL,A,Z,PFISH,IFISH)
C
C        THIS ROUTINE COMPUTES THE ASYMPTOTIC PARTS OF THE PARTITION
C        FUNCTIONS FOLLOWING
C           BASCHEK ET AL., ABH. HAMB. VIII, 26 (1966) IF IFISH = 0
C           FISCHEL AND SPARKS, AP. J. 164, 359 (1971) IF IFISH = 1
C           (APPROXIMATING THE ZETA FUNCTIONS BY INTEGRALS).
C
C        XL=QUANTUM NUMBER FOR THE FIRST LEVEL OF THE ASYMPTOTIC PART
C        H=QUANTUM NUMBER OF THE CUT (FOR IFISH=0)
C        A=DZ(FISCHEL AND SPARKS)=ALFA(BASCHEK ET AL.)
C        PFISH=P(FISCHEL AND SPARKS), ONLY NECESSARY IF IFISH = 1
C
C
C
      INCLUDE 'PREC'
      COMMON/UTPUT/ IREAD, IWRIT
C
      QAS=0.0
      IF(IFISH.LT.0) RETURN
C
C        WHICH TYPE
      IF(IFISH.GT.0)GO TO 1
C
C        BASCHEK ET AL.
      QAS=0.333333*(H*(H+1.)*(H+0.5)-XL*(XL+1.)*(XL+0.5)) +
     *                       A*(H-XL)+0.5*A*A*(H-XL)/(H*XL)
      RETURN
C
C        FISCHEL AND SPARKS
    1 P=PFISH*Z
C
C        FISCHEL AND SPARKS, EQ. (26)
      P2=P*P
      P3=P2*P
      IF(P.LE.XL)GO TO 2
      XLM1=XL-1.
      R2=XLM1*XLM1
      R3=R2*XLM1
      QAS=1.3333333*P3+0.5*P2+0.16666667*P+1.33333333*A*P-0.4*A*A/P-
     *0.33333333*R3-0.5*R2-0.16666667*XLM1-A*XLM1+0.5*A*A/XL
      RETURN
C
C        FISCHEL AND SPARKS, EQ. (27)
    2 AXL2=A/(XL*XL)
      QAS=P*P*P*P*(1.+AXL2*(0.33333333+0.1*AXL2))/XL
      RETURN
      END
*/
