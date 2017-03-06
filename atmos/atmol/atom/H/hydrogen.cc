#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "types.h"
#include "uts.h"
#include "const.h"

#include "atomcfg.h"

#include "hydrogen.h"

hydrogen::hydrogen(atmcfg *cfg,io_class &io_in):atom(cfg,io_in)
{
//  io_in.msg(IOL_INFO,"hydrogen::hydrogen: %s[%s]\n",cfg->name,cfg->id);
}

hydrogen::hydrogen(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in):atom(buf,offs,do_swap,io_in)
{
//  io_in.msg(IOL_INFO,"hydrogen::hydrogen: %lX %d\n",buf,offs);
}

hydrogen::~hydrogen(void)
{
}

// CHECK THIS! Had problem with lf and lt starting at 0, not at 1....
fp_t hydrogen::cbb(uint08_t zf,uint16_t lf,uint16_t lt,atmol**,uint16_t,fp_t Temp,fp_t ne,int32_t ix1,int32_t ix2,int32_t ix3)
{
// Collisional B-B rates from collisions with electrons
//
  if(lt<=lf) io.msg(IOL_ERROR|IOL_FATAL,"hydrogen::cbb: lowel level (%d) > upper level (%d)\n",lf,lt); // Upward collisions only
  if((lf==0)&&(lt==1)){       // Ly-alpha: Janev et al. (1987) pp. 18-21 and 257
    fp_t x=log(k*Temp/ev);    // polynomial variable x is in eV
// The 1S to 2P transitions
    fp_t C1[]={-2.814949375869E1,1.009828023274E1,-4.771961915818,1.467805963618,-2.979799374553E-1,
                3.861631407174E-2,-3.051685780771E-3,1.335472720988E-4,-2.476088392502E-6};
    fp_t RH1=C1[0]+x*(C1[1]+x*(C1[2]+x*(C1[3]+x*(C1[4]+x*(C1[5]+x*(C1[6]+x*(C1[7]+x*C1[8])))))));
// The 1S to 2S transitions
    fp_t C2[]={-2.833259375256E1,9.587356325603,-4.833579851041,1.415863373520,-2.537887918825E-1,
                2.800713977946E-2,-1.871408172571E-3,6.986668318407E-5,-1.123758504195E-6};
    fp_t RH2=C2[0]+x*(C2[1]+x*(C2[2]+x*(C2[3]+x*(C2[4]+x*(C2[5]+x*(C2[6]+x*(C2[7]+x*C2[8])))))));
    return ne*exp(RH1)+exp(RH2);
  }else{                      // all other transitions: Vriens and Smeets (1980), eqn. (17)
    fp_t EION=109678.764;
    fp_t ecm[]={0.000,82259.047,97492.279,102823.882,105291.646,106632.157,107440.441,107965.049,108324.718};
    fp_t ai=lf+1,aj=lt+1;
    fp_t ef=(lf<9)?ecm[lf]:EION*(1.0-1.0/(ai*ai));
    fp_t et=(lt<9)?ecm[lt]:EION*(1.0-1.0/(aj*aj));
// Johnson's (1972) formula for oscillator strength (see Janev et al. (1987), pp. 315-316)
    fp_t x=1.0-sqr(ai/aj),gmx;
    switch(lf){
      case(0): gmx=1.1330+(-0.4059+0.07014/x)/x; break;
      case(1): gmx=1.0785+(-0.2319+0.02947/x)/x; break;
      default:{ // >=2
        fp_t g0=   0.9935+(0.2328-0.1296/ai)/ai;
        fp_t g1= -(0.6282-(0.5598-0.5299/ai)/ai)/ai;
        fp_t g2=  (0.3887-(1.1810-1.4700/ai)/ai)/(ai*ai);
        gmx=g0+(g1+g2/x)/x;
      }
    }
    fp_t f=1.96028052*gmx*(ai/aj)/(x*x*x*ai*aj); // Oscillator strength - AGREES WITH VALUES IN JANEV ET AL. (1987), P. 315
// ef, et, EION and eng all in cm-1
    fp_t derg=et-ef;        // transisiton energy
    fp_t DION=EION-ef;      // ionization energy from the lower level=
// APN, BP and BPN of Vriens and Smeets (1980), eqns. (11), (13) and (12)
    fp_t amn=2.0*EION*f/derg;
    fp_t bm=(1.4*log(ai)-0.7-(0.51-(1.16-0.55/ai)/ai)/ai)/ai;
    fp_t bmn=4.0*sqr(EION/derg)/(aj*aj*aj)*(1.0+(4.0/3.0)*(DION/derg)+bm*sqr(DION/derg));
//
    fp_t eng=k*Temp/(h*c);  // ENERGY ASSOCIATED TO TEMPERATURE
    fp_t S=fabs(lf-lt);
// DELTA_PN and GAMMA_PN of Vriens and Smeets (1980), eqns. (18) and (19)
    fp_t DELMN=exp(-bmn/amn)+0.06*sqr(S/ai)/aj;
    fp_t GAMMN=EION*(3.0+11.0*sqr(S/ai))*log(1.0+ai*ai*ai*eng/EION)/(6.0+1.6*aj*S+0.3/(S*S)+0.8*aj*sqrt(aj/S)*fabs(S-0.6));
    return ne*1.6E-7*sqrt(ev/(h*c))*(amn*log(0.3*eng/EION+DELMN)+bmn)*sqrt(eng)*exp(-derg/eng)/(eng+GAMMN);
  }
  return 0.0;
}

#define H_Eion 109678.764       // Bashkin & Stoner (1975), P. 3

fp_t hydrogen::cbf(uint08_t zf,uint16_t lf,atmol**,uint16_t,fp_t Temp,fp_t ne,int32_t ix1,int32_t ix2,int32_t ix3)
{
// Collisional B-F rates from collisions with electrons
//
// Hydrogen mean energy levels taken from Bashkin and Stoner (1975), PP. 2-3
// E_ion=109678.764: Hydrogen ionization level taken from Bashkin and Stoner (1975), P. 3
// Collisional B-F ionization rates using polynomial fits by Vriens and Smeets (198), eqn. (8)
//
  fp_t y=(k*Temp)/ev;       // Thermal energy in ev
  switch(lf){
    case(0):{ // Ground level H1S TO H+
      fp_t x=log(y);  // log(Thermal energy in ev)
      fp_t C[]={-3.271396786375E1,1.353655609057E1,-5.739328757388,1.563154982022,-2.877056004391E-1,
                 3.482559773737E-2,-2.631976175590E-3,1.119543953861E-4,-2.039149852002E-6};
      return ne*exp(C[0]+x*(C[1]+x*(C[2]+x*(C[3]+x*(C[4]+x*(C[5]+x*(C[6]+x*(C[7]+x*C[8]))))))));
    }
    case(1):{ // First level H2S TO H+
      fp_t x=log(y);  // log(Thermal energy in ev)
      fp_t C[]={-1.973476726029E1,3.992702671457,-1.773436308973,5.331949621358E-1,-1.181042453190E-1,
                 1.763136575032E-2,-1.616005335321E-3,8.093908992682E-5,-1.686664454913E-6};
      return ne*exp(C[0]+x*(C[1]+x*(C[2]+x*(C[3]+x*(C[4]+x*(C[5]+x*(C[6]+x*(C[7]+x*C[8]))))))));
    }
    default:{ // all other cases
      fp_t ecm[]={0.000,82259.047,97492.279,102823.882,105291.646,106632.157,107440.441,107965.049,108324.718};
      fp_t eng=(lf<9)?ecm[lf]:H_Eion*(1.0-1.0/(fp_t)((lf+1)*(lf+1)));
      fp_t x=h*c*(H_Eion-eng)/(k*Temp);
      return ne*9.56E-6*exp(-x)/(pow(y,1.5)*(pow(x,2.33)+4.38*pow(x,1.72)+1.32*x)); // units?
    }
  }
  return 0.0;
}

/*
fp_t hydrogen::cfb(uint08_t zt,uint16_t lt,atmol**,uint16_t,fp_t Temp,fp_t ne,int32_t ix1,int32_t ix2,int32_t ix3)
{
  fp_t Cbf=0.0; // collisions with electrons
//
// NSTAR(IHI,K) / NSTAR(ILO,K) is the Boltzmann factor (gu/gl) exp(-dE/kT)
//
// pop[ix1][ix2][ix3].N[0];
// MULTI:
// CDN = NE(K) * CT * G(ILO) * SQRT(TEMP(K)) / G(IHI)
// CUP= CDN * NSTAR(IHI,K) / NSTAR(ILO,K)
// C(IHI,ILO,K) = CDN + C(IHI,ILO,K)
// C(ILO,IHI,K) = CUP + C(ILO,IHI,K)

// typically: electrons, then neutral hydrogen, then protons, HeI/HeII/HeIII if it's hot, then some other
// molecules (mostly charge transfer reactions) and that's all... 

// collisions with other species?
  return Cbf;
}*/

/*
C
C***********************************************************************
C
      FUNCTION H1BB (I,J,T)
C
C  DETERMINES H COLLISIONAL RATE COEFFICIENTS UP TO N=50
C
C CONVERSION FACTORS:
C
C      CM**-1 TO EV    HH * CC / EE
C      CM**-1 TO K     HH * CC / BK
C
C      K     TO HZ    BK / HH
C      CM**-1 TO HZ    CC
C
C CONSTANTS:
C
C      EE   ELECTRON CHARGE      1.602189E-12         
C      HH   PLANCK CONSTANT      (ERG S)
C      CC   VELOCITY OF LIGHT    (CM S**-1)
C      EM   ELECTRON REST MASS   (G)
C      UU   ATOMIC MASS UNIT     (G)
C      BK   BOLTZMANN CONSTANT   (ERG K**-1)
C      PI   3.14159265359
C
C
C INPUT:
C
C   I        LOWER LEVEL
C   J        UPPER LEVEL
C   T        TEMPERATURE (K)
C
C OUTPUT:
C
C   H1BB    UPWARDS COLLISIONAL RATE                  I TO J
C
C: H1BB   90-02-28  MODIFICATIONS: (MARTIN J. STIFT)
C:        REWRITTEN.  UP TO 50 BOUND LEVELS + ONE CONTINUUM LEVEL
C:
      INCLUDE 'PREC'
C
      INTEGER I,J,K
      DOUBLE PRECISION
     *        T,S,X,
     *        AI,AMN,AJ,BM,BMN,DELMN,DERG,ENG,
     *        DION,FIJ,GAMMN,GMX,G0,G1,G2,
     *        EION,EV(50),EE,HH,CC,EM,UU,BK,PI,TEV,
     *        B1S2P(9),B1S2S(9),RH1,RH2,H1BB
C
C POLYNOMIAL FITS FOR COLLISIONAL RATE COEFFICIENTS (LYMAN ALPHA)
C
C H(1S) TO H (2P)
C
      DATA B1S2P / -2.814949375869D+1,
     *              1.009828023274D+1,
     *             -4.771961915818D+0,
     *              1.467805963618D+0,
     *             -2.979799374553D-1,
     *              3.861631407174D-2,
     *             -3.051685780771D-3,
     *              1.335472720988D-4,
     *             -2.476088392502D-6/
C
C H(1S) TO H(2S)
C
      DATA B1S2S / -2.833259375256D+1,
     *              9.587356325603D+0,
     *             -4.833579851041D+0,
     *              1.415863373520D+0,
     *             -2.537887918825D-1,
     *              2.800713977946D-2,
     *             -1.871408172571D-3,
     *              6.986668318407D-5,
     *             -1.123758504195D-6/
C
C HYDROGEN MEAN ENERGY LEVELS TAKEN FROM BASHKIN & STONER (1975), PP. 2-3
C
      DATA EV /      0.000,
     *           82259.047,
     *           97492.279,
     *          102823.882,
     *          105291.646,
     *          106632.157,
     *          107440.441,
     *          107965.049,
     *          108324.718,
     *          41 * 0.000/
C
C HYDROGEN IONIZATION LEVEL TAKEN FROM BASHKIN & STONER (1975), P. 3
C
      DATA EION / 109678.764/
C
C MISCELLANOUS CONSTANTS
C
      DATA EE / 1.602189D-12/,
     *     HH / 6.626176D-27/,
     *     CC / 2.99792458D10/,
     *     EM / 9.109534D-28/,
     *     UU / 1.6605655D-24/,
     *     BK / 1.380662D-16/,
     *     PI / 3.14159265359/
C
C J MUST BE GREATER THAN I
C
      IF (J.LE.I) CALL STOP('H1BB: ERROR IN ORDER OF ENERGY LEVELS')
C
      AI = FLOAT(I)
      AJ = FLOAT(J)
C
C DETERMINE ENERGY LEVELS UP TO N=50
C
      DO 1 K = 10,50
         EV(K) = EION * (1. - 1. / FLOAT(K*K))
    1 CONTINUE
C
C EV, EION, ENG  IN   CM**-1
C
C ENG   ENERGY ASSOCIATED TO TEMPERATURE
C
      ENG = T * BK / HH / CC
C
C ENERGY DIFFERENCE AND IONISATION ENERGY FROM LOWER LEVEL
C
      DERG = EV(J) - EV(I)
      DION = EION  - EV(I)
C
C JOHNSON'S (1972) FORMULA FOR OSCILLATOR STRENGTH
C SEE JANEV ET AL. (1987), PP. 315-316
C
      X = 1. - (AI/AJ) * (AI/AJ)
C
      IF (I.EQ.1) THEN
         G0 =  1.1330
         G1 = -0.4059
         G2 =  0.07014
      ELSE IF (I.EQ.2) THEN
         G0 =  1.0785
         G1 = -0.2319
         G2 =  0.02947
      ELSE IF (I.GE.3) THEN
         G0 =   0.9935 + (0.2328 - 0.1296 / AI) / AI
         G1 = -(0.6282 - (0.5598 - 0.5299 / AI) / AI ) / AI
         G2 =  (0.3887 - (1.1810 - 1.4700 / AI) / AI ) / (AI * AI)
      ENDIF
C
      GMX = G0 + (G1 + G2 / X) / X
C
C OSCILLATOR STRENGTH - AGREES WITH VALUES IN JANEV ET AL. (1987), P. 315
C
      FIJ = 1.96028052 * GMX * (AI/AJ) / (AJ * AJ * X * X * X)
C
C APN, BP AND BPN OF VRIENS AND SMEETS (1980), EQNS. (11), (13) AND (12)
C
      AMN = 2. * EION * FIJ / DERG
C
      BM = ( 1.4 * LOG(AI) - 0.7 - (0.51 - (1.16 - 0.55 / AI) / 
     *     AI ) / AI ) / AI
C
      BMN = 4. * (EION/DERG) * (EION/DERG) / (AJ * AJ * AJ) * 
     *      (1. + (4./3.) * (DION/DERG) + 
     *       BM * (DION/DERG) * (DION/DERG) )
C
C DELTA_PN AND GAMMA_PN OF VRIENS AND SMEETS (1980), EQNS. (18) AND (19)
C
      S = ABS(AI - AJ)
      DELMN = EXP(-BMN/AMN) + 0.06 * (S/AI) * (S/AI) / AJ
      GAMMN = EION * (3. + 11. * (S/AI)*(S/AI)) *
     *        LOG(1. + AI * AI * AI * ENG / EION) /
     *        (6. + 1.6 * AJ * S + 0.3 / (S * S) + 
     *         0.8 * AJ * SQRT(AJ/S) * ABS(S-0.6))
C
C UPWARD COLLISIONAL RATE  (I TO J) : 
C VRIENS AND SMEETS (1980), EQN. (17), EXCEPT LYMAN ALPHA
C LYMAN ALPHA ACCORDING TO JANEV ET AL. (1987) PP. 18-21 AND 257
C
C TEV IS IN EV IN POLYNOMIAL EXPRESSIONS
C
      IF (I.EQ.1 .AND. J.EQ.2) THEN
C
         TEV = LOG(T * BK / EE)
C
         RH1 = B1S2P(1) + TEV * (B1S2P(2) + TEV * (B1S2P(3) +
     *                    TEV * (B1S2P(4) + TEV * (B1S2P(5) +
     *                    TEV * (B1S2P(6) + TEV * (B1S2P(7) +
     *                    TEV * (B1S2P(8) + TEV *  B1S2P(9) )))))))
         RH2 = B1S2S(1) + TEV * (B1S2S(2) + TEV * (B1S2S(3) +
     *                    TEV * (B1S2S(4) + TEV * (B1S2S(5) +
     *                    TEV * (B1S2S(6) + TEV * (B1S2S(7) +
     *                    TEV * (B1S2S(8) + TEV *  B1S2S(9) )))))))
C
         H1BB = EXP(RH1) + EXP(RH2)
C
      ELSE
C
         H1BB = 1.6D-7 * SQRT (EE / HH / CC) *
     *           (AMN * LOG(0.3D0 * ENG / EION + DELMN) + BMN) *
     *           SQRT(ENG) * EXP(-DERG/ENG) / (ENG + GAMMN)
C
      ENDIF
C
      RETURN
      END
C
C*************************************************************************
C
      FUNCTION H1BF(I,T)
C
C  DETERMINES H COLLISIONAL IONIZATION COEFFICIENTS UP TO N=50
C
C CONVERSION FACTORS:
C
C      CM**-1 TO EV    HH * CC / EE
C      CM**-1 TO K     HH * CC / BK
C
C      K     TO HZ    BK / HH
C      CM**-1 TO HZ    CC
C
C CONSTANTS:
C
C      EE   ELECTRON CHARGE      1.602189E-12         
C      HH   PLANCK CONSTANT      (ERG S)
C      CC   VELOCITY OF LIGHT    (CM S**-1)
C      EM   ELECTRON REST MASS   (G)
C      UU   ATOMIC MASS UNIT     (G)
C      BK   BOLTZMANN CONSTANT   (ERG K**-1)
C      PI   3.14159265359
C
C
C INPUT:
C
C   I        LOWER LEVEL
C   T        TEMPERATURE (K)
C
C OUTPUT:
C
C   H1BF      COLLISIONAL IONIZATION RATE               I TO K
C
C...................................................................
C
C: H1BF   90-02-28  MODIFICATIONS: (MARTIN J. STIFT)
C:        REWRITTEN.  UP TO 50 BOUND LEVELS + ONE CONTINUUM LEVEL
C:
      INCLUDE 'PREC'
C
      INTEGER I,K
      DOUBLE PRECISION
     *        T,RATIK,AI,ENG,DION,
     *        EION,EV(50),EE,HH,CC,EM,UU,BK,PI,TEV,
     *        B1SI(9),B2SI(9),H1BF
C
C H1S TO H+
C
      DATA B1SI /  -3.271396786375D+1,
     *              1.353655609057D+1,
     *             -5.739328757388D+0,
     *              1.563154982022D+0,
     *             -2.877056004391D-1,
     *              3.482559773737D-2,
     *             -2.631976175590D-3,
     *              1.119543953861D-4,
     *             -2.039149852002D-6/
C
C H2S TO H+
C
      DATA B2SI /  -1.973476726029D+1,
     *              3.992702671457D+0,
     *             -1.773436308973D+0,
     *              5.331949621358D-1,
     *             -1.181042453190D-1,
     *              1.763136575032D-2,
     *             -1.616005335321D-3,
     *              8.093908992682D-5,
     *             -1.686664454913D-6/
C
C HYDROGEN MEAN ENERGY LEVELS TAKEN FROM BASHKIN & STONER (1975), PP. 2-3
C
      DATA EV /      0.000,
     *           82259.047,
     *           97492.279,
     *          102823.882,
     *          105291.646,
     *          106632.157,
     *          107440.441,
     *          107965.049,
     *          108324.718,
     *          41 * 0.000/
C
C HYDROGEN IONIZATION LEVEL TAKEN FROM BASHKIN & STONER (1975), P. 3
C
      DATA EION / 109678.764/
C
C MISCELLANOUS CONSTANTS
C
      DATA EE / 1.602189E-12/,
     *     HH / 6.626176E-27/,
     *     CC / 2.99792458E10/,
     *     EM / 9.109534E-28/,
     *     UU / 1.6605655E-24/,
     *     BK / 1.380662E-16/,
     *     PI / 3.14159265359/
C
      AI = FLOAT(I)
C
C DETERMINE ENERGY LEVELS UP TO N=50
C
      DO 1 K = 10,50
         EV(K) = EION * (1. - 1. / FLOAT(K*K))
    1 CONTINUE
C
C DION AND ENG  IN   CM**-1
C
C ENG   ENERGY ASSOCIATED TO TEMPERATURE
C DION  IONISATION ENERGY FROM LOWER LEVEL
C
      ENG  = T * BK / HH / CC
      DION = EION  - EV(I)
C
C COLLISIONAL IONIZATION RATE  (I TO K) : VRIENS AND SMEETS (198), EQN. (8)
C
      IF (I.EQ.1) THEN
C
         TEV = LOG(T * BK / EE)
         RATIK = B1SI(1) + TEV * (B1SI(2) + TEV * (B1SI(3) +
     *                     TEV * (B1SI(4) + TEV * (B1SI(5) +
     *                     TEV * (B1SI(6) + TEV * (B1SI(7) +
     *                     TEV * (B1SI(8) + TEV *  B1SI(9) )))))))
         H1BF = EXP(RATIK)
C
      ELSE IF (I.EQ.2) THEN
C
         TEV = LOG(T * BK / EE)
         RATIK = B2SI(1) + TEV * (B2SI(2) + TEV * (B2SI(3) +
     *                     TEV * (B2SI(4) + TEV * (B2SI(5) +
     *                     TEV * (B2SI(6) + TEV * (B2SI(7) +
     *                     TEV * (B2SI(8) + TEV *  B2SI(9) )))))))
         H1BF = EXP(RATIK)
C
      ELSE
C
         H1BF = 9.56D-6 * EXP(-DION/ENG) / ENG / SQRT(ENG) *
     *     SQRT(EE / HH / CC) * (EE / HH / CC) / 
     *    ( (DION/ENG)**2.33 + 4.38*(DION/ENG)**1.72 + 1.32*(DION/ENG) )
C
      ENDIF
C
      RETURN
      END
C
C*************************************************************************
C
      SUBROUTINE HCOL
C
C  COLLISIONAL RATES FOR HYDROGEN.
C:
C: HCOL   90-02-28  MODIFCATIONS: (MARTIN J. STIFT)
C:        REWRITTEN.  UP TO 50 BOUND LEVELS + ONE CONTINUUM LEVEL
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CCONST'
C
      DOUBLE PRECISION T,H1BB,H1BF
C
C  CALCULATE C BOUND-BOUND AND BOUND-FREE
C  REFERENCE: VRIENS AND SMEETS, PHYS. REV. A, 22, 940 (1980)
C             JANEV ET AL. (1987) 
C                  ELEMENTARY PROCESSES IN HYDDROGEN-HELIUM PLASMAS
C
      DO 4 K = 1,NDEP
C
        T = TEMP(K)
C
        DO 2 J = 2,NK-1
          DO 1 I = 1,J-1
            C(I,J,K) = NE(K) * H1BB (I,J,T)
            C(J,I,K) = C(I,J,K) * NSTAR(I,K) / NSTAR(J,K)
    1     CONTINUE
    2   CONTINUE
C
        DO 3 I = 1,NK-1
          C(I,NK,K) = NE(K) * H1BF (I,T)
          C(NK,I,K) = C(I,NK,K) * NSTAR(I,K) / NSTAR(NK,K)
    3   CONTINUE
C
    4 CONTINUE
C
      RETURN
      END
*/
