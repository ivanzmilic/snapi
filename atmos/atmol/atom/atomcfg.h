#ifndef  __ATOMCFG_H__ // __ATOMCFG_H__
#define  __ATOMCFG_H__

#include "types.h"
#include "io.h"

struct tpfcfg{  // Traving et. al. config info
  fp_t g2,ee,ll;
  fp_t *a,*g;
  int n;
  int08_t asym;
public:
  tpfcfg(char *tdata,io_class&);
  ~tpfcfg(void);
};

struct pcfg{    // partition config
  char *pftype; // type of PF to use
  fp_t value;   // for constant PF
  fp_t g0;      // for Traving et.al. PF
//
  tpfcfg **pfc;
  int npfc;
public:
  pcfg(char *pdata,io_class&);
  ~pcfg(void);
};

struct bfcfg{   // level bound-free crossection config
  char *bftype; // type of crossection to use
  fp_t *l,*v;
  int n;
public:
  bfcfg(char *bfdata,io_class&);
  ~bfcfg(void);
};

struct crcfg{   // transition collisional rates
  char *crtype; // type of rate to use
  fp_t *t,*v;
  int n;
public:
  crcfg(char*,io_class&);
  ~crcfg(void);
};

struct lcfg{    // level config
  fp_t *A;      // transitions
  int nt;
  fp_t ee;      // ecxitation energy
  fp_t g_lande; // lande factor
  uint08_t g;   // ecxitation weight <- statistical weights
  uint08_t j_qn; // l quantum no.
  uint08_t l_qn; // j quantum no.
  int n;
  bfcfg *bf;    // bound-free crossection
  crcfg **cr;   // collisional rates
  int ncr;
public:
  lcfg(char *ldata,io_class&);
  ~lcfg(void);
};

struct icfg{    // ion config
  int08_t charge;
  fp_t eion;
//
  pcfg *pf;
//
  lcfg **level;
  int nl;
public:
  icfg(int,char*,io_class&);
  ~icfg(void);
};

struct atmcfg{  // atom config
  char *id,*name;
  int08_t z;
  fp_t mass,abund;
  int nlte; // if is in nlte
//
  icfg **ion;
  int ni;
public:
  atmcfg(char *adata,io_class&);
  ~atmcfg(void);
};

 
#endif                 // __ATOMCFG_H__
