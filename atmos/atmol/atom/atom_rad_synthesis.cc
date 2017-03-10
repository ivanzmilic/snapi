#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>

#include "types.h"
#include "uts.h"
#include "io.h"
#include "pack.h"
#include "struts.h"
#include "mem.h"
#include "profiles.h"
#include "const.h"
#include "atomcfg.h"
#include "mathtools.h"

#include "H/hydrogen.h"
#include "He/helium.h"

#include "atomcfg.h"
#include "partf.h"  // partition functions
#include "bfcs.h"   // bound-free crossections
#include "colr.h"   // collisional rates
#include "atom.h"

// A separate file versions of functions which compute opacity, emissivity, and their perturbation variants. 
// We compute opacity and emissivity separately 
// The file contains following methods: 
// int atom::op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** B, fp_t theta,fp_t phi,
//   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector);

int atom::op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** Bmag, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  // Input quantities are usual ones, T, Ne, Vlos, Vt, B are atmospheric quantities
  // theta,phi,lambda and nlambda give us angles and wavelengths for which to compute opacity/emissivity
  // contributions from given atom are added to op_vector and em_vector.

  boundfree_op_em_vector(T,Ne,Vlos,theta,phi,lambda,nlambda,op_vector,em_vector);
  boundbound_op_em_vector(T,Ne,Vlos,Vt,Bmag,theta,phi,lambda,nlambda,op_vector,em_vector);

  return 0;
}

int atom::boundfree_op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  return 0;
}

int atom::boundbound_op_em_vector(fp_t*** T,fp_t*** Ne,fp_t*** Vlos,fp_t*** Vt, fp_t**** Bmag, fp_t theta,fp_t phi,
   fp_t* lambda,int nlambda,fp_t ****** op_vector, fp_t ***** em_vector){

  return 0;
}

