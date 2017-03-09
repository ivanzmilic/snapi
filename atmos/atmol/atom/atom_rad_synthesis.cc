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

// A separate file for "synthesise" versions of functions which compute opacity, emissivity, and their perturbation variants. 
// The file contains following methods: 
// fp_t **** atom:: opacity_vector_synth(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t*,int);
// fp_t ***  atom::emissivity_vector_synth(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t*,int);

fp_t **** atom:: opacity_vector_synth(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t*,int){
  

  return 0;
}

fp_t ***  atom::emissivity_vector_synth(fp_t***,fp_t***,fp_t***,fp_t***, fp_t****, fp_t,fp_t,fp_t*,int){
  return 0;
}