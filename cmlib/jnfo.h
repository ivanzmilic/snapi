#ifndef __JNFO_H__ // __JNFO_H__
#define __JNFO_H__

#include <sys/types.h>
#include "types.h"
#include "io.h"

#include "atmos.h"
//#include "obs.h"
//#include "data.h"
//#include "inv.h"
// best make classes for each (e.g. atmos, obs, data, etc..)

struct jnfo{
  int na;
  class atmosphere **atmos;
  int nm;
  class model **models; // This is the model used to fit the data, or to synthesize from. Not necessary input. 
                       // if this is zero then we just synthesize from given atmosphere(s).

//
  int no;
  //class observable **obs;
//
  // These basically describe the observation so you do not have to specify observable per se. Probably MVN found it 
  // as well
  fp_t *az,*el,**lambda;
  int *nlambda;
  char **name;
  int *to_invert; ///whether to invert or not.
  int *xl,*xh,*yl,*yh;
  int *return_model; //whether to return model as the result of the inversion
  int *return_atmos; //whether to return full atmosphere as the result of the inversion
//
  uid_t uid;
  gid_t gid;
  char *uname;
//
  uint08_t cdcl,nmth,nsth,nsln; // chunk compression level, # master threads, # slave threads, # slaves


//
  jnfo(byte*,byte,io_class&);
  jnfo(void);
  ~jnfo(void);
  byte *compress(int&,int,byte,io_class&);
};

#endif             // __JNFO_H__
