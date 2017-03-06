#ifndef __CMDLN_H__   // __CMDLN_H__
#define __CMDLN_H__

#include "types.h"

struct optel{ // option element
  char **opt_id,**opt_def,****opt_map,*opt_hlp;
  char *opt_opt,*opt_type; 
};

class cmdln{
  optel *opt;
  char *path;
  char **arg_v;
public:
  cmdln(int &argc,char **argv,const char *lines[]);
  cmdln(int &argc,char **argv,const char *lines[],const char *mappings[]);
  ~cmdln(void);
  char *fullpath(void){ return path; }
  char *help(const char*);	     // return help string
  int get(const char*,int&);       // int
  int get(const char*,fp_t&);    // double
  int get(const char*,char*&);     // string
  int get(const char*);	     // switch
};

#endif                 // __CMDLN_H__
