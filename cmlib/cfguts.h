#ifndef __CFGUTS_H__ // __CFGUTS_H__
#define __CFGUTS_H__ 

#include "types.h"

char *clean_path(char *wd);
char *file2buf(FILE *f);
char **scope_sep(char *fc,const char *scope_type,io_class &io);
char **str_sep(const char *str,char sep);
char *get_arg(char *args,const char *key,int isflag);
void append(int *&in,int nl,int &nh,int num);
char *arg_test(char *data);

#endif                 // __CFGUTS_H__ 
