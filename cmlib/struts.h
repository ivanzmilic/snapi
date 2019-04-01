#ifndef __STRUTS_H__  // __STRUTS_H__
#define __STRUTS_H__

#include "types.h"
#include "io.h"

char *add(char *,char*);
char **append(char**,const char*);
char **append(char**,char**);
char **remove(char**,const char*);
char *dup_str(const char*);
char **dup_str(char**);
void del_str(char**);
int32_t cmp_str(const char*,const char*);
int08_t str_tmpl_id(const char*,const char*);
void swap(char*&,char*&);
void swap(char**&,char**&);
char **string2strings(const char*,io_class&);
char *strings2string(char**,io_class&);
char *roman(int);
int unroman(const char*);

#endif                // __STRUTS_H__
