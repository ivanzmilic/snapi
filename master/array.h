#ifndef __ARRAY_H__  // __ARRAY_H__
#define __ARRAY_H__

#include "job.h"
#include "slv.h"
#include "slvcfg.h"

void array_add(job_class*,job_class**&);
void array_del(job_class*,job_class**&);
void array_ins(job_class*,job_class**&,int);
void array_mv(job_class*,job_class**&,job_class**&);

void array_add(slv_class*,slv_class**&);
void array_del(slv_class*,slv_class**&);
void array_ins(slv_class*,slv_class**&,int);
void array_mv(slv_class*,slv_class**&,slv_class**&);

void array_add(struct chunk*,struct chunk**&);
void array_del(struct chunk*,struct chunk**&);
void array_ins(struct chunk*,struct chunk**&,int);
void array_mv(struct chunk*,struct chunk**&,struct chunk**&);

int append(int x,int *&p,int len);
int reduce(int *&p,int len,int pos);
int insert(int x,int *&p,int len,int pos);

int append(fp_t x,fp_t *&p,int len);
int reduce(fp_t *&p,int len,int pos);
int insert(fp_t x,fp_t *&p,int len,int pos);

void array_add(void*,void**&);
void array_del(void*,void**&);
void array_add(slvcfg*,slvcfg**&);
void array_del(slvcfg*,slvcfg**&);

#endif             // __ARRAY_H__
