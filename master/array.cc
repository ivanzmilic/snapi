#include <string.h>
#include "array.h"
#include "job.h"
#include "uts.h"
#include "slv.h"

void array_add(job_class *p,job_class **&v)
{
  int len=0;
  while(v[len]) ++len;
  job_class **tv=new job_class* [len+2];
  memcpy(tv,v,len*sizeof(job_class*));
  delete[] v;
  v=tv;
  v[len]=p;
  v[len+1]=0;
}

void array_del(job_class *p,job_class **&v)
{
  int idx=-1,len;
  for(len=0;v[len];++len) if(v[len]==p) idx=len;
  if(idx<0) return;
  job_class **tv=new job_class* [len];
  if(idx) memcpy(tv,v,idx*sizeof(job_class*));
  if(len>idx) memcpy(tv+idx,v+idx+1,(len-idx-1)*sizeof(job_class*));
  delete[] v;
  v=tv;
  v[--len]=0;
}

void array_ins(job_class *p,job_class **&v,int pos)
{
  int len=0; 
  while(v[len]) ++len;
  if(pos>len) pos=len; // if pos>len append...
  job_class **tv=new job_class* [len+2];
  if(pos) memcpy(tv,v,pos*sizeof(job_class*));
  if(len>pos) memcpy(tv+pos+1,v+pos,(len-pos)*sizeof(job_class*));
  delete[] v;
  v=tv;
  v[pos]=p;
  v[len+1]=0;
}

void array_mv(job_class *p,job_class **&s,job_class**&d)
{
  array_add(p,d);
  array_del(p,s);
}

void array_add(slv_class *p,slv_class **&v)
{
  int len=0;
  while(v[len]) ++len;
  slv_class **tv=new slv_class* [len+2];
  memcpy(tv,v,len*sizeof(slv_class*));
  delete[] v;
  v=tv;
  v[len]=p;
  v[len+1]=0;
}

void array_del(slv_class *p,slv_class **&v)
{
  int idx=-1,len;
  for(len=0;v[len];++len) if(v[len]==p) idx=len;
  if(idx<0) return;
  slv_class **tv=new slv_class* [len];
  if(idx) memcpy(tv,v,idx*sizeof(slv_class*));
  if(len>idx) memcpy(tv+idx,v+idx+1,(len-idx-1)*sizeof(slv_class*));
  delete[] v;
  v=tv;
  v[--len]=0;
}

void array_ins(slv_class *p,slv_class **&v,int pos)
{
  int len=0; 
  while(v[len]) ++len;
  if(pos>len) pos=len; // if pos>len append...
  slv_class **tv=new slv_class* [len+2];
  if(pos) memcpy(tv,v,pos*sizeof(slv_class*));
  if(len>pos) memcpy(tv+pos+1,v+pos,(len-pos)*sizeof(slv_class*));
  delete[] v;
  v=tv;
  v[pos]=p;
  v[len+1]=0;
}

void array_mv(slv_class *p,slv_class **&s,slv_class**&d)
{
  array_add(p,d);
  array_del(p,s);
}

void array_add(struct chunk *p,struct chunk **&v)
{
  int len=0;
  while(v[len]) ++len;
  struct chunk **tv=new struct chunk* [len+2];
  if(len) memcpy(tv,v,len*sizeof(struct chunk*));
  delete[] v;
  v=tv;
  v[len]=p;
  v[len+1]=0;
}

void array_del(struct chunk *p,struct chunk **&v)
{
  int idx=-1,len;
  for(len=0;v[len];++len) if(v[len]==p) idx=len;
  if(idx<0) return;
  struct chunk **tv=new struct chunk* [len];
  if(idx) memcpy(tv,v,idx*sizeof(struct chunk*));
  if(len>idx) memcpy(tv+idx,v+idx+1,(len-idx-1)*sizeof(struct chunk*));
  delete[] v;
  v=tv;
  v[--len]=0;
}

void array_ins(struct chunk *p,struct chunk **&v,int pos)
{
  int len=0; 
  while(v[len]) ++len;
  if(pos>len) pos=len; // if pos>len append...
  struct chunk **tv=new struct chunk* [len+2];
  if(pos) memcpy(tv,v,pos*sizeof(struct chunk*));
  if(len>pos) memcpy(tv+pos+1,v+pos,(len-pos)*sizeof(struct chunk*));
  delete[] v;
  v=tv;
  v[pos]=p;
  v[len+1]=0;
}

void array_mv(struct chunk *p,struct chunk **&s,struct chunk**&d)
{
  array_add(p,d);
  array_del(p,s);
}

int append(int x,int *&p,int len)
{
  int *tmp=new int [len+1]-1;
  if(p){
    if(len) memcpy(tmp+1,p+1,len*sizeof(int));
    delete[] (p+1);
  }
  p=tmp;
  p[len+1]=x;
  return len+1;
}

int reduce(int *&p,int len,int pos)
{
  pos=min(pos,len);  // if pos > len take last element
  int *tmp=new int [len-1]-1;
  if(pos>1) memcpy(tmp+1,p+1,(pos-1)*sizeof(int));
  if(len<pos) memcpy(tmp+pos,p+pos+1,(len-pos)*sizeof(int));
  delete[] (p+1);
  p=tmp;
  return len-1;
}

int insert(int x,int *&p,int len,int pos)
{
  pos=min(pos,len); // if pos > len insert as last element
  int *tmp=new int [len+1]-1;
  if(p){
    if(pos>1) memcpy(tmp+1,p+1,(pos-1)*sizeof(int));
    if(len<pos) memcpy(tmp+pos+1,p+pos,(len-pos)*sizeof(int));
    delete[] (p+1);
  }
  p=tmp;
  p[pos]=x;
  return len+1;
}

int append(fp_t x,fp_t *&p,int len)
{
  fp_t *tmp=new fp_t [len+1]-1;
  if(p){
    if(len) memcpy(tmp+1,p+1,len*sizeof(fp_t));
    delete[] (p+1);
  }
  p=tmp;
  p[len+1]=x;
  return len+1;
}

int reduce(fp_t *&p,int len,int pos)
{
  pos=min(pos,len);  // if pos > len take last element
  fp_t *tmp=new fp_t [len-1]-1;
  if(pos>1) memcpy(tmp+1,p+1,(pos-1)*sizeof(fp_t));
  if(len<pos) memcpy(tmp+pos,p+pos+1,(len-pos)*sizeof(fp_t));
  delete[] (p+1);
  p=tmp;
  return len-1;
}

int insert(fp_t x,fp_t *&p,int len,int pos)
{
  pos=min(pos,len); // if pos > len insert as last element
  fp_t *tmp=new fp_t [len+1]-1;
  if(p){
    if(pos>1) memcpy(tmp+1,p+1,(pos-1)*sizeof(fp_t));
    if(len<pos) memcpy(tmp+pos+1,p+pos,(len-pos)*sizeof(fp_t));
    delete[] (p+1);
  }
  p=tmp;
  p[pos]=x;
  return len+1;
}

void array_add(slvcfg *p,slvcfg **&v)
{
  int len=0;
  while(v[len]) ++len;
  slvcfg **tv=new slvcfg* [len+2];
  if(len) memcpy(tv,v,len*sizeof(slvcfg*));
  delete[] v;
  v=tv;
  v[len]=p;
  v[len+1]=0;
//  for(int i=0;v[i];++i) v[i]->update(v);
}

void array_del(slvcfg *p,slvcfg **&v)
{
  int idx=-1,len;
  for(len=0;v[len];++len) if(v[len]==p) idx=len;
  if(idx<0) return;
  slvcfg **tv=new slvcfg* [len];
  if(idx) memcpy(tv,v,idx*sizeof(slvcfg*));
  if(len>idx) memcpy(tv+idx,v+idx+1,(len-idx-1)*sizeof(slvcfg*));
  delete[] v;
  v=tv;
  v[--len]=0;
//  for(int i=0;v[i];++i) v[i]->update(v);
}

void array_add(void *p,void **&v)
{
  int len=0;
  while(v[len]) ++len;
  void **tv=new void* [len+2];
  memcpy(tv,v,len*sizeof(void*));
  delete[] v;
  v=tv;
  v[len]=p;
  v[len+1]=0;
}

void array_del(void *p,void **&v)
{
  int idx=-1,len;
  for(len=0;v[len];++len) if(v[len]==p) idx=len;
  if(idx<0) return;
  void **tv=new void* [len];
  if(idx) memcpy(tv,v,idx*sizeof(void*));
  if(len>idx) memcpy(tv+idx,v+idx+1,(len-idx-1)*sizeof(void*));
  delete[] v;
  v=tv;
  v[--len]=0;
}

