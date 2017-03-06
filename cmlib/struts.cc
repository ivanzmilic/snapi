#include <stdio.h>
#include <string.h>

#include "types.h"
#include "uts.h"

#include "struts.h"

char **append(char **a,const char *b)
{
  int la=0;
  if(a) while(a[la]) ++la;
//
  char **tmp=new char* [la+2];
//
  if(la){
    memcpy(tmp,a,la*sizeof(char*));
    delete[] a;
  }
  tmp[la]=strcpy(new char [strlen(b)+1],b);
  tmp[la+1]=0;
//
  return tmp;
}

char **append(char **a,char **b)
{
  int la=0,lb=0;
  if(a) while(a[la]) ++la;
  if(b) while(b[lb]) ++lb;
//
  char **tmp=new char* [la+lb+1];
//
  if(la) memcpy(tmp,a,la*sizeof(char*));
  if(lb) memcpy(tmp+la,b,lb*sizeof(char*));
  tmp[la+lb]=0;
//
  return tmp;
}

char **remove(char **a,const char *b)
{
  int l=0,i=-1;
  if(a){
    while(a[l]){
      if(!strcmp(a[l],b)) i=l;
      ++l;
    }
    if(i<0) return a;
    if(l){
      char **tmp=new char* [l];
      if(i) memcpy(tmp,a,i*sizeof(char*));
      memcpy(tmp+i,a+i+1,(l-i)*sizeof(char*)); // must copy the trailing null
      delete[] a;
      return tmp;
    }
  }
  return a;
}

char *dup_str(const char *s)
{
  return (s)?strcpy(new char [strlen(s)+1],s):0;
}

char **dup_str(char **s)
{
  if(s){
    int len=0;
    while(s[len]) ++len;
    char **rv=new char* [len+1];
    for(int i=0;i<len;++i) rv[i]=strcpy(new char [strlen(s[i])+1],s[i]);
    rv[len]=0;
    return rv;
  }
  return 0;
}

void del_str(char **s)
{
  if(s){
    for(int i=0;s[i];++i) delete[] s[i];
    delete[] s;
  }
}

int32_t cmp_str(const char *s1,const char *s2)
{
  if(s1&&s2) return strcmp(s1,s2);
  return (s1||s2);
}

int08_t str_tmpl_id(const char *id,const char *tmpl)
{
  int ch;
  if(sscanf(id,tmpl,&ch)==1){
    char s[strlen(id)+20];
    sprintf(s,tmpl,ch);
    if(!strcmp(s,id)) return (int08_t)ch;
  }
  return -1;
}

void swap(char *&a,char *&b)
{
  char *c=a;
  a=b;
  b=c;
}

void swap(char **&a,char **&b)
{
  char **c=a;
  a=b;
  b=c;
}

char **string2strings(const char *txt,io_class &io)
{
  int nl=0;
  const char *p=txt;
  char **data=0;
  while(const char *a=strchr(p,'"')){
    if(const char *b=strchr(a+1,'"')){
      int len=b-(a+1);
      char **tmp=new char* [nl+2];
      if(nl){
        memcpy(tmp,data,nl*sizeof(char*));
        delete[] data;
      }
      data=tmp;
      data[nl]=new char [len+1];
      memcpy(data[nl],a+1,len);
      data[nl][len]=0;
      data[++nl]=0;
      p=b+1;
    }else{
      io.msg(IOL_ERROR,"string2strings: parse error in string \"%s\"\n",nl);
      for(int l=0;data[l];++l) delete[] data[l];
      delete[] data;
      return 0;
    }
  }
  return data;
}

char *strings2string(char **txt,io_class &io)
{
  if(!txt) return 0;
  int sz=0;
  for(int i=0;txt[i];++i) sz+=strlen(txt[i])+3;
  char *p=new char [sz];
  sprintf(p,"\"%s\"",txt[0]);
  for(int i=1;txt[i];++i){
    strcat(p,",\"");
    strcat(p,txt[i]);
    strcat(p,"\"");
  }
  return p;
}

char *roman(int n)
{
  const char *singles[]={"","I","II","III","IV","V","VI","VII","VIII","IX"};
  const char *tens[]={"","X","XX","XXX","XL","L","LX","LXX","LXXX","XC"};
  const char *hundreds[]={"","C","CC","CCC","CD","D","DC","DCC","DCCC","CM"};
  const char *thousands[]={"","M","MM","MMM"};
  int th=(n%10000)/1000;
  int h=(n%1000)/100;
  int t=(n%100)/10;
  int s=n%10;
  char *rv=new char [strlen(thousands[th])+strlen(hundreds[h])+strlen(tens[t])+strlen(singles[s])+1];
  sprintf(rv,"%s%s%s%s",thousands[th],hundreds[h],tens[t],singles[s]);
  return rv;
}

int unroman(const char *s) // case insensitive
{
  int r=0,len=strlen(s);
  for(int i=0;i<len;++i){
    char chr=0xDF&(s[i]);
    char prv=0xDF&((i)?s[i-1]:0);
    char nxt=0xDF&(s[i+1]);
    switch(chr){
      case('I'):{
        switch(nxt){
          case('X'): r+=9; break;
          case('V'): r+=4; break;
          default: r+=1; break;
        }
        break;
      }
      case('V'):{
        if(prv!='I') r+=5;
        break;
      }
      case('X'):{
        if(prv!='I') // handled above
          switch(nxt){
            case('C'): r+=90; break;
            case('L'): r+=40; break;
            default: r+=10; break;
          }
        break;
      }
      case('L'):{
        if(prv!='X') r+=50;
        break;
      }
      case('C'):{
        if(prv!='X') // handled above
          switch(nxt){
            case('M'): r+=900; break;
            case('D'): r+=400; break;
            default: r+=100; break;
          }
        break;
      }
      case('D'):{
        if(prv!='C') r+=500;
        break;
      }
      case('M'):{
        if(prv!='C') r+=1000;
        break;
      }
      default: return -1; // error!
    } // switch
  } // loop
  return r;
}
