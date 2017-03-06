#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include "cmdln.h"
#include "types.h"

char *reducepath(char *path)
{
  char *p=strcpy(new char [strlen(path)+1],path);
  while(char *q=strstr(p,"/../")){
    int i=0;
    while(q[--i]!='/');
    q[++i]=0;
    q+=4;
    char *t=new char [strlen(p)+strlen(q)+1];
    sprintf(t,"%s%s",p,q);
    delete[] p;
    p=t;
  }
  return p;
}

int vrf_field(char *val,char type)
{
  switch(type){
    case('i'):{                     // int
      int a;
      return sscanf(val,"%d",&a);
    }
    case('d'):{                     // double
      double a;
      return sscanf(val,"%lf",&a);
    }
    case('s'): return (val)?strlen(val):1;  // string
  }
  return 0;
}

char **parse(char *line)
{
  int len=strlen(line);
  if(line[len-1]=='\n') line[--len]=0;
  for(int k=0,mode=1;k<=len;++k){
    if(line[k]=='"') mode=(!mode);
    if((line[k]==' ')&&(mode)) line[k]=0;
  }
  int m=0,*idx=new int [1],no=0;
  do{
    while((!line[m])&&(m<len)) ++m;
    if(m<len){
      int *t=new int [no+1];
      if(no) memcpy(t,idx,no*sizeof(int));
      delete[] idx;
      idx=t;
      idx[++no-1]=m;
      m+=strlen(line+m);
    }
  }while(++m<len);
//
  char esc[256];
  for(int k=0;k<=255;++k) esc[k]=k;
  esc['r']='\r';
  esc['t']='\t';
  esc['n']='\n';
  esc['f']='\f';
  esc['"']='"';
//
  char **tmp=new char* [no+1];
  for(int k=0;k<=no-1;++k){
    tmp[k]=strcpy(new char [strlen(line+idx[k])+1],line+idx[k]);
  // handle escaped characters
    char *q=tmp[k];
    while(q=strchr(q,'\\')){
      memmove(q,q+1,strlen(q));
      q[0]=esc[q[0]];
    }
  }
  delete[] idx;
  tmp[no]=0;
  return tmp;
}

char **split_fields(char *line,char sep)
{
  int len=strlen(line);
//
  char *p=line;
  while(p=strchr(p,sep)) (++p)[-1]=0;
//
  char **tmp=new char* [1];
  int no=0;
  p=line;
  while((p-line)<len){
    int sl=strlen(p)+1;
    if(sl>1){
      char **t=new char* [no+2];
      memcpy(t,tmp,no*sizeof(char*));
      delete[] tmp;
      tmp=t;
      tmp[no]=new char [sl];
      memcpy(tmp[no],p,sl);
      p+=sl;
      ++no;
    }
  }
  tmp[no]=0;
  return tmp;
}

void arg_extract(char *field,char *&def,char ***&map,char &type,char &opt)
{
  int k,ok=0;
  opt=0;
  for(k=0;!ok;++k)
    switch(field[k]){
      case('o'):{
        opt=1;
        break;
      }
      case('i'):
      case('d'):
      case('s'):{
        type=field[k];
	ok=1;
        break;
      }
      default:{
        exit(fprintf(stderr,"Unknown type identifier \'%c\'\n",field[k]));
        break;
      }
    }
  if(field[k]=='['){ // string -> numeric value mappings
    if(char *p=strchr(field+k,']')){
      int n=(p-field-k)-1,l;
      char *q=(char*)memcpy(new char [n+1],field+k+1,n);
      (p=q)[n]=0;
      for(l=0;p=strchr(p,',');++p) ++l;
// allocate matrix of string pointers
      map=new char** [l+2];                    // (0,l+1)
      map[0]=new char* [2*(l+1)];              // (0,l,0,2)
      for(int m=1;m<=l;++m) map[m]=map[m-1]+2; // set remaining pointers
      map[l+1]=0;                              // eot marker
//
      p=q;
      for(int m=0;m<=l;++m){
        char *s=strchr(p,'=')+1;
	s[-1]=0;
        char *t=strchr(s,',');
        (t=(t)?t+1:s+strlen(s)+1)[-1]=0;
        map[m][0]=strcpy(new char [strlen(p)+1],p);
        map[m][1]=strcpy(new char [strlen(s)+1],s);
	p=t;
      }
//
      delete[] q;
      k+=n+2;
    }else exit(fprintf(stderr,"No closing brackets found for mappings in \"%s\"\n",field));
  }else{
    map=new char** [1];
    map[0]=0;
  }
  if((field[k]=='"')&&(field[strlen(field+1)]=='"')){
    field[strlen(field+1)]=0;
    ++k;
  }
  def=strcpy(new char [strlen(field+k)+1],field+k);
  if(!strcmp(def,"<none>")){ // special keyword: no default value
    delete[] def;
    def=0;
  }
}

const char *types(char type)
{
  switch(type){
    case('i'): return "int";
    case('d'): return "double";
    case('s'): return "string";
  }
  return "unknown";
}

char **del_element(int n,char **argv)
{
  int argc;
  for(argc=0;argv[argc];++argc);
  char **tmpargs=new char* [argc];
  memcpy(tmpargs,argv,n*sizeof(char*));
  memcpy(tmpargs+n,argv+n+1,(argc-n-1)*sizeof(char*));
  tmpargs[argc-1]=0;
  delete[] argv;
  return tmpargs;
}

int get_opt(char **&argv,optel *opt)
{
  int n=0;
  while(argv[n]){
    int ok=0;
    for(int i=0;opt[i].opt_id;++i)
      for(int j=0;opt[i].opt_id[j];++j){
        if(!strcmp(opt[i].opt_id[j],argv[n])){                                  // got it!
          argv=del_element(n,argv);                                             // remove option
          if(opt[i].opt_opt[0]>=0){                                             // at least one argument
            for(int k=0;opt[i].opt_opt[k]>=0;++k){                              // assign arguments (if any)
              if(argv[n]){
	        if(opt[i].opt_def[k]) delete[] opt[i].opt_def[k];                                       // delete old value
                opt[i].opt_def[k]=strcpy(new char [strlen(argv[n])+1],argv[n]); // copy new value
                for(int m=0;opt[i].opt_map[k][m];++m)                           // for all mappings (if any)
                  if(!strcmp(opt[i].opt_map[k][m][0],opt[i].opt_def[k])){       // mapping found!
                    delete[] opt[i].opt_def[k];
                    opt[i].opt_def[k]=strcpy(new char [strlen(opt[i].opt_map[k][m][1])+1],opt[i].opt_map[k][m][1]);
                  }
                if(!vrf_field(opt[i].opt_def[k],opt[i].opt_type[k])){
                  const char *th[4]={"st","nd","rd","th"};
                  exit(fprintf(stderr,"Failed to convert %d%s command line argument \"%s\" for option \"%s\" to type \"%s\"\n%s\n",k+1,th[(k>3)?3:k],opt[i].opt_def[k],opt[i].opt_id[j],types(opt[i].opt_type[k]),opt[i].opt_hlp));
                }
                argv=del_element(n,argv);                                       // remove argument
              }else{
	        if(opt[i].opt_opt[k]==0) exit(fprintf(stderr,"Option \"%s\" requires more arguments\n",opt[i].opt_id[1]));
                opt[i].opt_def[k]=strcpy(new char [strlen("<none>")+1],"<none>");
              }
            }
          }else ++(opt[i].opt_type[0]);                                         // switch: increment
          j=0;
	  ok=1;                                                            // argv[n] found
          while(opt[i+1].opt_id) ++i;                                           // end loop over i
          while(opt[i].opt_id[j+1]) ++j;                                        // end loop over j
        }
      }
    if(!ok) ++n;
  }
  return n;
}

void map_line(char *s,char **mappings[])
{
  for(int k=0;mappings[k];++k)
    for(int l=1;mappings[k][l];++l)
      if(!strncmp(s,mappings[k][l],strlen(mappings[k][l]))){
        char *p=s+strlen(mappings[k][l]);   // end of keyword
        while(p[0]==' ') ++p;
        if(p[0]=='='){                      // got it: map this keyword to it's equivalent command line option
          char *q=s+strlen(mappings[k][0]); // end of command line option
          memmove(q,p,strlen(p)+1);
          memcpy(s,mappings[k][0],strlen(mappings[k][0]));
          q[0]=' ';
          return;
        }
      }
}

cmdln::cmdln(int &argc,char **argv_in,const char *lines_in[])
{
  cmdln(argc,argv_in,lines_in,0);
}

cmdln::cmdln(int &argc,char **argv_in,const char *lines_in[],const char *mappings_in[])
{
  char **lines=new char* [1];
  for(int nn=0;lines_in[nn];++nn){
    char **tmp=new char* [nn+2];
    memcpy(tmp,lines,nn*sizeof(char*));
    delete[] lines;
    lines=tmp;
    lines[nn]=strcpy(new char [strlen(lines_in[nn])+1],lines_in[nn]);
    lines[nn+1]=0;
  }
  char ***mappings=0;
  if(mappings_in){
    mappings=new char** [1];
    mappings[0]=0;
    for(int nn=0;mappings_in[nn];++nn){
      char ***tmp=new char** [nn+2];
      memcpy(tmp,mappings,nn*sizeof(char**));
      delete[] mappings;
      mappings=tmp;
      char *p;
      mappings[nn]=parse(p=strcpy(new char [strlen(mappings_in[nn])+1],mappings_in[nn]));
      delete[] p;
      mappings[nn+1]=0;
    }
  }
  opt=new optel [1];                                                  // start array
  for(int nn=0;lines[nn];++nn){
    optel *tmp=new optel [nn+2];                                      // extend array
    memcpy(tmp,opt,nn*sizeof(optel));
    delete[] opt;
    opt=tmp;
    opt[nn+1].opt_id=0;                                               // eot marker
//
// initialise new element
//
    char **fields=parse(lines[nn]);                                   // separate fields
    opt[nn].opt_id=split_fields(fields[0],'|');                       // separate aliases
    int k;
    for(k=0;fields[k];++k);                                           // find end of array
    if(fields[(k-=2)+1][0]=='h'){
      int len=strlen(fields[k+1]+1);
      memcpy(opt[nn].opt_hlp=new char [len+1],fields[k+1]+1,len+1);   // help field always last
      if((opt[nn].opt_hlp[0]=='"')&&(opt[nn].opt_hlp[len-1]=='"')) ((char*)memmove(opt[nn].opt_hlp,opt[nn].opt_hlp+1,len-1))[len-2]=0;
    }
    opt[nn].opt_def=new char* [k];
    opt[nn].opt_map=new char*** [k];
    opt[nn].opt_type=new char [k+1];
    opt[nn].opt_opt=new char [k+1];
    if(k>0){
      for(int l=0;l<=k-1;++l){
        arg_extract(fields[l+1],opt[nn].opt_def[l],opt[nn].opt_map[l],opt[nn].opt_type[l],opt[nn].opt_opt[l]);
        if(opt[nn].opt_def[l]){
          for(int m=0;opt[nn].opt_map[l][m];++m)                      // for all mappings (if any)
            if(!strcmp(opt[nn].opt_map[l][m][0],opt[nn].opt_def[l])){ // mapping found!
              delete[] opt[nn].opt_def[l];
              opt[nn].opt_def[l]=strcpy(new char [strlen(opt[nn].opt_map[l][m][1])+1],opt[nn].opt_map[l][m][1]);
            }
          if(!vrf_field(opt[nn].opt_def[l],opt[nn].opt_type[l])){
            const char *th[4]={"st","nd","rd","th"};
            exit(fprintf(stderr,"Error: Failed to convert default value for %d%s argument \"%s\" for option \"%s\" to type \"%s\"\n",l,th[(l>4)?3:l-1],opt[nn].opt_def[l],opt[nn].opt_id[0],types(opt[nn].opt_type[l])));
          }
        }
      }
    }
    opt[nn].opt_type[k]=0;                                      // flag: default status = clear
    opt[nn].opt_opt[k]=-1;
    for(int l=0;fields[l];++l) delete[] fields[l];
    delete[] fields;
//    for(int q=0;opt[nn].opt_id[q];++q) fprintf(stderr,"%s\n",opt[nn].opt_id[q]);
//    for(int q=0;opt[nn].opt_opt[q]>0;++q) fprintf(stderr,"\"%d\" \n",(int)(opt[nn].opt_type[q]));
  }
//
// store full path to the current process
//
  char *cwd=getcwd(NULL,10000);
  cwd[strlen(cwd)+1]=0;
  cwd[strlen(cwd)]='/';
  char *u=new char [strlen(cwd)+strlen(argv_in[0])+1];
  sprintf(u,"%s%s",cwd,argv_in[0]);
  free(cwd);
  path=reducepath(u);
  delete[] u;
//
// make a copy of the command line
//
  char **argv=new char* [argc+1];
  memcpy(argv,argv_in,argc*sizeof(char*));
  argv[argc]=0;
//
// find -cfg option in command line
//
/*
  optel *topt=new optel [2];
  topt[1].opt_id=0;
  for(int i=0;opt[i].opt_id;++i)
    for(int j=0;opt[i].opt_id[j];++j)
      if(!strcmp(opt[i].opt_id[j],"-cfg"))
        topt[0]=opt[i];
  int cfg=(argc-get_opt(argv,topt));                                  // extract option from the command line
  topt[0]=topt[1];                                                    // avoid deleting dynamically allocated elements of the original array
  delete[] topt;
//
// process the config file (if it exists)
//
  char *filename;
  if(get("-cfg",filename)){                                           // should be true if option in lines
    struct stat buf;
    if(stat(filename,&buf)){                                          // file does not exist
      if(cfg) fprintf(stderr,"Warning: config file \"%s\" not found\n",filename);
    }else{                                                            // file exists
      if(!cfg) fprintf(stderr,"Reading config file \"%s\"\n",filename);
      FILE *f=fopen(filename,"r");
      char s[10000];
      fgets(s,10000,f);
      int ln=1;
      do{
        if((s[0]!='#')&&(s[0]!=';')){
     	  map_line(s,mappings);
          char **fields=parse(s);				      // separate options
     	  get_opt(fields,opt);  				      // extract options
     	  for(int k=0;fields[k];++k) fprintf(stderr,"Warning: Ignoring unrecognised option \"%s\" in file \"%s\" line %d\n",fields[k],filename,ln);
     	  for(int k=0;fields[k];++k) delete fields[k];  	      // clean up
     	  delete[] fields;
        }
     	fgets(s,10000,f); // next line
     	++ln;
      }while(!feof(f));
      fclose(f);
    }
  }
*/
//
// process the rest of the command line
//
  if(mappings){
    for(int nn=0;mappings[nn];++nn){
      for(int mm=0;mappings[nn][mm];++mm) delete[] mappings[nn][mm];
      delete[] mappings[nn];
    }
    delete[] mappings;
  }
//
  get_opt(argv,opt);                                                  // extract left over options from the command line
  for(argc=0;argv[argc];++argc);                                      // update argc
  memcpy(argv_in,argv,argc*sizeof(char*));
  delete[] argv;
  for(int nn=0;lines[nn];++nn) delete[] lines[nn];
  delete[] lines;
}

cmdln::~cmdln(void)
{
  if(path) delete[] path;
  for(int i=0;opt[i].opt_id;++i){
    for(int j=0;opt[i].opt_id[j];++j) delete[] opt[i].opt_id[j];      // all id's
    delete[] opt[i].opt_id;
    for(int j=0;opt[i].opt_opt[j]>=0;++j){                            // all arguments
      if(opt[i].opt_def[j]) delete[] opt[i].opt_def[j];
      if(opt[i].opt_map[j]){
        for(int m=0;opt[i].opt_map[j][m];++m){                          // all mappings for each argument
          delete[] opt[i].opt_map[j][m][0];
          delete[] opt[i].opt_map[j][m][1];
        }
        delete[] opt[i].opt_map[j][0]; // data
        delete[] opt[i].opt_map[j]; // pointers
      }
      
    }
    delete[] opt[i].opt_def;
    delete[] opt[i].opt_map;
    delete[] opt[i].opt_hlp;
    delete[] opt[i].opt_opt;
    delete[] opt[i].opt_type;
  }
  delete[] opt;
}

char *cmdln::help(const char *opt_in)
{
  for(int i=0;opt[i].opt_id;++i)
    for(int j=0;opt[i].opt_id[j];++j)
      if(!strcmp(opt[i].opt_id[j],opt_in)) return opt[i].opt_hlp;
  exit(fprintf(stderr,"No help available for unrecognised option \"%s\"\n",opt_in));
  return NULL;
}

int cmdln::get(const char *opt_in,int &vo)
{
  for(int i=0;opt[i].opt_id;++i)
    for(int j=0;opt[i].opt_id[j];++j)
      if(!strcmp(opt[i].opt_id[j],opt_in)){
        if(opt[i].opt_opt[0]<0) return (vo=(int)(opt[i].opt_type[0])); // special treatment for flags
        if(opt[i].opt_type[1]!=0) exit(fprintf(stderr,"Option \"%s\" takes more than one argument\n",opt_in));
        return sscanf(opt[i].opt_def[0],"%d",&vo);
      }
  return 0;
}

int cmdln::get(const char *opt_in,fp_t &vo)
{
  for(int i=0;opt[i].opt_id;++i)
    for(int j=0;opt[i].opt_id[j];++j)
      if(!strcmp(opt[i].opt_id[j],opt_in)){
        if(opt[i].opt_type[1]!=0) exit(fprintf(stderr,"Option \"%s\" takes more than one argument\n",opt_in));
        const char *fmt=(sizeof(fp_t)>4)?"%lf":"%f";
        return (sizeof(fp_t)>4)?sscanf(opt[i].opt_def[0],fmt,&vo):sscanf(opt[i].opt_def[0],fmt,&vo);
      }
  return 0;
}

int cmdln::get(const char *opt_in,char *&vo)
{
  for(int i=0;opt[i].opt_id;++i)
    for(int j=0;opt[i].opt_id[j];++j)
      if(!strcmp(opt[i].opt_id[j],opt_in))
        return ((vo=opt[i].opt_def[0])!=0);                           // no default value: option not given
  return 0;
}

int cmdln::get(const char *opt_in)
{
  for(int i=0;opt[i].opt_id;++i)
    for(int j=0;opt[i].opt_id[j];++j)
      if(!strcmp(opt[i].opt_id[j],opt_in)) return 1;
  return 0;
}

