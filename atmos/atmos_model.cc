#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "types.h"
#include "io.h"
#include "mem.h"
#include "pack.h"
#include "atmos.h"
#include "profiles.h"
#include "const.h"

#define T_MIN 3400.0
#define T_MAX 12000.
#define VT_MIN 1E2
#define VT_MAX 15E5
#define VS_MIN -15E5
#define VS_MAX 15E5
#define B_MIN 1.0
#define B_MAX 1E4

model::model(){

  N_nodes_temp = 0;
  N_nodes_vt = 0;
  N_nodes_vs = 0;
  N_nodes_B = 0;
  N_nodes_theta = 0;
  N_nodes_phi = 0;

  temp_nodes_tau = 0;
  temp_nodes_temp = 0;
  vt_nodes_tau = 0;
  vt_nodes_vt = 0;
  vs_nodes_tau = 0;
  vs_nodes_vs = 0;
  B_nodes_tau = 0;
  B_nodes_B = 0;
  theta_nodes_tau = 0;
  theta_nodes_theta = 0;
  phi_nodes_tau = 0;
  phi_nodes_phi = 0;

  response_to_parameters = 0;
  tau_min = -5.0;
  tau_max = 1.0;
  N_depths = 41;

}

model::model(int N_nodes_temp_in, int N_nodes_vt_in, int N_nodes_vs_in, int N_nodes_B_in){

  // If numbers of nodes is > 0 then set the number of nodes and allocate memory for nodes
  // If it is zero, then set pointer to zero and that means we are reading the values from the file
  // If it is smaller then zero, return error
  
  if (N_nodes_temp_in){
  	N_nodes_temp = N_nodes_temp_in;
  	temp_nodes_tau = new fp_t [N_nodes_temp]-1;
  	temp_nodes_temp = new fp_t [N_nodes_temp]-1;
  }
  else if (N_nodes_temp_in == 0){
  	N_nodes_temp = 0;
  	temp_nodes_tau = 0;
  	temp_nodes_temp = 0;
  }
  else{
  	printf("model::cannot create model with negative number of nodes in temperature \n"); // TO BE REPLACED WITH IO.
  	exit(1);
  }

  if (N_nodes_vt_in){
  	N_nodes_vt = N_nodes_vt_in;
  	vt_nodes_tau = new fp_t [N_nodes_vt]-1;
  	vt_nodes_vt = new fp_t [N_nodes_vt]-1;
  }
  else if (N_nodes_vt_in == 0){
  	N_nodes_vt = 0;
  	vt_nodes_tau = 0;
  	vt_nodes_vt = 0;
  }
  else{
  	printf("model::cannot create model with negative number of nodes in vt \n");
  	exit(1);
  }

  if (N_nodes_vs_in){
  	N_nodes_vs = N_nodes_vs_in;
  	vs_nodes_tau = new fp_t [N_nodes_vs]-1; // just forget vx and vy... Just forget them man.
  	vs_nodes_vs = new fp_t [N_nodes_vs]-1;
  }
  else if (N_nodes_vs_in == 0){
  	N_nodes_vs = 0;
  	vs_nodes_tau = 0;
  	vs_nodes_vs = 0;
  }
  else{
  	printf("model::cannot create model with negative number of nodes in vs \n");
  	exit(1);
  }

  if (N_nodes_B_in){
  	N_nodes_B = N_nodes_B_in;
  	B_nodes_tau = new fp_t [N_nodes_B]-1;
  	B_nodes_B = new fp_t [N_nodes_B]-1;
  }
  else if (N_nodes_B_in == 0){
  	N_nodes_B = 0;
  	B_nodes_tau = 0;
  	B_nodes_B = 0;
  }
  else{
  	printf("model::cannot create model with negative number of nodes in B \n");
  	exit(1);
  } 
  N_nodes_theta = 0;
  N_nodes_phi = 0;
  response_to_parameters = 0;
  N_depths = 41;
  tau_min = -5.0;
  tau_max = 1.0;
}

model::model(int N_nodes_temp_in, int N_nodes_vt_in, int N_nodes_vs_in, int N_nodes_B_in, int N_nodes_theta_in, int N_nodes_phi_in){

  // If numbers of nodes is > 0 then set the number of nodes and allocate memory for nodes
  // If it is zero, then set pointer to zero and that means we are reading the values from the file
  // If it is smaller then zero, return error
  
  // If numbers of nodes is > 0 then set the number of nodes and allocate memory for nodes
  // If it is zero, then set pointer to zero and that means we are reading the values from the file
  // If it is smaller then zero, return error
  
  if (N_nodes_temp_in){
    N_nodes_temp = N_nodes_temp_in;
    temp_nodes_tau = new fp_t [N_nodes_temp]-1;
    temp_nodes_temp = new fp_t [N_nodes_temp]-1;
  }
  else if (N_nodes_temp_in == 0){
    N_nodes_temp = 0;
    temp_nodes_tau = 0;
    temp_nodes_temp = 0;
  }
  else{
    printf("model::cannot create model with negative number of nodes in temperature \n"); // TO BE REPLACED WITH IO.
    exit(1);
  }

  if (N_nodes_vt_in){
    N_nodes_vt = N_nodes_vt_in;
    vt_nodes_tau = new fp_t [N_nodes_vt]-1;
    vt_nodes_vt = new fp_t [N_nodes_vt]-1;
  }
  else if (N_nodes_vt_in == 0){
    N_nodes_vt = 0;
    vt_nodes_tau = 0;
    vt_nodes_vt = 0;
  }
  else{
    printf("model::cannot create model with negative number of nodes in vt \n");
    exit(1);
  }

  if (N_nodes_vs_in){
    N_nodes_vs = N_nodes_vs_in;
    vs_nodes_tau = new fp_t [N_nodes_vs]-1; // vx, vy, vz // Latter two probably not used
    vs_nodes_vs = new fp_t [N_nodes_vs]-1;
  }
  else if (N_nodes_vs_in == 0){
    N_nodes_vs = 0;
    vs_nodes_tau = 0;
    vs_nodes_vs = 0;
  }
  else{
    printf("model::cannot create model with negative number of nodes in vs \n");
    exit(1);
  }

  if (N_nodes_B_in){
    N_nodes_B = N_nodes_B_in;
    B_nodes_tau = new fp_t [N_nodes_B]-1;
    B_nodes_B = new fp_t [N_nodes_B]-1;
  }
  else if (N_nodes_B_in == 0){
    N_nodes_B = 0;
    B_nodes_tau = 0;
    B_nodes_B = 0;
  }
  else{
    printf("model::cannot create model with negative number of nodes in B \n");
    exit(1);
  }

  if (N_nodes_theta_in){
    N_nodes_theta = N_nodes_theta_in;
    theta_nodes_tau = new fp_t [N_nodes_theta]-1;
    theta_nodes_theta = new fp_t [N_nodes_theta]-1;
  }
  else if (N_nodes_theta_in == 0){
    N_nodes_theta = 0;
    theta_nodes_tau = 0;
    theta_nodes_theta = 0;
  }
  else{
    printf("model::cannot create model with negative number of nodes in theta \n");
    exit(1);
  }

  if (N_nodes_phi_in){
    N_nodes_phi = N_nodes_phi_in;
    phi_nodes_tau = new fp_t [N_nodes_phi]-1;
    phi_nodes_phi = new fp_t [N_nodes_phi]-1;
  }
  else if (N_nodes_phi_in == 0){
    N_nodes_phi = 0;
    phi_nodes_tau = 0;
    phi_nodes_phi = 0;
  }
  else{
    printf("model::cannot create model with negative number of nodes in phi \n");
    exit(1);
  }
  response_to_parameters = 0; 
  N_depths = 41;
  tau_min = -5.0;
  tau_max = 1.0;
}

model::model(mcfg *cfg,io_class &io_in){

  N_depths = 0;
  N_nodes_temp = N_nodes_vt = N_nodes_vs = N_nodes_B = N_nodes_theta = N_nodes_phi = 0;
  temp_nodes_tau = temp_nodes_temp = 0;
  vt_nodes_tau = vt_nodes_vt = 0;
  vs_nodes_tau = vs_nodes_vs = 0;
  B_nodes_tau = B_nodes_B = 0;
  theta_nodes_tau = theta_nodes_theta = 0;
  phi_nodes_tau = phi_nodes_phi = 0;
  tau_min = cfg->tau_min;
  tau_max = cfg->tau_max;
  for (int p=0;p<cfg->np;++p){
    if (strcmp(cfg->par[p]->id,"TEMP") == 0){
      N_nodes_temp = cfg->par[p]->n;
      temp_nodes_tau = new fp_t [N_nodes_temp]-1;
      temp_nodes_temp = new fp_t [N_nodes_temp]-1;
      memcpy(temp_nodes_tau+1,cfg->par[p]->tau,N_nodes_temp*sizeof(fp_t));
      memcpy(temp_nodes_temp+1,cfg->par[p]->value,N_nodes_temp*sizeof(fp_t));
      temp_reg_type = cfg->par[p]->reg_type;
      temp_reg_alpha = cfg->par[p]->reg_alpha;
    }
    else if (strcmp(cfg->par[p]->id,"VT") == 0){
      N_nodes_vt = cfg->par[p]->n;
      vt_nodes_tau = new fp_t [N_nodes_vt]-1;
      vt_nodes_vt = new fp_t [N_nodes_vt]-1;
      memcpy(vt_nodes_tau+1,cfg->par[p]->tau,N_nodes_vt*sizeof(fp_t));
      memcpy(vt_nodes_vt+1,cfg->par[p]->value,N_nodes_vt*sizeof(fp_t));
      vt_reg_type = cfg->par[p]->reg_type;
      vt_reg_alpha = cfg->par[p]->reg_alpha;
    }
    else if (strcmp(cfg->par[p]->id,"VS") == 0){
      N_nodes_vs = cfg->par[p]->n;
      vs_nodes_tau = new fp_t [N_nodes_vs]-1;
      vs_nodes_vs = new fp_t [N_nodes_vs]-1;
      memcpy(vs_nodes_tau+1,cfg->par[p]->tau,N_nodes_vs*sizeof(fp_t));
      memcpy(vs_nodes_vs+1,cfg->par[p]->value,N_nodes_vs*sizeof(fp_t));
      vs_reg_type = cfg->par[p]->reg_type;
      vs_reg_alpha = cfg->par[p]->reg_alpha;
    }
    else if (strcmp(cfg->par[p]->id,"B") == 0){
      N_nodes_B = cfg->par[p]->n;
      B_nodes_tau = new fp_t [N_nodes_B]-1;
      B_nodes_B = new fp_t [N_nodes_B]-1;
      memcpy(B_nodes_tau+1,cfg->par[p]->tau,N_nodes_B*sizeof(fp_t));
      memcpy(B_nodes_B+1,cfg->par[p]->value,N_nodes_B*sizeof(fp_t));
      B_reg_type = cfg->par[p]->reg_type;
      B_reg_alpha = cfg->par[p]->reg_alpha;
    }
    else if (strcmp(cfg->par[p]->id,"THETA") == 0){
      N_nodes_theta = cfg->par[p]->n;
      theta_nodes_tau = new fp_t [N_nodes_theta]-1;
      theta_nodes_theta = new fp_t [N_nodes_theta]-1;
      memcpy(theta_nodes_tau+1,cfg->par[p]->tau,N_nodes_theta*sizeof(fp_t));
      memcpy(theta_nodes_theta+1,cfg->par[p]->value,N_nodes_theta*sizeof(fp_t));
      theta_reg_type = cfg->par[p]->reg_type;
      theta_reg_alpha = cfg->par[p]->reg_alpha;
    }
    else if (strcmp(cfg->par[p]->id,"PHI") == 0){
      N_nodes_phi = cfg->par[p]->n;
      phi_nodes_tau = new fp_t [N_nodes_phi]-1;
      phi_nodes_phi = new fp_t [N_nodes_phi]-1;
      memcpy(phi_nodes_tau+1,cfg->par[p]->tau,N_nodes_phi*sizeof(fp_t));
      memcpy(phi_nodes_phi+1,cfg->par[p]->value,N_nodes_phi*sizeof(fp_t));
      phi_reg_type = cfg->par[p]->reg_type;
      phi_reg_alpha = cfg->par[p]->reg_alpha;
    }
  }
  response_to_parameters=0;
}

model::model(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in){
  offs+=unpack(buf+offs,do_swap,io_in);
}

model::~model(void){

  if(temp_nodes_tau)
    delete[](temp_nodes_tau+1);
  if (temp_nodes_temp)
    delete[](temp_nodes_temp+1);
  if (vt_nodes_tau)
    delete[](vt_nodes_tau+1);
  if (vt_nodes_vt)
    delete[](vt_nodes_vt+1);
  if (vs_nodes_tau)
    delete[](vs_nodes_tau+1);
  if (vs_nodes_vs)
    delete[](vs_nodes_vs+1);
  if (B_nodes_tau)
    delete[](B_nodes_tau+1);
  if (B_nodes_B)
  	delete[](B_nodes_B+1);
  if (theta_nodes_tau)
    delete[](theta_nodes_tau+1);
  if (theta_nodes_theta)
    delete[](theta_nodes_theta+1);
  if (phi_nodes_tau)
    delete[](phi_nodes_tau+1);
  if (phi_nodes_phi)
    delete[](phi_nodes_phi+1);
  int N_parameters = get_N_nodes_total();
  if(response_to_parameters)
    del_ft3dim(response_to_parameters,1,N_parameters,1,7,1,N_depths);

}

int32_t model::size(io_class &io_in)
{
  int32_t sz=7*sizeof(int);
  sz += 6*sizeof(int); // regularizaton types
  sz += 6*sizeof(fp_t); // regularization alphas
  sz += 2*sizeof(fp_t); // tau_min, tau_max
  sz+=2*N_nodes_temp*sizeof(fp_t);
  sz+=2*N_nodes_vt*sizeof(fp_t);
  sz+=2*N_nodes_vs*sizeof(fp_t);
  sz+=2*N_nodes_B*sizeof(fp_t);
  sz+=2*N_nodes_theta*sizeof(fp_t);
  sz+=2*N_nodes_phi*sizeof(fp_t);

  return sz;
}

int32_t model::pack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{

  int32_t offs=::pack(buf,N_depths,do_swap);
  offs+=::pack(buf+offs,tau_min,do_swap);
  offs+=::pack(buf+offs,tau_max,do_swap);
  int P[]={N_nodes_temp,N_nodes_vt,N_nodes_vs,N_nodes_B,N_nodes_theta,N_nodes_phi,-1};
  for (int i=0;P[i]>-1;++i) offs+=::pack(buf+offs,P[i],do_swap);
  
  if (N_nodes_temp){
    offs+=::pack(buf+offs,temp_nodes_tau,1,N_nodes_temp,do_swap);
    offs+=::pack(buf+offs,temp_nodes_temp,1,N_nodes_temp,do_swap);
  }
  offs+=::pack(buf+offs,temp_reg_type,do_swap);
  offs+=::pack(buf+offs,temp_reg_alpha,do_swap);
  if (N_nodes_vt){
    offs+=::pack(buf+offs,vt_nodes_tau,1,N_nodes_vt,do_swap);
    offs+=::pack(buf+offs,vt_nodes_vt,1,N_nodes_vt,do_swap);
  }
  offs+=::pack(buf+offs,vt_reg_type,do_swap);
  offs+=::pack(buf+offs,vt_reg_alpha,do_swap);
  if (N_nodes_vs){
    offs+=::pack(buf+offs,vs_nodes_tau,1,N_nodes_vs,do_swap);
    offs+=::pack(buf+offs,vs_nodes_vs,1,N_nodes_vs,do_swap);
  }
  offs+=::pack(buf+offs,vs_reg_type,do_swap);
  offs+=::pack(buf+offs,vs_reg_alpha,do_swap);
  if (N_nodes_B){
    offs+=::pack(buf+offs,B_nodes_tau,1,N_nodes_B,do_swap);
    offs+=::pack(buf+offs,B_nodes_B,1,N_nodes_B,do_swap);
  }
  offs+=::pack(buf+offs,B_reg_type,do_swap);
  offs+=::pack(buf+offs,B_reg_alpha,do_swap);
  if (N_nodes_theta){
    offs+=::pack(buf+offs,theta_nodes_tau,1,N_nodes_theta,do_swap);
    offs+=::pack(buf+offs,theta_nodes_theta,1,N_nodes_theta,do_swap);
  }
  offs+=::pack(buf+offs,theta_reg_type,do_swap);
  offs+=::pack(buf+offs,theta_reg_alpha,do_swap);
  if (N_nodes_phi){
    offs+=::pack(buf+offs,phi_nodes_tau,1,N_nodes_phi,do_swap);
    offs+=::pack(buf+offs,phi_nodes_phi,1,N_nodes_phi,do_swap);
  }
  offs+=::pack(buf+offs,phi_reg_type,do_swap);
  offs+=::pack(buf+offs,phi_reg_alpha,do_swap);

  return offs;
}

int32_t model::unpack(uint08_t *buf,uint08_t do_swap,io_class &io_in)
{

  int32_t offs=::unpack(buf,N_depths,do_swap);
  offs+=::unpack(buf+offs,tau_min,do_swap);
  offs+=::unpack(buf+offs,tau_max,do_swap);
  offs+=::unpack(buf+offs,N_nodes_temp,do_swap);
  offs+=::unpack(buf+offs,N_nodes_vt,do_swap);
  offs+=::unpack(buf+offs,N_nodes_vs,do_swap);
  offs+=::unpack(buf+offs,N_nodes_B,do_swap);
  offs+=::unpack(buf+offs,N_nodes_theta,do_swap);
  offs+=::unpack(buf+offs,N_nodes_phi,do_swap);
  temp_nodes_tau = 0;
  temp_nodes_temp = 0;
  vt_nodes_tau = 0;
  vt_nodes_vt = 0;
  vs_nodes_tau = 0;
  vs_nodes_vs = 0;
  B_nodes_tau = 0;
  B_nodes_B = 0;
  theta_nodes_tau = 0;
  theta_nodes_theta = 0;
  phi_nodes_tau = 0;
  phi_nodes_phi = 0;
  if (N_nodes_temp){
    temp_nodes_tau = new fp_t [N_nodes_temp] -1;
    temp_nodes_temp = new fp_t [N_nodes_temp] -1;
    offs+=::unpack(buf+offs,temp_nodes_tau,1,N_nodes_temp,do_swap);
    offs+=::unpack(buf+offs,temp_nodes_temp,1,N_nodes_temp,do_swap);
  }
  offs+=::unpack(buf+offs,temp_reg_type,do_swap);
  offs+=::unpack(buf+offs,temp_reg_alpha,do_swap);
  if (N_nodes_vt){
    vt_nodes_tau = new fp_t [N_nodes_vt] -1;
    vt_nodes_vt = new fp_t [N_nodes_vt] -1;
    offs+=::unpack(buf+offs,vt_nodes_tau,1,N_nodes_vt,do_swap);
    offs+=::unpack(buf+offs,vt_nodes_vt,1,N_nodes_vt,do_swap);
  }
  offs+=::unpack(buf+offs,vt_reg_type,do_swap);
  offs+=::unpack(buf+offs,vt_reg_alpha,do_swap);
  if (N_nodes_vs){
    vs_nodes_tau = new fp_t [N_nodes_vs] -1;
    vs_nodes_vs = new fp_t [N_nodes_vs] -1;
    offs+=::unpack(buf+offs,vs_nodes_tau,1,N_nodes_vs,do_swap);
    offs+=::unpack(buf+offs,vs_nodes_vs,1,N_nodes_vs,do_swap);
  }
  offs+=::unpack(buf+offs,vs_reg_type,do_swap);
  offs+=::unpack(buf+offs,vs_reg_alpha,do_swap);
  if (N_nodes_B){
    B_nodes_tau = new fp_t [N_nodes_B] -1;
    B_nodes_B = new fp_t [N_nodes_B] -1;
    offs+=::unpack(buf+offs,B_nodes_tau,1,N_nodes_B,do_swap);
    offs+=::unpack(buf+offs,B_nodes_B,1,N_nodes_B,do_swap);
  }
  offs+=::unpack(buf+offs,B_reg_type,do_swap);
  offs+=::unpack(buf+offs,B_reg_alpha,do_swap);
  if (N_nodes_theta){
    theta_nodes_tau = new fp_t [N_nodes_theta] -1;
    theta_nodes_theta = new fp_t [N_nodes_theta] -1;
    offs+=::unpack(buf+offs,theta_nodes_tau,1,N_nodes_theta,do_swap);
    offs+=::unpack(buf+offs,theta_nodes_theta,1,N_nodes_theta,do_swap);
  }
  offs+=::unpack(buf+offs,theta_reg_type,do_swap);
  offs+=::unpack(buf+offs,theta_reg_alpha,do_swap);
  if (N_nodes_phi){
    phi_nodes_tau = new fp_t [N_nodes_phi] -1;
    phi_nodes_phi = new fp_t [N_nodes_phi] -1;
    offs+=::unpack(buf+offs,phi_nodes_tau,1,N_nodes_phi,do_swap);
    offs+=::unpack(buf+offs,phi_nodes_phi,1,N_nodes_phi,do_swap);
  }  
  offs+=::unpack(buf+offs,phi_reg_type,do_swap);
  offs+=::unpack(buf+offs,phi_reg_alpha,do_swap);
//
  response_to_parameters = 0;
  return offs;
}

int model::set_temp_nodes(fp_t * tau_in, fp_t * temp_in){

	for (int i=1;i<=N_nodes_temp;++i){
		temp_nodes_tau[i] = tau_in[i];
		temp_nodes_temp[i] = temp_in[i];
	}
  return 0;
}

int model::get_N_nodes_temp(){
  return N_nodes_temp;
}

fp_t * model::get_temp_nodes_tau(){

  fp_t * copy = new fp_t [N_nodes_temp] - 1;
  memcpy(copy+1,temp_nodes_tau+1,N_nodes_temp*sizeof(fp_t));

  return copy;
}

fp_t * model::get_temp_nodes_temp(){

  fp_t * copy = new fp_t [N_nodes_temp] - 1;
  memcpy(copy+1,temp_nodes_temp+1,N_nodes_temp*sizeof(fp_t));

  return copy;
}

int model::set_vt_nodes(fp_t * tau_in, fp_t * vt_in){
  for (int i=1;i<=N_nodes_vt;++i){
    vt_nodes_tau[i] = tau_in[i];
    vt_nodes_vt[i] = vt_in[i];
  }
  return 0;
}

int model::get_N_nodes_vt(){
  return N_nodes_vt;
}

fp_t * model::get_vt_nodes_tau(){

  fp_t * copy = new fp_t [N_nodes_vt] - 1;
  memcpy(copy+1,vt_nodes_tau+1,N_nodes_vt*sizeof(fp_t));

  return copy;
}

fp_t * model::get_vt_nodes_vt(){
  fp_t * copy = new fp_t [N_nodes_vt] - 1;
  memcpy(copy+1,vt_nodes_vt+1,N_nodes_vt*sizeof(fp_t));

  return copy;
}

int model::set_vs_nodes(fp_t * tau_in, fp_t * vs_in){
  for (int i=1;i<=N_nodes_vs;++i){
    vs_nodes_tau[i] = tau_in[i];
    vs_nodes_vs[i] = vs_in[i];
  }
  return 0;
}

int model::get_N_nodes_vs(){
  return N_nodes_vs;
}

fp_t * model::get_vs_nodes_tau(){
  fp_t * copy = new fp_t [N_nodes_vs] - 1;
  memcpy(copy+1,vs_nodes_tau+1,N_nodes_vs*sizeof(fp_t));
  return copy;
}

fp_t * model::get_vs_nodes_vs(){
  fp_t * copy = new fp_t [N_nodes_vs] - 1;
  memcpy(copy+1,vs_nodes_vs+1,N_nodes_vs*sizeof(fp_t));
  return copy;
}

int model::set_B_nodes(fp_t * tau_in, fp_t * B_in){
  for (int i=1;i<=N_nodes_B;++i){
    B_nodes_tau[i] = tau_in[i];
    B_nodes_B[i] = B_in[i];
  }
  return 0;
}

int model::get_N_nodes_B(){
  return N_nodes_B;
}

fp_t * model::get_B_nodes_tau(){
  fp_t * copy = new fp_t [N_nodes_B] - 1;
  memcpy(copy+1,B_nodes_tau+1,N_nodes_B*sizeof(fp_t));
  return copy;
}

fp_t * model::get_B_nodes_B(){
  fp_t * copy = new fp_t [N_nodes_B] - 1;
  memcpy(copy+1,B_nodes_B+1,N_nodes_B*sizeof(fp_t));
  return copy;
}

int model::set_theta_nodes(fp_t * tau_in, fp_t * theta_in){
  for (int i=1;i<=N_nodes_theta;++i){
    theta_nodes_tau[i] = tau_in[i];
    theta_nodes_theta[i] = theta_in[i];
  }
  return 0;
}

int model::get_N_nodes_theta(){
  return N_nodes_theta;
}

fp_t * model::get_theta_nodes_tau(){
  fp_t * copy = new fp_t [N_nodes_theta] - 1;
  memcpy(copy+1,theta_nodes_tau+1,N_nodes_theta*sizeof(fp_t));
  return copy;
}

fp_t * model::get_theta_nodes_theta(){
  fp_t * copy = new fp_t [N_nodes_theta] - 1;
  memcpy(copy+1,theta_nodes_theta+1,N_nodes_theta*sizeof(fp_t));
  return copy;
}

int model::set_phi_nodes(fp_t * tau_in, fp_t * phi_in){
  for (int i=1;i<=N_nodes_phi;++i){
    phi_nodes_tau[i] = tau_in[i];
    phi_nodes_phi[i] = phi_in[i];
  }
  return 0;
}

int model::get_N_nodes_phi(){
  return N_nodes_phi;
}

fp_t * model::get_phi_nodes_tau(){
  fp_t * copy = new fp_t [N_nodes_phi] - 1;
  memcpy(copy+1,phi_nodes_tau+1,N_nodes_phi*sizeof(fp_t));
  return copy;
}

fp_t * model::get_phi_nodes_phi(){
  fp_t * copy = new fp_t [N_nodes_phi] - 1;
  memcpy(copy+1,phi_nodes_phi+1,N_nodes_phi*sizeof(fp_t));
  return copy;
}


int model::get_N_nodes_total(){
  return N_nodes_temp + N_nodes_vt + N_nodes_vs + N_nodes_B + N_nodes_theta + N_nodes_phi;
}

fp_t model::get_parameter(int i){

  if (i>0 && i<=get_N_nodes_total()){

    if (i<=N_nodes_temp){
      return temp_nodes_temp[i];
    }
    else if (i<=N_nodes_temp + N_nodes_vt){
      return vt_nodes_vt[i-N_nodes_temp]; 
    }
    else if (i<=N_nodes_temp + N_nodes_vt + N_nodes_vs){
      return vs_nodes_vs[i-N_nodes_temp-N_nodes_vt];
    }
    else if (i<=N_nodes_temp+N_nodes_vt + N_nodes_vs + N_nodes_B){
      return B_nodes_B[i-N_nodes_temp-N_nodes_vt-N_nodes_vs];
    }
    else if (i<=N_nodes_temp+N_nodes_vt + N_nodes_vs + N_nodes_B + N_nodes_theta){
      return theta_nodes_theta[i-N_nodes_temp-N_nodes_vt-N_nodes_vs-N_nodes_B];
    }
    else{ 
      return phi_nodes_phi[i-N_nodes_temp-N_nodes_vt-N_nodes_vs-N_nodes_B-N_nodes_theta];
    }
  }
  printf("Invalid input. \n");
  return -1;
}

int model::perturb_node_value(int parameter_no, fp_t perturbation){

  // Here we perturb model parameter with given number:

  if (parameter_no <= N_nodes_temp){
    temp_nodes_temp[parameter_no] += perturbation;
    return 0;
  }
  else if (parameter_no <= N_nodes_temp + N_nodes_vt){
    vt_nodes_vt[parameter_no - N_nodes_temp] += perturbation;
    return 0;
  }
  else if (parameter_no <= N_nodes_temp + N_nodes_vt + N_nodes_vs){
    vs_nodes_vs[parameter_no - N_nodes_temp - N_nodes_vt] += perturbation;
    return 0;
  }
  else if (parameter_no <= N_nodes_temp + N_nodes_vt + N_nodes_vs + N_nodes_B){
    B_nodes_B[parameter_no - N_nodes_temp - N_nodes_vt - N_nodes_vs] += perturbation;
    return 0;
  }
  else if (parameter_no <= N_nodes_temp + N_nodes_vt + N_nodes_vs + N_nodes_B + N_nodes_theta){
    theta_nodes_theta[parameter_no - N_nodes_temp - N_nodes_vt - N_nodes_vs - N_nodes_B] += perturbation;
    return 0;
  }
  else if (parameter_no <= N_nodes_temp + N_nodes_vt + N_nodes_vs + N_nodes_B + N_nodes_theta + N_nodes_phi){
    phi_nodes_phi[parameter_no - N_nodes_temp - N_nodes_vt - N_nodes_vs - N_nodes_B - N_nodes_theta] += perturbation;
    return 0;
  }
   
  return 0;

}

int model::which_parameter(int i){

  if (i<=N_nodes_temp)
    return 1;
  else if (i<=N_nodes_temp+N_nodes_vt)
    return 2;
  else if (i<=N_nodes_temp+N_nodes_vt+N_nodes_vs)
    return 3;
  else if (i<=N_nodes_temp+N_nodes_vt+N_nodes_vs+N_nodes_B)
    return 4;
  else if (i<=N_nodes_temp+N_nodes_vt+N_nodes_vs+N_nodes_B+N_nodes_theta)
    return 5;
  else if (i<=N_nodes_temp+N_nodes_vt+N_nodes_vs+N_nodes_B+N_nodes_theta+N_nodes_phi)
    return 6;
  
  return 0;
}

int model::print(){

  if (N_nodes_temp){
    printf("atmos::model : temp. nodes: \n");
    for (int i=1;i<=N_nodes_temp;++i)
      printf("%d %5.15e %5.15e \n", i, temp_nodes_tau[i], temp_nodes_temp[i]);
  }
  else 
    printf("atmos::mode : no temp. nodes \n");

  if (N_nodes_vt){
    printf("atmos::model : vt nodes: \n");
    for (int i=1;i<=N_nodes_vt;++i)
      printf("%d %5.15e %5.15e \n", i, vt_nodes_tau[i], vt_nodes_vt[i]);
  }
  else 
    printf("atmos::mode : no vt nodes \n");


  if (N_nodes_vs){
    printf("atmos::model : vs nodes: \n");
    for (int i=1;i<=N_nodes_vs;++i)
      printf("%d %5.15e %5.15e \n", i, vs_nodes_tau[i], vs_nodes_vs[i]);
  }
  else 
    printf("atmos::mode : no vs nodes \n");

  if (N_nodes_B){
    printf("atmos::model : B nodes: \n");
    for (int i=1;i<=N_nodes_B;++i)
      printf("%d %5.15e %5.15e \n", i, B_nodes_tau[i], B_nodes_B[i]);
  }
  else 
    printf("atmos::mode : no B nodes \n");

  if (N_nodes_theta){
    printf("atmos::model : theta nodes: \n");
    for (int i=1;i<=N_nodes_theta;++i)
      printf("%d %5.15e %5.15e \n", i, theta_nodes_tau[i], theta_nodes_theta[i]* 180.0/pi);
  }
  else 
    printf("atmos::mode : no theta nodes \n");

  if (N_nodes_phi){
    printf("atmos::model : phi nodes: \n");
    for (int i=1;i<=N_nodes_phi;++i)
      printf("%d %5.15e %5.15e \n", i, phi_nodes_tau[i], phi_nodes_phi[i]* 180.0/pi);
  }
  else 
    printf("atmos::mode : no phi nodes \n");

  return 0;
}

int model::print_to_file(FILE * output){

  if (N_nodes_temp){
    fprintf(output, "atmos::model : temp. nodes: \n");
    for (int i=1;i<=N_nodes_temp;++i)
      fprintf(output, "%d %5.15e %5.15e \n", i, temp_nodes_tau[i], temp_nodes_temp[i]);
  }
  else 
    fprintf(output, "atmos::mode : no temp. nodes \n");

  if (N_nodes_vt){
    fprintf(output, "atmos::model : vt nodes: \n");
    for (int i=1;i<=N_nodes_vt;++i)
      fprintf(output, "%d %5.15e %5.15e \n", i, vt_nodes_tau[i], vt_nodes_vt[i]);
  }
  else 
    printf("atmos::mode : no vt nodes \n");


  if (N_nodes_vs){
    fprintf(output,"atmos::model : vs nodes: \n");
    for (int i=1;i<=N_nodes_vs;++i)
      fprintf(output,"%d %5.15e %5.15e \n", i, vs_nodes_tau[i], vs_nodes_vs[i]);
  }
  else 
    fprintf(output,"atmos::mode : no vs nodes \n");

  if (N_nodes_B){
    fprintf(output,"atmos::model : B nodes: \n");
    for (int i=1;i<=N_nodes_B;++i)
      fprintf(output,"%d %5.15e %5.15e \n", i, B_nodes_tau[i], B_nodes_B[i]);
  }
  else 
    fprintf(output,"atmos::mode : no B nodes \n");

  if (N_nodes_theta){
    fprintf(output,"atmos::model : theta nodes: \n");
    for (int i=1;i<=N_nodes_theta;++i)
      fprintf(output,"%d %5.15e %5.15e \n", i, theta_nodes_tau[i], theta_nodes_theta[i] * 180.0/pi);
  }
  else 
    fprintf(output,"atmos::mode : no theta nodes \n");

  if (N_nodes_phi){
    fprintf(output,"atmos::model : phi nodes: \n");
    for (int i=1;i<=N_nodes_phi;++i)
      fprintf(output,"%d %5.15e %5.15e \n", i, phi_nodes_tau[i], phi_nodes_phi[i] * 180.0/pi);
  }
  else 
    fprintf(output,"atmos::mode : no phi nodes \n");

  return 0;
}

int model::set_response_to_parameters(fp_t *** responses_in, int Nd_in){

  N_depths = Nd_in;
  int N_parameters = N_nodes_temp + N_nodes_vt + N_nodes_vs + N_nodes_B + N_nodes_theta + N_nodes_phi;
  if (response_to_parameters==0){
    response_to_parameters = ft3dim(1,N_parameters,1,7,1,Nd_in);
  }
  memset(response_to_parameters[1][1]+1,0,N_parameters*7*N_depths*sizeof(fp_t));

  for (int p=1;p<=N_parameters;++p)
    for (int q=1;q<=7;++q)
      for (int x3k=1;x3k<=N_depths;++x3k)
        response_to_parameters[p][q][x3k] = responses_in[p][q][x3k];
  return 0;
}

int model::cpy_values_from(model* input){

  tau_min = input->get_tau_min();
  tau_max = input->get_tau_max();
  fp_t * tau, * value;
  tau = input->get_temp_nodes_tau();
  value = input->get_temp_nodes_temp();
  set_temp_nodes(tau,value);
  delete[](tau+1);delete[](value+1);
  tau = input->get_vt_nodes_tau();
  value = input->get_vt_nodes_vt();
  set_vt_nodes(tau,value);
  delete[](tau+1);delete[](value+1);
  tau = input->get_vs_nodes_tau();
  value = input->get_vs_nodes_vs();
  set_vs_nodes(tau,value);
  delete[](tau+1);delete[](value+1);
  tau = input->get_B_nodes_tau();
  value = input->get_B_nodes_B();
  set_B_nodes(tau,value);
  delete[](tau+1);delete[](value+1);
  tau = input->get_theta_nodes_tau();
  value = input->get_theta_nodes_theta();
  set_theta_nodes(tau,value);
  delete[](tau+1);delete[](value+1);
  tau = input->get_phi_nodes_tau();
  value = input->get_phi_nodes_phi();
  set_phi_nodes(tau,value);
  delete[](tau+1);delete[](value+1);
  return 0;
}

fp_t *** model::get_response_to_parameters(){

  int N_parameters = N_nodes_temp + N_nodes_vt + N_nodes_vs + N_nodes_B + N_nodes_theta + N_nodes_phi;
  fp_t *** response_to_parameters_cpy = ft3dim(1,N_parameters,1,7,1,N_depths);
  memcpy(response_to_parameters_cpy[1][1]+1,response_to_parameters[1][1]+1,N_parameters*7*N_depths*sizeof(fp_t));
  return response_to_parameters_cpy;
}

void model::set_temp_reg(int type, fp_t str){
  temp_reg_type = type;
  temp_reg_alpha = str;
}
void model::set_vt_reg(int type, fp_t str){
  vt_reg_type = type;
  vt_reg_alpha = str;
}
void model::set_vs_reg(int type, fp_t str){
  vs_reg_type = type;
  vs_reg_alpha = str;
}
void model::set_B_reg(int type, fp_t str){
  B_reg_type = type;
  B_reg_alpha = str;
}
void model::set_theta_reg(int type, fp_t str){
  theta_reg_type = type;
  theta_reg_alpha = str;
}
void model::set_phi_reg(int type, fp_t str){
  phi_reg_type = type;
  phi_reg_alpha = str;
}

int model::correct(fp_t * correction){

  // This function applies an array of corrections to the model.
  // We need to take care that corrections and not too severe. Which means very high or very low temperatures
  // Negative magnetic field, negative microturbulent velocities or phi and theta which are in some unreasonable range.
  int N_parameters = N_nodes_temp+N_nodes_vt+N_nodes_vs+N_nodes_B+N_nodes_theta+N_nodes_phi;

  // We need to backup node values so that we can properly set-up "correction" array, so that we can
  // later, if necessary, "revert" to the old model. Is this most convenient? Probably not, but it 
  // is the fastest, and does not require writing new models and things like that.

  // First make sure that corrections are not too big:
  /*for (int i=1;i<=N_nodes_temp;++i)
    if (fabs(correction[i]) > 1000.0) correction[i] = 1000.0 * fabs(correction[i])/correction[i];*/
  
  /*int vt_start = N_nodes_temp+1;
  int vt_stop = N_nodes_temp+N_nodes_vt;
  for (int i=vt_start;i<=vt_stop;++i)
    if (fabs(correction[i]) > 3E5) correction[i] = 3E5 * fabs(correction[i])/correction[i];

  int vs_start = N_nodes_temp+N_nodes_vt+1;
  int vs_stop = N_nodes_temp+N_nodes_vt+N_nodes_vs;
  for (int i=vs_start; i<=vs_stop; ++i)
    if (fabs(correction[i]) > 3E5) correction[i] = 3E5 * fabs(correction[i])/correction[i];*/
    
  /*int B_start = N_nodes_temp+N_nodes_vt+N_nodes_vs+1;
  int B_stop  = B_start-1+N_nodes_B;

  for (int i=B_start;i<=B_stop;++i)
    if (fabs(correction[i]) > 1000.0) correction[i] = 1000.0 * fabs(correction[i])/correction[i];
  */
  for (int i=1;i<=N_parameters;++i)
    perturb_node_value(i,correction[i]);
  /*
  // Check temperatures:
  for (int i=1;i<=N_nodes_temp;++i){
    if (temp_nodes_temp[i] < 3400.0) temp_nodes_temp[i] = 3400.0; // Lowest possible T
    if (temp_nodes_temp[i] > 10000.0) temp_nodes_temp[i] = 20000.0; // Highest possible T (?)
                                                                    // These quantities are very questionable
  }
  for (int i=1;i<=N_nodes_vt;++i){
    if (vt_nodes_vt[i] < 0) vt_nodes_vt[i] *= (-1.0);
    if (vt_nodes_vt[i] > 20E5) vt_nodes_vt[i] = 20E5; // highest possible vt
  }
  
  for (int i=1;i<=N_nodes_vs;++i)
    if (fabs(vs_nodes_vs[i])>20E5) vs_nodes_vs[i] = 20E5*fabs(vs_nodes_vs[i])/vs_nodes_vs[i];
  
  for (int i=1;i<=N_nodes_B;++i){
    if (B_nodes_B[i] < 0.0) B_nodes_B[i] = fabs(B_nodes_B[i]);
    if (B_nodes_B[i] > 10000.0) B_nodes_B[i] = 10000.0; // Highest possible B
  }
  */
  return 0;
}

int model::bracket_parameter_values(){

  for (int i=1;i<=N_nodes_temp;++i){
    if (temp_nodes_temp[i]<T_MIN) temp_nodes_temp[i] = T_MIN;
    if (temp_nodes_temp[i]>T_MAX) temp_nodes_temp[i] = T_MAX;
  }
  for (int i=1;i<=N_nodes_vt;++i){
    if (vt_nodes_vt[i]<VT_MIN) vt_nodes_vt[i] = VT_MIN;
    if (vt_nodes_vt[i]>VT_MAX) vt_nodes_vt[i] = VT_MAX;
  }
  for (int i=1;i<=N_nodes_vs;++i){
    if (vs_nodes_vs[i]<VS_MIN) vs_nodes_vs[i] = VS_MIN;
    if (vs_nodes_vs[i]>VS_MAX) vs_nodes_vs[i] = VS_MAX;
  }
  for (int i=1;i<=N_nodes_B;++i){
    if (B_nodes_B[i]<B_MIN) B_nodes_B[i] = B_MIN;
    if (B_nodes_B[i]>B_MAX) B_nodes_B[i] = B_MAX;
  }

  return 0;
}

int model::set_parameters(fp_t * par_input){

  // This functions sets the values of the parameters to the values specified in the input array
  int N_parameters = N_nodes_temp+N_nodes_vt+N_nodes_vs+N_nodes_B+N_nodes_theta+N_nodes_phi;
  
  for (int i=1;i<=N_parameters;++i){
    fp_t old = get_parameter(i);
    fp_t difference = par_input[i]-old;
    perturb_node_value(i,difference);
  }
  return 0;
}

int model::polish_angles(){
  for (int i=1;i<=N_nodes_theta;++i){
    fp_t temp = cos(theta_nodes_theta[i]);
    theta_nodes_theta[i] = acos(temp);
  }
  for (int i=1;i<=N_nodes_phi;++i){
    fp_t temp = tan(phi_nodes_phi[i]);
    phi_nodes_phi[i] = atan(temp);
  }
  return 0;
}

model * model_new(int N_nodes_temp_in, int N_nodes_vt_in, int N_nodes_vs_in, int N_nodes_B_in){
  return new model(N_nodes_temp_in, N_nodes_vt_in, N_nodes_vs_in, N_nodes_B_in);
}

model * model_new(int N_nodes_temp_in, int N_nodes_vt_in, int N_nodes_vs_in, int N_nodes_B_in, int N_nodes_theta_in, int N_nodes_phi_in){
 return new model(N_nodes_temp_in, N_nodes_vt_in, N_nodes_vs_in, N_nodes_B_in, N_nodes_theta_in, N_nodes_phi_in); 
}

model * model_new(mcfg *cfg,io_class &io_in){
  return new model(cfg,io_in);
}

model* model_new(uint08_t *buf,int32_t &offs,uint08_t do_swap,io_class &io_in){
  return new model(buf,offs,do_swap,io_in);
}

model * clone(model * input){

  model * output;
  int N_nodes_temp = input->get_N_nodes_temp();
  int N_nodes_vt = input->get_N_nodes_vt();
  int N_nodes_vs = input->get_N_nodes_vs();
  int N_nodes_B = input->get_N_nodes_B();
  int N_nodes_theta = input->get_N_nodes_theta();
  int N_nodes_phi = input->get_N_nodes_phi();

  output = new model(N_nodes_temp,N_nodes_vt,N_nodes_vs,N_nodes_B,N_nodes_theta,N_nodes_phi);
  output->set_tau_min(input->get_tau_min());
  output->set_tau_max(input->get_tau_max());

  fp_t * tau, * value;
  tau = input->get_temp_nodes_tau();
  value = input->get_temp_nodes_temp();
  output->set_temp_nodes(tau,value);
  output->set_temp_reg(input->get_temp_reg_type(),input->get_temp_reg_alpha());
  delete[](tau+1);delete[](value+1);
  tau = input->get_vt_nodes_tau();
  value = input->get_vt_nodes_vt();
  output->set_vt_nodes(tau,value);
  output->set_vt_reg(input->get_vt_reg_type(),input->get_vt_reg_alpha());
  delete[](tau+1);delete[](value+1);
  tau = input->get_vs_nodes_tau();
  value = input->get_vs_nodes_vs();
  output->set_vs_nodes(tau,value);
  output->set_vs_reg(input->get_vs_reg_type(),input->get_vs_reg_alpha());
  delete[](tau+1);delete[](value+1);
  tau = input->get_B_nodes_tau();
  value = input->get_B_nodes_B();
  output->set_B_nodes(tau,value);
  output->set_B_reg(input->get_B_reg_type(),input->get_B_reg_alpha());
  delete[](tau+1);delete[](value+1);
  tau = input->get_theta_nodes_tau();
  value = input->get_theta_nodes_theta();
  output->set_theta_nodes(tau,value);
  output->set_theta_reg(input->get_theta_reg_type(),input->get_theta_reg_alpha());
  delete[](tau+1);delete[](value+1);
  tau = input->get_phi_nodes_tau();
  value = input->get_phi_nodes_phi();
  output->set_phi_nodes(tau,value);
  output->set_phi_reg(input->get_phi_reg_type(),input->get_phi_reg_alpha());
  delete[](tau+1);delete[](value+1);

  return output;
}

// --------------------------------------------------------------------------------------------------------------

// Modelcube related properties:
modelcube::modelcube(){
  N_nodes_temp=N_nodes_vt=N_nodes_vs=N_nodes_B=N_nodes_theta=N_nodes_phi=0;
  nx=ny=0;
  data = 0;
  temp_nodes_tau = 0;
  vt_nodes_tau = 0;
  vs_nodes_tau = 0;
  theta_nodes_tau = 0;
  phi_nodes_tau = 0;
}

modelcube::modelcube(model * make_from, int nx_in, int ny_in){

  // Initialize
  N_nodes_temp=N_nodes_vt=N_nodes_vs=N_nodes_B=N_nodes_theta=N_nodes_phi=0;
  nx=ny=0;
  data = 0;
  temp_nodes_tau = 0;
  vt_nodes_tau = 0;
  vs_nodes_tau = 0;
  B_nodes_tau=0;
  theta_nodes_tau = 0;
  phi_nodes_tau = 0;
  // Set proper values (if we have them)
  N_nodes_temp=make_from->get_N_nodes_temp();
  N_nodes_vt=make_from->get_N_nodes_vt();
  N_nodes_vs=make_from->get_N_nodes_vs();
  N_nodes_B=make_from->get_N_nodes_B();
  N_nodes_theta=make_from->get_N_nodes_theta();
  N_nodes_phi=make_from->get_N_nodes_phi();

  N_parameters = make_from->get_N_nodes_total();

  fp_t * temp; // temporary array
  if (N_nodes_temp){
    temp = make_from->get_temp_nodes_tau();
    temp_nodes_tau = new fp_t [N_nodes_temp] - 1;
    for (int i=1;i<=N_nodes_temp;++i)
      temp_nodes_tau[i] = temp[i];
    delete[](temp+1);
  }
  if (N_nodes_vt){
    temp = make_from->get_vt_nodes_tau();
    vt_nodes_tau = new fp_t [N_nodes_vt] - 1;
    for (int i=1;i<=N_nodes_vt;++i)
      vt_nodes_tau[i] = temp[i];
    delete[](temp+1);
  }
  if (N_nodes_vs){
    temp = make_from->get_vs_nodes_tau();
    vs_nodes_tau = new fp_t [N_nodes_vs] - 1;
    for (int i=1;i<=N_nodes_vs;++i)
      vs_nodes_tau[i] = temp[i];
    delete[](temp+1);
  }
  if (N_nodes_B){
    temp = make_from->get_B_nodes_tau();
    B_nodes_tau = new fp_t [N_nodes_B] - 1;
    for (int i=1;i<=N_nodes_B;++i)
      B_nodes_tau[i] = temp[i];
    delete[](temp+1);
  }
  if (N_nodes_theta){
    temp = make_from->get_theta_nodes_tau();
    theta_nodes_tau = new fp_t [N_nodes_theta] - 1;
    for (int i=1;i<=N_nodes_theta;++i)
      theta_nodes_tau[i] = temp[i];
    delete[](temp+1);
  }
  if (N_nodes_phi){
    temp = make_from->get_vt_nodes_tau();
    phi_nodes_tau = new fp_t [N_nodes_phi] - 1;
    for (int i=1;i<=N_nodes_phi;++i)
      phi_nodes_tau[i] = temp[i];
    delete[](temp+1);
  }

  nx = nx_in;ny=ny_in;
  data=ft3dim(1,nx,1,ny,1,N_parameters);
  memset(data[1][1]+1,0,nx*ny*N_parameters*sizeof(fp_t));

}

modelcube::~modelcube(){
  if (temp_nodes_tau)
    delete[](temp_nodes_tau+1);
  if (vt_nodes_tau)
    delete[](vt_nodes_tau+1);
  if (vs_nodes_tau)
    delete[](vs_nodes_tau+1);
  if (B_nodes_tau)
    delete[](B_nodes_tau+1);
  if (theta_nodes_tau)
    delete[](theta_nodes_tau+1);
  if (phi_nodes_tau)
    delete[](phi_nodes_tau+1);
  if (data)
    del_ft3dim(data,1,nx,1,ny,1,N_parameters);

}

void modelcube::add_model(model * to_add, int i, int j){

  for (int kk=1;kk<=N_parameters;++kk){
    data[i][j][kk] = to_add->get_parameter(kk);
  }

}

void modelcube::simple_print(const char* printhere){

  FILE * output;
  output = fopen(printhere,"w");
  for (int i=1;i<=nx;++i)
    for (int j=1;j<=ny;++j){
      fprintf(output,"%d %d ", i, j);
      for (int kk=1;kk<=N_parameters;++kk)
        fprintf(output,"%f ", data[i][j][kk]);
      fprintf(output, "\n");
    }
  fclose(output);
}

fp_t *** modelcube::get_data(int &nx_in, int &ny_in, int &np_in){
  // Returns the copy of the datacube, which we write down then. 
  // This might get us in trouble with memory with larger cubes, but datacube is private so we 
  // cannot really return pointer, can we. Maybe we can and we should, not sure though.

  fp_t *** data_copy = ft3dim(1,nx,1,ny,1,N_parameters);
  memcpy(data_copy[1][1]+1,data[1][1]+1,nx*ny*N_parameters*sizeof(fp_t));
  nx_in = nx;
  ny_in = ny;
  np_in = N_parameters;

  return data_copy;
}

void modelcube::set_data(fp_t *** data_in){
  memcpy(data[1][1]+1,data_in[1][1]+1,nx*ny*N_parameters*sizeof(fp_t));
}