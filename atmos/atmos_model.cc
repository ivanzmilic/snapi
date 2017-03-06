#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "types.h"
#include "io.h"
#include "mem.h"
#include "atmos.h"
#include "profiles.h"
#include "const.h"

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

int model::get_parameter(int i){

  if (i>0 && i<=get_N_nodes_total()){

    if (i<=N_nodes_temp)
      return temp_nodes_temp[i];
    else if (i<=N_nodes_temp + N_nodes_vt)
      return vt_nodes_vt[i-N_nodes_temp]; 
    else if (i<=N_nodes_temp + N_nodes_vt + N_nodes_vs)
      return vs_nodes_vs[i-N_nodes_temp-N_nodes_vt];
    else if (i<=N_nodes_temp+N_nodes_vt + N_nodes_vs + N_nodes_B)
      return B_nodes_B[i-N_nodes_temp-N_nodes_vt-N_nodes_vs];
    else if (i<=N_nodes_temp+N_nodes_vt + N_nodes_vs + N_nodes_B + N_nodes_theta)
      return theta_nodes_theta[i-N_nodes_temp-N_nodes_vt-N_nodes_vs-N_nodes_B];
    else return phi_nodes_phi[i-N_nodes_temp-N_nodes_vt-N_nodes_vs-N_nodes_B-N_nodes_theta];
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

model * model_new(int N_nodes_temp_in, int N_nodes_vt_in, int N_nodes_vs_in, int N_nodes_B_in){
  return new model(N_nodes_temp_in, N_nodes_vt_in, N_nodes_vs_in, N_nodes_B_in);
}

model * model_new(int N_nodes_temp_in, int N_nodes_vt_in, int N_nodes_vs_in, int N_nodes_B_in, int N_nodes_theta_in, int N_nodes_phi_in){
 return new model(N_nodes_temp_in, N_nodes_vt_in, N_nodes_vs_in, N_nodes_B_in, N_nodes_theta_in, N_nodes_phi_in); 
}




