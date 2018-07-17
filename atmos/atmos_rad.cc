#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "const.h"
#include "types.h"
#include "io.h"
#include "mathtools.h"
#include "mem.h"

#include "atmol/atmol.h"

#include "atmos.h"

fp_t ***atmosphere::project(fp_t ***v1,fp_t ***v2,fp_t ***v3,fp_t theta,fp_t phi,int32_t ll1,int32_t ul1,int32_t ll2,int32_t ul2,int32_t ll3,int32_t ul3)
{
  int32_t nn=(ul1-ll1+1)*(ul2-ll2+1)*(ul3-ll3+1);
  fp_t ***vp=ft3dim(ll1,ul1,ll2,ul2,ll3,ul3);
  fp_t *p=vp[ll1][ll2]+ll3,*p1=v1[ll1][ll2]+ll3,*p2=v2[ll1][ll2]+ll3,*p3=v3[ll1][ll2]+ll3;
  fp_t w1=sin(phi)*sin(pi-theta),w2=cos(phi)*sin(pi-theta),w3=cos(pi-theta);
  for(int32_t i=0;i<nn;++i) p[i]=w1*p1[i]+w2*p2[i]+w3*p3[i];
  return vp;
}
  
fp_t ****atmosphere::transform(fp_t ***v1,fp_t ***v2,fp_t ***v3,fp_t theta,fp_t phi,int32_t ll1,int32_t ul1,int32_t ll2,int32_t ul2,int32_t ll3,int32_t ul3)
{
  int32_t nn=(ul1-ll1+1)*(ul2-ll2+1)*(ul3-ll3+1);
  fp_t ****vp=ft4dim(1,3,ll1,ul1,ll2,ul2,ll3,ul3);
  fp_t *p1=v1[ll1][ll2]+ll3,*p2=v2[ll1][ll2]+ll3,*p3=v3[ll1][ll2]+ll3;
  fp_t *mg=vp[1][ll1][ll2]+ll3,*in=vp[2][ll1][ll2]+ll3,*az=vp[3][ll1][ll2]+ll3;
  
  // Keep in mind that theta is actually pi-theta which we uaually use
  fp_t cp=cos(phi),sp=sin(phi),ct=cos(pi-theta),st=sin(pi-theta);
  
  for(int32_t i=0;i<nn;++i){
    mg[i]=sqrt(p1[i]*p1[i]+p2[i]*p2[i]+p3[i]*p3[i]);

    // Lets apply the rotation per definition:

    fp_t newB_x = ct*cp*p1[i] + sp*p2[i] - st*cp*p3[i];
    fp_t newB_y = -ct*sp*p1[i] + cp*p2[i] + st*sp*p3[i];
    fp_t newB_z = st*p1[i] + ct*p3[i];
    
    in[i] = acos(newB_z/mg[i]) * 180.0/pi;
    az[i] = atan(newB_y/newB_x)* 180.0/pi;

    //in[i]=90.0-atan2(tmp2,sqrt(tmp1*tmp1+tmp3*tmp3))*180.0/pi; // 0=vertical up
    //az[i]=atan2(tmp3,tmp1)*180.0/pi; // !! counter-clockwise=positive, as seen by the observer
  }
  return vp;
}

void atmosphere::transform_responses(fp_t **** responses, fp_t theta, fp_t phi, int l_low, int l_high){

  // So we now have the responses in the frame of the ray. We want to turn them back to the responses for theta and phi in the frame of the atmosphere. 

  fp_t cp=cos(phi),sp=sin(phi),ct=cos(pi-theta),st=sin(pi-theta);

  for (int x3k=x3l;x3k<=x3h;++x3k){ // For each depth

    fp_t newB_x = ct*cp*Bx[x1l][x2l][x3k] + sp*By[x1l][x2l][x3k] - st*cp*Bz[x1l][x2l][x3k];
    fp_t newB_y = -ct*sp*Bx[x1l][x2l][x3k] + cp*By[x1l][x2l][x3k] + st*sp*Bz[x1l][x2l][x3k];
    fp_t newB_z = st*Bx[x1l][x2l][x3k] + ct*Bz[x1l][x2l][x3k];
    fp_t B_mag = sqrt(Bx[x1l][x2l][x3k]*Bx[x1l][x2l][x3k] + By[x1l][x2l][x3k]*By[x1l][x2l][x3k] + Bz[x1l][x2l][x3k]*Bz[x1l][x2l][x3k]);
    
    fp_t d_inc_d_Bx = (st*B_mag - Bx[x1l][x2l][x3k]/B_mag*newB_z)/B_mag/B_mag * (-1.0)/(sqrt(1.0 - newB_z*newB_z/B_mag/B_mag));
    fp_t d_inc_d_By = (-By[x1l][x2l][x3k]/B_mag*newB_z)/B_mag/B_mag * (-1.0)/(sqrt(1.0 - newB_z*newB_z/B_mag/B_mag));
    fp_t d_inc_d_Bz = (ct*B_mag - Bz[x1l][x2l][x3k]/B_mag*newB_z)/B_mag/B_mag * (-1.0)/(sqrt(1.0 - newB_z*newB_z/B_mag/B_mag));

    fp_t d_az_d_Bx = (-ct*sp*newB_x - ct*cp*newB_y)/newB_x/newB_x / (1.0 + newB_y*newB_y/newB_x/newB_x);
    fp_t d_az_d_By = (cp*newB_x - sp*newB_y)/newB_x/newB_x / (1.0 + newB_y*newB_y/newB_x/newB_x);
    fp_t d_az_d_Bz = (st*sp*newB_x + st*cp*newB_y)/newB_x/newB_x / (1.0 + newB_y*newB_y/newB_x/newB_x);

    fp_t theta_B = acos(Bz[x1l][x2l][x3k]/B_mag);
    fp_t phi_B = atan(By[x1l][x2l][x3k]/Bx[x1l][x2l][x3k]); // These are in the atmosphere frame of the reference

    fp_t dBx_d_theta_B = B_mag*cos(theta_B)*cos(phi_B);
    fp_t dBx_d_phi_B = -B_mag*sin(theta_B)*sin(phi_B);
    fp_t dBy_d_theta_B = B_mag*cos(theta_B)*sin(phi_B);
    fp_t dBy_d_phi_B = B_mag*sin(theta_B)*cos(phi_B);
    fp_t dBz_d_theta_B = -B_mag*sin(theta_B);
    fp_t dBz_d_phi_B = 0.0;

    fp_t d_inc_d_theta_B = d_inc_d_Bx*dBx_d_theta_B + d_inc_d_By*dBy_d_theta_B + d_inc_d_Bz*dBz_d_theta_B;
    fp_t d_inc_d_phi_B = d_inc_d_Bx*dBx_d_phi_B + d_inc_d_By*dBy_d_phi_B + d_inc_d_Bz*dBz_d_phi_B;
    fp_t d_az_d_theta_B = d_az_d_Bx*dBx_d_theta_B + d_az_d_By*dBy_d_theta_B + d_az_d_Bz*dBz_d_theta_B;
    fp_t d_az_d_phi_B = d_az_d_Bx*dBx_d_phi_B + d_az_d_By*dBy_d_phi_B + d_az_d_Bz*dBz_d_phi_B;

    for (int l=l_low;l<=l_high;++l)
      for (int s=1;s<=4;++s){

        fp_t to_in = responses[6][x3k][l][s];
        fp_t to_az = responses[7][x3k][l][s];

        responses[4][x3k][l][s] *= cp; // While we are there, this one is also important. 

        responses[6][x3k][l][s] = to_in * d_inc_d_theta_B + to_az * d_az_d_theta_B;
        responses[7][x3k][l][s] = to_in * d_inc_d_phi_B + to_az * d_az_d_phi_B;
      }
  }  
}

//
// frequency-by-frequency versions
//

fp_t ***atmosphere::opacity(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in,
                            fp_t ****B,fp_t theta,fp_t phi,fp_t lambda)
{
  fp_t ***op=thomson_sc(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h); // electron scattering
  //memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  for(int a=0;a<natm;++a) op=add(atml[a]->opacity(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op,x1l,x1h,x2l,x2h,x3l,x3h);
  return op;
}

fp_t ***atmosphere::emissivity(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in,fp_t ****B,fp_t theta,fp_t phi,fp_t lambda)
{
  fp_t ***em=thomson_em(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h); // electron scattering
  //memset(em[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  for(int a=0;a<natm;++a) em=add(atml[a]->emissivity(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),em,x1l,x1h,x2l,x2h,x3l,x3h); 
  return em;
}

void atmosphere::normalize_to_referent_opacity(fp_t *** op, fp_t *** em){

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        op[x1i][x2i][x3i] /= op_referent[x1i][x2i][x3i];
        em[x1i][x2i][x3i] /= op_referent[x1i][x2i][x3i];
      }
}

void atmosphere::normalize_to_referent_opacity(fp_t ***** op, fp_t **** em){

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        for (int s=1;s<=4;++s){
          em[x1i][x2i][x3i][s] /= op_referent[x1i][x2i][x3i];
          for (int sp=1;sp<=4;++sp)
            op[x1i][x2i][x3i][s][sp] /= op_referent[x1i][x2i][x3i];
      }
}

void atmosphere::normalize_to_referent_opacity(fp_t ***** op, fp_t **** em,fp_t ******* op_pert, fp_t ****** em_pert){

  for (int p=1;p<=7;++p)
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i)
            for (int s=1;s<=4;++s){
              em_pert[p][x3k][x1i][x2i][x3i][s] = em_pert[p][x3k][x1i][x2i][x3i][s]/op_referent[x1i][x2i][x3i] 
                -em[x1i][x2i][x3i][s] * op_referent_derivative[p][x3k][x1i][x2i][x3i] / op_referent[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i];   
              for (int sp=1;sp<=4;++sp)
                op_pert[p][x3k][x1i][x2i][x3i][s][sp] = op_pert[p][x3k][x1i][x2i][x3i][s][sp]/op_referent[x1i][x2i][x3i]
                  -op[x1i][x2i][x3i][s][sp] * op_referent_derivative[p][x3k][x1i][x2i][x3i] / op_referent[x1i][x2i][x3i] / op_referent[x1i][x2i][x3i];
  }           
}


void atmosphere::de_normalize(fp_t *** op, fp_t *** em){
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        op[x1i][x2i][x3i] *= op_referent[x1i][x2i][x3i];
        em[x1i][x2i][x3i] *= op_referent[x1i][x2i][x3i];
      }
}

void atmosphere::de_normalize(fp_t ***** op, fp_t **** em){

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        for (int s=1;s<=4;++s){
          em[x1i][x2i][x3i][s] *= op_referent[x1i][x2i][x3i];
          for (int sp=1;sp<=4;++sp)
            op[x1i][x2i][x3i][s][sp] *= op_referent[x1i][x2i][x3i];
      }
}

fp_t ***atmosphere::opacity_lte(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in,
                            fp_t ****B,fp_t theta,fp_t phi,fp_t lambda)
{
  fp_t ***op=thomson_sc(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h); // electron scattering
  memset(op[x1l][x2l]+x3l,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  for(int a=0;a<natm;++a)
    if (!atml[a]->check_if_nlte()) 
      op=add(atml[a]->opacity(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op,x1l,x1h,x2l,x2h,x3l,x3h);
  return op;
}


fp_t atmosphere::opacity_continuum(fp_t T_in, fp_t Ne_in, fp_t lambda, int x1i, int x2i, int x3i){
  // An ad-hoc function which computes continuum opacity for given point, and given frequency
  fp_t op_cont = 6.65E-25 * Ne_in;
  // And then add f-f and b-f sources from all other atoms and molecules
  for (int a=0;a<natm;++a)
    if (!atml[a]->check_if_nlte()) 
      op_cont += atml[a]->opacity_continuum(T_in, Ne_in, lambda, x1i,x2i,x3i);
  return op_cont;
}

fp_t ** atmosphere::opacity_continuum_derivative(fp_t T_in, fp_t Ne_in, fp_t lambda, int x1i, int x2i, int x3i){

  fp_t ** op_cont_pert = ft2dim(1,7,x3l,x3h);
  memset(op_cont_pert[1]+x3l,0,7*(x3h-x3l+1)*sizeof(fp_t));
  for (int p=1;p<=7;++p)
      op_cont_pert[p][x3i] = Ne_lte_der[p][x1i][x2i][x3i] * 6.65E-25;
  
  for (int a=0;a<natm;++a)
    if (!atml[a]->check_if_nlte()) 
      op_cont_pert = add(atml[a]->opacity_continuum_pert(T_in, Ne_in, lambda, x1i,x2i,x3i), op_cont_pert,1,7,x3l,x3h);

  return op_cont_pert;
}



fp_t ***** atmosphere::opacity_pert(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in,fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){

  fp_t ***** op_pert = thomson_sc_pert(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);
  //memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  for (int a=0;a<natm;++a) op_pert = add(atml[a]->opacity_pert(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
  return op_pert;
}

fp_t ***** atmosphere::opacity_pert_lte(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in,fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){
  fp_t ***** op_pert = thomson_sc_pert(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);
  //memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  for (int a=0;a<natm;++a)
    if (!atml[a]->check_if_nlte()){
      op_pert = add(atml[a]->opacity_pert(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
    }
      
  return op_pert;
}

fp_t ***** atmosphere::emissivity_pert(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){

  fp_t ***** em_pert = thomson_em_pert(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);
  //memset(em_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  for (int a=0;a<natm;++a)em_pert = add(atml[a]->emissivity_pert(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),em_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);

  return em_pert;
}

fp_t ***** atmosphere::emissivity_pert_lte(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in,fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){
  fp_t ***** em_pert = thomson_em_pert(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);
  //memset(em_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t)); 
  for (int a=0;a<natm;++a)
    if (!atml[a]->check_if_nlte()){
      em_pert = add(atml[a]->emissivity_pert(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),em_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
    }
  return em_pert;
}


// -------------------------------------------------------------------------------------------------------------------------------------------------

fp_t ***** atmosphere::opacity_vector(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){
  
  fp_t ***** op = ft5dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  memset(op[x1l][x2l][x3l][1]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  fp_t *** thompson_scattering = thomson_sc(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        op[x1i][x2i][x3i][1][1] = op[x1i][x2i][x3i][2][2] = op[x1i][x2i][x3i][3][3] = op[x1i][x2i][x3i][4][4] = thompson_scattering[x1i][x2i][x3i];

  del_ft3dim(thompson_scattering,x1l,x1h,x2l,x3h,x3l,x3h);

  for (int a=0;a<natm;++a)
    add(atml[a]->opacity_vector(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);

  return op;
}

fp_t ****atmosphere::emissivity_vector(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){
  fp_t **** em_vector = ft4dim(x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  memset(em_vector[x1l][x2l][x3l]+1,0,(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));

  fp_t ***em_scalar=thomson_em(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h); // electron scattering
  
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i)
        em_vector[x1i][x2i][x3i][1] = em_scalar[x1i][x2i][x3i];

  del_ft3dim(em_scalar,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int a=0;a<natm;++a)
    add(atml[a]->emissivity_vector(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),em_vector,x1l,x1h,x2l,x2h,x3l,x3h,1,4);

  return em_vector;
}

// Then the same but perturbations:

fp_t ******* atmosphere::opacity_vector_pert(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){
  fp_t ******* op_vector_pert = ft7dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);
  memset(op_vector_pert[1][x3l][x1l][x2l][x3l][1]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*16*sizeof(fp_t));

  fp_t ***** op_scalar_pert = thomson_sc_pert(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int p=1;p<=7;++p)
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i){
            op_vector_pert[p][x3k][x1i][x2i][x3i][1][1] = op_vector_pert[p][x3k][x1i][x2i][x3i][2][2] = op_vector_pert[p][x3k][x1i][x2i][x3i][3][3] = 
              op_vector_pert[p][x3k][x1i][x2i][x3i][4][4] = op_scalar_pert[p][x3k][x1i][x2i][x3i];
  }

  del_ft5dim(op_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int a=0;a<natm;++a)
    add(atml[a]->opacity_vector_pert(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),op_vector_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4,1,4);

  return op_vector_pert;
}

fp_t ****** atmosphere::emissivity_vector_pert(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in, fp_t ****B,fp_t theta,fp_t phi,fp_t lambda){
  
  fp_t ****** em_vector_pert = ft6dim(1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);
  memset(em_vector_pert[1][x3l][x1l][x2l][x3l]+1,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*4*sizeof(fp_t));
  
  fp_t ***** em_scalar_pert = thomson_em_pert(Ne_in,lambda,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int p=1;p<=7;++p)
    for (int x3k=x3l;x3k<=x3h;++x3k)
      for (int x1i=x1l;x1i<=x1h;++x1i)
        for (int x2i=x2l;x2i<=x2h;++x2i)
          for (int x3i=x3l;x3i<=x3h;++x3i)
            em_vector_pert[p][x3k][x1i][x2i][x3i][1] = em_scalar_pert[p][x3k][x1i][x2i][x3i];

  del_ft5dim(em_scalar_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);

  for (int a=0;a<natm;++a)
    add(atml[a]->emissivity_vector_pert(T_in,Ne_in,Vlos,Vt_in, B, theta,phi,lambda),em_vector_pert,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h,1,4);


  return em_vector_pert;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------

fp_t ***atmosphere::thomson_sc(fp_t ***Ne_in,fp_t lambda,int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
// Thomson scattering on electrons: classic Thomson crossection times electron density
{
  fp_t ***op=ft3dim(x1l_in,x1h_in,x2l_in,x2h_in,x3l_in,x3h_in);
  fp_t sigma=6.65E-25; // As given in Stellar Atmospheres, ed by Milic
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i)
        op[x1i][x2i][x3i]=Ne_in[x1i][x2i][x3i]*sigma;
  return op;
}

fp_t ***** atmosphere::thomson_sc_pert(fp_t ***Ne_in,fp_t lambda,int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in){

  fp_t sigma=6.65E-25;
  fp_t ***** op_pert = ft5dim(1,7,x3l,x3h,x1l, x1h, x2l, x2h, x3l, x3h);
  memset(op_pert[1][x3l][x1l][x2l]+x3l, 0, 7*(x3h-x3l+1)*(x3h-x3l+1)*(x2h-x2l+1)*(x1h-x1l+1)*sizeof(fp_t));
  for(int param=1;param<=7;++param)
    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i)
          op_pert[param][x3i][x1i][x2i][x3i] = Ne_lte_der[param][x1i][x2i][x3i] * sigma;
  return op_pert;
}


fp_t ***atmosphere::thomson_em(fp_t ***Ne_in,fp_t lambda,int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in)
// Thomson scattering on electrons: classic Thomson crossection times electron density
{
  fp_t ***em=ft3dim(x1l_in,x1h_in,x2l_in,x2h_in,x3l_in,x3h_in);
  fp_t sigma=6.65E-25; // As given in Stellar Atmospheres, added by Milic
  for(int x1i=x1l;x1i<=x1h;++x1i)
    for(int x2i=x2l;x2i<=x2h;++x2i)
      for(int x3i=x3l;x3i<=x3h;++x3i){
        fp_t J = Planck_f(lambda, T[x1i][x2i][x3i]);
        em[x1i][x2i][x3i]=Ne_in[x1i][x2i][x3i]*sigma * J;
      }
  return em;
}

fp_t ***** atmosphere::thomson_em_pert(fp_t ***Ne_in,fp_t lambda,int32_t x1l_in,int32_t x1h_in,int32_t x2l_in,int32_t x2h_in,int32_t x3l_in,int32_t x3h_in){

  fp_t sigma=6.65E-25;
  fp_t ***** em_pert = ft5dim(1,7,x3l,x3h,x1l, x1h, x2l, x2h, x3l, x3h);
  memset(em_pert[1][x3l][x1l][x2l]+x3l, 0, 7*(x3h-x3l+1)*(x3h-x3l+1)*(x2h-x2l+1)*(x1h-x1l+1)*sizeof(fp_t));
  for (int x1i=x1l;x1i<=x1h;++x1i)
    for (int x2i=x2l;x2i<=x2h;++x2i)
      for (int x3i=x3l;x3i<=x3h;++x3i){
        fp_t J = Planck_f(lambda, T[x1i][x2i][x3i]);
        for (int param=1;param<=7;++param)
          em_pert[param][x3i][x1i][x2i][x3i]=Ne_lte_der[param][x1i][x2i][x3i]*sigma * J;
        em_pert[1][x3i][x1i][x2i][x3i] += Ne_in[x1i][x2i][x3i]*sigma * Planck_f_derivative(lambda ,T[x1i][x2i][x3i]);
      }

  return em_pert;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// =================================================================================================================================================================================
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// old:
// ===============================================================================================
fp_t *atmosphere::opacity(fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i){}
fp_t *atmosphere::emissivity(fp_t *lambda,int32_t nlambda,int32_t x1i,int32_t x2i,int32_t x3i){}
fp_t *atmosphere::thomson_op(fp_t ne,fp_t *lambda,int32_t nlambda){}
fp_t *atmosphere::thomson_em(fp_t Te,fp_t ne,fp_t *lambda,int32_t nlambda){}
// ===============================================================================================

// Fetching:
fp_t atmosphere::get_T(int x1i, int x2i, int x3i){
  return T[x1i][x2i][x3i];
}

fp_t * atmosphere::get_magnetic_field(int x1i, int x2i, int x3i){} // You cannot get transformed magnetic field.

//-------------------------------------------------------------------------------------------------------------------------------------------

// Debug function:

int atmosphere::op_em_pert_numerical_scalar(fp_t ***T_in,fp_t ***Ne_in,fp_t ***Vlos,fp_t ***Vt_in,
                            fp_t ****B,fp_t theta,fp_t phi,fp_t lambda, fp_t ***** op_pert, fp_t ***** em_pert){

  memset(op_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));
  memset(em_pert[1][x3l][x1l][x2l]+x3l,0,7*(x3h-x3l+1)*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));


  for (int x3k = x3l;x3k<=x3h;++x3k){

    // Temperature finite-difference perturbations

    T[x1l][x2l][x3k] += delta_T * 0.5;

    nltepops();

    for (int a = 0; a<natm; ++a){
      atml[a]->prof_setup();
      atml[a]->prof_init();
    }

    fp_t *** op = opacity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);
    fp_t *** em = emissivity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);

    for (int a = 0; a<natm; ++a)
      atml[a]->prof_clear();


    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          op_pert[1][x3k][x1i][x2i][x3i] = op[x1i][x2i][x3i];
          em_pert[1][x3k][x1i][x2i][x3i] = em[x1i][x2i][x3i];
        }

    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);

    T[x1l][x2l][x3k] -= delta_T;

    nltepops();

    for (int a = 0; a<natm; ++a){
      atml[a]->prof_setup();
      atml[a]->prof_init();
    }

    op = opacity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);
    em = emissivity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);

    for (int a = 0; a<natm; ++a)
      atml[a]->prof_clear();

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          op_pert[1][x3k][x1i][x2i][x3i] -= op[x1i][x2i][x3i];
          em_pert[1][x3k][x1i][x2i][x3i] -= em[x1i][x2i][x3i];
          op_pert[1][x3k][x1i][x2i][x3i] /= delta_T;
          em_pert[1][x3k][x1i][x2i][x3i] /= delta_T;
        }

    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);

    T[x1l][x2l][x3k] += 0.5 * delta_T;

    // Density finite difference perturbations:

    fp_t d_Nt = 1E-3 * Nt[x1l][x2l][x3k];
    Nt[x1l][x2l][x3k] += 0.5 * d_Nt;

    nltepops();

    for (int a = 0; a<natm; ++a){
      atml[a]->prof_setup();
      atml[a]->prof_init();
    }

    op = opacity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);
    em = emissivity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);

    for (int a = 0; a<natm; ++a)
      atml[a]->prof_clear();

    for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          op_pert[2][x3k][x1i][x2i][x3i] = op[x1i][x2i][x3i];
          em_pert[2][x3k][x1i][x2i][x3i] = em[x1i][x2i][x3i];
        }

    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);


    Nt[x1l][x2l][x3k] -= d_Nt;

    nltepops();

    for (int a = 0; a<natm; ++a){
      atml[a]->prof_setup();
      atml[a]->prof_init();
    }

    op = opacity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);
    em = emissivity(T_in, Ne_in, Vlos, Vt_in, B, theta, phi, lambda);

   for (int a = 0; a<natm; ++a)
      atml[a]->prof_clear();

   for (int x1i=x1l;x1i<=x1h;++x1i)
      for (int x2i=x2l;x2i<=x2h;++x2i)
        for (int x3i=x3l;x3i<=x3h;++x3i){
          op_pert[2][x3k][x1i][x2i][x3i] -= op[x1i][x2i][x3i];
          em_pert[2][x3k][x1i][x2i][x3i] -= em[x1i][x2i][x3i];
          op_pert[2][x3k][x1i][x2i][x3i] /= d_Nt;
          em_pert[2][x3k][x1i][x2i][x3i] /= d_Nt;
        }

    del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
    del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
    Nt[x1l][x2l][x3k] += 0.5 * d_Nt;

    nltepops();

    for (int a = 0; a<natm; ++a){
      atml[a]->prof_setup();
      atml[a]->prof_init();
    }

    
  }
  return 0;  
}

fp_t atmosphere::get_opacity_fudge(fp_t lambda){

  lambda *= 1E8;
  fp_t fudge=1.0;

  for (int l=1;l<=N_of-1;++l){
    if (lambda_of[l] < lambda && lambda_of[l+1] > lambda){
      fudge = (lambda-lambda_of[l+1])/(lambda_of[l]-lambda_of[l+1]) * value_of[l] +
        (lambda-lambda_of[l])/(lambda_of[l+1]-lambda_of[l]) * value_of[l+1];
    }
  }
  double lb_coeff = 2.1177*exp(-(lambda-2087.7)*(lambda-2087.7)/2.421E6) + 0.68738;
  if (lambda > 4500.0)
    lb_coeff=0.0;

  fudge*=(1.0+0.666667*lb_coeff);
  return fudge;
}



