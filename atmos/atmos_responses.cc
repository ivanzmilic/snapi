#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "types.h"
#include "io.h"
#include "mem.h"
#include "atmos.h"
#include "profiles.h"
#include "const.h"


int atmosphere::compute_nlte_population_responses(int lvl_of_approximation){
//
	// This is the function which should compute the response of all the relevant populations to the perturbation of all the possible physical
	// quantities at all the depths in the atmosphere. This is hard and interesting work.

	io.msg(IOL_INFO, "nlte population responses::computing level population responses using analytical approach with [%d]\n", lvl_of_approximation);

	
	// Some continuum wavelengths, not sure if they are needed:
	int32_t nlambda=0;
  fp_t *lambda;           // wavelengths for NLTE calculations, these are @ several pre-defined positions, for the continuum. 
  if (nlambda){
    lambda = new fp_t [nlambda];
    for (int l=0;l<nlambda;++l)
    	lambda[l] = 300E-7 + (2200E-7 - 300E-7) / (nlambda-1) * l;
  }
  if (nlambda == 1)
	lambda[0] = 300E-7;
	//Now we have some setup phase, where we pick wavelengths etc, that is, we try to cleverly sort out the wavelength grid required for our computations. For the moment it looks nice.
  for(int a=0;a<natm;++a) atml[a]->rtsetup(x1l,x1h,x2l,x2h,x3l,x3h); // initialize angular/wavelength redist/integration
  for(int a=0;a<natm;++a) lambda=atml[a]->getlambda(lambda,nlambda,T[x1l][x2l][x3h],Nt[x1l][x2l][x3h],Ne[x1l][x2l][x3h]); // compute wavelength grid for NLTE populations
	// This we can keep, why not.
	for(int a=0;a<natm;++a) atml[a]->rtinit();         // clear transition parameters for each species
	for(int a=0;a<natm;++a) if(atml[a]->check_if_nlte()) atml[a]->responses_init(); // clear response parameters for each species
	
	// First solve ones which are in lte:
	for (int a=0;a<natm;++a)
		if(!atml[a]->check_if_nlte()){
			atml[a]->responses_init();
			atml[a]->compute_lte_population_responses_analytical(T,Ne);
			//atml[a]->compute_lte_population_responses();
	}
	
	// Same ordering like in the case of NLTEpop
	if (nlambda ==0){
		io.msg(IOL_INFO, "atmosphere::nlt population responses : no nlte spieces. We are done here.\n");
		for(int a=0;a<natm;++a) atml[a]->rtclean(0,nlambda,x1l,x1h,x2l,x2h,x3l,x3h); // uninitialize angular/wavelength redist/integration
		return 0; 
	}

	// Here we will have to set-up something like a one loop of a NLTE solution:// radiative quantities
	fp_t ***S=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // Monochromatic intensity in given direction
	fp_t ** response_to_op = ft2dim(x3l, x3h, x3l, x3h); // Full lambda operator-like response of each intensity to each opacity
	fp_t ** response_to_em = ft2dim(x3l, x3h, x3l, x3h); // Same but for opacity

	// setup geometry specific angular quadrature, which is in 1D case gaussian quadrature, typically with 3 directions from -1,0 and viceversa.
	int ntp;
	fp_t *th,*ph;
	fp_t *bin=anglesetup(th,ph,ntp); // bin - integration weights
	                                   // th - theta
	                                   // ph - phi
	                                   // ntp - number of angles, now as this is 1D, unpolarized case, it is actually only number of theta points. In general case it is N_theta x N_phi

	io.msg(IOL_INFO,"atmosphere::nltepops: lambda=%E .. %E nlambda = %d \n",lambda[0],lambda[nlambda-1], nlambda);
	
	// Set up integrating weights for the wavelength
	fp_t * lambda_w = new fp_t [nlambda];
  	if (nlambda >= 2){
      lambda_w[0] = (lambda[1] - lambda[0]) * 0.5;
      lambda_w[nlambda - 1] = (lambda[nlambda-1] - lambda[nlambda-2]) * 0.5;
      for (int l = 1; l<nlambda-1; ++l)
        lambda_w[l] = (lambda[l+1] - lambda[l-1]) * 0.5;
  	}
    else if (nlambda)
      lambda_w[0] = 1.0;  

	// Perform the angular integration: 
	clock_t begin = clock();  
	for(int tp=1;tp<=ntp;++tp){

		fp_t ***Vr=project(Vx,Vy,Vz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h);   // LOS projected velocity
		fp_t ****B=transform(Bx,By,Bz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h); // transform to (B,Inc,Az) representation

		for (int a =0; a<natm; ++a)
			atml[a]->compute_profile_norm(th[tp], ph[tp], lambda, lambda_w, Vr, nlambda); // Normalize the profile as we need the norm

	    for(int l=0;l<nlambda;++l){
  	  
	        for (int a = 0; a<natm; ++a) // Initialize the profile which will be computed only in the opacity and then used from then on
	            atml[a]->prof_init();

	        fp_t ***op=opacity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]); // angle dependent in the scalar case -> NO, at least not at the moment.  
	        fp_t ***em=emissivity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
	        
	        fp_t ***** op_pert_lte = opacity_pert_lte(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
	        fp_t ***** em_pert_lte = emissivity_pert_lte(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
	        
	       	if (tau_grid) normalize_to_referent_opacity(op,em);

	        if(formal_with_responses(rt_grid, S,0,response_to_op,response_to_em,op,em,th[tp],ph[tp], boundary_condition_for_rt)){ // solution and approximate operator
	            io.msg(IOL_ERROR,"atmosphere::nltepops: for angle (%d), wavelength %d\n",tp,l);
	            del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft2dim(response_to_op,x3l, x3h, x3l, x3h);
	            del_ft2dim(response_to_em,x3l, x3h, x3l, x3h); 
	        }
	        if (tau_grid) de_normalize(op,em);
  				for(int a=0;a<natm;++a){ 
	        	atml[a]->add_response_contributions_new(S, response_to_op, response_to_em, op, em, lambda[l], lambda_w[l], th[tp], ph[tp], bin[tp], Vr,
	        	op_pert_lte, em_pert_lte); // give each species access to radiation field, that is, add the radiation field to the mean intensity        
	      }

	      del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h); // cannot be reused due to projections
	      del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
	      del_ft5dim(op_pert_lte,1,7, x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
				del_ft5dim(em_pert_lte,1,7, x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
	    }
	    del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
	    del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
	} 
	clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  //printf("Time spent on adding contributions = %f \n", time_spent); 
	
	delete[] lambda;
	// cleanup angular quadrature arrays
	delete[] (bin+1);
	delete[] (ph+1);
	delete[] (th+1);
	delete[]lambda_w;
	// 
	del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
	del_ft2dim(response_to_op,x3l, x3h, x3l, x3h); 
	del_ft2dim(response_to_em,x3l, x3h, x3l, x3h); 

	// Solve big linear system which gives the responses in the analytical way:
	for (int a = 0; a<natm; ++a)
		if(atml[a]->check_if_nlte())
	    	atml[a]->compute_nlte_population_responses();

	// Print to test:
	//for (int a=0;a<natm;++a)
		//atml[a]->print_population_responses("responses_analytical.txt", x3l, x3h);

	for (int a=0;a<natm;++a){
		atml[a]->rtinit();
	}
	for(int a=0;a<natm;++a) atml[a]->rtclean(ntp,nlambda,x1l,x1h,x2l,x2h,x3l,x3h); // uninitialize angular/wavelength redist/integration	
	io.msg(IOL_INFO, "nlte population responses::responses of the population to the pertubations have been computed\n");

}

void atmosphere::compute_anisotropy_responses(){

	// Here we try to compute anisotropy responses. The whole process is similar to computing the level responses except we now know the level responses.

	io.msg(IOL_INFO, "atmosphere::compute_anisotropy_responses : now computing simple anisotropy responses \n");

	// Here we will have to set-up something like a one loop of a NLTE solution:// radiative quantities
	fp_t ***S=ft3dim(x1l,x1h,x2l,x2h,x3l,x3h); // Monochromatic intensity in given direction
	fp_t ****dS = ft4dim(x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
	fp_t ** response_to_op = ft2dim(x3l, x3h, x3l, x3h); // Full lambda operator-like response of each intensity to each opacity
	fp_t ** response_to_em = ft2dim(x3l, x3h, x3l, x3h); // Same but for opacity

	
	// Some continuum wavelengths, not sure if they are needed:
	int32_t nlambda=5;
    fp_t *lambda=new fp_t [nlambda];           // wavelengths for NLTE calculations, these are @ several pre-defined positions, for the continuum. 
    for (int l=0;l<nlambda;++l)
    lambda[l] = 300E-7 + (2200E-7 - 300E-7) / (nlambda-1) * l;
	lambda[0] = 300E-7;
	// setup geometry specific angular quadrature, which is in 1D case gaussian quadrature, typically with 3 directions from -1,0 and viceversa.
	int ntp;
	fp_t *th,*ph;
	fp_t *bin=anglesetup(th,ph,ntp); // bin - integration weights
	                                   // th - theta
	                                   // ph - phi
	                                   // ntp - number of angles, now as this is 1D, unpolarized case, it is actually only number of theta points. In general case it is N_theta x N_phi

	for(int a=0;a<natm;++a) atml[a]->rtsetup(x1l,x1h,x2l,x2h,x3l,x3h); // initialize angular/wavelength redist/integration
	for(int a=0;a<natm;++a) lambda=atml[a]->getlambda(lambda,nlambda,T[x1l][x2l][x3h],Nt[x1l][x2l][x3h],Ne[x1l][x2l][x3h]); // compute wavelength grid for NLTE populations
	
	// Set up integrating weights for the wavelength
	fp_t * lambda_w = new fp_t [nlambda];
  	if (nlambda >= 2){
      lambda_w[0] = (lambda[1] - lambda[0]) * 0.5;
      lambda_w[nlambda - 1] = (lambda[nlambda-1] - lambda[nlambda-2]) * 0.5;
      for (int l = 1; l<nlambda-1; ++l)
        lambda_w[l] = (lambda[l+1] - lambda[l-1]) * 0.5;
  	}
    else if (nlambda)
      lambda_w[0] = 1.0;  

	// Perform the angular integration:    
	for(int tp=1;tp<=ntp;++tp){

		fp_t ***Vr=project(Vx,Vy,Vz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h);   // LOS projected velocity
		fp_t ****B=transform(Bx,By,Bz,th[tp],ph[tp],x1l,x1h,x2l,x2h,x3l,x3h); // transform to (B,Inc,Az) representation

		for (int a =0; a<natm; ++a)
			atml[a]->compute_profile_norm(th[tp], ph[tp], lambda, lambda_w, Vr, nlambda); // Normalize the profile as we need the norm

	    for(int l=0;l<nlambda;++l){
  	  
	        for (int a = 0; a<natm; ++a) // Initialize the profile which will be computed only in the opacity and then used from then on
	            atml[a]->prof_init();


	        fp_t ***op=opacity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]); // angle dependent in the scalar case -> NO, at least not at the moment.  
	        fp_t ***em=emissivity(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
	        
	        fp_t ***** op_pert = opacity_pert(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);
	        fp_t ***** em_pert = emissivity_pert(T,Ne,Vr,Vt,B,th[tp],ph[tp],lambda[l]);

	        if(formal_with_responses(x3, S,0,response_to_op,response_to_em,op,em,th[tp],ph[tp], boundary_condition_for_rt)){ // solution and approximate operator
	            io.msg(IOL_ERROR,"atmosphere::nltepops: for angle (%d), wavelength %d\n",tp,l);
	            del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
	            del_ft2dim(response_to_op,x3l, x3h, x3l, x3h);
	            del_ft2dim(response_to_em,x3l, x3h, x3l, x3h); 
	        }

	        // But what we need now is also a response of the intensity, explicitly, only for the first three parameters
	        for (int param=1;param<=3;++param){
      			formal_pert_analytical(dS, op, em, op_pert[param], em_pert[param], th[tp], ph[tp], boundary_condition_for_rt);
      			for (int a=0;a<natm;++a)
      			 	atml[a]->add_to_radiation_tensor_perturbation(S, response_to_op, response_to_em, op, em, lambda[l], lambda_w[l], th[tp], ph[tp], bin[tp], Vr,
	        	op_pert, em_pert);
      		}


	       	del_ft3dim(em,x1l,x1h,x2l,x2h,x3l,x3h); // cannot be reused due to projections
	        del_ft3dim(op,x1l,x1h,x2l,x2h,x3l,x3h);
	        del_ft5dim(op_pert,1,7, x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
			del_ft5dim(em_pert,1,7, x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
	    }
	    del_ft4dim(B,1,3,x1l,x1h,x2l,x2h,x3l,x3h);
	    del_ft3dim(Vr,x1l,x1h,x2l,x2h,x3l,x3h);
	}  
	
	delete[] lambda;
	// cleanup angular quadrature arrays
	delete[] (bin+1);
	delete[] (ph+1);
	delete[] (th+1);
	delete[]lambda_w;

	// 
	del_ft3dim(S,x1l,x1h,x2l,x2h,x3l,x3h);
	del_ft4dim(dS,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
	del_ft2dim(response_to_op,x3l, x3h, x3l, x3h); 
	del_ft2dim(response_to_em,x3l, x3h, x3l, x3h); 

	for (int a=0;a<natm;++a)
		if (atml[a]->check_if_nlte())
		atml[a]->print_radiation_field_tensor();

	io.msg(IOL_INFO, "atmosphere::compute_anisotropy_responses : finished computing simple anisotropy responses \n");


}

void atmosphere::compute_nlte_population_responses_numerical(int from, int to){
//
	// This does the same as the function above but using numerical quadratures. This is easy but long computation. However, we have to do it in order to be able to 
	// compare! 

	io.msg(IOL_INFO, "nlte population responses::computing level population responses using finite difference approach, from point %d to point %d\n", from, to);
	
	// You should reset this before beggingin:
	for(int a=0;a<natm;++a) atml[a]->responses_init(); 

	// Compute lte ones:
	for (int a=0;a<natm;++a)
		if (!atml[a]->check_if_nlte())
			atml[a]->compute_lte_population_responses();
	
	io.msg(IOL_INFO, "responses initialized \n");
		
	for (int x3i = from; x3i<=to; ++x3i){

		// Temperature response:
		fp_t T_step = delta_T; //

		T[x1l][x2l][x3i] += T_step * 0.5;
		nltepops();

		for (int a=0;a<natm;++a)
			atml[a]->add_pops_to_response(x3i,1);

		T[x1l][x2l][x3i] -= T_step;
		nltepops();

		for (int a=0;a<natm;++a)
			atml[a]->subtract_pops_from_response(x3i,1);

		for (int a=0;a<natm;++a)
			atml[a]->divide_responses_by_step(x3i,1, T_step);

		T[x1l][x2l][x3i] += T_step * 0.5;

		// Density response:
			
		fp_t Nt_step = 1E-3 * Nt[x1l][x2l][x3i];
		Nt[x1l][x2l][x3i] += Nt_step * 0.5;
		nltepops();

		for (int a=0;a<natm;++a)
			atml[a]->add_pops_to_response(x3i,2);

		Nt[x1l][x2l][x3i] -= Nt_step;
		nltepops();

		for (int a=0;a<natm;++a)
			atml[a]->subtract_pops_from_response(x3i,2);

		for (int a=0;a<natm;++a)
			atml[a]->divide_responses_by_step(x3i,2, Nt_step);

		Nt[x1l][x2l][x3i] += Nt_step * 0.5;

		// Vt response:
		Vt[x1l][x2l][x3i] += delta_vt;
		nltepops();

		for (int a=0;a<natm;++a)
			atml[a]->add_pops_to_response(x3i,3);

		Vt[x1l][x2l][x3i] -= delta_vt;
		nltepops();

		for (int a=0;a<natm;++a)
			atml[a]->subtract_pops_from_response(x3i,3);

		for (int a=0;a<natm;++a)
			atml[a]->divide_responses_by_step(x3i,3, delta_vt);

		Vt[x1l][x2l][x3i] += delta_vt * 0.5;
		
		io.msg(IOL_INFO, "computed population responses with respect to pertrubation at point # %d\n", x3i);
	}
	nltepops();

	for (int a=0;a<natm;++a)
		atml[a]->print_population_responses("responses_numerical.txt", from, to);
	
	io.msg(IOL_INFO, "nlte population responses::responses of the population to the pertubations have been computed using finite differences.\n");
}

// Same version but in taugrid:

void atmosphere::compute_nlte_population_responses_numerical_taugrid(int from, int to){
//
	// This does the same as the function above but using numerical quadratures. This is easy but long computation. However, we have to do it in order to be able to 
	// compare! 

	io.msg(IOL_INFO, "nlte population responses::computing level population responses using finite difference approach, from point %d to point %d\n", from, to);
	
	// You should reset this before beggingin:
	for(int a=0;a<natm;++a) atml[a]->responses_init(); 

	// Compute lte ones:
	for (int a=0;a<natm;++a)
		if (!atml[a]->check_if_nlte())
			atml[a]->compute_lte_population_responses();
	
	io.msg(IOL_INFO, "responses initialized \n");
		
	for (int x3i = from; x3i<=to; ++x3i){

		// Temperature response:
		fp_t T_step = delta_T; //

		T[x1l][x2l][x3i] += T_step * 0.5;
		nltepops_taugrid();

		for (int a=0;a<natm;++a)
			atml[a]->add_pops_to_response(x3i,1);

		T[x1l][x2l][x3i] -= T_step;
		nltepops_taugrid();

		for (int a=0;a<natm;++a)
			atml[a]->subtract_pops_from_response(x3i,1);

		for (int a=0;a<natm;++a)
			atml[a]->divide_responses_by_step(x3i,1, T_step);

		T[x1l][x2l][x3i] += T_step * 0.5;

		// Density response:
		fp_t Nt_step = 1E-3 * Nt[x1l][x2l][x3i];
		Nt[x1l][x2l][x3i] += Nt_step * 0.5;
		nltepops_taugrid();

		for (int a=0;a<natm;++a)
			atml[a]->add_pops_to_response(x3i,2);

		Nt[x1l][x2l][x3i] -= Nt_step;
		nltepops_taugrid();

		for (int a=0;a<natm;++a)
			atml[a]->subtract_pops_from_response(x3i,2);

		for (int a=0;a<natm;++a)
			atml[a]->divide_responses_by_step(x3i,2, Nt_step);

		Nt[x1l][x2l][x3i] += Nt_step * 0.5;

		// Vt response:
		Vt[x1l][x2l][x3i] += delta_vt;
		nltepops_taugrid();

		for (int a=0;a<natm;++a)
			atml[a]->add_pops_to_response(x3i,3);

		Vt[x1l][x2l][x3i] -= delta_vt;
		nltepops_taugrid();

		for (int a=0;a<natm;++a)
			atml[a]->subtract_pops_from_response(x3i,3);

		for (int a=0;a<natm;++a)
			atml[a]->divide_responses_by_step(x3i,3, delta_vt);

		Vt[x1l][x2l][x3i] += delta_vt * 0.5;

		
		io.msg(IOL_INFO, "computed population responses with respect to pertrubation at point # %d\n", x3i);
	}
	nltepops_taugrid();

	for (int a=0;a<natm;++a)
		atml[a]->print_population_responses("responses_numerical.txt", from, to);
	
	io.msg(IOL_INFO, "nlte population responses::responses of the population to the pertubations have been computed using finite differences.\n");
}


fp_t **** atmosphere::compute_intensity_response_numerical(int from, int to, fp_t theta, fp_t phi, fp_t lambda){
	return 0;
}

void atmosphere::ne_lte_derivatives(){

	memset(Ne_lte_der[1][x1l][x2l]+x3l,0,7*(x1h-x1l+1)*(x2h-x2l+1)*(x3h-x3l+1)*sizeof(fp_t));

	for (int x1i=x1l;x1i<=x1h;++x1i)
		for (int x2i=x2l;x2i<=x2h;++x2i)
			for (int x3i=x3l;x3i<=x3h;++x3i){

				// Temperature
				T[x1i][x2i][x3i] += delta_T * 0.5;
				execute_chemeq_for_point(x1i,x2i,x3i);
				Ne_lte_der[1][x1i][x2i][x3i] = Ne[x1i][x2i][x3i];

				for (int a=0;a<natm;++a)
					atml[a]->add_to_ion_responses(1,x3i, 1.0);

				T[x1i][x2i][x3i] -= delta_T;
				execute_chemeq_for_point(x1i,x2i,x3i);
				Ne_lte_der[1][x1i][x2i][x3i] -= Ne[x1i][x2i][x3i];
				Ne_lte_der[1][x1i][x2i][x3i] /= delta_T;

				for (int a=0;a<natm;++a){
					atml[a]->add_to_ion_responses(1,x3i, -1.0);
					atml[a]->divide_ion_responses(1,x3i,delta_T);
				}

				T[x1i][x2i][x3i] += delta_T * 0.5;
				
				fp_t Nt_step = delta_Nt_frac * Nt[x1i][x2i][x3i];

				Nt[x1i][x2i][x3i] += Nt_step * 0.5;
				execute_chemeq_for_point(x1i,x2i,x3i);
				Ne_lte_der[2][x1i][x2i][x3i] = Ne[x1i][x2i][x3i];

				for (int a=0;a<natm;++a)
					atml[a]->add_to_ion_responses(2,x3i, 1.0);

				Nt[x1i][x2i][x3i] -= Nt_step;
				execute_chemeq_for_point(x1i,x2i,x3i);
				Ne_lte_der[2][x1i][x2i][x3i] -= Ne[x1i][x2i][x3i];
				Ne_lte_der[2][x1i][x2i][x3i] /= Nt_step;

				for (int a=0;a<natm;++a){
					atml[a]->add_to_ion_responses(2,x3i, -1.0);
					atml[a]->divide_ion_responses(2,x3i,Nt_step);
				}

				Nt[x1i][x2i][x3i] += Nt_step * 0.5;
				execute_chemeq_for_point(x1i,x2i,x3i);
			}
}

fp_t atmosphere::get_neutral_H_derivative_lte(int param, int x1i, int x2i, int x3i){

	return atml[0]->get_population_response(param,x3i,x1i,x2i,x3i,0);
}

fp_t atmosphere::get_ne_lte_derivative(int param, int x1i, int x2i, int x3i){

	return Ne_lte_der[param][x1i][x2i][x3i];

}

void atmosphere::clear_ne_lte_derivatives(){

	del_ft4dim(Ne_lte_der,1,7,x1l,x1h,x2l,x2h,x3l,x3h);
}

void atmosphere::respsetup(){

	for(int a=0;a<natm;++a) 
		atml[a]->responses_setup();

	Ne_lte_der = ft4dim(1,7,x1l,x1h,x2l,x2h,x3l,x3h);
	memset(Ne_lte_der[1][x1l][x2l]+x3l,0,7*(x1h-x1l+1)*(x2h-x1l+1)*(x3h-x3l+1)*sizeof(fp_t));

}

void atmosphere::respclean(){
  for (int a = 0; a<natm; ++a){
    atml[a]->responses_clear();
  }
  clear_ne_lte_derivatives();
}

void atmosphere::delete_op_referent_derivative(){
	del_ft5dim(op_referent_derivative,1,7,x3l,x3h,x1l,x1h,x2l,x2h,x3l,x3h);
}


// Obsolete:
void atmosphere::compute_nlte_population_responses_taugrid(int lvl_of_approximation){}

