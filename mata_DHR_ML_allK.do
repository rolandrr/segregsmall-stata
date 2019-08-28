version 14.2
/******************************************************************************/
/* This do-file creates the different mata structures and functions used to 
perform ML for allK (random K case cf. Appendix B.1. of DHR) */
/******************************************************************************/


/******************************************************************************/
/* structure : results_prob_K_hat */
/* this structure contains the estimation related to the distribution of K 
The choice whether to include the units with K = 1 is done in a higher level
command defining the database db used in the estimation
Here all functions apply to a given db */
/******************************************************************************/

capture mata: mata drop results_prob_K_hat()
mata:
mata set matastrict on
struct results_prob_K_hat {
	real scalar Kbar, expectation_K_hat, nb_obs_allK
	real colvector probabilities
}
end

capture mata: mata drop _View_results_prob_K_hat()
mata:
mata set matastrict on
void _View_results_prob_K_hat(struct results_prob_K_hat scalar prob_K_hat) {
	printf("structure results_prob_K_hat - 4 elements \n")
	printf("Kbar = ")
	prob_K_hat.Kbar
	printf("expectation_K_hat = ")
	prob_K_hat.expectation_K_hat
	printf("nb_obs_allK = ")
	prob_K_hat.nb_obs_allK
	printf("probabilities = \n")
	prob_K_hat.probabilities
}
end

capture mata: mata drop _Get_results_prob_K_hat_colv()
mata:
mata set matastrict on
real colvector _Get_results_prob_K_hat_colv(struct results_prob_K_hat scalar prob_K_hat) {
	real colvector res_prob_K_hat_colv
	res_prob_K_hat_colv = J(prob_K_hat.Kbar + 3, 1, .)
	res_prob_K_hat_colv[1] = prob_K_hat.Kbar
	res_prob_K_hat_colv[2] = prob_K_hat.nb_obs_allK
	res_prob_K_hat_colv[3] = prob_K_hat.expectation_K_hat
	res_prob_K_hat_colv[4..prob_K_hat.Kbar + 3] = prob_K_hat.probabilities
	return(res_prob_K_hat_colv)
}
end

capture mata: mata drop Construct_results_prob_K_hat()
mata:
mata set matastrict on
struct results_prob_K_hat scalar Construct_results_prob_K_hat(real matrix db) {
	struct results_prob_K_hat scalar prob_K_hat
	
	prob_K_hat.Kbar = rows(db)
	prob_K_hat.nb_obs_allK = quadsum(db[1..prob_K_hat.Kbar, 2], 1)
	prob_K_hat.probabilities = db[1..prob_K_hat.Kbar, 2] :/ prob_K_hat.nb_obs_allK
	prob_K_hat.expectation_K_hat = (1..prob_K_hat.Kbar) * prob_K_hat.probabilities
	
	return(prob_K_hat)
}
end


/******************************************************************************/
/* structure : raw_ML_results_allK */
/* this structure concatenates the one-fixed-K structure raw_ML_results
it performs ML_estimation for each K */
/******************************************************************************/

capture mata: mata drop raw_ML_results_allK()
mata:
mata set matastrict on
struct raw_ML_results_allK {
	real scalar Kbar, nb_obs_allK
	struct raw_ML_results colvector colv_ML
	real scalar I_all_Fpk_constrained
	real colvector index_K_with_nb_obs_positive, m1k_hat
}
end

capture mata: mata drop _View_I_constrained_allK()
mata:
mata set matastrict on
void _View_I_constrained_allK(struct raw_ML_results_allK scalar ML_allK) {
	real matrix res_I_constrained
	real scalar k
	res_I_constrained = J(ML_allK.Kbar, 2, .)
	for(k = 1; k <= ML_allK.Kbar; k++) {
		res_I_constrained[k,1] = k
		res_I_constrained[k,2] = ML_allK.colv_ML[k].I_constrained_distribution
	}
	printf("Stack vector of I_constrained_distribution for each K = \n")
	res_I_constrained
	printf("I_all_Fpk_constrained = ")
	ML_allK.I_all_Fpk_constrained
}
end

capture mata: mata drop Get_I_constrained_from_ML_allK()
mata:
mata set matastrict on
real colvector Get_I_constrained_from_ML_allK(struct raw_ML_results_allK scalar ML_allK) {
	real colvector I_constrained_KbyK
	real scalar k
	I_constrained_KbyK = J(ML_allK.Kbar, 1, .)
	for(k = 1; k <= ML_allK.Kbar; k++) {
		I_constrained_KbyK[k] = ML_allK.colv_ML[k].I_constrained_distribution
	}
	return(I_constrained_KbyK)
}
end

capture mata: mata drop _View_K_with_nb_obs_positive()
mata:
mata set matastrict on
void _View_K_with_nb_obs_positive(struct raw_ML_results_allK scalar ML_allK) {
	printf("index_K_with_nb_obs_positive (from structure raw_ML_results_allK) = \n")
	ML_allK.index_K_with_nb_obs_positive
}
end

capture mata: mata drop _View_m1k_hat()
mata:
mata set matastrict on
void _View_m1k_hat(struct raw_ML_results_allK scalar ML_allK) {
	printf("m1k_hat (from structure raw_ML_results_allK) = \n")
	ML_allK.m1k_hat
}
end

capture mata: mata drop _Get_m1k_hat()
mata:
mata set matastrict on
real colvector _Get_m1k_hat(struct raw_ML_results_allK scalar ML_allK) {
	return(ML_allK.m1k_hat)
}
end

capture mata: mata drop _Get_m1k_hat_nb_obs_positive()
mata:
mata set matastrict on
real colvector _Get_m1k_hat_nb_obs_positive(struct raw_ML_results_allK scalar ML_allK) {
	return(ML_allK.m1k_hat[ML_allK.index_K_with_nb_obs_positive])
}
end

capture mata: mata drop _View_raw_ML_results_from_allK()
mata:
mata set matastrict on
void _View_raw_ML_results_from_allK(struct raw_ML_results_allK scalar ML_allK, ///
									real scalar K) {
	_View_raw_ML_results(ML_allK.colv_ML[K])																
}
end	

capture mata: mata drop _Get_xy_hat_from_allK()
mata:
mata set matastrict on
real matrix _Get_xy_hat_from_allK(struct raw_ML_results_allK scalar ML_allK, real scalar K){
	return(ML_allK.colv_ML[K].CML.xy_hat)
}
end

capture mata: mata drop _Get_m_tilde_oneK()
mata:
mata set matastrict on
real colvector _Get_m_tilde_oneK(	struct raw_ML_results_allK scalar ML_allK,
									real scalar K) {
	return(_Get_m_tilde_raw_ML_results(ML_allK.colv_ML[K]))						
}
end

capture mata: mata drop _Get_m_tilde_from_ML_allK()
mata:
mata set matastrict on
real colvector _Get_m_tilde_from_ML_allK(	struct raw_ML_results_allK scalar ML_allK, ///
											real scalar K) {
	return(ML_allK.colv_ML[K].UML.m_tilde)
}
end

capture mata: mata drop Estimation_ML_allK()
mata:
mata set matastrict on
struct raw_ML_results_allK scalar Estimation_ML_allK(	real matrix db, ///
														struct results_prob_K_hat scalar prob_K_hat, ///
														real scalar I_trace_estimation, ///
														real scalar I_trace_estimation_type, ///
														| struct options_optimization_CML optionsoptCML_arg) {
	struct raw_ML_results_allK scalar ML_allK
	
	external optionsoptCML_default; struct options_optimization_CML scalar optionsoptCML
	optionsoptCML = ((args() == 5) ? optionsoptCML_arg : optionsoptCML_default)	
	
	ML_allK.Kbar = prob_K_hat.Kbar
	ML_allK.nb_obs_allK = prob_K_hat.nb_obs_allK
	
	real scalar k, mod_for_space, index_K
	real colvector colv_I_constrained
	
	ML_allK.colv_ML = raw_ML_results(prob_K_hat.Kbar)
	ML_allK.m1k_hat = J(prob_K_hat.Kbar, 1, .)
	colv_I_constrained = J(prob_K_hat.Kbar, 1, .)
	ML_allK.index_K_with_nb_obs_positive = selectindex(db[., 2])

	if (I_trace_estimation) {
		displayas("text")
		printf("Estimation - current unit size analyzed (out of %f distinct sizes):\n", length(ML_allK.index_K_with_nb_obs_positive))
		displayflush()
	}
	
	if (I_trace_estimation) mod_for_space = 40
	if (I_trace_estimation_type) mod_for_space = 30
		
	for(index_K = 1; index_K <= length(ML_allK.index_K_with_nb_obs_positive); index_K++) { /* loop over k (begin) */
		
		if ((I_trace_estimation)|(I_trace_estimation_type)) { /* trace estimation (begin) */
			displayas("text")
			if (!(mod(index_K,10))){
				if ((!(mod(index_K,mod_for_space)))&(I_trace_estimation)){ /* case to put a line back */
					printf("%f \n", index_K)
					displayflush()
				}
				else {
					printf("+")
					displayflush()
				}
			}
			else {
				printf(".")
				displayflush()
			}
		} /* trace estimation (end) */
		
		k = ML_allK.index_K_with_nb_obs_positive[index_K]
		ML_allK.colv_ML[k] = Estimation_ML(db[k, 1..k+3], optionsoptCML)		
		colv_I_constrained[k] = ML_allK.colv_ML[k].I_constrained_distribution
		ML_allK.m1k_hat[k] = ML_allK.colv_ML[k].CML.m1_hat		
	} /* loop over k (end) */
	
	/* OLD
	for(k = 1; k <= prob_K_hat.Kbar; k++) { /* loop over k, with or without observations (begin) */
		
		if ((I_trace_estimation) | (I_trace_estimation_type)) { /* trace estimation */
			displayas("text")
			if (!(mod(k,10))){
				if (!(mod(k,mod_for_space))) { /* case to put a line back */
					if (I_trace_estimation) {
						printf("%f \n", k)
						displayflush()
					}
					else {
						printf("%f \n {space 45}", k)
						displayflush()
					}
				}	
				else { /* case to put a "+" */
					printf("+")
					displayflush()
				}
			}
			else {
				printf(".")
				displayflush()
			}	
		} /* trace estimation (end) */	
				
		ML_allK.colv_ML[k] = Estimation_ML(db[k, 1..k+3], optionsoptCML)		
		colv_I_constrained[k] = ML_allK.colv_ML[k].I_constrained_distribution
		ML_allK.m1k_hat[k] = ML_allK.colv_ML[k].CML.m1_hat	
				
		/*if (prob_K_hat.Kbar > 20){
			if ((!mod(k, 5)) & (k>=20)) {
				printf("end of optimization for K ="); k
			}
		
		}*/
	} /* loop over k, with or without observations (end) */
	*/
		
	if ((I_trace_estimation) | (I_trace_estimation_type)) {
		printf("\n")
		displayflush()
	}
	
	ML_allK.I_all_Fpk_constrained = (sum(colv_I_constrained[ML_allK.index_K_with_nb_obs_positive]) :== length(ML_allK.index_K_with_nb_obs_positive))

	return(ML_allK)
}
end
