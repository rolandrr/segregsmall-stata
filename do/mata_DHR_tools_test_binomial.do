/******************************************************************************/
/* This do-file creates the different mata structures and functions used to 
perform the test of the binomial mixture model - DHR Appendix C.1 + extension
for test across K 
Following discussion with Xavier and Roland (2018/09/07), only the test across
K (the test K by K can be done using the if and in options) */
/******************************************************************************/

/******************************************************************************/
/* Computation of the test statistics and p-value */
/******************************************************************************/

/* Compute_LR_oneK()
compute loglikelihood ratio test statistic LR for a given K - 
cf. formula of LR_n, Annexe C.1 Test of the binomial mixture model 
Warning: it uses the particular structure of data_Nks */
capture mata: mata drop Compute_LR_oneK()
mata:
mata set matastrict on
real scalar Compute_LR_oneK(real rowvector data_Nks, struct raw_ML_results scalar raw_ML){
	if (raw_ML.I_constrained_distribution) { /* constrained case (begin) */	
		return(2 * ((data_Nks[3..length(data_Nks)])[raw_ML.init_esti.list_Nks_Xm1_positive_obs] * ln((complete_P(raw_ML.UML.P_tilde)[raw_ML.init_esti.list_Nks_Xm1_positive_obs]) :/ (complete_P(raw_ML.CML.P_hat)[raw_ML.init_esti.list_Nks_Xm1_positive_obs]))))
	} /* constrained case (begin) */
	else { /* unconstrained case (begin) */
		return(0)
	} /* unconstrained case (end) */
}
end
/* OLD
real colvector P_completed_hat_with_obs, P_completed_tilde_with_obs
real rowvector data_Nks_with_obs
P_completed_hat_with_obs = complete_P(raw_ML.CML.P_hat)[raw_ML.init_esti.list_Nks_Xm1_positive_obs]
P_completed_tilde_with_obs = complete_P(raw_ML.UML.P_tilde)[raw_ML.init_esti.list_Nks_Xm1_positive_obs]
data_Nks_with_obs = (data_Nks[3..length(data_Nks)])[raw_ML.init_esti.list_Nks_Xm1_positive_obs]
*/
		
/* Compute test statistic for the test across K
weighted sum of the LR_one K K by K */
capture mata: mata drop Compute_sum_LR_allK()
mata:
mata set matastrict on
real scalar Compute_sum_LR_allK(real matrix db, struct results_prob_K_hat scalar prob_K_hat, ///
								struct raw_ML_results_allK scalar ML_allK){
	real colvector LRs_KbyK
	real scalar index_K, K
	LRs_KbyK = J(length(ML_allK.index_K_with_nb_obs_positive), 1, .)
	for(index_K = 1; index_K <= length(ML_allK.index_K_with_nb_obs_positive); index_K++){
		K = ML_allK.index_K_with_nb_obs_positive[index_K]
		if (ML_allK.colv_ML[K].I_constrained_distribution){
			LRs_KbyK[index_K] = Compute_LR_oneK(Extract_data_Nks_from_db(db, K), ML_allK.colv_ML[K])
		}
		else {
			LRs_KbyK[index_K] = 0
		}
	}
	return((prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive]') * LRs_KbyK)
}
end

capture mata: mata drop Extract_data_Nks_from_db()
mata:
mata set matastrict on
real rowvector Extract_data_Nks_from_db(real matrix db, real scalar k_chosen){
	return(db[k_chosen,1..(k_chosen+3)])
}
end

/* Wrap-up of Compute_sum_LR_allK from db */
capture mata: mata drop Compute_sum_LR_allK_from_db()
mata:
mata set matastrict on
real scalar Compute_sum_LR_allK_from_db(real matrix db, ///
										real scalar I_trace_estimation, ///
										struct options_optimization_CML scalar optionsoptCML){																
	struct results_prob_K_hat scalar prob_K_hat	
	struct raw_ML_results_allK scalar ML_allK
	prob_K_hat = Construct_results_prob_K_hat(db)
	ML_allK = Estimation_ML_allK(db, prob_K_hat, I_trace_estimation, 0, optionsoptCML)
	return(Compute_sum_LR_allK(db, prob_K_hat, ML_allK))
}
end

/* my_empirical_p_value()
compute p_value from an empirical distribution : colvector of realization
of bootstrap statistics and the observed statistics */
capture mata: mata drop my_empirical_p_value()
mata:
mata set matastrict on
real scalar my_empirical_p_value(	real scalar stat_observed, ///
									real colvector stats_bootstrap) {
	return(quadsum((stats_bootstrap :> stat_observed),1) / length(stats_bootstrap))												
}
end		


/******************************************************************************/
/* Wrap-up function to perform the test starting only from db and prob_K_hat
used for other method (beta, ct) and in the case assumption Kpind */
/******************************************************************************/
capture mata: mata drop Test_binomial_acrossK()
mata:
mata set matastrict on
real rowvector Test_binomial_acrossK(real matrix db, ///
									struct results_prob_K_hat scalar prob_K_hat, ///
									real scalar nb_bootstrap_repetition, ///
									struct options_optimization_CML scalar optionsoptCML, ///
									real scalar nb_bootstrap_distance){
	
	real rowvector test_binomial_results; test_binomial_results = J(1,2,.)
	struct raw_ML_results_allK scalar ML_allK
	real scalar I_geq1_constrained_distribution, b
	real colvector store_bootstrap_LR_allK
	real matrix db_b
	struct struct_In_Pb_allK scalar	In_Pb_allK
	
	displayas("text"); printf("\n*** Additional test of the binomial model ***\n"); displayflush()
	ML_allK = Estimation_ML_allK(db, prob_K_hat, 1, 0, optionsoptCML)
	I_geq1_constrained_distribution = (sum(Get_I_constrained_from_ML_allK(ML_allK)[ML_allK.index_K_with_nb_obs_positive]) > 0)
	
	if (I_geq1_constrained_distribution) { /* case where at least one of the distribution F_p^k is constrained (begin) */
	
		test_binomial_results[1] = Compute_sum_LR_allK(db, prob_K_hat, ML_allK)
	
		displayas("text"); printf("Preparation of bootstrap -\n"); displayflush()
		In_Pb_allK = Define_In_Pb_allK(prob_K_hat, ML_allK, optionsoptCML, nb_bootstrap_distance)
		store_bootstrap_LR_allK = J(nb_bootstrap_repetition, 1, .)
	
		displayas("text"); printf("Bootstrap - current bootstrap iteration (out of %f):\n", nb_bootstrap_repetition); displayflush()
		for(b = 1; b <= nb_bootstrap_repetition; b++){ /* loop of iterations of bootstrap (begin) */
			/* display trace bootstrap (begin) */
			displayas("text")
			if (!(mod(b,10))){
				if (!(mod(b,50))) {
					printf("%f \n", b)
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
			/* display trace bootstrap (end) */
			db_b = Draw_db_bootstrap(prob_K_hat, ML_allK, In_Pb_allK)
			store_bootstrap_LR_allK[b] = Compute_sum_LR_allK_from_db(db_b, 0, optionsoptCML)
		} /* loop of iterations of bootstrap (end) */		
	
		test_binomial_results[2] = my_empirical_p_value(test_binomial_results[1], store_bootstrap_LR_allK)
	
	} /* case where at least one of the distribution F_p^k is constrained (end) */
	
	else { /* case where all the distributions F_p^k are unconstrained (begin) */
		test_binomial_results[1] = .
		test_binomial_results[1] = 1
	} /* case where all the distributions F_p^k are unconstrained (end) */
	
	return(test_binomial_results)
}
end									

capture mata: mata drop Test_binomial_acrossK_from_db()
mata:
mata set matastrict on
real rowvector Test_binomial_acrossK_from_db(real matrix db, ///									
									real scalar nb_bootstrap_repetition, ///
									struct options_optimization_CML scalar optionsoptCML, ///
									real scalar nb_bootstrap_distance){
	struct results_prob_K_hat scalar prob_K_hat
	prob_K_hat = Construct_results_prob_K_hat(db)
	return(Test_binomial_acrossK(db, prob_K_hat, nb_bootstrap_repetition, optionsoptCML, nb_bootstrap_distance)						)									
}
end
