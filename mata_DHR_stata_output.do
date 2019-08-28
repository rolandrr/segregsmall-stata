version 14.2
/******************************************************************************/
/* This do-files define the wrap-up functions that are called from the 
ado-file defining the final command */
/******************************************************************************/

/******************************************************************************/
/* Tool fonction to select the bounds (quantities of interest) from the 
matrix of estimates for bootstrap */
/******************************************************************************/
capture mata: mata drop Get_bounds_for_bootstrap_HR()
mata:
mata set matastrict on
real rowvector Get_bounds_for_bootstrap_HR(real matrix bounds_estimates){
	return((bounds_estimates[1,(3,4)],bounds_estimates[2,(3,4)],bounds_estimates[3,(3,4)],bounds_estimates[4,(3,4)],bounds_estimates[5,(3,4)],bounds_estimates[6,(3,4)],bounds_estimates[7,(3,4)],bounds_estimates[8,(3,4)]))
}
end

/******************************************************************************/
/* Function to determinate the option nb_bootstrap_distance
as a function of nb_bootstrap_repetition */
/******************************************************************************/
capture mata: mata drop Define_nb_bootstrap_distance()
mata:
mata set matastrict on
real scalar Define_nb_bootstrap_distance(real scalar nb_bootstrap_repetition){
	real scalar nb_bootstrap_distance
	if (nb_bootstrap_repetition < 25) nb_bootstrap_distance = 10
	else if (nb_bootstrap_repetition < 50) nb_bootstrap_distance = 25
	else if (nb_bootstrap_repetition < 100) nb_bootstrap_distance = 50
	else if (nb_bootstrap_repetition < 300) nb_bootstrap_distance = 100
	else if (nb_bootstrap_repetition < 500) nb_bootstrap_distance = 150
	else if (nb_bootstrap_repetition < 1000) nb_bootstrap_distance = floor(0.5 * nb_bootstrap_repetition)
	else nb_bootstrap_distance = floor(0.75 * nb_bootstrap_repetition)
	return(nb_bootstrap_distance)
}
end

/******************************************************************************/
/* Unconditional analysis */
/******************************************************************************/
capture mata: mata drop St_bounds_ci_np_DTACW_uncond()
mata:
mata set matastrict on
struct Results_eci scalar St_bounds_ci_np_DTACW_uncond(struct Struct_db_uncond scalar struct_db_uncond, ///
													real scalar b_atkinson, ///
													real scalar I_hyp_independenceKp, ///
													real scalar nb_bootstrap_repetition, ///
													real scalar I_noinference, ///
													real scalar I_testbinomial, ///
													| real scalar alpha){
	real scalar nb_arg_option; nb_arg_option = 7
	struct Results_eci results_eci													
	struct options_optimization_cib scalar optionsoptcib
	
	/* technical options for optimization (begin) */
	/* technical options relative to optimization: set default 
	Elements of options_opt_default: nb_max_iter (in moptimize), ptol, vtol, nrtol, nb_candidate_cr 
	nb_bootstrap_distance: number of boostrap iteration used to compute k_n, hence d_n
	and thus to decide whether I_n is equal to 1 or 0, to determine which type of boostrap
	and type of confidence intervals are used */
	struct options_optimization scalar options_opt_default
	options_opt_default = Construct_options_optimization(150, 10^(-8), 10^(-9), 10^(-7), 400)
	real scalar nb_bootstrap_distance
	nb_bootstrap_distance = Define_nb_bootstrap_distance(nb_bootstrap_repetition)
	/* technical options for optimization (end) */
	// nb_bootstrap_distance_default
	
	optionsoptcib = Cons_options_optimization_cib(options_opt_default, Cons_options_ci_bootstrap(nb_bootstrap_distance, nb_bootstrap_repetition))
	real matrix estimates_ci, bounds_DTACW, bounds_DTACW_b, db_b  
	struct results_prob_K_hat scalar prob_K_hat, prob_K_hat_b
	struct raw_ML_results_allK scalar ML_allK, ML_allK_b
	struct struct_In_Pb_allK scalar	In_Pb_allK
	struct raw_ML_results_Kpind scalar raw_ML_Kpind
	struct struct_In_mb_Kpind scalar In_mb_Kpind
	real scalar nb_obs, nb_col, b
	real scalar k, index_line, nb_support_points, nb_maximal_support_points, index_K
	real scalar i_r /* index for row in estimates_ci */
	real rowvector i_c_b /* index for the column in store_bootstrap_bounds */
	real matrix store_infodistributionofp
	real scalar I_fixed_K
	real scalar I_geq1_constrained_distribution, LR_allK_data
	real colvector store_bootstrap_LR_allK
		
	real matrix store_bootstrap_bounds
	/* definition of the columns of store_bootstrap_bounds
	segregation index, unit or individual level, first low, second up
	store_bootstrap_bounds[column,] = (1,1,1,1,2,2,2,2,3+b_atkinson,3+b_atkinson,3+b_atkinson,3+b_atkinson,4,4,4,4)
	store_bootstrap_bounds[column,] = (0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1) */
			
	results_eci = Results_eci()
	
	/* Definition of info_eci (begin) */
	results_eci.info_eci.I_method_np = 1
	results_eci.info_eci.I_method_beta = 0
	results_eci.info_eci.I_method_ct = 0
	results_eci.info_eci.I_conditional = 0
	results_eci.info_eci.I_unit_level_characteristic = .
	results_eci.info_eci.b_atkinson = b_atkinson
	results_eci.info_eci.I_noinference = I_noinference
	results_eci.info_eci.nb_bootstrap_repetition = nb_bootstrap_repetition
	results_eci.info_eci.specified_alpha = (args() == nb_arg_option ? alpha : .)
	results_eci.info_eci.I_deltamethod = .
	results_eci.info_eci.I_hyp_independenceKp = I_hyp_independenceKp
	results_eci.info_eci.I_testbinomial = I_testbinomial
	results_eci.info_eci.nb_ct_repetition = .
	/* Definition of info_eci (end) */
	
	if (I_noinference) {
		displayas("text"); printf("\n*** Estimation (without inference) ***\n"); displayflush()
	}
	else {
		displayas("text"); printf("\n*** Estimation and inference ***\n"); displayflush()
	}
	
	prob_K_hat = Construct_results_prob_K_hat(struct_db_uncond.db) /* information about K - in both K independent of p or not */
	I_fixed_K = (struct_db_uncond.info_data.nb_K_positive_obs == 1) /* for consistency of the results, as sometimes small numerical differences in the two optimization, if only one K, take CML */	
	nb_obs = prob_K_hat.nb_obs_allK /* number of observations (i.e. units) */
	
	if (I_noinference){
		estimates_ci = J(8, 4, .)
	}
	else {
		nb_col = ((args() == nb_arg_option) ? 13 : 11)
		estimates_ci = J(8, nb_col, .) /* for output matrix bounds and ci */
	}
	
	if ((I_hyp_independenceKp) & (!I_fixed_K)){ /* hypothesis: independence K and p (begin) */
		
		/* Estimation of the bounds (begin) */
		displayas("text"); printf("Estimation - K and p assumed independent: units are merged (maximal size = %f)\n", prob_K_hat.Kbar); displayflush()
		raw_ML_Kpind = Estimation_ML_Kpind(struct_db_uncond.db,  optionsoptcib.optimization.CML)
		results_eci.uncond.I_all_Fpk_constrained = raw_ML_Kpind.I_constrained_distribution
		bounds_DTACW = Bounds_DTACW_Kpind(b_atkinson, raw_ML_Kpind, optionsoptcib.optimization)
		estimates_ci[,1..4] = bounds_DTACW
		/* Estimation of the bounds (end) */
		
		if (!I_noinference){ /* bootstrap Kpind case (begin) */
		
			displayas("text"); printf("Preparation of bootstrap -\n"); displayflush()
			In_mb_Kpind = Define_In_mb(prob_K_hat, raw_ML_Kpind, optionsoptcib.optimization.CML, optionsoptcib.cib.nb_bootstrap_distance)
			
			store_bootstrap_bounds = J(nb_bootstrap_repetition, 16, .)
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
				db_b = Draw_db_bootstrap_Kpind(prob_K_hat, In_mb_Kpind.mb)
				bounds_DTACW_b = Bounds_DTACW_Kpind_from_db(b_atkinson, db_b, 0, optionsoptcib.optimization)
				store_bootstrap_bounds[b,] = Get_bounds_for_bootstrap_HR(bounds_DTACW_b)
			} /* loop of iterations of bootstrap (end) */
			printf("\n"); displayflush()
			
			for(i_r = 1; i_r <= rows(estimates_ci); i_r++){
				i_c_b = (2*i_r-1, 2*i_r)
				estimates_ci[i_r,5] = Define_In_ci(estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], nb_obs, raw_ML_Kpind.I_constrained_distribution)
				estimates_ci[i_r,6..7] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], 0.01)
				estimates_ci[i_r,8..9] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], 0.05)
				estimates_ci[i_r,10..11] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], 0.10)
				if (args() == nb_arg_option) estimates_ci[i_r,12..13] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], alpha)
			}
		} /* bootstrap Kpind case (end) */
	
		/* info about the distribution of p Kpind case (begin) */
		if (raw_ML_Kpind.I_constrained_distribution){ /* case constrained (begin) */
			nb_support_points = cols(raw_ML_Kpind.CML.xy_hat)
			store_infodistributionofp = J(3, max((6, nb_support_points)), .)
			store_infodistributionofp[1,(1..5)] = (prob_K_hat.Kbar, struct_db_uncond.info_data.nb_units_studied, 1, 1, nb_support_points)
			store_infodistributionofp[(2,3),(1..nb_support_points)] = raw_ML_Kpind.CML.xy_hat
		} /* case constrained (end) */
		else { /* case unconstrained (begin) */
			store_infodistributionofp = J(3, 6, .)
			store_infodistributionofp[1,(1..4)] = (prob_K_hat.Kbar, struct_db_uncond.info_data.nb_units_studied, 1, 0)
		} /* case unconstrained (end) */
		results_eci.uncond.store_infodistributionofp_np = store_infodistributionofp
		/* info about the distribution of p Kpind case (end) */
	
		if ((I_testbinomial)&(nb_bootstrap_repetition>0)) { /* test binomial assumption Kpind case (begin) */
			results_eci.uncond.test_binomial_results = Test_binomial_acrossK(struct_db_uncond.db, prob_K_hat, nb_bootstrap_repetition, optionsoptcib.optimization.CML, optionsoptcib.cib.nb_bootstrap_distance)	
			printf("\n"); displayflush()
		} /* test binomial assumption Kpind case (end) */
	
	} /* hypothesis: independence K and p (end) */
	
	else { /* no hypothesis about dependence between K and p (begin) */
		/* NB: by default, if inference, computation and storage of the test binomial as it is the same bootstrap procedure
		hence, marginal additional cost to to the binomial test */
		
		/* Estimation of the bounds (begin) */
		ML_allK = Estimation_ML_allK(struct_db_uncond.db, prob_K_hat, 1, 0, optionsoptcib.optimization.CML)
		results_eci.uncond.I_all_Fpk_constrained = ML_allK.I_all_Fpk_constrained
		nb_obs = ML_allK.nb_obs_allK
		bounds_DTACW = Bounds_DTACW(b_atkinson, prob_K_hat, ML_allK, optionsoptcib.optimization)
		estimates_ci[,1..4] = bounds_DTACW
				
		/* Estimation of the bounds (end) */
	
		if (!I_noinference){ /* inference general case case (begin) */
		
			/* Computation of the test statistic from data (begin) */
			results_eci.uncond.test_binomial_results = J(1, 2, .)
			I_geq1_constrained_distribution = (sum(Get_I_constrained_from_ML_allK(ML_allK)[ML_allK.index_K_with_nb_obs_positive]) > 0)
			/* 	= 1 when at least one of the distribution F_p^k is contrained 
				= 0 if all distributions F_p^k are unconstrained, in this case, automatically accept H0 and p-value set to 1 by convention */
			if (I_geq1_constrained_distribution) LR_allK_data = Compute_sum_LR_allK(struct_db_uncond.db, prob_K_hat, ML_allK)
			/* Computation of the test statistic from data (end) */	
			
			displayas("text"); printf("Preparation of bootstrap -\n"); displayflush()
			In_Pb_allK = Define_In_Pb_allK(prob_K_hat, ML_allK, optionsoptcib.optimization.CML, optionsoptcib.cib.nb_bootstrap_distance)
			
			if (I_geq1_constrained_distribution) store_bootstrap_LR_allK = J(nb_bootstrap_repetition, 1, .)
			store_bootstrap_bounds = J(nb_bootstrap_repetition, 16, .)
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
				prob_K_hat_b = Construct_results_prob_K_hat(db_b)
				ML_allK_b = Estimation_ML_allK(db_b, prob_K_hat_b, 0, 0, optionsoptcib.optimization.CML)
				bounds_DTACW_b = Bounds_DTACW(b_atkinson, prob_K_hat_b, ML_allK_b, optionsoptcib.optimization)
				store_bootstrap_bounds[b,] = Get_bounds_for_bootstrap_HR(bounds_DTACW_b)			
				if (I_geq1_constrained_distribution) store_bootstrap_LR_allK[b] = Compute_sum_LR_allK(db_b, prob_K_hat_b, ML_allK_b)
			} /* loop of iterations of bootstrap (end) */
			printf("\n"); displayflush()
							
			for(i_r = 1; i_r <= rows(estimates_ci); i_r++){
				i_c_b = (2*i_r-1, 2*i_r)
				estimates_ci[i_r,5] = Define_In_ci(estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], nb_obs, ML_allK.I_all_Fpk_constrained)
				estimates_ci[i_r,6..7] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], 0.01)
				estimates_ci[i_r,8..9] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], 0.05)
				estimates_ci[i_r,10..11] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], 0.10)
				if (args() == nb_arg_option) estimates_ci[i_r,12..13] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], store_bootstrap_bounds[,i_c_b], alpha)
			}
			
			if (I_geq1_constrained_distribution) {
				results_eci.uncond.test_binomial_results[1] = LR_allK_data
				results_eci.uncond.test_binomial_results[2] = my_empirical_p_value(LR_allK_data, store_bootstrap_LR_allK)
			}
			else {
				results_eci.uncond.test_binomial_results[1] = .
				results_eci.uncond.test_binomial_results[2] = 1 /* all F_p^k are unconstrained, p-value set to 1 by convention */
			}
				
		} /* inference general case case (end) */
	
		/* if weird user demand: no inference, no Kpind but required binomial test */
		if ((I_noinference)&(I_testbinomial)&(nb_bootstrap_repetition>0)) {
			results_eci.uncond.test_binomial_results = Test_binomial_acrossK(struct_db_uncond.db, prob_K_hat, nb_bootstrap_repetition, optionsoptcib.optimization.CML, optionsoptcib.cib.nb_bootstrap_distance)	
			printf("\n"); displayflush()
		}	

		/* info about the distribution of p general case (begin) */
		store_infodistributionofp = J(prob_K_hat.Kbar*3, max((6,floor(prob_K_hat.Kbar/2)+2)), .)
		nb_maximal_support_points = 0		
		for(k = 1; k <= prob_K_hat.Kbar; k++) {	
			index_line = 3*(k-1)+1		
			store_infodistributionofp[index_line,(1,2)] = (k, struct_db_uncond.db[k,2])
		}
		for(index_K = 1; index_K <= struct_db_uncond.info_data.nb_K_positive_obs; index_K++){
			k = struct_db_uncond.info_data.list_K_positive_obs[index_K]
			index_line = 3*(k-1)+1
			store_infodistributionofp[index_line,(3,4)] = (ML_allK.colv_ML[k].init_esti.nb_obs/struct_db_uncond.info_data.nb_units_studied, ML_allK.colv_ML[k].I_constrained_distribution)
			if (ML_allK.colv_ML[k].I_constrained_distribution) {
				nb_support_points = cols(ML_allK.colv_ML[k].CML.xy_hat)
				nb_maximal_support_points = max((nb_maximal_support_points, nb_support_points))
				store_infodistributionofp[index_line,5] = nb_support_points
				store_infodistributionofp[(index_line+1,index_line+2),1..nb_support_points] = ML_allK.colv_ML[k].CML.xy_hat
			}		
		}
		results_eci.uncond.store_infodistributionofp_np = store_infodistributionofp[,(1..max((6,nb_maximal_support_points)))]
		/* info about the distribution of p general case (end) */	
		
	} /* no hypothesis about dependence between K and p (end) */
	
	results_eci.uncond.estimates_ci = estimates_ci

/*OLD: for check of the bootstrap distribution */
/*
	real scalar k
	for (k = 1; k <= rows(db); k++){
		displayas("text")
		printf("For k = %f: ", k)
		displayflush()
		In_Pb_allK.colv_In_Pb[k].In
		printf("\n")
		displayflush()
	}
	results_HR.store_bootstrap_bounds = store_bootstrap_bounds
*/
	return(results_eci)
}
end	


/******************************************************************************/
/* Conditional analysis */
/******************************************************************************/
capture mata: mata drop St_bounds_ci_np_DTACW_cond()
mata:
mata set matastrict on
struct Results_eci scalar St_bounds_ci_np_DTACW_cond(struct Struct_db_cond scalar struct_db_cond, ///
														real scalar b_atkinson, ///
														real scalar I_hyp_independenceKp, ///
														real scalar nb_bootstrap_repetition, ///
														| real scalar alpha){
	real scalar nb_arg_option; nb_arg_option = 5
	struct Results_eci scalar results_eci; results_eci = Results_eci()
	
	struct options_optimization_cib scalar optionsoptcib
	
	/* technical options for optimization (begin) */
	/* technical options relative to optimization: set default 
	Elements of options_opt_default: nb_max_iter (in moptimize), ptol, vtol, nrtol, nb_candidate_cr 
	nb_bootstrap_distance: number of boostrap iteration used to compute k_n, hence d_n
	and thus to decide whether I_n is equal to 1 or 0, to determine which type of boostrap
	and type of confidence intervals are used */
	struct options_optimization scalar options_opt_default
	options_opt_default = Construct_options_optimization(150, 10^(-8), 10^(-9), 10^(-7), 400)
	real scalar nb_bootstrap_distance
	nb_bootstrap_distance = Define_nb_bootstrap_distance(nb_bootstrap_repetition)
	/* technical options for optimization (end) */
	
	optionsoptcib = Cons_options_optimization_cib(options_opt_default, Cons_options_ci_bootstrap(nb_bootstrap_distance, nb_bootstrap_repetition))
	real scalar nb_col, type, b, nb_obs	
	struct Struct_store_bootstrap scalar struct_store_bootstrap
	struct Struct_db_cond scalar struct_db_cond_boot
	struct Results_cond scalar results_cond_boot
	real scalar i_r /* index for row in estimates_ci */
	real rowvector i_c_b /* index for the column in store_bootstrap_bounds */
	real matrix estimates_ci
		
	/* Definition of info_eci (begin) */
	results_eci.info_eci.I_method_np = 1
	results_eci.info_eci.I_method_beta = 0
	results_eci.info_eci.I_method_ct = 0
	results_eci.info_eci.I_conditional = 1
	results_eci.info_eci.I_unit_level_characteristic = struct_db_cond.I_unit_level_characteristic
	results_eci.info_eci.b_atkinson = b_atkinson
	results_eci.info_eci.I_noinference = (nb_bootstrap_repetition == 0)
	results_eci.info_eci.nb_bootstrap_repetition = nb_bootstrap_repetition
	results_eci.info_eci.specified_alpha = (args() == nb_arg_option ? alpha : .)
	results_eci.info_eci.I_deltamethod = .
	results_eci.info_eci.I_hyp_independenceKp = I_hyp_independenceKp
	results_eci.info_eci.I_testbinomial = .
	results_eci.info_eci.nb_ct_repetition = .
	/* Definition of info_eci (end) */
			
	if (nb_bootstrap_repetition == 0){ /* only estimation without inference by bootstrap (begin) */
	
		displayas("text"); printf("\n*** Estimation (without inference) ***\n"); displayflush()
		results_eci.cond = Bounds_DTACW_from_struct_db_cond(struct_db_cond, b_atkinson, I_hyp_independenceKp, 0, optionsoptcib.cib.nb_bootstrap_distance, 1, optionsoptcib.optimization)
		
	} /* only estimation without inference by bootstrap (end) */
	
	
	else { /* estimation and inference by bootstrap (begin) */
	
		displayas("text"); printf("\n*** Estimation and inference ***\n"); displayflush()
	
		/* estimation and preparing bootstrap (begin) */
		results_eci.cond = Bounds_DTACW_from_struct_db_cond(struct_db_cond, b_atkinson, I_hyp_independenceKp, 1, optionsoptcib.cib.nb_bootstrap_distance, 1, optionsoptcib.optimization)
		struct_store_bootstrap = Struct_store_bootstrap()
		struct_store_bootstrap.per_type = Struct_results(struct_db_cond.nb_types)
		for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types (begin) */
			struct_store_bootstrap.per_type[type].res = J(nb_bootstrap_repetition, 16, .)
		} /* loop over types (end) */
		struct_store_bootstrap.aggregated = J(nb_bootstrap_repetition, 16, .)
		/* estimation and preparing bootstrap (end) */
		
		/* bootstrap (begin) */
		displayas("text"); printf("Bootstrap - current bootstrap iteration (out of %f):\n", nb_bootstrap_repetition); displayflush()
		for(b = 1; b <= nb_bootstrap_repetition; b++){ /* loop over bootstrap iterations (begin) */
		
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
			
			if (I_hyp_independenceKp){ /* hypothesis: independence K and p (begin) */
				struct_db_cond_boot = Draw_Struct_db_cond_bootstrap_i(struct_db_cond, results_eci.cond.In_mb_Kpind_pertype)
			} /* hypothesis: independence K and p (end) */
			else{ /* no hypothesis about dependence between K and p (begin) */
				struct_db_cond_boot = Draw_Struct_db_cond_bootstrap(struct_db_cond, results_eci.cond.In_Pb_allK_pertype)
			} /* no hypothesis about dependence between K and p (end) */

			results_cond_boot = Bounds_DTACW_from_struct_db_cond(struct_db_cond_boot, b_atkinson, I_hyp_independenceKp, 0, optionsoptcib.cib.nb_bootstrap_distance, 0, optionsoptcib.optimization)
			for(type = 1; type <= struct_db_cond.nb_types; type++){
				struct_store_bootstrap.per_type[type].res[b,] = Get_bounds_for_bootstrap_HR(results_cond_boot.estimates_ci_per_type[type].estimates_ci)	
			}
			struct_store_bootstrap.aggregated[b,] = Get_bounds_for_bootstrap_HR(results_cond_boot.estimates_ci_aggregated)

		} /* loop over bootstrap iterations (end) */
		printf("\n"); displayflush()
		/* bootstrap (end) */

		nb_col = ((args() == nb_arg_option) ? 13 : 11)
		
		for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types for bounds and ci (begin) */
			estimates_ci = J(8, nb_col, .)
			estimates_ci[,1..4] = results_eci.cond.estimates_ci_per_type[type].estimates_ci
			nb_obs = struct_db_cond.nb_units_studied_per_type[type]	
			for(i_r = 1; i_r <= rows(estimates_ci); i_r++){
				i_c_b = (2*i_r-1, 2*i_r)
				estimates_ci[i_r,5] = Define_In_ci(estimates_ci[i_r,3..4], struct_store_bootstrap.per_type[type].res[,i_c_b], nb_obs, results_eci.cond.I_all_Fpk_constrained_per_type[type])
				estimates_ci[i_r,6..7] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.per_type[type].res[,i_c_b], 0.01)
				estimates_ci[i_r,8..9] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.per_type[type].res[,i_c_b], 0.05)
				estimates_ci[i_r,10..11] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.per_type[type].res[,i_c_b], 0.10)
				if (args() == nb_arg_option) estimates_ci[i_r,12..13] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.per_type[type].res[,i_c_b], alpha)
			}
			results_eci.cond.estimates_ci_per_type[type].estimates_ci = estimates_ci
		} /* loop over types for bounds and ci (end) */
		
		/* bounds and ci for aggregated index (begin) */
		estimates_ci = J(8, nb_col, .)
		estimates_ci[,1..4] = results_eci.cond.estimates_ci_aggregated
		
		/* nb_obs = results_HR_cond.results_cond.nb_units ?
		Ok pour unit-level characteristic 
		Quid pour individual-level ? Plutot nb_obs = sum(results_HR_cond.results_cond.nb_units_per_type)
		Ne joue pas dans IC interior ni IC boundary car auto-normalisation des racine(n) ou n est le nombre d unites 
		Mais va jouer que dans la definition de I_n et donc choix du type d IC (interior or boundary) par contre non ?
		Oui, vu avec Xavier et Roland (point Skype 2018/10/23, on prend la somme des cellules K x W dans le cas individual-level covariate
		i.e. sum(struct_db_cond.nb_units_studied_per_type)
		dans le cas unit-level covariate, cela correspond bien a la somme des unites des differents types et donc au nombre total d'unites */
		nb_obs = sum(struct_db_cond.nb_units_studied_per_type)
		for(i_r = 1; i_r <= rows(estimates_ci); i_r++){
			i_c_b = (2*i_r-1, 2*i_r)
			estimates_ci[i_r,5] = Define_In_ci(estimates_ci[i_r,3..4], struct_store_bootstrap.aggregated[,i_c_b], nb_obs, results_eci.cond.I_all_type_all_Fpk_constrained)
			estimates_ci[i_r,6..7] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.aggregated[,i_c_b], 0.01)
			estimates_ci[i_r,8..9] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.aggregated[,i_c_b], 0.05)
			estimates_ci[i_r,10..11] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.aggregated[,i_c_b], 0.10)
			if (args() == nb_arg_option) estimates_ci[i_r,12..13] = CI_1(estimates_ci[i_r,5], estimates_ci[i_r,3..4], struct_store_bootstrap.aggregated[,i_c_b], alpha)
		}
		results_eci.cond.estimates_ci_aggregated = estimates_ci
		/* bounds and ci for aggregated index (end) */
		
	} /* estimation and inference by bootstrap (end) */
	
	return(results_eci)
}
end
