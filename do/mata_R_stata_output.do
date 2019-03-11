/* Wrap-up functions for estimation and inference */

/******************************************************************************/
/* Unconditional analysis */
/******************************************************************************/
capture mata: mata drop St_esti_ci_beta_DTACWG_uncond()
mata:
mata set matastrict on
struct Results_eci scalar St_esti_ci_beta_DTACWG_uncond(struct Struct_db_uncond scalar struct_db_uncond, ///
													real scalar b_atkinson, ///
													real scalar I_hyp_independenceKp, ///
													real scalar nb_bootstrap_repetition, ///
													real scalar I_noinference, ///
													real scalar I_deltamethod, ///
													real scalar I_testbinomial, ///
													| real scalar alpha){

	real scalar nb_arg_option; nb_arg_option = 8
	
	struct Results_eci results_eci	
	results_eci = Results_eci()
	
	/* technical options for optimization (begin) */
	struct options_optimization_CML optionsoptML_beta_default
	optionsoptML_beta_default = Cons_options_optimization_CML(150, 10^(-8), 10^(-9), 10^(-7))
	struct options_optimization_CML optionsoptCML_default
	optionsoptCML_default = Cons_options_optimization_CML(150, 10^(-8), 10^(-9), 10^(-7))
	real scalar nb_bootstrap_distance_default
	nb_bootstrap_distance_default = 50
	/* technical options for optimization (end) */
	
	real scalar I_fixed_K, b, nb_K_positive_obs_b, index_seg_weight
	real matrix db_b
	struct Info_data scalar info_data_b
		
	real matrix store_bootstrap_estimates
	/* definition of the columns of store_bootstrap_estimate
	couple of estimated index, first unit-level, second individual-level
	1-2nd columns: Duncan; 3-4th: Theil; 5-6th: Atkinson; 7-8th: Coworker; 9-10th: Gini; 
	
	segregation index, unit or individual level, estimates
	store_bootstrap_bounds[column,] = (1,1,1,1,2,2,2,2,3+b_atkinson,3+b_atkinson,3+b_atkinson,3+b_atkinson,4,4,4,4)
	store_bootstrap_bounds[column,] = (0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1) */
	

	/* Definition of info_eci (begin) */
	results_eci.info_eci.I_method_np = 0
	results_eci.info_eci.I_method_beta = 1
	results_eci.info_eci.I_method_ct = 0
	results_eci.info_eci.I_conditional = 0
	results_eci.info_eci.I_unit_level_characteristic = .
	results_eci.info_eci.b_atkinson = b_atkinson
	results_eci.info_eci.I_noinference = I_noinference
	results_eci.info_eci.nb_bootstrap_repetition = nb_bootstrap_repetition
	results_eci.info_eci.specified_alpha = (args() == nb_arg_option ? alpha : .)
	results_eci.info_eci.I_deltamethod = I_deltamethod
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
		
	I_fixed_K = (struct_db_uncond.info_data.nb_K_positive_obs == 1) /* for consistency of the results, as sometimes small numerical differences in the two optimization, if only one K, take the optimization */	
	
	if ((I_hyp_independenceKp) & (!I_fixed_K)){ /* hypothesis: independence K and p (begin) */
		
		displayas("text"); printf("Estimation - K and p assumed independent: units are merged (maximal size = %f)\n", struct_db_uncond.info_data.Kbar); displayflush()
		
		if (I_noinference){ /* Kpind noinference (begin) */
			results_eci.uncond = Estimates_DTACWG_beta_Kpind_runc(b_atkinson, struct_db_uncond.db, struct_db_uncond.info_data, 1, optionsoptML_beta_default)
		} /* Kpind noinference (end) */
	
		else { /* Kpind inference (begin) */
		
			if (I_deltamethod) { /* Kpind inference by deltamethod (begin) */
							
				if (args() == nb_arg_option) results_eci.uncond = Estici_dm_DTACWG_beta_Kpind_runc(b_atkinson, struct_db_uncond.db, struct_db_uncond.info_data, 1, optionsoptML_beta_default, alpha)
				else results_eci.uncond = Estici_dm_DTACWG_beta_Kpind_runc(b_atkinson, struct_db_uncond.db, struct_db_uncond.info_data, 1, optionsoptML_beta_default)
							
			} /* Kpind inference by deltamethod (end) */
			
			else { /* Kpind inference by bootstrap (begin) */
			
				results_eci.uncond = Estimates_DTACWG_beta_Kpind_runc(b_atkinson, struct_db_uncond.db, struct_db_uncond.info_data, 1, optionsoptML_beta_default)
				results_eci.uncond.estimates_ci = (results_eci.uncond.estimates_ci, J(rows(results_eci.uncond.estimates_ci), ((args() == nb_arg_option) ? 8 : 6), .))
			
				store_bootstrap_estimates = J(nb_bootstrap_repetition, 10, .)
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
					db_b = Draw_db_bootstrap_unit_level(struct_db_uncond.db)
					nb_K_positive_obs_b = length(selectindex(db_b[,2]))
					store_bootstrap_estimates[b,] = Estimates_DTACWG_beta_Kpind_eci(b_atkinson, db_b, nb_K_positive_obs_b, 1, optionsoptML_beta_default, 1)
				} /* loop of iterations of bootstrap (end) */
				printf("\n"); displayflush()
				
				for(index_seg_weight = 1; index_seg_weight <= 10; index_seg_weight++){
					results_eci.uncond.estimates_ci[index_seg_weight,(4,5)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], 0.01)
					results_eci.uncond.estimates_ci[index_seg_weight,(6,7)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], 0.05)
					results_eci.uncond.estimates_ci[index_seg_weight,(8,9)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], 0.10)
					if (args() == nb_arg_option) results_eci.uncond.estimates_ci[index_seg_weight,(10,11)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], alpha)
				}		
			} /* Kpind inference by bootstrap (end) */
				
		} /* Kpind inference (end) */
		
	
	} /* hypothesis: independence K and p (end) */
	
	else { /* no hypothesis about dependence between K and p (begin) */
		/* in this case, inference only by bootstrap */
		
		results_eci.uncond = Estimates_DTACWG_beta_resuncond(b_atkinson, struct_db_uncond.db, struct_db_uncond.info_data, 1, optionsoptML_beta_default)
			
		if (!I_noinference){ /* inference general case case (begin) */
			
			results_eci.uncond.estimates_ci = (results_eci.uncond.estimates_ci, J(rows(results_eci.uncond.estimates_ci), ((args() == nb_arg_option) ? 8 : 6), .))
			store_bootstrap_estimates = J(nb_bootstrap_repetition, 10, .)
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
				db_b = Draw_db_bootstrap_unit_level(struct_db_uncond.db)
				info_data_b = Info_data()
				info_data_b.list_K_positive_obs = selectindex(db_b[,2])
				info_data_b.nb_K_positive_obs = length(info_data_b.list_K_positive_obs)
				store_bootstrap_estimates[b,] = Estimates_DTACWG_beta_eci(b_atkinson, db_b, info_data_b, 1, optionsoptML_beta_default, 1, 0)
			} /* loop of iterations of bootstrap (end) */
			printf("\n"); displayflush()
			
			for(index_seg_weight = 1; index_seg_weight <= 10; index_seg_weight++){
				results_eci.uncond.estimates_ci[index_seg_weight,(4,5)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], 0.01)
				results_eci.uncond.estimates_ci[index_seg_weight,(6,7)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], 0.05)
				results_eci.uncond.estimates_ci[index_seg_weight,(8,9)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], 0.10)
				if (args() == nb_arg_option) results_eci.uncond.estimates_ci[index_seg_weight,(10,11)] = CI_by_percentile_bootstrap(store_bootstrap_estimates[,index_seg_weight], alpha)
			}		
		
		} /* inference general case case (end) */
	
	} /* no hypothesis about dependence between K and p (end) */
	
	if ((I_testbinomial)&(nb_bootstrap_repetition>0)) { /* If required: binomial test (begin) */
		results_eci.uncond.test_binomial_results = Test_binomial_acrossK_from_db(struct_db_uncond.db, nb_bootstrap_repetition, optionsoptCML_default, nb_bootstrap_distance_default)	
		printf("\n"); displayflush()
	} /* If required: binomial test (end) */	
		
	return(results_eci)
}
end													


/******************************************************************************/
/* Conditional analysis */
/******************************************************************************/

capture mata: mata drop Get_estimates_for_bootstrap_R()
mata:
mata set matastrict on
real rowvector Get_estimates_for_bootstrap_R(real matrix estimates){
	return(estimates[,3]')
}
end


capture mata: mata drop St_esti_ci_beta_DTACWG_cond()
mata:
mata set matastrict on
struct Results_eci scalar St_esti_ci_beta_DTACWG_cond(	struct Struct_db_cond scalar struct_db_cond, ///
														real scalar b_atkinson, ///
														real scalar I_hyp_independenceKp, ///
														real scalar nb_bootstrap_repetition, ///
														| real scalar alpha){

	real scalar nb_arg_option; nb_arg_option = 5
	struct Results_eci results_eci;	results_eci = Results_eci()
	
	/* technical options for optimization (begin) */
	struct options_optimization_CML optionsoptML_beta_default
	optionsoptML_beta_default = Cons_options_optimization_CML(150, 10^(-8), 10^(-9), 10^(-7))
	/* technical options for optimization (end) */
	
	real scalar nb_col, type, b, index_seg_weight
	struct Struct_store_bootstrap scalar struct_store_bootstrap
	struct Struct_db_cond scalar struct_db_cond_boot
	struct Results_cond scalar results_cond_boot	
	
	/* Definition of info_eci (begin) */
	results_eci.info_eci.I_method_np = 0
	results_eci.info_eci.I_method_beta = 1
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
		results_eci.cond = Esti_DTACWG_beta_from_structcond(struct_db_cond, b_atkinson, I_hyp_independenceKp, optionsoptML_beta_default, 1)
		
	} /* only estimation without inference by bootstrap (end) */

	else { /* estimation and inference by bootstrap (begin) */ 

		displayas("text"); printf("\n*** Estimation and inference ***\n"); displayflush()
		
		/* estimation and preparing bootstrap (begin) */
		results_eci.cond = Esti_DTACWG_beta_from_structcond(struct_db_cond, b_atkinson, I_hyp_independenceKp, optionsoptML_beta_default, 1)
		struct_store_bootstrap = Struct_store_bootstrap()
		struct_store_bootstrap.per_type = Struct_results(struct_db_cond.nb_types)
		for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types (begin) */
			struct_store_bootstrap.per_type[type].res = J(nb_bootstrap_repetition, 10, .)
		} /* loop over types (end) */
		struct_store_bootstrap.aggregated = J(nb_bootstrap_repetition, 10, .)
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
	
			if (struct_db_cond.I_unit_level_characteristic) struct_db_cond_boot = Draw_struct_db_condunit_unit(struct_db_cond) /* unit-level characteristic */
			else struct_db_cond_boot = struct_db_cond_boot = Draw_struct_db_condindi_unit(struct_db_cond) /* individual-level characteristic */
				
			results_cond_boot = Esti_DTACWG_beta_from_structcond(struct_db_cond_boot, b_atkinson, I_hyp_independenceKp, optionsoptML_beta_default, 0)
			
			for(type = 1; type <= struct_db_cond.nb_types; type++){
				struct_store_bootstrap.per_type[type].res[b,] = Get_estimates_for_bootstrap_R(results_cond_boot.estimates_ci_per_type[type].estimates_ci)	
			}
			struct_store_bootstrap.aggregated[b,] = Get_estimates_for_bootstrap_R(results_cond_boot.estimates_ci_aggregated)

		} /* loop over bootstrap iterations (end) */
		printf("\n"); displayflush()
		/* bootstrap (end) */
	
		nb_col = ((args() == nb_arg_option) ? 8 : 6)
	
		for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types for estimates and ci (begin) */
			results_eci.cond.estimates_ci_per_type[type].estimates_ci = (results_eci.cond.estimates_ci_per_type[type].estimates_ci, J(rows(results_eci.cond.estimates_ci_per_type[type].estimates_ci), nb_col, .))
			for(index_seg_weight = 1; index_seg_weight <= 10; index_seg_weight++){
				results_eci.cond.estimates_ci_per_type[type].estimates_ci[index_seg_weight,(4,5)] = CI_by_percentile_bootstrap(struct_store_bootstrap.per_type[type].res[,index_seg_weight], 0.01)
				results_eci.cond.estimates_ci_per_type[type].estimates_ci[index_seg_weight,(6,7)] = CI_by_percentile_bootstrap(struct_store_bootstrap.per_type[type].res[,index_seg_weight], 0.05)
				results_eci.cond.estimates_ci_per_type[type].estimates_ci[index_seg_weight,(8,9)] = CI_by_percentile_bootstrap(struct_store_bootstrap.per_type[type].res[,index_seg_weight], 0.10)
				if (args() == nb_arg_option) results_eci.cond.estimates_ci_per_type[type].estimates_ci[index_seg_weight,(10,11)] = CI_by_percentile_bootstrap(struct_store_bootstrap.per_type[type].res[,index_seg_weight], alpha)
			}	
		} /* loop over types for estimates and ci (end) */
	
		/* estimates and ci for aggregated indices (begin) */
		results_eci.cond.estimates_ci_aggregated = (results_eci.cond.estimates_ci_aggregated, J(rows(results_eci.cond.estimates_ci_aggregated), nb_col, .))
		for(index_seg_weight = 1; index_seg_weight <= 10; index_seg_weight++){
			results_eci.cond.estimates_ci_aggregated[index_seg_weight, (4,5)] = CI_by_percentile_bootstrap(struct_store_bootstrap.aggregated[,index_seg_weight], 0.01)
			results_eci.cond.estimates_ci_aggregated[index_seg_weight, (6,7)] = CI_by_percentile_bootstrap(struct_store_bootstrap.aggregated[,index_seg_weight], 0.05)
			results_eci.cond.estimates_ci_aggregated[index_seg_weight, (8,9)] = CI_by_percentile_bootstrap(struct_store_bootstrap.aggregated[,index_seg_weight], 0.10)
			if (args() == nb_arg_option) results_eci.cond.estimates_ci_aggregated[index_seg_weight, (10,11)] = CI_by_percentile_bootstrap(struct_store_bootstrap.aggregated[,index_seg_weight], alpha)
		}
		/* estimates and ci for aggregated indices (end) */	
	
	} /* estimation and inference by bootstrap (end) */
	
	return(results_eci)
}
end



