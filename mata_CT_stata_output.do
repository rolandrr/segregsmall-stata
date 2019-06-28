/* Wrap-up functions for estimation and inference */

/******************************************************************************/
/* Unconditional analysis */
/******************************************************************************/

capture mata: mata drop St_esti_ct_DTACWG_uncond()
mata:
mata set matastrict on
struct Results_eci scalar St_esti_ct_DTACWG_uncond(struct Struct_db_uncond scalar struct_db_uncond, ///
													real scalar b_atkinson, ///
													real scalar nb_ct_repetition, ///
													real scalar nb_bootstrap_repetition, ///
													real scalar I_testbinomial, ///
													| real scalar alpha){
	real scalar nb_arg_option; nb_arg_option = 6
	
	struct Results_eci results_eci	
	results_eci = Results_eci()
	real scalar b, prop_minority_hat_random, index_seg, multiple_for_trace, nb_col_estimates_ci
	real matrix db_random_allocation, store_indices_random_allocation
	
	/* technical options for optimization (begin) */
	struct options_optimization_CML optionsoptCML_default
	optionsoptCML_default = Cons_options_optimization_CML(150, 10^(-8), 10^(-9), 10^(-7))
	real scalar nb_bootstrap_distance_default
	nb_bootstrap_distance_default = 50
	/* technical options for optimization (end) */
			
	/* Definition of info_eci (begin) */
	results_eci.info_eci.I_method_np = 0
	results_eci.info_eci.I_method_beta = 0
	results_eci.info_eci.I_method_ct = 1
	results_eci.info_eci.I_conditional = 0
	results_eci.info_eci.I_unit_level_characteristic = .
	results_eci.info_eci.b_atkinson = b_atkinson
	results_eci.info_eci.I_noinference = .
	results_eci.info_eci.nb_bootstrap_repetition = nb_bootstrap_repetition
	results_eci.info_eci.specified_alpha = (args() == nb_arg_option ? alpha : .)
	results_eci.info_eci.I_deltamethod = .
	results_eci.info_eci.I_hyp_independenceKp = .
	results_eci.info_eci.I_testbinomial = I_testbinomial
	results_eci.info_eci.nb_ct_repetition = nb_ct_repetition
	/* Definition of info_eci (end) */
	
	displayas("text"); printf("\n*** Estimation and correction ***\n"); displayflush()
	
	nb_col_estimates_ci = ((args() == nb_arg_option) ? 14 : 12)
	results_eci.uncond.estimates_ci = J(5, nb_col_estimates_ci, .)
	results_eci.uncond.estimates_ci[,1] = (1, 2, 3+b_atkinson, 4, 5)'
	/* first column specifies the different rows i.e. the different segregation indices, 1=Duncan, 2=Theil, 3.b=Atkinson(b), 4=Coworker, 5=Gini
	second column:	"naive" proportion-based estimate
	third column:	mean of the index under random allocation
	fourth column:	Carrington-Troske corrected index 
	fifth column:	standard deviation of the index under random allocation
	sixth column:	Cortesen Falk and Cohen standardized score 
	seventh to end columns: quantiles of the index generated under random allocation 
	with quantiles 0.01, 0.05, 0.10, 0.90, 0.95, 0.99, alpha, 1-alpha (if a specified alpha)*/
	
	/* proportion_based indices (begin) */
	results_eci.uncond.estimates_ci[1,2] = Prop_Duncan_from_struct_db_unc(struct_db_uncond)
	results_eci.uncond.estimates_ci[2,2] = Prop_Theil_from_struct_db_unc(struct_db_uncond)
	results_eci.uncond.estimates_ci[3,2] = Prop_Atkinson_from_struct_db_unc(struct_db_uncond, b_atkinson)
	results_eci.uncond.estimates_ci[4,2] = Prop_Coworker_from_struct_db_unc(struct_db_uncond)
	results_eci.uncond.estimates_ci[5,2] = Prop_Gini_from_struct_db_unc(struct_db_uncond)
	/* proportion_based indices (end) */
	
	/* estimation under random allocation (begin) */
	if (nb_ct_repetition > 300) multiple_for_trace = 50
	else multiple_for_trace = 10
	displayas("text"); printf("CT correction - current random allocation iteration x%f (out of %f):\n", multiple_for_trace, nb_ct_repetition); displayflush()
	store_indices_random_allocation = J(nb_ct_repetition, 5, .)
	for(b = 1; b <= nb_ct_repetition; b++) {
	
		/* trace ct correction (begin) */
		displayas("text")
		if (!(mod(b,multiple_for_trace))){
			if (!(mod(b,10*multiple_for_trace))){
				printf("+")
				displayflush()
			}
			else {
				printf(".")
				displayflush()
			}
		}
		/* trace ct correction (end) */
	
		db_random_allocation = Construct_db_random_allocation(struct_db_uncond.db, struct_db_uncond.info_data.Kbar, struct_db_uncond.info_data.prop_minority_hat)
		prop_minority_hat_random = Compute_Pi_from_scratch(db_random_allocation) /* number of individuals is the same by construction, but with the binomial approximation, the number of minority individuals may change however */
		store_indices_random_allocation[b,1] = Prop_Duncan(db_random_allocation, prop_minority_hat_random, struct_db_uncond.info_data.nb_individuals)
		store_indices_random_allocation[b,2] = Prop_Theil(db_random_allocation, prop_minority_hat_random, struct_db_uncond.info_data.nb_individuals)
		store_indices_random_allocation[b,3] = Prop_Atkinson(db_random_allocation, prop_minority_hat_random, struct_db_uncond.info_data.nb_individuals, b_atkinson)
		store_indices_random_allocation[b,4] = Prop_Coworker(db_random_allocation, prop_minority_hat_random, struct_db_uncond.info_data.nb_individuals)
		store_indices_random_allocation[b,5] = Prop_Gini(db_random_allocation, prop_minority_hat_random, struct_db_uncond.info_data.nb_individuals)
	}
	printf("\n"); displayflush()
	/* estimation under random allocation (end) */
	
	/* expected mean and standard deviation under random allocation and CT and CFC corrections + quantile (begin) */
	for(index_seg = 1; index_seg <= 5; index_seg++){
		results_eci.uncond.estimates_ci[index_seg,3] = mean(store_indices_random_allocation[,index_seg])
		results_eci.uncond.estimates_ci[index_seg,5] = sqrt(quadvariance(store_indices_random_allocation[,index_seg]))
		results_eci.uncond.estimates_ci[index_seg,4] = Compute_CT_correction(results_eci.uncond.estimates_ci[index_seg,2], results_eci.uncond.estimates_ci[index_seg,3])
		results_eci.uncond.estimates_ci[index_seg,6] = Compute_CFC_standard_score(results_eci.uncond.estimates_ci[index_seg,2], results_eci.uncond.estimates_ci[index_seg,3], results_eci.uncond.estimates_ci[index_seg,5])
		/* quantiles of the segregation index generated under random allocation (begin) */
		results_eci.uncond.estimates_ci[index_seg,7] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.01)
		results_eci.uncond.estimates_ci[index_seg,8] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.05)
		results_eci.uncond.estimates_ci[index_seg,9] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.10)
		results_eci.uncond.estimates_ci[index_seg,10] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.90)
		results_eci.uncond.estimates_ci[index_seg,11] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.95)
		results_eci.uncond.estimates_ci[index_seg,12] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.99)
		if (args() == nb_arg_option){
			results_eci.uncond.estimates_ci[index_seg,13] = my_empirical_quantile(store_indices_random_allocation[,index_seg], alpha)
			results_eci.uncond.estimates_ci[index_seg,14] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 1-alpha)
		}
		/* quantiles of the segregation index generated under random allocation (end) */	
	}
	/* expected mean and standard deviation under random allocation and CT and CFC corrections + quantile (end) */
	
	if ((I_testbinomial)&(nb_bootstrap_repetition>0)) { /* If required: binomial test (begin) */
			results_eci.uncond.test_binomial_results = Test_binomial_acrossK_from_db(struct_db_uncond.db, nb_bootstrap_repetition, optionsoptCML_default, nb_bootstrap_distance_default)	
	} /* If required: binomial test (end) */	
	
	return(results_eci)
}
end													

/* for test only */
capture mata: mata drop _Get_ct_indices()
mata:
mata set matastrict on
real colvector _Get_ct_indices(struct Results_eci results_eci){
	return(results_eci.uncond.estimates_ci[,4])
}
end


/******************************************************************************/
/* Conditional analysis */
/******************************************************************************/
capture mata: mata drop St_esti_ct_DTACWG_cond()
mata:
mata set matastrict on
struct Results_eci scalar St_esti_ct_DTACWG_cond(	struct Struct_db_cond scalar struct_db_cond, ///
													real scalar b_atkinson, ///
													real scalar nb_ct_repetition, ///
													| real scalar alpha){
	
	real scalar nb_arg_option; nb_arg_option = 4
	
	struct Results_eci results_eci;	results_eci = Results_eci()
	real scalar type
	
	/* Definition of info_eci (begin) */
	results_eci.info_eci.I_method_np = 0
	results_eci.info_eci.I_method_beta = 0
	results_eci.info_eci.I_method_ct = 1
	results_eci.info_eci.I_conditional = 1
	results_eci.info_eci.I_unit_level_characteristic = struct_db_cond.I_unit_level_characteristic
	results_eci.info_eci.b_atkinson = b_atkinson
	results_eci.info_eci.I_noinference = .
	results_eci.info_eci.nb_bootstrap_repetition = .
	results_eci.info_eci.specified_alpha = (args() == nb_arg_option ? alpha : .)
	results_eci.info_eci.I_deltamethod = .
	results_eci.info_eci.I_hyp_independenceKp = .
	results_eci.info_eci.I_testbinomial = .
	results_eci.info_eci.nb_ct_repetition = nb_ct_repetition
	/* Definition of info_eci (end) */

	results_eci.cond.estimates_ci_per_type = Struct_eci(struct_db_cond.nb_types)
	
	displayas("text"); printf("\n*** Estimation and correction ***\n"); displayflush()
	
	displayas("text")
	printf("Estimation and correction - current type analyzed (out of %f distinct types):", struct_db_cond.nb_types) 
	displayflush()
	
	/* Estimation and CT correction type by type (begin) */
	for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types (begin) */
		printf("\n - %f, ", type)
		if (args() == nb_arg_option) 	results_eci.cond.estimates_ci_per_type[type].estimates_ci = Esti_ct_DTACWG_for_one_type(struct_db_cond.db_per_type[type].db, struct_db_cond.info_data_per_type[type], b_atkinson, nb_ct_repetition, alpha)
		else 							results_eci.cond.estimates_ci_per_type[type].estimates_ci = Esti_ct_DTACWG_for_one_type(struct_db_cond.db_per_type[type].db, struct_db_cond.info_data_per_type[type], b_atkinson, nb_ct_repetition)
	} /* loop over types (end) */
	printf("\n"); displayflush()
	/* Estimation and CT correction type by type (end) */
	
	/* Estimation of the aggregated index across type (begin) */
	results_eci.cond.estimates_ci_aggregated = Aggregation_across_types_CT(results_eci.cond, struct_db_cond.type_probabilities) /* both for Kpind or not, since aggregated index just; both for individual- or unit-level characteristics since already taken into account in the type_probabilities */
	/* Estimation of the aggregated index across type (end) */
	
	return(results_eci)													
}
end
