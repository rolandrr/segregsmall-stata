

capture mata: mata drop Esti_DTACWG_beta_from_structcond()
mata:
mata set matastrict on
struct Results_cond scalar Esti_DTACWG_beta_from_structcond(struct Struct_db_cond scalar struct_db_cond, ///
															real scalar b_atkinson, ///
															real scalar I_hyp_independenceKp,
															struct options_optimization options_opt, ///
															real scalar I_trace_estimation){
	struct Results_cond scalar results_cond
	results_cond = Results_cond()
	real scalar type
	results_cond.estimates_ci_per_type = Struct_eci(struct_db_cond.nb_types)


	if (I_trace_estimation) {
		displayas("text")
		printf("Estimation - current type analyzed (out of %f distinct types):\n", struct_db_cond.nb_types) 
		displayflush()
	}
	
	/* Estimation of the index type by type (begin) */
	for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types (begin) */
		
		if (I_trace_estimation) {
			displayas("text")
			if (I_hyp_independenceKp) {
				printf(" - %f, K and p assumed independent: units of different sizes are merged (maximal size = %f)", type, struct_db_cond.info_data_per_type[type].Kbar)
			}
			else {
				printf(" - %f, current unit size analyzed (out of %f distinct sizes): ", type, struct_db_cond.info_data_per_type[type].nb_K_positive_obs)
			}	
			displayflush()
		}
		
		if ((I_hyp_independenceKp) & ((struct_db_cond.info_data_per_type[type].nb_K_positive_obs > 1))) { /* hypothesis: independence K and p (begin) */
			/* Kpind and at least two different K ; just for the estimation for consistency of the output for the estimated bounds
			not necessary for bootstrap, it is the same with some numerical approximations only */
			results_cond.estimates_ci_per_type[type].estimates_ci = Estimates_DTACWG_beta_Kpind_eci(b_atkinson, struct_db_cond.db_per_type[type].db, struct_db_cond.info_data_per_type[type].nb_K_positive_obs, 1, options_opt, 0)
		} /* hypothesis: independence K and p (end) */
		
		else { /* no hypothesis about dependence between K and p (begin) */			
			results_cond.estimates_ci_per_type[type].estimates_ci = Estimates_DTACWG_beta_eci(b_atkinson, struct_db_cond.db_per_type[type].db, struct_db_cond.info_data_per_type[type], 1, options_opt, 0, I_trace_estimation)
		} /* no hypothesis about dependence between K and p (end) */
	
		if (I_trace_estimation) {
			displayflush()
			printf("\n")
		}
		
	} /* loop over types (end) */
	/* Estimation of the index type by type (end) */
	
	/* Estimation of the aggregated index across type (begin) */
	results_cond.estimates_ci_aggregated = Aggregation_across_types(results_cond, struct_db_cond.type_probabilities) /* both for Kpind or not, since aggregated index just; both for individual- or unit-level characteristics since already taken into account in the type_probabilities */
	/* Estimation of the aggregated index across type (end) */
	
	return(results_cond)
}
end
