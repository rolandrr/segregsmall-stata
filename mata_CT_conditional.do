version 14.2
/* functions used for conditional analyses with CT method */

/* same matrix output as results_eci.uncond.estimates_ci in mata_CT_stata_output
but just for one type
function to be used for each type in conditional analyses */
capture mata: mata drop Esti_ct_DTACWG_for_one_type()
mata:
mata set matastrict on
real matrix Esti_ct_DTACWG_for_one_type(real matrix db_one_type, ///
										struct Info_data scalar info_data_one_type, ///
										real scalar b_atkinson, ///
										real scalar nb_ct_repetition, ///
										| real scalar alpha){
	
	real scalar nb_arg_option; nb_arg_option = 5

	real matrix esti_ct_DTACWG, store_indices_random_allocation, db_random_allocation
	real scalar b, prop_minority_hat_random, index_seg, multiple_for_trace, nb_col_esti_ct_DTACWG
	
	nb_col_esti_ct_DTACWG = ((args() == nb_arg_option) ? 14 : 12)
	esti_ct_DTACWG = J(5, nb_col_esti_ct_DTACWG, .)
	esti_ct_DTACWG[,1] = (1, 2, 3+b_atkinson, 4, 5)'

	/* proportion_based indices (begin) */
	esti_ct_DTACWG[1,2] = Prop_Duncan(db_one_type, info_data_one_type.prop_minority_hat, info_data_one_type.nb_individuals)
	esti_ct_DTACWG[2,2] = Prop_Theil(db_one_type, info_data_one_type.prop_minority_hat, info_data_one_type.nb_individuals)
	esti_ct_DTACWG[3,2] = Prop_Atkinson(db_one_type, info_data_one_type.prop_minority_hat, info_data_one_type.nb_individuals, b_atkinson)
	esti_ct_DTACWG[4,2] = Prop_Coworker(db_one_type, info_data_one_type.prop_minority_hat, info_data_one_type.nb_individuals)
	esti_ct_DTACWG[5,2] = Prop_Gini(db_one_type, info_data_one_type.prop_minority_hat, info_data_one_type.nb_individuals)
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
	
		db_random_allocation = Construct_db_random_allocation(db_one_type, info_data_one_type.Kbar, info_data_one_type.prop_minority_hat)
		prop_minority_hat_random = Compute_Pi_from_scratch(db_random_allocation) /* number of individuals is the same by construction, but with the binomial approximation, the number of minority individuals may change however */
		store_indices_random_allocation[b,1] = Prop_Duncan(db_random_allocation, prop_minority_hat_random, info_data_one_type.nb_individuals)
		store_indices_random_allocation[b,2] = Prop_Theil(db_random_allocation, prop_minority_hat_random, info_data_one_type.nb_individuals)
		store_indices_random_allocation[b,3] = Prop_Atkinson(db_random_allocation, prop_minority_hat_random, info_data_one_type.nb_individuals, b_atkinson)
		store_indices_random_allocation[b,4] = Prop_Coworker(db_random_allocation, prop_minority_hat_random, info_data_one_type.nb_individuals)
		store_indices_random_allocation[b,5] = Prop_Gini(db_random_allocation, prop_minority_hat_random, info_data_one_type.nb_individuals)
	}
	/* estimation under random allocation (end) */
	
	/* expected mean and standard deviation under random allocation and CT and CFC corrections  + quantiles (begin) */
	for(index_seg = 1; index_seg <= 5; index_seg++){
		esti_ct_DTACWG[index_seg,3] = mean(store_indices_random_allocation[,index_seg])
		esti_ct_DTACWG[index_seg,5] = sqrt(quadvariance(store_indices_random_allocation[,index_seg]))
		esti_ct_DTACWG[index_seg,4] = Compute_CT_correction(esti_ct_DTACWG[index_seg,2],esti_ct_DTACWG[index_seg,3])
		esti_ct_DTACWG[index_seg,6] = Compute_CFC_standard_score(esti_ct_DTACWG[index_seg,2], esti_ct_DTACWG[index_seg,3], esti_ct_DTACWG[index_seg,5])
		/* quantiles of the segregation index generated under random allocation (begin) */
		esti_ct_DTACWG[index_seg,7] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.01)
		esti_ct_DTACWG[index_seg,8] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.05)
		esti_ct_DTACWG[index_seg,9] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.10)
		esti_ct_DTACWG[index_seg,10] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.90)
		esti_ct_DTACWG[index_seg,11] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.95)
		esti_ct_DTACWG[index_seg,12] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 0.99)
		if (args() == nb_arg_option){
			esti_ct_DTACWG[index_seg,13] = my_empirical_quantile(store_indices_random_allocation[,index_seg], alpha)
			esti_ct_DTACWG[index_seg,14] = my_empirical_quantile(store_indices_random_allocation[,index_seg], 1-alpha)
		}
		/* quantiles of the segregation index generated under random allocation + quantiles (end) */	
		
	}
	/* expected mean and standard deviation under random allocation and CT and CFC corrections (end) */

	return(esti_ct_DTACWG)
}
end

capture mata: mata drop Aggregation_across_types_CT()
mata:
mata set matastrict on
real matrix Aggregation_across_types_CT(struct Results_cond scalar results_cond, ///
										real colvector type_probabilities){
	
	real matrix aggregation, temp_matrix
	real scalar nb_rows, nb_cols
	real scalar index, type, nb_types
	
	nb_types = length(type_probabilities)
	nb_rows = rows(results_cond.estimates_ci_per_type[1].estimates_ci)
	nb_cols = cols(results_cond.estimates_ci_per_type[1].estimates_ci)
	aggregation = J(nb_rows, nb_cols, .)
	aggregation[,1] = results_cond.estimates_ci_per_type[1].estimates_ci[,1] /* the first column is for the characterization of the index : segregation index */
	
	for(index = 1; index <= nb_rows; index++){
		temp_matrix = J(nb_types, length((2..nb_cols)), .)
		for(type = 1; type <= nb_types; type++){
			temp_matrix[type,] = results_cond.estimates_ci_per_type[type].estimates_ci[index,(2..nb_cols)]
		}
		aggregation[index,(2..nb_cols)] = quadcolsum(type_probabilities :* temp_matrix, 1)
	}
	
	return(aggregation)								
}
end		
