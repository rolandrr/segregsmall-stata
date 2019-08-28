version 14.2
/******************************************************************************/
/* this do-files implements the functions to include covariates
(cf. Appendix B.2 of DHR) */
/******************************************************************************/

/******************************************************************************/
/* structure for results of conditional segregation */
/******************************************************************************/

capture mata: mata drop Struct_results()
mata:
mata set matastrict on
struct Struct_results {
	real matrix res
}
end

capture mata: mata drop Struct_store_bootstrap()
mata:
mata set matastrict on
struct Struct_store_bootstrap {
	struct Struct_results rowvector per_type
	real matrix aggregated
}
end

capture mata: mata drop _View_Results_cond_In_Pb()
mata:
mata set matastrict on
void _View_Results_cond_In_Pb(struct Results_cond scalar results_cond, real scalar type, real scalar K){
	printf("In_Pb_allK_pertype from struct Results_cond	\n")
	printf("For type: \n")
	type
	printf("For K: \n")
	K
	_View_struct_In_Pb_allK(results_cond.In_Pb_allK_pertype[type], K)
}
end

capture mata: mata drop _View_Results_cond_In_mb()
mata:
mata set matastrict on
void _View_Results_cond_In_mb(struct Results_cond scalar results_cond, real scalar type){
	printf("In_mb_Kpind_pertype from struct Results_cond	\n")
	printf("For type: \n")
	type
	_View_struct_In_mb_Kpind(results_cond.In_mb_Kpind_pertype[type])
}
end


/******************************************************************************/
/* Estimates conditional and aggregated */
/******************************************************************************/

capture mata: mata drop Aggregation_across_types()
mata:
mata set matastrict on
real matrix Aggregation_across_types(struct Results_cond scalar results_cond, ///
									real colvector type_probabilities){
	
	real matrix aggregation, temp_matrix
	real scalar nb_rows, nb_cols
	real scalar index, type, nb_types
	
	nb_types = length(type_probabilities)
	nb_rows = rows(results_cond.estimates_ci_per_type[1].estimates_ci)
	nb_cols = cols(results_cond.estimates_ci_per_type[1].estimates_ci)
	aggregation = J(nb_rows, nb_cols, .)
	aggregation[,(1,2)] = results_cond.estimates_ci_per_type[1].estimates_ci[,(1,2)] /* the two first columns are for the characterization of the index : segregation index + weights (unit or individual) */
	
	for(index = 1; index <= nb_rows; index++){
		temp_matrix = J(nb_types, length((3..nb_cols)), .)
		for(type = 1; type <= nb_types; type++){
			temp_matrix[type,] = results_cond.estimates_ci_per_type[type].estimates_ci[index,(3..nb_cols)]
		}
		aggregation[index,(3..nb_cols)] = quadcolsum(type_probabilities :* temp_matrix, 1)
	}
	
	return(aggregation)								
}
end									

capture mata: mata drop Bounds_DTACW_from_struct_db_cond()
mata:
mata set matastrict on
struct Results_cond scalar Bounds_DTACW_from_struct_db_cond(struct Struct_db_cond scalar struct_db_cond, ///
															real scalar b_atkinson, ///
															real scalar I_hyp_independenceKp,
															real scalar I_get_info_for_bootstrapping, ///
															real scalar nb_bootstrap_distance, ///
															real scalar I_trace_estimation, ///
															struct options_optimization scalar options_opt){
	struct Results_cond scalar results_cond
	results_cond = Results_cond()
	real scalar type
	results_cond.estimates_ci_per_type = Struct_eci(struct_db_cond.nb_types)
	results_cond.I_all_Fpk_constrained_per_type = J(struct_db_cond.nb_types, 1, .)
	
	struct raw_ML_results_Kpind colvector raw_ML_Kpind_per_type
	raw_ML_Kpind_per_type = raw_ML_results_Kpind(struct_db_cond.nb_types)
	
	struct results_prob_K_hat colvector prob_K_hat_per_type
	prob_K_hat_per_type = results_prob_K_hat(struct_db_cond.nb_types)
	struct raw_ML_results_allK colvector ML_allK_per_type
	ML_allK_per_type = raw_ML_results_allK(struct_db_cond.nb_types)

	/* Estimation of the bounds (begin) */
	if (I_trace_estimation) {
		displayas("text")
		printf("Estimation - current type analyzed (out of %f):\n", struct_db_cond.nb_types) 
		displayflush()
	}
	for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types (begin) */
		
		if (I_trace_estimation) {
			displayas("text")
			if (I_hyp_independenceKp) {
				printf(" - %f, ", type)
				displayflush()
			}
			else {
				printf(" - %f, current unit size analyzed (out of %f distinct sizes): ", type, struct_db_cond.info_data_per_type[type].nb_K_positive_obs)
				displayflush()
			}	
		}

		if ((I_hyp_independenceKp) & ((struct_db_cond.info_data_per_type[type].nb_K_positive_obs > 1))) { /* hypothesis: independence K and p (begin) */
			/* Kpind and at least two different K ; just for the estimation for consistency of the output for the estimated bounds
			not necessary for bootstrap, it is the same with some numerical approximations only */		
			if (I_trace_estimation) {
				displayas("text")
				printf("K and p assumed independent: units are merged (maximal size = %f)\n", struct_db_cond.info_data_per_type[type].Kbar)
				displayflush()
			}
			raw_ML_Kpind_per_type[type] = Estimation_ML_Kpind(struct_db_cond.db_per_type[type].db, options_opt.CML)
			results_cond.estimates_ci_per_type[type].estimates_ci = Bounds_DTACW_Kpind(b_atkinson, raw_ML_Kpind_per_type[type], options_opt)
			results_cond.I_all_Fpk_constrained_per_type[type] = raw_ML_Kpind_per_type[type].I_constrained_distribution
			/*OLD:results_cond.estimates_ci_per_type[type].estimates_ci = Bounds_DTACW_Kpind_from_db(b_atkinson, struct_db_cond.db_per_type[type].db, I_trace_estimation, options_opt)*/	
		} /* hypothesis: independence K and p (end) */
		
		else { /* no hypothesis about dependence between K and p (begin) */
			prob_K_hat_per_type[type] = Construct_results_prob_K_hat(struct_db_cond.db_per_type[type].db)
			ML_allK_per_type[type] = Estimation_ML_allK(struct_db_cond.db_per_type[type].db, prob_K_hat_per_type[type], 0, I_trace_estimation, options_opt.CML)
			results_cond.estimates_ci_per_type[type].estimates_ci = Bounds_DTACW(b_atkinson, prob_K_hat_per_type[type], ML_allK_per_type[type], options_opt)
			results_cond.I_all_Fpk_constrained_per_type[type] = ML_allK_per_type[type].I_all_Fpk_constrained
			/*OLD:results_cond.estimates_ci_per_type[type].estimates_ci = Bounds_DTACW_from_db(b_atkinson, struct_db_cond.db_per_type[type].db, I_trace_estimation, options_opt)*/
		} /* no hypothesis about dependence between K and p (end) */
		
	} /* loop over types (end) */
	
	results_cond.estimates_ci_aggregated = Aggregation_across_types(results_cond, struct_db_cond.type_probabilities) /* both for Kpind or not, since aggregated index just; both for individual- or unit-level characteristics since already taken into account in the type_probabilities */
	results_cond.I_all_type_all_Fpk_constrained = (sum(results_cond.I_all_Fpk_constrained_per_type) == struct_db_cond.nb_types)
	/* Estimation of the bounds (end) */
	
	if (I_get_info_for_bootstrapping){ /* Get information for bootstrapping (begin) */

		if (I_hyp_independenceKp){ /* hypothesis: independence K and p (begin) */
			results_cond.In_mb_Kpind_pertype = struct_In_mb_Kpind(struct_db_cond.nb_types)
		} /* hypothesis: independence K and p (end) */
		else { /* no hypothesis about dependence between K and p (begin) */
			results_cond.In_Pb_allK_pertype = struct_In_Pb_allK(struct_db_cond.nb_types)
		} /* no hypothesis about dependence between K and p (end) */	
		
		displayas("text")
		printf("Preparation of bootstrap - current type analyzed (out of %f):\n", struct_db_cond.nb_types) 
		displayflush()	
	
		for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types (begin) */
			displayas("text")
			if (!(mod(type,10))){
				if (!(mod(type,50))) {
					printf("%f \n", type)
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
				prob_K_hat_per_type[type] = Construct_results_prob_K_hat(struct_db_cond.db_per_type[type].db)
				results_cond.In_mb_Kpind_pertype[type] = Define_In_mb(prob_K_hat_per_type[type], raw_ML_Kpind_per_type[type], options_opt.CML, nb_bootstrap_distance)			
				/*OLD:results_cond.In_mb_Kpind_pertype[type] = Define_In_mb_from_db(struct_db_cond.db_per_type[type].db, options_opt.CML, nb_bootstrap_distance)*/
			} /* hypothesis: independence K and p (end) */
			
			else { /* no hypothesis about dependence between K and p (begin) */
				results_cond.In_Pb_allK_pertype[type] = Define_In_Pb_allK(prob_K_hat_per_type[type], ML_allK_per_type[type], options_opt.CML, nb_bootstrap_distance)
				/*OLD:results_cond.In_Pb_allK_pertype[type] = Define_In_Pb_allK_from_db(struct_db_cond.db_per_type[type].db, options_opt.CML, nb_bootstrap_distance)*/
			} /* no hypothesis about dependence between K and p (end) */		
			
		} /* loop over types (end) */
		printf("\n")
		
	} /* Get information for bootstrapping (end) */
		
	return(results_cond)
}
end	


/******************************************************************************/
/* Bootstrap */
/******************************************************************************/

/* perform the bootstrap for the conditional analysis with unit-level covariates 
draw in the empirical distribution of (Z, K) i.e. in the empirical distribution
of individual-level type and unit size */
capture mata: mata drop Draw_boot_store_ZK_nb_units()
mata:
mata set matastrict on
real matrix Draw_boot_store_ZK_nb_units(real matrix store_ZK_nb_units, real scalar nb_observations){
	real matrix boot_store_ZK_nb_units
	real colvector draws_number
	real scalar index, zk

	draws_number = floor(nb_observations:*runiform(nb_observations,1):+1)

	boot_store_ZK_nb_units = store_ZK_nb_units
	boot_store_ZK_nb_units[,3] = J(rows(boot_store_ZK_nb_units), 1, .)

	index = 0
	for(zk = 1; zk <= rows(store_ZK_nb_units); zk++){
		if (store_ZK_nb_units[zk,3] > 0){
			boot_store_ZK_nb_units[zk,3] = sum((draws_number :>= (index+1)) :* (draws_number :<= (index+store_ZK_nb_units[zk,3])))
			index = index + store_ZK_nb_units[zk,3]
		}
		else {
			boot_store_ZK_nb_units[zk,3] = 0
		}
	}
	
	return(boot_store_ZK_nb_units)
}
end

/* perform the bootstrap at the level of the unit holding constant
the composition of the unit - hence the use of this matrix store_WK_nb_units */
capture mata: mata drop Draw_boot_store_WK_nb_units()
mata:
mata set matastrict on
real matrix Draw_boot_store_WK_nb_units(real matrix store_WK_nb_units){
	real matrix boot_store_WK_nb_units
	real colvector draws_number
	real scalar index, wk
	real scalar nb_types, nb_observations
	
	nb_types = cols(store_WK_nb_units) - 2
	
	nb_observations = quadsum(store_WK_nb_units[,(nb_types+1)],1)
	draws_number = floor(nb_observations:*runiform(nb_observations,1):+1)
	
	boot_store_WK_nb_units = store_WK_nb_units
	boot_store_WK_nb_units[,nb_types+1] = J(rows(boot_store_WK_nb_units), 1, .)
	
	index = 0
	for(wk = 1; wk <= rows(store_WK_nb_units); wk++){
		if (store_WK_nb_units[wk,nb_types+1] > 0){
			boot_store_WK_nb_units[wk,nb_types+1] = sum((draws_number :>= (index+1)) :* (draws_number :<= (index+store_WK_nb_units[wk,nb_types+1])))
			index = index + store_WK_nb_units[wk,nb_types+1]
		}
		else {
			boot_store_WK_nb_units[wk,nb_types+1] = 0
		}
	}

	return(boot_store_WK_nb_units)
}
end

capture mata: mata drop Complete_db_bootstrap()
mata:
mata set matastrict on
real matrix Complete_db_bootstrap(real matrix draws_K_nb_units, ///
								struct struct_In_Pb_allK scalar	In_Pb_allK){
	real matrix db_b
	db_b = J(rows(draws_K_nb_units), rows(draws_K_nb_units)+3, -99)
	
	real colvector list_K_0_obs, list_K_positive_obs
	real scalar i, K
	
	db_b[,(1,2)] = draws_K_nb_units
	
	/* fill K with no observations */
	list_K_0_obs = selectindex(draws_K_nb_units[,2] :== 0)
	if (length(list_K_0_obs)){
		for(i = 1; i <= length(list_K_0_obs); i++) {
			K = list_K_0_obs[i]
			db_b[K, 3..K+3] = J(1, K+1, 0)
		}
	}
	
	/* fill K with observations */
	list_K_positive_obs = selectindex(draws_K_nb_units[,2])
	for(i = 1; i <= length(list_K_positive_obs); i++) {
		K = list_K_positive_obs[i]
		db_b[K, 3..K+3] = my_rndmultinomial(draws_K_nb_units[K,2], In_Pb_allK.colv_In_Pb[K].Pb)
	}
	
	return(db_b)
}
end								

capture mata: mata drop Complete_db_bootstrap_Kpind()
mata:
mata set matastrict on
real matrix Complete_db_bootstrap_Kpind(real matrix draws_K_nb_units, ///
								struct struct_In_mb_Kpind scalar In_mb){
	real matrix db_b
	db_b = J(rows(draws_K_nb_units), rows(draws_K_nb_units)+3, -99)
	
	real colvector list_K_0_obs, list_K_positive_obs
	real colvector P_hat_b_k
	real scalar i, K
	
	db_b[,(1,2)] = draws_K_nb_units
	
	/* fill K with no observations */
	list_K_0_obs = selectindex(draws_K_nb_units[,2] :== 0)
	if (length(list_K_0_obs)){
		for(i = 1; i <= length(list_K_0_obs); i++) {
			K = list_K_0_obs[i]
			db_b[K, 3..K+3] = J(1, K+1, 0)
		}
	}
	
	/* fill K with observations */
	list_K_positive_obs = selectindex(draws_K_nb_units[,2])
	for(i = 1; i <= length(list_K_positive_obs); i++) {
		K = list_K_positive_obs[i]
		/* Definition of the P_hat_b_k (cf. appendix D.4 DHR) */
		P_hat_b_k = complete_P(Def_Q(K) * In_mb.mb[(1..K)])
		db_b[K, 3..K+3] = my_rndmultinomial(draws_K_nb_units[K,2], P_hat_b_k)
	}
	
	return(db_b)
}
end	

capture mata: mata drop Draw_Struct_db_cond_bootstrap()
mata:
mata set matastrict on
struct Struct_db_cond scalar Draw_Struct_db_cond_bootstrap(struct Struct_db_cond scalar struct_db_cond, ///
														struct struct_In_Pb_allK vector In_Pb_allK_pertype){
	
	struct Struct_db_cond scalar struct_db_cond_boot
	real scalar type, K, Kpertype_bar
	real matrix draws_K_nb_units, temp_Kpertype_nb_units
	struct_db_cond_boot = Struct_db_cond()
	struct_db_cond_boot.I_unit_level_characteristic = struct_db_cond.I_unit_level_characteristic
	struct_db_cond_boot.I_excludingsingletonpertype = struct_db_cond.I_excludingsingletonpertype
	struct_db_cond_boot.nb_types = struct_db_cond.nb_types

	struct_db_cond_boot.db_per_type = Struct_db(struct_db_cond_boot.nb_types)
	struct_db_cond_boot.info_data_per_type = Info_data(struct_db_cond_boot.nb_types)
	
	if (struct_db_cond_boot.I_unit_level_characteristic) { /* unit-level characteristic (begin) */
		
		struct_db_cond_boot.store_ZK_nb_units = Draw_boot_store_ZK_nb_units(struct_db_cond.store_ZK_nb_units, struct_db_cond.info_data_uncond.nb_units_studied)
		struct_db_cond_boot.info_data_uncond.nb_units_studied = quadsum(struct_db_cond_boot.store_ZK_nb_units[,3],1)
		struct_db_cond_boot.info_data_uncond.nb_individuals = quadsum(struct_db_cond_boot.store_ZK_nb_units[,2] :* struct_db_cond_boot.store_ZK_nb_units[,3],1)	
		struct_db_cond_boot.type_frequencies = Cons_type_frequencies_unit(struct_db_cond_boot)
		struct_db_cond_boot.type_probabilities = Cons_type_probabilities_unit(struct_db_cond_boot)
		struct_db_cond_boot.nb_units_studied_per_type = struct_db_cond_boot.type_frequencies
		
		for(type = 1; type <= struct_db_cond_boot.nb_types; type++){ /* loop over types (begin) */
			draws_K_nb_units = select(struct_db_cond_boot.store_ZK_nb_units[,(2,3)], (struct_db_cond_boot.store_ZK_nb_units[,1] :== type))
			struct_db_cond_boot.db_per_type[type].db = Complete_db_bootstrap(draws_K_nb_units, In_Pb_allK_pertype[type])
			struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs = selectindex(struct_db_cond_boot.db_per_type[type].db[,2])
			struct_db_cond_boot.info_data_per_type[type].nb_K_positive_obs = length(struct_db_cond.info_data_per_type[type].list_K_positive_obs)
		} /* loop over types (end) */	
		
	} /* unit-level characteristic (end) */
	
	else { /* individual-level characteristic (begin) */

		struct_db_cond_boot.store_WK_nb_units = Draw_boot_store_WK_nb_units(struct_db_cond.store_WK_nb_units) 
		
		for(type = 1; type <= struct_db_cond_boot.nb_types; type++){ /* loop over types (begin) */	
			Kpertype_bar = max(struct_db_cond_boot.store_WK_nb_units[,type])
			temp_Kpertype_nb_units = select(struct_db_cond_boot.store_WK_nb_units[,(type,(struct_db_cond_boot.nb_types+1))], struct_db_cond_boot.store_WK_nb_units[,type])
			draws_K_nb_units = J(Kpertype_bar, 2, .)	
			for(K = 1; K <= Kpertype_bar; K++){
				if (sum(temp_Kpertype_nb_units[,1] :== K)) {
					draws_K_nb_units[K,] = (K, quadsum(select(temp_Kpertype_nb_units[,2], (temp_Kpertype_nb_units[,1] :== K)),1))				
				}
				else {
					draws_K_nb_units[K,] = (K,0)
				}
			}
			struct_db_cond_boot.db_per_type[type].db = Complete_db_bootstrap(draws_K_nb_units, In_Pb_allK_pertype[type])
			struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs = selectindex(struct_db_cond_boot.db_per_type[type].db[,2])
			struct_db_cond_boot.info_data_per_type[type].nb_K_positive_obs = length(struct_db_cond.info_data_per_type[type].list_K_positive_obs)
		} /* loop over types (end) */
		
		struct_db_cond_boot.nb_units_studied_per_type = Cons_nb_units_per_type_ind(struct_db_cond_boot)
		struct_db_cond_boot.type_frequencies = Cons_type_frequencies_ind(struct_db_cond_boot)
		struct_db_cond_boot.type_probabilities = Cons_type_probabilities_ind(struct_db_cond_boot)
	
	} /* individual-level characteristic (end) */
	
	return(struct_db_cond_boot)
}
end														

/* _i (for avoiding too long name : Kpind case */
capture mata: mata drop Draw_Struct_db_cond_bootstrap_i()
mata:
mata set matastrict on
struct Struct_db_cond scalar Draw_Struct_db_cond_bootstrap_i(struct Struct_db_cond scalar struct_db_cond, ///
																struct struct_In_mb_Kpind vector In_mb_Kpind_pertype){
	
	struct Struct_db_cond scalar struct_db_cond_boot
	real scalar type, K, Kpertype_bar
	real matrix draws_K_nb_units, temp_Kpertype_nb_units
	struct_db_cond_boot = Struct_db_cond()
	struct_db_cond_boot.I_unit_level_characteristic = struct_db_cond.I_unit_level_characteristic
	struct_db_cond_boot.I_excludingsingletonpertype = struct_db_cond.I_excludingsingletonpertype
	struct_db_cond_boot.nb_types = struct_db_cond.nb_types

	struct_db_cond_boot.db_per_type = Struct_db(struct_db_cond_boot.nb_types)
	struct_db_cond_boot.info_data_per_type = Info_data(struct_db_cond_boot.nb_types)
	
	if (struct_db_cond_boot.I_unit_level_characteristic) { /* unit-level characteristic (begin) */
		
		struct_db_cond_boot.store_ZK_nb_units = Draw_boot_store_ZK_nb_units(struct_db_cond.store_ZK_nb_units, struct_db_cond.info_data_uncond.nb_units_studied)
		struct_db_cond_boot.info_data_uncond.nb_units_studied = quadsum(struct_db_cond_boot.store_ZK_nb_units[,3],1)
		struct_db_cond_boot.info_data_uncond.nb_individuals = quadsum(struct_db_cond_boot.store_ZK_nb_units[,2] :* struct_db_cond_boot.store_ZK_nb_units[,3],1)	
		struct_db_cond_boot.type_frequencies = Cons_type_frequencies_unit(struct_db_cond_boot)
		struct_db_cond_boot.type_probabilities = Cons_type_probabilities_unit(struct_db_cond_boot)
		struct_db_cond_boot.nb_units_studied_per_type = struct_db_cond_boot.type_frequencies
		
		for(type = 1; type <= struct_db_cond_boot.nb_types; type++){ /* loop over types (begin) */
			draws_K_nb_units = select(struct_db_cond_boot.store_ZK_nb_units[,(2,3)], (struct_db_cond_boot.store_ZK_nb_units[,1] :== type))
			struct_db_cond_boot.db_per_type[type].db = Complete_db_bootstrap_Kpind(draws_K_nb_units, In_mb_Kpind_pertype[type])
			struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs = selectindex(struct_db_cond_boot.db_per_type[type].db[,2])
			struct_db_cond_boot.info_data_per_type[type].nb_K_positive_obs = length(struct_db_cond.info_data_per_type[type].list_K_positive_obs)
		} /* loop over types (end) */	
		
	} /* unit-level characteristic (end) */
	
	else { /* individual-level characteristic (begin) */

		struct_db_cond_boot.store_WK_nb_units = Draw_boot_store_WK_nb_units(struct_db_cond.store_WK_nb_units) 
				
		for(type = 1; type <= struct_db_cond_boot.nb_types; type++){ /* loop over types (begin) */	
			Kpertype_bar = max(struct_db_cond_boot.store_WK_nb_units[,type])
			temp_Kpertype_nb_units = select(struct_db_cond_boot.store_WK_nb_units[,(type,(struct_db_cond_boot.nb_types+1))], struct_db_cond_boot.store_WK_nb_units[,type])
			draws_K_nb_units = J(Kpertype_bar, 2, .)	
			for(K = 1; K <= Kpertype_bar; K++){
				if (sum(temp_Kpertype_nb_units[,1] :== K)) {
					draws_K_nb_units[K,] = (K, quadsum(select(temp_Kpertype_nb_units[,2], (temp_Kpertype_nb_units[,1] :== K)),1))				
				}
				else {
					draws_K_nb_units[K,] = (K,0)
				}
			}
			struct_db_cond_boot.db_per_type[type].db = Complete_db_bootstrap_Kpind(draws_K_nb_units, In_mb_Kpind_pertype[type])
			struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs = selectindex(struct_db_cond_boot.db_per_type[type].db[,2])
			struct_db_cond_boot.info_data_per_type[type].nb_K_positive_obs = length(struct_db_cond.info_data_per_type[type].list_K_positive_obs)
			
		} /* loop over types (end) */
		
		struct_db_cond_boot.nb_units_studied_per_type = Cons_nb_units_per_type_ind(struct_db_cond_boot)
		struct_db_cond_boot.type_frequencies = Cons_type_frequencies_ind(struct_db_cond_boot)
		struct_db_cond_boot.type_probabilities = Cons_type_probabilities_ind(struct_db_cond_boot)
	
	} /* individual-level characteristic (end) */
	
	return(struct_db_cond_boot)
}
end	
