/******************************************************************************/
/* This do-file creates structures and functions used later to compute
confidence interval (ci) for Duncan, Theil, Coworker, Atkinson for allK i.e. 
with K random and assumption K and p independent */
/******************************************************************************/

capture mata: mata drop Get_m1hat_Kpind()
mata:
mata set matastrict on
real scalar Get_m1hat_Kpind(struct raw_ML_results_Kpind scalar raw_ML){
	if (raw_ML.I_constrained_distribution){
		return(raw_ML.CML.m1_hat)
	}
	else {
		return(raw_ML.UML.m_tilde[1])
	}
}
end

capture mata: mata drop Theil_oneK_Kpind()
mata:
mata set matastrict on 
real rowvector Theil_oneK_Kpind(struct raw_ML_results_Kpind scalar raw_ML){
	real rowvector int_h_Fpk 
	real scalar m1_hat
	m1_hat = Get_m1hat_Kpind(raw_ML)
	int_h_Fpk = Def_int_h_Fpk_Theil_Kpind(raw_ML, m1_hat)
	return(compute_nu_Theil(int_h_Fpk[1], m1_hat), compute_nu_Theil(int_h_Fpk[2], m1_hat))
}
end

capture mata: mata drop struct_In_mb_Kpind()
mata:
mata set matastrict on 
struct struct_In_mb_Kpind{
	real scalar In
	real colvector mb
}
end

capture mata: mata drop _View_struct_In_mb_Kpind()
mata:
mata set matastrict on
void _View_struct_In_mb_Kpind(struct struct_In_mb_Kpind scalar In_mb){
	printf("Structure struct_In_mb_Kpind \n")
	printf("In =\n")
	In_mb.In
	printf("mb =\n")
	In_mb.mb
}
end

capture mata: mata drop Draw_db_bootstrap_Kpind()
mata:
mata set matastrict on
real matrix Draw_db_bootstrap_Kpind(struct results_prob_K_hat scalar prob_K_hat, ///
									real colvector mb){

	real matrix db_b
	db_b = J(prob_K_hat.Kbar, prob_K_hat.Kbar+3, -99)
	real colvector draws_K, list_K_0_obs, list_K_positive_obs
	real colvector P_hat_b_k
	real scalar i, K
	
	/* Draw from K */
	db_b[.,1] = (1..prob_K_hat.Kbar)'
	draws_K = Draw_K_b(prob_K_hat)
	db_b[.,2] = draws_K
									
	/* fill K with no observations */
	list_K_0_obs = selectindex(draws_K :== 0)
	if (length(list_K_0_obs)){
		for(i = 1; i <= length(list_K_0_obs); i++) {
			K = list_K_0_obs[i]
			db_b[K, 3..K+3] = J(1, K+1, 0)
		}
	}
	
	/* fill K with observations */
	list_K_positive_obs = selectindex(draws_K)
	for(i = 1; i <= length(list_K_positive_obs); i++) {
		K = list_K_positive_obs[i]
		/* Definition of the P_hat_b_k (cf. appendix D.4 DHR) */
		P_hat_b_k = complete_P(Def_Q(K) * mb[(1..K)])
		db_b[K, 3..K+3] = my_rndmultinomial(draws_K[K], P_hat_b_k)
	}
									
	return(db_b)													
}
end

capture mata: mata drop Define_In_mb()
mata:
mata set matastrict on
struct struct_In_mb_Kpind scalar Define_In_mb(struct results_prob_K_hat scalar prob_K_hat, ///
												struct raw_ML_results_Kpind scalar raw_ML_Kpind, ///
												struct options_optimization_CML scalar optionsoptCML, ///
												real scalar nb_bootstrap_distance){

	struct struct_In_mb_Kpind scalar In_mb_Kpind
	struct raw_ML_results_Kpind scalar raw_ML_Kpind_b
	real matrix db_b
	real scalar delta_hat, b, distance
	real rowvector thetas, thetas_b
	real colvector deltas_b, mb_temporary
		
	In_mb_Kpind = struct_In_mb_Kpind()
	
	if (raw_ML_Kpind.I_constrained_distribution){ /* constrained case (begin) */
		In_mb_Kpind.In = 1
		In_mb_Kpind.mb = Vector_first_moments(raw_ML_Kpind.CML.xy_hat, prob_K_hat.Kbar)
	} /* constrained case (end) */
	
	else { /* unconstrained case (begin) */
		
		thetas = Theil_oneK_Kpind(raw_ML_Kpind)
		delta_hat = thetas[2] - thetas[1]
		
		deltas_b = J(nb_bootstrap_distance, 1, .)
		mb_temporary = Vector_first_moments(raw_ML_Kpind.CML.xy_hat, prob_K_hat.Kbar)
		for(b = 1; b <= nb_bootstrap_distance; b++){
			db_b = Draw_db_bootstrap_Kpind(prob_K_hat, mb_temporary)
			raw_ML_Kpind_b = Estimation_ML_Kpind(db_b, optionsoptCML)
			thetas_b = Theil_oneK_Kpind(raw_ML_Kpind_b)
			deltas_b[b] = thetas_b[2] - thetas_b[1]
		}
		
		distance = delta_hat / (sqrt(2*log(log(raw_ML_Kpind.init_esti.nb_obs))) * sqrt(variance(deltas_b)))
		/* condition_In_1 = (sum((deltas_b :== 0)) == length(deltas_b))
		=> in this case distance = 0 since variance = 0, then it is the case In = 1 */	
		
		if (distance <= 1){
			In_mb_Kpind.In = 1
			In_mb_Kpind.mb = Project_frontier_M(raw_ML_Kpind.UML.m_tilde)
		}	
		else {
			In_mb_Kpind.In = 0
			In_mb_Kpind.mb = raw_ML_Kpind.UML.m_tilde
		}	
	} /* unconstrained case (end) */
	
	return(In_mb_Kpind)
}
end

/*OLD
capture mata: mata drop Define_In_mb_from_db()
mata:
mata set matastrict on
struct struct_In_mb_Kpind scalar Define_In_mb_from_db(real matrix db, ///
													struct options_optimization_CML scalar optionsoptCML, ///
													real scalar nb_bootstrap_distance){
	struct results_prob_K_hat scalar prob_K_hat
	struct raw_ML_results_Kpind scalar raw_ML_Kpind
	prob_K_hat = Construct_results_prob_K_hat(db)
	raw_ML_Kpind = Estimation_ML_Kpind(db, optionsoptCML)
	return(Define_In_mb(prob_K_hat, raw_ML_Kpind, optionsoptCML, nb_bootstrap_distance))													
}
end
*/
