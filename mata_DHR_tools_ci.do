version 14.2
/******************************************************************************/
/* This do-file creates structures and functions used later to compute
confidence interval (ci) for Duncan, Theil, Coworker, Atkinson for allK i.e. 
with K random */
/******************************************************************************/

/******************************************************************************/
/* Inputs for bootstrap */
/******************************************************************************/

capture mata: mata drop struct_In_Pb()
mata:
mata set matastrict on
struct struct_In_Pb {
	real colvector Pb
	real scalar In
}
end

/* complete_P() 
when vector of probabilities for X = 1..K, construct the complete vector i.e. 
also with probability that X = 0 */
capture mata: mata drop complete_P()
mata:
mata set matastrict on
real colvector complete_P(real colvector P_incomplete) {
	real scalar sum_proba, min_proba
	
	sum_proba = quadsum(P_incomplete,1)
	min_proba = min(P_incomplete)
	if ((sum_proba <= 1) & (min_proba >= 0)) {
		return(((1 - sum_proba) \ P_incomplete))
	}
	else {
		errprintf("The vector to be completed is not an adequate vector of probabilities - function complete_P()\n")
		exit()
	}
}
end

/* Get first l moments from a discrete distribution */
capture mata: mata drop Vector_first_moments()
mata:
mata set matastrict on
real colvector Vector_first_moments(real matrix xy, real scalar l){
	real colvector moments_xy
	real scalar p
	moments_xy = J(l, 1, .)
	for(p = 1; p <= l; p++){
		moments_xy[p] = Moment_discrete_distribution(xy, p)
	}
	return(moments_xy)
}
end

/* distance euclidian two vectors */
capture mata: mata drop Distance_euclidian()
mata:
mata set matastrict on
real scalar Distance_euclidian(real colvector m1, real colvector m2){
	real colvector diff
	if ((length(m1) == length(m2))){
		diff = m1 - m2
		return(sqrt((diff')*diff))
	}
	else {
		errprintf("Different dimensions in Distance_euclidian() - exit \n")
		exit()
	}
}
end


/* Project_frontier_M */
capture mata: mata drop Project_frontier_M()
mata:
mata set matastrict on
real colvector Project_frontier_M(real colvector m_in_int_M){

	real colvector m_projected, m_line, m_curve, m_pr0, m_pr1
	real scalar m1_line, m2_line, m1_curve, m2_curve, appro_slope
	real scalar K, tol_frontier, threshold, nb_ite_max, iteration, add_slope_curve
	struct principal_representations scalar pr
	
	/* options for manual projection in case K = 2 */
	tol_frontier = 10^(-9)
	threshold = 10^(-5)
	nb_ite_max = 0.3*(1/threshold)
	add_slope_curve = -0.5
	
	K = length(m_in_int_M)
	
	if (K == 1) { /* case K = 1 */
		/* extreme case, frontier of M has just {0} and {1}
		pick the closet */
		m_projected = ((m_in_int_M > 0.5) ? 1 : 0)
	}
	else if (K == 2) { /* case K = 2 */
		/* manual projection to frontier of M */
	
		/* project on line y = x */
		m1_line = m_in_int_M[1]
		m2_line = m_in_int_M[2]
		iteration = 0
		while (((m1_line - m2_line) >= tol_frontier) & (iteration <= nb_ite_max)){
			iteration = iteration + 1
			m2_line = m2_line + threshold
			m1_line = m1_line - threshold
		}
		m_line = (m1_line \ m2_line)
		
		/* projet on curve y = x^2 */
		m1_curve = m_in_int_M[1]
		m2_curve = m_in_int_M[2]
		appro_slope = 1/(2*(m1_curve*(1+add_slope_curve)))
		iteration = 0
		while (((m2_curve - m1_curve^2) >= tol_frontier) & (iteration <= nb_ite_max)){
			iteration = iteration + 1
			m2_curve = m2_curve - appro_slope*threshold
			m1_curve = m1_curve + threshold
		}
		m_curve = (m1_curve \ m2_curve)
		
		/* choose the closest among this two possibility */
		if (Distance_euclidian(m_in_int_M, m_line) >= Distance_euclidian(m_in_int_M, m_curve)){
			m_projected = m_curve
		}
		else {
			m_projected = m_line
		}	
	}
	else { /* case K > 2 */
		/* use principal representation to get something on frontier of M
		and only difference in the last moment, not sure actually to do the projection
		but as K increases, difference between m_in_int_M and m_projected decreases
		and above all it seems difficult to do better with guarantee and in reasonable
		amount of time in Mata */
		
		pr = Principal_representations(m_in_int_M[1..(K-1)])
		m_pr0 = Vector_first_moments(Get_principal_representation(pr,0), K)
		m_pr1 = Vector_first_moments(Get_principal_representation(pr,1), K)
		if (Distance_euclidian(m_in_int_M, m_pr0) >= Distance_euclidian(m_in_int_M, m_pr1)){
			m_projected = m_pr1
		}
		else {
			m_projected = m_pr0
		}
	}
	
	return(m_projected)
}
end

/* Get P_b : P_bootstrap and I_n for oneK
In Matlab code, the distance is based on the Theil bounds 
Choice to be determined here */
capture mata: mata drop Define_In_Pb_oneK()
mata:
mata set matastrict on
struct struct_In_Pb scalar Define_In_Pb_oneK(struct raw_ML_results scalar raw_ML, ///
											struct options_optimization_CML scalar optionsoptCML, ///
											real scalar nb_bootstrap_distance){
			
	struct struct_In_Pb scalar In_Pb
	struct raw_ML_results scalar raw_ML_b
	real rowvector thetas, thetas_b, data_Nks_b
	real scalar delta_hat, b, distance
	real colvector deltas_b

	if (raw_ML.I_constrained_distribution) { /* constrained case */
		In_Pb.In = 1
		In_Pb.Pb = complete_P(raw_ML.CML.P_hat)
	}
	else { /* unconstrained case */
	
		thetas = Theil_oneK(raw_ML)
		/*thetas = Duncan_oneK(raw_ML)*/
		/*thetas = Atkinson_oneK(raw_ML, 0.5)*/
		delta_hat = thetas[2] - thetas[1]
		
		deltas_b = J(nb_bootstrap_distance, 1, .)
		for(b = 1; b <= nb_bootstrap_distance; b++){
			data_Nks_b = (raw_ML.init_esti.K, raw_ML.init_esti.nb_obs, my_rndmultinomial(raw_ML.init_esti.nb_obs, complete_P(raw_ML.CML.P_hat)))
			raw_ML_b = Estimation_ML(data_Nks_b, optionsoptCML)
			thetas_b = Theil_oneK(raw_ML_b)
			/*thetas_b = Duncan_oneK(raw_ML_b)*/
			/*thetas_b = Atkinson_oneK(raw_ML_b, 0.5)*/
			deltas_b[b] = thetas_b[2] - thetas_b[1]
		}
		distance = delta_hat / (sqrt(2*log(log(raw_ML.init_esti.nb_obs))) * sqrt(variance(deltas_b)))
		/* condition_In_1 = (sum((deltas_b :== 0)) == length(deltas_b))
		=> in this case distance = 0 since variance = 0, then it is the case In = 1 */
		
		if (distance <= 1){
			In_Pb.In = 1
			In_Pb.Pb = complete_P(Def_Q(raw_ML.init_esti.K) * Project_frontier_M(raw_ML.UML.m_tilde))
		}
		else {
			In_Pb.In = 0
			In_Pb.Pb = complete_P(raw_ML.UML.P_tilde)
		}
	}
	return(In_Pb)
}
end

capture mata: mata drop _View_struct_In_Pb()
mata:
mata set matastrict on
void _View_struct_In_Pb(struct struct_In_Pb scalar In_Pb){
	printf("Struct struct_In_Pb - 2 elements \n")
	printf("Pb = \n")
	In_Pb.Pb
	printf("In = \n")
	In_Pb.In
}
end
	

	
/******************************************************************************/
/* Empirical distribution of K */
/******************************************************************************/

/* Draw K in its empirical distribution */
capture mata: mata drop Draw_K_b()
mata:
mata set matastrict on
real colvector Draw_K_b(struct results_prob_K_hat scalar prob_K_hat){
	return(my_rndmultinomial(prob_K_hat.nb_obs_allK, prob_K_hat.probabilities)')
}
end

/* my_rndmultinomial() - random multinomial generator */
capture mata: mata drop my_rndmultinomial()
mata:
mata set matastrict on
real rowvector my_rndmultinomial(real scalar nb_draws, real colvector P) {
	real rowvector draws, result
	real scalar i
	draws = rdiscrete(1, nb_draws, P) 
	result = J(1, length(P), .)
	for(i = 1; i <= length(P); i++) {
		result[i] = quadsum((draws :== i),1)
	}
	return(result)
}
end
	
	
/******************************************************************************/
/* Iteration of bootstrap */
/******************************************************************************/

capture mata: mata drop struct_In_Pb_allK()
mata:
mata set matastrict on
struct struct_In_Pb_allK {
	struct struct_In_Pb colvector colv_In_Pb
}
end

capture mata: mata drop Define_In_Pb_allK()
mata:
mata set matastrict on
struct struct_In_Pb_allK scalar Define_In_Pb_allK(struct results_prob_K_hat scalar prob_K_hat, ///
												struct raw_ML_results_allK scalar ML_allK, ///
												struct options_optimization_CML scalar optionsoptCML, ///
												real scalar nb_bootstrap_distance){

	struct struct_In_Pb_allK scalar	In_Pb_allK										
	real scalar index, K
	
	In_Pb_allK.colv_In_Pb = struct_In_Pb(prob_K_hat.Kbar)
	
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		In_Pb_allK.colv_In_Pb[K] = Define_In_Pb_oneK(ML_allK.colv_ML[K], optionsoptCML, nb_bootstrap_distance)
	}
	
	return(In_Pb_allK)
}
end

/*OLD
capture mata: mata drop Define_In_Pb_allK_from_db()
mata:
mata set matastrict on
struct struct_In_Pb_allK scalar Define_In_Pb_allK_from_db(real matrix db, ///
														struct options_optimization_CML scalar optionsoptCML, ///
														real scalar nb_bootstrap_distance){

	struct results_prob_K_hat scalar prob_K_hat
	struct raw_ML_results_allK scalar ML_allK
	prob_K_hat = Construct_results_prob_K_hat(db)
	ML_allK = Estimation_ML_allK(db, prob_K_hat, 0, 0, optionsoptCML)	

	return(Define_In_Pb_allK(prob_K_hat, ML_allK, optionsoptCML, nb_bootstrap_distance))
}
end
*/

capture mata: mata drop _View_struct_In_Pb_allK()
mata:
mata set matastrict on
void _View_struct_In_Pb_allK(struct struct_In_Pb_allK scalar In_Pb_allK	, real scalar K){
	_View_struct_In_Pb(In_Pb_allK.colv_In_Pb[K])
}
end

capture mata: mata drop Draw_db_bootstrap()
mata:
mata set matastrict on
real matrix Draw_db_bootstrap(struct results_prob_K_hat scalar prob_K_hat, ///
							struct raw_ML_results_allK scalar ML_allK, ///
							struct struct_In_Pb_allK scalar	In_Pb_allK) {	
	real matrix db_b
	db_b = J(ML_allK.Kbar, ML_allK.Kbar+3, -99)
	
	real colvector draws_K, list_K_0_obs, list_K_positive_obs
	real scalar i, K
	
	/* Draw from K */
	db_b[.,1] = (1..ML_allK.Kbar)'
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
		db_b[K, 3..K+3] = my_rndmultinomial(draws_K[K], In_Pb_allK.colv_In_Pb[K].Pb)
	}
	
	return(db_b)
}
end


/******************************************************************************/
/* options for bootstrap */
/******************************************************************************/

capture mata: mata drop options_ci_bootstrap()
mata:
mata set matastrict on
struct options_ci_bootstrap
{
	real scalar nb_bootstrap_distance
	real scalar nb_bootstrap_repetition
}
end

capture mata: mata drop Cons_options_ci_bootstrap()
mata:
mata set matastrict on
struct options_ci_bootstrap scalar Cons_options_ci_bootstrap(real scalar nb_bootstrap_distance_arg, ///
															real scalar nb_bootstrap_repetition_arg){
	struct options_ci_bootstrap scalar optionscib
	optionscib.nb_bootstrap_distance = nb_bootstrap_distance_arg
	optionscib.nb_bootstrap_repetition = nb_bootstrap_repetition_arg
	return(optionscib)
}
end


/******************************************************************************/
/* options for optimization and bootstrap */
/******************************************************************************/

/* structure options_optimization_cib,
to concatenate options_optimization and options_ci_bootstrap */
capture mata: mata drop options_optimization_cib()
mata:
mata set matastrict on
struct options_optimization_cib 
{
	struct options_optimization scalar optimization
	struct options_ci_bootstrap scalar cib
}
end

capture mata: mata drop Cons_options_optimization_cib()
mata:
mata set matastrict on
struct options_optimization_cib scalar Cons_options_optimization_cib(struct options_optimization scalar optimization_arg, ///
														struct options_ci_bootstrap scalar cib_arg){
	struct options_optimization_cib scalar optionsoptcib
	optionsoptcib.optimization = optimization_arg
	optionsoptcib.cib = cib_arg
	return(optionsoptcib)
}
end


/******************************************************************************/
/* Define I_n (for allK, for oneK or for K and p assumed independent)
specific to each index */
/******************************************************************************/
/* for a given index i.e. type and unit or individual level, compute I_n
that determine whether we use the boundary or interior case CI */
capture mata: mata drop Define_In_ci()
mata:
mata set matastrict on
real scalar Define_In_ci(real rowvector bounds_hat, real matrix bounds_bootstrap, real scalar nb_obs, real scalar I_all_Fpk_constrained){
	/* arguments (following DHR notation, part 3.2):
	bounds_hat[1] = theta_low_hat
	bounds_hat[2] = theta_up_hat 
	bounds_bootstrap[,1] = vector of the bootstrap theta_low_hat
	bounds_bootstrap[,2] = vector of the bootstrap theta_up_hat 
	I_all_Fpk_constrained = indicator whether all the Fpk are constrained, in this case the bounds are equal and I_n = 1 
	value :
	indicator I_n as defined in the paper */
	
/* NB: OLD : version with k_n as defined in the paper but problem, cf. correction avec Xavier July 2018 and my paper version of the paper
modified 2018/10/12, may change previous results, in particular about some problem of inference for the Duncan perhaps? */
/*
	real scalar In_ci_allK
	In_ci_allK = (nb_obs*(bounds_hat[2]-bounds_hat[1])*sqrt(variance(bounds_bootstrap[,2] - bounds_bootstrap[,1])) / sqrt(2*log(log(nb_obs))) <= 1)
	return(In_ci_allK)
*/	
	real scalar distance_dn, variance_delta_hat_bootstrap
	
	if (I_all_Fpk_constrained){ /* constrained case (begin) */
		return(1)
	} /* constrained case (end) */
	
	else { /* unconstrained case (begin) */
		if (bounds_hat[1] == bounds_hat[2]){ 
			return(1)
		} 
		else { 
			variance_delta_hat_bootstrap = variance((bounds_bootstrap[,2] :- bounds_bootstrap[,1]))
			if (variance_delta_hat_bootstrap == 0){
				return(1)
			}
			else {
				distance_dn = (bounds_hat[2]-bounds_hat[1]) / sqrt(2*log(log(nb_obs))*variance_delta_hat_bootstrap)
				return((distance_dn <= 1))
			}
		}	
	} /* unconstrained case (end) */
}
end

/******************************************************************************/
/* CI_ : compute the different CI, 
return a rowvector: first column = CI_low, second = CI_up*/
/******************************************************************************/

/* my_empirical_quantile()
get the quantile tau from an empirical distribution 
i.e. simply the tau-th % number in ascending order */
capture mata: mata drop my_empirical_quantile()
mata:
mata set matastrict on
real scalar my_empirical_quantile(real colvector X, real scalar tau) {
	if ((tau >= 1) | (tau <= 0)) {
		errprintf("The quantile entered is not valid (outside (0,1)) - function my_empirical_quantile() ")
		exit()
	}
	else {
		_sort(X,1)
		return(X[ceil(tau * length(X))])
	}	
}
end

/* OLD: Unnecessary since auto-normalization of the number of units in CI_interior
capture mata: mata drop CI_interior_with_rootn()
mata:
mata set matastrict on
real rowvector CI_interior_with_rootn(	real rowvector bounds_hat, real matrix bounds_bootstrap, ///
										real scalar nb_obs, real scalar alpha){
	real scalar CI_low, CI_up						
	real colvector T_low_stars, T_up_stars
	real scalar quantile_T_low_stars, quantile_T_up_stars
	
	
	T_low_stars = sqrt(nb_obs) :* (bounds_bootstrap[,1] :- bounds_hat[1])
	T_up_stars = sqrt(nb_obs) :* (bounds_bootstrap[,2] :- bounds_hat[2])
	quantile_T_low_stars = my_empirical_quantile(T_low_stars, 1 - alpha)
	quantile_T_up_stars = my_empirical_quantile(T_up_stars, alpha)
	
	CI_low = bounds_hat[1] - (quantile_T_low_stars / sqrt(nb_obs))
	CI_up = bounds_hat[2] - (quantile_T_up_stars / sqrt(nb_obs))
	
	return((CI_low, CI_up))
}
end
*/

capture mata: mata drop CI_interior()
mata:
mata set matastrict on
real rowvector CI_interior(real rowvector bounds_hat, real matrix bounds_bootstrap, real scalar alpha){
	real scalar CI_low, CI_up						
	real colvector T_low_stars_without_rootn, T_up_stars_without_rootn
	
	T_low_stars_without_rootn = (bounds_bootstrap[,1] :- bounds_hat[1])
	T_up_stars_without_rootn = (bounds_bootstrap[,2] :- bounds_hat[2])	
	CI_low = max((0, bounds_hat[1] - my_empirical_quantile(T_low_stars_without_rootn, 1 - alpha)))
	CI_up = min((1, bounds_hat[2] - my_empirical_quantile(T_up_stars_without_rootn, alpha)))
	return((CI_low, CI_up))
}
end

/* OLD: Unnecessary since auto-normalization of the number of units in CI_interior
capture mata: mata drop CI_boundary_with_rootn()
mata:
mata set matastrict on
real rowvector CI_boundary_with_rootn( real rowvector bounds_hat, real matrix bounds_bootstrap, ///
							real scalar nb_obs, real scalar alpha){
	real scalar CI_low, CI_up
	real colvector T_s_stars
	real scalar quantile_T_s_stars

	T_s_stars = sqrt(nb_obs) :* abs(((bounds_bootstrap[,1]:+bounds_bootstrap[,2]):/2) :- ((bounds_hat[1]+bounds_hat[2])/2))
	quantile_T_s_stars = my_empirical_quantile(T_s_stars, 1 - alpha)
	CI_low = bounds_hat[1] - quantile_T_s_stars / sqrt(nb_obs)
	CI_up = bounds_hat[2] + quantile_T_s_stars / sqrt(nb_obs)
	
	return((CI_low, CI_up))
}
end
*/

capture mata: mata drop CI_boundary()
mata:
mata set matastrict on
real rowvector CI_boundary(real rowvector bounds_hat, real matrix bounds_bootstrap, real scalar alpha){
	real scalar CI_low, CI_up
	real colvector T_s_stars_without_rootn
	real scalar quantile_T_s_stars

	T_s_stars_without_rootn = abs(((bounds_bootstrap[,1]:+bounds_bootstrap[,2]):/2) :- ((bounds_hat[1]+bounds_hat[2])/2))
	quantile_T_s_stars = my_empirical_quantile(T_s_stars_without_rootn, 1 - alpha)
	CI_low = max((0, bounds_hat[1] - quantile_T_s_stars))
	CI_up = min((1, bounds_hat[2] + quantile_T_s_stars))
	return((CI_low, CI_up))
}
end

/* for a given index, estimated bounds, bootstrap bounds, alpha and I_n,
compute CI_1_(1-alpha)
CI_1 to follow the notation of DHR 
CI_2 (cf. Appendix B.2 with K random is treated similarly with bootstrapping
the aggregated index */
capture mata: mata drop CI_1()
mata:
mata set matastrict on
real rowvector CI_1(real scalar I_n, real rowvector bounds_hat, ///
					real matrix bounds_bootstrap, real scalar alpha){			
	if (I_n) return(CI_boundary(bounds_hat, bounds_bootstrap, alpha))
	else return(CI_interior(bounds_hat, bounds_bootstrap, alpha))
}
end	
