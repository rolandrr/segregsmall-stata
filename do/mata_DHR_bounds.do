/******************************************************************************/
/* This do-file performs bounds estimation for the DHR method for
Duncan, Theil, Coworker, Atkinson for allK i.e. with K random (it also works
with a fixed K)
The general architecture is to use expression of type B.3-B.4 of DHR Appendix
or equivalently the decomposition of the quantity of interest from 
Assumption 2.1 and first estimate the brick : int h(x, m_{01])dF_p^k(x) 
denoted here int_h_Fpk */
/******************************************************************************/

/******************************************************************************/
/* Functions nu() with notation of DHR */
/******************************************************************************/
/* compute nu(u,v) for the Duncan and other index */

capture mata: mata drop compute_nu_Duncan()
mata
mata set matastrict on
real scalar compute_nu_Duncan(real scalar u, real scalar v) {
	if (u == 0) return(0)
	else if ((v == 0)|(v == 1)) return(0)
	else return(u / (2 * v * (1-v)))
}
end

capture mata: mata drop compute_nu_Theil()
mata
mata set matastrict on
real scalar compute_nu_Theil(real scalar u, real scalar v) {
	if ((v == 0)|(v == 1)) return(0)
	else return(1 - (u / (xlnx(v) + xlnx(1-v))))
}
end

/* xlnx() - extension of the function in 0 by continuity */
capture mata: mata drop xlnx()
mata:
mata set matastrict on
real matrix xlnx(real matrix x){
	return(x :* ln(x :+ (x :== 0)))
}
end	

capture mata: mata drop compute_nu_Atkinson()
mata
mata set matastrict on
real scalar compute_nu_Atkinson(real scalar u, real scalar v, real scalar b) {
	if ((v == 0)|(v == 1)) return(0)
	else {
		if ((b >= 1) | (b <= 0)) {
			errprintf("Parameter b of Atkinson is invalid - function compute_nu_atkinson()")
			return(.)
		}
		else {
			return(1 - (v^(-b/(1-b)) * (1 - v)^(-1) * u^(1/(1-b))))	
		}
	}
}
end

capture mata: mata drop compute_nu_Coworker()
mata
mata set matastrict on
real scalar compute_nu_Coworker(real scalar u, real scalar v) {
	if (u == 0) return(0)
	else if ((v == 0)|(v == 1)) return(0)
	else return(u / (v - (v^2)))
}
end


/******************************************************************************/
/* Functions int_h_Fp with notation of DHR */
/******************************************************************************/
/* compute int_h_Fp for the different segregation indices
when the distribution Fp is discrete and given by a matrix xy
with first row = support points and second row = associated masses */

capture mata: mata drop compute_int_h_Fp_dis_Duncan()
mata:
mata set matastrict on
real scalar compute_int_h_Fp_dis_Duncan(real matrix xy, real scalar m){
	return(abs(xy[1,.] :- m) * (xy[2,.]'))
}
end

capture mata: mata drop compute_int_h_Fp_dis_Theil()
mata:
mata set matastrict on
real scalar compute_int_h_Fp_dis_Theil(real matrix xy){
	return(((xy[1,.] :* ln(xy[1,.] :+ (xy[1,.] :== 0))) + ((1:-xy[1,.]) :* ln((1:-xy[1,.]) :+ ((1:-xy[1,.]) :== 0)))) * (xy[2,.]'))
}
end

capture mata: mata drop compute_int_h_Fp_dis_Atkinson()
mata:
mata set matastrict on
real scalar compute_int_h_Fp_dis_Atkinson(real matrix xy, real scalar b){
	if ((b >= 1) | (b <= 0)) {
		errprintf("Parameter b of Atkinson is invalid - function compute_int_h_Fp_dis_Atkinson()")
		return(.)
	}
	else {
		return((xy[1,.]:^b :* (1 :- xy[1,.]):^(1-b)) * (xy[2,.]'))
	}
}
end

capture mata: mata drop compute_int_h_Fp_dis_Coworker()
mata:
mata set matastrict on
real scalar compute_int_h_Fp_dis_Coworker(real matrix xy, real scalar m){
	return(m^2 + Moment_discrete_distribution(xy, 2) - 2 * m * Moment_discrete_distribution(xy, 1))
}
end


/******************************************************************************/
/* functions for m01_hat, Kp_hat, K_hat */
/******************************************************************************/

capture mata: mata drop Def_m01_hat_allK()
mata:
mata set matastrict on
real scalar Def_m01_hat_allK(	struct results_prob_K_hat scalar prob_K_hat, ///
								struct raw_ML_results_allK scalar ML_allK) {
	return((prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive])' * ML_allK.m1k_hat[ML_allK.index_K_with_nb_obs_positive])						
}
end	

capture mata: mata drop Def_expectation_Kp_hat()
mata:
mata set matastrict on
real scalar Def_expectation_Kp_hat(	struct results_prob_K_hat scalar prob_K_hat, ///
									struct raw_ML_results_allK scalar ML_allK) {
	return(((1..prob_K_hat.Kbar)[ML_allK.index_K_with_nb_obs_positive] :* (prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive]')) * ML_allK.m1k_hat[ML_allK.index_K_with_nb_obs_positive])
}
end			

capture mata: mata drop Def_expectation_K_hat()
mata:
mata set matastrict on
real scalar Def_expectation_K_hat(	struct results_prob_K_hat scalar prob_K_hat, ///
									struct raw_ML_results_allK scalar ML_allK) {
	return((1..prob_K_hat.Kbar)[ML_allK.index_K_with_nb_obs_positive] * prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive])
}
end


/******************************************************************************/
/* Quantities int_h_Fpk
whether it is weigthed or unweighted depends on the argument m 
oneK, K by K */
/******************************************************************************/

/* K per K */
capture mata: mata drop Def_int_h_Fpk_Duncan()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Duncan(struct raw_ML_results scalar raw_ML, ///
									real scalar m, real scalar nb_candidates_cr){
		
	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	
	if (raw_ML.init_esti.nb_obs) { /* case : with observations */
		if (raw_ML.I_constrained_distribution) { /* case : constrained distribution */
			int_h_Fpk[1] = compute_int_h_Fp_dis_Duncan(raw_ML.CML.xy_hat, m)
			int_h_Fpk[2] = int_h_Fpk[1]
		}
		else { /* case : unconstrained distribution */
			if(raw_ML.init_esti.K == 1) { /* case : K = 1 */
				int_h_Fpk = (0, 2*m*(1-m)) /* to have at the end 0 and 1 for the bound */
			}
			else { /* case : K > 1 */
				int_h_Fpk = int_h_Fp_hats_A21_un_Duncan_cr(raw_ML.UML.m_tilde, m, nb_candidates_cr)
			}
		}
	}
	return(int_h_Fpk)
}
end

/* K per K */
/* using second part of Theorem 2.2. - verified by Xavier that ok for the traditional definition
of the Theil index (for symetry); also checked for the Atkinson later */
capture mata: mata drop Def_int_h_Fpk_Theil()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Theil(struct raw_ML_results scalar raw_ML, real scalar m){

	struct principal_representations scalar pr
	real scalar value_pr0, value_pr1
	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	
	if (raw_ML.init_esti.nb_obs) { /* case : with observations */
		if (raw_ML.I_constrained_distribution) { /* case : constrained distribution */
			int_h_Fpk[1] = compute_int_h_Fp_dis_Theil(raw_ML.CML.xy_hat)
			int_h_Fpk[2] = int_h_Fpk[1]
		}
		else { /* case : unconstrained distribution */
			if(raw_ML.init_esti.K == 1) { /* case : K = 1 */
				int_h_Fpk = (xlnx(m) + xlnx(1-m), 0) /* to have at the end 0 and 1 for the bound */
			}
			else { /* case : K > 1 */
				pr = Principal_representations(raw_ML.UML.m_tilde)
				value_pr0 = compute_int_h_Fp_dis_Theil(Get_principal_representation(pr,0))
				value_pr1 = compute_int_h_Fp_dis_Theil(Get_principal_representation(pr,1))
				int_h_Fpk = (min((value_pr0, value_pr1)), max((value_pr0, value_pr1)))
				/* this order since the derivative of nu_Theil w.r.t. u (cf. DHR notation)
				is positive for any v in [0,1] */
			}
		}
	}
	return(int_h_Fpk)
}
end

/* K per K */
capture mata: mata drop Def_int_h_Fpk_Atkinson()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Atkinson(struct raw_ML_results scalar raw_ML, real scalar m, real scalar b){
	if ((b >= 0) & (b <= 1)){
	
		struct principal_representations scalar pr
		real scalar value_pr0, value_pr1
		real rowvector int_h_Fpk 
		int_h_Fpk = J(1, 2, .)
		
		if (raw_ML.init_esti.nb_obs) { /* case : with observations */
			if (raw_ML.I_constrained_distribution) { /* case : constrained distribution */
				int_h_Fpk[1] = compute_int_h_Fp_dis_Atkinson(raw_ML.CML.xy_hat, b)
				int_h_Fpk[2] = int_h_Fpk[1]
			}
			else { /* case : unconstrained distribution */
				if(raw_ML.init_esti.K == 1) { /* case : K = 1 */
					int_h_Fpk = ((m^(-b/(1-b))/(1-m))^((-1)*(1-b)), 0) /* to have at the end 0 and 1 for the bound */
				}
				else { /* case : K > 1 */
					pr = Principal_representations(raw_ML.UML.m_tilde)
					value_pr0 = compute_int_h_Fp_dis_Atkinson(Get_principal_representation(pr,0), b)
					value_pr1 = compute_int_h_Fp_dis_Atkinson(Get_principal_representation(pr,1), b)
					int_h_Fpk = (max((value_pr0, value_pr1)), min((value_pr0, value_pr1)))
					/* this order since the derivative of nu_Atkinson w.r.t. u (cf. DHR notation)
					is negative for any v in [0,1] */
				}
			}
		}
		return(int_h_Fpk)
	}
	else {
		errprintf("Invalid argument b in function Def_int_h_Fpk_Atkinson() \n")
		exit()
	}
}
end	

/* K per K */	
capture mata: mata drop Def_int_h_Fpk_Coworker()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Coworker(struct raw_ML_results scalar raw_ML, real scalar m){

	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	
	if (raw_ML.init_esti.nb_obs) { /* case : with observations */	
	
		if (raw_ML.I_constrained_distribution) { /* case : constrained distribution */
			int_h_Fpk[1] = compute_int_h_Fp_dis_Coworker(raw_ML.CML.xy_hat, m)
			int_h_Fpk[2] = int_h_Fpk[1]
		}
		else { /* case : unconstrained distribution */
			if(raw_ML.init_esti.K == 1) { /* case : K = 1 */
				int_h_Fpk = (0, m - m^2) /* to have at the end 0 and 1 for the bound */
			}
			else { /* case : K > 1 */
				int_h_Fpk[1] = m^2 + raw_ML.UML.m_tilde[2] - 2*m*raw_ML.UML.m_tilde[1]
				int_h_Fpk[2] = int_h_Fpk[1]
			}
		}
	}
	return(int_h_Fpk)
}
end


/******************************************************************************/
/* Bounds K by K */
/******************************************************************************/

capture mata: mata drop Get_m1hat_oneK()
mata:
mata set matastrict on
real scalar Get_m1hat_oneK(struct raw_ML_results scalar raw_ML){
	if (raw_ML.I_constrained_distribution) return(raw_ML.CML.m1_hat)
	else return(raw_ML.UML.m_tilde[1])
}
end

capture mata: mata drop Duncan_oneK()
mata:
mata set matastrict on
real rowvector Duncan_oneK(struct raw_ML_results scalar raw_ML, real scalar nb_candidates_cr){
	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	real scalar m_01_k_hat
	m_01_k_hat = Get_m1hat_oneK(raw_ML)
	int_h_Fpk = Def_int_h_Fpk_Duncan(raw_ML, m_01_k_hat, nb_candidates_cr)
	return(compute_nu_Duncan(int_h_Fpk[1], m_01_k_hat), compute_nu_Duncan(int_h_Fpk[2], m_01_k_hat))
}
end

capture mata: mata drop Theil_oneK()
mata:
mata set matastrict on 
real rowvector Theil_oneK(struct raw_ML_results scalar raw_ML){
	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	real scalar m_01_k_hat
	m_01_k_hat = Get_m1hat_oneK(raw_ML)
	int_h_Fpk = Def_int_h_Fpk_Theil(raw_ML, m_01_k_hat)
	return(compute_nu_Theil(int_h_Fpk[1], m_01_k_hat), compute_nu_Theil(int_h_Fpk[2], m_01_k_hat))
}
end

capture mata: mata drop Atkinson_oneK()
mata:
mata set matastrict on 
real rowvector Atkinson_oneK(struct raw_ML_results scalar raw_ML, real scalar b){
	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	real scalar m_01_k_hat
	m_01_k_hat = Get_m1hat_oneK(raw_ML)
	int_h_Fpk = Def_int_h_Fpk_Atkinson(raw_ML, m_01_k_hat, b)
	return(compute_nu_Atkinson(int_h_Fpk[1], m_01_k_hat, b), compute_nu_Atkinson(int_h_Fpk[2], m_01_k_hat, b))
}
end

capture mata: mata drop Coworker_oneK()
mata:
mata set matastrict on 
real rowvector Coworker_oneK(struct raw_ML_results scalar raw_ML){
	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	real scalar m_01_k_hat
	m_01_k_hat = Get_m1hat_oneK(raw_ML)
	int_h_Fpk = Def_int_h_Fpk_Coworker(raw_ML, m_01_k_hat)
	return(compute_nu_Coworker(int_h_Fpk[1], m_01_k_hat), compute_nu_Coworker(int_h_Fpk[2], m_01_k_hat))
}
end



/******************************************************************************/
/* Bounds
	implicitly for allK but it works for a database with only one fixed K
	the costly part computationally speaking is the estimation of CML
	then the computation of the different indices is instantaneous
	hence the computation of the four index and also both weighted
	and unweighted version */
/******************************************************************************/

/* structure option optimization, to concatenate options_optimization_CML and 
options_opti_21 simply */
capture mata: mata drop options_optimization()
mata:
mata set matastrict on
struct options_optimization
{
	struct options_optimization_CML scalar CML
	real scalar nb_candidates_cr
}
end

capture mata: mata drop _View_options_optimization()
mata:
mata set matastrict on
void _View_options_optimization(struct options_optimization scalar options_opt){
	printf("Struct options_optimization - 2 sub-structures \n")
	_View_options_optimization_CML(options_opt.CML)
	printf("nb_candidates_cr \n")
	options_opt.nb_candidates_cr
}
end

capture mata: mata drop Construct_options_optimization()
mata:
mata set matastrict on
struct options_optimization scalar Construct_options_optimization(
	real scalar nb_max_iter_CML, ///
	real scalar ptol, real scalar vtol, real scalar nrtol, ///
	real scalar nb_candidates_cr){

	struct options_optimization scalar options_opt
	options_opt.CML = Cons_options_optimization_CML(nb_max_iter_CML, ptol, vtol, nrtol)
	options_opt.nb_candidates_cr = nb_candidates_cr
	return(options_opt)
}
end

/* output of Bounds_DTACW* functions is a matrix :
first column is a numerotation of the index : 
1 = Duncan, 2 = Theil, 3.b = Atkinson(b), 4 = Coworker
second column is an indicator of weighted i.e. individual-level
0 = unweighted that is unit level index
1 = weighted that is individual level index
the third and fourth column is the estimated lower and upper bounds */

capture mata: mata drop Bounds_DTACW()
mata:
mata set matastrict on
real matrix Bounds_DTACW(real scalar b_atkinson, ///
						struct results_prob_K_hat scalar prob_K_hat, ///
						struct raw_ML_results_allK scalar ML_allK, ///
						struct options_optimization options_opt){
	
	real matrix bounds_DTACW, store_estimated_bounds_K_by_K
	real rowvector weights_unit_level, weights_indi_level
	real scalar index_K, K
	
	/* definition of the weights for the aggregated index for K random (begin) */
	weights_unit_level = (prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive])'
	weights_indi_level = (1..prob_K_hat.Kbar)[ML_allK.index_K_with_nb_obs_positive] :* ((prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive])') :/ Def_expectation_K_hat(prob_K_hat, ML_allK)
	/* definition of the weights for the aggregated index for K random (end) */
	
	/* definition of the estimated bounds K by K (begin) */
	store_estimated_bounds_K_by_K = J(length(ML_allK.index_K_with_nb_obs_positive), 8, .)
	for(index_K = 1; index_K <= length(ML_allK.index_K_with_nb_obs_positive); index_K++){
		K = ML_allK.index_K_with_nb_obs_positive[index_K]
		store_estimated_bounds_K_by_K[index_K, (1,2)] = Duncan_oneK(ML_allK.colv_ML[K], options_opt.nb_candidates_cr)
		store_estimated_bounds_K_by_K[index_K, (3,4)] = Theil_oneK(ML_allK.colv_ML[K])
		store_estimated_bounds_K_by_K[index_K, (5,6)] = Atkinson_oneK(ML_allK.colv_ML[K], b_atkinson)
		store_estimated_bounds_K_by_K[index_K, (7,8)] = Coworker_oneK(ML_allK.colv_ML[K])
	}
	/* definition of the estimated bounds K by K (end) */
	
	/* definition of the aggregated estimated bounds (begin) */
	bounds_DTACW = J(8, 4, .)
	bounds_DTACW[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4)'
	bounds_DTACW[,2] = (0,1,0,1,0,1,0,1)'
	
	bounds_DTACW[1,3] = weights_unit_level * store_estimated_bounds_K_by_K[,1]
	bounds_DTACW[1,4] = weights_unit_level * store_estimated_bounds_K_by_K[,2]
	bounds_DTACW[2,3] = weights_indi_level * store_estimated_bounds_K_by_K[,1]
	bounds_DTACW[2,4] = weights_indi_level * store_estimated_bounds_K_by_K[,2]
	
	bounds_DTACW[3,3] = weights_unit_level * store_estimated_bounds_K_by_K[,3]
	bounds_DTACW[3,4] = weights_unit_level * store_estimated_bounds_K_by_K[,4]
	bounds_DTACW[4,3] = weights_indi_level * store_estimated_bounds_K_by_K[,3]
	bounds_DTACW[4,4] = weights_indi_level * store_estimated_bounds_K_by_K[,4]
	
	bounds_DTACW[5,3] = weights_unit_level * store_estimated_bounds_K_by_K[,5]
	bounds_DTACW[5,4] = weights_unit_level * store_estimated_bounds_K_by_K[,6]
	bounds_DTACW[6,3] = weights_indi_level * store_estimated_bounds_K_by_K[,5]
	bounds_DTACW[6,4] = weights_indi_level * store_estimated_bounds_K_by_K[,6]
	
	bounds_DTACW[7,3] = weights_unit_level * store_estimated_bounds_K_by_K[,7]
	bounds_DTACW[7,4] = weights_unit_level * store_estimated_bounds_K_by_K[,8]
	bounds_DTACW[8,3] = weights_indi_level * store_estimated_bounds_K_by_K[,7]
	bounds_DTACW[8,4] = weights_indi_level * store_estimated_bounds_K_by_K[,8]	
	/* definition of the aggregated estimated bounds (end) */
	
	return(bounds_DTACW)
}
end

capture mata: mata drop Bounds_DTACW_from_db()
mata:
mata set matastrict on
real matrix Bounds_DTACW_from_db(real scalar b_atkinson, real matrix db, ///
						real scalar I_trace_estimation_type, ///
						struct options_optimization options_opt){
	struct results_prob_K_hat scalar prob_K_hat
	struct raw_ML_results_allK scalar ML_allK
	prob_K_hat = Construct_results_prob_K_hat(db)
	ML_allK = Estimation_ML_allK(db, prob_K_hat, 0, I_trace_estimation_type, options_opt.CML)
	return(Bounds_DTACW(b_atkinson, prob_K_hat, ML_allK, options_opt))
}
end

/*
/* OLD: with original definition of the parameters of interest in case K random
as defined in Equations (B.1) and (B.2) of DHR, Appendix B.1 Random unit size */
capture mata: mata drop Bounds_DTACW_Eq_B1_B2()
mata:
mata set matastrict on
real matrix Bounds_DTACW_Eq_B1_B2(real scalar b_atkinson, ///
								struct results_prob_K_hat scalar prob_K_hat, ///
								struct raw_ML_results_allK scalar ML_allK, ///
								struct options_optimization options_opt){
			
	real matrix bounds_DTACW
	bounds_DTACW = J(8, 4, .)
	real scalar exp_K_hat, exp_Kp_hat, m01_hat, m_indiv_level
	real scalar index, K
	real matrix store_int_h_Fpk
	real rowvector positive_prob_K_hat, weights_K_indiv_level
	
	bounds_DTACW[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4)'
	bounds_DTACW[,2] = (0,1,0,1,0,1,0,1)'
	
	/* m01_hat, exp_Kp_hat, exp_K_hat, positive_prob_K_hat */
	exp_K_hat = Def_expectation_K_hat(prob_K_hat, ML_allK)
	exp_Kp_hat = Def_expectation_Kp_hat(prob_K_hat, ML_allK)
	m01_hat = Def_m01_hat_allK(prob_K_hat, ML_allK)
	m_indiv_level = exp_Kp_hat / exp_K_hat
	positive_prob_K_hat = (prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive])'
	weights_K_indiv_level = (1..prob_K_hat.Kbar)[ML_allK.index_K_with_nb_obs_positive] :* ((prob_K_hat.probabilities[ML_allK.index_K_with_nb_obs_positive])') :/ exp_K_hat
		
	/* Duncan unit level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Duncan(ML_allK.colv_ML[K], m01_hat, options_opt.nb_candidates_cr)
	}
	bounds_DTACW[1,3] = compute_nu_Duncan(positive_prob_K_hat * store_int_h_Fpk[,1], m01_hat)
	bounds_DTACW[1,4] = compute_nu_Duncan(positive_prob_K_hat * store_int_h_Fpk[,2], m01_hat)
		
	/* Duncan individual level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Duncan(ML_allK.colv_ML[K], m_indiv_level, options_opt.nb_candidates_cr)
	}
	bounds_DTACW[2,3] = compute_nu_Duncan(weights_K_indiv_level * store_int_h_Fpk[,1], m_indiv_level)
	bounds_DTACW[2,4] = compute_nu_Duncan(weights_K_indiv_level * store_int_h_Fpk[,2], m_indiv_level)
	
	/* Theil unit level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Theil(ML_allK.colv_ML[K], m01_hat)
	}
	bounds_DTACW[3,3] = compute_nu_Theil(positive_prob_K_hat * store_int_h_Fpk[,1], m01_hat)
	bounds_DTACW[3,4] = compute_nu_Theil(positive_prob_K_hat * store_int_h_Fpk[,2], m01_hat)
	
	/* Theil individual level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Theil(ML_allK.colv_ML[K], m_indiv_level)
	}
	bounds_DTACW[4,3] = compute_nu_Theil(weights_K_indiv_level * store_int_h_Fpk[,1], m_indiv_level)
	bounds_DTACW[4,4] = compute_nu_Theil(weights_K_indiv_level * store_int_h_Fpk[,2], m_indiv_level)
	
	/* Atkinson(b=0.5) unit level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Atkinson(ML_allK.colv_ML[K], m01_hat, b_atkinson)
	}
	bounds_DTACW[5,3] = compute_nu_Atkinson(positive_prob_K_hat * store_int_h_Fpk[,1], m01_hat, b_atkinson)
	bounds_DTACW[5,4] = compute_nu_Atkinson(positive_prob_K_hat * store_int_h_Fpk[,2], m01_hat, b_atkinson)
	
	/* Atkinson(b=0.5) individual level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Atkinson(ML_allK.colv_ML[K], m_indiv_level, b_atkinson)
	}
	bounds_DTACW[6,3] = compute_nu_Atkinson(weights_K_indiv_level * store_int_h_Fpk[,1], m_indiv_level, b_atkinson)
	bounds_DTACW[6,4] = compute_nu_Atkinson(weights_K_indiv_level * store_int_h_Fpk[,2], m_indiv_level, b_atkinson)
	
	/* Coworker unit level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Coworker(ML_allK.colv_ML[K], m01_hat)
	}
	bounds_DTACW[7,3] = compute_nu_Coworker(positive_prob_K_hat * store_int_h_Fpk[,1], m01_hat)
	bounds_DTACW[7,4] = compute_nu_Coworker(positive_prob_K_hat * store_int_h_Fpk[,2], m01_hat)
	
	/* Coworker individual level */
	store_int_h_Fpk = J(length(ML_allK.index_K_with_nb_obs_positive), 2, .)
	for(index = 1; index <= length(ML_allK.index_K_with_nb_obs_positive); index++){
		K = ML_allK.index_K_with_nb_obs_positive[index]
		store_int_h_Fpk[index,] =  Def_int_h_Fpk_Coworker(ML_allK.colv_ML[K], m_indiv_level)
	}
	bounds_DTACW[8,3] = compute_nu_Coworker(weights_K_indiv_level * store_int_h_Fpk[,1], m_indiv_level)
	bounds_DTACW[8,4] = compute_nu_Coworker(weights_K_indiv_level * store_int_h_Fpk[,2], m_indiv_level)
	
	return(bounds_DTACW)
}
end

capture mata: mata drop Bounds_DTACW_Eq_B1_B2_from_db()
mata:
mata set matastrict on
real matrix Bounds_DTACW_Eq_B1_B2_from_db(real scalar b_atkinson, real matrix db, ///
										real scalar I_trace_estimation_type, ///
										struct options_optimization options_opt){
	struct results_prob_K_hat scalar prob_K_hat
	struct raw_ML_results_allK scalar ML_allK
	prob_K_hat = Construct_results_prob_K_hat(db)
	ML_allK = Estimation_ML_allK(db, prob_K_hat, 0, I_trace_estimation_type, options_opt.CML)
	return(Bounds_DTACW_Eq_B1_B2(b_atkinson, prob_K_hat, ML_allK, options_opt))
}
end
*/
