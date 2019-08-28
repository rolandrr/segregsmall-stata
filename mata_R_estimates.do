version 14.2
/******************************************************************************/
/* this do-file implements Beta method of Rathelot 2012 (R) for the case 
components of the mixture c = 1 (at least for the moment) - more precisely
this do-file does the estimation of the indices 
cf. other mata_R_* do-files */
/******************************************************************************/

/******************************************************************************/
/* Compute segregation indices from the estimated coefficients */
/******************************************************************************/

/* Functions that compute the segregation indices based on the estimated 
coefficients of the mixture of beta-binomial distribution */

/* Duncan_beta */
/* Compute Duncan index based on the estimate of the ML of the beta binomial
model */
capture mata: mata drop Duncan_beta()
mata:
mata set matastrict on		
real scalar Duncan_beta(real rowvector parameters) {

	real scalar theta, c
	real scalar alpha, beta, mu
	real rowvector alphas, betas, lambdas, mujs
	
	/* case : c = 1 */
	if (length(parameters) == 2) {
		alpha = parameters[1]
		beta = parameters[2]
		mu = alpha / (alpha + beta)
		/*theta = ((alpha + beta)^2 / (2*alpha*beta)) * (mu * (ibeta(alpha, beta, mu) - ibeta(alpha+1, beta, mu)) + (1-mu) * (ibeta(alpha, beta+1, mu) - ibeta(alpha, beta, mu)))	
		below is an equivalent simpler expression for the Duncan index*/
		theta = (ibeta(alpha, beta, mu) - ibeta(alpha+1, beta, mu)) / (1-mu)		
	}
	/* case : c > 1 */
	else {
		c = length(parameters)/3
		alphas = parameters[1..c]
		betas = parameters[(c+1)..(2*c)]
		lambdas = parameters[(2*c+1)..(3*c)]
		mujs = alphas :/ (alphas :+ betas)
		mu = lambdas * mujs'
		/*theta = (lambdas * ((1-mu) :+ ((2*mu-1):*ibeta(alphas, betas, mu)) - (mujs:*ibeta(alphas:+1, betas, mu)) - ((1:-mujs):*(1:- ibeta(alphas, betas:+1, mu))))') / (2*mu*(1-mu))
		below is an equivalent simpler expression for the Duncan index*/
		theta = (mu*(lambdas * ibeta(alphas, betas, mu)') - (lambdas :* mujs) * ibeta(alphas:+1, betas, mu)') / (mu*(1-mu))
	}
	return(theta)
}
end	

/* Theil_beta */
/* Compute Theil index based on the estimate of the ML of the beta binomial
model */
capture mata: mata drop Theil_beta()
mata:
mata set matastrict on
real scalar Theil_beta(real rowvector parameters) {

	real scalar theta, c
	real scalar alpha, beta, mu
	real rowvector alphas, betas, lambdas, mujs
	
	/* case : c = 1 */
	if (length(parameters) == 2) {
		alpha = parameters[1]
		beta = parameters[2]
		mu = alpha / (alpha + beta)
		theta = 1 - ((mu*digamma(alpha+1) + (1-mu)*digamma(beta+1) - digamma(alpha+beta+1)) / (xlnx(mu) + xlnx(1-mu)))
	}
	/* case : c > 1 */
	else {
		c = length(parameters)/3
		alphas = parameters[1..c]
		betas = parameters[(c+1)..(2*c)]
		lambdas = parameters[(2*c+1)..(3*c)]
		mujs = alphas :/ (alphas :+ betas)
		mu = lambdas * mujs'		
		theta = 1 - (lambdas * (mujs:*digamma(alphas:+1) + (1:-mujs):*digamma(betas:+1) - digamma(alphas:+betas:+1))') / (xlnx(mu) + xlnx(1-mu))
	}
	return(theta)
}
end

/* Atkinson_beta */
/* Compute Atkinson index based on the estimate of the ML of the beta binomial
model */
capture mata: mata drop Atkinson_beta()
mata:
mata set matastrict on
real scalar Atkinson_beta(real rowvector parameters, real scalar b) {

	real scalar theta, c
	real scalar alpha, beta, mu, int_temp
	real rowvector alphas, betas, lambdas, mujs, ints
	
	/* case : c = 1 */
	if (length(parameters) == 2) {
		alpha = parameters[1]
		beta = parameters[2]
		mu = alpha / (alpha + beta)
		int_temp = exp(lngamma(alpha+b) + lngamma(beta+1-b) - lngamma(alpha) - lngamma(beta)) / (alpha + beta)
		theta = 1 - ((mu^(-b/(1-b))) / (1 - mu)) * (int_temp^(1/(1-b)))
	}
	/* case : c > 1 */
	else {
		c = length(parameters)/3
		alphas = parameters[1..c]
		betas = parameters[(c+1)..(2*c)]
		lambdas = parameters[(2*c+1)..(3*c)]
		mujs = alphas :/ (alphas :+ betas)
		mu = lambdas * mujs'
		ints = exp(lngamma(alphas:+b) + lngamma(betas:+1:-b) - lngamma(alphas) - lngamma(betas)) :/ (alphas:+betas)
		theta = 1 - ((mu^(-b/(1-b))) / (1 - mu)) * ((lambdas * ints')^(1/(1-b)))
	}
	return(theta)
}
end
/* previous expression work 
	theta = 1 - (((mu^(-b/(1-b))) / (1 - mu)) * ((fun_beta(alpha+b, beta-b+1) / fun_beta(alpha, beta))^(1/(1-b))))
*/	

/* Coworker_beta */
/* Compute Coworker index based on the estimate of the ML of the beta binomial
model */
capture mata: mata drop Coworker_beta()
mata:
mata set matastrict on
real scalar Coworker_beta(real rowvector parameters) {

	real scalar theta, c
	real scalar alpha, beta, mu, mu2
	real rowvector alphas, betas, lambdas, mujs, mu2js

	/* case : c = 1 */
	if (length(parameters) == 2) {
		alpha = parameters[1]
		beta = parameters[2]
		theta = 1 / (alpha + beta + 1)
	}
	
	/* case : c > 1 */
	else {
		c = length(parameters)/3	
		alphas = parameters[1..c]
		betas = parameters[(c+1)..(2*c)]
		lambdas = parameters[(2*c+1)..(3*c)]
		mujs = alphas :/ (alphas :+ betas)
		mu = lambdas * mujs'
		mu2js = (alphas :* (alphas:+1)) :/ ((alphas :+ betas) :* (alphas :+ betas :+ 1))
		mu2 = lambdas * mu2js'
		theta = (mu2 - mu^2) / (mu - mu^2)
	}
	return(theta)
}
end
	
/* Gini_beta */
/* Compute Gini index based on the estimate of the ML of the beta binomial
model */

/* Approximation by the trapezoidal rule */
/* alpha and beta are the two parameters of the beta
nb_points is the number of points used in the trapezoidal approximation 
for c = 1 (ceq1) */
capture mata: mata drop Approx_int_Fp2_beta_ceq1()
mata:
mata set matastrict on
real scalar Approx_int_Fp2_beta_ceq1(real scalar alpha, real scalar beta, real scalar nb_points){
	real colvector evaluations
	evaluations = ibeta(alpha, beta, rangen(0, 1, nb_points)):^2
	return((quadsum(evaluations[2..(nb_points-1)], 1)/(nb_points-1)) + (evaluations[1] + evaluations[nb_points])/(2*(nb_points-1)))
}
end

/* idem as previous function with c strictly greater than 1 (csq1)
and for the integrand defined in equation (6) of Rathelot (2012) */
capture mata: mata drop Approx_int_beta_csg1()
mata:
mata set matastrict on
real scalar Approx_int_beta_csg1(real rowvector parameters, real scalar nb_points){
	real colvector evaluations
	real colvector points
	real scalar i
	points = rangen(0, 1, nb_points)
	evaluations = J(nb_points, 1, .)
	for(i = 1; i <= nb_points; i++){
		evaluations[i] = Function_integrand_equation6(points[i], parameters)
	}
	return((quadsum(evaluations[2..(nb_points-1)], 1)/(nb_points-1)) + (evaluations[1] + evaluations[nb_points])/(2*(nb_points-1)))
}
end

capture mata: mata drop Function_f_R_section3()
mata:
mata set matastrict on
real scalar Function_f_R_section3(real scalar p, real rowvector parameters){
	real scalar alpha, beta, c, j
	real rowvector alphas, betas, lambdas
	real colvector store_fractions
	if (length(parameters) == 2) { /* case c = 1 (begin) */
		alpha = parameters[1]
		beta = parameters[2]
		return(p^(alpha-1) * (1-p)^(beta-1) / fun_beta(alpha, beta))
	} /* case c = 1 (end) */
	else { /* case : c > 1 (begin) */
		c = length(parameters)/3	
		alphas = parameters[1..c]
		betas = parameters[(c+1)..(2*c)]
		lambdas = parameters[(2*c+1)..(3*c)]
		store_fractions = J(c, 1, .)
		for(j = 1; j <= c; j++){
			store_fractions[j] = p^(alphas[j]-1) * (1-p)^(betas[j]-1) / fun_beta(alphas[j], betas[j])
		}
		return(lambdas * store_fractions)
	} /* case : c > 1 (end) */
}
end

capture mata: mata drop Function_integrand_equation6()
mata:
mata set matastrict on 
real scalar Function_integrand_equation6(real scalar p, real rowvector parameters){

	real scalar alpha, beta, c, j 
	real rowvector alphas, betas, lambdas
	real colvector store_part1, store_part2
	
	if (length(parameters) == 2) { /* case c = 1 (begin) */
		alpha = parameters[1]
		beta = parameters[2]
		return(((beta/(alpha+beta)) * Function_f_R_section3(p, (alpha, beta+1))) * ((alpha/(alpha+beta)) * ibeta(alpha+1, beta, p)))
	} /* case c = 1 (end) */
	
	else { /* case : c > 1 (begin) */
		c = length(parameters)/3	
		alphas = parameters[1..c]
		betas = parameters[(c+1)..(2*c)]
		lambdas = parameters[(2*c+1)..(3*c)]
		store_part1 = J(c, 1, .)
		store_part2 = J(c, 1, .)
		for(j = 1; j <= c; j++){
			store_part1[j] = (betas[j] / (alphas[j] + betas[j])) * Function_f_R_section3(p, (alphas[j], betas[j]+1))
			store_part2[j] = (alphas[j] / (alphas[j] + betas[j])) * ibeta(alphas[j]+1, betas[j], p)
		}
		return((lambdas * store_part1) * (lambdas * store_part2))
	} /* case : c > 1 (end) */
}
end

capture mata: mata drop Gini_beta()
mata:
mata set matastrict on 
real scalar Gini_beta(real rowvector parameters){

	real scalar mu, c
	real rowvector alphas, betas, lambdas, mujs
	real scalar nb_points, gini 
	/* to perform trapezoidal approximation 
	of int_{0}^p F_p^2(t) dt, that appears in Gini index for c = 1
	and of the quantity derived in equation (6) of Rathelot (2012) for c > 1 */
	nb_points = 10^4
	
	if (length(parameters) == 2) { /* case c = 1 (begin) */
		mu = parameters[1] / (parameters[1] + parameters[2])
		gini = (1 - mu - Approx_int_Fp2_beta_ceq1(parameters[1], parameters[2], nb_points)) / (mu * (1-mu))
	} /* case c = 1 (end) */
	
	else { /* case : c > 1 (begin) */
		c = length(parameters)/3	
		alphas = parameters[1..c]
		betas = parameters[(c+1)..(2*c)]
		lambdas = parameters[(2*c+1)..(3*c)]
		mujs = alphas :/ (alphas :+ betas)
		mu = lambdas * mujs'
		gini = 1 - 2*Approx_int_beta_csg1(parameters, nb_points)/(mu*(1-mu))
	} /* case : c > 1 (end) */
	
	if (gini > 1) gini = 1
	if (gini < 0) gini = 0
	
	return(gini)
}
end

/******************************************************************************/
/* Estimates DTACW with method beta (Rathelot 2012) */
/******************************************************************************/

/* same output matrix as in Bounds_DTACW with DHR method 
with only one colum for the estimation since one estimate (versus bounds)
with two additional rows for the Gini index (index = 5) */

capture mata: mata drop Def_expec_K_hat_from_id_pKh()
mata:
mata set matastrict on
real scalar Def_expec_K_hat_from_id_pKh(struct Info_data scalar info_data, ///
										struct results_prob_K_hat scalar prob_K_hat){
	return((1..prob_K_hat.Kbar)[info_data.list_K_positive_obs] * prob_K_hat.probabilities[info_data.list_K_positive_obs])								
}
end

capture mata: mata drop Get_db_for_oneK()
mata:
mata set matastrict on
real matrix Get_db_for_oneK(real matrix db, real scalar K){
	real matrix db_for_oneK
	real scalar k
	db_for_oneK = db[1..K,1..(K+3)]
	for(k = 1; k < K; k++){
		db_for_oneK[k,(2..(k+3))] = J(1, (k+2), 0)
	}
	return(db_for_oneK)
}
end

capture mata: mata drop Estimates_DTACWG_beta_resuncond()
mata:
mata set matastrict on
struct Results_uncond scalar Estimates_DTACWG_beta_resuncond(	real scalar b_atkinson, ///
																real matrix db, ///
																struct Info_data scalar info_data, ///
																real scalar c, ///
																struct options_optimization_CML optionsoptCML){
	
	real matrix store_estimates_K_by_K
	real rowvector weights_unit_level, weights_indi_level, parameters_hat, index_columns
	real scalar index_K, K, index_index, nb_parameters, mod_for_space
	struct results_prob_K_hat scalar prob_K_hat
	struct Results_uncond scalar results_uncond
	results_uncond = Results_uncond()

	displayas("text")
	printf("Estimation - current unit size analyzed (out of %f distinct sizes):\n", info_data.nb_K_positive_obs)
	displayflush()
	
	/* definition of the weights for the aggregated index for K random (begin) */
	prob_K_hat = Construct_results_prob_K_hat(db) /* estimation of prob_K_hat */
	weights_unit_level = (prob_K_hat.probabilities[info_data.list_K_positive_obs])'
	weights_indi_level = (1..prob_K_hat.Kbar)[info_data.list_K_positive_obs] :* ((prob_K_hat.probabilities[info_data.list_K_positive_obs])') :/ Def_expec_K_hat_from_id_pKh(info_data, prob_K_hat)
	/* definition of the weights for the aggregated index for K random (end) */
	
	/* matrix for store_infodistributionofp_beta (begin) */
	if (c == 1) nb_parameters = 2
	else nb_parameters = 3*c
	results_uncond.store_infodistributionofp_beta = J(info_data.nb_K_positive_obs, (4+nb_parameters), .)
	/* matrix for store_infodistributionofp_beta (end) */
	
	mod_for_space = 40
	
	/* definition of the estimates K by K (begin) */
	store_estimates_K_by_K = J(info_data.nb_K_positive_obs, 5, .)
	for(index_K = 1; index_K <= info_data.nb_K_positive_obs; index_K++){ /* loop over unit size K (begin) */
	
		/* follow estimation (begin) */
		displayas("text")
		if (!(mod(index_K,10))){
			if (!(mod(index_K,mod_for_space))) { /* case to put a line back */
				printf("%f \n", index_K)
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
		/* follow estimation (end) */
		
		K = info_data.list_K_positive_obs[index_K]
		parameters_hat = beta_ML(Get_db_for_oneK(db, K), 1, c, optionsoptCML)[1,.]
		store_estimates_K_by_K[index_K, 1] = Duncan_beta(parameters_hat)
		store_estimates_K_by_K[index_K, 2] = Theil_beta(parameters_hat)
		store_estimates_K_by_K[index_K, 3] = Atkinson_beta(parameters_hat, b_atkinson)
		store_estimates_K_by_K[index_K, 4] = Coworker_beta(parameters_hat)
		store_estimates_K_by_K[index_K, 5] = Gini_beta(parameters_hat)
		results_uncond.store_infodistributionofp_beta[index_K,] = (K, db[K,2], db[K,2]/info_data.nb_units_studied, c, parameters_hat) 
	} /* loop over unit size K (end) */	
	/* definition of the estimates K by K (end) */
	printf("\n"); displayflush() /* follow estimation */
	
	/* definition of the aggregated estimated bounds (begin) */
	results_uncond.estimates_ci = J(10, 3, .)
	results_uncond.estimates_ci[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4,5,5)'
	results_uncond.estimates_ci[,2] = (0,1,0,1,0,1,0,1,0,1)'
	
	for(index_index = 1; index_index <= 5; index_index++){
		index_columns = (2*index_index - 1, 2*index_index)
		results_uncond.estimates_ci[index_columns[1], 3] = weights_unit_level * store_estimates_K_by_K[,index_index]
		results_uncond.estimates_ci[index_columns[2], 3] = weights_indi_level * store_estimates_K_by_K[,index_index]
	}
	/* definition of the aggregated estimated bounds (end) */
	
	return(results_uncond)							
}
end	

capture mata: mata drop Estimates_DTACWG_beta_eci()
mata:
mata set matastrict on
real matrix Estimates_DTACWG_beta_eci(	real scalar b_atkinson, ///
										real matrix db, ///
										struct Info_data scalar info_data, ///
										real scalar c, ///
										struct options_optimization_CML optionsoptCML, ///
										real scalar I_output_for_bootstrap, ///
										real scalar I_trace_estimation){
	
	real matrix estimates_DTACWG_beta, store_estimates_K_by_K
	real rowvector weights_unit_level, weights_indi_level, parameters_hat, index_columns
	real scalar index_K, K, index_index
	struct results_prob_K_hat scalar prob_K_hat
	
	/* definition of the weights for the aggregated index for K random (begin) */
	prob_K_hat = Construct_results_prob_K_hat(db) /* estimation of prob_K_hat */
	weights_unit_level = (prob_K_hat.probabilities[info_data.list_K_positive_obs])'
	weights_indi_level = (1..prob_K_hat.Kbar)[info_data.list_K_positive_obs] :* ((prob_K_hat.probabilities[info_data.list_K_positive_obs])') :/ Def_expec_K_hat_from_id_pKh(info_data, prob_K_hat)
	/* definition of the weights for the aggregated index for K random (end) */
	
	/* definition of the estimates K by K (begin) */
	
	store_estimates_K_by_K = J(info_data.nb_K_positive_obs, 5, .)
	for(index_K = 1; index_K <= info_data.nb_K_positive_obs; index_K++){ /* loop over unit size K (begin) */
		
		
		if (I_trace_estimation) { /* trace estimation (begin) */
			displayas("text")
			if (!(mod(index_K,10))){
				printf("+")
				displayflush()
			}
			else {
				printf(".")
				displayflush()
			}		
		} /* trace estimation (end) */
		
		K = info_data.list_K_positive_obs[index_K]	
		parameters_hat = beta_ML(Get_db_for_oneK(db, K), 1, c, optionsoptCML)[1,.]		
		store_estimates_K_by_K[index_K, 1] = Duncan_beta(parameters_hat)
		store_estimates_K_by_K[index_K, 2] = Theil_beta(parameters_hat)
		store_estimates_K_by_K[index_K, 3] = Atkinson_beta(parameters_hat, b_atkinson)
		store_estimates_K_by_K[index_K, 4] = Coworker_beta(parameters_hat)
		store_estimates_K_by_K[index_K, 5] = Gini_beta(parameters_hat)
		
	} /* loop over unit size K (end) */	
	
	/* definition of the estimates K by K (end) */
	
	/* definition of the aggregated estimated bounds (begin) */
	estimates_DTACWG_beta = J(10, 3, .)
	estimates_DTACWG_beta[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4,5,5)'
	estimates_DTACWG_beta[,2] = (0,1,0,1,0,1,0,1,0,1)'
	
	for(index_index = 1; index_index <= 5; index_index++){
		index_columns = (2*index_index - 1, 2*index_index)
		estimates_DTACWG_beta[index_columns[1], 3] = weights_unit_level * store_estimates_K_by_K[,index_index]
		estimates_DTACWG_beta[index_columns[2], 3] = weights_indi_level * store_estimates_K_by_K[,index_index]
	}
	/* definition of the aggregated estimated bounds (end) */
	
	if (I_output_for_bootstrap) return((estimates_DTACWG_beta[,3]'))
	else return(estimates_DTACWG_beta)						
}
end										

/* _runc : result_uncond as an output */
capture mata: mata drop Estimates_DTACWG_beta_Kpind_runc()
mata:
mata set matastrict on
struct Results_uncond scalar Estimates_DTACWG_beta_Kpind_runc(	real scalar b_atkinson, ///
																real matrix db, ///
																struct Info_data scalar info_data, ///
																real scalar c, ///
																struct options_optimization_CML optionsoptCML){
	struct Results_uncond scalar results_uncond
	results_uncond = Results_uncond()												
	real rowvector parameters_hat
	
	parameters_hat = beta_ML(db, info_data.nb_K_positive_obs, c, optionsoptCML)[1,.]
		
	results_uncond.estimates_ci = J(10, 3, .)
	results_uncond.estimates_ci[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4,5,5)'
	results_uncond.estimates_ci[,2] = (0,1,0,1,0,1,0,1,0,1)'
	results_uncond.estimates_ci[1,3] = Duncan_beta(parameters_hat)
	results_uncond.estimates_ci[2,3] = results_uncond.estimates_ci[1,3]
	results_uncond.estimates_ci[3,3] = Theil_beta(parameters_hat)
	results_uncond.estimates_ci[4,3] = results_uncond.estimates_ci[3,3]
	results_uncond.estimates_ci[5,3] = Atkinson_beta(parameters_hat, b_atkinson)
	results_uncond.estimates_ci[6,3] = results_uncond.estimates_ci[5,3]
	results_uncond.estimates_ci[7,3] = Coworker_beta(parameters_hat)
	results_uncond.estimates_ci[8,3] = results_uncond.estimates_ci[7,3]
	results_uncond.estimates_ci[9,3] = Gini_beta(parameters_hat)
	results_uncond.estimates_ci[10,3] = results_uncond.estimates_ci[9,3]
	
	results_uncond.store_infodistributionofp_beta = (info_data.Kbar, info_data.nb_units_studied, 1, c, parameters_hat)
	
	return(results_uncond)
}
end

/* _eci : estimates_ci : return directly a matrix 10 x 3 of the estimates */
capture mata: mata drop Estimates_DTACWG_beta_Kpind_eci()
mata:
mata set matastrict on
real matrix Estimates_DTACWG_beta_Kpind_eci(real scalar b_atkinson, ///
											real matrix db, ///
											real scalar nb_K_positive_obs, ///
											real scalar c, ///
											struct options_optimization_CML optionsoptCML, ///
											real scalar I_output_for_bootstrap){
	real matrix estimates_DTACWG_beta
	real rowvector parameters_hat
	
	parameters_hat = beta_ML(db, nb_K_positive_obs, c, optionsoptCML)[1,.]
	
	estimates_DTACWG_beta = J(10, 3, .)
	estimates_DTACWG_beta[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4,5,5)'
	estimates_DTACWG_beta[,2] = (0,1,0,1,0,1,0,1,0,1)'
	estimates_DTACWG_beta[1,3] = Duncan_beta(parameters_hat)
	estimates_DTACWG_beta[2,3] = estimates_DTACWG_beta[1,3]
	estimates_DTACWG_beta[3,3] = Theil_beta(parameters_hat)
	estimates_DTACWG_beta[4,3] = estimates_DTACWG_beta[3,3]
	estimates_DTACWG_beta[5,3] = Atkinson_beta(parameters_hat, b_atkinson)
	estimates_DTACWG_beta[6,3] = estimates_DTACWG_beta[5,3]
	estimates_DTACWG_beta[7,3] = Coworker_beta(parameters_hat)
	estimates_DTACWG_beta[8,3] = estimates_DTACWG_beta[7,3]
	estimates_DTACWG_beta[9,3] = Gini_beta(parameters_hat)
	estimates_DTACWG_beta[10,3] = estimates_DTACWG_beta[9,3]
	
	if (I_output_for_bootstrap) return((estimates_DTACWG_beta[,3]'))
	else return(estimates_DTACWG_beta)
}
end
