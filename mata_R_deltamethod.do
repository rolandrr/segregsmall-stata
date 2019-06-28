/******************************************************************************/
/* this do-file implements inference (confidence intervals) for the Beta method 
of Rathelot 2012 (R), or at least for the number of components of the 
mixture c = 1 for the moment 
R 2012 suggests two possible methods for inference :
bootstrap
delta method
This do-files does the delta-method part */
/******************************************************************************/

/******************************************************************************/
/* Get numerical derivatives 
for some cases (example coworker with c = 1), there is a simple analytical 
derivative but simpler to use numerical derivatives for most of cases */
/******************************************************************************/

/* eval_Duncan_beta */
/* evaluator for computing numerical derivative */
capture mata: mata drop D
capture mata: mata drop eval_Duncan_beta()
mata:
mata set matastrict on		
void eval_Duncan_beta(real rowvector parameters, theta) {
	theta = Duncan_beta(parameters)
}
end

/* get_d_Duncan_beta */
/* return jacobian for Duncan_beta evaluated at `parameters'*/
capture mata: mata drop get_d_Duncan_beta()
mata:
mata set matastrict on
real rowvector get_d_Duncan_beta(real rowvector parameters) {
	transmorphic D
	D = deriv_init()
	deriv_init_evaluator(D, &eval_Duncan_beta())
	deriv_init_evaluatortype(D, "d")
	deriv_init_params(D, parameters)
	return(deriv(D, 1))
}
end

/* eval_Theil_beta */
/* evaluator for computing numerical derivative */
capture mata: mata drop D
capture mata: mata drop eval_Theil_beta()
mata:
mata set matastrict on
void eval_Theil_beta(real rowvector parameters, theta) {
	theta = Theil_beta(parameters)
}
end

/* get_d_Theil_beta */
/* return jacobian for Theil_beta evaluated at `parameters'*/
capture mata: mata drop get_d_Theil_beta()
mata:
mata set matastrict on
real rowvector get_d_Theil_beta(real rowvector parameters) {
	transmorphic D
	D = deriv_init()
	deriv_init_evaluator(D, &eval_Theil_beta())
	deriv_init_evaluatortype(D, "d")
	deriv_init_params(D, parameters)
	return(deriv(D, 1))
}
end

/* eval_Atkinson_beta */
/* evaluator for computing numerical derivative */
capture mata: mata drop D
capture mata: mata drop eval_Atkinson_beta()
mata:
mata set matastrict on
void eval_Atkinson_beta(real rowvector parameters, real scalar b, theta) {
	theta = Atkinson_beta(parameters, b)
}
end

/* get_d_Atkinson_beta */
/* return jacobian for Atkinson_beta evaluated at `parameters' and with b parameter `b'*/
capture mata: mata drop get_d_Atkinson_beta()
mata:
mata set matastrict on
real rowvector get_d_Atkinson_beta(real rowvector parameters, real scalar b) {
	transmorphic D
	D = deriv_init()
	deriv_init_evaluator(D, &eval_Atkinson_beta())
	deriv_init_evaluatortype(D, "d")
	deriv_init_params(D, parameters)
	deriv_init_argument(D, 1, b)
	return(deriv(D, 1))
}
end

/* eval_Coworker_beta */
/* evaluator for computing numerical derivative (for c > 1 ; c = 1, 
easy analytical solution */
capture mata: mata drop D
capture mata: mata drop eval_Coworker_beta()
mata:
mata set matastrict on
void eval_Coworker_beta(real rowvector parameters, theta) {
	theta = Coworker_beta(parameters)
}
end

/* get_d_Coworker_beta */
/* return jacobian for Coworker_beta evaluated at `parameters'*/
capture mata: mata drop get_d_Coworker_beta()
mata:
mata set matastrict on
real rowvector get_d_Coworker_beta(real rowvector parameters) {
	
	real scalar deriv_analytical
	transmorphic D
	
	/* case : c = 1 */
	if (length(parameters) == 2) {
		deriv_analytical = -1 / (quadsum(parameters,1) + 1)^2
		return((deriv_analytical, deriv_analytical))
	}
	
	/* case : c > 1 */
	else {
		D = deriv_init()
		deriv_init_evaluator(D, &eval_Coworker_beta())
		deriv_init_evaluatortype(D, "d")
		deriv_init_params(D, parameters)
		return(deriv(D, 1))
	}
}
end

/* eval_Gini_beta */
/* evaluator for computing numerical derivative */
capture mata: mata drop D
capture mata: mata drop eval_Gini_beta()
mata:
mata set matastrict on		
void eval_Gini_beta(real rowvector parameters, theta) {
	theta = Gini_beta(parameters)
}
end

/* get_d_Gini_beta */
/* return jacobian for Gini_beta evaluated at `parameters'*/
capture mata: mata drop get_d_Gini_beta()
mata:
mata set matastrict on
real rowvector get_d_Gini_beta(real rowvector parameters) {
	transmorphic D
	D = deriv_init()
	deriv_init_evaluator(D, &eval_Gini_beta())
	deriv_init_evaluatortype(D, "d")
	deriv_init_params(D, parameters)
	return(deriv(D, 1))
}
end

	
/******************************************************************************/
/* construct confidence-interval with delta-method */
/******************************************************************************/
/* only for unconditional case and K and p assumed independent */

/* NB: there is a first transformation from the parameters of optimize
of beta_ML => take the exponential in case c = 1 for the alpha and beta
parameters 
this initial step has to be taken into account in the delta method */

capture mata: mata drop CI_from_asymptotic_normality()
mata:
mata set matastrict on
real rowvector CI_from_asymptotic_normality(estimator, nb_obs, standard_error, alpha) {
	real scalar demi_length
	demi_length = invnormal(1-(alpha/2))*standard_error/sqrt(nb_obs)
	return((estimator - demi_length, estimator + demi_length))
}
end

/* _runc : result_uncond as an output */
capture mata: mata drop Estici_dm_DTACWG_beta_Kpind_runc()
mata:
mata set matastrict on
struct Results_uncond scalar Estici_dm_DTACWG_beta_Kpind_runc(	real scalar b_atkinson, ///
																real matrix db, ///
																struct Info_data scalar info_data, ///
																real scalar c, ///
																struct options_optimization_CML optionsoptCML, ///
																| real scalar alpha){
	real scalar nb_arg_option; nb_arg_option = 6
	
	struct Results_uncond scalar results_uncond
	results_uncond = Results_uncond()		
	real matrix ML_results, variance_parameter_tilde_hat
	real rowvector parameters_hat, jacobian
	real colvector store_standard_error_hat
	real scalar index_seg, index_row_1, index_row_2
	
	ML_results = beta_ML(db, info_data.nb_K_positive_obs, c, optionsoptCML)
	parameters_hat = ML_results[1,]
	variance_parameter_tilde_hat = ML_results[(2..(length(parameters_hat)+1)),]
	
	results_uncond.estimates_ci = J(10, ((args() == nb_arg_option) ? 11 : 9), .)
	store_standard_error_hat = J(5, 1, .)
	
	/* Duncan */
	results_uncond.estimates_ci[1,3] = Duncan_beta(parameters_hat)
	jacobian = parameters_hat :* get_d_Duncan_beta(parameters_hat)
	store_standard_error_hat[1] = sqrt(jacobian * variance_parameter_tilde_hat * (jacobian'))
	/* Theil */
	results_uncond.estimates_ci[3,3] = Theil_beta(parameters_hat)
	jacobian = parameters_hat :* get_d_Theil_beta(parameters_hat)
	store_standard_error_hat[2] = sqrt(jacobian * variance_parameter_tilde_hat * (jacobian'))
	/* Atkinson */
	results_uncond.estimates_ci[5,3] = Atkinson_beta(parameters_hat, b_atkinson)
	jacobian = parameters_hat :* get_d_Atkinson_beta(parameters_hat, b_atkinson)
	store_standard_error_hat[3] = sqrt(jacobian * variance_parameter_tilde_hat * (jacobian'))
	/* Coworker */
	results_uncond.estimates_ci[7,3] = Coworker_beta(parameters_hat)
	jacobian = parameters_hat :* get_d_Coworker_beta(parameters_hat)
	store_standard_error_hat[4] = sqrt(jacobian * variance_parameter_tilde_hat * (jacobian'))
	/* Gini */
	results_uncond.estimates_ci[9,3] = Gini_beta(parameters_hat)
	jacobian = parameters_hat :* get_d_Gini_beta(parameters_hat)
	store_standard_error_hat[5] = sqrt(jacobian * variance_parameter_tilde_hat * (jacobian'))
	
	for(index_seg = 1; index_seg <= 5; index_seg++){
		index_row_1 = 2*index_seg-1
		index_row_2 = 2*index_seg
		results_uncond.estimates_ci[index_row_1,(4,5)] = CI_from_asymptotic_normality(results_uncond.estimates_ci[index_row_1,3], info_data.nb_units_studied, store_standard_error_hat[index_seg], 0.01)
		results_uncond.estimates_ci[index_row_1,(6,7)] = CI_from_asymptotic_normality(results_uncond.estimates_ci[index_row_1,3], info_data.nb_units_studied, store_standard_error_hat[index_seg], 0.05)
		results_uncond.estimates_ci[index_row_1,(8,9)] = CI_from_asymptotic_normality(results_uncond.estimates_ci[index_row_1,3], info_data.nb_units_studied, store_standard_error_hat[index_seg], 0.10)
		if (args() == nb_arg_option) results_uncond.estimates_ci[index_row_1,(10,11)] = CI_from_asymptotic_normality(results_uncond.estimates_ci[index_row_1,3], info_data.nb_units_studied, store_standard_error_hat[index_seg], alpha)
		results_uncond.estimates_ci[index_row_2,] = results_uncond.estimates_ci[index_row_1,]
	}
	
	results_uncond.estimates_ci[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4,5,5)'
	results_uncond.estimates_ci[,2] = (0,1,0,1,0,1,0,1,0,1)'
	
	results_uncond.store_infodistributionofp_beta = (info_data.Kbar, info_data.nb_units_studied, 1, c, parameters_hat)

	return(results_uncond)
}
end
