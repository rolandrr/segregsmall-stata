version 14.2
/******************************************************************************/
/* this do-file implements part of Beta method of Rathelot 2012 (R) for the case 
components of the mixture c = 1 (at least for the moment); it does the 
Maximum Likelihood Estimation part
cf. other mata_R_* do-files */
/******************************************************************************/

/******************************************************************************/
/* miscellanous functions */
/******************************************************************************/

/* fonction fun_beta() */
/* older version : return((gamma(a) * gamma(b)) / gamma(a+b))
possible numerical issues for large arguments => . since gamma(too large) = .
whereas the correct value of the beta function is not .
Attempt to correct : use built-in function lngamma() and take exponential of 
the logarithm */
capture mata: mata drop fun_beta()
mata:
mata set matastrict on
real scalar fun_beta(real scalar a, real scalar b) {
	return(exp(lngamma(a) + lngamma(b) - lngamma(a+b)))
}
end

/* fun_gamma_deriv1() */
capture mata: mata drop fun_gamma_deriv1()
mata:
mata set matastrict on
real scalar fun_gamma_deriv1(real scalar x) {
	return(digamma(x)*gamma(x))
}
end

/* fun_gamma_deriv2() */
capture mata: mata drop fun_gamma_deriv2()
mata:
mata set matastrict on
real scalar fun_gamma_deriv2(real scalar x) {
	return(gamma(x) * (trigamma(x) + (digamma(x))^2))
}
end

/* fun_my_beta_deriv1() */
/* Warning: defined for the exponential transformation of tilde argument into
the argument as real positive numbers */
capture mata: mata drop fun_my_beta_deriv1()
mata:
mata set matastrict on
real scalar fun_my_beta_deriv1(	real scalar deriv, ///
								real scalar primitiv, ///
								real scalar constant) {
	real scalar sum_arg
	sum_arg = primitiv + constant
	return(((gamma(constant) * deriv) / ((gamma(sum_arg))^2)) * (fun_gamma_deriv1(primitiv)*gamma(sum_arg) - gamma(primitiv)*fun_gamma_deriv1(sum_arg)))
}
end		

/* trying to improve numerical overflow issue but weird behaviour of optimizer with that..
	real scalar sum_arg, sca_temporary
	sum_arg = primitiv + constant
	sca_temporary = digamma(primitiv) - digamma(sum_arg)
	return(sign(sca_temporary) * exp(lngamma(constant) + ln(deriv) - lngamma(sum_arg) + lngamma(primitiv) + ln(abs(sca_temporary))))
*/					
/* older version
real scalar sum_arg
sum_arg = primitiv + constant
return(((gamma(constant) * deriv) / ((gamma(sum_arg))^2)) * (fun_gamma_deriv1(primitiv)*gamma(sum_arg) - gamma(primitiv)*fun_gamma_deriv1(sum_arg)))
*/

/* fun_my_beta_deriv2() */
/* Warning: defined for the exponential transformation of tilde argument into
the argument as real positive numbers */
capture mata: mata drop fun_my_beta_deriv2()
mata:
mata set matastrict on
real scalar fun_my_beta_deriv2(	real scalar deriv, ///
								real scalar primitiv, ///
								real scalar constant) {						
	real scalar sum_arg
	sum_arg = primitiv + constant						
	return((gamma(constant) * deriv / gamma(sum_arg)^2) * ((fun_gamma_deriv1(primitiv)*gamma(sum_arg) - gamma(primitiv)*fun_gamma_deriv1(sum_arg)) * (1 - 2*deriv*fun_gamma_deriv1(sum_arg)/gamma(sum_arg)) + deriv * (fun_gamma_deriv2(primitiv)*gamma(sum_arg) - gamma(primitiv)*fun_gamma_deriv2(sum_arg))))
}
end

/* trying to improve numerical overflow issue but weird behaviour of optimizer with that..
	real scalar sum_arg, psi_primitiv, psi_sum, sca_temporary
	sum_arg = primitiv + constant
	psi_primitiv = digamma(primitiv)
	psi_sum = digamma(sum_arg)
	sca_temporary = (psi_primitiv - psi_sum)*(1 + deriv*(psi_primitiv - psi_sum)) + deriv*(trigamma(primitiv) - trigamma(sum_arg))
	return(sign(sca_temporary) * exp(lngamma(constant) + lngamma(deriv) + lngamma(primitiv) - lngamma(sum_arg) + ln(abs(sca_temporary))))
*/
/* older version 
before trying to improve numerical overflow issues 
	return((gamma(constant) * deriv / gamma(sum_arg)^2) * ((fun_gamma_deriv1(primitiv)*gamma(sum_arg) - gamma(primitiv)*fun_gamma_deriv1(sum_arg)) * (1 - 2*deriv*fun_gamma_deriv1(sum_arg)/gamma(sum_arg)) + deriv * (fun_gamma_deriv2(primitiv)*gamma(sum_arg) - gamma(primitiv)*fun_gamma_deriv2(sum_arg))))
*/

/* fun_my_beta_deriv12() */
/* Warning: defined for the exponential transformation of tilde argument into
the argument as real positive numbers 
cross-derivative */
capture mata: mata drop fun_my_beta_deriv12()
mata:
mata set matastrict on
real scalar fun_my_beta_deriv12(real scalar deriv_1, ///
								real scalar primitiv_1, ///
								real scalar constant_1, ///
								real scalar deriv_2) {
	real scalar sum_arg
	sum_arg = primitiv_1 + constant_1
	
	real scalar sca_temporary_1,sca_temporary_2 
	sca_temporary_1 = fun_gamma_deriv1(primitiv_1)*gamma(sum_arg) - gamma(primitiv_1)*fun_gamma_deriv1(sum_arg)
	sca_temporary_2 = fun_gamma_deriv1(sum_arg)*fun_gamma_deriv1(primitiv_1) - gamma(primitiv_1)*fun_gamma_deriv2(sum_arg)
	
	return((deriv_1 * deriv_2 / gamma(sum_arg)^2) * ((sca_temporary_1 * (fun_gamma_deriv1(constant_1) - 2 * fun_gamma_deriv1(sum_arg) * gamma(constant_1) / gamma(sum_arg))) + sca_temporary_2 * gamma(constant_1)))
}
end								

/* Def_matrix_betas() */
/* matrix_betas : for size of unit k = 1..Kbar (k = m in Rathelot) that are the rows
	and number of minority individual x = 0..k (x = n in Rathelot) that are the columns, for x > k, entry is set to 0
	for given parameters alpha and beta
	entry is the latest fraction of Rathelot 2012 equation (5) */
capture mata: mata drop Def_matrix_betas()
mata:
mata set matastrict on
real matrix Def_matrix_betas(real scalar Kbar, real scalar alpha, real scalar beta, real scalar lambda) {

	real matrix matrix_betas
	real scalar k, x
	
	matrix_betas = J(Kbar, Kbar+1, 0)
	
	for(k = 1; k <= Kbar; k++) {
		for(x = 0; x <= Kbar; x++) {
			if (x <= k) {
				matrix_betas[k, x+1] = lambda * (fun_beta(alpha + x, beta + k - x) / fun_beta(alpha, beta))
			}
		}
	}
	return(matrix_betas)
}
end

/* my_prod */
/* NB : only for positive numbers (for ln), ok here in the application */
capture mata: mata drop my_prod()
mata:
mata set matastrict on
real scalar my_prod(real vector my_vector) {
	return(exp(quadsum(ln(my_vector), 0)))
}
end

capture mata: mata drop Get_Exp_K_from_db()
mata:
mata set matastrict on
real scalar Get_Exp_K_from_db(real matrix db) {
	return((db[.,1]' * db[.,2]) / sum(db[.,2]))
}
end

capture mata: mata drop Get_Exp_X_from_db()
mata:
mata set matastrict on
real scalar Get_Exp_X_from_db(real matrix db, real scalar p) {
	real scalar i, j
	real matrix db_temp
	db_temp = db
	for (i = 1; i <= rows(db); i++) {
		for(j = 5; j <= cols(db); j++) {
			if (j >= (i+4)) {
				db_temp[i,j] = .
			}
		}
	}
	return((quadcolsum(db_temp[.,3..cols(db)], 0) * ((0..rows(db)):^p)') / sum(db[.,2]))
}
end
/* db_temp utilise car Mata semble modifie en dehors de l'environnement local
de la fonction la valeur de db : Attention !! */


/******************************************************************************/
/* eval_parametric_ML() */
/* evaluator for optimization of ML */
/******************************************************************************/
capture mata: mata drop R
capture mata: mata drop eval_parametric_ML()
mata:
mata set matastrict on
void eval_parametric_ML(real scalar todo, ///
						real rowvector parameters_tilde, ///
						real matrix db, ///
						val, grad, hess) {
	
	real scalar alpha, beta, D
	real scalar Kbar, c
	real scalar k, x, j, j1, j2, sca_temporary_1, sca_temporary_2, sca_temporary_3
	real matrix matrix_betas, mat_temporary, mat_A
	real rowvector alphas, betas, lambdas_tilde, lambdas
	real matrix matrix_derivatives_arg1, matrix_derivatives_arg2, matrix_derivatives_arg11, matrix_derivatives_arg12, matrix_derivatives_arg22
	
	Kbar = rows(db)
	
	/* case : 1-component mixture beta = beta distribution simply */
	if (length(parameters_tilde) == 2) {
	
		/* Computation of value - case c = 1 */
		alpha = exp(parameters_tilde[1])
		beta = exp(parameters_tilde[2])
		/* matrix_betas : for size of unit k = 1..Kbar (k = m in Rathelot) that are the rows
		and number of minority individual x = 0..k (x = n in Rathelot) that are the columns, for x > k, entry is set to 0
		for given parameters alpha and beta
		entry is the latest fraction of Rathelot 2012 equation (5) */
		matrix_betas = J(Kbar, Kbar+1, 0)
		for(k = 1; k <= Kbar; k++) {
			for(x = 0; x <= Kbar; x++) {
				if (x <= k) {
					matrix_betas[k, x+1] = fun_beta(alpha + x, beta + k - x) / fun_beta(alpha, beta)
				}
			}
		}
		
		val = quadsum(db[.,3..cols(db)] :* log(matrix_betas :+ (matrix_betas :== 0)), 1)
	
		/* Computation of gradient - case c = 1 */
		if (todo >= 1) {
		
			matrix_derivatives_arg1 = J(Kbar, Kbar+1, 0)
			matrix_derivatives_arg2 = J(Kbar, Kbar+1, 0)
			for(k = 1; k <= Kbar; k++) {
				for(x = 0; x <= Kbar; x++) {
					if (x <= k) {
						matrix_derivatives_arg1[k, x+1] = alpha * (digamma(alpha+x) - digamma(alpha+beta+k) - digamma(alpha) + digamma(alpha+beta))
						matrix_derivatives_arg2[k, x+1] = beta * (digamma(beta+k-x) - digamma(alpha+beta+k) - digamma(beta) + digamma(alpha+beta))
					}
				}
			}
		
			grad[1] = quadsum(db[.,3..cols(db)] :* matrix_derivatives_arg1, 1)
			grad[2] = quadsum(db[.,3..cols(db)] :* matrix_derivatives_arg2, 1)
			
			/* Computation of hessian - case c = 1 */
			if (todo == 2) {
			
				matrix_derivatives_arg11 = J(Kbar, Kbar+1, 0)
				matrix_derivatives_arg22 = J(Kbar, Kbar+1, 0)
				matrix_derivatives_arg12 = J(Kbar, Kbar+1, 0)
				for(k = 1; k <= Kbar; k++) {
					for(x = 0; x <= Kbar; x++) {
						if (x <= k) {
						
							matrix_derivatives_arg11[k, x+1] = matrix_derivatives_arg1[k, x+1] + (alpha^2) * (trigamma(alpha+x) + trigamma(alpha+beta) - trigamma(alpha+beta+k) - trigamma(alpha))
							matrix_derivatives_arg22[k, x+1] = matrix_derivatives_arg2[k, x+1] + (beta^2) * (trigamma(beta+k-x) + trigamma(alpha+beta) - trigamma(alpha+beta+k) - trigamma(beta))
							matrix_derivatives_arg12[k, x+1] = alpha * beta * (trigamma(alpha+beta) - trigamma(alpha+beta+k))
														
						}
					}
				}
				
				hess[1,1] = quadsum(db[.,3..cols(db)] :* matrix_derivatives_arg11, 1)
				hess[2,2] = quadsum(db[.,3..cols(db)] :* matrix_derivatives_arg22, 1)
				hess[2,1] = quadsum(db[.,3..cols(db)] :* matrix_derivatives_arg12, 1)
				_makesymmetric(hess)
			}
			
		}
	}
	
	/* case : c-component mixture beta, c > 1 */
	if (length(parameters_tilde) >= 3) {
	
		/* Computation of value - case c > 1 */
		c = (length(parameters_tilde)+1)/3
		alphas = exp(parameters_tilde[1..c])
		betas = exp(parameters_tilde[(c+1)..(2*c)])
		lambdas_tilde = parameters_tilde[(2*c+1)..(3*c-1)]
		D = 1 + quadsum(fun_logistic(lambdas_tilde), 1) /* scalar D = 1 + sum_{j=1}^{c-1} Logistic(lambda_tilde_j) */
		lambdas = (fun_logistic(lambdas_tilde), 1) :/ D
			
		matrix_betas = J(Kbar, Kbar+1, 0)
		for(j = 1; j <= c; j++) {
			matrix_betas = matrix_betas :+ Def_matrix_betas(Kbar, alphas[j],  betas[j], lambdas[j])
		}
		
		val = quadsum(db[.,3..cols(db)] :* log(matrix_betas :+ (matrix_betas :== 0)), 1)
	
		/* Computation of gradient - case c > 1 */
		if (todo >= 1) {
		
			/* derivative w.r.t. alpha_j, j = 1..c */
			for(j = 1; j <= c; j++) {
		
				mat_temporary = J(Kbar, Kbar+1, 0)
				for(k = 1; k <= Kbar; k++) {
					for(x = 0; x <= Kbar; x++) {
						if (x <= k) {
							mat_temporary[k, x+1] = (lambdas[j] / (fun_beta(alphas[j], betas[j])^2)) * ((fun_my_beta_deriv1(alphas[j], alphas[j]+x, betas[j]+k-x)*fun_beta(alphas[j], betas[j])) - (fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(alphas[j], alphas[j], betas[j]))) / matrix_betas[k, x+1]
						}
					}
				}
				grad[j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
			}
			
			/* derivative w.r.t. beta_j, j = 1..c */
			for(j = 1; j <= c; j++) {
		
				mat_temporary = J(Kbar, Kbar+1, 0)
				for(k = 1; k <= Kbar; k++) {
					for(x = 0; x <= Kbar; x++) {
						if (x <= k) {
							mat_temporary[k, x+1] = (lambdas[j] / (fun_beta(alphas[j], betas[j])^2)) * ((fun_my_beta_deriv1(betas[j], betas[j]+k-x, alphas[j]+x)*fun_beta(alphas[j], betas[j])) - (fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(betas[j], betas[j], alphas[j]))) / matrix_betas[k, x+1]
						}
					}
				}
				grad[c+j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
			}
			
			/* derivative w.r.t. lambda_j, j = 1..c-1 */
			for(j = 1; j <= c-1; j++) {
			
				mat_temporary = J(Kbar, Kbar+1, 0)
				for(k = 1; k <= Kbar; k++) {
					for(x = 0; x <= Kbar; x++) {
						if (x <= k) {
							mat_temporary[k, x+1] = (fun_logistic_deriv1(lambdas_tilde[j]) / D) * ((fun_beta(alphas[j]+x, betas[j]+k-x) / fun_beta(alphas[j], betas[j])) - matrix_betas[k, x+1]) / matrix_betas[k, x+1]	
						}
					}
				}
				grad[(2*c)+j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
			}
			
			/* Computation of hessian - case c > 1 */
			if (todo == 2) {
			
				/* hessian alpha_j1 alpha_j2 j1 != j2, j1, j2 = 1..c */
				for(j1 = 2; j1 <= c; j1++) {
					for(j2 = 1; j2 < j1; j2++) {	
						mat_temporary = J(Kbar, Kbar+1, 0)
						for(k = 1; k <= Kbar; k++) {
							for(x = 0; x <= Kbar; x++) {
								if (x <= k) {
									mat_temporary[k, x+1] = (lambdas[j1] / (fun_beta(alphas[j1], betas[j1])^2)) * ((fun_my_beta_deriv1(alphas[j1], alphas[j1]+x, betas[j1]+k-x)*fun_beta(alphas[j1], betas[j1])) - (fun_beta(alphas[j1]+x, betas[j1]+k-x)*fun_my_beta_deriv1(alphas[j1], alphas[j1], betas[j1]))) * (-lambdas[j2] / (matrix_betas[k, x+1]^2 * fun_beta(alphas[j2], betas[j2])^2)) * (fun_my_beta_deriv1(alphas[j2], alphas[j2]+x, betas[j2]+k-x)*fun_beta(alphas[j2], betas[j2]) - fun_beta(alphas[j2]+x, betas[j2]+k-x)*fun_my_beta_deriv1(alphas[j2], alphas[j2], betas[j2]))					
								}
							}
						}
						hess[j1,j2] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
					}
				}
				
				/* hessian alpha_j1 beta_j2 j1 != j2, j1, j2 = 1..c */
				for(j1 = 1; j1 <= c; j1++) {
					for(j2 = 1; j2 <= c; j2++) {	
						if (j1 != j2) {
							mat_temporary = J(Kbar, Kbar+1, 0)
							for(k = 1; k <= Kbar; k++) {
								for(x = 0; x <= Kbar; x++) {
									if (x <= k) {
										mat_temporary[k, x+1] = (lambdas[j1] / (fun_beta(alphas[j1], betas[j1])^2)) * ((fun_my_beta_deriv1(alphas[j1], alphas[j1]+x, betas[j1]+k-x)*fun_beta(alphas[j1], betas[j1])) - (fun_beta(alphas[j1]+x, betas[j1]+k-x)*fun_my_beta_deriv1(alphas[j1], alphas[j1], betas[j1]))) * (-lambdas[j2] / (matrix_betas[k, x+1]^2 * fun_beta(alphas[j2], betas[j2])^2)) * (fun_my_beta_deriv1(betas[j2], betas[j2]+k-x, alphas[j2]+x)*fun_beta(alphas[j2], betas[j2]) - fun_beta(alphas[j2]+x, betas[j2]+k-x)*fun_my_beta_deriv1(betas[j2], betas[j2], alphas[j2]))					
									}
								}
							}
							hess[c+j2,j1] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
						}
					}
				}
				
				/* hessian alpha_j1 alpha_j2 j1 = j2 = j, j = 1..c */
				for(j = 1; j <= c; j++) {
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {
								sca_temporary_1 = (fun_my_beta_deriv1(alphas[j], alphas[j]+x, betas[j]+k-x)*fun_beta(alphas[j], betas[j])) - (fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(alphas[j], alphas[j], betas[j]))
								sca_temporary_2 = (fun_my_beta_deriv2(alphas[j], alphas[j]+x, betas[j]+k-x)*fun_beta(alphas[j], betas[j])) - (fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv2(alphas[j], alphas[j], betas[j]))
								mat_temporary[k, x+1] = (lambdas[j] / (fun_beta(alphas[j], betas[j])^2 * matrix_betas[k, x+1])) * (sca_temporary_2 - (2 * sca_temporary_1 * fun_my_beta_deriv1(alphas[j], alphas[j], betas[j]) / fun_beta(alphas[j], betas[j])) - (lambdas[j] * sca_temporary_1^2 / (fun_beta(alphas[j], betas[j])^2 * matrix_betas[k, x+1])))
							}
						}
					}
					hess[j,j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
				}
				
				/* hessian alpha_j1 beta_j2 j1 = j2 = j, j = 1..c */
				for(j = 1; j <= c; j++) {
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {						
								sca_temporary_1 = (fun_my_beta_deriv1(alphas[j], alphas[j]+x, betas[j]+k-x)*fun_beta(alphas[j], betas[j])) - (fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(alphas[j], alphas[j], betas[j]))
								sca_temporary_2 = fun_my_beta_deriv12(alphas[j], alphas[j]+x,  betas[j]+k-x, betas[j])*fun_beta(alphas[j], betas[j]) + fun_my_beta_deriv1(alphas[j], alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(betas[j], betas[j], alphas[j]) - fun_my_beta_deriv1(betas[j], betas[j]+k-x, alphas[j]+x)*fun_my_beta_deriv1(alphas[j], alphas[j], betas[j]) - fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv12(alphas[j], alphas[j], betas[j], betas[j])
								sca_temporary_3 = fun_my_beta_deriv1(betas[j], betas[j]+k-x, alphas[j]+x)*fun_beta(alphas[j], betas[j]) - fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(betas[j], betas[j], alphas[j])
								mat_temporary[k, x+1] = (lambdas[j] / (fun_beta(alphas[j], betas[j])^2 * matrix_betas[k, x+1])) * (sca_temporary_2 - (2 * sca_temporary_1 * fun_my_beta_deriv1(betas[j], betas[j], alphas[j]) / fun_beta(alphas[j], betas[j])) - (lambdas[j] * sca_temporary_1 * sca_temporary_3 / (fun_beta(alphas[j], betas[j])^2 * matrix_betas[k, x+1])))
							}
						}
					}
					hess[c+j,j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
				}
				
				/* hessian beta_j1 beta_j2, j1 != j2 j1, j2 = 1..c */
				for(j1 = 2; j1 <= c; j1++) {
					for(j2 = 1; j2 < j1; j2++) {
						mat_temporary = J(Kbar, Kbar+1, 0)
						for(k = 1; k <= Kbar; k++) {
							for(x = 0; x <= Kbar; x++) {
								if (x <= k) {
									sca_temporary_1 = (fun_my_beta_deriv1(betas[j1], betas[j1]+k-x, alphas[j1]+x)*fun_beta(alphas[j1], betas[j1]) - fun_beta(alphas[j1]+x, betas[j1]+k-x)*fun_my_beta_deriv1(betas[j1], betas[j1], alphas[j1])) / (fun_beta(alphas[j1], betas[j1])^2)			
									sca_temporary_2 = (fun_my_beta_deriv1(betas[j2], betas[j2]+k-x, alphas[j2]+x)*fun_beta(alphas[j2], betas[j2]) - fun_beta(alphas[j2]+x, betas[j2]+k-x)*fun_my_beta_deriv1(betas[j2], betas[j2], alphas[j2])) / (fun_beta(alphas[j2], betas[j2])^2)
									mat_temporary[k, x+1] = lambdas[j1] * sca_temporary_1 * (-lambdas[j2]) * sca_temporary_2 / (matrix_betas[k, x+1]^2)
								}
							}
						}
						hess[c+j1,c+j2] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
					}			
				}
				
				/* hessian beta_j1 beta_j2, j1 = j2 = j, j = 1..c */
				for(j = 1; j <= c; j++) {
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {							
								sca_temporary_1 = fun_my_beta_deriv1(betas[j], betas[j]+k-x, alphas[j]+x)*fun_beta(alphas[j], betas[j]) - fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(betas[j], betas[j], alphas[j])								
								sca_temporary_2 = fun_my_beta_deriv2(betas[j], betas[j]+k-x, alphas[j]+x)*fun_beta(alphas[j], betas[j]) - fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv2(betas[j], betas[j], alphas[j])							
								mat_temporary[k, x+1] = (lambdas[j] / (fun_beta(alphas[j], betas[j])^2 * matrix_betas[k, x+1])) * (sca_temporary_2 - (2 * sca_temporary_1 * fun_my_beta_deriv1(betas[j], betas[j], alphas[j]) / fun_beta(alphas[j], betas[j])) - (lambdas[j] * sca_temporary_1^2 / (fun_beta(alphas[j], betas[j])^2 * matrix_betas[k, x+1])))							
							}
						}
					}
					hess[c+j,c+j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)
				}	
				
				
				mat_A = D :* matrix_betas
				
				/* hessian alpha_j1 lambda_j2, j1 = j2 = j, j = 1..c-1 */
				for(j = 1; j <= c-1; j++) {
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {	
								sca_temporary_1 = fun_my_beta_deriv1(alphas[j], alphas[j]+x, betas[j]+k-x)*fun_beta(alphas[j], betas[j]) - fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(alphas[j], alphas[j], betas[j])
								mat_temporary[k, x+1] = (sca_temporary_1 * fun_logistic_deriv1(lambdas_tilde[j]) / (mat_A[k, x+1] * fun_beta(alphas[j], betas[j])^2)) * (1 - (fun_logistic(lambdas_tilde[j]) * fun_beta(alphas[j]+x, betas[j]+k-x) / (mat_A[k, x+1] * fun_beta(alphas[j], betas[j]))))
							}
						}
					}
					hess[2*c+j,j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
				}
				
				/* hessian beta_j1 lambda_j2, j1 = j2 = j, j = 1..c-1 */
				for(j = 1; j <= c-1; j++) {
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {	
								sca_temporary_1 = fun_my_beta_deriv1(betas[j], betas[j]+k-x, alphas[j]+x)*fun_beta(alphas[j], betas[j]) - fun_beta(alphas[j]+x, betas[j]+k-x)*fun_my_beta_deriv1(betas[j], betas[j], alphas[j])
								mat_temporary[k, x+1] = (sca_temporary_1 * fun_logistic_deriv1(lambdas_tilde[j]) / (mat_A[k, x+1] * fun_beta(alphas[j], betas[j])^2)) * (1 - (fun_logistic(lambdas_tilde[j]) * fun_beta(alphas[j]+x, betas[j]+k-x) / (mat_A[k, x+1] * fun_beta(alphas[j], betas[j]))))
							}
						}
					}
					hess[2*c+j,c+j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
				}
				
				/* hessian alpha_j1 lambda_j2, j1 != j2, j1 = 1..c (+ special case j1 = c), j2 = 1..c-1 */
				for(j1 = 1; j1 <= c-1; j1++) {
					for(j2 = 1; j2 <= c-1; j2++) {
						if(j1 != j2) {
							mat_temporary = J(Kbar, Kbar+1, 0)
							for(k = 1; k <= Kbar; k++) {
								for(x = 0; x <= Kbar; x++) {
									if (x <= k) {	
										sca_temporary_1 = fun_my_beta_deriv1(alphas[j1], alphas[j1]+x, betas[j1]+k-x)*fun_beta(alphas[j1], betas[j1]) - fun_beta(alphas[j1]+x, betas[j1]+k-x)*fun_my_beta_deriv1(alphas[j1], alphas[j1], betas[j1])
										mat_temporary[k, x+1] = (sca_temporary_1 * fun_logistic(lambdas_tilde[j1]) / (fun_beta(alphas[j1], betas[j1])^2)) * (-fun_logistic_deriv1(lambdas_tilde[j2])) * fun_beta(alphas[j2]+x, betas[j2]+k-x) / (fun_beta(alphas[j2], betas[j2]) * mat_A[k, x+1]^2)
									}
								}
							}
							hess[2*c+j2, j1] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
						}
					}
				}
				for(j2 = 1; j2 <= c-1; j2++) {
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {	
								sca_temporary_1 = fun_my_beta_deriv1(alphas[j1], alphas[j1]+x, betas[j1]+k-x)*fun_beta(alphas[j1], betas[j1]) - fun_beta(alphas[j1]+x, betas[j1]+k-x)*fun_my_beta_deriv1(alphas[j1], alphas[j1], betas[j1])
								mat_temporary[k, x+1] = (sca_temporary_1 / (fun_beta(alphas[j1], betas[j1])^2)) * (-fun_logistic_deriv1(lambdas_tilde[j2])) * fun_beta(alphas[j2]+x, betas[j2]+k-x) / (fun_beta(alphas[j2], betas[j2]) * mat_A[k, x+1]^2)
							}
						}
					}
					hess[2*c+j2, c] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
				}
				
				/* hessian beta_j1 lambda_j2, j1 != j2, j1 = 1..c (+ special case j1 = c), j2 = 1..c-1 */
				for(j1 = 1; j1 <= c-1; j1++) {
					for(j2 = 1; j2 <= c-1; j2++) {
						if(j1 != j2) {
							mat_temporary = J(Kbar, Kbar+1, 0)
							for(k = 1; k <= Kbar; k++) {
								for(x = 0; x <= Kbar; x++) {
									if (x <= k) {	
										sca_temporary_1 = fun_my_beta_deriv1(betas[j1], betas[j1]+k-x, alphas[j1]+x)*fun_beta(alphas[j1], betas[j1]) - fun_beta(alphas[j1]+x, betas[j1]+k-x)*fun_my_beta_deriv1(betas[j1], betas[j1], alphas[j1])
										mat_temporary[k, x+1] = (sca_temporary_1 * fun_logistic(lambdas_tilde[j1]) / (fun_beta(alphas[j1], betas[j1])^2)) * (-fun_logistic_deriv1(lambdas_tilde[j2])) * fun_beta(alphas[j2]+x, betas[j2]+k-x) / (fun_beta(alphas[j2], betas[j2]) * mat_A[k, x+1]^2)
									}
								}
							}
							hess[2*c+j2, c+j1] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
						}
					}
				}
				for(j2 = 1; j2 <= c-1; j2++) {
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {	
								sca_temporary_1 = fun_my_beta_deriv1(betas[j1], betas[j1]+k-x, alphas[j1]+x)*fun_beta(alphas[j1], betas[j1]) - fun_beta(alphas[j1]+x, betas[j1]+k-x)*fun_my_beta_deriv1(betas[j1], betas[j1], alphas[j1])
								mat_temporary[k, x+1] = (sca_temporary_1 / (fun_beta(alphas[j1], betas[j1])^2)) * (-fun_logistic_deriv1(lambdas_tilde[j2])) * fun_beta(alphas[j2]+x, betas[j2]+k-x) / (fun_beta(alphas[j2], betas[j2]) * mat_A[k, x+1]^2)
							}
						}
					}
					hess[2*c+j2, 2*c] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
				}
				
				/* hessian lambda_j1 lambda_j2, j1 != j2, j1, j2 = 1..c-1*/
				for(j1 = 2; j1 <= c-1; j1++) {
					for(j2 = 1; j2 <= c-1; j2++) {
						if(j1 != j2) {
							mat_temporary = J(Kbar, Kbar+1, 0)
							for(k = 1; k <= Kbar; k++) {
								for(x = 0; x <= Kbar; x++) {
									if (x <= k) {	
										sca_temporary_1 = fun_beta(alphas[j1]+x, betas[j1]+k-x) / fun_beta(alphas[j1], betas[j1])
										sca_temporary_2 = fun_beta(alphas[j2]+x, betas[j2]+k-x) / fun_beta(alphas[j2], betas[j2])
										mat_temporary[k, x+1] = fun_logistic_deriv1(lambdas_tilde[j1]) * fun_logistic_deriv1(lambdas_tilde[j2]) * ((sca_temporary_1 * sca_temporary_2 / (mat_A[k, x+1]^2)) + (1/(D^2)))										
									}
								}
							}
							hess[2*c+j1, 2*c+j2] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
						}
					}
				}
				
				/* hessian lambda_j1 lambda_j2, j1 = j2 = j, j = 1..c-1*/
				for(j = 1; j <= c-1; j++) {
				
					mat_temporary = J(Kbar, Kbar+1, 0)
					for(k = 1; k <= Kbar; k++) {
						for(x = 0; x <= Kbar; x++) {
							if (x <= k) {
								sca_temporary_1 = fun_beta(alphas[j]+x, betas[j]+k-x) / fun_beta(alphas[j], betas[j])
								mat_temporary[k, x+1] = fun_logistic_deriv2(lambdas_tilde[j]) * ((sca_temporary_1 / mat_A[k, x+1]) - (1/D)) + (fun_logistic_deriv1(lambdas_tilde[j])^2) * ((1/D^2) - ((sca_temporary_1 / mat_A[k, x+1])^2))								
							}
						}
					}
					hess[2*c+j, 2*c+j] = quadsum(db[.,3..cols(db)] :* mat_temporary, 1)	
				}
				
				_makesymmetric(hess)
				
			
			} /* end todo == 2 */
		} /* end todo >= 1 */
	} /* length(parameters_tilde) >= 3 */
}
end


/******************************************************************************/
/* Method of Moments (MM) */
/* Used for defining initial values in the ML optimization */
/******************************************************************************/

capture mata: mata drop MM_beta_estimate_oneK()
mata:
mata set matastrict on
real rowvector MM_beta_estimate_oneK(real matrix db, real scalar c) {

	real scalar K, nb_obs, m1, m2, denominator
	real rowvector parameters_hat, db_oneK

	db_oneK = Extract_data_Nks_from_db(db, rows(db))
	
	/* case : c = 1 */
	if (c == 1) {

		/* source https://en.wikipedia.org/wiki/Beta-binomial_distribution#Method_of_moments */
		/* Problem if denominator = 0, return missing values */
		K = db_oneK[1]
		nb_obs = db_oneK[2]
		m1 = (db_oneK[3..K+3] * (0..K)') / nb_obs
		m2 = (db_oneK[3..K+3] * ((0..K):^2)') / nb_obs
		parameters_hat = J(1, 2, .)
		denominator = K*((m2/m1)-m1-1) + m1
		parameters_hat[1] = (K*m1 - m2) / denominator	
		parameters_hat[2] = ((K-m1)*(K - (m2/m1))) / denominator
	}
	
	/* case : c > 1 */
	else {
		errprintf("Error: c > 1")
		exit()
	}
	
	return(parameters_hat)
}
end

capture mata: mata drop MM_beta_estimate_allK()
mata:
mata set matastrict on
real rowvector MM_beta_estimate_allK(real matrix db, real scalar c) {

	real scalar EK, m1, m2, denominator
	real rowvector parameters_hat

	/* case : c = 1 */
	if (c == 1) {
	
		/* source https://en.wikipedia.org/wiki/Beta-binomial_distribution#Method_of_moments */
		/* Problem if denominator = 0, return missing values */
		EK = Get_Exp_K_from_db(db)
		m1 = Get_Exp_X_from_db(db, 1)
		m2 = Get_Exp_X_from_db(db, 2)
		parameters_hat = J(1, 2, .)
		denominator = EK*((m2/m1)-m1-1) + m1
		parameters_hat[1] = (EK*m1 - m2) / denominator	
		parameters_hat[2] = ((EK-m1)*(EK - (m2/m1))) / denominator
	}
	
	/* case : c > 1 */
	else {
		errprintf("Error: c > 1")
		exit()
	}
	
	return(parameters_hat)
}
end


/******************************************************************************/
/* Transform_beta_parameters */
/* from the tilde transform of the parameter obtained as output of optimize,
get the transformed one */
/******************************************************************************/
capture mata: mata drop Transform_beta_parameters()
mata:
mata set matastrict on
real rowvector Transform_beta_parameters(real rowvector parameters_tilde) {

	real rowvector parameters, lambdas_tilde
	real scalar c, D
	
	/* case : c = 1 */
	if (length(parameters_tilde) == 2) {
		parameters = exp(parameters_tilde)
	}
	
	/* case : c > 1 */
	else {
		c = (length(parameters_tilde)+1)/3
		parameters = J(1, 3*c, .)
		parameters[1..(2*c)] = exp(parameters_tilde[1..(2*c)])
		lambdas_tilde = parameters_tilde[(2*c+1)..(3*c-1)]
		D = 1 + quadsum(fun_logistic(lambdas_tilde), 1) /* scalar D = 1 + sum_{j=1}^{c-1} Logistic(lambda_tilde_j) */
		parameters[(2*c+1)..(3*c)] = (fun_logistic(lambdas_tilde), 1) :/ D
	}
	
	return(parameters)
}
end


/******************************************************************************/
/* beta_ML() */
/* get the estimated parameters and estimated variance
NB: the first row of the returned matrix is the estimated parameters
NB : the following lines are the matrix of estimated variance */
/******************************************************************************/
capture mata: mata drop beta_ML()
mata:
mata set matastrict on
real matrix beta_ML(	real matrix db, ///
						real scalar nb_K_positive_obs, ///
						real scalar c, ///
						struct options_optimization_CML optionsoptCML) {
	
	transmorphic R
	real rowvector parameters_tilde, parameters_MM_hat
	real scalar code_optimize
	real rowvector parameters_hat
	real matrix variance_hat
	real scalar epsilon_ln
	epsilon_ln = 10^(-7)
	real scalar lower_bound, upper_bound
	
	/* ML estimation of the parameters */

	if ((rows(db) == 1)&(db[1,1] == 1)&(db[1,3] == db[1,4])) { /* very particular case where ML impossible (begin) */	
		/* particular case where EMV is impossible (cf. details in manuscript complement R 2012 paper 
		return the estimation of a Uniform([0,1]) i.e. alpha = beta = 1	*/
		parameters_hat = (1,1)
		variance_hat = (1,0)\(0,1)
	} /* very particular case where ML impossible (end) */
	
	else { /* general case (begin) */
	
		if (c == 1) {
		
			if (nb_K_positive_obs > 1) 	parameters_MM_hat = MM_beta_estimate_allK(db, 1)
			else parameters_MM_hat = MM_beta_estimate_oneK(db, 1)
					
			if (sum(parameters_MM_hat :== .)) { /* in case missing values in parameters_MM_hat (begin) */
					parameters_tilde = (0,0)
			} /* in case missing values in parameters_MM_hat (end) */		
			else {			
				/* in case the MM estimates yields negative values => to avoid missing 
				values for the optimisation starting values : just a little smaller than 0 for the ln 
				=> 0 for the exponential */
				if (sum((parameters_MM_hat :<= epsilon_ln)) > 0) {
					parameters_MM_hat[selectindex(parameters_MM_hat :<= epsilon_ln)] = J(1, length(selectindex(parameters_MM_hat :<= epsilon_ln)), epsilon_ln)
				}
				parameters_tilde = ln(parameters_MM_hat)
			}
						
			R = optimize_init()
			optimize_init_which(R, "max")
			optimize_init_evaluator(R, &eval_parametric_ML())
			optimize_init_evaluatortype(R, "d2")
			optimize_init_technique(R, "nr")
			optimize_init_conv_maxiter(R, optionsoptCML.nb_max_iter)
			optimize_init_params(R, parameters_tilde)
			optimize_init_argument(R, 1, db)
			optimize_init_conv_ptol(R, optionsoptCML.ptol)
			optimize_init_conv_vtol(R, optionsoptCML.vtol)
			optimize_init_conv_nrtol(R, optionsoptCML.nrtol)
			optimize_init_verbose(R, 0)
			optimize_init_tracelevel(R, "none")
			optimize_init_conv_warning(R, "off")

			code_optimize = _optimize(R)
			if (code_optimize == 0){
				parameters_hat = Transform_beta_parameters(optimize_result_params(R))
				variance_hat = optimize_result_V_oim(R) 
				/* Remark: this is variance of the parameters_tilde before transformation, it is taken into account in Delta-Method OK */
			}
			else {
				errprintf("Error in optimization - function beta_ML() \n")
				exit()
			}
		}
		else {
			errprintf("Error c > 1 - not done yet beta_ML()")
			exit()
		}
		
		/* control for extreme values of parameters_hat, set a minimum and a maximum
		to avoid numerical under or overflow when computing the segregation indices */
		lower_bound = 10^(-4)
		upper_bound = 10^(4)
		if (sum((parameters_hat :< lower_bound))){
			parameters_hat[selectindex(parameters_hat :< lower_bound)] = J(1, length(selectindex(parameters_hat :< lower_bound)), lower_bound)
		}
		if (sum((parameters_hat :> upper_bound))){
			parameters_hat[selectindex(parameters_hat :> upper_bound)] = J(1, length(selectindex(parameters_hat :> upper_bound)), upper_bound)
		}
		
	} /* general case (end) */
	
	return((parameters_hat \ variance_hat))
}
end
