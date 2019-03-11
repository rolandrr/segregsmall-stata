/* This do-files express the bounds for the case K and p independent
In this setting, cf. Appendix, the bounds are estimated as in the cae with
a single unit size, hence no difference between unweighted (unit level)
or weighted (individual level) indices
Hence, directly the use of these function for computing the bounds 
when K and p are independent 
instead of trying to reuse the function and structure of the random K case */

/******************************************************************************/
/* Use some functions from CML K random / K by K  and adapt other for Kpind */
/******************************************************************************/

capture mata: mata drop Def_int_h_Fpk_Duncan_Kpind()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Duncan_Kpind(struct raw_ML_results_Kpind scalar raw_ML, ///
									real scalar m, real scalar nb_candidates_cr){
		
	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	
	if (raw_ML.init_esti.nb_obs) { /* case : with observations */
		if (raw_ML.I_constrained_distribution) { /* case : constrained distribution */
			int_h_Fpk[1] = compute_int_h_Fp_dis_Duncan(raw_ML.CML.xy_hat, m)
			int_h_Fpk[2] = int_h_Fpk[1]
		}
		else { /* case : unconstrained distribution */
			if(raw_ML.init_esti.Kbar == 1) { /* case : K = 1 */
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

capture mata: mata drop Def_int_h_Fpk_Theil_Kpind()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Theil_Kpind(struct raw_ML_results_Kpind scalar raw_ML, real scalar m){

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
			if(raw_ML.init_esti.Kbar == 1) { /* case : K = 1 */
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

capture mata: mata drop Def_int_h_Fpk_Atkinson_Kpind()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Atkinson_Kpind(struct raw_ML_results_Kpind scalar raw_ML, real scalar m, real scalar b){

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
			if(raw_ML.init_esti.Kbar == 1) { /* case : K = 1 */
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
end	

capture mata: mata drop Def_int_h_Fpk_Coworker_Kpind()
mata:
mata set matastrict on
real rowvector Def_int_h_Fpk_Coworker_Kpind(struct raw_ML_results_Kpind scalar raw_ML, real scalar m){

	real rowvector int_h_Fpk 
	int_h_Fpk = J(1, 2, .)
	
	if (raw_ML.init_esti.nb_obs) { /* case : with observations */	
	
		if (raw_ML.I_constrained_distribution) { /* case : constrained distribution */
			int_h_Fpk[1] = compute_int_h_Fp_dis_Coworker(raw_ML.CML.xy_hat, m)
			int_h_Fpk[2] = int_h_Fpk[1]
		}
		else { /* case : unconstrained distribution */
			if(raw_ML.init_esti.Kbar == 1) { /* case : K = 1 */
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
/* Bounds_DTACW_Kpind() */
/******************************************************************************/
/* in the setting K and p independent, unit and individual level
weights indices coincide ; nonetheless, we return the same output
cf. matrix of Bounds_DTACW() for ease of use in the Stata output
but only has to be computed once, and no just one K here 
so can simplify the different sums over K of the case with K random */
capture mata: mata drop Bounds_DTACW_Kpind()
mata:
mata set matastrict on
real matrix Bounds_DTACW_Kpind(real scalar b_atkinson, ///
							struct raw_ML_results_Kpind scalar raw_ML_Kpind, ///
							struct options_optimization scalar options_opt){

	real matrix bounds_DTACW
	real rowvector int_h_Fpk
	real scalar m01_hat
	
	bounds_DTACW = J(8, 4, .)
	bounds_DTACW[,1] = (1,1,2,2,3+b_atkinson,3+b_atkinson,4,4)'
	bounds_DTACW[,2] = (0,1,0,1,0,1,0,1)'
	
	m01_hat = (raw_ML_Kpind.I_constrained_distribution ? raw_ML_Kpind.CML.m1_hat : raw_ML_Kpind.UML.m_tilde[1])
	
	/* Duncan */
	int_h_Fpk = Def_int_h_Fpk_Duncan_Kpind(raw_ML_Kpind, m01_hat, options_opt.nb_candidates_cr)
	bounds_DTACW[1,3] = compute_nu_Duncan(int_h_Fpk[1], m01_hat)
	bounds_DTACW[1,4] = compute_nu_Duncan(int_h_Fpk[2], m01_hat)
	bounds_DTACW[2,(3,4)] = bounds_DTACW[1,(3,4)]
	
	/* Theil */
	int_h_Fpk = Def_int_h_Fpk_Theil(raw_ML_Kpind, m01_hat)
	bounds_DTACW[3,3] = compute_nu_Theil(int_h_Fpk[1], m01_hat)
	bounds_DTACW[3,4] = compute_nu_Theil(int_h_Fpk[2], m01_hat)
	bounds_DTACW[4,(3,4)] = bounds_DTACW[3,(3,4)]
	
	/* Atkinson */
	int_h_Fpk = Def_int_h_Fpk_Atkinson(raw_ML_Kpind, m01_hat, b_atkinson)
	bounds_DTACW[5,3] = compute_nu_Atkinson(int_h_Fpk[1], m01_hat, b_atkinson)
	bounds_DTACW[5,4] = compute_nu_Atkinson(int_h_Fpk[2], m01_hat, b_atkinson)
	bounds_DTACW[6,(3,4)] = bounds_DTACW[5,(3,4)]
	
	/* Coworker */
	int_h_Fpk = Def_int_h_Fpk_Coworker(raw_ML_Kpind, m01_hat)
	bounds_DTACW[7,3] = compute_nu_Coworker(int_h_Fpk[1], m01_hat)
	bounds_DTACW[7,4] = compute_nu_Coworker(int_h_Fpk[2], m01_hat)
	bounds_DTACW[8,(3,4)] = bounds_DTACW[7,(3,4)]
	
	return(bounds_DTACW)
}
end

capture mata: mata drop Bounds_DTACW_Kpind_from_db()
mata:
mata set matastrict on
real matrix Bounds_DTACW_Kpind_from_db(real scalar b_atkinson, real matrix db, ///
							real scalar I_trace_estimation_type, ///
							struct options_optimization options_opt){

	struct raw_ML_results_Kpind raw_ML_Kpind
	if (I_trace_estimation_type){
		displayas("text")
		printf("K and p assumed independent: units are merged (maximal size = %f)\n", rows(db))
		displayflush()
	}
	raw_ML_Kpind = Estimation_ML_Kpind(db, options_opt.CML)
	return(Bounds_DTACW_Kpind(b_atkinson, raw_ML_Kpind, options_opt))						
}
end
