version 14.2
/* This do-file performs the optimization in the case of assumption 2.1 only,
Actually, following comparisons with Matlab optimization, additional theoretical
results, no numerical optimization is performed (constraints appear to be too
hard to integrate in Mata optimizer) ; a search is made along the feasible
canonical representations */

/* return min and max denoted eta up eta low in page 2 of Appendix (B.1)
mu is the vector of moment to match, corresponds to m^k in the paper
m corresponds to m_{01} in the paper 
nb_candidates_cr is to define the degree of discretization - default will be 400 */
capture mata: mata drop int_h_Fp_hats_A21_un_Duncan_cr()
mata:
mata set matastrict on
real rowvector int_h_Fp_hats_A21_un_Duncan_cr(real colvector mu, real scalar m, real scalar nb_candidates_cr) {
																			
	struct principal_representations scalar pr, pr_next
	real matrix int_h_Fp_candidates
	real colvector grid_next_moment
	real scalar next_moment_min, next_moment_max, i
	
	int_h_Fp_candidates = J(nb_candidates_cr-1, 2, .)
	pr = Principal_representations(mu)
	
	/* int_h_Fp with the initial principal representations */
	int_h_Fp_candidates[1,1] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr, 0), m)
	int_h_Fp_candidates[1,2] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr, 1), m)
	
	/* minimum and maximum moment of order length(mu)+1 */
	next_moment_min = Moment_discrete_distribution(Get_principal_representation(pr, 0), length(mu)+1)
	next_moment_max = Moment_discrete_distribution(Get_principal_representation(pr, 1), length(mu)+1)
	grid_next_moment = rangen(next_moment_min, next_moment_max, nb_candidates_cr)
	
	for(i = 2; i <= nb_candidates_cr-1; i++){
		pr_next = Principal_representations((mu\grid_next_moment[i]))
		int_h_Fp_candidates[i,1] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr_next, 0), m)
		int_h_Fp_candidates[i,2] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr_next, 1), m)
	}
	
	return((min(int_h_Fp_candidates), max(int_h_Fp_candidates)))
	/* this order since the derivative of nu_Duncan w.r.t. u 
	(cf. DHR notation) is positive for any v in [0,1] */
}
end
