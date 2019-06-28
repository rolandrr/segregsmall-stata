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


/******************************************************************************/
/* Exploration - old for test */
/******************************************************************************/
/*	

capture mata: mata drop Res_Duncan_A21_test()
mata:
mata set matastrict on
struct Res_Duncan_A21_test{
	real scalar max, min
	real matrix min_representation, max_representation
}
end

capture mata: mata drop _View_Res_Duncan_A21_test()
mata:
mata set matastrict on
void _View_Res_Duncan_A21_test(struct Res_Duncan_A21_test scalar res_Duncan_A21_test){
	printf("max")
	res_Duncan_A21_test.max
	printf("min")
	res_Duncan_A21_test.min
	printf("min_representation")
	res_Duncan_A21_test.min_representation
	printf("max_representation")
	res_Duncan_A21_test.max_representation
}
end

capture mata: mata drop _Get_pr_Res_Duncan_A21_test()
mata:
mata set matastrict on
real matrix _Get_pr_Res_Duncan_A21_test(struct Res_Duncan_A21_test scalar res_Duncan_A21_test, real scalar choice){
	if (choice == 1) return(res_Duncan_A21_test.min)
	if (choice == 2) return(res_Duncan_A21_test.max)
	if (choice == 3) return(res_Duncan_A21_test.min_representation)
	if (choice == 4) return(res_Duncan_A21_test.max_representation)
}
end

/* Approximate max and min for the Duncan by testing principal representations
one order higher that is varying in a single dimensional grid that is 
the power of higher order (comprised between the two principal
representation of the targeted vector of moments */
capture mata: mata drop int_h_Fp_hats_A21_uncons_Duncan2()
mata:
mata set matastrict on
struct Res_Duncan_A21_test scalar int_h_Fp_hats_A21_uncons_Duncan2(real colvector mu, real scalar m, ///
												| real scalar nb_candidates_arg) {

													
	struct Res_Duncan_A21_test scalar res_Duncan_A21_test 
	
	external nb_candidates_default
	real scalar nb_candidates
	nb_candidates = ((args() == 3) ? nb_candidates_arg : nb_candidates_default)		
																						
	struct principal_representations scalar pr, pr_next
	real matrix int_h_Fp_candidates
	real colvector grid_next_moment
	real scalar next_moment_min, next_moment_max, i
	
	int_h_Fp_candidates = J(nb_candidates-1, 2, .)
	pr = Principal_representations(mu)
	
	/* int_h_Fp with the initial principal representations */
	int_h_Fp_candidates[1,1] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr, 0), m)
	int_h_Fp_candidates[1,2] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr, 1), m)
	
	/* minimum and maximum moment of order length(mu)+1 */
	next_moment_min = Moment_discrete_distribution(Get_principal_representation(pr, 0), length(mu)+1)
	next_moment_max = Moment_discrete_distribution(Get_principal_representation(pr, 1), length(mu)+1)
	grid_next_moment = rangen(next_moment_min, next_moment_max, nb_candidates)
	
	for(i = 2; i <= nb_candidates-1; i++){
	pr_next = Principal_representations((mu\grid_next_moment[i]))
/*
Get_principal_representation(pr_next, 0)
Get_principal_representation(pr_next, 1)		
*/
		int_h_Fp_candidates[i,1] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr_next, 0), m)
		int_h_Fp_candidates[i,2] = compute_int_h_Fp_dis_Duncan(Get_principal_representation(pr_next, 1), m)
	}
/*
compute_nu_Duncan(min(int_h_Fp_candidates), m)
compute_nu_Duncan(max(int_h_Fp_candidates), m)
*/

/*	return((min(int_h_Fp_candidates), max(int_h_Fp_candidates)))*/
	/* this order since the derivative of nu_Duncan w.r.t. u 
	(cf. DHR notation) is positive for any v in [0,1] */

real scalar i_min_pr0, i_min_pr1, i_max_pr0, i_max_pr1
real matrix w
real scalar I_min_in_pr0, I_max_in_pr0
real scalar i_min, i_max
struct principal_representations scalar pr_min, pr_max
real matrix min_representation, max_representation

	i_min_pr0 = .
	i_min_pr1 = .
	w = .
	minindex(int_h_Fp_candidates[,1], 1, i_min_pr0, w)
	minindex(int_h_Fp_candidates[,2], 1, i_min_pr1, w)
	I_min_in_pr0 = (int_h_Fp_candidates[i_min_pr0,1] < int_h_Fp_candidates[i_min_pr1,2])
	i_min = (I_min_in_pr0 ? i_min_pr0 : i_min_pr1)
	pr_min = Principal_representations((mu\grid_next_moment[i_min]))
	min_representation = Get_principal_representation(pr_min, 1-I_min_in_pr0)
/*
printf("ceci est la representation qui atteint le min \n")	
min_representation	
*/	
	i_max_pr0 = .
	i_max_pr1 = .
	w = .
	maxindex(int_h_Fp_candidates[,1], 1, i_max_pr0, w)
	maxindex(int_h_Fp_candidates[,2], 1, i_max_pr1, w)
	I_max_in_pr0 = (int_h_Fp_candidates[i_max_pr0,1] > int_h_Fp_candidates[i_max_pr1,2])
	i_max = (I_max_in_pr0 ? i_max_pr0 : i_max_pr1)
	pr_max = Principal_representations((mu\grid_next_moment[i_max]))
	max_representation = Get_principal_representation(pr_max, 1-I_max_in_pr0)
/*	
printf("ceci est la representation qui atteint le max \n")	
max_representation		
*/	
printf("i_min puis i_max I_min_in_pr0 I_max_in_pr0 pour voir si viendrait meme point sur la ligne \n")
i_min
i_max
I_min_in_pr0
I_max_in_pr0


res_Duncan_A21_test.min = compute_nu_Duncan(min(int_h_Fp_candidates), m)
res_Duncan_A21_test.max = compute_nu_Duncan(max(int_h_Fp_candidates), m)
res_Duncan_A21_test.min_representation = min_representation
res_Duncan_A21_test.max_representation = max_representation
	
	/*
	return((compute_nu_Duncan(min(int_h_Fp_candidates), m), compute_nu_Duncan(max(int_h_Fp_candidates), m)))
	*/
	return(res_Duncan_A21_test)
}
end

*/

