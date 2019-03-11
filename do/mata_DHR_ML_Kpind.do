/* across K - likelihood summing up the different contribution K by K */

/******************************************************************************/
/* structure : init_estimation_Kpind */
/* Some basic statistics (initialization of estimation) */
/******************************************************************************/

/* + use structure results_prob_K_hat for empirical distribution of K for bootstrap */

capture mata: mata drop init_estimation_Kpind()
mata:
mata set matastrict on
struct init_estimation_Kpind {
	real colvector list_K_positive_obs
	real scalar nb_K_positive_obs
	real Kbar, Lbar, I_Kbar_odd
	real nb_obs
}
end

capture mata: mata drop Construct_init_estimation_Kpind()
mata:
mata set matastrict on
struct init_estimation_Kpind scalar Construct_init_estimation_Kpind(real matrix db){
	struct init_estimation_Kpind scalar init_esti_Kpind
	init_esti_Kpind.list_K_positive_obs = selectindex(db[,2])
	init_esti_Kpind.nb_K_positive_obs = length(init_esti_Kpind.list_K_positive_obs)
	init_esti_Kpind.Kbar = max(init_esti_Kpind.list_K_positive_obs)
	init_esti_Kpind.Lbar = floor((init_esti_Kpind.Kbar+1)/2)
	init_esti_Kpind.I_Kbar_odd = init_esti_Kpind.Kbar - 2*floor(init_esti_Kpind.Kbar/2)
	init_esti_Kpind.nb_obs = quadsum(db[,2],1)	
	return(init_esti_Kpind)
}
end

capture mata: mata drop _View_init_estimation_Kpind()
mata:
mata set matastrict on
void _View_init_estimation_Kpind(struct init_estimation_Kpind scalar init_esti_Kpind){
	printf("Struct init_estimation_Kpind \n")
	printf("list_K_positive_obs =\n")
	init_esti_Kpind.list_K_positive_obs
	printf("nb_K_positive_obs =")
	init_esti_Kpind.nb_K_positive_obs
	printf("Kbar =")
	init_esti_Kpind.Kbar
	printf("I_Kbar_odd =")
	init_esti_Kpind.I_Kbar_odd
	printf("Lbar =")
	init_esti_Kpind.Lbar
	printf("nb_obs =")
	init_esti_Kpind.nb_obs
}
end


/******************************************************************************/
/* structure UML_results_Kpind */
/* Equivalent of UML K per K */
/******************************************************************************/

capture mata: mata drop UML_results_Kpind()
mata:
mata set matastrict on
struct UML_results_Kpind {
	real colvector P_tilde, m_tilde
	real scalar I_m_tilde_in_M_tol
	real scalar K_max_info_m_tilde
}
end

capture mata: mata drop _View_UML_results_Kpind()
mata:
mata set matastrict on
void _View_UML_results_Kpind(struct UML_results_Kpind scalar UML_Kpind){
	printf("Structure UML_results_Kpind \n")
	printf("P_tilde = \n")
	UML_Kpind.P_tilde
	printf("m_tilde = \n")
	UML_Kpind.m_tilde
	printf("I_m_tilde_in_M_tol =")
	UML_Kpind.I_m_tilde_in_M_tol
	printf("K_max_info_m_tilde =")
	UML_Kpind.K_max_info_m_tilde
}
end

capture mata: mata drop Construct_UML_results_Kpind()
mata:
mata set matastrict on
struct UML_results_Kpind scalar Construct_UML_results_Kpind(real matrix db, struct init_estimation_Kpind scalar init_esti_Kpind){
	struct UML_results_Kpind scalar UML_Kpind
	real colvector m_tilde_Kpind_EqD2
	real scalar tolerance
	
	tolerance = 0
	UML_Kpind = UML_results_Kpind()

	/* Case with observations */
	if (init_esti_Kpind.nb_obs > 0){
		/* Test with Equation D2 if works */
		m_tilde_Kpind_EqD2 = Get_m_tilde_Kpind_EqD2(db, init_esti_Kpind)
		if (I_valid_mu(m_tilde_Kpind_EqD2)){
			UML_Kpind.m_tilde = m_tilde_Kpind_EqD2
			UML_Kpind.K_max_info_m_tilde = rows(UML_Kpind.m_tilde)
			UML_Kpind.I_m_tilde_in_M_tol = Test_mu_in_M(UML_Kpind.m_tilde, tolerance)
			UML_Kpind.P_tilde = . /* a priori, P_tilde not necessary in the case K and p independent to make estimation and CI */
		}
		else {
			UML_Kpind = Cons_UML_Kpind_from_UML_perK(db, init_esti_Kpind)
		}	
	}
	/* Case without observations */
	else {
		UML_Kpind.m_tilde = .
		UML_Kpind.P_tilde = .
		UML_Kpind.K_max_info_m_tilde = .
		UML_Kpind.I_m_tilde_in_M_tol = .
	}
	return(UML_Kpind)
}
end

/***************************************/
/* Using previous structure UML K by K */
/***************************************/

/* About the construction of m_tilde */
/* for any k, m_tilde of size k is the vector k x 1 whose x element is: E[p^x|K=k]
If K and p are independent, E[p^x|K=k] = E[p^x] for any x and any k
Idea:
(i) use function Construct_UML_results() from the K by K case to compute for each
k m_tilde
(ii) for each x, we have 1+(k-x) estimation of m_tilde(x)
(iii) under hypothesis Kpind, the same just variability with estimation, take 
a weighted average by the number of observations by units of size K
to have the idea of granting more weights to more precise (with more observations)
estimation of m_tilde
(iv) to get a vector m_tilde that will be used for initial starting values of CML
method */

/* OLD problem si pas info utiliser K_max_info_m_tilde */
/*
capture mata: mata drop Get_m_tilde_Kpind_meanUMLperK_o()
mata:
mata set matastrict on
real colvector Get_m_tilde_Kpind_meanUMLperK_o(real matrix db, struct init_estimation_Kpind scalar init_esti_Kpind){
	real matrix store_m_tilde_K_by_K
	real scalar index, K, order
	real rowvector list_not_na, list_K_high_enough
	real colvector m_tilde_Kpind
	store_m_tilde_K_by_K = J(init_esti_Kpind.Kbar, init_esti_Kpind.nb_K_positive_obs, .)
	for(index = 1; index <= init_esti_Kpind.nb_K_positive_obs; index++){
		K = init_esti_Kpind.list_K_positive_obs[index]
		store_m_tilde_K_by_K[1..K, index] = _Get_m_tilde_from_data(db[K, (1..(K+3))])
	}	
	m_tilde_Kpind = J(init_esti_Kpind.Kbar, 1, .)
	for(order = 1; order <= init_esti_Kpind.Kbar; order++){
		list_not_na = selectindex(store_m_tilde_K_by_K[order,] :!= .)
		if (length(list_not_na)) {
			m_tilde_Kpind[order] = ((db[init_esti_Kpind.list_K_positive_obs[list_not_na],2]') :/ quadsum((db[init_esti_Kpind.list_K_positive_obs[list_not_na],2]'),1)) * (store_m_tilde_K_by_K[order,list_not_na]')
		}
	}
	/* Check if any problem in m_tilde_Kpind
	if any take the m_tilde up to moment Lbar (the highest moment necessary to determine starting points 
	if again impossible, return . and it will be random starting points in the CML */
	if (sum(m_tilde_Kpind :== .)){
		list_K_high_enough = selectindex(init_esti_Kpind.list_K_positive_obs :>= init_esti_Kpind.Lbar)
		if (length(list_K_high_enough)){
			K = min(init_esti_Kpind.list_K_positive_obs[selectindex(init_esti_Kpind.list_K_positive_obs :>= init_esti_Kpind.Lbar)])
			m_tilde_Kpind = _Get_m_tilde_from_data(db[K, (1..(K+3))])
		}
		if (sum(m_tilde_Kpind :== .)){
			m_tilde_Kpind = .
		}
	}
	return(m_tilde_Kpind)
}
end
*/

capture mata: mata drop Cons_UML_Kpind_from_UML_perK()
mata:
mata set matastrict on
struct UML_results_Kpind scalar Cons_UML_Kpind_from_UML_perK(real matrix db, struct init_estimation_Kpind scalar init_esti_Kpind){
	struct UML_results_Kpind scalar UML_Kpind	
	struct UML_results colvector UML_perK
	real scalar index_K, K, order
	real colvector list_K_max_info_m_tilde_perK, list_not_na
	real matrix store_m_tilde_K_by_K
	
	UML_Kpind = UML_results_Kpind()
	
	/* Save UML K per K, including m_tilde */
	UML_perK = UML_results(init_esti_Kpind.nb_K_positive_obs)
	for(index_K = 1; index_K <= init_esti_Kpind.nb_K_positive_obs; index_K++){
		K = init_esti_Kpind.list_K_positive_obs[index_K]
		UML_perK[index_K] = Cons_UML_results_from_dataNks(db[K, (1..(K+3))])
	}

	/* Get the K_max_info_m_tilde across K */
	list_K_max_info_m_tilde_perK = J(init_esti_Kpind.nb_K_positive_obs, 1, .)
	for(index_K = 1; index_K <= init_esti_Kpind.nb_K_positive_obs; index_K++){
		list_K_max_info_m_tilde_perK[index_K] = UML_perK[index_K].K_max_info_m_tilde
	}
	UML_Kpind.K_max_info_m_tilde = max(list_K_max_info_m_tilde_perK)

	/* stock the different m_tilde for aggregating them */
	store_m_tilde_K_by_K = J(UML_Kpind.K_max_info_m_tilde, init_esti_Kpind.nb_K_positive_obs, .)
	for(index_K = 1; index_K <= init_esti_Kpind.nb_K_positive_obs; index_K++){
		K = init_esti_Kpind.list_K_positive_obs[index_K]
		store_m_tilde_K_by_K[(1..length(UML_perK[index_K].m_tilde)), index_K] = UML_perK[index_K].m_tilde
	}

	UML_Kpind.m_tilde = J(UML_Kpind.K_max_info_m_tilde, 1, .)
	for(order = 1; order <= UML_Kpind.K_max_info_m_tilde; order++){
		list_not_na = selectindex(store_m_tilde_K_by_K[order,] :!= .)
		if (length(list_not_na)) {
			UML_Kpind.m_tilde[order] = ((db[init_esti_Kpind.list_K_positive_obs[list_not_na],2]') :/ quadsum((db[init_esti_Kpind.list_K_positive_obs[list_not_na],2]'),1)) * (store_m_tilde_K_by_K[order,list_not_na]')
		}
	}
	
	/* UML_Kpind.I_m_tilde_in_M_tol */
	real scalar tolerance; tolerance = 0
	UML_Kpind.I_m_tilde_in_M_tol = Test_mu_in_M(UML_Kpind.m_tilde, tolerance)
	
	/* UML_Kpind.P_tilde */
	UML_Kpind.P_tilde = . /* a priori, P_tilde not necessary in the case K and p independent to make estimation and CI */
	
	return(UML_Kpind)
}
end

/*******************************************/
/* Using directly Equation D2 of DHR paper */
/*******************************************/
	
capture mata: mata drop Def_Q_stacked()
mata:
mata set matastrict on
real matrix Def_Q_stacked(struct init_estimation_Kpind scalar init_esti_Kpind){
	real matrix Q_stacked
	real scalar k, index_k, index_row
	Q_stacked = J(sum(init_esti_Kpind.list_K_positive_obs,1), init_esti_Kpind.Kbar, 0)
	index_row = 1
	for(index_k = 1; index_k <= init_esti_Kpind.nb_K_positive_obs; index_k++){
		k = init_esti_Kpind.list_K_positive_obs[index_k]
		Q_stacked[(index_row..(index_row+k-1)), (1..k)] = Def_Q(k)
		index_row = index_row + k
	}
	return(Q_stacked)
}
end
	
capture mata: mata drop Get_P_tilde_stacked()
mata:
mata set matastrict on
real colvector Get_P_tilde_stacked(real matrix db, struct init_estimation_Kpind scalar init_esti_Kpind){
	real colvector P_tilde_stacked
	real scalar k, index_k, index_row
	P_tilde_stacked = J(sum(init_esti_Kpind.list_K_positive_obs,1), 1, 0)
	index_row = 1
	for(index_k = 1; index_k <= init_esti_Kpind.nb_K_positive_obs; index_k++){
		k = init_esti_Kpind.list_K_positive_obs[index_k]
		P_tilde_stacked[(index_row..(index_row+k-1))] = Def_P_tilde(db[k,(1..k+3)])
		index_row = index_row + k
	}
	return(P_tilde_stacked)	
}
end	

capture mata: mata drop Get_m_tilde_Kpind_EqD2()
mata:
mata set matastrict on
real colvector Get_m_tilde_Kpind_EqD2(real matrix db, struct init_estimation_Kpind scalar init_esti_Kpind){
	real matrix Q_stacked
	real colvector P_tilde_stacked
	Q_stacked = Def_Q_stacked(init_esti_Kpind)
	P_tilde_stacked = Get_P_tilde_stacked(db, init_esti_Kpind)
	return(cholsolve(Q_stacked' * Q_stacked, Q_stacked' * P_tilde_stacked))
}
end


/******************************************************************************/
/* CML in case K and p are assumed independent (Appendix D.4) */
/* Follows the general structure of Get_xy_hat_CML(), in which optimization
is solved K by K, here only once with the whole likelihood across K 
So can be used in the K by K case also with an input matrix db that has only
one row for the case K by K without any assumption about the dependence between
K and p
In some trials, algorithms are quite close; enable coherence issue if use with
only one K, with or without option independence between K and p */
/******************************************************************************/

capture mata: mata drop CML_results_Kpind()
mata:
mata set matastrict on
struct CML_results_Kpind {
	real matrix xy_hat
	real colvector P_hat
	real scalar m1_hat /* m1_hat : first moment of the distribution xy_hat */	
}
end

capture mata: mata drop _View_CML_results_Kpind()
mata:
mata set matastrict on
void _View_CML_results_Kpind(struct CML_results_Kpind scalar CML_Kpind) {
	printf("Struct CML_results_Kpind - 3 elements \n")
	printf("xy_hat = \n")
	CML_Kpind.xy_hat
	printf("P_hat = \n")
	CML_Kpind.P_hat
	printf("m1_hat = \n")
	CML_Kpind.m1_hat
}
end

capture mata: mata drop Construct_CML_results_Kpind()
mata:
mata set matastrict on
struct CML_results_Kpind scalar Construct_CML_results_Kpind(real matrix db, ///
												struct init_estimation_Kpind scalar init_esti_Kpind, ///
												struct UML_results_Kpind scalar UML_Kpind, ///
												struct options_optimization_CML optionsoptCML) {

	struct CML_results_Kpind scalar CML_Kpind
	CML_Kpind = CML_results_Kpind()
	if (init_esti_Kpind.nb_obs > 0){
		CML_Kpind.xy_hat = Get_xy_hat_CML_K_p_independent(db, init_esti_Kpind.list_K_positive_obs, init_esti_Kpind.nb_K_positive_obs, init_esti_Kpind.Lbar, UML_Kpind.m_tilde, UML_Kpind.K_max_info_m_tilde, optionsoptCML)
		CML_Kpind.m1_hat = Moment_discrete_distribution(CML_Kpind.xy_hat, 1)
	}
	else {
		CML_Kpind.xy_hat = .
		CML_Kpind.m1_hat = .
	}
	return(CML_Kpind)
}
end

/******************************************************************************/
/* raw_ML (wrap-up) in case K and p are assumed independent (Appendix D.4) */
/******************************************************************************/

capture mata: mata drop raw_ML_results_Kpind()
mata:
mata set matastrict on
struct raw_ML_results_Kpind {
	struct init_estimation_Kpind scalar init_esti
	struct UML_results_Kpind scalar UML
	struct CML_results_Kpind scalar CML
	real scalar I_constrained_distribution
}
end

capture mata: mata drop _View_raw_ML_results_Kpind()
mata:
mata set matastrict on
void _View_raw_ML_results_Kpind(struct raw_ML_results_Kpind scalar raw_ML_Kpind) {
	printf("raw_ML_results_Kpind - 4 elements : 3 structures + 1 scalar \n")
	_View_init_estimation_Kpind(raw_ML_Kpind.init_esti)
	_View_UML_results_Kpind(raw_ML_Kpind.UML)
	_View_CML_results_Kpind(raw_ML_Kpind.CML)
	printf("I_constrained_distribution = ")
	raw_ML_Kpind.I_constrained_distribution
}
end

capture mata: mata drop Estimation_ML_Kpind()
mata:
mata set matastrict on
struct raw_ML_results_Kpind scalar Estimation_ML_Kpind(real matrix db, struct options_optimization_CML optionsoptCML){

	struct raw_ML_results_Kpind scalar raw_ML_Kpind
	raw_ML_Kpind = raw_ML_results_Kpind()
	raw_ML_Kpind.init_esti = Construct_init_estimation_Kpind(db)
	raw_ML_Kpind.UML = Construct_UML_results_Kpind(db, raw_ML_Kpind.init_esti)
	raw_ML_Kpind.CML = Construct_CML_results_Kpind(db, raw_ML_Kpind.init_esti, raw_ML_Kpind.UML, optionsoptCML) 
	if ((raw_ML_Kpind.UML.m_tilde[1] == 0) | (raw_ML_Kpind.UML.m_tilde[1] == 1)){
		raw_ML_Kpind.I_constrained_distribution = 1
	}
	else {
		raw_ML_Kpind.I_constrained_distribution = (!raw_ML_Kpind.UML.I_m_tilde_in_M_tol)
	}
	return(raw_ML_Kpind)
}
end


/******************************************************************************/
/* CML optimization */
/******************************************************************************/

/* Get_xy_hat_CML_K_p_independent() */
capture mata: mata drop Get_xy_hat_CML_K_p_independent()
mata:
mata set matastrict on
real matrix Get_xy_hat_CML_K_p_independent(real matrix db, ///
										real colvector list_K_positive_obs, ///
										real scalar nb_K_positive_obs, ///
										real scalar L, ///
										real colvector m_tilde, ///
										real scalar K_max_info_m_tilde, ///
										struct options_optimization_CML optionsoptCML) {

	/* Declaration variable (begin) */
	real scalar index, K
	real colvector count_units_X_eq_K	
	real matrix results_opti, results_opti_converged, store_opti, w
	real scalar nb_supp, remaining_nb_supp, nb_supp_upper_bound, init_technique, xB, mass_xB, i, nb_supp_hat
	real rowvector xtiytim10, xm1, x, ym1, y, ytildem1, x_draw, y_draw, xy_opti_max, x_hat_CML, y_hat_CML
	real scalar nb_use_alea_fail_opti, step_B_while_2, crit_continue_while_1, crit_continue_while_2, crit_continue_while_3
	real scalar ll_max_converged, ll_max_all, ll_max
	real colvector index_max_nbsupp_converged, index_max_nbsupp_all
	real colvector index_max_ll_store_opti
	transmorphic S
	/* Declaration variable (end) */


	/********************************************/
	/* Degenerate cases */
	/********************************************/
	
	/* If only K takes only one value in data, units have all the same size 
	get back to case K by K ; check and deviation at an higher level 
	when choosing the function used 
	here: assumes to have at least 2 modalities for K with >= 1 units for each */
	count_units_X_eq_K = J(nb_K_positive_obs, 1, .)
	for(index = 1; index <= nb_K_positive_obs; index++){
		K = list_K_positive_obs[index]
		count_units_X_eq_K[index] = db[K, K+3]
	}
	/* if only observations = units with zero minorities: X_i = 0 for all i */
	if (sum((db[list_K_positive_obs,3] :== db[list_K_positive_obs,2])) == nb_K_positive_obs){
		x_hat_CML = 0
		y_hat_CML = 1
	}
	/* if only observations = units with full minority individuals: X_i = K_i for all i */
	else if (sum(count_units_X_eq_K :== db[list_K_positive_obs,2]) == nb_K_positive_obs) {
		x_hat_CML = 1
		y_hat_CML = 1
	}
	
	/********************************************/
	/* Non-degenerated cases */
	/********************************************/
	
	/* Constrained Maximum Likelihood (CML) optimization */
	else {
	
		/***************************************************/
		/* Various options and parameters for optimization */
		/***************************************************/
		real scalar Lplus1;					Lplus1 = L + 1
		results_opti = J(Lplus1, 2*Lplus1+6, .)
		
		real scalar epsilon_h;				epsilon_h = 10^(-6)
		real scalar nb_h_grid; 				nb_h_grid = 30
		real scalar epsilon_x;				epsilon_x = 10^(-6)
		real scalar epsilon_y; 				epsilon_y = 10^(-5)
		real scalar threshold_xm1_0;		threshold_xm1_0 = 0.05
		real scalar threshold_xm1_1;		threshold_xm1_1 = 0.1
	
		real scalar nb0_alea_fail_opti;					nb0_alea_fail_opti = 15
		real scalar nb_max_iter_alea_fail;				nb_max_iter_alea_fail = 25
		real scalar nb0_alea_without_m_tilde;			nb0_alea_without_m_tilde = 15
		real scalar nb_max_iter_alea_without_m_tilde;	nb_max_iter_alea_without_m_tilde = 25
		real scalar nb_max_iter_conv_continue;			nb_max_iter_conv_continue = 35
		real scalar B_auto_continue;					B_auto_continue = 2
		real scalar max_while1_use_alea_fail_opti;		max_while1_use_alea_fail_opti = 3
		real scalar threshold_reldif_ll_max_while_3;	threshold_reldif_ll_max_while_3 = 10^(-6)
		real scalar threshold_reldif_max_ll_conv_all;	threshold_reldif_max_ll_conv_all = 10^(-9)
		real scalar step_1;								step_1 = 2
		real scalar step_2;								step_2 = 4
		real scalar step_3;								step_3 = 6
		real scalar step_4;								step_4 = 8		
	
		/******************************/
		/* begin first "while" boucle */
		/******************************/
		nb_supp = 0
		nb_use_alea_fail_opti = 0
		crit_continue_while_1 = 1
		while ((nb_supp <= Lplus1-1) & (crit_continue_while_1)) {
		
			nb_supp = nb_supp + 1
/*	
printf("Entree boucle while 1 pour B = \n")	
nb_supp			
*/
			/*************************/
			/* Particular case B = 1 */
			/*************************/
			if (nb_supp == 1) {

				S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
				optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
				optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
				optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
				optimize_init_params(S, get_xtiytim10(m_tilde[1]))

				/* in case of failure of the optimization */
				if (_optimize(S) != 0) {
					results_opti[nb_supp,(cols(results_opti)-5)] = 1
					continue
				}
				
				if (optimize_result_value(S) >= 0) {
					/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
					l'ajout de l'option missing mais au cas ou 
					on met une valeur . a la logvraisemblance (pour ne pas la retenir
					a la fin - le max par defaut ignore les valeurs manquantes .)
					et un indicateur de degenere */
					results_opti[nb_supp,(cols(results_opti)-5)] = 1
					results_opti[nb_supp,(cols(results_opti)-2)] = .
				}
				else {
					results_opti[nb_supp,1] = fun_logistic(optimize_result_params(S))
					results_opti[nb_supp,Lplus1+1] = 1
					results_opti[nb_supp,(cols(results_opti)-5)] = 0
					results_opti[nb_supp,(cols(results_opti)-4)] = nb_supp
					results_opti[nb_supp,(cols(results_opti)-3)] = optimize_result_value0(S)
					results_opti[nb_supp,(cols(results_opti)-2)] = optimize_result_value(S)
					results_opti[nb_supp,cols(results_opti)-1] = optimize_result_iterations(S)
					results_opti[nb_supp,cols(results_opti)] = optimize_result_converged(S)
				}
	
			}
			
			/***************/
			/* Case B > 1 */
			/***************/
			else {
		
				/* verification qu'on a bien les infos sur m_tilde pour definir les xtiytim10
				sinon, on passe a la partie avec des points de depart aleatoires */
				if (nb_supp <= K_max_info_m_tilde + 1) {
/*	
printf("OK, info sur m_tilde pour nb_supp = \n")
nb_supp	
*/			
					/* test avec les deux choix possibles pour initialisation (init) */
					store_opti = J(2, 2*nb_supp+6, .)

					for (init_technique = 1; init_technique <= 2; init_technique++) {
					
						/* avec x1 partant a 0 */
						if (init_technique == 1) {
							xtiytim10 = get_xtiytim10(get_xy_0_s0(nb_supp, m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
						/* avec x1 partant de h */
						else {
							xtiytim10 = get_xtiytim10(get_xy_0_sh(nb_supp, m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
					
						/* Cas ou probleme lors du calcul de x0 et y0 avec des . (lors inversion matrice) 
						on prend alors un point de depart provenant des points obtenus precedemment */                        
						if (sum(xtiytim10 :== .) > 0) { 
/*					
errprintf("recours a l'alea (xtiytim10 avec des .) \n")
*/
							xm1 = results_opti[nb_supp-1,1..nb_supp-1]
							if (min(xm1) >= threshold_xm1_0) {
									/* on prend un point proche de 0 s'il n'y en a pas deja */
									xB = (threshold_xm1_0 - epsilon_x) * runiform(1,1) + epsilon_x
							}
							else {
									/* sinon, on prend un point un peu plus grand que le dernier */
									xB = max(xm1) + threshold_xm1_1 * runiform(1,1)
							}
							x = (xm1, xB)   
							/* on attribue a ce point une probabilite arbitraire mass_xB */
							ym1 = results_opti[nb_supp-1, Lplus1+1..Lplus1+nb_supp-1]
							mass_xB = 1/nb_supp
							y = (ym1,mass_xB) :/ quadsum((ym1,mass_xB),1)
							xtiytim10 = get_xtiytim10((x\y))                        
						}
						
						S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
						optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
						optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
						optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
						optimize_init_params(S, xtiytim10)			
					
						if (_optimize(S) != 0) {
								/* prend les cas ou arret de la procedure car par exemple
								error 4 : Hessian is not negative semidefinite
								error 26: missing parameter valies not allowed */
								store_opti[init_technique,(cols(store_opti)-5)] = 1
								continue
						}

						if (optimize_result_value(S) >= 0) {
							/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
							l'ajout de l'option missing mais au cas ou 
							on met une valeur . a la logvraisemblance (pour ne pas la retenir
							a la fin - le max par defaut ignore les valeurs manquantes .)
							et un indicateur de degenere ou alors un autre probleme */
							store_opti[init_technique,(cols(store_opti)-5)] = 1
							store_opti[init_technique,(cols(store_opti)-2)] = .
						}
						else {
							/* cas normal */
							store_opti[init_technique,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
							ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
							store_opti[init_technique,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
							store_opti[init_technique,(cols(store_opti)-5)] = 0
							store_opti[init_technique,(cols(store_opti)-4)] = nb_supp
							store_opti[init_technique,(cols(store_opti)-3)] = optimize_result_value0(S)
							store_opti[init_technique,(cols(store_opti)-2)] = optimize_result_value(S)
							store_opti[init_technique,cols(store_opti)-1] = optimize_result_iterations(S)
							store_opti[init_technique,cols(store_opti)] = optimize_result_converged(S)
						}
					}
								
					/* en cas d'echec des deux optimisations, on en relance nb0_alea aleatoires */							
					if (sum(store_opti[.,(cols(store_opti)-5)]) == rows(store_opti)) {
/*	
errprintf("recours a l'alea (echec des optimisations avec xtiytim10) \n")	
*/			
						nb_use_alea_fail_opti = nb_use_alea_fail_opti + 1
						store_opti = J(nb0_alea_fail_opti, 2*nb_supp+6, .)
						for (i = 1; i <= nb0_alea_fail_opti; i++) {
							
							x_draw = runiform(1,nb_supp); y_draw = runiform(1,nb_supp); y_draw = y_draw :/ quadsum(y_draw,1); xtiytim10 = get_xtiytim10((x_draw \ y_draw))
							
							S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, nb_max_iter_alea_fail); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
							optimize_init_params(S, xtiytim10)	
						
							if (_optimize(S) != 0) {
									/* prend les cas ou arret de la procedure car par exemple
									error 4 : Hessian is not negative semidefinite
									error 26: missing parameter valies not allowed */
									store_opti[i,(cols(store_opti)-5)] = 1
									continue
							}
							
							if (optimize_result_value(S) >= 0) {
									/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
									l'ajout de l'option missing mais au cas ou 
									on met une valeur . a la logvraisemblance (pour ne pas la retenir
									a la fin - le max par defaut ignore les valeurs manquantes .)
									et un indicateur de degenere */
									store_opti[i,(cols(store_opti)-5)] = 1
									store_opti[i,(cols(store_opti)-2)] = .
							}
							else {
									/* cas normal */
									store_opti[i,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
									ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
									store_opti[i,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
									store_opti[i,(cols(store_opti)-5)] = 0
									store_opti[i,(cols(store_opti)-4)] = nb_supp
									store_opti[i,(cols(store_opti)-3)] = optimize_result_value0(S)
									store_opti[i,(cols(store_opti)-2)] = optimize_result_value(S)
									store_opti[i,cols(store_opti)-1] = optimize_result_iterations(S)
									store_opti[i,cols(store_opti)] = optimize_result_converged(S)
							}
						}
					}
				}
			
				/* pas d'info sur m_tilde, on part de points aleatoires */
				else {
/*				
errprintf("recours a l'alea (pas assez d'infos sur m_tilde, depart de points aleatoires), pour nb_supp = \n")
nb_supp
*/
					store_opti = J(nb0_alea_without_m_tilde, 2*nb_supp+6, .)
					
					for (i = 1; i <= nb0_alea_without_m_tilde; i++) {
						
						x_draw = runiform(1,nb_supp); y_draw = runiform(1,nb_supp); y_draw = y_draw :/ quadsum(y_draw,1); xtiytim10 = get_xtiytim10((x_draw \ y_draw))
						
						S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
						optimize_init_conv_maxiter(S, nb_max_iter_alea_without_m_tilde); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
						optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
						optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
						optimize_init_params(S, xtiytim10)
												
						if (_optimize(S) != 0) {
							/* prend les cas ou arret de la procedure car par exemple
							error 4 : Hessian is not negative semidefinite
							error 26: missing parameter valies not allowed */
							store_opti[i,(cols(store_opti)-5)] = 1
							continue
						}
						
						if (optimize_result_value(S) >= 0) {
							/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
							l'ajout de l'option missing mais au cas ou 
							on met une valeur . a la logvraisemblance (pour ne pas la retenir
							a la fin - le max par defaut ignore les valeurs manquantes .)
							et un indicateur de degenere */
							store_opti[i,(cols(store_opti)-5)] = 1
							store_opti[i,(cols(store_opti)-2)] = .
						}
						else {
							/* cas normal */
							store_opti[i,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
							ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
							store_opti[i,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
							store_opti[i,(cols(store_opti)-5)] = 0
							store_opti[i,(cols(store_opti)-4)] = nb_supp
							store_opti[i,(cols(store_opti)-3)] = optimize_result_value0(S)
							store_opti[i,(cols(store_opti)-2)] = optimize_result_value(S)
							store_opti[i,cols(store_opti)-1] = optimize_result_iterations(S)
							store_opti[i,cols(store_opti)] = optimize_result_converged(S)
						}
					}
				}
			
				/* selection de la maximum vraisemblance entre les deux points initiaux ou les points aleatoires 
				NB: max par defaut ne prend pas en consideration les valeurs manquantes ; different en comparaison
				immediate attention */
				w = .
				index_max_ll_store_opti = .
				maxindex(store_opti[.,(cols(store_opti)-2)], 1, index_max_ll_store_opti, w)
							
				if (length(index_max_ll_store_opti) <= 0) {
					/* que des . dans store_opti, aucune n'a marche : on essaye le suivant */
					continue
				}
				else {
					results_opti[nb_supp, 1..nb_supp] = store_opti[index_max_ll_store_opti[1], 1..nb_supp]
					results_opti[nb_supp, Lplus1+1..Lplus1+nb_supp] = store_opti[index_max_ll_store_opti[1], nb_supp+1..2*nb_supp]
					results_opti[nb_supp,(cols(results_opti)-5)] = 0
					results_opti[nb_supp,(cols(results_opti)-4)] = nb_supp
					results_opti[nb_supp,(cols(results_opti)-3)..cols(results_opti)] = store_opti[index_max_ll_store_opti[1],(cols(store_opti)-3)..cols(store_opti)]
				}
			} /* end case B > 1 */
		
			/* Critere poursuite boucle while 1 */
/*	
printf("avant critere poursuite boucle while 1, marque 5, valeur de B ici = \n")
nb_supp
*/
			if (nb_supp == 1) {
				ll_max = .
			}
			else {
/* printf("marque 6 entree dans le else \n") */
				ll_max = max(results_opti[1..nb_supp-1, (cols(store_opti)-2)])
			}
			
			crit_continue_while_1 = (nb_supp <= B_auto_continue) | ///
									(results_opti[nb_supp, cols(results_opti)] & (results_opti[nb_supp, cols(results_opti)-1] <= nb_max_iter_conv_continue)) | ///
									((results_opti[nb_supp, (cols(results_opti)-2)] != .) & (results_opti[nb_supp, (cols(results_opti)-2)] > ll_max) & (nb_use_alea_fail_opti <= max_while1_use_alea_fail_opti))
		}
		/****************************/
		/* end first "while" boucle */
		/****************************/
		
		remaining_nb_supp = Lplus1 - nb_supp
		
		if ((nb_supp <= 4) | (Lplus1 <= 6)) {
			step_B_while_2 = min((remaining_nb_supp, step_1))
		}
		else if ((nb_supp <= 6) | (Lplus1 <= 12)) {
			step_B_while_2 = min((remaining_nb_supp, step_2))
		}
		else if ((nb_supp <= 8) | (Lplus1 <= 14)) {
			step_B_while_2 = min((remaining_nb_supp, step_3))
		}
		else {
			step_B_while_2 = step_4
		}
/*
printf("marque 7, voici le step_B_while_2 = \n")
step_B_while_2
*/
		if (remaining_nb_supp > 0) {
		
			/*******************************/
			/* begin second "while" boucle */
			/*******************************/
			
			crit_continue_while_2 = 1
			nb_supp = nb_supp + step_B_while_2
			nb_use_alea_fail_opti = 0
			while((nb_supp <= Lplus1) & crit_continue_while_2) {
/*
printf("marque 8 entree boucle while 2 pour B = \n")	
nb_supp
*/
				/* verification qu'on a bien les infos sur m_tilde pour definir les xtiytim10
				sinon, on passe a la partie avec des points de depart aleatoires */
				if (nb_supp <= K_max_info_m_tilde + 1) {
/*		
printf("OK, info sur m_tilde pour nb_supp = \n")
nb_supp	
*/
					/* test avec les deux choix possibles pour initialisation (init) */
					store_opti = J(2, 2*nb_supp+6, .)

					for (init_technique = 1; init_technique <= 2; init_technique++) {
					
						/* avec x1 partant a 0 */
						if (init_technique == 1) {
							xtiytim10 = get_xtiytim10(get_xy_0_s0(nb_supp, m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
						/* avec x1 partant de h */
						else {
							xtiytim10 = get_xtiytim10(get_xy_0_sh(nb_supp, m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
					
						/* Cas ou probleme lors du calcul de x0 et y0 avec des . (lors inversion matrice) 
						on prend alors un point de depart provenant des points obtenus precedemment */                        
						if (sum(xtiytim10 :== .) > 0) { 
/*					
errprintf("recours a l'alea (xtiytim10 avec des .) \n")
*/
								xm1 = results_opti[nb_supp-1,1..nb_supp-1]
								if (min(xm1) >= threshold_xm1_0) {
										/* on prend un point proche de 0 s'il n'y en a pas deja */
										xB = (threshold_xm1_0 - epsilon_x) * runiform(1,1) + epsilon_x
								}
								else {
										/* sinon, on prend un point un peu plus grand que le dernier */
										xB = max(xm1) + threshold_xm1_1 * runiform(1,1)
								}
								x = (xm1, xB)   
								/* on attribue a ce point une probabilite arbitraire mass_xB */
								ym1 = results_opti[nb_supp-1, Lplus1+1..Lplus1+nb_supp-1]
								mass_xB = 1/nb_supp
								y = (ym1,mass_xB) :/ quadsum((ym1,mass_xB),1)
								xtiytim10 = get_xtiytim10((x\y))                        
							}

							S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
							optimize_init_params(S, xtiytim10)	
						
							if (_optimize(S) != 0) {
									/* prend les cas ou arret de la procedure car par exemple
									error 4 : Hessian is not negative semidefinite
									error 26: missing parameter valies not allowed */
									store_opti[init_technique,(cols(store_opti)-5)] = 1
									continue
							}

							if (optimize_result_value(S) >= 0) {
								/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
								l'ajout de l'option missing mais au cas ou 
								on met une valeur . a la logvraisemblance (pour ne pas la retenir
								a la fin - le max par defaut ignore les valeurs manquantes .)
								et un indicateur de degenere ou alors un autre probleme */
								store_opti[init_technique,(cols(store_opti)-5)] = 1
								store_opti[init_technique,(cols(store_opti)-2)] = .
							}
							else {
								/* cas normal */
								store_opti[init_technique,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
								ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
								store_opti[init_technique,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
								store_opti[init_technique,(cols(store_opti)-5)] = 0
								store_opti[init_technique,(cols(store_opti)-4)] = nb_supp
								store_opti[init_technique,(cols(store_opti)-3)] = optimize_result_value0(S)
								store_opti[init_technique,(cols(store_opti)-2)] = optimize_result_value(S)
								store_opti[init_technique,cols(store_opti)-1] = optimize_result_iterations(S)
								store_opti[init_technique,cols(store_opti)] = optimize_result_converged(S)
							}
						}
									
						/* en cas d'echec des deux optimisations, on en relance nb0_alea aleatoires */							
						if (sum(store_opti[.,(cols(store_opti)-5)]) == rows(store_opti)) {
/*				
errprintf("recours a l'alea (echec des optimisations avec xtiytim10) \n")	
*/					
							nb_use_alea_fail_opti = nb_use_alea_fail_opti + 1
							store_opti = J(nb0_alea_fail_opti, 2*nb_supp+6, .)
							for (i = 1; i <= nb0_alea_fail_opti; i++) {
								
								x_draw = runiform(1,nb_supp); y_draw = runiform(1,nb_supp); y_draw = y_draw :/ quadsum(y_draw,1); xtiytim10 = get_xtiytim10((x_draw \ y_draw))
								
								S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
								optimize_init_conv_maxiter(S, nb_max_iter_alea_fail); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
								optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
								optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
								optimize_init_params(S, xtiytim10)
															
								if (_optimize(S) != 0) {
										/* prend les cas ou arret de la procedure car par exemple
										error 4 : Hessian is not negative semidefinite
										error 26: missing parameter valies not allowed */
										store_opti[i,(cols(store_opti)-5)] = 1
										continue
								}
								
								if (optimize_result_value(S) >= 0) {
										/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
										l'ajout de l'option missing mais au cas ou 
										on met une valeur . a la logvraisemblance (pour ne pas la retenir
										a la fin - le max par defaut ignore les valeurs manquantes .)
										et un indicateur de degenere */
										store_opti[i,(cols(store_opti)-5)] = 1
										store_opti[i,(cols(store_opti)-2)] = .
								}
								else {
										/* cas normal */
										store_opti[i,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
										ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
										store_opti[i,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
										store_opti[i,(cols(store_opti)-5)] = 0
										store_opti[i,(cols(store_opti)-4)] = nb_supp
										store_opti[i,(cols(store_opti)-3)] = optimize_result_value0(S)
										store_opti[i,(cols(store_opti)-2)] = optimize_result_value(S)
										store_opti[i,cols(store_opti)-1] = optimize_result_iterations(S)
										store_opti[i,cols(store_opti)] = optimize_result_converged(S)
								}
							}
						}
					}
				
					/* pas assez d'info sur m_tilde, on part de points aleatoires */
					else {
/*					
errprintf("recours a l'alea (pas assez d'infor sur m_tilde, depart de points aleatoires), pour nb_supp = \n")
nb_supp
*/
						store_opti = J(nb0_alea_without_m_tilde, 2*nb_supp+6, .)
						
						for (i = 1; i <= nb0_alea_without_m_tilde; i++) {
							
							x_draw = runiform(1,nb_supp); y_draw = runiform(1,nb_supp); y_draw = y_draw :/ quadsum(y_draw,1); xtiytim10 = get_xtiytim10((x_draw \ y_draw))
							
							S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, nb_max_iter_alea_without_m_tilde); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
							optimize_init_params(S, xtiytim10)
														
							if (_optimize(S) != 0) {
								/* prend les cas ou arret de la procedure car par exemple
								error 4 : Hessian is not negative semidefinite
								error 26: missing parameter valies not allowed */
								store_opti[i,(cols(store_opti)-5)] = 1
								continue
							}
							
							if (optimize_result_value(S) >= 0) {
								/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
								l'ajout de l'option missing mais au cas ou 
								on met une valeur . a la logvraisemblance (pour ne pas la retenir
								a la fin - le max par defaut ignore les valeurs manquantes .)
								et un indicateur de degenere */
								store_opti[i,(cols(store_opti)-5)] = 1
								store_opti[i,(cols(store_opti)-2)] = .
							}
							else {
								/* cas normal */
								store_opti[i,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
								ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
								store_opti[i,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
								store_opti[i,(cols(store_opti)-5)] = 0
								store_opti[i,(cols(store_opti)-4)] = nb_supp
								store_opti[i,(cols(store_opti)-3)] = optimize_result_value0(S)
								store_opti[i,(cols(store_opti)-2)] = optimize_result_value(S)
								store_opti[i,cols(store_opti)-1] = optimize_result_iterations(S)
								store_opti[i,cols(store_opti)] = optimize_result_converged(S)
							}
						}
					}
				
					/* selection de la maximum vraisemblance entre les deux points initiaux ou les points aleatoires 
					NB: max par defaut ne prend pas en consideration les valeurs manquantes ; different en comparaison
					immediate attention */
					w = .; index_max_ll_store_opti = .; maxindex(store_opti[.,(cols(store_opti)-2)], 1, index_max_ll_store_opti, w)
								
					if (length(index_max_ll_store_opti) <= 0) {
						/* que des . dans store_opti, aucune n'a marche : on essaye le suivant */
						continue
					}
					else {
						results_opti[nb_supp, 1..nb_supp] = store_opti[index_max_ll_store_opti[1], 1..nb_supp]
						results_opti[nb_supp, Lplus1+1..Lplus1+nb_supp] = store_opti[index_max_ll_store_opti[1], nb_supp+1..2*nb_supp]
						results_opti[nb_supp,(cols(results_opti)-5)] = 0
						results_opti[nb_supp,(cols(results_opti)-4)] = nb_supp
						results_opti[nb_supp,(cols(results_opti)-3)..cols(results_opti)] = store_opti[index_max_ll_store_opti[1],(cols(store_opti)-3)..cols(store_opti)]
					}
				
				
				/* critere poursuite boucle while 2 */
				ll_max = max(results_opti[1..nb_supp-1, (cols(store_opti)-2)])
				crit_continue_while_2 = ((results_opti[nb_supp, (cols(store_opti)-2)] != .) & (results_opti[nb_supp, (cols(store_opti)-2)] > ll_max))
			
				if (crit_continue_while_2) {
					remaining_nb_supp = Lplus1 - nb_supp
					nb_supp = nb_supp + min((remaining_nb_supp+1, step_B_while_2))
				}
				
			} 
			/*******************************/
			/* end second "while" boucle */
			/*******************************/
		
			remaining_nb_supp = Lplus1 - nb_supp
/*
printf("marque 10 end of boucle while 2, value of remaining_nb_supp = \n")
remaining_nb_supp
*/			
			if (remaining_nb_supp > 0) {
			
				nb_supp_upper_bound = nb_supp
				nb_supp = nb_supp_upper_bound - (step_B_while_2 - 1)
				crit_continue_while_3 = 1
				nb_use_alea_fail_opti = 0
				
				/******************************/
				/* begin third "while" boucle */
				/******************************/
				while((nb_supp <= min((Lplus1, nb_supp_upper_bound))) & crit_continue_while_3) {
/*
printf("marque 11 entree boucle while 3, nb_supp = \n")				
nb_supp
*/	
					/* verification qu'on a bien les infos sur m_tilde pour definir les xtiytim10
					sinon, on passe a la partie avec des points de depart aleatoires */
					if (nb_supp <= K_max_info_m_tilde + 1) {
/*
printf("OK, info sur m_tilde pour nb_supp = \n")
nb_supp	
*/
					/* test avec les deux choix possibles pour initialisation (init) */
					store_opti = J(2, 2*nb_supp+6, .)

					for (init_technique = 1; init_technique <= 2; init_technique++) {
					
						/* avec x1 partant a 0 */
						if (init_technique == 1) {
							xtiytim10 = get_xtiytim10(get_xy_0_s0(nb_supp, m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
						/* avec x1 partant de h */
						else {
							xtiytim10 = get_xtiytim10(get_xy_0_sh(nb_supp, m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
					
						/* Cas ou probleme lors du calcul de x0 et y0 avec des . (lors inversion matrice) 
						on prend alors un point de depart provenant des points obtenus precedemment */                        
						if (sum(xtiytim10 :== .) > 0) { 
/*				
errprintf("recours a l'alea (xtiytim10 avec des .) \n")
*/
								xm1 = results_opti[nb_supp-1,1..nb_supp-1]
								if (min(xm1) >= threshold_xm1_0) {
										/* on prend un point proche de 0 s'il n'y en a pas deja */
										xB = (threshold_xm1_0 - epsilon_x) * runiform(1,1) + epsilon_x
								}
								else {
										/* sinon, on prend un point un peu plus grand que le dernier */
										xB = max(xm1) + threshold_xm1_1 * runiform(1,1)
								}
								x = (xm1, xB)   
								/* on attribue a ce point une probabilite arbitraire mass_xB */
								ym1 = results_opti[nb_supp-1, Lplus1+1..Lplus1+nb_supp-1]
								mass_xB = 1/nb_supp
								y = (ym1,mass_xB) :/ quadsum((ym1,mass_xB),1)
								xtiytim10 = get_xtiytim10((x\y))                        
							}
							
							S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
							optimize_init_params(S, xtiytim10)
													
							if (_optimize(S) != 0) {
									/* prend les cas ou arret de la procedure car par exemple
									error 4 : Hessian is not negative semidefinite
									error 26: missing parameter valies not allowed */
									store_opti[init_technique,(cols(store_opti)-5)] = 1
									continue
							}

							if (optimize_result_value(S) >= 0) {
								/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
								l'ajout de l'option missing mais au cas ou 
								on met une valeur . a la logvraisemblance (pour ne pas la retenir
								a la fin - le max par defaut ignore les valeurs manquantes .)
								et un indicateur de degenere ou alors un autre probleme */
								store_opti[init_technique,(cols(store_opti)-5)] = 1
								store_opti[init_technique,(cols(store_opti)-2)] = .
							}
							else {
								/* cas normal */
								store_opti[init_technique,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
								ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
								store_opti[init_technique,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
								store_opti[init_technique,(cols(store_opti)-5)] = 0
								store_opti[init_technique,(cols(store_opti)-4)] = nb_supp
								store_opti[init_technique,(cols(store_opti)-3)] = optimize_result_value0(S)
								store_opti[init_technique,(cols(store_opti)-2)] = optimize_result_value(S)
								store_opti[init_technique,cols(store_opti)-1] = optimize_result_iterations(S)
								store_opti[init_technique,cols(store_opti)] = optimize_result_converged(S)
							}
						}
									
						/* en cas d'echec des deux optimisations, on en relance nb0_alea aleatoires */							
						if (sum(store_opti[.,(cols(store_opti)-5)]) == rows(store_opti)) {
/*				
errprintf("recours a l'alea (echec des optimisations avec xtiytim10) \n")	
*/				
							nb_use_alea_fail_opti = nb_use_alea_fail_opti + 1
							store_opti = J(nb0_alea_fail_opti, 2*nb_supp+6, .)
							for (i = 1; i <= nb0_alea_fail_opti; i++) {
								
								x_draw = runiform(1,nb_supp); y_draw = runiform(1,nb_supp); y_draw = y_draw :/ quadsum(y_draw,1); xtiytim10 = get_xtiytim10((x_draw \ y_draw))
								
								S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
								optimize_init_conv_maxiter(S, nb_max_iter_alea_fail); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
								optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
								optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
								optimize_init_params(S, xtiytim10)

								if (_optimize(S) != 0) {
										/* prend les cas ou arret de la procedure car par exemple
										error 4 : Hessian is not negative semidefinite
										error 26: missing parameter valies not allowed */
										store_opti[i,(cols(store_opti)-5)] = 1
										continue
								}
								
								if (optimize_result_value(S) >= 0) {
										/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
										l'ajout de l'option missing mais au cas ou 
										on met une valeur . a la logvraisemblance (pour ne pas la retenir
										a la fin - le max par defaut ignore les valeurs manquantes .)
										et un indicateur de degenere */
										store_opti[i,(cols(store_opti)-5)] = 1
										store_opti[i,(cols(store_opti)-2)] = .
								}
								else {
										/* cas normal */
										store_opti[i,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
										ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
										store_opti[i,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
										store_opti[i,(cols(store_opti)-5)] = 0
										store_opti[i,(cols(store_opti)-4)] = nb_supp
										store_opti[i,(cols(store_opti)-3)] = optimize_result_value0(S)
										store_opti[i,(cols(store_opti)-2)] = optimize_result_value(S)
										store_opti[i,cols(store_opti)-1] = optimize_result_iterations(S)
										store_opti[i,cols(store_opti)] = optimize_result_converged(S)
								}
							}
						}
					}
				
					/* pas assez d'info sur m_tilde, on part de points aleatoires */
					else {
/*				
errprintf("recours a l'alea (pas assez d'infor sur m_tilde, depart de points aleatoires), pour nb_supp = \n")
nb_supp
*/
						store_opti = J(nb0_alea_without_m_tilde, 2*nb_supp+6, .)
						
						for (i = 1; i <= nb0_alea_without_m_tilde; i++) {
							
							x_draw = runiform(1,nb_supp); y_draw = runiform(1,nb_supp); y_draw = y_draw :/ quadsum(y_draw,1); xtiytim10 = get_xtiytim10((x_draw \ y_draw))
							
							S = optimize_init(); optimize_init_evaluator(S, &eval_CML_indep_Kp()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, nb_max_iter_alea_without_m_tilde); optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							optimize_init_argument(S, 1, db); optimize_init_argument(S, 2, list_K_positive_obs); optimize_init_argument(S, 3, nb_K_positive_obs)
							optimize_init_params(S, xtiytim10)
														
							if (_optimize(S) != 0) {
								/* prend les cas ou arret de la procedure car par exemple
								error 4 : Hessian is not negative semidefinite
								error 26: missing parameter valies not allowed */
								store_opti[i,(cols(store_opti)-5)] = 1
								continue
							}
							
							if (optimize_result_value(S) >= 0) {
								/* cas degenere - les arg max mis a 0 ou 1 (a priori n'arrive plus avec
								l'ajout de l'option missing mais au cas ou 
								on met une valeur . a la logvraisemblance (pour ne pas la retenir
								a la fin - le max par defaut ignore les valeurs manquantes .)
								et un indicateur de degenere */
								store_opti[i,(cols(store_opti)-5)] = 1
								store_opti[i,(cols(store_opti)-2)] = .
							}
							else {
								/* cas normal */
								store_opti[i,1..nb_supp] = fun_logistic(optimize_result_params(S)[1..nb_supp])
								ytildem1 = optimize_result_params(S)[nb_supp+1..2*nb_supp-1]
								store_opti[i,nb_supp+1..2*nb_supp] = (fun_logistic(ytildem1),1) :/ (1 + quadsum(fun_logistic(ytildem1),1))
								store_opti[i,(cols(store_opti)-5)] = 0
								store_opti[i,(cols(store_opti)-4)] = nb_supp
								store_opti[i,(cols(store_opti)-3)] = optimize_result_value0(S)
								store_opti[i,(cols(store_opti)-2)] = optimize_result_value(S)
								store_opti[i,cols(store_opti)-1] = optimize_result_iterations(S)
								store_opti[i,cols(store_opti)] = optimize_result_converged(S)
							}
						}
					}
				
					/* selection de la maximum vraisemblance entre les deux points initiaux ou les points aleatoires 
					NB: max par defaut ne prend pas en consideration les valeurs manquantes ; different en comparaison
					immediate attention */
					w = .; index_max_ll_store_opti = .; maxindex(store_opti[.,(cols(store_opti)-2)], 1, index_max_ll_store_opti, w)
								
					if (length(index_max_ll_store_opti) <= 0) {
						/* que des . dans store_opti, aucune n'a marche : on essaye le suivant */
						continue
					}
					else {
						results_opti[nb_supp, 1..nb_supp] = store_opti[index_max_ll_store_opti[1], 1..nb_supp]
						results_opti[nb_supp, Lplus1+1..Lplus1+nb_supp] = store_opti[index_max_ll_store_opti[1], nb_supp+1..2*nb_supp]
						results_opti[nb_supp,(cols(results_opti)-5)] = 0
						results_opti[nb_supp,(cols(results_opti)-4)] = nb_supp
						results_opti[nb_supp,(cols(results_opti)-3)..cols(results_opti)] = store_opti[index_max_ll_store_opti[1],(cols(store_opti)-3)..cols(store_opti)]
					}
				
				
					/* critere poursuite boucle while 3 */
					ll_max = max(results_opti[1..nb_supp-1, (cols(store_opti)-2)])
					crit_continue_while_3 = ((results_opti[nb_supp, (cols(store_opti)-2)] != .) & (results_opti[nb_supp, (cols(store_opti)-2)] > ll_max) & (nb_use_alea_fail_opti == 0)) | ///
											((results_opti[nb_supp, (cols(store_opti)-2)] != .) & (results_opti[nb_supp, (cols(store_opti)-2)] > ll_max) & (reldif(results_opti[nb_supp, (cols(store_opti)-2)], ll_max) > threshold_reldif_ll_max_while_3))
				
					if (crit_continue_while_3) {
						nb_supp = nb_supp + 1
					}
				
				}
				/****************************/
				/* end third "while" boucle */
				/****************************/
				
			} /* end if (remaining_nb_supp > 0) apres boucle while 2 */
		
		} /* end if (remaining_nb_supp > 0) apres boucle while 1*/
		
		
		/**************************************************************/
		/* Get final result that is xy_hat argmax of the optimization */
		/**************************************************************/
			
		/* max ll amongst the converged */
		results_opti_converged = select(results_opti, results_opti[.,cols(results_opti)])
		ll_max_converged = max(results_opti_converged[.,(cols(results_opti_converged)-2)])
		index_max_nbsupp_converged = .; w = .; maxindex(results_opti_converged[.,(cols(results_opti_converged)-2)], 1, index_max_nbsupp_converged, w)
		
		/* max ll amongst all */
		ll_max_all = max(results_opti[.,(cols(results_opti)-2)])
		index_max_nbsupp_all = .; w = .; maxindex(results_opti[.,(cols(results_opti)-2)], 1, index_max_nbsupp_all, w)
		
		if (sum(results_opti[index_max_nbsupp_all, cols(results_opti)])) {
			/* if converged amongs the arg max over all, choose the smallest B amongs converged */
			xy_opti_max = results_opti_converged[index_max_nbsupp_converged[1], .]
		}
		else {
			if (reldif(ll_max_converged, ll_max_all) <= threshold_reldif_max_ll_conv_all) {
				/* if converged enough close to not converged, choose the converged */	
				xy_opti_max = results_opti_converged[index_max_nbsupp_converged[1], .]	
			}
			else {
				xy_opti_max = results_opti[index_max_nbsupp_all[1], .]
			}
		}
		
		nb_supp_hat = xy_opti_max[(cols(results_opti)-4)]
		x_hat_CML = xy_opti_max[1..nb_supp_hat]
		y_hat_CML = xy_opti_max[Lplus1+1..Lplus1+nb_supp_hat]
	
	} /* end if cas non degenere */


	return (clean_xy(x_hat_CML \ y_hat_CML))
}
end
	
/******************************************************************************/
/* Subroutines for CML in case Kpind */	 
/******************************************************************************/

/* in addition to the evaluator, it uses subroutines from the K per K case */

/* Evaluator for the optimization in case K p independent, across K*/
/*******************************************************************/
/* Embranchements separes par les if selon la valeur de todo 
pour eviter de calculer forcement val, grad et hess
meme si en theorie, on precise le type d'evaluateur avec l'option d2
en pratique, la version v2 utilisee etait deux a trois fois plus lente
que la version v1 ou on ne calcule pas a chaque fois val grad et hess
Pour cela la version v1 calculait plusieurs fois les memes quantites
donc manque d efficacite 
d ou cette version 3 et les grands embranchements if 
*/

capture mata: mata drop S
capture mata: mata drop eval_CML_indep_Kp()
mata:
mata set matastrict on
void eval_CML_indep_Kp( real scalar todo, ///
						real rowvector xtildeytildem1, ///
						real matrix db, ///
						real colvector list_K_positive_obs, ///
						real scalar nb_K_positive_obs, ///
						val, grad, hess){
			
	/* Declaration variable (begin) */
	real scalar K
	real scalar index
	real colvector store_contribution_logL_by_K
	real matrix store_contribution_grad_by_K, cumulated_value_hess_K_by_K
	real scalar nb_supp, D
	real colvector Nk, sumjnbsuppyxk1mxKmk
	real rowvector xtilde, x, ytildem1, y
	real matrix xjk, onemxjKmk, xjkm1, onemxjKmkm1, kxjkm1onemxjKmk, KmkxjkonemxjKmkm1
	real matrix A, B, xjkm2, onemxjKmkm2, Aprime
	real rowvector W, Wpprime 
	real scalar j, j1, j2
	real colvector store_contribution_h_by_K_Beq1
	/* Declaration variable (end) */

	nb_supp = (length(xtildeytildem1)+1)/2

	if (todo == 0){ /* if todo = 0 (begin) */ 
	
		store_contribution_logL_by_K = J(nb_K_positive_obs, 1, .)
	
		if (nb_supp > 1) { /* if nb_supp > 1 (begin) */
			
			xtilde = xtildeytildem1[1..nb_supp]
			x = fun_logistic(xtilde)
			ytildem1 = xtildeytildem1[nb_supp+1..2*nb_supp-1]
			D = 1 + quadsum(fun_logistic(ytildem1),1) /* scalar D = 1 + sum_{i=1}^{nb_supp-1} Logistic(ytilde_i) */
			y = (fun_logistic(ytildem1),1) :/ D

			/* Computation of value (log-likelihood to maximize) */
			for(index = 1; index <= nb_K_positive_obs; index++){ /* loop over K (begin) */
				K = list_K_positive_obs[index]
				Nk = db[K, (3..(K+3))]'	
				/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
				xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
				/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
				onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
				/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
				sumjnbsuppyxk1mxKmk = quadrowsum(xjk :* y :* onemxjKmk, 1)
				/* val (last formula of Lemma 3.1) */
				store_contribution_logL_by_K[index] = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
			} /* loop over K (end) */
			
		} /* if nb_supp > 1 (end) */
		
		else { /* if nb_supp = 1 (begin) */	
			
			xtilde = xtildeytildem1[1]
			x = fun_logistic(xtilde)
			
			/* Computation of value (log-likelihood to maximize) */
			for(index = 1; index <= nb_K_positive_obs; index++){ /* loop over K (begin) */
				K = list_K_positive_obs[index]
				Nk = db[K, (3..(K+3))]'	
				/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
				xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
				/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
				onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
				/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
				sumjnbsuppyxk1mxKmk = xjk :* onemxjKmk
				/* val (last formula of Lemma 3.1) */
				store_contribution_logL_by_K[index] = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
			} /* loop over K (end) */
		} /* if nb_supp = 1 (end) */	
		
		val = quadsum(store_contribution_logL_by_K, 1)
		
	} /* if todo = 0 (end) */ 
	
	else if(todo == 1) { /* if todo = 1 (begin) */ 
		
		store_contribution_logL_by_K = J(nb_K_positive_obs, 1, .)
		store_contribution_grad_by_K = J(nb_K_positive_obs, 2*nb_supp-1, .)
		
		if (nb_supp > 1) { /* if nb_supp > 1 (begin) */

			xtilde = xtildeytildem1[1..nb_supp]
			x = fun_logistic(xtilde)
			ytildem1 = xtildeytildem1[nb_supp+1..2*nb_supp-1]
			D = 1 + quadsum(fun_logistic(ytildem1),1) /* scalar D = 1 + sum_{i=1}^{nb_supp-1} Logistic(ytilde_i) */
			y = (fun_logistic(ytildem1),1) :/ D		
		
			for(index = 1; index <= nb_K_positive_obs; index++){ /* loop over K (begin) */
				K = list_K_positive_obs[index]
				Nk = db[K, (3..(K+3))]'
		
				/* part computed in the computation of val (begin) */
				/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
				xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
				/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
				onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
				/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
				sumjnbsuppyxk1mxKmk = quadrowsum(xjk :* y :* onemxjKmk, 1)
				/* val (last formula of Lemma 3.1) */			
				store_contribution_logL_by_K[index] = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
				/* part computed in the computation of val (end) */
				
				/* part computed in the computation of grad (begin) */
				/* matrix xjkm1 : (x_j)^(k-1), row = k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 0 to have 0 and not . (0^(-1)) row k = 0 */
				xjkm1 = (J(K+1,1,1)*x):^((-1..K-1)'*J(1,nb_supp,1))
				if (sum((x:<=0)) > 0) {
					xjkm1[1,selectindex((x:<=0))] = J(1,length(selectindex((x:<=0))),0)
				}
				/* matrix onemxjKmkm1 : (1-x_j)^(K-k-1), row: k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 1 to have 0 and not . (0^(-1)), row k = K */
				onemxjKmkm1 = (J(K+1,1,1)*(1:-x)):^((K-1..-1)'*J(1,nb_supp,1))
				if (sum((x:>=1)) > 0) {
					onemxjKmkm1[rows(onemxjKmkm1),selectindex((x:>=1))] = J(1,length(selectindex((x:>=1))),0)
				}
				kxjkm1onemxjKmk = (0..K)' :* xjkm1 :* onemxjKmk
				KmkxjkonemxjKmkm1 = (K..0)' :* xjk :*  onemxjKmkm1
				A = y :* (kxjkm1onemxjKmk :- KmkxjkonemxjKmkm1)
				B = xjk :* onemxjKmk
				
				for (j = 1; j <= nb_supp-1; j++) {
				
					/* Derivative with respect to xtilde_j */
					store_contribution_grad_by_K[index, j] = fun_logistic_deriv1(xtilde[j]) * quadsum((A[.,j] :/ sumjnbsuppyxk1mxKmk) :* Nk,1)
						
					/* Derivative with respect to ytilde_j for j = 1..nb_supp-1 in case nb_supp >= 2 */
					store_contribution_grad_by_K[index, j+nb_supp] = (fun_logistic_deriv1(ytildem1[j]) / D) * quadsum(Nk :* (B[.,j] :- sumjnbsuppyxk1mxKmk) :/ sumjnbsuppyxk1mxKmk,1)
				}
				
				/* Derivative with respect to xtilde_j, for j = nb_supp */	
				store_contribution_grad_by_K[index, nb_supp] = fun_logistic_deriv1(xtilde[nb_supp]) * quadsum((A[.,nb_supp] :/ sumjnbsuppyxk1mxKmk) :* Nk,1)
				
				/* part computed in the computation of grad (end) */
			} /* loop over K (end) */
		} /* if nb_supp > 1 (end) */
		
		else { /* if nb_supp = 1 (begin) */
	
			xtilde = xtildeytildem1[1]
			x = fun_logistic(xtilde)
		
			for(index = 1; index <= nb_K_positive_obs; index++){ /* loop over K (begin) */
				K = list_K_positive_obs[index]
				Nk = db[K, (3..(K+3))]'
				
				/* part computed in the computation of val (begin) */
				/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
				xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
				/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
				onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
				/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
				sumjnbsuppyxk1mxKmk = xjk :* onemxjKmk
				/* val (last formula of Lemma 3.1) */
				store_contribution_logL_by_K[index] = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
				/* part computed in the computation of val (end) */
				
				/* part computed in the computation of grad (begin) */
				/* matrix xjkm1 : (x_j)^(k-1), row = k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 0 to have 0 and not . (0^(-1)) row k = 0 */
				xjkm1 = (J(K+1,1,1)*x):^((-1..K-1)'*J(1,nb_supp,1))
				if (sum((x:<=0)) > 0) {
					xjkm1[1,selectindex((x:<=0))] = J(1,length(selectindex((x:<=0))),0)
				}
				/* matrix onemxjKmkm1 : (1-x_j)^(K-k-1), row: k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 1 to have 0 and not . (0^(-1)), row k = K */
				onemxjKmkm1 = (J(K+1,1,1)*(1:-x)):^((K-1..-1)'*J(1,nb_supp,1))
				if (sum((x:>=1)) > 0) {
					onemxjKmkm1[rows(onemxjKmkm1),selectindex((x:>=1))] = J(1,length(selectindex((x:>=1))),0)
				}
				kxjkm1onemxjKmk = (0..K)' :* xjkm1 :* onemxjKmk
				KmkxjkonemxjKmkm1 = (K..0)' :* xjk :*  onemxjKmkm1
				A = kxjkm1onemxjKmk :- KmkxjkonemxjKmkm1
				/* Derivative with respect to xtilde (univariate function if nb_supp = 1) */			
				store_contribution_grad_by_K[index,1] = fun_logistic_deriv1(xtilde) * quadsum((A :/ sumjnbsuppyxk1mxKmk) :* Nk,1)
				/* part computed in the computation of grad (end) */
				
			} /* loop over K (end) */
		} /* if nb_supp = 1 (end) */
		
		val = quadsum(store_contribution_logL_by_K, 1)
		grad = quadcolsum(store_contribution_grad_by_K, 1)
	
	} /* if todo = 1 (end) */ 
	
	else { /* if todo = 2 (begin) */ 
	
		store_contribution_logL_by_K = J(nb_K_positive_obs, 1, .)
		store_contribution_grad_by_K = J(nb_K_positive_obs, 2*nb_supp-1, .)
	
		if (nb_supp > 1) { /* if nb_supp > 1 (begin) */

			xtilde = xtildeytildem1[1..nb_supp]
			x = fun_logistic(xtilde)
			ytildem1 = xtildeytildem1[nb_supp+1..2*nb_supp-1]
			D = 1 + quadsum(fun_logistic(ytildem1),1) /* scalar D = 1 + sum_{i=1}^{nb_supp-1} Logistic(ytilde_i) */
			y = (fun_logistic(ytildem1),1) :/ D
			
			cumulated_value_hess_K_by_K = J(2*nb_supp-1, 2*nb_supp-1, 0)
	
			for(index = 1; index <= nb_K_positive_obs; index++){ /* loop over K (begin) */
				K = list_K_positive_obs[index]
				Nk = db[K, (3..(K+3))]'
			
				/* part computed in the computation of val (begin) */
				/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
				xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
				/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
				onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
				/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
				sumjnbsuppyxk1mxKmk = quadrowsum(xjk :* y :* onemxjKmk, 1)
				/* val (last formula of Lemma 3.1) */
				store_contribution_logL_by_K[index] = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
				/* part computed in the computation of val (end) */
				
				/* part computed in the computation of grad (begin) */
				/* matrix xjkm1 : (x_j)^(k-1), row = k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 0 to have 0 and not . (0^(-1)) row k = 0 */
				xjkm1 = (J(K+1,1,1)*x):^((-1..K-1)'*J(1,nb_supp,1))
				if (sum((x:<=0)) > 0) {
					xjkm1[1,selectindex((x:<=0))] = J(1,length(selectindex((x:<=0))),0)
				}
				/* matrix onemxjKmkm1 : (1-x_j)^(K-k-1), row: k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 1 to have 0 and not . (0^(-1)), row k = K */
				onemxjKmkm1 = (J(K+1,1,1)*(1:-x)):^((K-1..-1)'*J(1,nb_supp,1))
				if (sum((x:>=1)) > 0) {
					onemxjKmkm1[rows(onemxjKmkm1),selectindex((x:>=1))] = J(1,length(selectindex((x:>=1))),0)
				}
				kxjkm1onemxjKmk = (0..K)' :* xjkm1 :* onemxjKmk
				KmkxjkonemxjKmkm1 = (K..0)' :* xjk :*  onemxjKmkm1
				A = y :* (kxjkm1onemxjKmk :- KmkxjkonemxjKmkm1)
				B = xjk :* onemxjKmk
				
				for (j = 1; j <= nb_supp-1; j++) {
			
					/* Derivative with respect to xtilde_j */
					store_contribution_grad_by_K[index, j] = fun_logistic_deriv1(xtilde[j]) * quadsum((A[.,j] :/ sumjnbsuppyxk1mxKmk) :* Nk,1)
						
					/* Derivative with respect to ytilde_j for j = 1..nb_supp-1 in case nb_supp >= 2 */
					store_contribution_grad_by_K[index, j+nb_supp] = (fun_logistic_deriv1(ytildem1[j]) / D) * quadsum(Nk :* (B[.,j] :- sumjnbsuppyxk1mxKmk) :/ sumjnbsuppyxk1mxKmk,1)
				}
				
				/* Derivative with respect to xtilde_j, for j = nb_supp */	
				store_contribution_grad_by_K[index, nb_supp] = fun_logistic_deriv1(xtilde[nb_supp]) * quadsum((A[.,nb_supp] :/ sumjnbsuppyxk1mxKmk) :* Nk,1)				
				/* part computed in the computation of grad (end) */
				
				/* part computed in the computation of hess (begin) */
				/* matrix xjkm2 : (x_j)^(k-2) - row: k = 0..K, column: j = 1..nb_supp
				with specific moodification in case x = 0 to have 0 and not . (0^(-1) and 0^(-2)), rows k = 0..1 */
				xjkm2 = (J(K+1,1,1)*x):^((-2..K-2)'*J(1,nb_supp,1))
				if (sum((x:<=0)) > 0) {
					xjkm2[1..2,selectindex((x:<=0))] = J(2,length(selectindex((x:<=0))),0)
				}
				
				/* matrix onemxjKmkm2 : (1-x_j)^(K-k-2) - row: k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 1 to have 0 and not . (0^(-1) and 0^(-2)), rows k = K..K-1 */
				onemxjKmkm2 = (J(K+1,1,1)*(1:-x)):^((K-2..-2)'*J(1,nb_supp,1))
				if (sum((x:>=1)) > 0) {
					onemxjKmkm2[rows(onemxjKmkm2)..rows(onemxjKmkm2)-1, selectindex((x:>=1))] = J(2,length(selectindex((x:>=1))),0)
				}
				
				/* matrix Aprime : derivative of matrix A with respect to x_j */
				/*	part1Aprime = (0..K)' :* ( ((-1..K-1)' :* xjkm2 :* onemxjKmk) :- ((K..0)' :*  xjkm1 :* onemxjKmkm1))
					part2Aprime = (K..0)' :* (((0..K)' :* xjkm1 :* onemxjKmkm1) :- ((K-1..-1)' :* xjk :* onemxjKmkm2)) */
				Aprime =  y :* ( ((0..K)' :* (((-1..K-1)' :* xjkm2 :* onemxjKmk) :- ((K..0)' :*  xjkm1 :* onemxjKmkm1))) :- ((K..0)' :* (((0..K)' :* xjkm1 :* onemxjKmkm1) :- ((K-1..-1)' :* xjk :* onemxjKmkm2))) )
				
				/* rowvector W : column: j = 1..nb_supp, the entry j is equal to :
				sum_{k=0}^K N_k * ( B_{k,j} - Sigma_k) / Sigma_k
				where Sigma_k = sumjnbsuppyxk1mxKmk[k] */	
				W = quadcolsum(Nk :* (B :- sumjnbsuppyxk1mxKmk) :/ sumjnbsuppyxk1mxKmk,1)
				
				/* rowvector Wpprime (pseudo-prime) : column: j = 1..nb_supp, the entry j is equal to :
				sum_{k=0}^K N_k * B_{k,j} * ( B_{k,j} - Sigma_k) / (Sigma_k)^2 */
				Wpprime = quadcolsum(Nk :* B :* ((B :- sumjnbsuppyxk1mxKmk) :/ (sumjnbsuppyxk1mxKmk:^2)),1)

				for (j1 = 1; j1 <= nb_supp; j1++) {	
					for (j2 = 1; j2 <= nb_supp; j2++) {
					
						if (j1 == j2) {
							
							/* Hessian with respect to xtilde_j1 and xtilde_j2 with j1 = j2 */
							cumulated_value_hess_K_by_K[j1,j2] = cumulated_value_hess_K_by_K[j1,j2] + (fun_logistic_deriv2(xtilde[j1]) * quadsum(Nk :* A[.,j1] :/ sumjnbsuppyxk1mxKmk,1)) + ((fun_logistic_deriv1(xtilde[j1])^2) * quadsum(Nk :* ((Aprime[.,j1] :* sumjnbsuppyxk1mxKmk) :- (A[.,j1]:^2)) :/ (sumjnbsuppyxk1mxKmk:^2),1))
							
							if (j2 < nb_supp) {	
							
								/* Hessian with respect to ytilde_j2 and xtilde_j1 with j1 = j2	*/
								cumulated_value_hess_K_by_K[nb_supp+j2,j1] = cumulated_value_hess_K_by_K[nb_supp+j2,j1] + (fun_logistic_deriv1(ytildem1[j1]) * fun_logistic_deriv1(xtilde[j1]) / D) * quadsum(Nk :* (((kxjkm1onemxjKmk :- KmkxjkonemxjKmkm1)[.,j1]) :/ (sumjnbsuppyxk1mxKmk:^2)) :* (sumjnbsuppyxk1mxKmk :- (B[.,j1] :* y[j1])),1)
							
								/* Hessian with respect to ytilde_j1 and ytilde_j2 with j1 = j2 */
								cumulated_value_hess_K_by_K[nb_supp+j1,nb_supp+j2] = cumulated_value_hess_K_by_K[nb_supp+j1,nb_supp+j2] + (W[j1] * fun_logistic_deriv2(ytildem1[j1]) / D) - (((fun_logistic_deriv1(ytildem1[j1]) / D)^2) * (W[j1] + Wpprime[j1]))
							
							}
						}
						
						else if (j2 < j1) {
						
							/* Hessian with respect to xtilde_j1 and xtilde_j2 with j1 != j2 */
							cumulated_value_hess_K_by_K[j1,j2] = cumulated_value_hess_K_by_K[j1,j2] - fun_logistic_deriv1(xtilde[j1]) * fun_logistic_deriv1(xtilde[j2]) * quadsum(Nk :* A[.,j1] :* A[.,j2] :/ (sumjnbsuppyxk1mxKmk:^2),1)
							
							/* Hessian with respect to ytilde_j2 and xtilde_j1, with j1 != j2 */
							cumulated_value_hess_K_by_K[nb_supp+j2,j1] = cumulated_value_hess_K_by_K[nb_supp+j2,j1] + (-fun_logistic_deriv1(ytildem1[j2]) * fun_logistic_deriv1(xtilde[j1]) / D) * quadsum(Nk :* B[.,j2] :* A[.,j1] :/ (sumjnbsuppyxk1mxKmk:^2),1)
							
							if (j1 < nb_supp) {
							
								/* Hessian with respect to ytilde_j1 and ytilde_j2 with j1 != j2 */
								cumulated_value_hess_K_by_K[nb_supp+j1,nb_supp+j2] = cumulated_value_hess_K_by_K[nb_supp+j1,nb_supp+j2] + (-fun_logistic_deriv1(ytildem1[j1]) * fun_logistic_deriv1(ytildem1[j2]) / (D^2)) * (W[j1] + quadsum(Nk :* B[.,j1] :* (B[.,j2] :- sumjnbsuppyxk1mxKmk) :/ (sumjnbsuppyxk1mxKmk:^2),1))
								
							}
						}
						
						else {
							/*(l2 > l1)*/ 
							
							if (j2 < nb_supp) {
							
								/* Hessian with respect to ytilde_j2 and xtilde_j1, with j1 != j2 */
								cumulated_value_hess_K_by_K[nb_supp+j2,j1] = cumulated_value_hess_K_by_K[nb_supp+j2,j1] + (-fun_logistic_deriv1(ytildem1[j2]) * fun_logistic_deriv1(xtilde[j1]) / D) * quadsum(Nk :* B[.,j2] :* A[.,j1] :/ (sumjnbsuppyxk1mxKmk:^2),1)
							
							}
						}
					}
				}				
			/* part computed in the computation of hess (end) */
			} /* loop over K (end) */
			
			hess = makesymmetric(cumulated_value_hess_K_by_K)
		} /* if nb_supp > 1 (end) */
		
		else { /* if nb_supp = 1 (begin) */
		
			xtilde = xtildeytildem1[1]
			x = fun_logistic(xtilde)		
		
			store_contribution_h_by_K_Beq1 = J(nb_K_positive_obs,1,.)
		
			for(index = 1; index <= nb_K_positive_obs; index++){ /* loop over K (begin) */
				K = list_K_positive_obs[index]
				Nk = db[K, (3..(K+3))]'
			
				/* part computed in the computation of val (begin) */
				/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
				xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
				/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
				onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
				/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
				sumjnbsuppyxk1mxKmk = xjk :* onemxjKmk
				/* val (last formula of Lemma 3.1) */
				store_contribution_logL_by_K[index] = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
				/* part computed in the computation of val (end) */
				
				/* part computed in the computation of grad (begin) */
				/* matrix xjkm1 : (x_j)^(k-1), row = k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 0 to have 0 and not . (0^(-1)) row k = 0 */
				xjkm1 = (J(K+1,1,1)*x):^((-1..K-1)'*J(1,nb_supp,1))
				if (sum((x:<=0)) > 0) {
					xjkm1[1,selectindex((x:<=0))] = J(1,length(selectindex((x:<=0))),0)
				}
				/* matrix onemxjKmkm1 : (1-x_j)^(K-k-1), row: k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 1 to have 0 and not . (0^(-1)), row k = K */
				onemxjKmkm1 = (J(K+1,1,1)*(1:-x)):^((K-1..-1)'*J(1,nb_supp,1))
				if (sum((x:>=1)) > 0) {
					onemxjKmkm1[rows(onemxjKmkm1),selectindex((x:>=1))] = J(1,length(selectindex((x:>=1))),0)
				}
				kxjkm1onemxjKmk = (0..K)' :* xjkm1 :* onemxjKmk
				KmkxjkonemxjKmkm1 = (K..0)' :* xjk :*  onemxjKmkm1
				A = kxjkm1onemxjKmk :- KmkxjkonemxjKmkm1
				/* Derivative with respect to xtilde (univariate function if nb_supp = 1) */			
				store_contribution_grad_by_K[index,1] = fun_logistic_deriv1(xtilde) * quadsum((A :/ sumjnbsuppyxk1mxKmk) :* Nk,1)
				/* part computed in the computation of grad (end) */
			
				/* part computed in the computation of hess (begin) */
				/* matrix xjkm2 : (x_j)^(k-2), row: k = 0..K, column: j = 1..nb_supp
				with specific moodification in case x = 0 to have 0 and not . (0^(-1) and 0^(-2)), rows k = 0..1 */
				xjkm2 = (J(K+1,1,1)*x):^((-2..K-2)'*J(1,nb_supp,1))
				if (sum((x:<=0)) > 0) {
					xjkm2[1..2,selectindex((x:<=0))] = J(2,length(selectindex((x:<=0))),0)
				}
				/* matrix onemxjKmkm2 : (1-x_j)^(K-k-2), row: k = 0..K, column: j = 1..nb_supp
				with specific modification in case x = 1 to have 0 and not . (0^(-1) and 0^(-2)), rows k = K..K-1 */
				onemxjKmkm2 = (J(K+1,1,1)*(1:-x)):^((K-2..-2)'*J(1,nb_supp,1))
				if (sum((x:>=1)) > 0) {
					onemxjKmkm2[rows(onemxjKmkm2)..rows(onemxjKmkm2)-1, selectindex((x:>=1))] = J(2,length(selectindex((x:>=1))),0)
				}
				/* matrix Aprime - derivative of matrix A with respect to x_j */
				/*	part1Aprime = (0..K)' :* ( ((-1..K-1)' :* xjkm2 :* onemxjKmk) :- ((K..0)' :*  xjkm1 :* onemxjKmkm1))
					part2Aprime = (K..0)' :* (((0..K)' :* xjkm1 :* onemxjKmkm1) :- ((K-1..-1)' :* xjk :* onemxjKmkm2)) */
				Aprime = ((0..K)' :* (((-1..K-1)' :* xjkm2 :* onemxjKmk) :- ((K..0)' :*  xjkm1 :* onemxjKmkm1))) :- ((K..0)' :* (((0..K)' :* xjkm1 :* onemxjKmkm1) :- ((K-1..-1)' :* xjk :* onemxjKmkm2)))
				store_contribution_h_by_K_Beq1[index] = (fun_logistic_deriv2(xtilde) * quadsum(Nk :* A :/ sumjnbsuppyxk1mxKmk,1)) + ((fun_logistic_deriv1(xtilde)^2) * quadsum(Nk :* ((Aprime :* sumjnbsuppyxk1mxKmk) :- (A:^2)) :/ (sumjnbsuppyxk1mxKmk:^2),1))
				/* part computed in the computation of hess (end) */
			} /* loop over K (end) */		
		
			hess = quadsum(store_contribution_h_by_K_Beq1, 1)
		} /* if nb_supp = 1 (end) */
	
		val = quadsum(store_contribution_logL_by_K, 1)
		grad = quadcolsum(store_contribution_grad_by_K, 1)
	
	} /* if todo = 2 (end) */
	
}
end


