/******************************************************************************/
/* This do-file creates the different mata structures and functions used to 
perform Maximum Likelihood (ML) estimation, both unconstrained and constrained
for the method of D'Haultfoeuille and Rathelot 2017 (DHR). 
If any, subfunctions used in a function are displayed just after the definition
of the latter.
oneK means that the ML is done in this file only for a fixed known K (setting of
DHR main text).
By contrast, the suffix "allK" will indicate that the corresponding structure or
functions act for the aggregate segregation indices when the size of units K is 
allowed to vary (cf. appendix B.1 in DHR 2017). */
/******************************************************************************/


/******************************************************************************/
/* structure : init_estimation */
/* Some basic statistics (initialization of estimation) */
/******************************************************************************/

capture mata: mata drop init_estimation()
mata:
mata set matastrict on
struct init_estimation {
	real scalar K, L, I_K_odd, nb_obs
	real vector list_Nks_Xm1_positive_obs /* list des X - 1 (Xm1) with observations */
}
end
/* NB: throughout, variables named I_xxx refers to the indicator of xxx and
takes values 0 (if xxx is false) or 1 (if xxx is true) */

capture mata: mata drop Construct_init_estimation()
mata:
mata set matastrict on
struct init_estimation scalar Construct_init_estimation(real rowvector data_Nks) {
	struct init_estimation scalar init_esti
	init_esti.K = data_Nks[1]
	init_esti.L = floor((init_esti.K+1)/2)
	init_esti.I_K_odd = init_esti.K - 2*floor(init_esti.K/2)
	init_esti.nb_obs = data_Nks[2]
	init_esti.list_Nks_Xm1_positive_obs = selectindex(data_Nks[3..length(data_Nks)])
	return(init_esti)
}
end

capture mata: mata drop _View_init_estimation()
mata:
mata set matastrict on
void _View_init_estimation(struct init_estimation scalar init_esti) {
	printf("Struct init_estimation - 5 elements \n")
	printf("K = ")
	init_esti.K
	printf("L = ")
	init_esti.L
	printf("I_K_odd = ")
	init_esti.I_K_odd
	printf("nb_obs = ")
	init_esti.nb_obs
}
end
/* NB: throughout, functions whose name is _xxx will indicate functions of 
internal use, be it to visualize structures (for debugging or check) or 
check some parts of the code */


/******************************************************************************/
/* structure : UML_results */
/* UML_results contains the results of the Unconstrained Maximum Likelihood 
estimation */
/******************************************************************************/

capture mata: mata drop UML_results()
mata:
mata set matastrict on
struct UML_results {
	real colvector P_tilde, m_tilde
	real scalar I_m_tilde_in_M_tol
	real scalar K_max_info_m_tilde
}
end

capture mata: mata drop Construct_UML_results()
mata:
mata set matastrict on
struct UML_results scalar Construct_UML_results(real rowvector data_Nks, ///
												struct init_estimation scalar init_esti) {
	struct UML_results scalar UML
	real colvector m_tilde_temp
	real scalar k, crit_continue_UML_m_tilde										
	
	/* case when there is some observations */
	if (init_esti.nb_obs > 0) {
		UML.P_tilde = Def_P_tilde(data_Nks)
		UML.m_tilde = lusolve(Def_Q(init_esti.K), UML.P_tilde)
		UML.K_max_info_m_tilde = init_esti.K
		
		/* Check for the computation of m_tilde using lusolve()
		Even if matrix Q is theoretically invertible whatever its size as an
		upper triangular matrix with non-zero diagonal elements, when K is large
		and thus Q (K x K matrix) is of large size, there may be numerical issues
		It is not really important in the sense that when K is large, basically
		we will be in the constrained case and therefore m_hat instead of m_tilde
		will be used. However, m_tilde is used to search for good starting values
		for the CML. Hence, the following computations to have at least the
		K_max_info_m_tilde first moments of m_tilde 
		Remark: qrsolve() was not used as experience reveal it yields not consistent
		results here (consistent with respect to the fact that m being the moments
		of a distribution in [0,1], they should be in [0,1] too 
		The different if conditions depending on how large K are here to speed up
		the procedure by starting from higher K for higher values of K
		Typically, numerical difficulties happen around K = 40 */
		
		if (!I_valid_mu(UML.m_tilde)) {
		
			/* case K <= 15 */
			if (init_esti.K <= 15) {
				k = 0
				crit_continue_UML_m_tilde = 1
				while ((k <= init_esti.K) & crit_continue_UML_m_tilde) {					
					k = k + 1
					m_tilde_temp = lusolve(Def_Q(k), UML.P_tilde[1..k])
					crit_continue_UML_m_tilde = I_valid_mu(m_tilde_temp)					
				}
			}
			/* case K > 15 & <= 30 */
			else if (init_esti.K <= 30) {
				k = 15
				crit_continue_UML_m_tilde = 1
				while ((k <= init_esti.K) & crit_continue_UML_m_tilde) {
					k = k + 1
					m_tilde_temp = lusolve(Def_Q(k), UML.P_tilde[1..k])
					crit_continue_UML_m_tilde = I_valid_mu(m_tilde_temp)	
				}
			}
			/* case K > 30 */
			else {
				k = 30
				crit_continue_UML_m_tilde = 1
				while ((k <= init_esti.K) & crit_continue_UML_m_tilde) {	
					k = k + 1
					m_tilde_temp = lusolve(Def_Q(k), UML.P_tilde[1..k])
					crit_continue_UML_m_tilde = I_valid_mu(m_tilde_temp)
				}
			}
				
			UML.K_max_info_m_tilde = k - 1
			UML.m_tilde = lusolve(Def_Q(UML.K_max_info_m_tilde), UML.P_tilde[1..UML.K_max_info_m_tilde])
			
			/* in case of problem, start from 1, takes more time but safer */
			if (!I_valid_mu(UML.m_tilde)) {
				k = 0
				crit_continue_UML_m_tilde = 1
				while ((k <= init_esti.K) & crit_continue_UML_m_tilde) {
					k = k + 1
					m_tilde_temp = lusolve(Def_Q(k), UML.P_tilde[1..k])
					crit_continue_UML_m_tilde = I_valid_mu(m_tilde_temp)
				}
				UML.K_max_info_m_tilde = k - 1
				UML.m_tilde = lusolve(Def_Q(UML.K_max_info_m_tilde), UML.P_tilde[1..UML.K_max_info_m_tilde])
			}
		}
		
		/* Check if the vector m_tilde belongs to the moment space M 
		to determine whether we are in the constrained or unconstrained case
		cf. infra - tolerance fixed = 0 by default */
		real scalar tolerance
		tolerance = 0
		UML.I_m_tilde_in_M_tol = Test_mu_in_M(UML.m_tilde, tolerance)
	}
	
	/* case when there is no observation */
	else {
		UML.P_tilde = J(init_esti.K, 1, .)
		UML.m_tilde = J(init_esti.K, 1, .)
		UML.I_m_tilde_in_M_tol = .
		UML.K_max_info_m_tilde = .	
	}
	return(UML)
}
end

capture mata: mata drop Cons_UML_results_from_dataNks()
mata:
mata set matastrict on
struct UML_results scalar Cons_UML_results_from_dataNks(real rowvector data_Nks){
	return(Construct_UML_results(data_Nks, Construct_init_estimation(data_Nks)))
}
end

/* Def_P_tilde() */
/* Definition of P_tilde (unconstrained maximum likelihood estimator) */
capture mata: mata drop Def_P_tilde()
mata:
mata set matastrict on
real colvector Def_P_tilde(real rowvector data_Nks) {	
	return((data_Nks[4..length(data_Nks)] :/ data_Nks[2])')
}
end
/* NB: Def_P_tilde() eminently uses the particular form of data data_Nks
cf. adofile db_from_unit_level.ado that constructs the database for details */

/* Def_Q() */
/* Creation of matrix Q for a given real scalar K */
capture mata: mata drop Def_Q()
mata:
mata set matastrict on
real matrix Def_Q(real scalar K) {
	real matrix Q
	real scalar i, j
	Q = J(K, K, 0)
	for (j = 1; j <= K; j++) {
		for (i = 1; i <= j; i++) {
			Q[i,j] = comb(K, j)*comb(j, i)*(-1)^(j-i)
		}
	}
	return(Q)
}
end

/* I_valid_mu() */
/* check whether mu can be a defined (without .) vector of moment of a 
distribution in [0,1] */
capture mata: mata drop I_valid_mu()
mata:
mata set matastrict on
real scalar I_valid_mu(real colvector mu) {
	return(!(sum(mu :== .) | (sum(mu :> 1)) | (sum(mu :< 0))))
}
end	

/* Test_mu_in_M()
Test that an input given colvector of reals belongs to the moment space M (truncated 
Hausdorff problem) ; conditions from Krein and Nudel'man, article DHR Appendix D.1. 
Proposition D.1 ; tol = tolerance for numerical issues in computing eigenvalues
in UML, fixed to 0 as the default option - choice made after several empirical
tests */
capture mata: mata drop Test_mu_in_M()
mata:
mata set matastrict on
real scalar Test_mu_in_M(real colvector mu, real scalar tol) {
	real scalar K, I_K_odd, L, i, j
	real scalar I_mu_in_M
	real matrix B, A, mC /* consistent with DHR's notation, mC = minus C */
	K = length(mu)
	I_K_odd = K - 2*floor(K/2)
	L = floor((K+1)/2)
	
	if (I_K_odd) { 	/* Case: K odd */
	
		/* test B_mu semi-definite positive */
		B = construct_B(L, mu)	
		if (min(symeigenvalues(B)) >= -tol) {
			/* test A_mu - B_mu semi-definite positive */
			A = construct_A(L, mu)		
			I_mu_in_M = (min(symeigenvalues(A-B)) >= -tol)
		}
		else {
			I_mu_in_M = 0
		}
	}
	
	else { /* Case: K even */
	
		/* test A_mu semi-definite positive */
		A = construct_A(L+1, mu)
		if (min(symeigenvalues(A)) >= -tol) {
			/* test - C_mu (mC) semi-definite positive */
			mC = J(L, L, 0)
			for (i = 1; i <= L; i++) {
				for (j = 1; j <= L; j++) {
					mC[i,j] = - mu[i+j] + mu[i+j-1]
				}
			}
			I_mu_in_M = (min(symeigenvalues(mC)) >= -tol)
		}
		else {
			I_mu_in_M = 0
		}
	}
	return(I_mu_in_M)
}
end

/* construct_A() */
/* Construct matrix A as defined in HR section 2.3 */
capture mata: mata drop construct_A()
mata:
mata set matastrict on
real matrix construct_A(real scalar size, real colvector mu) {
	real scalar i, j
	real matrix A; A = J(size, size, 1)
	for (j = 2; j <= size; j++) {
		A[1,j] = mu[1+j-2]
	}
	for (i = 2; i <= size; i++) {
		for (j = 1; j <= size; j++) {
			A[i,j] = mu[i+j-2]
		}
	}
	return(A)
}
end

/* construct_B() */
/* Construct matrix B as defined in HR section 2.3 */
capture mata: mata drop construct_B()
mata:
mata set matastrict on
real matrix construct_B(real scalar size, real colvector mu) {
	real scalar i, j
	real matrix B; 	B = J(size, size, .)
	for (i = 1; i <= size; i++) {
		for (j = 1; j <= size; j++) {
			B[i,j] = mu[i+j-1]
		}
	}
	return(B)
}
end

/* construct_C() */
/* Construct matrix C as defined in HR section 2.3 */
capture mata: mata drop construct_C()
mata:
mata set matastrict on
real matrix construct_C(real scalar size, real colvector mu) {
	real scalar i, j
	real matrix C; 	C = J(size, size, .)
	for (i = 1; i <= size; i++) {
		for (j = 1; j <= size; j++) {
			C[i,j] = mu[i+j] - mu[i+j-1]
		}
	}
	return(C)
}
end

capture mata: mata drop _View_UML_results()
mata:
mata set matastrict on
void _View_UML_results(struct UML_results scalar UML) {
	printf("Struct UML_results - 6 elements \n")
	printf("P_tilde = \n")
	UML.P_tilde
	printf("m_tilde = \n")
	UML.m_tilde
	printf("I_m_tilde_in_M_tol = ")
	UML.I_m_tilde_in_M_tol
	printf("K_max_info_m_tilde = ")
	UML.K_max_info_m_tilde
}
end

capture mata: mata drop _Get_UML_results()
mata:
mata set matastrict on
real matrix _Get_UML_results(struct UML_results scalar UML) {
	real matrix matrix_UML
	matrix_UML = J(length(UML.P_tilde), 4, .)
	matrix_UML[1,1] = UML.K_max_info_m_tilde
	matrix_UML[1,2] = UML.I_m_tilde_in_M_tol
	matrix_UML[.,3] = UML.P_tilde
	matrix_UML[.,4] = UML.m_tilde
	return(matrix_UML)
}
end

capture mata: mata drop _Get_mu_tilde()
mata:
mata set matastrict on
real colvector _Get_mu_tilde(struct UML_results scalar UML) {
	return(UML.m_tilde)
}
end

capture mata: mata drop _Get_m_tilde_from_data()
mata:
mata set matastrict on
real colvector _Get_m_tilde_from_data(real rowvector data_Nks){
	struct init_estimation scalar init_esti
	struct UML_results scalar UML
	init_esti = Construct_init_estimation(data_Nks)
	UML = Construct_UML_results(data_Nks, init_esti)
	return(UML.m_tilde)
}
end


/******************************************************************************/
/* structure : options_optimization_CML */
/* structure that contains the technical options for the mata optimize()
function to perform the Constrained Maximum Likelihood (CML) */
/******************************************************************************/

capture mata: mata drop options_optimization_CML()
mata:
mata set matastrict on
struct options_optimization_CML {
	real scalar nb_max_iter
	real scalar ptol, vtol, nrtol
}
end

capture mata: mata drop Cons_options_optimization_CML()
mata:
mata set matastrict on
struct options_optimization_CML scalar Cons_options_optimization_CML(	real scalar nb_max_iter, ///
																		real scalar ptol, ///
																		real scalar vtol, ///
																		real scalar nrtol) {
	struct options_optimization_CML scalar optionsoptCML
	optionsoptCML.nb_max_iter = nb_max_iter
	optionsoptCML.ptol = ptol
	optionsoptCML.vtol = vtol
	optionsoptCML.nrtol = nrtol
	return(optionsoptCML)
}
end

capture mata: mata drop _View_options_optimization_CML()
mata:
mata set matastrict on
void _View_options_optimization_CML(struct options_optimization_CML scalar optionsoptCML) {
	printf("Struct options_optimization_CML - 4 elements \n")
	printf("nb_max_iter = ")
	optionsoptCML.nb_max_iter
	printf("ptol = ")
	optionsoptCML.ptol
	printf("vtol = ")
	optionsoptCML.vtol
	printf("nrtol = ")
	optionsoptCML.nrtol
}
end


/******************************************************************************/
/* structure : CML_results */
/* CML_results contains the results of the Constrained Maximum Likelihood 
estimation - i.e. the results of the maximization program of Lemma 3.1 in DHR */
/******************************************************************************/

capture mata: mata drop CML_results()
mata:
mata set matastrict on
struct CML_results {
	real matrix xy_hat /* cf. infra Moment_discrete_distribution() for explanation */
	real colvector P_hat
	real scalar m1_hat /* m1_hat : first moment of the distribution xy_hat */	
}
end

capture mata: mata drop Construct_CML_results()
mata:
mata set matastrict on
struct CML_results scalar Construct_CML_results(real rowvector data_Nks, ///
												struct init_estimation scalar init_esti, ///
												struct UML_results scalar UML, ///
												struct options_optimization_CML optionsoptCML) {
	struct CML_results scalar CML
	CML = CML_results()
	if (init_esti.nb_obs > 0) {
		CML.xy_hat = Get_xy_hat_CML(data_Nks, init_esti, UML, optionsoptCML)
		CML.P_hat = Get_P_hat_CML(CML.xy_hat, init_esti.K)
		CML.m1_hat = Moment_discrete_distribution(CML.xy_hat, 1)
	}
	else {
		CML.xy_hat = J(2, 1, .)
		CML.P_hat = J(init_esti.K, 1, .)
		CML.m1_hat = .
	}
	return(CML)
}
end

capture mata: mata drop _View_CML_results()
mata:
mata set matastrict on
void _View_CML_results(struct CML_results scalar CML) {
	printf("Struct CML_results - 3 elements \n")
	printf("xy_hat = \n")
	CML.xy_hat
	printf("P_hat = \n")
	CML.P_hat
	printf("m1_hat = \n")
	CML.m1_hat
}
end

/* Moment_discrete_distribution() */
/* compute moment of order k's from a given discrete distribution
the distribution being represented by a matrix xy
whose first line (x) are the support points aka roots
and second line (y) are the corresponding weights aka masses */
capture mata: mata drop Moment_discrete_distribution()
mata:
mata set matastrict on
real scalar Moment_discrete_distribution(real matrix xy, real scalar k) {
	return((xy[1,.]:^k) * (xy[2,.]'))
}
end

/* Get_P_hat_CML() */
/* get P_hat from xy_hat using first formula of DHR Lemma 3.1 */
capture mata: mata drop Get_P_hat_CML()
mata:
mata set matastrict on
real colvector Get_P_hat_CML(real matrix xy_hat, real scalar K) {

	real rowvector x_hat, y_hat
	real scalar nb_supp_hat
	real matrix xjk_hat, onemxjKmk_hat
	real colvector sumjnbsuppyxk1mxKmk_hat, combkK
	
	x_hat = xy_hat[1,.]
	y_hat = xy_hat[2,.]
	nb_supp_hat = cols(xy_hat)
	
	/* matrix xjk_hat : ({x_hat}_j)^k, row: k = 1..K, column: j = 1..nb_supp_hat */
	xjk_hat = (J(K,1,1)*x_hat):^((1..K)'*J(1,nb_supp_hat,1))

	/* matrix onemxjKmk_hat : (1 - {x_hat}_j)^(K-k), row: k = 1..K, column: j = 1..nb_supp_hat */
	onemxjKmk_hat = (J(K,1,1)*(1:-x_hat)):^((K-1..0)'*J(1,nb_supp_hat,1))
	
	/* colvector sumjnbsuppyxk1mxKmk_hat : sum_{j=1}^{nb_supp_hat} {y_hat}_j ({x_hat}_j)^k (1-{x_hat}_j)^(K-k), row: k = 1..K */
	sumjnbsuppyxk1mxKmk_hat = quadrowsum(y_hat :* xjk_hat :* onemxjKmk_hat, 1)
	
	/* colvector combkK : choose k out of K, row: k = 1..K */
	combkK = comb(K,(1..K))'

	/* return colvector P_hat : P_hat[k] = comb(K,k) * sum_{j=1}^{nb_supp_hat} {y_hat}_j ({x_hat}_j)^k (1-{x_hat}_j)^(K-k),
	row: k = 1..K */ 
	return((combkK :* sumjnbsuppyxk1mxKmk_hat))
}	
end


/******************************************************************************/
/* structure : raw_ML_results */
/* raw_ML_results contains the results of the UML and CML 
(as previously for a fixed single K - setting of DHR main text */
/******************************************************************************/

capture mata: mata drop raw_ML_results()
mata:
mata set matastrict on
struct raw_ML_results {
	struct init_estimation scalar init_esti
	struct UML_results scalar UML
	struct CML_results scalar CML
	real scalar I_constrained_distribution
}
end

capture mata: mata drop _View_raw_ML_results()
mata:
mata set matastrict on
void _View_raw_ML_results(struct raw_ML_results scalar raw_ML) {
	printf("raw_ML_results - 4 elements : 3 structures + 1 scalar \n")
	_View_init_estimation(raw_ML.init_esti)
	_View_UML_results(raw_ML.UML)
	_View_CML_results(raw_ML.CML)
	printf("I_constrained_distribution = ")
	raw_ML.I_constrained_distribution
}
end

capture mata: mata drop _Get_m_tilde_raw_ML_results()
mata:
mata set matastrict on
real colvector _Get_m_tilde_raw_ML_results(struct raw_ML_results scalar raw_ML) {	
	return(_Get_mu_tilde(raw_ML.UML))
}
end

capture mata: mata drop _Get_xy_hat_raw_ML_results()
mata:
mata set matastrict on
real matrix _Get_xy_hat_raw_ML_results(struct raw_ML_results scalar raw_ML) {
	return(raw_ML.CML.xy_hat)
}
end

capture mata: mata drop _Get_xy_hat_dataNks()
mata:
mata set matastrict on
real matrix _Get_xy_hat_dataNks(real rowvector data_Nks, struct options_optimization_CML optionsoptCML) {
	struct raw_ML_results scalar raw_ML
	raw_ML = Estimation_ML(data_Nks, optionsoptCML)
	return(raw_ML.CML.xy_hat)
}
end

capture mata: mata drop Construct_raw_ML_results()
mata:
mata set matastrict on
struct raw_ML_results scalar Construct_raw_ML_results(	struct init_estimation scalar init_esti, ///
														struct UML_results scalar UML, ///
														struct CML_results scalar CML) {
	struct raw_ML_results scalar raw_ML
	raw_ML.init_esti = init_esti
	raw_ML.UML = UML
	raw_ML.CML = CML

	if (init_esti.nb_obs > 0) {
		if ((UML.m_tilde[1] == 0) | (UML.m_tilde[1] == 1)) {
			raw_ML.I_constrained_distribution = 1
		}
		else {
			raw_ML.I_constrained_distribution = (!UML.I_m_tilde_in_M_tol)
			/* DHR Appendix D.1 - test whether the distribution of p F_p is uniquely
			determined by its moment i.e. m_hat in the frontier of moment space M,
			which is implied by m_tilde not in the moment space M 
			(the test implied by DHR Proposition D.1 can be made with a certain tolerance
			default option is tolerance = 0 - cf. supra in the constructor of UML_results) */
		}
	}
	else {
		raw_ML.I_constrained_distribution = .
	}
	return(raw_ML)
}
end


/* Estimation_ML() */
/* Estimation_ML() performs the entire ML estimation in DHR, namely the first
init_esti step, UML and CML ; it yields as output a structure raw_ML_results */
capture mata: mata drop Estimation_ML()
mata:
mata set matastrict on
struct raw_ML_results scalar Estimation_ML(	real rowvector data_Nks, ///
											| struct options_optimization_CML optionsoptCML_arg) {
	external optionsoptCML_default
	struct options_optimization_CML scalar optionsoptCML
	optionsoptCML = ((args() == 2) ? optionsoptCML_arg : optionsoptCML_default)	
	
	struct init_estimation scalar init_esti
	struct UML_results scalar UML
	struct CML_results scalar CML
	init_esti = Construct_init_estimation(data_Nks)
	UML = Construct_UML_results(data_Nks, init_esti)
	CML = Construct_CML_results(data_Nks, init_esti, UML, optionsoptCML)
	return(Construct_raw_ML_results(init_esti, UML, CML))
}
end


/******************************************************************************/
/* Constrained Maximum Likelihood */
/* Functions that solve the CML problem (HR Lemma 3.1) */
/******************************************************************************/

/* OLD : idea to wrap-up (wu) with the only function for Kpind */
/*
capture mata: mata drop	Get_xy_hat_CML_wu()
mata:
mata set matastrict on
real matrix Get_xy_hat_CML_wu(	real rowvector data_Nks, ///
								struct init_estimation scalar init_esti, ///
								struct UML_results scalar UML, /// 
								struct options_optimization_CML optionsoptCML) {

	/* wrap-up function to use the Get_xy_hat_CML_K_p_independent()
	with only one K, as if K fixed so ok */
	real matrix db_for_one_K
	real scalar k
	db_for_one_K = J(init_esti.K, init_esti.K+3, 0)
	db_for_one_K[,1] = (1..init_esti.K)'
	db_for_one_K[init_esti.K,] = data_Nks
	if(init_esti.K > 1){
		for(k = 1; k < init_esti.K; k++){
			db_for_one_K[k,((k+4)..(init_esti.K+3))] = J(1, (init_esti.K-k), -99)
		}
	}
	return(Get_xy_hat_CML_K_p_independent(db_for_one_K, init_esti.K, 1, init_esti.L, UML.m_tilde, UML.K_max_info_m_tilde, optionsoptCML))
}
end	
*/	

/******************************************************************************/
/* Get_xy_hat_CML() 
get xy_hat from constrained maximum likelihood optimization - Lemma 3.1 */
capture mata: mata drop	Get_xy_hat_CML()
mata:
mata set matastrict on
real matrix Get_xy_hat_CML(	real rowvector data_Nks, ///
							struct init_estimation scalar init_esti, ///
							struct UML_results scalar UML, /// 
							struct options_optimization_CML scalar optionsoptCML) {

	
	real matrix results_opti, results_opti_converged, store_opti, w
	real scalar nb_supp, remaining_nb_supp, nb_supp_upper_bound, init_technique, xB, mass_xB, i, nb_supp_hat
	real rowvector xtiytim10, xm1, x, ym1, y, ytildem1, x_draw, y_draw, xy_opti_max, x_hat_CML, y_hat_CML
	real scalar nb_use_alea_fail_opti, step_B_while_2, crit_continue_while_1, crit_continue_while_2, crit_continue_while_3
	real scalar ll_max_converged, ll_max_all, ll_max
	real colvector index_max_nbsupp_converged, index_max_nbsupp_all
	real colvector index_max_ll_store_opti
	transmorphic S
		
	/********************************************/
	/* Degenerate cases */
	/********************************************/
	
	/* if only observations such that for all i, X_i = 0 */
	if (data_Nks[3] == data_Nks[2]) { 
		x_hat_CML = 0
		y_hat_CML = 1
	}
	
	/* if only observations such that for all i, X_i = K */
	else if (data_Nks[init_esti.K + 3] == data_Nks[2]) {
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
		real scalar Lplus1;					Lplus1 = init_esti.L + 1
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

				S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
				optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_params(S, get_xtiytim10(UML.m_tilde[1])); optimize_init_argument(S, 1, data_Nks)
				optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)	
				optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")

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
				if (nb_supp <= UML.K_max_info_m_tilde + 1) {
/*	
errprintf("OK, info sur m_tilde pour nb_supp = \n")
nb_supp	
*/
					/* test avec les deux choix possibles pour initialisation (init) */
					store_opti = J(2, 2*nb_supp+6, .)

					for (init_technique = 1; init_technique <= 2; init_technique++) {
					
						/* avec x1 partant a 0 */
						if (init_technique == 1) {
							xtiytim10 = get_xtiytim10(get_xy_0_s0(nb_supp, UML.m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
						/* avec x1 partant de h */
						else {
							xtiytim10 = get_xtiytim10(get_xy_0_sh(nb_supp, UML.m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
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

						S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
						optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
						optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)		
						optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
					
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
							
							S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, nb_max_iter_alea_fail); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
							optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)		
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
						
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
						
						S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
						optimize_init_conv_maxiter(S, nb_max_iter_alea_without_m_tilde); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
						optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)			
						optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
						
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
				if (nb_supp <= UML.K_max_info_m_tilde + 1) {
/*		
errprintf("OK, info sur m_tilde pour nb_supp = \n")
nb_supp	
*/
					/* test avec les deux choix possibles pour initialisation (init) */
					store_opti = J(2, 2*nb_supp+6, .)

					for (init_technique = 1; init_technique <= 2; init_technique++) {
					
						/* avec x1 partant a 0 */
						if (init_technique == 1) {
							xtiytim10 = get_xtiytim10(get_xy_0_s0(nb_supp, UML.m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
						/* avec x1 partant de h */
						else {
							xtiytim10 = get_xtiytim10(get_xy_0_sh(nb_supp, UML.m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
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

							S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
							optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)		
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
						
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
								
								S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
								optimize_init_conv_maxiter(S, nb_max_iter_alea_fail); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
								optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)		
								optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							
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
							
							S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, nb_max_iter_alea_without_m_tilde); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
							optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)			
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							
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
					if (nb_supp <= UML.K_max_info_m_tilde + 1) {
/*	
errprintf("OK, info sur m_tilde pour nb_supp = \n")
nb_supp	
*/
					/* test avec les deux choix possibles pour initialisation (init) */
					store_opti = J(2, 2*nb_supp+6, .)

					for (init_technique = 1; init_technique <= 2; init_technique++) {
					
						/* avec x1 partant a 0 */
						if (init_technique == 1) {
							xtiytim10 = get_xtiytim10(get_xy_0_s0(nb_supp, UML.m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
						}
						/* avec x1 partant de h */
						else {
							xtiytim10 = get_xtiytim10(get_xy_0_sh(nb_supp, UML.m_tilde, epsilon_h, nb_h_grid, epsilon_x, epsilon_y))
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

							S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, optionsoptCML.nb_max_iter); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
							optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)		
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
						
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
								
								S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
								optimize_init_conv_maxiter(S, nb_max_iter_alea_fail); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
								optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)		
								optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							
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
							
							S = optimize_init(); optimize_init_evaluator(S, &eval_xtiytim1_formulebas_mis_CML()); optimize_init_evaluatortype(S, "d2"); optimize_init_technique(S, "nr")
							optimize_init_conv_maxiter(S, nb_max_iter_alea_without_m_tilde); optimize_init_params(S, xtiytim10); optimize_init_argument(S, 1, data_Nks)
							optimize_init_conv_ptol(S, optionsoptCML.ptol); optimize_init_conv_vtol(S, optionsoptCML.vtol); optimize_init_conv_nrtol(S, optionsoptCML.nrtol)			
							optimize_init_verbose(S, 0); optimize_init_tracelevel(S, "none"); optimize_init_conv_warning(S, "off")
							
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
/* eval_xtiytim1_formulebas_mis_CML()
Evaluator function for optimization in Mata */
capture mata: mata drop S
capture mata: mata drop eval_xtiytim1_formulebas_mis_CML()
mata:
mata set matastrict on
void eval_xtiytim1_formulebas_mis_CML(	real scalar todo, ///
										real rowvector xtildeytildem1, ///
										real rowvector data_Nks, ///
										val, grad, hess) {

	real scalar K, nb_supp
	real colvector Nk, sumjnbsuppyxk1mxKmk
	real rowvector xtilde, x, ytildem1, y
	real matrix xjk, onemxjKmk, xjkm1, onemxjKmkm1, kxjkm1onemxjKmk, KmkxjkonemxjKmkm1
	real matrix A, B, D, xjkm2, onemxjKmkm2, Aprime
	real rowvector W, Wpprime 
	real scalar j, j1, j2
	
	K = data_Nks[1]
	Nk = data_Nks[3..length(data_Nks)]' /* colvector Nk : N_k, row: k = 0..K */
	nb_supp = (length(xtildeytildem1)+1)/2
	
	/* if nb_supp > 1 */
	if (nb_supp > 1) {
	
		xtilde = xtildeytildem1[1..nb_supp]
		x = fun_logistic(xtilde)
		ytildem1 = xtildeytildem1[nb_supp+1..2*nb_supp-1]
		D = 1 + quadsum(fun_logistic(ytildem1),1) /* scalar D = 1 + sum_{i=1}^{nb_supp-1} Logistic(ytilde_i) */
		y = (fun_logistic(ytildem1),1) :/ D 

		/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
		xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
	
		/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
		onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
		
		/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
		sumjnbsuppyxk1mxKmk = quadrowsum(xjk :* y :* onemxjKmk, 1)
		
		/* val (last formula of Lemma 3.1) */
		val = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
		
		if (todo >= 1) {
			/* Computation of gradient */
			
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
				grad[j] = fun_logistic_deriv1(xtilde[j]) * quadsum((A[.,j] :/ sumjnbsuppyxk1mxKmk) :* Nk,1)
					
				/* Derivative with respect to ytilde_j for j = 1..nb_supp-1 in case nb_supp >= 2 */
				grad[j+nb_supp] = (fun_logistic_deriv1(ytildem1[j]) / D) * quadsum(Nk :* (B[.,j] :- sumjnbsuppyxk1mxKmk) :/ sumjnbsuppyxk1mxKmk,1)
			}
			
			/* Derivative with respect to xtilde_j, for j = nb_supp */	
			grad[nb_supp] = fun_logistic_deriv1(xtilde[nb_supp]) * quadsum((A[.,nb_supp] :/ sumjnbsuppyxk1mxKmk) :* Nk,1)		
		
			if (todo == 2) {
				/* Computation of hessian */
				
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
							hess[j1,j2] = (fun_logistic_deriv2(xtilde[j1]) * quadsum(Nk :* A[.,j1] :/ sumjnbsuppyxk1mxKmk,1)) ///
											+ ((fun_logistic_deriv1(xtilde[j1])^2) * quadsum(Nk :* ((Aprime[.,j1] :* sumjnbsuppyxk1mxKmk) :- (A[.,j1]:^2)) :/ (sumjnbsuppyxk1mxKmk:^2),1))
							
							if (j2 < nb_supp) {	
							
								/* Hessian with respect to ytilde_j2 and xtilde_j1 with j1 = j2	*/
								hess[nb_supp+j2,j1] = (fun_logistic_deriv1(ytildem1[j1]) * fun_logistic_deriv1(xtilde[j1]) / D) * ///
														quadsum(Nk :* (((kxjkm1onemxjKmk :- KmkxjkonemxjKmkm1)[.,j1]) :/ (sumjnbsuppyxk1mxKmk:^2)) :* (sumjnbsuppyxk1mxKmk :- (B[.,j1] :* y[j1])),1)
							
								/* Hessian with respect to ytilde_j1 and ytilde_j2 with j1 = j2 */
								hess[nb_supp+j1,nb_supp+j2] = (W[j1] * fun_logistic_deriv2(ytildem1[j1]) / D) - (((fun_logistic_deriv1(ytildem1[j1]) / D)^2) * (W[j1] + Wpprime[j1]))
							
							}
						}
						
						else if (j2 < j1) {
						
							/* Hessian with respect to xtilde_j1 and xtilde_j2 with j1 != j2 */
							hess[j1,j2] = -fun_logistic_deriv1(xtilde[j1]) * fun_logistic_deriv1(xtilde[j2]) * quadsum(Nk :* A[.,j1] :* A[.,j2] :/ (sumjnbsuppyxk1mxKmk:^2),1)
							
							/* Hessian with respect to ytilde_j2 and xtilde_j1, with j1 != j2 */
							hess[nb_supp+j2,j1] = (-fun_logistic_deriv1(ytildem1[j2]) * fun_logistic_deriv1(xtilde[j1]) / D) * quadsum(Nk :* B[.,j2] :* A[.,j1] :/ (sumjnbsuppyxk1mxKmk:^2),1)
							
							if (j1 < nb_supp) {
							
								/* Hessian with respect to ytilde_j1 and ytilde_j2 with j1 != j2 */						
								hess[nb_supp+j1,nb_supp+j2] = (-fun_logistic_deriv1(ytildem1[j1]) * fun_logistic_deriv1(ytildem1[j2]) / (D^2)) * ///
																( W[j1] + quadsum(Nk :* B[.,j1] :* (B[.,j2] :- sumjnbsuppyxk1mxKmk) :/ (sumjnbsuppyxk1mxKmk:^2),1) )
								
							}
						}
						
						else {
							/*(l2 > l1)*/ 
							
							if (j2 < nb_supp) {
							
								/* Hessian with respect to ytilde_j2 and xtilde_j1, with j1 != j2 */
								hess[nb_supp+j2,j1] = (-fun_logistic_deriv1(ytildem1[j2]) * fun_logistic_deriv1(xtilde[j1]) / D) * quadsum(Nk :* B[.,j2] :* A[.,j1] :/ (sumjnbsuppyxk1mxKmk:^2),1)
							
							}
						}
					}
				}
				_makesymmetric(hess)
			}
		}
	}
	/* if nb_supp = 1 - just one argument : the point of a Dirac mass */
	else {
	
		xtilde = xtildeytildem1[1]
		x = fun_logistic(xtilde)
	
		/* computation of val */
		
		/* matrix xjk : (x_j)^k, row: k = 0..K, column: j = 1..nb_supp */
		xjk = (J(K+1,1,1)*x):^((0..K)'*J(1,nb_supp,1))
		
		/* matrix onemxjKmk : (1 - x_j)^(K-k), row: k = 0..K, column: j = 1..nb_supp */
		onemxjKmk = (J(K+1,1,1)*(1:-x)):^((K..0)'*J(1,nb_supp,1))
		
		/* colvector sumjnbsuppyxk1mxKmk : sum_{j=1}^{nb_supp} y_j (x_j)^k (1-x_j)^(K-k), row: k = 0..K */
		sumjnbsuppyxk1mxKmk = xjk :* onemxjKmk
		
		/* val (last formula of Lemma 3.1) */
		val = quadsum(ln(sumjnbsuppyxk1mxKmk) :* Nk,1)
	
		if (todo >= 1) {
			/* Computation of gradient */
			
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
			grad = fun_logistic_deriv1(xtilde) * quadsum((A :/ sumjnbsuppyxk1mxKmk) :* Nk,1)
			
			if (todo == 2) {
				/* Computation of hessian */
				
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
				
				hess = (fun_logistic_deriv2(xtilde) * quadsum(Nk :* A :/ sumjnbsuppyxk1mxKmk,1)) +  ((fun_logistic_deriv1(xtilde)^2) * quadsum(Nk :* ((Aprime :* sumjnbsuppyxk1mxKmk) :- (A:^2)) :/ (sumjnbsuppyxk1mxKmk:^2),1))
			}
		}
	}
}							
end


/******************************************************************************/
/* get_xy_0_sh() 
Define the initial starting point xy_0 for optimization
first row of returned matrix: roots (starting from h), second row: masses */
capture mata: mata drop get_xy_0_sh()
mata:
mata set matastrict on
real matrix get_xy_0_sh(	real scalar nb_supp, real colvector m_tilde, ///
							real scalar epsilon_h, real scalar nb_h_grid, ///
							real scalar epsilon_x, real scalar epsilon_y) {
	
	real matrix V, V_inv, store_y, index_max, w
	real colvector h_grid, m_mod_init, y
	real scalar i
	real rowvector x

	/* Computation of initial x and y that fit best the unconstrained moment
	vector m_tilde (starting from h for support points) */
	V = Vandermonde((1..nb_supp)')'
	V_inv = luinv(V)
	h_grid = rangen(epsilon_h, 1/nb_supp, nb_h_grid)
	store_y = J(length(h_grid), nb_supp+1, .)
	for (i = 1; i <= length(h_grid); i++) {
		m_mod_init = (1 \ m_tilde[1..nb_supp-1]) :/ J(nb_supp,1,h_grid[i]):^((0..nb_supp-1)')
		y = V_inv * m_mod_init
		store_y[i, .] = (y',min(y))
	}
	index_max = .
	w = .
	maxindex(store_y[.,nb_supp+1], 1, index_max, w)
	x = (1..nb_supp) :* h_grid[index_max[1]]
	y = store_y[index_max[1], 1..nb_supp]
	
	/* Modification of x and y to respect support points in (0,1) and probabilities */
	
	/* Modification of x */
	if (sum(x :<= epsilon_x) > 0) {
        x[selectindex(x :<= epsilon_x)] = J(1,length(selectindex(x :<= epsilon_x)),epsilon_x)
	}
	if (sum(x :>= 1-epsilon_x) > 0) {
        x[selectindex(x :>= 1-epsilon_x)] = J(1,length(selectindex(x :>= 1-epsilon_x)),1-epsilon_x)
	}

	/* Modification of y */
	if (sum(y :<= epsilon_y) > 0) {
        y[selectindex(y :<= epsilon_y)] = J(1,length(selectindex(y :<= epsilon_y)),epsilon_y)
	}
	y = y :/ quadsum(y,1)

	return((x\y))	
}
end


/******************************************************************************/
/* get_xy_0_s0() 
Define the initial starting point xy_0 for optimization
first row of returned matrix: roots (starting from 0), second row: masses */
capture mata: mata drop get_xy_0_s0()
mata:
mata set matastrict on
real matrix get_xy_0_s0(	real scalar nb_supp, real colvector m_tilde, ///
							real scalar epsilon_h, real scalar nb_h_grid, ///
							real scalar epsilon_x, real scalar epsilon_y) {

	real matrix V, V_inv, store_y, index_max, w
	real colvector h_grid, m_mod_init, y
	real scalar i
	real rowvector x

	/* Computation of initial x and y that fit best the unconstrained moment
	vector m_tilde (starting from h for support points) */
	V = Vandermonde((0..nb_supp-1)')'
	V_inv = luinv(V)
	h_grid = rangen(epsilon_h, 1/(nb_supp-1), nb_h_grid)
	store_y = J(length(h_grid), nb_supp+1, .)
	for (i = 1; i <= length(h_grid); i++) {
		m_mod_init = (1 \ m_tilde[1..nb_supp-1]) :/ J(nb_supp,1,h_grid[i]):^((0..nb_supp-1)')
		y = V_inv * m_mod_init
		store_y[i, .] = (y',min(y))
	}
	index_max = .
	w = .
	maxindex(store_y[.,nb_supp+1], 1, index_max, w)
	x = (0..nb_supp-1) :* h_grid[index_max[1]]
	y = store_y[index_max[1], 1..nb_supp]
	
	/* Modification of x and y to respect support points in (0,1) and probabilities */
	
	/* Modification of x */
	if (sum(x :<= epsilon_x) > 0) {
        x[selectindex(x :<= epsilon_x)] = J(1,length(selectindex(x :<= epsilon_x)),epsilon_x)
	}
	if (sum(x :>= 1-epsilon_x) > 0) {
        x[selectindex(x :>= 1-epsilon_x)] = J(1,length(selectindex(x :>= 1-epsilon_x)),1-epsilon_x)
	}

	/* Modification of y */
	if (sum(y :<= epsilon_y) > 0) {
        y[selectindex(y :<= epsilon_y)] = J(1,length(selectindex(y :<= epsilon_y)),epsilon_y)
	}
	y = y :/ quadsum(y,1)

	return((x\y))	
}
end

/******************************************************************************/
/* get_xtiytim10() 
get xtilde ytilde initial point to feed the evaluator and the optimization from 
the original xy_0
tilde corresponds to the logistic transformation for the x (roots) and a
logistic transformation plus a normalization (sum = 1) for the y (masses) */
capture mata: mata drop get_xtiytim10()
mata:
mata set matastrict on
real rowvector get_xtiytim10(real matrix xy0) {

	real rowvector x, y, xtilde, ym1, ytildem1, indexm1
	real scalar nb_supp
	real matrix index_max_y, w
	
	nb_supp = cols(xy0)
	
	if (nb_supp > 1){

		x = xy0[1,.]
		y = xy0[2,.]

		index_max_y = .
		w = .
		maxindex(y, 1, index_max_y, w)
		indexm1 = selectindex((1..nb_supp) :!= index_max_y[1])
		
		x = (x[indexm1], x[index_max_y[1]])
		xtilde = fun_logit(x)

		ym1 = (1 - y[index_max_y[1]]) :* y[indexm1] :/ quadsum(y[indexm1],1)
		ytildem1 = fun_logit(lusolve((ym1' * J(1, nb_supp-1, 1) - I(nb_supp-1)), -ym1'))'
		
		return((xtilde,ytildem1))
	}
	if (nb_supp == 1) {
		
		if (rows(xy0) == 1){
			xtilde = fun_logit(xy0)
		}
		else {
			xtilde = fun_logit(xy0[1])
		}
		return(xtilde)
	}
}
end

/******************************************************************************/
/* fun_logistic() 
usual logistic function */
capture mata: mata drop fun_logistic()
mata:
mata set matastrict on
numeric matrix fun_logistic(numeric matrix X) {
	return(1 :/ (1 :+ exp(-X)))
}
end

/* fun_logit() 
usual logit function = inverse of logistic function */
capture mata: mata drop fun_logit()
mata:
mata set matastrict on
numeric matrix fun_logit(numeric matrix X) {
	return(ln(X:/(1 :- X)))
}
end

/* fun_logistic_deriv1() 
first derivative of logistic function */
capture mata: mata drop fun_logistic_deriv1()
mata:
mata set matastrict on
numeric matrix fun_logistic_deriv1(numeric matrix X) {
	return(exp(X) :/ ((1 :+ exp(X)):^2))
}
end

/* fun_logistic_deriv2() 
second derivative of logistic function */
capture mata: mata drop fun_logistic_deriv2()
mata:
mata set matastrict on
numeric matrix fun_logistic_deriv2(numeric matrix X) {
	return((exp(X) :* (1 :- exp(X))) :/ ((1 :+ exp(X)):^3))	
}
end

/******************************************************************************/
/* clean_xy() */
capture mata: mata drop clean_xy()
mata:
mata set matastrict on
real matrix clean_xy(real matrix xy_arg) {
	
	real matrix xy
	real matrix index_asc_order; index_asc_order = .
	real matrix w; w = .
	
	xy = select(xy_arg, xy_arg[2,.])
	minindex(xy[1,.], length(xy[1,.]), index_asc_order, w)
	xy[1,.] = xy[1, index_asc_order]
	xy[2,.] = xy[2, index_asc_order]
	
	return(xy)
}
end
