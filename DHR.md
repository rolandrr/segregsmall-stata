# Method of D'Haultfoeuille and Rathelot (DHR) 

## mata_DHR_ML_oneK
creates mata structures and functions to perform Maximum Likelihood (ML) K by K (oneK)

struct init_estimation {
	real scalar K, L, I_K_odd, nb_obs
	real vector list_Nks_Xm1_positive_obs
}
Construct_init_estimation()
_View_init_estimation()

struct UML_results {
	real colvector P_tilde, m_tilde
	real scalar I_m_tilde_in_M_tol
	real scalar K_max_info_m_tilde
}
Construct_UML_results()
Cons_UML_results_from_dataNks()
Def_P_tilde()
Def_Q()
I_valid_mu()
Test_mu_in_M()
construct_A()
construct_B()
construct_C()
_View_UML_results()
_Get_UML_results()
_Get_mu_tilde()
_Get_m_tilde_from_data()

options_optimization_CML {
	real scalar nb_max_iter
	real scalar ptol, vtol, nrtol
}
Cons_options_optimization_CML()
_View_options_optimization_CML()

struct CML_results {
	real matrix xy_hat /* cf. infra Moment_discrete_distribution() for explanation */
	real colvector P_hat
	real scalar m1_hat /* m1_hat : first moment of the distribution xy_hat */	
}
Construct_CML_results()
_View_CML_results()
Moment_discrete_distribution()
Get_P_hat_CML()

struct raw_ML_results {
	struct init_estimation scalar init_esti
	struct UML_results scalar UML
	struct CML_results scalar CML
	real scalar I_constrained_distribution
}
_View_raw_ML_results()
_Get_m_tilde_raw_ML_results()
_Get_xy_hat_raw_ML_results()
_Get_xy_hat_dataNks()
Construct_raw_ML_results()
Estimation_ML()

Get_xy_hat_CML()
eval_xtiytim1_formulebas_mis_CML()
get_xy_0_sh()
get_xy_0_s0()
get_xtiytim10()
fun_logistic()
fun_logit()
fun_logistic_deriv1()
fun_logistic_deriv2()
clean_xy()


## mata_DHR_ML_allK
creates mata structure and functions to perform ML in the case K random

struct results_prob_K_hat {
	real scalar Kbar, expectation_K_hat, nb_obs_allK
	real colvector probabilities
}
_View_results_prob_K_hat()
_Get_results_prob_K_hat_colv()
Construct_results_prob_K_hat()

struct raw_ML_results_allK {
	real scalar Kbar, nb_obs_allK
	struct raw_ML_results colvector colv_ML
	real scalar I_constrained_distribution_allK
	real colvector index_K_with_nb_obs_positive, m1k_hat
}
_View_I_constrained_allK()
Get_I_constrained_from_ML_allK()
_View_K_with_nb_obs_positive()
_View_m1k_hat()
_Get_m1k_hat()
_Get_m1k_hat_nb_obs_positive()
_View_raw_ML_results_from_allK()
_Get_xy_hat_from_allK()
_Get_m_tilde_oneK()
_Get_m_tilde_from_ML_allK()
Estimation_ML_allK()


## mata_DHR_principal_representation
creates the structure principal_representations (pr) and commpute pr

principal_representations
_View_principal_representations()
Get_principal_representation()
Principal_representations()
get_xy_from_x_and_mu()


## mata_DHR_optibycr_A21
performs (approximate) optimization to obtain the bounds in the case of assumption 2.1 only (Duncan); approximate since it uses a search through feasible canonical representations (cr) only (but ok in practice and far better than Mata optimizer for which the constraints appear too difficult), and some theoretical justifications in progress

int_h_Fp_hats_A21_un_Duncan_cr()


## mata_DHR_bounds
performs bounds estimation for Duncan, Theil, Coworker and Atkinson for random K (and thus also for single fixed K)

compute_nu_Duncan()
compute_nu_Theil()
xlnx()
compute_nu_Atkinson()
compute_nu_Coworker()

compute_int_h_Fp_dis_Duncan()
compute_int_h_Fp_dis_Theil()
compute_int_h_Fp_dis_Atkinson()
compute_int_h_Fp_dis_Coworker()

Def_m01_hat_allK()
Def_expectation_Kp_hat()
Def_expectation_K_hat()

Def_int_h_Fpk_Duncan()
Def_int_h_Fpk_Theil()
Def_int_h_Fpk_Atkinson()
Def_int_h_Fpk_Coworker()

Get_m1hat_oneK()
Duncan_oneK()
Theil_oneK()
Atkinson_oneK()
Coworker_oneK()

options_optimization {
	struct options_optimization_CML scalar CML
	real scalar nb_candidates_cr
}
_View_options_optimization()
Construct_options_optimization()
options_opt_default = Construct_options_optimization(200, 10^(-8), 10^(-9), 10^(-7), 400)

Bounds_DTACW()
Bounds_DTACW_from_db()

Bounds_DTACW_Eq_B1_B2()
Bounds_DTACW_Eq_B1_B2_from_db()


## mata_DHR_ML_Kpind
creates the structure and performs Maximum Likelihood estimation when assumption K and p independent (Kpind)

struct init_estimation_Kpind {
	real colvector list_K_positive_obs
	real scalar nb_K_positive_obs
	real Kbar, Lbar, I_Kbar_odd
	real nb_obs
}
Construct_init_estimation_Kpind()
_View_init_estimation_Kpind()

struct UML_results_Kpind {
	real colvector P_tilde, m_tilde
	real scalar I_m_tilde_in_M_tol
	real scalar K_max_info_m_tilde
}
_View_UML_results_Kpind()
Construct_UML_results_Kpind()
Cons_UML_Kpind_from_UML_perK()
Def_Q_stacked()
Get_P_tilde_stacked()
Get_m_tilde_Kpind_EqD2()

struct CML_results_Kpind {
	real matrix xy_hat
	real colvector P_hat
	real scalar m1_hat /* m1_hat : first moment of the distribution xy_hat */	
}
_View_CML_results_Kpind()
Construct_CML_results_Kpind()

struct raw_ML_results_Kpind {
	struct init_estimation_Kpind scalar init_esti
	struct UML_results_Kpind scalar UML
	struct CML_results_Kpind scalar CML
	real scalar I_constrained_distribution
}
_View_raw_ML_results_Kpind()
Estimation_ML_Kpind()

Get_xy_hat_CML_K_p_independent()
eval_CML_indep_Kp()


## mata_DHR_bounds_Kpind
estimates the bounds when assumption Kpind

Def_int_h_Fpk_Duncan_Kpind()
Def_int_h_Fpk_Theil_Kpind()
Def_int_h_Fpk_Atkinson_Kpind()
Def_int_h_Fpk_Coworker_Kpind()

Bounds_DTACW_Kpind()
Bounds_DTACW_Kpind_from_db()


## mata_DHR_tools_ci
creates structures and functions used to compute confidence intervals in K random unconditional

struct struct_In_Pb {
	real colvector Pb
	real scalar In
}
complete_P()
Vector_first_moments()
Distance_euclidian()
Project_frontier_M() /* use principal representation to get something on frontier of M and only difference in the last moment, not sure actually to do the projection but as K increases, difference between m_in_int_M and m_projected decreases and above all it seems difficult to do better with guarantee and in reasonable amount of time in Mata */
Define_In_Pb_oneK() /* use for oneK only one index - the Theil - for all index used, similar to Matlab code */
_View_struct_In_Pb()

Draw_K_b()
my_rndmultinomial()

struct struct_In_Pb_allK {
	struct struct_In_Pb colvector colv_In_Pb
}
Define_In_Pb_allK()
Define_In_Pb_allK_from_db()
_View_struct_In_Pb_allK()
Draw_db_bootstrap()

options_ci_bootstrap {
	real scalar nb_bootstrap_distance /* for computing V bootstrap de Delta et definir I_n */
	real scalar nb_bootstrap_repetition /* number of bootstrap repetition for inference */
}
Cons_options_ci_bootstrap()

options_optimization_cib {
	struct options_optimization scalar optimization
	struct options_ci_bootstrap scalar cib
}
Cons_options_optimization_cib()

Define_In_ci() /* define I_n for allK aggregated measure, specific to each index, to determine whether boundary or interior type confidence interval is used */

my_empirical_quantile()
CI_interior()
CI_boundary()
CI_1()


## mata_DHR_tools_ci_Kpind
creates structures and functions used to compute confidence intervals in K random, Kpind, unconditional

Get_m1hat_Kpind()
Theil_oneK_Kpind()

struct struct_In_mb_Kpind{
	real scalar In
	real colvector mb
}
_View_struct_In_mb_Kpind()

Draw_db_bootstrap_Kpind()

Define_In_mb()
Define_In_mb_from_db()


## mata_DHR_conditional
defines structures and functions to perform conditional analysis

Struct_results {
	real matrix res
}
Struct_store_bootstrap {
	struct Struct_results rowvector per_type
	real matrix aggregated
}
_View_Results_cond_In_Pb()
_View_Results_cond_In_mb()

Aggregation_across_types()

Bounds_DTACW_from_struct_db_cond()

Draw_boot_store_ZK_nb_units()
Draw_boot_store_WK_nb_units()

Complete_db_bootstrap()
Complete_db_bootstrap_Kpind()

Draw_Struct_db_cond_bootstrap()
Draw_Struct_db_cond_bootstrap_i()


## mata_DHR_tools_test_binomial
defines functions necessary to perform the test of the binomial assumption

Compute_LR_oneK()
Compute_sum_LR_allK()
Extract_data_Nks_from_db()
Compute_sum_LR_allK_from_db()

my_empirical_p_value()

Test_binomial_acrossK()
Test_binomial_acrossK_from_db()


## mata_DHR_stata_output
defines the wrap-up functions that are called from the ado-file defining the final command, both for unconditional and conditional analyses

Get_bounds_for_bootstrap_HR()

Define_nb_bootstrap_distance

St_bounds_ci_np_DTACW_uncond()

St_bounds_ci_np_DTACW_cond()

