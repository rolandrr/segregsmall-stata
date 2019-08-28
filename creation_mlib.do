version 14.2
/* this do-files creates the mata functions used for the package segregsmall 
and save and gather them in a mata library lsegregsmall.mlib */

/******************************************************************************/
/* first launch the do-files that create the mata structures and functions */
/******************************************************************************/
/* NB: the order of the files matter as some structures use previous ones */
/* directory where are the .do that create the mata structures and functions */
clear all
cd "C:/Temp/Dropbox/Segregation_polarisation/Code Stata XDH-RR/Current/beta_version_v1/do_files_for_mata"
//cd "/Users/lucasgirard/Dropbox/Segregation_polarisation/Code Stata XDH-RR/Current/beta_version_v1/do_files_for_mata"

/* data base manipulation and output */
do mata_output
do mata_tools_database
/* dhr */
do mata_DHR_ML_oneK
do mata_DHR_ML_allK
do mata_DHR_principal_representation
do mata_DHR_optibycr_A21
do mata_DHR_bounds
do mata_DHR_ML_Kpind
do mata_DHR_bounds_Kpind
do mata_DHR_tools_ci_Kpind
do mata_DHR_tools_ci
do mata_DHR_conditional
do mata_DHR_tools_test_binomial
do mata_DHR_stata_output
/* rbeta */
do mata_R_ML
do mata_R_estimates
do mata_R_tools_ci
do mata_R_deltamethod
do mata_R_conditional
do mata_R_stata_output
/* ct */
do mata_CT_proportion_estimates
do mata_CT_correction_ct
do mata_CT_conditional
do mata_CT_stata_output


/******************************************************************************/
/* Creation of the mlib */
/******************************************************************************/
/* directory where to store the .mlib files (later to be put online to be shared) */
cd "C:/Temp/Dropbox/Segregation_polarisation/Code Stata XDH-RR/Current/beta_version_v1/files_to_distribute"
//cd "/Users/lucasgirard/Dropbox/Segregation_polarisation/Code Stata XDH-RR/Current/beta_version_v1/files_to_distribute"

/* lsegregsmall_dataoutput (begin) */
mata: mata mlib create lsegregsmall, replace

local names_functions_in_lib ///
Info_data() /// data manipulation and output (begin)
Struct_db_uncond() ///
Struct_db() ///
Struct_db_cond() ///
Info_estimation_inference() ///
Results_uncond() ///
Struct_eci() ///
Results_cond() ///
Results_eci() ///
_View_Info_data() ///
_View_Struct_db_uncond() ///
_View_Struct_db_cond() ///
_View_Info_estimation_inference() ///
_View_Results_uncond() ///
_View_Results_cond() ///
Get_nb_units_studied() ///
Compute_Pi_from_scratch() ///
Compute_S() ///
Compute_S_X() ///
Compute_Pi() ///
Cons_store_ZK_nb_units() ///
Cons_type_frequencies_unit() ///
Cons_type_probabilities_unit() ///
Cons_store_WK_nb_units() ///
Cons_nb_units_per_type_ind() ///
Cons_type_frequencies_ind() ///
Cons_type_probabilities_ind() ///
Cons_summary_info_data_per_type() ///
Corr_info_data_uncond_condind() ///
Cons_nb_singleton_cells_all_type() /// data manipulation and output (end)
init_estimation()  /// DHR (begin)
Construct_init_estimation() ///
_View_init_estimation() ///
UML_results() ///
Construct_UML_results() ///
Cons_UML_results_from_dataNks() ///
Def_P_tilde() ///
Def_Q() ///
I_valid_mu() ///
Test_mu_in_M() ///
construct_A() ///
construct_B() ///
construct_C() ///
_View_UML_results() ///
_Get_UML_results() ///
_Get_mu_tilde() ///
_Get_m_tilde_from_data() ///
options_optimization_CML() ///
Cons_options_optimization_CML() ///
_View_options_optimization_CML() ///
CML_results() ///
Construct_CML_results() ///
_View_CML_results() ///
Moment_discrete_distribution() ///
Get_P_hat_CML() ///
raw_ML_results() ///
_View_raw_ML_results() ///
_Get_m_tilde_raw_ML_results() ///
_Get_xy_hat_raw_ML_results() ///
_Get_xy_hat_dataNks() ///
Construct_raw_ML_results() ///
Estimation_ML() ///
Get_xy_hat_CML() ///
eval_xtiytim1_formulebas_mis_CML() ///
get_xy_0_sh() ///
get_xy_0_s0() ///
get_xtiytim10() ///
fun_logistic() ///
fun_logit() ///
fun_logistic_deriv1() ///
fun_logistic_deriv2() ///
clean_xy() ///
results_prob_K_hat() ///
_View_results_prob_K_hat() ///
_Get_results_prob_K_hat_colv() ///
Construct_results_prob_K_hat() ///
raw_ML_results_allK() ///
_View_I_constrained_allK() ///
Get_I_constrained_from_ML_allK() ///
_View_K_with_nb_obs_positive() ///
_View_m1k_hat() ///
_Get_m1k_hat() ///
_Get_m1k_hat_nb_obs_positive() ///
_View_raw_ML_results_from_allK() ///
_Get_xy_hat_from_allK() ///
_Get_m_tilde_oneK() ///
_Get_m_tilde_from_ML_allK() ///
Estimation_ML_allK() ///
principal_representations() ///
_View_principal_representations() ///
Get_principal_representation() ///
Principal_representations() ///
get_xy_from_x_and_mu() ///
int_h_Fp_hats_A21_un_Duncan_cr() ///
compute_nu_Duncan() ///
compute_nu_Theil() ///
xlnx() ///
compute_nu_Atkinson() ///
compute_nu_Coworker() ///
compute_int_h_Fp_dis_Duncan() ///
compute_int_h_Fp_dis_Theil() ///
compute_int_h_Fp_dis_Atkinson() ///
compute_int_h_Fp_dis_Coworker() ///
Def_m01_hat_allK() ///
Def_expectation_Kp_hat() ///
Def_expectation_K_hat() ///
Def_int_h_Fpk_Duncan() ///
Def_int_h_Fpk_Theil() ///
Def_int_h_Fpk_Atkinson() ///
Def_int_h_Fpk_Coworker() ///
Get_m1hat_oneK() ///
Duncan_oneK() ///
Theil_oneK() ///
Atkinson_oneK() ///
Coworker_oneK() ///
options_optimization() ///
_View_options_optimization() ///
Construct_options_optimization() ///
Bounds_DTACW() ///
Bounds_DTACW_from_db() ///
init_estimation_Kpind() ///
Construct_init_estimation_Kpind() ///
_View_init_estimation_Kpind() ///
UML_results_Kpind() ///
_View_UML_results_Kpind() ///
Construct_UML_results_Kpind() ///
Cons_UML_Kpind_from_UML_perK() ///
Def_Q_stacked() ///
Get_P_tilde_stacked() ///
Get_m_tilde_Kpind_EqD2() ///
CML_results_Kpind() ///
_View_CML_results_Kpind() ///
Construct_CML_results_Kpind() ///
raw_ML_results_Kpind() ///
_View_raw_ML_results_Kpind() ///
Estimation_ML_Kpind() ///
Get_xy_hat_CML_K_p_independent() ///
eval_CML_indep_Kp() ///
Def_int_h_Fpk_Duncan_Kpind() ///
Def_int_h_Fpk_Theil_Kpind() ///
Def_int_h_Fpk_Atkinson_Kpind() ///
Def_int_h_Fpk_Coworker_Kpind() ///
Bounds_DTACW_Kpind() ///
Bounds_DTACW_Kpind_from_db() ///
struct_In_Pb() ///
complete_P() ///
Vector_first_moments() ///
Distance_euclidian() ///
Project_frontier_M() ///
Define_In_Pb_oneK() ///
_View_struct_In_Pb() ///
Draw_K_b() ///
my_rndmultinomial() ///
struct_In_Pb_allK() ///
Define_In_Pb_allK() ///
_View_struct_In_Pb_allK() ///
Draw_db_bootstrap() ///
options_ci_bootstrap() ///
Cons_options_ci_bootstrap() ///
options_optimization_cib() ///
Cons_options_optimization_cib() ///
Define_In_ci() ///
my_empirical_quantile() ///
CI_interior() ///
CI_boundary() ///
CI_1() ///
Get_m1hat_Kpind() ///
Theil_oneK_Kpind() ///
struct_In_mb_Kpind() ///
_View_struct_In_mb_Kpind() ///
Draw_db_bootstrap_Kpind() ///
Define_In_mb() ///
Struct_results() ///
Struct_store_bootstrap() ///
_View_Results_cond_In_Pb() ///
_View_Results_cond_In_mb() ///
Aggregation_across_types() ///
Bounds_DTACW_from_struct_db_cond() ///
Draw_boot_store_ZK_nb_units() ///
Draw_boot_store_WK_nb_units() ///
Complete_db_bootstrap() ///
Complete_db_bootstrap_Kpind() ///
Draw_Struct_db_cond_bootstrap() ///
Draw_Struct_db_cond_bootstrap_i() ///
Compute_LR_oneK() ///
Compute_sum_LR_allK() ///
Extract_data_Nks_from_db() ///
Compute_sum_LR_allK_from_db() ///
my_empirical_p_value() ///
Test_binomial_acrossK() ///
Test_binomial_acrossK_from_db() ///
Get_bounds_for_bootstrap_HR() ///
Define_nb_bootstrap_distance() ///
St_bounds_ci_np_DTACW_uncond() ///
St_bounds_ci_np_DTACW_cond() /// DHR (end)
fun_beta() /// Rbeta (begin)
fun_gamma_deriv1() ///
fun_gamma_deriv2() ///
fun_my_beta_deriv1() ///
fun_my_beta_deriv2() ///
fun_my_beta_deriv12() ///
Def_matrix_betas() ///
my_prod() ///
Get_Exp_K_from_db() ///
Get_Exp_X_from_db() ///
eval_parametric_ML() ///
MM_beta_estimate_oneK() ///
MM_beta_estimate_allK() ///
Transform_beta_parameters() ///
beta_ML() ///
Duncan_beta() ///
Theil_beta() ///
Atkinson_beta() ///
Coworker_beta() ///
Approx_int_Fp2_beta_ceq1() ///
Approx_int_beta_csg1() ///
Function_f_R_section3() ///
Function_integrand_equation6() ///
Gini_beta() ///
Def_expec_K_hat_from_id_pKh() ///
Get_db_for_oneK() ///
Estimates_DTACWG_beta_resuncond() ///
Estimates_DTACWG_beta_eci() ///
Estimates_DTACWG_beta_Kpind_runc() ///
Estimates_DTACWG_beta_Kpind_eci() ///
Draw_db_bootstrap_unit_level() ///
Draw_struct_db_condunit_unit() ///
Draw_boot_store_WXK_nb_units() ///
Draw_struct_db_condindi_unit() ///
CI_by_percentile_bootstrap() ///
eval_Duncan_beta() ///
get_d_Duncan_beta() ///
eval_Theil_beta() ///
get_d_Theil_beta() ///
eval_Atkinson_beta() ///
get_d_Atkinson_beta() ///
eval_Coworker_beta() ///
get_d_Coworker_beta() ///
eval_Gini_beta() ///
get_d_Gini_beta() ///
CI_from_asymptotic_normality() ///
Estici_dm_DTACWG_beta_Kpind_runc() ///
Esti_DTACWG_beta_from_structcond() ///
St_esti_ci_beta_DTACWG_uncond() ///
Get_estimates_for_bootstrap_R() ///
St_esti_ci_beta_DTACWG_cond() /// Rbeta (end)
Prop_Duncan() /// CT (begin)
Prop_Duncan_from_struct_db_unc() ///
compute_xlog2x() ///
Prop_Theil() ///
Prop_Theil_from_struct_db_unc() ///
Prop_Atkinson() ///
Prop_Atkinson_from_struct_db_unc() ///
Prop_Coworker() ///
Prop_Coworker_from_struct_db_unc() ///
Prop_Gini() ///
Prop_Gini_from_struct_db_unc() ///
Construct_db_random_allocation() ///
count_rbinomial() ///
Compute_CT_correction() ///
Compute_CFC_standard_score() ///
Esti_ct_DTACWG_for_one_type() ///
Aggregation_across_types_CT() ///
St_esti_ct_DTACWG_uncond() ///
_Get_ct_indices() ///
St_esti_ct_DTACWG_cond() /* CT (end) */

/* add each function in the library (begin) */
foreach name_function in `names_functions_in_lib' {
	di "Function to be added: `name_function'"
	mata: mata mlib add  lsegregsmall `name_function'
}
/* add each function in the library (end) */	
