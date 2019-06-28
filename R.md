# Method Rathelot (R, or beta, Rbeta) 

## mata_R_ML
does the Maximum Likelihood Estimation of R (for c = 1) only assumption of a Beta distribution for F_p

fun_beta()
fun_gamma_deriv1()
fun_gamma_deriv2()
fun_my_beta_deriv1()
fun_my_beta_deriv2()
fun_my_beta_deriv12()
Def_matrix_betas()
my_prod()

Get_Exp_K_from_db()
Get_Exp_X_from_db()

eval_parametric_ML()

MM_beta_estimate_oneK()
MM_beta_estimate_allK()

Transform_beta_parameters()

beta_ML()


## mata_R_estimates
does the estimation of the indices from the coefficients of beta mixture model estimated in mata_R_ML

Duncan_beta()
Theil_beta()
Atkinson_beta()
Coworker_beta()
/* sub-functions for computation of the Gini */
Approx_int_Fp2_beta_ceq1()
Approx_int_beta_csg1()
Function_f_R_section3()
Function_integrand_equation6()
Gini_beta()

Def_expec_K_hat_from_id_pKh()
Get_db_for_oneK()

Estimates_DTACWG_beta_resuncond()
Estimates_DTACWG_beta_eci()

Estimates_DTACWG_beta_Kpind_runc()
Estimates_DTACWG_beta_Kpind_eci()


## mata_R_tools_ci
does the bootstrap sample for unconditional case, conditional case with unit-level covariates or individual-level covariates. For all, it is a classical bootstrap performed at the level of the units; the difference in the code come from the differences structures of the data structure "struct_db_cond" and it uses only Mata to do the bootstrap, no need to return to the Stata database at unit-level

Draw_db_bootstrap_unit_level()

Draw_struct_db_condunit_unit()

Draw_boot_store_WXK_nb_units()
Draw_struct_db_condindi_unit()

CI_by_percentile_bootstrap()


## mata_R_deltamethod
estimation and confidence interval by delta-method (case Kpind, unconditional)

eval_Duncan_beta()
get_d_Duncan_beta()
eval_Theil_beta()
get_d_Theil_beta()
eval_Atkinson_beta()
get_d_Atkinson_beta()
eval_Coworker_beta()
get_d_Coworker_beta()
eval_Gini_beta()
get_d_Gini_beta()

CI_from_asymptotic_normality()

Estici_dm_DTACWG_beta_Kpind_runc()


## mata_R_conditional
estimation from struct_db_cond for conditional analyses

Esti_DTACWG_beta_from_structcond()


## mata_R_stata_output
wrap-up functions for Stata output

St_esti_ci_beta_DTACWG_uncond()

Get_estimates_for_bootstrap_R()
St_esti_ci_beta_DTACWG_cond()

