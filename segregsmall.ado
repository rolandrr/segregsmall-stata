capture program drop segregsmall
program define segregsmall, eclass

	version 13
	
	/**************************/
	/* Syntax command (begin) */
	syntax varlist(numeric) [if] [in], ///
		Format(string) Method(string) /// format of database and method
		[atkinson(real 0.5)] /// parameter for the Atkinson index
		[WITHsingle] [EXCLUdingsinglepertype] /// treatment of singletons
		[INDEPendencekp] /// estimation and inference when K and p are supposed independent
		[REPBootstrap(integer 200)] [Level(real 0)] [NOCI] /// options for inference
		[repct(integer 50)] /// number of repetitions for CT (Carrington-Troske correction)
		[TESTBinomial] /// test of binomial assumption
		[CONDItional(string)] // for conditional analysis
	/* Syntax command (end) */
	/************************/
	
	/*************************/
	/* Option checks (begin) */
	/* format */
	if (("`format'" != "unit") & ("`format'" != "indiv")){
		display as error "Option format incorrectly specified: it must be either unit or indiv"
		exit 198 
	}

	/* method */
	if (("`method'" != "np") & ("`method'" != "beta") & ("`method'" != "ct")){
		display as error "Option method incorrectly specified: it must be one of np, beta, ct"
		exit 198
	}

	/* atkinson */
	if ((`atkinson' <= 0) | (`atkinson' >= 1)){
			display as error "Option atkinson incorrectly specified: it must be a real in (0,1)"
			exit 198
	}

	/* conditional */
	if ("`conditional'" == "") local I_conditional 0
	else local I_conditional 1
	if (`I_conditional') {
		if (("`conditional'" != "unit") & ("`conditional'" != "indiv")){
			display as error "Option conditional incorrectly specified: it must be either unit or indiv"
			exit 198
		}
		if (("`conditional'" == "indiv") & ("`format'" != "indiv")){
			display as error "Option conditional(indiv) requires format(indiv) databases"
			exit 184
		}
	}
	
	/* independencekp */
	if ("`independencekp'" == "") local I_independencekp 0
	else local I_independencekp 1

	/* testbinomial */
	if ("`testbinomial'" == "") local I_testbinomial 0
	else local I_testbinomial 1
	if ((`I_conditional') & (`I_testbinomial')){
		display as error "Options conditional and testbinomial cannot be combined"
		exit 184
	} 

	/* withsingle and excludingsinglepertype */
	if ("`withsingle'" == "") local I_withsingle 0
	else local I_withsingle 1
	if ("`excludingsinglepertype'" == "") local I_excludingsinglepertype 0
	else local I_excludingsinglepertype 1	
	if ((`I_withsingle') & (`I_excludingsinglepertype')){
		display as error "Options withsingle and excludingsinglepertype cannot be combined"
		exit 184
	}
	if (((`I_excludingsinglepertype') & (!`I_conditional')) | ((`I_excludingsinglepertype') & ("`conditional'" == "unit"))){
		display as error "Option excludingsinglepertype is only available with option conditional(indiv)"
		exit 184
	}

	/* repbootstrap */
	if (`repbootstrap' <= 1){
		display as error "Option repbootstrap incorrectly specified: it must be a positive integer, > 1"
		exit 198
	}

	/* repct */
	if (`repct' <= 1){
		display as error "Option repct incorrectly specified: it must be a positive integer, > 1"
		exit 198
	}

	/* level */
	if (`level' == 0) local I_specified_alpha 0
	else local I_specified_alpha 1
	if ((`I_specified_alpha') & ((`level' < 0) | (`level' >= 1))){
			display as error "Option level incorrectly specified: it must be a real scalar in (0,1)"
			exit 198
	}
	local alpha = 1 - `level'
	if (`I_specified_alpha') local alpha_text = round(100*(1-`alpha'), 0.1)
	else local alpha_text 95
	
	/* noinference */		
	if ("`noci'" == "") local I_noinference 0
	else local I_noinference 1		
	if ((`I_noinference') & (`I_specified_alpha')){
		display as error "Option level and noci cannot be combined"
		exit 184
	} 

	/* irrelevant options with method(ct) */
	if ("`method'" == "ct"){
		if(`I_noinference'){
			display as error "Option noci is irrelevant with method(ct)"
			exit 184
		}
		if(`I_independencekp'){
			display as error "Option independencekp is irrelevant with method(ct)"
			exit 184
		}
	}
	/* Option checks (end) */
	/***********************/

	/******************************************************/
	/* Selecton of observations [if] [in] options (begin) */
	preserve
	marksample touse, novarlist 
	/* novarlist option: no automatic drop of missing values, cf. checks in 
	ado-files db_from_* */
	quietly: keep if `touse'
	/* Selecton of observations [if] [in] options (end) */
	/****************************************************/

	/************************************/
	/* estimation and inference (begin) */
	display as text "" /* to leave some space, but less that display _newline */
	display as text "*** Construction of relevant databases for the analysis ***"		
	tempname results

	if (`I_conditional'){ /* Conditional analysis (begin) */

		/* import data (begin) */
		if ("`conditional'" == "unit") { /* unit-level characteristic (begin) */
			if (`I_withsingle') db_from_`format'_level_cond_unit `varlist', withsingleton
			else db_from_`format'_level_cond_unit `varlist'
		} /* unit-level characteristic (end) */
		else { /* indiv-level characteristic (begin) */
	
			if ((!`I_withsingle') & (!`I_excludingsinglepertype')) 		db_from_indiv_level_cond_indiv `varlist'
			else if ((`I_withsingle') & (!`I_excludingsinglepertype')) 	db_from_indiv_level_cond_indiv `varlist', withsingleton
			else 														db_from_indiv_level_cond_indiv `varlist', excludingsingletonpertype	
		} /* indiv-level characteristic (end) */
		/* import data (end) */
		
		if (`I_noinference') local repbootstrap_for_cond 0
		else local repbootstrap_for_cond = `repbootstrap'
		
		if ("`method'" == "np"){ /* Method np (begin) */
			if(`I_specified_alpha'){
				mata: `results' = St_bounds_ci_np_DTACW_cond(struct_db_cond, `atkinson', `I_independencekp', `repbootstrap_for_cond', `alpha')
			}
			else {
				mata: `results' = St_bounds_ci_np_DTACW_cond(struct_db_cond, `atkinson', `I_independencekp', `repbootstrap_for_cond')
			}
		} /* Method np (end) */
		
		else if ("`method'" == "beta"){ /* Method beta (begin) */
			if(`I_specified_alpha'){
				mata: `results' = St_esti_ci_beta_DTACWG_cond(struct_db_cond, `atkinson', `I_independencekp', `repbootstrap_for_cond', `alpha')
			}
			else {
				mata: `results' = St_esti_ci_beta_DTACWG_cond(struct_db_cond, `atkinson', `I_independencekp', `repbootstrap_for_cond')
			}
		} /* Method beta (end) */
		
		else { /* Method ct (begin) */
			if(`I_specified_alpha'){
				mata: `results' = St_esti_ct_DTACWG_cond(struct_db_cond, `atkinson', `repct', `alpha')
			}
			else {
				mata: `results' = St_esti_ct_DTACWG_cond(struct_db_cond, `atkinson', `repct')
			}			
		} /* Method ct (end) */
		
	} /* Conditional analysis (end) */

	else { /* Unconditional analysis (begin) */

		/* import data (begin) */
		if (`I_withsingle') db_from_`format'_level `varlist', withsingleton
		else db_from_`format'_level `varlist'		 
		/* import data (end) */

		if("`method'" == "np"){ /* Method np (begin) */
			if(`I_specified_alpha'){
				mata: `results' = St_bounds_ci_np_DTACW_uncond(struct_db_uncond, `atkinson', `I_independencekp', `repbootstrap', `I_noinference', `I_testbinomial', `alpha')
			}
			else {
				mata: `results' = St_bounds_ci_np_DTACW_uncond(struct_db_uncond, `atkinson', `I_independencekp', `repbootstrap', `I_noinference', `I_testbinomial')
			}
		} /* Method np (end) */
		
		else if ("`method'" == "beta"){ /* Method beta (begin) */
			if(`I_specified_alpha'){
				mata: `results' = St_esti_ci_beta_DTACWG_uncond(struct_db_uncond, `atkinson', `I_independencekp', `repbootstrap', `I_noinference', 0, `I_testbinomial', `alpha')
			}
			else {
				mata: `results' = St_esti_ci_beta_DTACWG_uncond(struct_db_uncond, `atkinson', `I_independencekp', `repbootstrap', `I_noinference', 0, `I_testbinomial')
			}
		} /* Method beta (end) */
		
		else { /* Method ct (begin) */
			if(`I_specified_alpha'){
				mata: `results' = St_esti_ct_DTACWG_uncond(struct_db_uncond, `atkinson', `repct', `repbootstrap', `I_testbinomial', `alpha')
			}
			else {
				mata: `results' = St_esti_ct_DTACWG_uncond(struct_db_uncond, `atkinson', `repct', `repbootstrap', `I_testbinomial')
			}
		} /* Method ct (end) */

	} /* Unconditional analysis (end) */
	/* estimation and inference (end) */
	/**********************************/

	/******************/
	/* Output (begin) */
	Output_`method'_method 	`I_conditional' `results' `I_specified_alpha' `alpha_text' `I_noinference' `I_independencekp' `I_testbinomial' `repbootstrap' `I_withsingle' `I_excludingsinglepertype' `repct'
	ereturn local cmd_arguments "`0'"
	ereturn local cmd "segregsmall"
	/* Output (end) */
	/****************/

end


/******************************************************************************/
/* LOCAL SUBROUTINES FOR OUTPUT */
/******************************************************************************/

/******************************************************************************/
/* Sub-programs for the command measureseg */
/******************************************************************************/

/******************************************************************************/
/* Output_info_data_uncond */
/* transfer into Stata objects (for return later) general information about data 
in unconditional analyses */
capture program drop Output_info_data_uncond
program define Output_info_data_uncond, eclass

	mata: st_numscalar("I_withsingle", struct_db_uncond.info_data.I_withsingleton)
	mata: st_numscalar("nb_units_total", struct_db_uncond.info_data.nb_units_total)
	mata: st_numscalar("nb_units_singleton", struct_db_uncond.info_data.nb_units_singleton)
	mata: st_numscalar("nb_units_studied", struct_db_uncond.info_data.nb_units_studied)
	mata: st_numscalar("nb_individuals", struct_db_uncond.info_data.nb_individuals)
	mata: st_numscalar("nb_minority_individuals", struct_db_uncond.info_data.nb_minority_individuals)
	mata: st_numscalar("prop_minority_hat", struct_db_uncond.info_data.prop_minority_hat)
	mata: st_numscalar("Kbar", struct_db_uncond.info_data.Kbar)
	mata: st_matrix("list_K_positive_obs", struct_db_uncond.info_data.list_K_positive_obs)
	matrix colnames list_K_positive_obs = K
	mata: st_numscalar("nb_K_positive_obs", struct_db_uncond.info_data.nb_K_positive_obs)
	
	ereturn scalar I_withsingle = I_withsingle
	ereturn scalar nb_units_total = nb_units_total
	ereturn scalar nb_units_single = nb_units_singleton
	ereturn scalar nb_units_studied = nb_units_studied
	ereturn scalar nb_individuals = nb_individuals
	ereturn scalar nb_minority_individuals = nb_minority_individuals
	ereturn scalar prop_minority_hat = prop_minority_hat
	ereturn scalar K_max = Kbar
	ereturn matrix list_K_with_obs = list_K_positive_obs
	ereturn scalar nb_K_with_obs = nb_K_positive_obs
	
end	
/******************************************************************************/

/******************************************************************************/
/* Output_info_data_cond */
/* transfer into Stata objects (for return later) general information about data 
in conditional analyses */
capture program drop Output_info_data_cond
program define Output_info_data_cond, eclass

	/* aggregated information (begin) */
	mata: st_numscalar("nb_types", struct_db_cond.nb_types)
	mata: st_numscalar("I_unit_level_characteristic", struct_db_cond.I_unit_level_characteristic)
	mata: st_numscalar("I_excludingsinglepertype", struct_db_cond.I_excludingsingletonpertype)
	mata: st_matrix("type_frequencies", struct_db_cond.type_frequencies)
	mata: st_matrix("type_probabilities", struct_db_cond.type_probabilities)
	mata: st_matrix("nb_units_studied_per_type", struct_db_cond.nb_units_studied_per_type)
	mata: st_numscalar("I_withsingle", struct_db_cond.info_data_uncond.I_withsingleton)
	mata: st_numscalar("nb_units_total", struct_db_cond.info_data_uncond.nb_units_total)
	mata: st_numscalar("nb_units_singleton", struct_db_cond.info_data_uncond.nb_units_singleton)
	mata: st_numscalar("nb_units_studied", struct_db_cond.info_data_uncond.nb_units_studied)
	mata: st_numscalar("nb_individuals", struct_db_cond.info_data_uncond.nb_individuals)
	mata: st_numscalar("nb_minority_individuals", struct_db_cond.info_data_uncond.nb_minority_individuals)
	mata: st_numscalar("prop_minority_hat", struct_db_cond.info_data_uncond.prop_minority_hat)
	mata: st_numscalar("Kbar", struct_db_cond.info_data_uncond.Kbar)
	mata: st_matrix("list_K_positive_obs", struct_db_cond.info_data_uncond.list_K_positive_obs)
	matrix colnames list_K_positive_obs = K
	mata: st_numscalar("nb_K_positive_obs", struct_db_cond.info_data_uncond.nb_K_positive_obs)	
	/* aggregated information (end) */
	
	/* information type by type (begin) */
	mata: st_matrix("summary_info_data_per_type", struct_db_cond.summary_info_data_per_type)
	
	local rownames_type 
	forvalues type = 1/`=nb_types' {
		local rownames_type `rownames_type' type_`type'
	}
	
	matrix rownames type_frequencies = `rownames_type'
	matrix rownames type_probabilities = `rownames_type'
	matrix rownames nb_units_studied_per_type = `rownames_type'
	matrix rownames summary_info_data_per_type = `rownames_type'
	
	if (I_unit_level_characteristic){ /* unit-level characteristic (begin) */
		matrix colnames type_frequencies = nb_s_unit
		matrix colnames type_probabilities = prop_s_unit
		matrix colnames nb_units_studied_per_type = nb_s_unit
		matrix colnames summary_info_data_per_type = I_w_single nb_t_unit nb_t_single nb_s_unit nb_ind nb_mino_ind prop_mino K_max nb_K_w_obs
		ereturn matrix nb_units_studied_per_type = nb_units_studied_per_type
	} /* unit-level characteristic (end) */
	else { /* individual-level characteristic (begin) */
		matrix colnames type_frequencies = nb_ind
		matrix colnames type_probabilities = prop_ind
		matrix colnames nb_units_studied_per_type = nb_s_cells
		matrix colnames summary_info_data_per_type = I_w_singlec nb_t_cell nb_t_singlec nb_s_cell nb_ind nb_mino_ind prop_mino K_max nb_K_w_obs
		ereturn matrix nb_cells_studied_per_type = nb_units_studied_per_type
		mata: st_numscalar("nb_cells_studied_sum_across_type", sum(struct_db_cond.nb_units_studied_per_type))
		ereturn scalar nb_cells_studied_sum_across_type = nb_cells_studied_sum_across_type
		mata: st_numscalar("nb_single_cells_sum_across_type", sum(struct_db_cond.nb_singleton_cells_all_type))
		ereturn scalar nb_single_cells_sum_across_type = nb_single_cells_sum_across_type
	} /* individual-level characteristic (end) */
	/* information type by type (end) */
	
	ereturn scalar nb_types = nb_types
	ereturn scalar I_unit_level_characteristic = I_unit_level_characteristic
	ereturn scalar I_withsingle = I_withsingle
	ereturn scalar I_excludingsinglepertype = I_excludingsinglepertype
	ereturn matrix type_frequencies = type_frequencies
	ereturn matrix type_probabilities = type_probabilities
	ereturn scalar nb_units_total = nb_units_total
	ereturn scalar nb_units_single = nb_units_singleton
	ereturn scalar nb_units_studied = nb_units_studied
	ereturn scalar nb_individuals = nb_individuals
	ereturn scalar nb_minority_individuals = nb_minority_individuals
	ereturn scalar prop_minority_hat = prop_minority_hat
	ereturn scalar K_max = Kbar
	ereturn matrix list_K_with_obs = list_K_positive_obs
	ereturn scalar nb_K_with_obs = nb_K_positive_obs
	ereturn matrix summary_info_data_per_type = summary_info_data_per_type
	
end	
/******************************************************************************/

/******************************************************************************/
/* Output_info_estimation_inference */
/* transfer into Stata objects (for return later) user's choice as regards
estimation and inference */
capture program drop Output_info_estimation_inference
program define Output_info_estimation_inference, eclass
	args results
	mata: st_numscalar("I_method_np", `results'.info_eci.I_method_np)
	mata: st_numscalar("I_method_beta", `results'.info_eci.I_method_beta)
	mata: st_numscalar("I_method_ct", `results'.info_eci.I_method_ct)
	mata: st_numscalar("I_conditional", `results'.info_eci.I_conditional)
	mata: st_numscalar("I_unit_level_characteristic", `results'.info_eci.I_unit_level_characteristic)
	mata: st_numscalar("b_atkinson", `results'.info_eci.b_atkinson)
	mata: st_numscalar("I_noinference", `results'.info_eci.I_noinference)
	mata: st_numscalar("nb_bootstrap_repetition", `results'.info_eci.nb_bootstrap_repetition)
	mata: st_numscalar("specified_alpha", `results'.info_eci.specified_alpha)
	mata: st_numscalar("I_hyp_independenceKp", `results'.info_eci.I_hyp_independenceKp)
	mata: st_numscalar("I_testbinomial", `results'.info_eci.I_testbinomial)
	mata: st_numscalar("nb_ct_repetition", `results'.info_eci.nb_ct_repetition)
	
	ereturn scalar I_method_np = I_method_np
	ereturn scalar I_method_beta = I_method_beta
	ereturn scalar I_method_ct = I_method_ct
	ereturn scalar I_conditional = I_conditional
	ereturn scalar I_unit_level_characteristic = I_unit_level_characteristic
	ereturn scalar b_atkinson = b_atkinson
	ereturn scalar I_noci = I_noinference
	ereturn scalar nb_bootstrap_repetition = nb_bootstrap_repetition
	ereturn scalar specified_level = 1 - specified_alpha
	ereturn scalar I_hyp_independenceKp = I_hyp_independenceKp
	ereturn scalar I_testbinomial = I_testbinomial
	ereturn scalar nb_ct_repetition = nb_ct_repetition
	
end	
/******************************************************************************/

/******************************************************************************/
capture program drop Display_test_binomial
program define Display_test_binomial
	args results repbootstrap
	
	mata: st_numscalar("stat_test", `results'.uncond.test_binomial_results[1])
	mata: st_numscalar("p_value", `results'.uncond.test_binomial_results[2])
		
	display as text "" /* to leave some space, but less that display _newline */
	display as text "Test of binomial assumption (H0: conditional binomial distribution):"
	display as text "(distribution under the null obtained by bootstrap, " as input "`repbootstrap' repetitions" as text ")"

	local col_indent_first_col 12
	local col_vert_line 30
	local col_stat_test 35
	local col_p_value 65
	
	display as text _col(`col_indent_first_col') "Result" ///
		_col(`col_vert_line') "{c |}" ///
		_col(`col_stat_test') "value of test statistic" ///
		_col(`col_p_value') "p-value"
	display as text "{hline `=`col_vert_line'-1'}{c +}{hline 43}"
	
	display as text _col(`col_vert_line') "{c |}" ///
			as text _col(`=`col_stat_test'+7') as result stat_test ///
			as text _col(`=`col_p_value'+2') as result p_value
	
end
/******************************************************************************/

/******************************************************************************/
capture program drop Display_information_data
program define Display_information_data

	args I_conditional I_noinference I_independencekp repbootstrap I_withsingle I_excludingsinglepertype I_method_np_or_beta repct
	
	if (`I_conditional') { /* conditional analysis (begin) */
			if (I_unit_level_characteristic) { /* unit-level characteristics (begin) */
				display as text "Conditional analysis with unit-level covariates (" as result nb_types as text " distinct types)"
				display as text "Number of units studied in the analysis (sum across types): " as result nb_units_studied
				if (nb_units_singleton > 1) local accord_unit units
				else local accord_unit unit 
				if (`I_withsingle') 	display as text "(" as result nb_units_singleton as text " `accord_unit' with a single individual are" as input " included " as text "in the analysis)"
				else 					display as text "(" as result nb_units_singleton as text " `accord_unit' with a single individual are" as input " excluded " as text "from the analysis)"	
			} /* unit-level characteristics (end) */
			
			else { /* individual-level characteristics (begin) */
				display as text "Conditional analysis with individual-level covariates (" as result nb_types as text " distinct types)"
				display as text "Number of units studied in the analysis (sum across types): " as result nb_units_studied
				if (nb_units_singleton > 1) local accord_unit units
				else local accord_unit unit 
				if (`I_withsingle') 	display as text "(" as result nb_units_singleton as text " `accord_unit' with a single individual are" as input " included " as text "in the analysis)"
				else 					display as text "(" as result nb_units_singleton as text " `accord_unit' with a single individual are" as input " excluded " as text "from the analysis)"	
				display as text "It results (sum across types) in " as result nb_cells_studied_sum_across_type as text " cells unit x individual-type studied in the analysis"		
				if (nb_single_cells_sum_across_type > 1) local accord_cell cells
				else local accord_cell cell				
				if (`I_excludingsinglepertype') display as text "(" as result nb_single_cells_sum_across_type as text " `accord_cell' unit x individual-type with a single individual are" as input " excluded " as text "from the analysis)"	
				else 							display as text "(" as result nb_single_cells_sum_across_type as text " `accord_cell' unit x individual-type with a single individual are" as input " included " as text "in the analysis)"
						
			} /* individual-level characteristics (end) */
			
		} /* conditional analysis (end) */
		
		else { /* unconditional analysis (begin) */
			display as text "Unconditional analysis"
			display as text "Number of units studied in the analysis: " as result nb_units_studied
			if (nb_units_singleton > 1) local accord_unit units
			else local accord_unit unit 
			if (`I_withsingle') 	display as text "(" as result nb_units_singleton as text " `accord_unit' with a single individual are" as input " included " as text "in the analysis)"
			else 					display as text "(" as result nb_units_singleton as text " `accord_unit' with a single individual are" as input " excluded " as text "from the analysis)"	
		} /* unconditional analysis (end) */

		display as text "Number of individuals studied: " as result nb_individuals
		display as text "Proportion of minority (or reference) group: " as result %-4.0g prop_minority_hat

		if (`I_method_np_or_beta'){
			if (`I_independencekp') display as text "Assumption on dependence between K and p for estimation and inference: " as input "independence"
			else 					display as text "Assumption on dependence between K and p for estimation and inference: " as input "none"
			
			if (`I_noinference') 	display as text "Inference: " as input "none, only estimation"
			else 					display as text "Inference: " as input " by bootstrap, `repbootstrap' repetitions"
		}
		else {
			display as text "No inference in the proportion-based indices framework"
			display as text "CT-correction is made using " as input "`repct' draws " as text "under random allocation"
		}
		
		display as text "" /* to leave some space, but less than display _newline */
		
		if (`I_conditional') { /* conditional analysis (begin) */
			display as text "Conditional aggregated segregation indices:"
			display as text "(per type indices are available in e(estimates_ci_type_#) matrices)"
		} /* conditional analysis (end) */
		else { /* unconditional analysis (begin) */		
			display as text "Unconditional segregation indices:"
		} /* unconditional analysis (end) */

end
/******************************************************************************/

/******************************************************************************/
capture program drop Output_np_method
program define Output_np_method, eclass

	args I_conditional results I_specified_alpha alpha_text I_noinference I_independencekp I_testbinomial repbootstrap I_withsingle I_excludingsinglepertype repct

	ereturn clear
	Output_info_estimation_inference `results' /* save and ereturn user's choices for estimation and inference */
	
	if (`I_conditional') { /* conditional case (begin) */
	
		Output_info_data_cond /* general information about data */
		
		/* I_constrained_case aggregated (begin) */
		mata: st_numscalar("I_constrained_case", `results'.cond.I_all_type_all_Fpk_constrained) 
		ereturn scalar I_constrained_case = I_constrained_case
		/* I_constrained_case aggregated (end) */
		
		/* I_constrained_case per type (begin) */
		mata: st_matrix("I_constrained_case_per_type", `results'.cond.I_all_Fpk_constrained_per_type)
		local rownames_type 
		forvalues type = 1/`=nb_types' {
			local rownames_type `rownames_type' type_`type'
		}
		matrix rownames I_constrained_case_per_type = `rownames_type'
		matrix colnames I_constrained_case_per_type = I_cons_case
		ereturn matrix I_constrained_case_per_type = I_constrained_case_per_type
		/* I_constrained_case per type (end) */
		
		/* matrix of estimates of the bounds and confidence intervals aggregated (begin) */
		mata: st_matrix("estimates_ci_aggregated", `results'.cond.estimates_ci_aggregated)
		matrix rownames estimates_ci_aggregated = Duncan_unit Duncan_ind Theil_unit Theil_ind Atkinson_unit Atkinson_ind Coworker_unit Coworker_ind
		if (`I_noinference') { /* no inference case (begin) */
			matrix colnames estimates_ci_aggregated = index I_indiv_weight bounds_hat:low bounds_hat:up
		} /* no inference case (end) */		
		else { /* inference case (begin) */
			if (`I_specified_alpha') { /* specified alpha case (begin) */
				matrix colnames estimates_ci_aggregated = index I_indiv_weight bounds_hat:low bounds_hat:up ///
					I_ciboundary ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up ci`alpha_text':low ci`alpha_text':up	
				local col_ci_l 12
				local col_ci_u 13
			} /* specified alpha case (end) */
			else { /* no specified alpha case (begin) */
				matrix colnames estimates_ci_aggregated = index I_indiv_weight bounds_hat:low bounds_hat:up ///
					I_ciboundary ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up
				local col_ci_l 8
				local col_ci_u 9	
			} /* no specified alpha case (end) */
		} /* inference case (end) */
		ereturn matrix estimates_ci_aggregated = estimates_ci_aggregated
		/* matrix of estimates of the bounds and confidence intervals aggregated (end) */
		
		/* elements of the matrix that will be displayed in Stata log (begin) */
		forvalues index = 1/8 {
			mata: st_numscalar("bounds_l_`index'", `results'.cond.estimates_ci_aggregated[`index',3])
			mata: st_numscalar("bounds_u_`index'", `results'.cond.estimates_ci_aggregated[`index',4])
			if (`I_noinference') {
				mata: st_numscalar("ci_l_`index'", .)
				mata: st_numscalar("ci_u_`index'", .)
			}
			else {
				mata: st_numscalar("ci_l_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_ci_l'])
				mata: st_numscalar("ci_u_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_ci_u'])
			}
		}	
		/* elements of the matrix that will be displayed in Stata log (end) */
	
		/* matrix infodistributionofp_np: not available with conditional analyses */
		
		/* matrix test_binomial_results: not available with conditional analyses */
		
		/* matrices of estimates of the bounds and confidence interval per type (begin) */
		forvalues type = 1/`=nb_types' { /* loop over types (begin) */
			mata: st_matrix("estimates_ci_type_`type'", `results'.cond.estimates_ci_per_type[`type'].estimates_ci)
			matrix rownames estimates_ci_type_`type' = Duncan_unit Duncan_ind Theil_unit Theil_ind Atkinson_unit Atkinson_ind Coworker_unit Coworker_ind
			if (`I_noinference') { /* no inference case (begin) */
				matrix colnames estimates_ci_type_`type' = index I_indiv_weight bounds_hat:low bounds_hat:up
			} /* no inference case (end) */
			else { /* inference case (begin) */
				if (`I_specified_alpha') { /* specified alpha case (begin) */
					matrix colnames estimates_ci_type_`type' = index I_indiv_weight bounds_hat:low bounds_hat:up ///
						I_ciboundary ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up ci`alpha_text':low ci`alpha_text':up	
				} /* specified alpha case (end) */
				else { /* no specified alpha case (begin) */
					matrix colnames estimates_ci_type_`type' = index I_indiv_weight bounds_hat:low bounds_hat:up ///
						I_ciboundary ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up
				} /* no specified alpha case (end) */
			} /* inference case (end) */
			ereturn matrix estimates_ci_type_`type' = estimates_ci_type_`type'
		} /* loop over types (end) */	
		/* matrices of estimates of the bounds and confidence interval per type (end) */
		
	} /* conditional case (end) */
	
	else { /* unconditional case (begin) */
		
		Output_info_data_uncond /* general information about data */
	
		/* I_constrained_case (begin) */
		mata: st_numscalar("I_constrained_case", `results'.uncond.I_all_Fpk_constrained) 
		ereturn scalar I_constrained_case = I_constrained_case
		/* I_constrained_case (end) */
	
		/* matrix of estimates of the bounds and confidence intervals (begin) */
		mata: st_matrix("estimates_ci", `results'.uncond.estimates_ci)
		matrix rownames estimates_ci = Duncan_unit Duncan_ind Theil_unit Theil_ind Atkinson_unit Atkinson_ind Coworker_unit Coworker_ind
		if (`I_noinference') { /* no inference case (begin) */
			matrix colnames estimates_ci = index I_indiv_weight bounds_hat:low bounds_hat:up
		} /* no inference case (end) */
		else { /* inference case (begin) */
			if (`I_specified_alpha') { /* specified alpha case (begin) */
				matrix colnames estimates_ci = index I_indiv_weight bounds_hat:low bounds_hat:up ///
					I_ciboundary ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up ci`alpha_text':low ci`alpha_text':up	
				local col_ci_l 12
				local col_ci_u 13
			} /* specified alpha case (end) */
			else { /* no specified alpha case (begin) */
				matrix colnames estimates_ci = index I_indiv_weight bounds_hat:low bounds_hat:up ///
					I_ciboundary ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up
				local col_ci_l 8
				local col_ci_u 9	
			} /* no specified alpha case (end) */
		} /* inference case (end) */
		ereturn matrix estimates_ci = estimates_ci
		/* matrix of estimates of the bounds and confidence intervals (end) */
		
		/* elements of the matrix that will be displayed in Stata log (begin) */
		forvalues index = 1/8 {
			mata: st_numscalar("bounds_l_`index'", `results'.uncond.estimates_ci[`index',3])
			mata: st_numscalar("bounds_u_`index'", `results'.uncond.estimates_ci[`index',4])
			if (`I_noinference') {
				mata: st_numscalar("ci_l_`index'", .)
				mata: st_numscalar("ci_u_`index'", .)
			}
			else {
				mata: st_numscalar("ci_l_`index'", `results'.uncond.estimates_ci[`index',`col_ci_l'])
				mata: st_numscalar("ci_u_`index'", `results'.uncond.estimates_ci[`index',`col_ci_u'])
			}
		}	
		/* elements of the matrix that will be displayed in Stata log (end) */
	
		/* matrix infodistributionofp_np (begin) */
		mata: st_matrix("info_distribution_of_p", `results'.uncond.store_infodistributionofp_np)
		if (`I_independencekp') { /* case K and p assumed independent (begin) */
			matrix rownames info_distribution_of_p = info_data p_uncond:nodes p_uncond:masses
			matrix colnames info_distribution_of_p = K_max nb_units prop_unit I_cons_Fp nb_nodes xy
		} /* case K and p assumed independent (end) */
		else { /* case K and p not assumed independent (begin) */
			local rownames
			forvalues k = 1/`=Kbar'{
				local rownames `rownames' case_K=`k':info p_cond_K=`k':nodes p_cond_K=`k':masses
			}
			matrix rownames info_distribution_of_p = `rownames'	
			matrix colnames info_distribution_of_p = K nb_units prop_unit I_cons_Fpk nb_nodes xy		
		} /* case K and p not assumed independent (end) */
		ereturn matrix info_distribution_of_p = info_distribution_of_p
		/* matrix infodistributionofp_np (end) */
		
		/* matrix test_binomial_results when computed (begin) */
		if ((`I_independencekp'&`I_testbinomial'&(`repbootstrap'>0)) | ((!`I_independencekp')&(((!`I_noinference')&(`repbootstrap'>0)) | ((`I_noinference')&(`I_testbinomial')&(`repbootstrap'>0))))) {
			mata: st_matrix("test_binomial_results", `results'.uncond.test_binomial_results)
			matrix colnames test_binomial_results = stat_test p_value
			matrix rownames test_binomial_results = test_result
			ereturn matrix test_binomial_results = test_binomial_results
		}
		/* matrix test_binomial_results when computed (end) */

	} /* unconditional case (end) */
	
	
	/**********************/
	/* log output (begin) */
	/**********************/
	
	display as text "" /* to leave some space, but less that display _newline */
	display as text "Bounds for segregation indices using non parametric (np) method:"
	display as text "{hline 65}"
	
	Display_information_data `I_conditional' `I_noinference' `I_independencekp' `repbootstrap' `I_withsingle' `I_excludingsinglepertype' 1 `repct'
		
	local col_indent_first_col 4
	local col_weight 16
	local col_vert_line 30
	local col_bounds_l 33
	local col_bounds_u 48
	local col_conf_int 63
	local width_format 7
	display as text _col(`col_indent_first_col') "Index" _col(`col_weight') "Weight-level" _col(`col_vert_line') "{c |}" ///
		_col(`col_bounds_l') "Lower bound" _col(`col_bounds_u') "Upper bound" _col(`col_conf_int') "[`alpha_text'% Conf. Interval]"
	display as text "{hline `=`col_vert_line'-1'}{c +}{hline 53}"
	local name_1 Duncan
	local name_2 Duncan
	local name_3 Theil
	local name_4 Theil
	local name_5 Atkinson
	local name_6 Atkinson
	local name_7 Coworker
	local name_8 Coworker
	
	if (`I_conditional'){
		if (I_unit_level_characteristic) local I_show_unit_level 1
		else local I_show_unit_level 0
	}
	else local I_show_unit_level 1
	
	forvalues index = 1/8 {
		if ((`index' == 1) | (`index' == 3) | (`index' == 5) | (`index' == 7)) {			
			if (`I_show_unit_level') {
				Display_line_np `name_`index'' "unit" bounds_l_`index' bounds_u_`index' ci_l_`index' ci_u_`index' ///
					`col_indent_first_col' `col_weight' `col_vert_line' `col_bounds_l' `col_bounds_u' `col_conf_int' `width_format' `I_noinference'
			}
		}
		else {
			Display_line_np `name_`index'' "individual" bounds_l_`index' bounds_u_`index' ci_l_`index' ci_u_`index' ///
				`col_indent_first_col' `col_weight' `col_vert_line' `col_bounds_l' `col_bounds_u' `col_conf_int' `width_format' `I_noinference'
		}
	}
	
	if (`I_testbinomial') {
		Display_test_binomial `results' `repbootstrap'
	}
	
	/********************/
	/* log output (end) */
	/********************/
	
end
/******************************************************************************/

/******************************************************************************/
capture program drop Display_line_np
program define Display_line_np
	args name_index name_weight bounds_l bounds_u ci_low ci_up ///
		col_indent_first_col col_weight col_vert_line col_bounds_l col_bounds_u col_conf_int width_format I_noinference
	if(`I_noinference') local indent_conf_int 5
	else local indent_conf_int 2
	local space_center_result 3
	display as text _col(`col_indent_first_col') "`name_index'"  _col(`col_weight') "`name_weight'" _col(`col_vert_line') "{c |}" ///
		as text _col(`=`col_bounds_l'+`space_center_result'') as result %-`width_format'.0g `bounds_l' ///
		as text _col(`=`col_bounds_u'+`space_center_result'') as result %-`width_format'.0g `bounds_u' ///
		as text _col(`=`col_conf_int'+`indent_conf_int'') as result %-`width_format'.0g `ci_low' _skip(2) as result %-`width_format'.0g `ci_up'
end
/******************************************************************************/

/******************************************************************************/
capture program drop Output_beta_method
program define Output_beta_method, eclass

	args I_conditional results I_specified_alpha alpha_text I_noinference I_independencekp I_testbinomial repbootstrap I_withsingle I_excludingsinglepertype repct
 	
	ereturn clear
	Output_info_estimation_inference `results' /* save and ereturn user's choices for estimation and inference */

	local rownames_estimates_ci "Duncan_unit Duncan_ind Theil_unit Theil_ind Atkinson_unit Atkinson_ind Coworker_unit Coworker_ind Gini_unit Gini_ind"
	local nb_seg_weight 10
	local col_index_hat_in_mata 3
	local col_ci_l_with_alpha 10
	local col_ci_u_with_alpha 11
	local col_ci_l_without_alpha 6
	local col_ci_u_without_alpha 7
	
	if (`I_conditional') { /* conditional case (begin) */

		Output_info_data_cond /* general information about data */

		/* matrix of estimates of the indices and confidence intervals aggregated (begin) */
		mata: st_matrix("estimates_ci_aggregated", `results'.cond.estimates_ci_aggregated)		
		matrix rownames estimates_ci_aggregated = `rownames_estimates_ci'
		if (`I_noinference') { /* no inference case (begin) */
			matrix colnames estimates_ci_aggregated = index I_indiv_weight index_hat
		} /* no inference case (end) */		
		else { /* inference case (begin) */	
			if (`I_specified_alpha') { /* specified alpha case (begin) */
				matrix colnames estimates_ci_aggregated = index I_indiv_weight index_hat ///
					ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up ci`alpha_text':low ci`alpha_text':up	
				local col_ci_l = `col_ci_l_with_alpha'
				local col_ci_u = `col_ci_u_with_alpha'
			} /* specified alpha case (end) */
			else { /* no specified alpha case (begin) */		
				matrix colnames estimates_ci_aggregated = index I_indiv_weight index_hat ///
					ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up			
				local col_ci_l = `col_ci_l_without_alpha'
				local col_ci_u = `col_ci_u_without_alpha'			
			} /* no specified alpha case (end) */
		} /* inference case (end) */
		ereturn matrix estimates_ci_aggregated = estimates_ci_aggregated	
		/* matrix of estimates of the bounds and confidence interval aggregated (end) */
	
		/* elements of the matrix that will be displayed in Stata log (begin) */
		forvalues index = 1/`nb_seg_weight' {
			mata: st_numscalar("index_hat_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_index_hat_in_mata'])
			if (`I_noinference') {
				mata: st_numscalar("ci_l_`index'", .)
				mata: st_numscalar("ci_u_`index'", .)
			}
			else {
				mata: st_numscalar("ci_l_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_ci_l'])
				mata: st_numscalar("ci_u_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_ci_u'])
			}
		}	
		/* elements of the matrix that will be displayed in Stata log (end) */

		/* matrix info_distribution_of_p: not available with conditional analyses */
		
		/* matrix test_binomial_results: not available with conditional analyses */
		
		/* matrices of estimates of the bounds and confidence interval per type (begin) */
		forvalues type = 1/`=nb_types' { /* loop over types (begin) */
			mata: st_matrix("estimates_ci_type_`type'", `results'.cond.estimates_ci_per_type[`type'].estimates_ci)
			matrix rownames estimates_ci_type_`type' = `rownames_estimates_ci'
			if (`I_noinference') { /* no inference case (begin) */
				matrix colnames estimates_ci_type_`type' = index I_indiv_weight index_hat
			} /* no inference case (end) */
			else { /* inference case (begin) */
				if (`I_specified_alpha') { /* specified alpha case (begin) */
					matrix colnames estimates_ci_type_`type' = index I_indiv_weight index_hat ///
						ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up ci`alpha_text':low ci`alpha_text':up	
				} /* specified alpha case (end) */
				else { /* no specified alpha case (begin) */
					matrix colnames estimates_ci_type_`type' = index I_indiv_weight index_hat ///
						ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up
				} /* no specified alpha case (end) */
			} /* inference case (end) */
			ereturn matrix estimates_ci_type_`type' = estimates_ci_type_`type'
		} /* loop over types (end) */	
		/* matrices of estimates of the bounds and confidence interval per type (end) */
		
	} /* conditional case (end) */

	else { /* unconditional case (begin) */
		
		Output_info_data_uncond /* general information about data */
	
		/* matrix of estimates (point) and confidence intervals (begin) */
		mata: st_matrix("estimates_ci", `results'.uncond.estimates_ci)
		matrix rownames estimates_ci = `rownames_estimates_ci'
		if (`I_noinference') { /* no inference case (begin) */
			matrix colnames estimates_ci = index I_indiv_weight index_hat
		} /* no inference case (end) */
		else { /* inference case (begin) */
			if (`I_specified_alpha') { /* specified alpha case (begin) */
				matrix colnames estimates_ci = index I_indiv_weight index_hat ///
					ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up ci`alpha_text':low ci`alpha_text':up			
				local col_ci_l = `col_ci_l_with_alpha'		
				local col_ci_u = `col_ci_u_with_alpha'
			} /* specified alpha case (end) */
			else { /* no specified alpha case (begin) */
				matrix colnames estimates_ci = index I_indiv_weight index_hat ///
					ci99:low ci99:up ci95:low ci95:up ci90:low ci90:up
				local col_ci_l = `col_ci_l_without_alpha'
				local col_ci_u = `col_ci_u_without_alpha'
			} /* no specified alpha case (end) */
		} /* inference case (end) */
		ereturn matrix estimates_ci = estimates_ci
		/* matrix of estimates (point) and confidence intervals (end) */
		
		/* elements of the matrix that will be displayed in Stata log (begin) */
		forvalues index = 1/`nb_seg_weight' {
			mata: st_numscalar("index_hat_`index'", `results'.uncond.estimates_ci[`index',`col_index_hat_in_mata'])
			if (`I_noinference') {
				mata: st_numscalar("ci_l_`index'", .)
				mata: st_numscalar("ci_u_`index'", .)
			}
			else {
				mata: st_numscalar("ci_l_`index'", `results'.uncond.estimates_ci[`index',`col_ci_l'])
				mata: st_numscalar("ci_u_`index'", `results'.uncond.estimates_ci[`index',`col_ci_u'])
			}
		}	
		/* elements of the matrix that will be displayed in Stata log (end) */

		/* matrix store_infodistributionofp_beta (begin) */
		mata: st_matrix("info_distribution_of_p", `results'.uncond.store_infodistributionofp_beta)
		if (`I_independencekp') { /* case K and p assumed independent (begin) */
			matrix rownames info_distribution_of_p = info_p
			matrix colnames info_distribution_of_p = K_max nb_units prop_unit nb_comp alpha beta /* only for nb_comp = c = 1 */
		} /* case K and p assumed independent (end) */
		else { /* case K and p not assumed independent (begin) */
			matrix rownames info_distribution_of_p = info_p_|_K
			matrix colnames info_distribution_of_p = K nb_units prop_unit nb_comp alpha beta /* only for nb_comp = c = 1 */
		} /* case K and p not assumed independent (end) */
		ereturn matrix info_distribution_of_p = info_distribution_of_p
		/* matrix infodistributionofp_np (end) */
				
		/* matrix test_binomial_results when computed (begin) */
		if (`I_testbinomial'&(`repbootstrap'>0)) {
			mata: st_matrix("test_binomial_results", `results'.uncond.test_binomial_results)
			matrix colnames test_binomial_results = stat_test p_value
			matrix rownames test_binomial_results = test_result
			ereturn matrix test_binomial_results = test_binomial_results
		}
		/* matrix test_binomial_results when computed (end) */

	} /* unconditional case (end) */
	
	
	/**********************/
	/* log output (begin) */
	/**********************/

	display as text "" /* to leave some space, but less that display _newline */
	display as text "Estimates for segregation indices using parametric (beta) method:"
	display as text "{hline 66}"
	
	Display_information_data `I_conditional' `I_noinference' `I_independencekp' `repbootstrap' `I_withsingle' `I_excludingsinglepertype' 1 `repct'
	
	local col_indent_first_col 4
	local col_weight 16
	local col_vert_line 30
	local col_index_hat 33	
	local col_conf_int 51
	local width_format 7
	display as text _col(`col_indent_first_col') "Index" _col(`col_weight') "Weight-level" _col(`col_vert_line') "{c |}" ///
		_col(`col_index_hat') "Point estimate" ///
		_col(`col_conf_int') "[`alpha_text'% Conf. Interval]"
	display as text "{hline `=`col_vert_line'-1'}{c +}{hline 41}" //53
	local name_1 Duncan
	local name_2 Duncan
	local name_3 Theil
	local name_4 Theil
	local name_5 Atkinson
	local name_6 Atkinson
	local name_7 Coworker
	local name_8 Coworker
	local name_9 Gini
	local name_10 Gini
	
	if (`I_conditional'){
		if (I_unit_level_characteristic) local I_show_unit_level 1
		else local I_show_unit_level 0
	}
	else local I_show_unit_level 1
	
	forvalues index = 1/`nb_seg_weight' {
		if ((`index' == 1) | (`index' == 3) | (`index' == 5) | (`index' == 7) | (`index' == 9)) {			
			if (`I_show_unit_level') {
				Display_line_beta `name_`index'' "unit" index_hat_`index' ci_l_`index' ci_u_`index' ///
					`col_indent_first_col' `col_weight' `col_vert_line' `col_index_hat' `col_conf_int' `width_format' `I_noinference'
			}
		}
		else {
			Display_line_beta `name_`index'' "individual" index_hat_`index' ci_l_`index' ci_u_`index' ///
				`col_indent_first_col' `col_weight' `col_vert_line' `col_index_hat' `col_conf_int' `width_format' `I_noinference'
		}
	}
	
	if (`I_testbinomial') {
		Display_test_binomial `results' `repbootstrap'
	}


	/********************/
	/* log output (end) */
	/********************/
	

end
/******************************************************************************/

/******************************************************************************/
capture program drop Display_line_beta
program define Display_line_beta
	args name_index name_weight index_hat ci_low ci_up ///
		col_indent_first_col col_weight col_vert_line col_index_hat col_conf_int width_format I_noinference
	if(`I_noinference') local indent_conf_int 5
	else local indent_conf_int 2
	local space_center_result 4
	display as text _col(`col_indent_first_col') "`name_index'"  _col(`col_weight') "`name_weight'" _col(`col_vert_line') "{c |}" ///
		as text _col(`=`col_index_hat'+`space_center_result'') as result %-`width_format'.0g `index_hat' ///
		as text _col(`=`col_conf_int'+`indent_conf_int'') as result %-`width_format'.0g `ci_low' _skip(2) as result %-`width_format'.0g `ci_up'
end
/******************************************************************************/

/******************************************************************************/
capture program drop Output_ct_method
program define Output_ct_method, eclass

	args I_conditional results I_specified_alpha alpha_text I_noinference I_independencekp I_testbinomial repbootstrap I_withsingle I_excludingsinglepertype repct
	
	ereturn clear
	Output_info_estimation_inference `results' /* save and ereturn user's choices for estimation and inference */

	local rownames_estimates_ci "Duncan Theil Atkinson Coworker Gini"
	if (`I_specified_alpha') {
		local text_perc = 100 - `alpha_text'
		local colnames_estimates_ci "index prop_based mean_und_ra CT_corr sd_under_ra CFC_score perc_1_ra perc_5_ra perc_10_ra perc_90_ra perc_95_ra perc_99_ra perc_`text_perc'_ra perc_`alpha_text'_ra"
	}
	else {
		local colnames_estimates_ci "index prop_based mean_und_ra CT_corr sd_under_ra CFC_score perc_1_ra perc_5_ra perc_10_ra perc_90_ra perc_95_ra perc_99_ra"
	}
	
	local nb_seg_weight 5
	local col_index_hat_in_mata 3
	local col_ci_l_with_alpha 10
	local col_ci_u_with_alpha 11
	local col_ci_l_without_alpha 6
	local col_ci_u_without_alpha 7
	
	if (`I_conditional') { /* conditional case (begin) */

		Output_info_data_cond /* general information about data */

		/* matrix of estimates, data, under random allocation and corrections aggregated (begin) */
		mata: st_matrix("estimates_ci_aggregated", `results'.cond.estimates_ci_aggregated)		
		matrix rownames estimates_ci_aggregated = `rownames_estimates_ci'
		matrix colnames estimates_ci_aggregated = `colnames_estimates_ci'
		ereturn matrix estimates_ci_aggregated = estimates_ci_aggregated
		/* matrix of estimates, data, under random allocation and corrections aggregated (begin) */
			
		/* elements of the matrix that will be displayed in Stata log (begin) */
		local col_index_prop_based 2
		local col_index_expected_under_ra 3
		local col_index_corr_CT 4
		forvalues index = 1/`nb_seg_weight' {
			mata: st_numscalar("index_prop_based_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_index_prop_based'])
			mata: st_numscalar("index_expected_under_ra_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_index_expected_under_ra'])
			mata: st_numscalar("index_corr_CT_`index'", `results'.cond.estimates_ci_aggregated[`index',`col_index_corr_CT'])
		}	
		/* elements of the matrix that will be displayed in Stata log (end) */
	
		/* matrix of estimates, data, under random allocation and corrections per type (begin) */
		forvalues type = 1/`=nb_types' { /* loop over types (begin) */
			mata: st_matrix("estimates_ci_type_`type'", `results'.cond.estimates_ci_per_type[`type'].estimates_ci)
			matrix rownames estimates_ci_type_`type' = `rownames_estimates_ci'
			matrix colnames estimates_ci_type_`type' = `colnames_estimates_ci'
			ereturn matrix estimates_ci_type_`type' = estimates_ci_type_`type'
		} /* loop over types (end) */			
		/* matrix of estimates, data, under random allocation and corrections per type (end) */
		
	} /* conditional case (end) */

	else { /* unconditional case (begin) */
		
		Output_info_data_uncond /* general information about data */
	
		/* matrix of estimates, data, under random allocation and corrections (begin) */
		mata: st_matrix("estimates_ci", `results'.uncond.estimates_ci)
		matrix rownames estimates_ci = `rownames_estimates_ci'
		matrix colnames estimates_ci = `colnames_estimates_ci'
		ereturn matrix estimates_ci = estimates_ci
		/* matrix of estimates, data, under random allocation and corrections (end) */
		
		/* elements of the matrix that will be displayed in Stata log (begin) */
		local col_index_prop_based 2
		local col_index_expected_under_ra 3
		local col_index_corr_CT 4
		forvalues index = 1/`nb_seg_weight' {
			mata: st_numscalar("index_prop_based_`index'", `results'.uncond.estimates_ci[`index',`col_index_prop_based'])
			mata: st_numscalar("index_expected_under_ra_`index'", `results'.uncond.estimates_ci[`index',`col_index_expected_under_ra'])
			mata: st_numscalar("index_corr_CT_`index'", `results'.uncond.estimates_ci[`index',`col_index_corr_CT'])
		}	
		/* elements of the matrix that will be displayed in Stata log (end) */
				
		/* matrix test_binomial_results when computed (begin) */
		if (`I_testbinomial'&(`repbootstrap'>0)) {
			mata: st_matrix("test_binomial_results", `results'.uncond.test_binomial_results)
			matrix colnames test_binomial_results = stat_test p_value
			matrix rownames test_binomial_results = test_result
			ereturn matrix test_binomial_results = test_binomial_results
		}
		/* matrix test_binomial_results when computed (end) */

	} /* unconditional case (end) */
	
	
	/**********************/
	/* log output (begin) */
	/**********************/

	display as text "" /* to leave some space, but less that display _newline */
	display as text "Estimates for segregation indices using proportion-based and CT correction (ct) method:"
	display as text "{hline 88}"
	
	Display_information_data `I_conditional' `I_noinference' `I_independencekp' `repbootstrap' `I_withsingle' `I_excludingsinglepertype' 0 `repct'
	
	local col_indent_first_col 4
	local col_weight 16
	local col_vert_line 30
	local col_proportion_based 33	
	local col_mean_under_ra 53
	local col_CT_corrected 83
	local width_format 7
	display as text _col(`col_indent_first_col') "Index" _col(`col_weight') "Weight-level" _col(`col_vert_line') "{c |}" ///
		_col(`col_proportion_based') "Proportion-based" ///
		_col(`col_mean_under_ra') "Expected under rand. alloc." ///
		_col(`col_CT_corrected') "CT-corrected"
	display as text "{hline `=`col_vert_line'-1'}{c +}{hline 65}"
	local name_1 Duncan
	local name_2 Theil
	local name_3 Atkinson
	local name_4 Coworker
	local name_5 Gini
	
	forvalues index = 1/`nb_seg_weight' {	
		Display_line_ct `name_`index'' "n.a." index_prop_based_`index' index_expected_under_ra_`index' index_corr_CT_`index' ///
			`col_indent_first_col' `col_weight' `col_vert_line' `col_proportion_based' `col_mean_under_ra' `col_CT_corrected' ///
			`width_format'
	}
	
	if (`I_testbinomial') {
		Display_test_binomial `results' `repbootstrap'
	}

	/********************/
	/* log output (end) */
	/********************/

end
/******************************************************************************/

/******************************************************************************/
capture program drop Display_line_ct
program define Display_line_ct
	args name_index name_weight index_prop_based index_expected_under_ra index_corr_CT ///
		col_indent_first_col col_weight col_vert_line col_proportion_based col_mean_under_ra col_CT_corrected ///
		width_format
	local space_center_weight 5
	local space_center_result_1 4
	local space_center_result_2 9
	local space_center_result_3 3
	display as text _col(`col_indent_first_col') "`name_index'"  _col(`=`col_weight'+`space_center_weight'') "`name_weight'" _col(`col_vert_line') "{c |}" ///
		as text _col(`=`col_proportion_based'+`space_center_result_1'') as result %-`width_format'.0g `index_prop_based' ///
		as text _col(`=`col_mean_under_ra'+`space_center_result_2'') as result %-`width_format'.0g `index_expected_under_ra' ///
		as text _col(`=`col_CT_corrected'+`space_center_result_3'') as result %-`width_format'.0g `index_corr_CT'
end
/******************************************************************************/
