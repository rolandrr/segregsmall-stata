/* structures used to save information about data and the results of
estimation and inference */

/******************************************************************************/
/* Structures to save information from data and estimation / inference */
/******************************************************************************/

capture mata: mata drop Info_data()
mata:
mata set matastrict on
struct Info_data {
	real scalar I_withsingleton
	real scalar nb_units_total, nb_units_singleton, nb_units_studied
	/* nb_units_total and nb_units_singleton keep the same values whatever I_withsingleton
	nb_units_total: total number of units in the data
	nb_units_singleton: total number of singleton units (size = 1) in the data
	nb_units_studied = number of units in the analysis
	nb_units_studied = nb_units_total - (1-I_withsingleton) * nb_units_singleton */
	real scalar nb_individuals, nb_minority_individuals, prop_minority_hat
	/* NB: over the population of unit studied (depending on option I_withsingleton) 
	nb_individuals = number of individuals studied i.e. within the population of studied units, 
	nb_minority_individuals = idem minority individuals
	prop_minority_hat = nb_minority_individuals / nb_individuals */
	real scalar Kbar
	/* largest unit size in the data */
	real colvector list_K_positive_obs
	/* list of sizes of the units present in data */
	real scalar nb_K_positive_obs
	/* length of this list: number of distinct sizes of the units in data */
}
end

capture mata: mata drop Struct_db_uncond()
mata:
mata set matastrict on
struct Struct_db_uncond {
	struct Info_data scalar info_data 
	real matrix db /* triangle matrix of data, unconditional analysis */
}
end

capture mata: mata drop Struct_db()
mata:
mata set matastrict on
struct Struct_db {
	real matrix db /* to stock the triangle matrix of data type by type in conditional analysis */
}
end

capture mata: mata drop Struct_db_cond()
mata:
mata set matastrict on
struct Struct_db_cond {
	struct Info_data scalar info_data_uncond 
	/* neglecting type, information about data as aggregated, as for unconditional case 
	but attention for individual-level conditional analyses, requires correction for 
	number of individuals, in case of option excludingsinglepertype */
	real scalar nb_types
	real scalar I_unit_level_characteristic
	real scalar I_excludingsingletonpertype
	struct Struct_db colvector db_per_type
	struct Info_data colvector info_data_per_type 
	/* type by type ; for cells (unit x individual-type) in the individual level covariate case
	in particular info_data_per_type[type].
	nb_units_total and nb_units_singleton are the same whatever the chosen options
	I_excludingsingletonpertype and info_data_uncond.I_withsingleton
	but not the case for nb_units_studied, it depends on those options */
	real matrix store_ZK_nb_units
	real matrix store_WK_nb_units
	real matrix store_WXK_nb_units
	real colvector type_frequencies, type_probabilities, nb_units_studied_per_type
	/* same names for both unit-level or individual-level analysis
	but difference in the meaning of course :
	type_frequencies : nunber of studied units of each type in unit-level analysis
	type_frequencies : nunber of individuals (within studied units) of each type in individual-level analysis
	type_probabilities : type_frequencies expressed in proportion of sum(type_frequencies)
	nb_units_studied_per_type : idem to type_frequencies in unit-level analysis
	nb_units_studied_per_type : a unit is played by a cell unit size x individual covariate in individual-level analysis */
	real matrix summary_info_data_per_type
	/* summary matrix for ereturn of the Stata command
	per type indexed by the row, some elements of info_data_per_type[type]
	column 1: I_withsingleton
	column 2: nb_units_total
	column 3: nb_units_singleton
	column 4: nb_units_studied
	column 5: nb_individuals
	column 6: nb_minority_individuals
	column 7: prop_minority_hat 
	column 8: Kbar 
	column 9: nb_K_positive_obs */	
	real scalar nb_singleton_cells_all_type
	/* sum over type of the singleton cells (unit x individual-type) 
	only relevant for individual-level characteristics
	for unit-level characteristics, the equivalent is simply info_data_uncond.nb_units_singleton
	(as there is not this notion of cells) */
}
end

capture mata: mata drop Info_estimation_inference()
mata:
mata set matastrict on
struct Info_estimation_inference {
	real scalar I_method_np, I_method_beta, I_method_ct
	real scalar I_conditional
	real scalar I_unit_level_characteristic
	real scalar b_atkinson
	real scalar I_noinference
	real scalar nb_bootstrap_repetition
	real scalar specified_alpha
	real scalar I_deltamethod
	real scalar I_hyp_independenceKp
	real scalar I_testbinomial
	real scalar nb_ct_repetition
}
end

capture mata: mata drop Results_uncond()
mata:
mata set matastrict on
struct Results_uncond {
	real matrix estimates_ci
	real rowvector test_binomial_results /* first element: value of the test statistic, second: p-value */
	real matrix store_infodistributionofp_np /* useful for DHR only */
	real matrix store_infodistributionofp_beta /* useful for R only */
	real scalar I_all_Fpk_constrained /* useful for DHR only */
}
end

capture mata: mata drop Struct_eci()
mata:
mata set matastrict on
struct Struct_eci {
	real matrix estimates_ci /* for the result type by type in conditional analysis */
}
end

capture mata: mata drop Results_cond()
mata:
mata set matastrict on
struct Results_cond {
	struct Struct_eci rowvector estimates_ci_per_type /* conditional index per type */
	real matrix estimates_ci_aggregated /* aggregated conditional index */
	real colvector I_all_Fpk_constrained_per_type /* useful for DHR only */
	real scalar I_all_type_all_Fpk_constrained /* useful for DHR only */
	struct struct_In_Pb_allK colvector In_Pb_allK_pertype /* useful for DHR only */
	struct struct_In_mb_Kpind colvector In_mb_Kpind_pertype /* useful for DHR only */
}
end

capture mata: mata drop Results_eci()
mata:
mata set matastrict on
struct Results_eci {
	struct Info_estimation_inference scalar info_eci
	struct Results_uncond scalar uncond
	struct Results_cond scalar cond
}
end


/******************************************************************************/
/* _View function for debugging and coding checks */
/******************************************************************************/

capture mata: mata drop _View_Info_data()
mata:
mata set matastrict on
void _View_Info_data(struct Info_data scalar info_data){
	printf("struct Info_data \n")
	printf("I_withsingleton =")
	info_data.I_withsingleton
	printf("nb_units_total =")
	info_data.nb_units_total
	printf("nb_units_singleton =")
	info_data.nb_units_singleton
	printf("nb_units_studied =")
	info_data.nb_units_studied
	printf("nb_individuals =")
	info_data.nb_individuals
	printf("nb_minority_individuals =")
	info_data.nb_minority_individuals
	printf("prop_minority_hat =")
	info_data.prop_minority_hat
	printf("Kbar =")
	info_data.Kbar
	printf("list_K_positive_obs =\n")
	info_data.list_K_positive_obs
	printf("nb_K_positive_obs =")
	info_data.nb_K_positive_obs
}
end

capture mata: mata drop _View_Struct_db_uncond()
mata:
mata set matastrict on
void _View_Struct_db_uncond(struct Struct_db_uncond scalar struct_db_uncond){
	printf("struct Struct_db_uncond \n")
	_View_Info_data(struct_db_uncond.info_data)
	printf("db = \n")
	struct_db_uncond.db
}
end

capture mata: mata drop _View_Struct_db_cond()
mata:
mata set matastrict on
void _View_Struct_db_cond(struct Struct_db_cond scalar struct_db_cond){
	printf("struct Struct_db_cond \n")
	printf("neglecting type, aggregated unconditional info on data \n")
	_View_Info_data(struct_db_cond.info_data_uncond)
	printf("nb_types =")
	struct_db_cond.nb_types
	printf("I_unit_level_characteristic =")
	struct_db_cond.I_unit_level_characteristic
	printf("I_excludingsingletonpertype =")
	struct_db_cond.I_excludingsingletonpertype
	real scalar type
	for(type = 1; type <= struct_db_cond.nb_types; type++){
		printf("Info_data by type, for type = ")
		type
		_View_Info_data(struct_db_cond.info_data_per_type[type])
	}
	printf("store_ZK_nb_units = \n")
	struct_db_cond.store_ZK_nb_units
	printf("store_WK_nb_units = \n")
	struct_db_cond.store_WK_nb_units
	printf("store_WXK_nb_units = \n")
	struct_db_cond.store_WXK_nb_units
	printf("type_frequencies = \n")
	struct_db_cond.type_frequencies
	printf("type_probabilities = \n")
	struct_db_cond.type_probabilities
	printf("nb_units_studied_per_type = \n")
	struct_db_cond.nb_units_studied_per_type
	for(type = 1; type <= struct_db_cond.nb_types; type++){
		printf("db by type, for type = ")
		type
		struct_db_cond.db_per_type[type].db
	}
	printf("summary_info_data_per_type = \n")
	struct_db_cond.summary_info_data_per_type
	printf("nb_singleton_cells_all_type =")
	struct_db_cond.nb_singleton_cells_all_type
}
end

capture mata: mata drop _View_Info_estimation_inference()
mata:
mata set matastrict on
void _View_Info_estimation_inference(struct Info_estimation_inference scalar info_eci){
	printf("struct Info_estimation_inference \n")
	printf("I_method_np =")
	info_eci.I_method_np
	printf("I_method_beta =")
	info_eci.I_method_beta
	printf("I_method_ct =")
	info_eci.I_method_ct
	printf("I_conditional =")
	info_eci.I_conditional
	printf("I_unit_level_characteristic =")
	info_eci.I_unit_level_characteristic
	printf("b_atkinson =")
	info_eci.b_atkinson
	printf("I_noinference =")
	info_eci.I_noinference
	printf("nb_bootstrap_repetition =")
	info_eci.nb_bootstrap_repetition
	printf("specified_alpha =")
	info_eci.specified_alpha
	printf("I_deltamethod =")
	info_eci.I_deltamethod
	printf("I_hyp_independenceKp =")
	info_eci.I_hyp_independenceKp
	printf("I_testbinomial =")
	info_eci.I_testbinomial
	printf("nb_ct_repetition =")
	info_eci.nb_ct_repetition
}
end

capture mata: mata drop _View_Results_uncond()
mata:
mata set matastrict on
void _View_Results_uncond(struct Results_uncond scalar results_uncond){
	printf("struct Results_uncond \n")
	printf("estimates_ci = \n")
	results_uncond.estimates_ci
	printf("test_binomial_results = \n")
	results_uncond.test_binomial_results
	printf("store_infodistributionofp_np = \n")
	results_uncond.store_infodistributionofp_np
	printf("store_infodistributionofp_beta = \n")
	results_uncond.store_infodistributionofp_beta
}
end

capture mata: mata drop _View_Results_cond()
mata:
mata set matastrict on
void _View_Results_cond(struct Results_cond scalar results_cond, real scalar nb_types){
	printf("struct Results_cond \n")
	real scalar type
	for(type = 1; type <= nb_types; type++){
		printf("estimates_ci_per_type, for type =")
		type
		results_cond.estimates_ci_per_type[type].estimates_ci
	}
	printf("estimates_ci_aggregated = \n")
	results_cond.estimates_ci_aggregated
	printf("I_all_Fpk_constrained_per_type = \n")
	results_cond.I_all_Fpk_constrained_per_type
	printf("I_all_type_all_Fpk_constrained =")
	results_cond.I_all_type_all_Fpk_constrained
}
end
