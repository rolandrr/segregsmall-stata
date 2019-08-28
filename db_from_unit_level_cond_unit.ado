capture program drop db_from_unit_level_cond_unit
program db_from_unit_level_cond_unit
	
	version 14.2
	
	/* NB : the function is a void function in a way but creates 
	a mata structure used for conditional segregation
	named struct_db_cond */
	
	/* Parsing arguments */
	syntax varlist(min=3 max=3 numeric), [withsingleton]
	tokenize `varlist', parse(" ")
	local K `1' 
	local X `2'
	local Z `3' /* Z is unit-level characteristic, following notations of DHR */
	
	/* Check for missing data => entail bug when construct the stata matrix db
	and preferable that the user checks her or his data */
	quietly: count if missing(`K')
	if (r(N) > 0) {
		display as error "There are missing values in variable `K' (size of the units): estimation cannot be done"
		exit
	}
	quietly: count if missing(`X')
	if (r(N) > 0) {
		display as error "There are missing values in variable `X' (number of minority individuals): estimation cannot be done"
		exit
	}
	quietly: count if missing(`Z')
	if (r(N) > 0) {
		display as error "There are missing values in variable `Z' (type of the units): estimation cannot be done"
		exit
	}
	
	/* Check that the values (different than the format which can be, although some memory loss
	a double for instance) of `K' are strictly positive integers, and of `X' are null or positive integers */
	
	quietly: summarize `X'
	if (r(min) < 0){
		display as error "There are invalid values in variable `X' (number of minority individuals): strictly lower than 0"
		exit
	}
	
	quietly: summarize `K'
	if (r(min) < 1){
		display as error "There are invalid values in variable `K' (size of the units): strictly lower than 1"
		exit
	}
	if ((r(min) == 1) & (r(max) == 1) & ("`withsingleton'" == "")){
		display as error "All the units are of size 1 and the option withsingleton is not specified: there is no data"
		exit
	}
	
	tempvar int_K int_X diff_int_K diff_int_X
	gen `int_K' = int(`K')
	gen `int_X' = int(`X')
	egen `diff_int_K' = diff(`K' `int_K')
	egen `diff_int_X' = diff(`X' `int_X')
	summarize `diff_int_K', meanonly
	if (r(mean)) {
		display as error "There are invalid values in variable `K' (size of the units): non integer values"
		exit
	}
	summarize `diff_int_X', meanonly
	if (r(mean)) {
		display as error "There are invalid values in variable `X' (number of minority individuals): non integer values"
		exit
	}
	
	/* Check that for a given unit (row), there is no more minority individuals 
	in a unit than its size */
	tempvar I_K_leq_X
	gen `I_K_leq_X' = (`K' >= `X')
	summarize `I_K_leq_X', meanonly
	if (r(mean) != 1) {
		display as error "There are invalid values in variables `K' and `X': for some units, `X' > `K'"
		exit
	}
	
	/* Check for values of Z - categorical variable in {1,2,3,...,\bar{Z}
	and definition of nb_types */
	preserve 
		quietly: duplicates drop `Z', force 
		sort `Z'
		tempvar steps_Z 
		quietly: generate `steps_Z' = `Z'[_n+1] - `Z'[_n]
		quietly: replace `steps_Z' = 1 if _n == _N
		summarize `steps_Z', meanonly
		local mean_steps_Z = r(mean)
		quietly: tabulate `Z'
		local nb_distinct_values_Z = r(r)
		quietly: summarize `Z'
		if ((r(min) != 1) | (r(max)<2) | (`mean_steps_Z' != 1) | (`nb_distinct_values_Z' != r(max))) {
			display as error "Variable `Z' (type of the units) is invalid: it should be a unit-level categorical variable with modalities (at least two) numbered 1,2,3..."
			exit
		}
		local Z_bar = r(max)
	restore
	
	/* Instanciation of struct_db_cond and definition of its attributes */
	mata: struct_db_cond = Struct_db_cond()
	mata: struct_db_cond.I_unit_level_characteristic = 1
	mata: struct_db_cond.nb_types = `Z_bar'
	if ("`withsingleton'" == "") db_from_unit_level `K' `X'
	else db_from_unit_level `K' `X', withsingleton
	mata: struct_db_cond.info_data_uncond = struct_db_uncond.info_data
	mata: struct_db_cond.I_excludingsingletonpertype = . /* option available only for individual covariates analysis */
	
	/* Definition of database for each subsample of units such that Z = z 
	for that, we use the ado-file db_from_unit_level */
	mata: struct_db_cond.db_per_type = Struct_db(struct_db_cond.nb_types)
	mata: struct_db_cond.info_data_per_type = Info_data(struct_db_cond.nb_types)
	forvalues type = 1/`Z_bar' {
		preserve
			quietly: keep if `Z' == `type'
			if ("`withsingleton'" == "") db_from_unit_level `K' `X'
			else db_from_unit_level `K' `X', withsingleton
		restore	
		mata: struct_db_cond.db_per_type[`type'].db = struct_db_uncond.db
		mata: struct_db_cond.info_data_per_type[`type'] = struct_db_uncond.info_data
	}
	
	mata: struct_db_cond.store_ZK_nb_units = Cons_store_ZK_nb_units(struct_db_cond)
	mata: struct_db_cond.type_frequencies = Cons_type_frequencies_unit(struct_db_cond)
	mata: struct_db_cond.type_probabilities = Cons_type_probabilities_unit(struct_db_cond)
	mata: struct_db_cond.nb_units_studied_per_type = struct_db_cond.type_frequencies /* idem with unit-level covariates */
	mata: struct_db_cond.summary_info_data_per_type = Cons_summary_info_data_per_type(struct_db_cond)
				
end
