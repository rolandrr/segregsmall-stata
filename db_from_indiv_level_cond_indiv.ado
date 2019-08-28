capture program drop db_from_indiv_level_cond_indiv
program db_from_indiv_level_cond_indiv
	
	version 14.2

	/* NB : the function is a void function in a way but creates 
	a mata structure used for conditional segregation 
	named struct_db_cond */
	
	/* Parsing arguments */
	syntax varlist(min=3 max=3 numeric), [withsingleton] [excludingsingletonpertype]
	tokenize `varlist', parse(" ")
	local id_unit `1'
	local I_minority `2'
	local W `3'
	
	/* Check for missing data => entail bug when construct the stata matrix db
	and preferable that the user checks her or his data */
	quietly: count if missing(`id_unit')
	if (r(N) > 0) {
		display as error "There are missing values in variable `id_unit' (identifier of the unit an individual belongs to): estimation cannot be done"
		exit
	}
	quietly: count if missing(`I_minority')
	if (r(N) > 0) {
		display as error "There are missing values in variable `I_minority' (indicator of minority individual): estimation cannot be done"
		exit
	}
	quietly: count if missing(`W')
	if (r(N) > 0) {
		display as error "There are missing values in variable `W' (type of the individuals): estimation cannot be done"
		exit
	}
	
	/* Check conformity of variables I_minority */
	quietly: tabulate `I_minority'
	local nb_distinct_values_I_minority = r(r)
	quietly: summarize `I_minority'
	if ((r(min) != 0) | (r(max) != 1) | (`nb_distinct_values_I_minority' != 2)){
		display as error "There are invalid values in variable `I_minority' (indicator of minority individual): it is not a binary 0-1 variable"	
		exit
	}	
	
	/* Check for values of W - categorical variable in {1,2,3,...,\bar{W}}
	and definition of nb_types */
	preserve 
		quietly: duplicates drop `W', force 
		sort `W'
		tempvar steps_W 
		quietly: generate `steps_W' = `W'[_n+1] - `W'[_n]
		quietly: replace `steps_W' = 1 if _n == _N
		summarize `steps_W', meanonly
		local mean_steps_W = r(mean)
		quietly: tabulate `W'
		local nb_distinct_values_W = r(r)
		quietly: summarize `W'
		if ((r(min) != 1) | (r(max)<2) | (`mean_steps_W' != 1) | (`nb_distinct_values_W' != r(max))){
			display as error "Variable `W' (type of the individuals) is invalid: it should be an individual-level categorical variable with modalities (at least two) numbered 1,2,3..."
			exit
		}
		local W_bar = r(max)
	restore
	
	/* Instanciation of struct_db_cond and definition of its attributes */
	mata: struct_db_cond = Struct_db_cond()
	mata: struct_db_cond.nb_types = `W_bar'
	mata: struct_db_cond.I_unit_level_characteristic = 0
	mata: struct_db_cond.I_excludingsingletonpertype = (("`excludingsingletonpertype'" != "") ? 1 : 0)
	if ("`withsingleton'" == "") db_from_indiv_level `id_unit' `I_minority'
	else db_from_indiv_level `id_unit' `I_minority', withsingleton
	mata: struct_db_cond.info_data_uncond = struct_db_uncond.info_data

	/*
	/* Older version: Definition of attribute store_WK_nb_units (begin) */
	preserve
		quietly{		
			keep `id_unit' `W'
			bysort `id_unit': egen count_unit = count(`id_unit')
			bysort `id_unit' `W': egen count_unit_type = count(`id_unit')
			bysort `id_unit' `W': keep if _n == 1
			if ("`withsingleton'" == "") drop if count_unit == 1
			drop count_unit
			if ("`excludingsingletonpertype'" != "") drop if count_unit_type == 1
		
			reshape wide count_unit_type, i(`id_unit') j(`W')
			
			local list_type
			forvalues type = 1/`W_bar' {
				replace count_unit_type`type' = 0 if count_unit_type`type' == .
				local list_type `list_type' count_unit_type`type'
			}
			tempvar nb_unit
			bysort `list_type': egen `nb_unit' = count(`id_unit')
			bysort `list_type': keep if _n == 1
			drop `id_unit'	
		}
		mata: struct_db_cond.store_WK_nb_units = Cons_store_WK_nb_units(st_data(.,.))
	restore
	/* Definition of attribute store_WK_nb_units (end) */
	*/
	
	/* Definition of attribute store_WK_nb_units and store_WXK_nb_units (begin) */		
	preserve
	quietly{
		keep `id_unit' `W' `I_minority' 
		bysort `id_unit': egen per_unit_K = count(`id_unit')
		bysort `id_unit' `W': egen per_unit_type_K = count(`id_unit')
		bysort `id_unit' `W': egen per_unit_type_X = total(`I_minority')
		bysort `id_unit' `W': keep if _n == 1
		if ("`withsingleton'" == "") drop if per_unit_K == 1
		drop per_unit_K `I_minority'
		if ("`excludingsingletonpertype'" != "") drop if per_unit_type_K == 1
			
		reshape wide per_unit_type_K per_unit_type_X, i(`id_unit') j(`W')

		local list_type_KX
		local list_type_K
		forvalues type = 1/`W_bar' {
			replace per_unit_type_K`type' = 0 if per_unit_type_K`type' == .
			replace per_unit_type_X`type' = 0 if per_unit_type_X`type' == .
			local list_type_KX `list_type_KX' per_unit_type_K`type' per_unit_type_X`type'
			local list_type_K `list_type_K' per_unit_type_K`type'
		}
		
		bysort `list_type_KX': egen nb_unit_KX = count(`id_unit')
		bysort `list_type_K': egen nb_unit_K = count(`id_unit')
		
		bysort `list_type_KX': keep if _n == 1
		drop `id_unit'
		
		mata: struct_db_cond.store_WXK_nb_units = st_data(.,(1..(cols(st_data(.,.))-1)))
		
		drop per_unit_type_X*
		drop nb_unit_KX
				
		bysort `list_type_K': keep if _n == 1
		
		mata: struct_db_cond.store_WK_nb_units = Cons_store_WK_nb_units(st_data(.,.))
	}	
	restore	
	/* Definition of attribute store_WK_nb_units and store_WXK_nb_units (end) */
	
	
	/* Definition of database for each subsample of occurrences such that type ind W = w 
	for that, we use the ado-file db_from_unit_level */
	mata: struct_db_cond.db_per_type = Struct_db(struct_db_cond.nb_types)
	mata: struct_db_cond.info_data_per_type = Info_data(struct_db_cond.nb_types)
	
	tempvar K
	bysort `id_unit': egen `K' = count(`I_minority')
	forvalues w = 1/`W_bar' {
		preserve
			quietly{
				keep if `W' == `w'
				if ("`withsingleton'" == "") drop if `K' == 1				
				tempvar K_per_type X_per_type
				bysort `id_unit': egen `K_per_type' = count(`I_minority')
				bysort `id_unit': egen `X_per_type' = total(`I_minority')
				bysort `id_unit': keep if _n == 1
				capture isid `id_unit'
				if (_rc){
					display as error "Variable `id_unit' (identifier of the unit an individual belongs to) is invalid: it does not uniquely define the units"
					exit 459
				}
				if ("`excludingsingletonpertype'" == "") db_from_unit_level `K_per_type' `X_per_type', withsingleton
				else db_from_unit_level `K_per_type' `X_per_type'
			}			
		restore
		mata: struct_db_cond.db_per_type[`w'].db = struct_db_uncond.db
		mata: struct_db_cond.info_data_per_type[`w'] = struct_db_uncond.info_data
	}
	
	/* Definition of different attributes */
	/* for check:
	if ("`withsingleton'" == "") mata: struct_db_cond.nb_units = sum(select(struct_db_cond.store_WK_nb_units[,5], struct_db_cond.store_WK_nb_units[,6] :> 1))
	else mata: struct_db_cond.nb_units = sum(struct_db_cond.store_WK_nb_units[,5])
	*/
	mata: struct_db_cond.nb_units_studied_per_type = Cons_nb_units_per_type_ind(struct_db_cond)
	mata: struct_db_cond.type_frequencies = Cons_type_frequencies_ind(struct_db_cond)
	mata: struct_db_cond.type_probabilities = Cons_type_probabilities_ind(struct_db_cond)
	mata: struct_db_cond.summary_info_data_per_type = Cons_summary_info_data_per_type(struct_db_cond)
	mata: struct_db_cond.info_data_uncond = Corr_info_data_uncond_condind(struct_db_cond)
	mata: struct_db_cond.nb_singleton_cells_all_type = Cons_nb_singleton_cells_all_type(struct_db_cond)
		
end
