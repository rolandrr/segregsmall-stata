capture program drop db_from_indiv_level
program db_from_indiv_level
	
	version 14.2
		
	/* Parsing arguments */
	syntax varlist(min=2 max=2 numeric), [withsingleton]
	tokenize `varlist', parse(" ")
	local id_unit `1'
	local I_minority `2'
	
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
	
	/* Check conformity of variables I_minority */
	quietly: tabulate `I_minority'
	local nb_distinct_values_I_minority = r(r)
	quietly: summarize `I_minority'
	if ((r(min) != 0) | (r(max) != 1) | (`nb_distinct_values_I_minority' != 2)){
		display as error "There are invalid values in variable `I_minority' (indicator of minority individual): it is not a binary 0-1 variable"	
		exit
	}
	
	/* Generate variables K and X, get unit-level database and call db_from_unit_level */
	preserve
		tempvar K X
		bysort `id_unit': egen `K' = count(`I_minority')
		bysort `id_unit': egen `X' = total(`I_minority')
		quietly: bysort `id_unit': keep if _n == 1
		capture isid `id_unit'
		if (_rc){
			display as error "Variable `id_unit' (identifier of the unit an individual belongs to) is invalid: it does not uniquely define the units"
			exit 459
		}
		if ("`withsingleton'" == "") db_from_unit_level `K' `X'
		else db_from_unit_level `K' `X', withsingleton		

end
