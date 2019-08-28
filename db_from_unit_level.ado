capture program drop db_from_unit_level
program db_from_unit_level
	
	version 14.2
	
	/* NB : the function is a void function in a way but creates 
	a mata structure called: struct_db_uncond */
	
	/* Parsing arguments */
	syntax varlist(min=2 max=2 numeric), [withsingleton]
	tokenize `varlist', parse(" ")
	local K `1' 
	local X `2'

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
	
	if ("`withsingleton'" == "") local I_withsingleton = 0
	else local I_withsingleton = 1
	
	if ((r(min) == 1) & (r(max) == 1) & (!`I_withsingleton')){
		display as error "All the units are of size 1 and the option withsingleton is not specified: there is no data"
		exit
	}
	local K_min = r(min)
	local K_max = r(max)
	local K_count = `K_max' - `K_min' + 1
	
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
		
	quietly: count if `K' == 1
	local nb_units_singleton = r(N)
	
	preserve

		if (!`I_withsingleton'){
			quietly: drop if `K' == 1
		}

		/* See for limit of matsize - help Stata 14.1.2 
		A CHANGER PEUT-ETRE POUR MODIFIER LA MATRICE UTILISEE AVEC SEULEMENT
		LES MODALITES REALISEES DE K
		*/
	
		local db_matsize = `K_max'+3
		if (`db_matsize' >= 400)&(`db_matsize' <= 11000) {
			display as text "Some units have more than 400 individuals, need to adjust Stata matsize default"
			display as error "Warning: with such unit size, the estimation could be quite long"
			quietly: set matsize `db_matsize'
			display as text "Matsize adjusted"
		}
		else if `db_matsize' > 110000 {
			display as error "Some units have more than 11,000 individuals"
			display as error "With our data structure, it exceeds the maximal matsize of Stata and estimation cannot be done"
			display as error "With such unit size, you may prefer other methods;"
			display as error "or perhaps abstract from exceptionnally large units and withdraw them from the analysis"
			exit
		}
			
		matrix define db = J(`K_max', `K_max'+3, 0)
		
		/* first colum of db : K size of unit */
		forvalues i = 1 / `K_max' {
			local j =`i'
			matrix define db[`j',1] = `j' 
		}
		
		/* 	fill the matrix : 
			entry [i,j]
			i = 1,...,`K_max'
			j = 3,...,`K_max'
			: number of units with size K = i and X = j-3 
			entry [i,2], i = 1,...,`K_max', number of units of size i in the dataset
			-99 for impossible observation instead of 0 i.e. when X > K
		*/

		tempvar count_K count_XK
		bysort `K': egen `count_K' = count(`X')
		bysort `K' `X': egen `count_XK' = count(`X')	
		quietly: duplicates drop `K' `X', force
		forvalues i = 1 / `=_N' {
			matrix define db[`K'[`i'],`X'[`i']+3] = `count_XK'[`i']
		}
		quietly: duplicates drop `K', force
		forvalues i = 1 / `=_N' {
			matrix define db[`K'[`i'],2] = `count_K'[`i']
		}		
		
		forvalues i = 1 / `=`K_max'-1' {
			*local j_start = `i'+4
			*local j_end = `K_max'+3
			forvalues j = `=`i'+4' / `=`K_max'+3' {
				matrix define db[`i',`j'] = -99
			}
		}
	
	/* define directly a mata structure object struct_db_uncond */
	mata: struct_db_uncond = Struct_db_uncond()	
	
	/* Definition of the matrix db */
	mata: struct_db_uncond.db = st_matrix("db")
	
	/* Defintion of the struct info_data */
	mata: struct_db_uncond.info_data.I_withsingleton = `I_withsingleton'
	mata: struct_db_uncond.info_data.nb_units_studied = Get_nb_units_studied(struct_db_uncond.db)
	if (`I_withsingleton'){
		mata: struct_db_uncond.info_data.nb_units_singleton = struct_db_uncond.db[1,2]
		mata: struct_db_uncond.info_data.nb_units_total = struct_db_uncond.info_data.nb_units_studied
	}
	else {
		mata: struct_db_uncond.info_data.nb_units_singleton = `nb_units_singleton'
		mata: struct_db_uncond.info_data.nb_units_total = struct_db_uncond.info_data.nb_units_studied + `nb_units_singleton'
	}
	mata: struct_db_uncond.info_data.nb_individuals = Compute_S(struct_db_uncond.db)
	mata: struct_db_uncond.info_data.nb_minority_individuals = Compute_S_X(struct_db_uncond.db)
	mata: struct_db_uncond.info_data.prop_minority_hat = Compute_Pi(struct_db_uncond.info_data)
	mata: struct_db_uncond.info_data.list_K_positive_obs = selectindex(struct_db_uncond.db[,2])
	mata: struct_db_uncond.info_data.Kbar = max(struct_db_uncond.info_data.list_K_positive_obs)
	mata: struct_db_uncond.info_data.nb_K_positive_obs = length(struct_db_uncond.info_data.list_K_positive_obs)
	
end
