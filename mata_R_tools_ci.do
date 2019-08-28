version 14.2
/* This do-file defines the structure and functions used to perform
inference in the R beta method */

/* draw_db_bootstrap_unit_level */
/* from an original matrix of data db, perform bootstrap (usual bootstrap with 
replacement at the level of the observations, namely here at the level of 
the units) 
It is used to draw from the units for unconditional Rbeta method */
capture mata: mata drop Draw_db_bootstrap_unit_level()
mata:
mata set matastrict on
real matrix Draw_db_bootstrap_unit_level(real matrix db) {

	real scalar nb_observations, Kbar, index, k, x
	real colvector draws_number
	real matrix db_bootstrap_R
	
	nb_observations = quadsum(db[,2], 1)
	Kbar = rows(db)
	
	/* draw the number of the observations that will be in the bootstrap sample */
	draws_number = floor(nb_observations:*runiform(nb_observations,1):+1) 
	
	/* construct database from bootstrap sample */
	db_bootstrap_R = J(Kbar, Kbar+3, 0)
	db_bootstrap_R[,1] = (1..Kbar)'

	index = 0
	for(k = 1; k <= Kbar; k++) {
		for(x = 0; x <= k; x++) {
			if (db[k,x+3] > 0) {
				db_bootstrap_R[k,x+3] = sum((draws_number :>= (index+1)) :* (draws_number :<= (index+db[k,x+3])))
				index = index+db[k,x+3]
			}
		}
	}
	db_bootstrap_R[,2] = rowsum(db_bootstrap_R[.,3..(Kbar+3)])
	for(k = 1; k <= (Kbar-1); k++) {
		db_bootstrap_R[k, ((k+4)..(Kbar+3))] = J(1,Kbar-k, -99)
	}
	
	return(db_bootstrap_R)
}
end

/* draw_struct_db_condunit_unit */
/* Identique au precedent mais dans le cas conditionnel avec unit-level 
caracteristiques : tire dans l'ensemble des unites, en melangeant taille
des unites, type des unites et leur composition en X / K 
Utilise pour bootstrap dans Rbeta method dans le cas conditional unit-level */
capture mata: mata drop Draw_struct_db_condunit_unit()
mata:
mata set matastrict on
struct Struct_db_cond scalar Draw_struct_db_condunit_unit(struct Struct_db_cond scalar struct_db_cond){

	struct Struct_db_cond scalar struct_db_cond_boot
	real scalar nb_observations, index, type, Kbar_type_temp, k, x
	real colvector draws_number
	
	struct_db_cond_boot = Struct_db_cond()
	struct_db_cond_boot.nb_types = struct_db_cond.nb_types 
	/* si pas assez d'observations d'un certain type de telle sorte que le nombre de type a l'issue du bootstrap pourrait changer
	pas forcement de souci, juste une ponderation de 0 pour ce type non represente dans le bootstrap 
	[mais signifie que peut etre pas assez de donnees pour avoir cette finesse de distinction entre les types] */
	
	struct_db_cond_boot.db_per_type = Struct_db(struct_db_cond_boot.nb_types)
	struct_db_cond_boot.info_data_per_type = Info_data(struct_db_cond_boot.nb_types)
	struct_db_cond_boot.nb_units_studied_per_type = J(struct_db_cond_boot.nb_types, 1, .)
	
	nb_observations = struct_db_cond.info_data_uncond.nb_units_studied
	draws_number = floor(nb_observations:*runiform(nb_observations,1):+1) 

	index = 0
	for(type = 1; type <= struct_db_cond_boot.nb_types; type++){ /* loop over types (begin) */

		Kbar_type_temp = struct_db_cond.info_data_per_type[type].Kbar /* temp (for temporary) since it may be that no unit of size Kbar of type = type are drawn */
		struct_db_cond_boot.db_per_type[type].db = J(Kbar_type_temp, Kbar_type_temp+3, 0)
		struct_db_cond_boot.db_per_type[type].db[,1] = (1..Kbar_type_temp)'

		/* draw of the units for struct_db_cond.db_per_type[type].db (begin) */
		for(k = 1; k <= Kbar_type_temp; k++){
			for(x = 0; x <= k; x++){
				if (struct_db_cond.db_per_type[type].db[k,x+3] > 0) {
						struct_db_cond_boot.db_per_type[type].db[k,x+3] = sum((draws_number :>= (index+1)) :* (draws_number :<= (index+struct_db_cond.db_per_type[type].db[k,x+3])))
						index = index + struct_db_cond.db_per_type[type].db[k,x+3]
				}
			}
		}
		/* draw of the units for struct_db_cond.db_per_type[type].db (end) */
		
		/* clean of db_per_type[type].db and attributes of info_data_per_type[type] (begin) */
		struct_db_cond_boot.db_per_type[type].db[,2] = rowsum(struct_db_cond_boot.db_per_type[type].db[,3..(Kbar_type_temp+3)])
		for(k = 1; k <= (Kbar_type_temp-1); k++){
			struct_db_cond_boot.db_per_type[type].db[k, ((k+4)..(Kbar_type_temp+3))] = J(1, Kbar_type_temp-k, -99)
		}
		struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs = selectindex(struct_db_cond_boot.db_per_type[type].db[,2])
		struct_db_cond_boot.info_data_per_type[type].Kbar = max(struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs)
		struct_db_cond_boot.info_data_per_type[type].nb_K_positive_obs = length(struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs)
		struct_db_cond_boot.db_per_type[type].db = struct_db_cond_boot.db_per_type[type].db[(1..struct_db_cond_boot.info_data_per_type[type].Kbar), (1..(struct_db_cond_boot.info_data_per_type[type].Kbar+3))]
		struct_db_cond_boot.nb_units_studied_per_type[type] = sum(struct_db_cond_boot.db_per_type[type].db[,2], 1)			
		/* clean of db_per_type[type].db and attributes of info_data_per_type[type] (end) */
		
	} /* loop over types (end) */
	
	struct_db_cond_boot.type_probabilities = struct_db_cond_boot.nb_units_studied_per_type :/ sum(struct_db_cond_boot.nb_units_studied_per_type, 1)
	
	return(struct_db_cond_boot)
}
end

capture mata: mata drop Draw_boot_store_WXK_nb_units()
mata:
mata set matastrict on
real matrix Draw_boot_store_WXK_nb_units(real matrix store_WXK_nb_units){
	
	real matrix boot_store_WXK_nb_units
	real scalar nb_observations, index, wxk
	real colvector draws_number
	
	nb_observations = quadsum(store_WXK_nb_units[,cols(store_WXK_nb_units)], 1)
	draws_number = floor(nb_observations:*runiform(nb_observations,1):+1)
	
	boot_store_WXK_nb_units = store_WXK_nb_units
	boot_store_WXK_nb_units[,cols(boot_store_WXK_nb_units)] = J(rows(boot_store_WXK_nb_units), 1, .)
	
	index = 0
	for(wxk = 1; wxk <= rows(boot_store_WXK_nb_units); wxk++){
		if (store_WXK_nb_units[wxk, cols(store_WXK_nb_units)] > 0){
			boot_store_WXK_nb_units[wxk, cols(boot_store_WXK_nb_units)] = sum((draws_number :>= (index+1)) :* (draws_number :<= (index+store_WXK_nb_units[wxk, cols(store_WXK_nb_units)])))
			index = index + store_WXK_nb_units[wxk, cols(store_WXK_nb_units)]
		}
		else {
			boot_store_WXK_nb_units[wxk, cols(boot_store_WXK_nb_units)] = 0
		}
	}

	return(boot_store_WXK_nb_units)
}
end

/* Bootstrap at the level of the units (keeping therefore fixed the compositon of
units in terms of the different individual types (defined by individual covariates)
Used for boostrap of Rbeta method in the case conditional analyses with individual
level covariates */
capture mata: mata drop Draw_struct_db_condindi_unit()
mata:
mata set matastrict on
struct Struct_db_cond scalar Draw_struct_db_condindi_unit(struct Struct_db_cond scalar struct_db_cond){

	struct Struct_db_cond scalar struct_db_cond_boot
	real scalar index_col_nbunits, type, index_col_K, index_col_X, k, x
	
	struct_db_cond_boot = Struct_db_cond()
	struct_db_cond_boot.nb_types = struct_db_cond.nb_types
	/* si pas assez d'observations d'un certain type de telle sorte que le nombre de type a l'issue du bootstrap pourrait changer
	pas forcement de souci, juste une ponderation de 0 pour ce type non represente dans le bootstrap 
	[mais signifie que peut etre pas assez de donnees pour avoir cette finesse de distinction entre les types] */
	
	struct_db_cond_boot.store_WXK_nb_units = Draw_boot_store_WXK_nb_units(struct_db_cond.store_WXK_nb_units)
	struct_db_cond_boot.db_per_type = Struct_db(struct_db_cond_boot.nb_types)
	struct_db_cond_boot.info_data_per_type = Info_data(struct_db_cond_boot.nb_types)
	struct_db_cond_boot.nb_units_studied_per_type = J(struct_db_cond_boot.nb_types, 1, .)

	index_col_nbunits = 2*struct_db_cond_boot.nb_types+1
	for(type = 1; type <= struct_db_cond_boot.nb_types; type++){ /* loop over type (begin) */
		index_col_K = 2*type-1
		index_col_X = 2*type	
		struct_db_cond_boot.info_data_per_type[type].Kbar = max(select(struct_db_cond_boot.store_WXK_nb_units[,index_col_K], (struct_db_cond_boot.store_WXK_nb_units[,index_col_nbunits] :> 0)))
		
		/* construction of the triangle of data db for type = type */
		struct_db_cond_boot.db_per_type[type].db = J(struct_db_cond_boot.info_data_per_type[type].Kbar, struct_db_cond_boot.info_data_per_type[type].Kbar+3, 0)
		for(k = 1; k <= struct_db_cond_boot.info_data_per_type[type].Kbar; k++){
			for(x = 0; x <= k; x++){
				struct_db_cond_boot.db_per_type[type].db[k,(x+3)] = sum((select(struct_db_cond_boot.store_WXK_nb_units[,index_col_nbunits], ((struct_db_cond_boot.store_WXK_nb_units[,index_col_K] :== k):*(struct_db_cond_boot.store_WXK_nb_units[,index_col_X] :== x)))))
			}
		}
		struct_db_cond_boot.db_per_type[type].db[,2] = rowsum(struct_db_cond_boot.db_per_type[type].db, 1)
		struct_db_cond_boot.db_per_type[type].db[,1] = (1..struct_db_cond_boot.info_data_per_type[type].Kbar)'
		for(k = 1; k <= struct_db_cond_boot.info_data_per_type[type].Kbar-1; k++){
			struct_db_cond_boot.db_per_type[type].db[k, ((k+4)..(struct_db_cond_boot.info_data_per_type[type].Kbar+3))] = J(1, struct_db_cond_boot.info_data_per_type[type].Kbar- k, -99)
		}
		
		struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs = selectindex(struct_db_cond_boot.db_per_type[type].db[,2])
		struct_db_cond_boot.info_data_per_type[type].Kbar = max(struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs)
		struct_db_cond_boot.info_data_per_type[type].nb_K_positive_obs = length(struct_db_cond_boot.info_data_per_type[type].list_K_positive_obs)
		struct_db_cond_boot.nb_units_studied_per_type[type] = sum(struct_db_cond_boot.db_per_type[type].db[,2], 1)		
			
	} /* loop over type (end) */
	
	struct_db_cond_boot.type_probabilities = struct_db_cond_boot.nb_units_studied_per_type :/ sum(struct_db_cond_boot.nb_units_studied_per_type, 1)
	
	return(struct_db_cond_boot)
}
end	

/* CI_by_percentile_bootstrap() */
/* Input : vector of bootstrap realization of the statistic of interest,
level alpha of the confidence interval
The use of functions ceil and floor guarantees conservativeness of the 
procedure */
capture mata: mata drop CI_by_percentile_bootstrap()
mata:
mata set matastrict on
real rowvector CI_by_percentile_bootstrap (	real vector realized_stat_bootstrap, ///
											real scalar alpha) {
	real scalar CI_low, CI_up										
	_sort(realized_stat_bootstrap, 1)
	CI_low = max((0, realized_stat_bootstrap[max((1,floor(length(realized_stat_bootstrap) * (alpha/2))))]))
	CI_up = min((1, realized_stat_bootstrap[ceil(length(realized_stat_bootstrap) * (1-(alpha/2)))]))
	return((CI_low, CI_up))										
}
end
