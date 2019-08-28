version 14.2
/* functions to manipulate data */

/******************************************************************************/
/* Different functions to compute quantities of interest from dataset db */
/******************************************************************************/

capture mata: mata drop Get_nb_units_studied()
mata:
mata set matastrict on
real scalar Get_nb_units_studied(real matrix db){
	return(quadsum(db[,2],1))
}
end

/* Compute_Pi_from_scratch() : from a given dataset db 
the choice of db determines whether we consider or not i : K_i = 1 in the data
Pi = overall proportion of minority individuals in the sample
It is an estimate in case of sample data */
capture mata: mata drop Compute_Pi_from_scratch()
mata:
mata set matastrict on
real scalar Compute_Pi_from_scratch(real matrix db) {
	return(Compute_S_X(db) / Compute_S(db))
}
end						

/* S = total number of individuals in the sample */
capture mata: mata drop Compute_S()
mata:
mata set matastrict on
real scalar Compute_S(real matrix db) {
	return(db[,1]' * db[,2])
}
end

/* S_X = total number of minority individuals in the sample */
capture mata: mata drop Compute_S_X()
mata:
mata set matastrict on
real scalar Compute_S_X(real matrix db) {
	real scalar S_X
	real scalar x
	S_X = 0
	for(x = 4; x <= cols(db); x++) {
		S_X = S_X + (x-3) * quadsum(db[(x-3)..rows(db),x],1)
	}
	return(S_X)		
}
end

/* Compute_Pi() : use S and S_X already computed in struct Info_data  */
capture mata: mata drop Compute_Pi()
mata:
mata set matastrict on
real scalar Compute_Pi(struct Info_data scalar info_data){
	return(info_data.nb_minority_individuals / info_data.nb_individuals)
}
end

/* structure of the object store_ZK_nb_units
matrix with three columns
each row is a type x size of a unit 
first column indicates the type (Z)
second column indicates the size of the unit (K)
third column indicates the number of such units of type = Z and size = K in the data */
capture mata: mata drop Cons_store_ZK_nb_units()
mata:
mata set matastrict on
real matrix Cons_store_ZK_nb_units(struct Struct_db_cond scalar struct_db_cond){
	real matrix store_ZK_nb_units
	real scalar type
	for(type = 1; type <= struct_db_cond.nb_types; type++){
		if(type == 1) store_ZK_nb_units = (J(rows(struct_db_cond.db_per_type[type].db),1,type),struct_db_cond.db_per_type[type].db[,(1,2)])
		else store_ZK_nb_units = (store_ZK_nb_units \ (J(rows(struct_db_cond.db_per_type[type].db),1,type),struct_db_cond.db_per_type[type].db[,(1,2)]))
	}
	return(store_ZK_nb_units)
}
end

capture mata: mata drop Cons_type_frequencies_unit()
mata:
mata set matastrict on
real colvector Cons_type_frequencies_unit(struct Struct_db_cond scalar struct_db_cond){
	real colvector type_frequencies_unit
	real scalar type
	type_frequencies_unit = J(struct_db_cond.nb_types, 1, .)
	for(type = 1; type <= struct_db_cond.nb_types; type++){
		type_frequencies_unit[type] = quadsum(select(struct_db_cond.store_ZK_nb_units[,3], struct_db_cond.store_ZK_nb_units[,1] :== type),1)
	}
	return(type_frequencies_unit)
}
end

capture mata: mata drop Cons_type_probabilities_unit()
mata:
mata set matastrict on
real colvector Cons_type_probabilities_unit(struct Struct_db_cond scalar struct_db_cond){
	return(struct_db_cond.type_frequencies :/ struct_db_cond.info_data_uncond.nb_units_studied)
}
end

/* structure of the object store_WK_nb_units
matrix 
each row is a given composition of units formed by individuals from different types
1st column : number of individuals of type 1 in this unit
2nd column : number of individuals of type 2 in this unit
....
nb_type-th column :  number of individuals of type nb_type in this unit
nb_type+1-th column : number of units with such composition 
nb_type+2-th column : number of individuals in the unit with such composition,
that is simply equal to the sum of the row from 1st to nb_type-th column

Example nb_type = 3

store_WK_nb_units : 
0 2 7 5 9
1 2 1 2 4
etc....

first row means : 
in the data, there are 
* 5 units that contain 
- 0 individuals of type 1
- 2 individuals of type 2
- 7 individuals of type 3
- and thus 9 individuals in total (aggregating types)
* 2 units that contain 
- 1 individuals of type 1
- 2 individuals of type 2
- 1 individuals of type 3
- and thus 4 individuals in total (aggregating types)

option withsingleton:
without (default): suppress units with only 1 individual
with: keep those units

option excludingsingletonpertype
with: do not take into accounts cells type x unit with only 1 individual,
as if set to 0 
hence it changes the number of units of a given composition

for instance if without option excludingsingletonpertype, store_WK_nb_units : 
....
0 0 2 2 2
....
0 1 2 2 3
....
1 0 2 2 3
....
1 1 2 1 4
....

those four lines become in store_WK_nb_units with option excludingsingletonpertype
only one line : 
0 0 2 7 2

after in the analysis and computation of index per type (before aggregating), 
the link between the identifier of the units and the type of individuals
belonging to this units is forgotten

*/

/* structure of the object store_WXK_nb_units

store_WXK_nb_units is similar to store_WK_nb_units
except that it distinguishes the composition of the unit both in terms 
of individual types and in terms of minority or majority individuals
columns are to be considered by pairs, except the last one
1st column: number of individuals of type 1 minority or majority in this unit
2nd column: number of individuals of type 1 minority in this unit
...

(2*nb_type-1)-th column: number of individuals of type nb_type minority or majority in this unit
(2*nb_type)-th column:number of individuals of type nb_type minority in this unit
nb_type+1-th column : number of units with such composition 

It is used to perform the boostrap at the unit-level (Rbeta method) keeping 
constant the composition of the units for conditional analyses with 
individual-level covariates 
*/

capture mata: mata drop Cons_store_WK_nb_units()
mata:
mata set matastrict on
real matrix Cons_store_WK_nb_units(real matrix store_WK_nb_units_withoutK){
	real colvector K
	K = quadrowsum(store_WK_nb_units_withoutK[,(1..(cols(store_WK_nb_units_withoutK)-1))],1)
	return((store_WK_nb_units_withoutK, K))
}
end

capture mata: mata drop Cons_nb_units_per_type_ind()
mata:
mata set matastrict on
real colvector Cons_nb_units_per_type_ind(struct Struct_db_cond scalar struct_db_cond){
	real colvector nb_units_per_type
	real scalar type
	nb_units_per_type = J(struct_db_cond.nb_types, 1, .)
	for(type = 1; type <= struct_db_cond.nb_types; type++){
		nb_units_per_type[type] = quadsum(struct_db_cond.db_per_type[type].db[,2],1)
	}
	return(nb_units_per_type)
}
end

capture mata: mata drop Cons_type_frequencies_ind()
mata:
mata set matastrict on
real colvector Cons_type_frequencies_ind(struct Struct_db_cond scalar struct_db_cond){
	real colvector type_frequencies_ind
	real scalar type
	type_frequencies_ind = J(struct_db_cond.nb_types, 1, .)
	for(type = 1; type <= struct_db_cond.nb_types; type++){
		type_frequencies_ind[type] = Compute_S(struct_db_cond.db_per_type[type].db)
	}
	return(type_frequencies_ind)
}
end

capture mata: mata drop Cons_type_probabilities_ind()
mata:
mata set matastrict on
real colvector Cons_type_probabilities_ind(struct Struct_db_cond scalar struct_db_cond){
	return((struct_db_cond.type_frequencies :/ quadsum(struct_db_cond.type_frequencies,1)))
}
end

capture mata: mata drop Cons_summary_info_data_per_type()
mata:
mata set matastrict on
real matrix Cons_summary_info_data_per_type(struct Struct_db_cond scalar struct_db_cond){
	real matrix summary_info_data_per_type
	real scalar type
	summary_info_data_per_type = J(struct_db_cond.nb_types, 9, .)
	for(type = 1; type <= struct_db_cond.nb_types; type++){
		summary_info_data_per_type[type,1] = struct_db_cond.info_data_per_type[type].I_withsingleton
		summary_info_data_per_type[type,2] = struct_db_cond.info_data_per_type[type].nb_units_total
		summary_info_data_per_type[type,3] = struct_db_cond.info_data_per_type[type].nb_units_singleton
		summary_info_data_per_type[type,4] = struct_db_cond.info_data_per_type[type].nb_units_studied
		summary_info_data_per_type[type,5] = struct_db_cond.info_data_per_type[type].nb_individuals
		summary_info_data_per_type[type,6] = struct_db_cond.info_data_per_type[type].nb_minority_individuals
		summary_info_data_per_type[type,7] = struct_db_cond.info_data_per_type[type].prop_minority_hat
		summary_info_data_per_type[type,8] = struct_db_cond.info_data_per_type[type].Kbar
		summary_info_data_per_type[type,9] = struct_db_cond.info_data_per_type[type].nb_K_positive_obs
	}
	return(summary_info_data_per_type)
}
end

/* for conditional analyses with individual-level covariates
nb_individuals studied may be different with the option excluding, hence this
correction (Corr) 
NB: it uses struct_db_cond.summary_info_data_per_type */
capture mata: mata drop Corr_info_data_uncond_condind()
mata:
mata set matastrict on
struct Info_data scalar Corr_info_data_uncond_condind(struct Struct_db_cond scalar struct_db_cond){

	struct Info_data scalar info_data_uncond
	real scalar type
	
	info_data_uncond = Info_data()
	info_data_uncond = struct_db_cond.info_data_uncond 	/* take the current info_data_uncond and just correct some parts */
	
	info_data_uncond.nb_individuals = sum(struct_db_cond.summary_info_data_per_type[,5], 1)
	info_data_uncond.nb_minority_individuals = sum(struct_db_cond.summary_info_data_per_type[,6], 1)
	info_data_uncond.prop_minority_hat = info_data_uncond.nb_minority_individuals / info_data_uncond.nb_individuals

	return(info_data_uncond)
}
end	
/* OLD: without using struct_db_cond.summary_info_data_per_type
info_data_uncond.nb_individuals = 0
info_data_uncond.nb_minority_individuals = 0
for(type = 1; type <= struct_db_cond.nb_types; type++){ /* loop over types (begin) */
	info_data_uncond.nb_individuals = info_data_uncond.nb_individuals + struct_db_cond.info_data_per_type[type].nb_individuals
	info_data_uncond.nb_minority_individuals = info_data_uncond.nb_minority_individuals + struct_db_cond.info_data_per_type[type].nb_minority_individuals
} /* loop over types (end) */
info_data_uncond.prop_minority_hat = info_data_uncond.nb_minority_individuals / info_data_uncond.nb_individuals
*/

/* NB: it uses struct_db_cond.summary_info_data_per_type */
capture mata: mata drop Cons_nb_singleton_cells_all_type()
mata:
mata set matastrict on
real scalar Cons_nb_singleton_cells_all_type(struct Struct_db_cond scalar struct_db_cond){
	return(sum(struct_db_cond.summary_info_data_per_type[,3], 1))
}
end
