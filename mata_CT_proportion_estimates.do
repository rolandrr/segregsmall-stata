version 14.2
/* this do-files creates the functions to do the proportion-based indices */

/* Input: 
- triangle of data db (the choice of db determines whether we consider or not i : K_i = 1 in the data)
- prop_minority_hat = overall proportion of minority individuals in the sample studied
- nb_individuals = total number of individuals in the sample studied
- nb_minority_individuals = total number of minority individuals in the sample studied
The function computes the indices based on the definition of Massey and Denton or 
equivalently James and Taeuber
It uses the structure of db and hopefully vectorization of the code to avoid loop
over units
*/


/* Duncan */
capture mata: mata drop Prop_Duncan()
mata:
mata set matastrict on
real scalar Prop_Duncan(real matrix db, real scalar prop_minority_hat, real scalar nb_individuals) {

	real scalar k
	
	real scalar matrix_proportions
	matrix_proportions = J(rows(db), rows(db)+1, 0)
	for(k = 1; k <= rows(db); k++) {
		matrix_proportions[k..rows(db),k+1] = (k :/ (k..rows(db)))'
	}

	real matrix matrix_coefficients
	matrix_coefficients = abs(matrix_proportions :- prop_minority_hat) :* db[,1]
	for(k = 1; k <= rows(matrix_coefficients)-1; k++) {
		matrix_coefficients[k, (2+k)..cols(matrix_coefficients)] = J(1, cols(matrix_coefficients)-(k+1), 0)
	}
		
	return(quadsum(matrix_coefficients :* db[.,3..cols(db)], 1) / (2*nb_individuals*prop_minority_hat*(1-prop_minority_hat)))
}
end

capture mata: mata drop Prop_Duncan_from_struct_db_unc()
mata:
mata set matastrict on
real scalar Prop_Duncan_from_struct_db_unc(struct Struct_db_uncond scalar struct_db_uncond){
	return(Prop_Duncan(struct_db_uncond.db, struct_db_uncond.info_data.prop_minority_hat, struct_db_uncond.info_data.nb_individuals))
}
end


/* Theil */
capture mata: mata drop compute_xlog2x()
mata:
mata set matastrict on
real matrix compute_xlog2x(real matrix x) {
	return(x :* (ln(x :+ (x :== 0)) :/ ln(2)))
}
end

capture mata: mata drop Prop_Theil()
mata:
mata set matastrict on
real scalar Prop_Theil(real matrix db, real scalar prop_minority_hat, real scalar nb_individuals) {

	real scalar k
	
	real scalar matrix_proportions
	matrix_proportions = J(rows(db), rows(db)+1, 0)
	for(k = 1; k <= rows(db); k++) {
		matrix_proportions[k..rows(db),k+1] = (k :/ (k..rows(db)))'
	}

	real scalar E /* Entropy of total population (in fact it is minus entropy but
	it will be similar for the entropy by unit so the sign cancels out */
	E = compute_xlog2x(prop_minority_hat) + compute_xlog2x(1-prop_minority_hat)
	
	real matrix matrix_coefficients
	matrix_coefficients = (E :- (compute_xlog2x(matrix_proportions) + compute_xlog2x(1 :- matrix_proportions))) :* db[,1]
	for(k = 1; k <= rows(matrix_coefficients)-1; k++) {
		matrix_coefficients[k, (2+k)..cols(matrix_coefficients)] = J(1, cols(matrix_coefficients)-(k+1), 0)
	}
	
	return(quadsum(matrix_coefficients :* db[.,3..cols(db)], 1) / (nb_individuals*E))
}
end

capture mata: mata drop Prop_Theil_from_struct_db_unc()
mata:
mata set matastrict on
real scalar Prop_Theil_from_struct_db_unc(struct Struct_db_uncond scalar struct_db_uncond){
	return(Prop_Theil(struct_db_uncond.db, struct_db_uncond.info_data.prop_minority_hat, struct_db_uncond.info_data.nb_individuals))
}
end


/* Atkinson */
capture mata: mata drop Prop_Atkinson()
mata:
mata set matastrict on
real scalar Prop_Atkinson(real matrix db, real scalar prop_minority_hat, real scalar nb_individuals, real scalar b) {

	real scalar k
	
	real scalar matrix_proportions
	matrix_proportions = J(rows(db), rows(db)+1, 0)
	for(k = 1; k <= rows(db); k++) {
		matrix_proportions[k..rows(db),k+1] = (k :/ (k..rows(db)))'
	}

	real matrix matrix_coefficients
	matrix_coefficients = (matrix_proportions :^ b) :* ((1 :- matrix_proportions):^(1-b)) :* db[,1]
	for(k = 1; k <= rows(matrix_coefficients)-1; k++) {
		matrix_coefficients[k, (2+k)..cols(matrix_coefficients)] = J(1, cols(matrix_coefficients)-(k+1), 0)
	}
		
	return(1 - ((prop_minority_hat/(1-prop_minority_hat)) * ((quadsum(matrix_coefficients :* db[.,3..cols(db)], 1) / (nb_individuals*prop_minority_hat))^(1/(1-b)))))
}
end
	
capture mata: mata drop Prop_Atkinson_from_struct_db_unc()
mata:
mata set matastrict on
real scalar Prop_Atkinson_from_struct_db_unc(struct Struct_db_uncond scalar struct_db_uncond, real scalar b){
	return(Prop_Atkinson(struct_db_uncond.db, struct_db_uncond.info_data.prop_minority_hat, struct_db_uncond.info_data.nb_individuals, b))
}
end	
	
	
/* Coworker */	
capture mata: mata drop Prop_Coworker()
mata:
mata set matastrict on
real scalar Prop_Coworker(real matrix db, real scalar prop_minority_hat, real scalar nb_individuals) {

	real scalar k
	
	real scalar matrix_proportions
	matrix_proportions = J(rows(db), rows(db)+1, 0)
	for(k = 1; k <= rows(db); k++) {
		matrix_proportions[k..rows(db),k+1] = (k :/ (k..rows(db)))'
	}

	real matrix matrix_coefficients
	matrix_coefficients = ((matrix_proportions :- prop_minority_hat):^2) :* db[,1]
	for(k = 1; k <= rows(matrix_coefficients)-1; k++) {
		matrix_coefficients[k, (2+k)..cols(matrix_coefficients)] = J(1, cols(matrix_coefficients)-(k+1), 0)
	}
	
	return(quadsum(matrix_coefficients :* db[.,3..cols(db)], 1) / (nb_individuals*prop_minority_hat*(1-prop_minority_hat)))
}
end

capture mata: mata drop Prop_Coworker_from_struct_db_unc()
mata:
mata set matastrict on
real scalar Prop_Coworker_from_struct_db_unc(struct Struct_db_uncond scalar struct_db_uncond){
	return(Prop_Coworker(struct_db_uncond.db, struct_db_uncond.info_data.prop_minority_hat, struct_db_uncond.info_data.nb_individuals))
}
end


/* Gini */
capture mata: mata drop Prop_Gini()
mata:
mata set matastrict on
real scalar Prop_Gini(real matrix db, real scalar prop_minority_hat, real scalar nb_individuals) {
	
	real matrix matrix_coefficients_i, matrix_coefficients_j
	real scalar Kbar
	real scalar k_i, x_i, k_j, x_j

	Kbar = rows(db)
	
	matrix_coefficients_i = J(Kbar, Kbar+1, 0)
	for(k_i = 1; k_i <= Kbar; k_i++) {
		for(x_i = 0; x_i <= k_i; x_i++){
			matrix_coefficients_j = J(Kbar, Kbar+1, 0)
			for(k_j = 1; k_j <= Kbar; k_j++){
				for(x_j = 0; x_j <= k_j; x_j++){
					matrix_coefficients_j[k_j, x_j+1] = k_i * k_j * abs(x_i/k_i - x_j/k_j)
				}
			}
			matrix_coefficients_i[k_i, x_i+1] = sum(matrix_coefficients_j :* db[,(3..(Kbar+3))],1) / nb_individuals /* division here by nb_individuals to try to avoid numerical overflow */
		}
	}

	return(sum(matrix_coefficients_i :* db[,(3..(Kbar+3))],1) / (2*nb_individuals*prop_minority_hat*(1-prop_minority_hat)))
}
end

capture mata: mata drop Prop_Gini_from_struct_db_unc()
mata:
mata set matastrict on
real scalar Prop_Gini_from_struct_db_unc(struct Struct_db_uncond scalar struct_db_uncond){
	return(Prop_Gini(struct_db_uncond.db, struct_db_uncond.info_data.prop_minority_hat, struct_db_uncond.info_data.nb_individuals))
}
end
