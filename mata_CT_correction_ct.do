version 14.2
/* Functions to do CT-correction */

/* Construct_db_random_allocation() */
/* construct a database under the hypothesis of randomness
cf. CFC page 633 and Appendix :
- given the number of units n
- given the size K_i of each unit
- given the sample overall minority proportion Pi
- assume that each unit of size K_i is a random sample from the total population
We will use the binomial approximation of the hypergeometric distribution
following CFC remarks (cf. mathematical appendix, CFC page 635
and CT remarks (first paragraph of section 1, page 403) */
capture mata: mata drop Construct_db_random_allocation()
mata:
mata set matastrict on
real matrix Construct_db_random_allocation(	real matrix db, ///
											real scalar Kbar, ///
											real scalar prop_minority_hat){
	real matrix db_random
	real scalar k
	
	db_random = J(Kbar, Kbar+3, -99)
	db_random[,(1,2)] = db[.,(1,2)]
	for (k = 1; k <= Kbar; k++) {
		if (db_random[k,2] == 0) 	db_random[k,3..(k+3)] = J(1, k+1, 0)
		else 						db_random[k,3..(k+3)] = count_rbinomial(rbinomial(1, db_random[k,2], k, prop_minority_hat), k)
	}

	/* in case of doubt about validity of binomial approximation */
	/*if (Kbar >= 0.1 * S) {
		errprintf("Warning : binomial approximation of hypergeometric doubtful for the Carrington-Troske correction")
	}*/
	
	return(db_random)
}
end								

/* try to do something faster but same computation time approximately
the main cost is the count_rbinomial and drawing of the binomial,
not the loop over unit size k */
/* OLD
capture mata: mata drop Construct_db_random_allocation_2()
mata:
mata set matastrict on
real matrix Construct_db_random_allocation_2(	real matrix db, ///
												real colvector list_K_positive_obs, ///
												real scalar nb_K_positive_obs, ///
												real scalar prop_minority_hat){
	real matrix db_random
	real scalar index_k, k
	
	db_random = db
	for(index_k = 1; index_k <= nb_K_positive_obs; index_k++) {
		k = list_K_positive_obs[index_k]
		db_random[k,3..(k+3)] = count_rbinomial(rbinomial(1, db_random[k,2], k, prop_minority_hat), k)
	}

	return(db_random)
}
end								
*/

capture mata: mata drop count_rbinomial()
mata:
mata set matastrict on
real rowvector count_rbinomial(	real rowvector draws, ///
								real scalar k) {
	real scalar i
	real rowvector count_rbinomial
	count_rbinomial = J(1, k+1, .)
	for (i = 0; i <= k; i++) {
		count_rbinomial[i+1] = sum(draws :== i)
	}
	return(count_rbinomial)
}
end								

/* Test of another method to be faster,
but turns out to be longer */
/* OLD:
capture mata: mata drop count_rbinomial_v2()
mata:
mata set matastrict on
real rowvector count_rbinomial_v2(	real rowvector draws, ///
									real scalar k) {
	real rowvector accumulate_Nk
	real scalar i
	accumulate_Nk = J(1, k+1, 0)
	for(i = 1; i <= length(draws); i++){
			accumulate_Nk[draws[i]+1] = accumulate_Nk[draws[i]+1] + 1
	}
	return(accumulate_Nk)			
}
end
*/

/* Compute_CT_correction() */
/* compute CT correction for a given proportion based segregation index 
realized in data: index_realized_in_data
and a given estimate of index under random allocation: index_random_allo */
capture mata: mata drop Compute_CT_correction()
mata:
mata set matastrict on
real scalar Compute_CT_correction(	real scalar index_realized_in_data, ///
									real scalar index_random_allo) {
							
	if (index_realized_in_data >= index_random_allo) {
		return((index_realized_in_data - index_random_allo) / (1 - index_random_allo))
	}
	else {
		return((index_realized_in_data - index_random_allo) / index_random_allo)
	}
}
end		

/* Compute_CFC_standard_score() */
/* compute Cortese Falk and Cohen standard score for a given proportion based 
segregation index realized in data: index_realized_in_data
and mean and standard deviation of indices under random allocation */
capture mata: mata drop Compute_CFC_standard_score()
mata:
mata set matastrict on
real scalar Compute_CFC_standard_score(	real scalar index_realized_in_data, ///
										real scalar index_random_allo_mean, ///
										real scalar index_random_allo_sd){										
	return((index_realized_in_data - index_random_allo_mean) / index_random_allo_sd)
}
end					
