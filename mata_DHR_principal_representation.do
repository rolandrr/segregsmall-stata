version 14.2
/******************************************************************************/
/* This do-file creates structures and functions for defining principal
representations
This is used in the bounds with assumption 2.2. */
/******************************************************************************/

/* structure : principal_representations
Constructor : Principal_representations()
compute lower and upper principal representations */

capture mata: mata drop principal_representations()
mata:
mata set matastrict on
struct principal_representations
{
	real matrix lower, upper
}
end
/* as in the rest of the code, a discrete distribution is represented by a 
2 row matrix, first row = support points, second row = associated masses,
the number of columns being the number of roots aka support points */

capture mata: mata drop _View_principal_representations()
mata:
mata set matastrict on
void _View_principal_representations(struct principal_representations scalar pr) {
	printf("Structure principal_representations - 2 elements \n ")
	printf("lower = \n ")
	pr.lower
	printf("upper = \n ")
	pr.upper
}
end

capture mata: mata drop Get_principal_representation()
mata:
mata set matastrict on
real matrix Get_principal_representation(	struct principal_representations scalar pr, ///
											real scalar chosen_pr) {
	/* upper principal representation for any true value of chosen_pr in particular for 1 => 1 
	as presence of this root determined if upper or not; lower representation for chosen_pr = 0 */
	if (chosen_pr) {
		return(pr.upper)
	}
	else {
		return(pr.lower)
	}								
}
end											

capture mata: mata drop Principal_representations()
mata:
mata set matastrict on
struct principal_representations scalar Principal_representations(real colvector mu) {

	struct principal_representations scalar pr
	real scalar K; K = length(mu)
	real scalar L; L = floor((K+1)/2)
	real scalar I_K_odd; I_K_odd = K - 2*floor(K/2)
	real matrix A, B, C 
	real colvector a_low, a_up, x_low, x_up


	/* case K odd */
	if (I_K_odd) {
	
		/* special case K = 1 NB: mu[1] = m1 different de 0 et 1 a priori 
		sinon cas contraint et pas d'interet a rechercher des pr */
		if (K == 1) {
		
			pr.lower = (mu \ 1)
			pr.upper = ((0,1) \ ((1 - mu), mu))
		
		}
		/* case K odd > 1 */
		else {
		
			A = construct_A(L, mu)
			C = construct_C(L-1, mu)
		
			/* lower principal representation */
			a_low = qrsolve(A, -mu[L..K])
			x_low = polyroots((a_low \ 1)')
			/* upper principal representation */
			a_up = qrsolve(C, mu[L..K-1] :- mu[L+1..K])
			x_up = (0, polyroots((a_up \ 1)'), 1)
			
		}
	}
	/* case K even */
	else {
	
		A = construct_A(L, mu)
		B = construct_B(L, mu)
		
		/* lower principal representation */
		a_low = qrsolve(B, -mu[L+1..K])
		x_low = (0, polyroots((a_low \ 1)'))
		
		/* upper principal representation */
		a_up = qrsolve(B-A, mu[L..K-1] :- mu[L+1..K])
		x_up = (polyroots((a_up \ 1)'), 1)
	
	}

	if (K > 1) {
		if ((isrealvalues(x_low)) & (isrealvalues(x_up)) & (length(x_low) >= 1) & (length(x_up) >= 1) &  (!(sum(x_low :== .) + sum(x_up :== .)))) {
			pr.lower = clean_xy(get_xy_from_x_and_mu(x_low, mu))
			pr.upper = clean_xy(get_xy_from_x_and_mu(x_up, mu))
		}
		else {
			errprintf("WARNING : no principal representation support points complex or . or nothing - function Principal_representations() \")
			pr.lower = (. \ .)
			pr.upper = (. \ .)
		}
	}
	
	return(pr)	
}
end


/* get_xy_from_x_and_mu()
Input: a rowvector of support points and the associated vector of moments -
starting from moment of order 1 and whose size is at least number of support 
points - 1 */
capture mata: mata drop get_xy_from_x_and_mu()
mata:
mata set matastrict on
real matrix get_xy_from_x_and_mu(numeric rowvector x, numeric colvector mu) {
	
	numeric matrix V
	real vector index_y_positive
	numeric colvector y
	real scalar crit_continue_while, nb_inversions
	real scalar threshold_positive_y; 	threshold_positive_y = 10^(-8)
	real scalar nb_inversions_max;		nb_inversions_max = 7
	real scalar threshold_x;			threshold_x = 10^(-6)
	real scalar tolerance_sum_y;		tolerance_sum_y = 10^(-9)

/* ajout pour opti, a voir pour principal representations */
	if (sum(x :== .)) {
		return((0\1))
	}
	
	if (isrealvalues(x)) {
		x = Re(x)
	}
	else {
		errprintf("EXIT: x contains non real values - function get_xy_from_x_and_mu()")
		exit()
	}
	
	
	/* keep only support points in [0,1] and round to 0 or 1 if happens */
	x = x[selectindex((x :>=0) :* (x :<= 1))]	
	if (sum(x :<= threshold_x)) x[selectindex(x :<= threshold_x)] = J(1, length(selectindex(x :<= threshold_x)), 0)
	if (sum(x :>= 1 - threshold_x)) x[selectindex(x :>= 1 - threshold_x)] = J(1, length(selectindex(x :>= 1 - threshold_x)), 1)

	/* select if same value after rounding at 0 or 1 */
	if(sum(x :== 0) > 1){
		x = x[(selectindex(x:>0), selectindex((x:==0))[1])]
	}
	if(sum(x :== 1) > 1){
		x = x[(selectindex(x:<1), selectindex((x:==1))[1])]
	}

	if (length(x) == 1) {
		y = 1
	}
	else { /* case length(x) > 1 */
	
		V = Vandermonde(x')'
		y = lusolve(V, (1 \ mu[1..(length(x)-1)]))
				
		if (isrealvalues(y)) { /* case y contains only real values */
			y = Re(y)
		
			if (sum(y :== .)) {
			
				errprintf("WARNING: lusolve of get_xy_from_x_and_mu() returns . - try with qrsolve()")
				y = qrsolve(V, (1 \ mu[1..(length(x)-1)]))
				
				if (sum(y :== .)) {
					errprintf("EXIT: qrsolve of get_xy_from_x_and_mu() return . - exit() the function")
					exit()
				}
				
				if (isrealvalues(y)) {
					y = Re(y)
				}
				else {
					errprintf("EXIT: y contains non real values - function get_xy_from_x_and_mu()")
					exit()
				}
			}

			crit_continue_while = 1
			nb_inversions = 0
			while ((crit_continue_while) & (nb_inversions <= nb_inversions_max)) {
			
				nb_inversions = nb_inversions + 1
				
				if (sum(y :< threshold_positive_y)) {
			
					index_y_positive = selectindex(y :>= threshold_positive_y)
					x = x[index_y_positive]
					
					if (length(x) == 1) {
						y = 1
						crit_continue_while = 0
					}
					
					else {

						V = Vandermonde(x')'
						y = lusolve(V, (1 \ mu[1..(length(x)-1)]))
						
						if (isrealvalues(y)) {
							y = Re(y)
						}
						else {
							errprintf("EXIT: y contains non real values - function get_xy_from_x_and_mu()")
							exit()
						}
						
						if (sum(y :== .)) {
							errprintf("WARNING: lusolve of get_xy_from_x_and_mu() returns . - try with qrsolve()")
							y = qrsolve(V, (1 \ mu[1..(length(x)-1)]))
							
							if (isrealvalues(y)) {
							y = Re(y)
							}
							else {
								errprintf("EXIT: y contains non real values - function get_xy_from_x_and_mu()")
								exit()
							}
						}
						crit_continue_while = 1
					}
				}
				else {
					crit_continue_while = 0
				}
			} /* end of boucle while */
					
			/* Ajout d'une tolerance pour bien un vecteur de probabilites */
			if (reldif(sum(y),1) > tolerance_sum_y) {
				errprintf("EXIT: y does not sum to 1 - function get_xy_from_x_and_mu()")
				exit()
			}
			
			
			if (((sum(y :< 0))) | ((sum(y :> 1)))) {
				errprintf("EXIT: y is not an adequate vector of probabilities (< 0, > 1) - function get_xy_from_x_and_mu()")
				exit()
			}
			/*
			if ((sum(x :< 0)) + (sum(x :> 1))) {
				if (sum(y[selectindex((x :< 0) + (x :> 1))]) > threshold_positive_y) {
					errprintf("EXIT : positive mass given to support points x outside [0,1] - function get_xy_from_x_and_mu()")
					exit()
				}
			}
			*/
		} /* end case : y contains only real values */
		else {
			errprintf("EXIT: y contains non real values - function get_xy_from_x_and_mu()")
			exit()
		}
	} /* end case length(x) > 1 */
		
	return(x \ (y'))
}
end
