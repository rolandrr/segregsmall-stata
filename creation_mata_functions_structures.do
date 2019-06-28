/* NB: the order of the files matter as some structures use previous ones */

/**************************/
/* data base manipulation */
/**************************/


/*******/
/* DHR */
/*******/
/* case without restriction as regards dependency beetween K and p */
do mata_DHR_ML_oneK
do mata_DHR_ML_allK
do mata_DHR_principal_representation
do mata_DHR_optibycr_A21
do mata_DHR_bounds
/* case with assumed independence between K and p */
do mata_DHR_ML_Kpind
do mata_DHR_bounds_Kpind
/* inference and conditional case */
do mata_DHR_tools_ci_Kpind
do mata_DHR_tools_ci
do mata_DHR_conditional
/* test binomial */
do mata_DHR_tools_test_binomial
/* Stata output functions (wrap-up) */
do mata_DHR_stata_output

/*********/
/* Rbeta */
/*********/


/*********************/
/* Carrington-Troske */
/*********************/

