# Do folder

## creation_mlib.do
Is the do file used to create the mlib

## The other do-files
Are used to create the mata 

### Structures for information about data and results 

- mata_output
- mata_tools_database
- db_from_unit_level
- db_from_indiv_level
- db_from_unit_level_cond_unit
- db_from_indiv_level_cond_unit
- db_from_indiv_level_cond_indiv

### DHR method 

#### case without restriction as regards dependency beetween K and p 
- mata_DHR_ML_oneK
- mata_DHR_ML_allK
- mata_DHR_principal_representation
- mata_DHR_optibycr_A21
- mata_DHR_bounds

#### case with assumed independence between K and p 
- mata_DHR_ML_Kpind
- mata_DHR_bounds_Kpind

#### inference and conditional case 
- mata_DHR_tools_ci_Kpind
- mata_DHR_tools_ci
- mata_DHR_conditional

#### test binomial 
- mata_DHR_tools_test_binomial

#### Stata output functions (wrap-up) 
- mata_DHR_stata_output

### R beta method 
- mata_R_ML
- mata_R_estimates
- mata_R_tools_ci
- mata_R_deltamethod
- mata_R_conditional
- mata_R_stata_output

### CT method 
- mata_CT_proportion_estimates
- mata_CT_correction_ct
- mata_CT_conditional
- mata_CT_stata_output