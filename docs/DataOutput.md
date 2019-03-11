#   Data manipulation and output

## mata_output
Info_data()
Struct_db_uncond()
Struct_db()
Struct_db_cond()
Info_estimation_inference()
Results_uncond()
Struct_eci()
Results_cond()
Results_eci()
_View_Info_data()
_View_Struct_db_uncond()
_View_Struct_db_cond()
_View_Info_estimation_inference()
_View_Results_uncond()
_View_Results_cond()


## mata_tools_database
Get_nb_units_studied()

Compute_Pi_from_scratch()
Compute_S()
Compute_S_X()
Compute_Pi()

Cons_store_ZK_nb_units()
Cons_type_frequencies_unit()
Cons_type_probabilities_unit()

Cons_store_WK_nb_units()
Cons_nb_units_per_type_ind()
Cons_type_frequencies_ind()
Cons_type_probabilities_ind()

Cons_summary_info_data_per_type()
Corr_info_data_uncond_condind()
Cons_nb_singleton_cells_all_type()


## db_from_unit_level.ado
Check for missing values in K or X, if any exit
Check for integer >= 1 in K, if any exit
Check for integer >= 0 in X, if any exit
Check that for each unit, X <= K, if any exit
=> struct_db_uncond


## db_from_indiv_level.ado
Check for missing values in id_unit or I_minority, if any exit
Check that I_minority is a binary indicator 0-1 variable, if not exit
Check that id_unit uniquely identifies the units, if not exit
It uses db_from_unit_level
=> struct_db_uncond


## db_from_unit_level_cond_unit.ado
Check for missing values in K, X, or Z, if any exit
Check for integer >= 1 in K, if any exit
Check for integer >= 0 in X, if any exit
Check that for each unit, X <= K, if any exit
(même si après utilise db_from_unit_level, mieux de tester et d'arrêter dès le début le cas échéant)
Check that Z is categorical variable in {1,2,3,...} with at least two modalities
It uses db_from_unit_level
=> struct_db_cond


## db_from_indiv_level_cond_unit.ado
Check for missing values in id_unit, I_minority, or Z, if any exit
Check that I_minority is a binary indicator 0-1 variable, if not exit
Check that Z does not vary within units, in this case exit
Check that id_unit uniquely identifies the units, if not exit
It uses db_from_unit_level_cond_unit
=> struct_db_cond


## db_from_indiv_level_cond_indiv.ado
Check for missing values in id_unit, I_minority, or W, if any exit
Check that I_minority is a binary indicator 0-1 variable, if not exit
Check that W is categorical variable in {1,2,3,...} with at least two modalities
It uses db_from_ind_level for aggregated part of struct_db_cond
Check that id_unit uniquely identifies the units, if not exit
and then db_from_unit_level
=> struct_db_cond
