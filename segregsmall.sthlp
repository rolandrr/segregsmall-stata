{smcl}
{* *! version 1.1.0 28aug2019}{...}
{viewerjumpto "Syntax" "segregsmall##syntax"}{...}
{viewerjumpto "Description" "segregsmall##description"}{...}
{viewerjumpto "Options" "segregsmall##options"}{...}
{viewerjumpto "Remarks" "segregsmall##remarks"}{...}
{viewerjumpto "Examples" "segregsmall##examples"}{...}
{viewerjumpto "Stored results" "segregsmall##results"}{...}
{viewerjumpto "Authors" "segregsmall##authors"}{...}
{viewerjumpto "References" "segregsmall##references"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{bf:segregsmall} {hline 2}}Estimation of segregation indices in small-unit settings{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:segregsmall}
{varlist}
{ifin}
{cmd:,} {opt f:ormat}{cmd:(}{it:format}{cmd:)}
{opt m:ethod}{cmd:(}{it:method}{cmd:)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{opt f:ormat(format)}}indicates database's format; {it:format} must be either {it:unit} or {it:indiv}{p_end}
{synopt :{opt m:ethod(method)}}specifies the method used to compute indices; {it:method} must be either {it:np}, {it:beta}, or {it:ct}{p_end}

{syntab:Conditional}
{synopt :{opt condi:tional(conditional)}}performs conditional segregation analysis; {it:conditional} must be either {it:unit} or {it:indiv} depending on the level of the characteristics{p_end}

{syntab:Population}
{synopt :{opt with:single}}includes in the analysis units with only one individual; by default, they are dropped as uninformative about segregation{p_end}
{synopt :{opt exclu:dingsinglepertype}}excludes from the analysis cells (unit x type) with only one individual; by default, they are included; only available in conditional analysis with individual characteristics{p_end}

{syntab:Assumption}
{synopt :{opt indep:endencekp}}assumes independence between units' sizes and the probabilities, within each unit, an individual belongs to the minority or reference group{p_end}

{syntab:Inference}
{synopt :{opt repb:ootstrap(#)}}sets the number of bootstrap iterations for inference; default is {cmd:repbootstrap(200)}{p_end}
{synopt :{opt l:evel(#)}}sets level in (0,1) for additional confidence intervals; by default usual 0.9, 0.95, and 0.99 confidence intervals are saved{p_end}
{synopt :{opt noci}}restricts the command to estimation (no confidence intervals){p_end}

{syntab:Model test}
{synopt :{opt testb:inomial}}performs a test of the binomial assumption that underlies the statistical model; in {it:np} method, it is done by default with inference{p_end}

{syntab:CT method}
{synopt :{opt repct(#)}}sets the number of repetitions used to compute random allocation in method {it:ct}; default is {cmd:repct(50)}{p_end}

{syntab:Index}
{synopt :{opt atkinson(#)}}sets the parameter (often denoted "b" or "beta") for the Akinson index; default and recommended is {cmd:atkinson(0.5)}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:segregsmall} estimates classical indices (Duncan, Theil, Atkinson, Coworker, and Gini) 
to measure segregation, namely the systematic relative concentration of a minority group 
(e.g. foreigners) in some units only.
Units can be, for instance, firms, neighbordhoods, geographical areas, or classrooms.

{pstd}
Three methods can be used: {it:np} for D'Haultfoeuille and Rathelot (2017),
{it:beta} for Rathelot (2012),
and {it:ct} for Carrington and Troske (1997).
Those methods deal with the small-unit bias (overestimation of segregation
in settings where units are small, that is are composed of few individuals).
With {it:np} and {it:beta} methods, inference on the segregation indices
can be done by bootstrap.

{pstd}
Conditional segregation indices can be computed:
they measure "net" or "residual" segregation taking into 
account other characteristics (either of individuals, or of units) 
that may influence the allocation of individuals into units. 


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt format} indicates the format of the database in memory that is used for the estimation.
The argument {it:format} of the option must be either {it:unit} or {it:indiv}.
Combined with {opt conditional}, they determine the relevant variables to be passed in {varlist}.
For unconditional analysis:

{pmore}
{varlist} needs to be {it: K X} for {hi:{it:unit-level}} databases:
an observation is a unit, the variable {it:K} is the number of individuals,
and the variable {it:X} is the number of minority individuals in the unit.
The values of {it:K} (respectively {it:X}) need to be positive (resp. positive or null)
integers; for each unit, {it:X} must be lower or equal to {it:K}.

{pmore}
{varlist} needs to be {it:id_unit I_minority} for {hi:{it:individual-level}} databases:
an observation is an individual, the variable {it:id_unit} is the identifier of the unit,
the variable {I_minority} has to be a dummy variable equal to 1 if the individual belongs
to the minority group, 0 otherwise.

{phang}
{opt method} specifies the method used to compute the segregation indices.
The argument {it:method} of the option must be one of: {it:np}, {it:beta}, {it:ct}.
{it:np} and {it:beta} methods, for a given index, say the Duncan,
estimate the same parameter of interest. 
The latter, a segregation index, is defined as a function of the distribution
of a random variable p.
Each unit has a realization of that random variable, which is the probability that 
an arbitrary individual within the unit belongs to the minority group.
The main assumptions regarding data-generating-process that enable identification,
estimation, and inference are the following:
(i) the units are i.i.d.,
(ii) for each unit, we assume that conditional on the realization of p and 
the unit's size K, the number of minority individuals X is distributed according
to a Binomial(K, p).
The parameter of interest estimated by {it:ct} method is different.
In this case, the index is defined as an inequality index on the list, across units,
of the empirical proportion X/K of minority individuals.

{pmore}
{hi:{it:np}} stands for non-parametric and implements the method of
D'Haultfoeuille, X. and Rathelot, R., 2017.
Measuring segregation on small units: A partial identification analysis.
Quantitative Economics, 8(1), pp.39-73.
In short, the distribution of p is not restricted with that method.
The {it:np} method allows to estimate and make inference on four classical segregation indices:
Duncan (also known as Dissimilarity), Theil, Atkinson, and Coworker (also known as Isolation).

{pmore}
{hi:{it:beta}} implements the method of 
Rathelot, R., 2012. 
Measuring segregation when units are small: a parametric approach. 
Journal of Business & Economic Statistics, 30(4), pp.546-553.
The method assumes a Beta distributions for the random variable p,
hence its name.
The {it:beta} method allows to estimate and make inference on the four previous 
segregation indices plus the Gini index.

{pmore}
{hi:{it:ct}} implements the method of
Carrington, W.J. and Troske, K.R., 1997.
On measuring segregation in samples with small units.
Journal of Business & Economic Statistics, 15(4), pp.402-409.
The basic idea of the method is to correct the direct index obtained
as an inequality index on the proportion X/K across units.
The more variation in the proportion X/K across units, the higher such index.
Yet, due to small-unit considerations, the lower the unit size K, the more variation.
Hence a second step which consists in a correction of the index obtained.
The correction is made by comparing the actual value of the inequality index in data
with the expected value of the index when individuals are allocated randomly into units,
keeping fixed the structure, namely the units and their sizes.

{dlgtab:Conditional}

{phang}
{opt conditional} estimates conditional segragation indices.
In a nutshell, conditional segregation indices aim at measuring the relative 
concentration of minority individuals in some units only taking into account
characteristics that may influence the allocation of individuals into units.
The analogy can be made with linear models in which covariates are introduced to
account for the variance of the outcome variable.
The idea of conditional segregation analysis is to compute indices that reflect
the residual or net level of segregation while the contribution of covariates to
segregation is removed.
Covariates can be defined at the unit level or at the finer level of a position
or an individual.
Consider for instance firm segregation with firms as units and foreigners
as minority individuals. 
A unit-level characteristic can be the sector of the firm,
whereas an individual-level characteristic can be skilled / unskilled position
or worker.
The argument {it:conditional} of the option must be either {it:unit} or {it:indiv}.
The choice depends on the level at which are defined the characteristics 
considered in the conditional analysis.

{pmore}
With {hi:{it:unit-level}} covariates, {varlist} needs to be: {it: K X Z} for unit-level
databases, {it: id_unit I_minority Z} for individual-level dabatases.
The additional {it:Z} is a categorical variable that indicates the type of each unit.
{it:Z} needs to take values in 2,3,4,....
It can be constructed from a natural discrete variable, quantiles of a continuous
characteristics, or interesections of various features of units.

{pmore}
With {hi:{it:individual-level}} (or position-level) covariates, individual-level
databases are necessary and {varlist} needs to be {it: id_unit I_minority W}.
The additional {it:W} is a categorical variable that indicates the type of each 
individual. {it:W} needs to take values in 2,3,4,....
Again, it can be defined from one discrete variable, quantiles of a continuous
attribute or interesections of various characteristics defined at the individual
or position-level.

{pmore}
In both cases, the aggregated conditional index is defined and estimated as 
a convex combination of indices restricted to each type.
For unit-level covariates, we consider the different subsamples of units defined
by type Z = z.
For individual-level covariates, we consider the different subsamples of individuals
defined by type W = w.
Note that, in the latter case, the small-unit bias is compounded.
Indeed, units are divided into cells that are defined as an intersection of a unit
and an individual type.

{dlgtab:Population}

{phang}
{opt withsingle} By default, single units i.e. with only one individual,
are excluded from the analysis.
Indeed, they are uninformative about segregation since the empirical proportion
of minority individuals for such units is either 0 or 1, that is the range of 
possible values for the empirical proportion and for the underlying probability p.
The option forces to include the units with only one individual in the analysis.
The option is available both for unconditional or conditional analyses, 
either at individual or unit level.

{phang}
{opt excludingsinglepertype} The option is only available for conditional analyses
with individual-level covariates.
In this case, an index is computed for each type W = w on the subsample of
individuals with type w. Within this subsample, a "unit"
(as in the unconditional case) corresponds to a cell defined as the intersection
of a unit and an individual type.
Even if single units are excluded at the start, some of the cells can be composed
of only one individual (it will be all the more frequent as the number of individual
types is large and the units' sizes are small), and thus would be uninformative
as regards segregation.
They are included by default in the analysis but the option enables to exclude them.

{dlgtab:Assumption}

{phang}
{opt independencekp} assumes independence between the unit size {it:K} and
the probability {it:p} that an arbitrary individual in a unit belongs to the
minority.
In method {it:ct}, following the original method, the issue of possible
dependencies between {it:K} and {it:p} is not addressed since the object {it:p}
is not considered: indices are based on the proportion X/K.
Therefore, the option is irrelevant and not available with method {it:ct}.
In methods {it:np} and {it:beta}, by default, the links between those two random
variables are left unrestricted.
Hence, each method first computes one index by unit size and then 
aggregate them into the final segregation index.
The option forces independence between {it:K} and {it:p} and enables to 
perform estimation in a single step, gathering units of all sizes.
The independance between {it:K} and {it:p} is an assumption and should be used
only when reliable.

{dlgtab:Inference}

{phang}
{opt repbootstrap} specifies the number fo bootstrap repetitions used to compute
asymptotic confidence intervals for methods {it:np} and {it:beta}.
Following the original paper, confidence intervals are not available for method
{it:ct}. The default number is 200.
The confidence intervals reported are for the different segregation indices.

{phang}
{opt level} sets a personalized level, in (0,1) for confidence intervals with {it:np} 
and {it:beta} methods; sets personalized percentile of the segregation indices computed 
under random allocation with {it:ct} method.
{it:np} and {it:beta} methods: 
By default (cf. Stored results), confidence intervals (CI) at the usual 90%, 95%, and
99% level are computed and saved, and the 95% level CI is displayed in the log.
With the option, the later three CIs are still computed and stored but an
additional CI at the specified level is computed, saved, and is the one 
displayed in the log.
Note that computing and saving several confidence intervals is negligible 
as regards computation time compared to the bootstrap procedure.
{it:ct} method: By default (cf. Stored results), the 1%, 5%, 10%, 90%, 95%, and 99%
empirical percentiles of the recpt(#) draws of indices under random allocation
are reported. The higher percentiles enable to test the null hypothesis
that the data are consistent with random allocation (by comparing the 
proportion-based index obtained in the data and the relevant quantiles).
With the option, the command also saves the 1-{opt level(#)}-th and {opt level(#)}-th 
empirical quantiles to perform the randomization test at a specified level.

{phang}
{opt noci} restricts the command to estimation: only estimation is performed
and no confidence intervals are reported.
The option is motivated by the computational cost of inference by bootstrap,
which can be substantial when the units' size are large
(above approximately 30, 40 individuals)
especially with {it:np} method.
Besides, due to the construction of conditional aggregated indices, the computational cost
is globally multiplied by the number of types included in the conditional analysis.

{dlgtab:Model test}

{phang}
{opt testbinomial} enables to test the main hypothesis of the model for the 
theoretical guarantees of {it:np} and {it:beta} methods, namely that,
for each unit, conditional on the realization of p and 
the unit's size K, the number of minority individuals X is distributed according
to a Binomial(K, p).
The test relies on the techniques of method {it:np}
when the links between {it:K} and {it:p} are left unrestricted.
It uses the same boostrap procedure.
The option is available only for the unconditional case.
In conditional analysis, the test can be performed separately for each type
using restrictions through the {ifin} option.

{pmore}
Consequently, with {it:np} method and inference without 
option {opt independencekp}, the test of the 
binomial assumption is made by default as the additional computation cost is
negligible. It is stored (cf. Stored results) but not shown in the log 
by default. In this case, the option only displays the result of the test.

{pmore}
On the contrary, for any other configurations
(with other methods, without inference, with assumption of independence
between {it:K} and {it:p}), 
in addition to estimation (and possible inference), the option performs the
test relying on the {it:np} methodology.

{pmore}
In both cases, the number of bootstrap repetitions used for the test
is the same as the one specified by option {opt repbootstrap}.
Although left possible, it is not recommanded to use the option
outside method {it:np} and with inference.
Indeed, at the same cost, in addition to the test, the user could get 
the estimated segregation indices by {it:np}.

{dlgtab:CT method}

{phang}
{opt repct} sets the number of repetitions (draws) used in the {it:ct}
method to compute the expected value of the index under random
allocation of individuals in units.
The default value is 50.
It might seem low but it happens that,
empirically in many settings,
a larger number of repetitions leads to very similar estimates but
at an higher computational cost.

{dlgtab:Index}

{phang}
{opt atkinson} sets the parameter for the Atkinson index.
Compared to the other four classical indices considered
(Duncan, Theil, Coworker, Gini), the Atkinson index is parametrized
by a value in (0,1), often denoted "b" or "beta" in the segregation or
income inequality literatures in which that index is used.
The default value is 0.5.
Although at own risk, the option allows to specify another parameter.
It is not recommanded since the default value is the only one
that ensures a desirable property of the Atkinson segregation index, 
namely symmetry.
That property says that relabelling the two groups does not change the
measure of segregation.


{marker remarks}{...}
{title:Remarks}

{pstd}
Previous descriptions and options provide the basic information required
to use {cmd:segregsmall} command.
However, we strongly encourage users to read the three papers (cf. References)
that define the methods
implemented by the command for in-depth understanding.


{marker examples}{...}
{title:Examples}

{pstd}Setup TODO INCLUDE A COMMAND TO LOAD A TEST DATABASE AT UNIT-LEVEL
INCLUDED AS ANCILLARY FILES IN THE PACKAGE{p_end}
{phang2}{cmd:. TO DO}{p_end}

{pstd}Compute unconditional segregation indices using {it:np} method,
inference is done by boostrap using default 200 bootstrap iterations,
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X, format(unit) method(np)}{p_end}

{pstd}Compute unconditional segregation indices using {it:beta} method,
inference is done by boostrap using 150 bootstrap iterations,
independence is assumed between {it:K} and {it:p},
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X, format(unit) method(beta) independencekp repbootstrap(150)}{p_end}

{pstd}Compute unconditional segregation indices using {it:ct} method,
the CT correction is done using default 50 repetitions,
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X, format(unit) method(ct)}{p_end}

{pstd}Compute unconditional segregation indices using {it:ct} method,
the CT correction is done using 100 repetitions,
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X, format(unit) method(ct) repct(100)}{p_end}

{pstd}Compute unconditional segregation indices using {it:beta} method,
without inference (only estimation),
including single units (with one individual) in the analysis,
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X, format(unit) method(beta) withsingle noci}{p_end}

{pstd}Compute unconditional segregation indices using {it:np} method,
and display in the output the result of the test of the binomial assumption,
inference is done by boostrap using default 200 bootstrap iterations
that are also used to do the test,
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X, format(unit) method(np) testbinomial}{p_end}

{pstd}Compute unconditional segregation indices using {it:beta} method,
compute, save and display in the output confidence interval at 99%,
inference is done by boostrap using 500 bootstrap iterations,
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X, format(unit) method(beta) repbootstrap(500) level(0.99)}{p_end}

{pstd}Compute conditional segregation indices using {it:np} method
with unit-level covariates,
inference is done by boostrap using default 200 bootstrap iterations,
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X Z, format(unit) method(np) conditional(unit)}{p_end}

{pstd}Compute conditional segregation indices using {it:beta} method
with unit-level covariates,
without inference (only estimation),
including single units (with one individual) in the analysis,
independence is assumed between {it:K} and {it:p},
database used is at unit-level format{p_end}
{phang2}{cmd:. segregsmall K X Z, format(unit) method(beta) conditional(unit) withsingle noci independencekp}{p_end}

{pstd}Setup TODO INCLUDE A COMMAND TO LOAD A TEST DATABASE AT INDIVIDUAL-LEVEL
INCLUDED AS ANCILLARY FILES IN THE PACKAGE{p_end}
{phang2}{cmd:. TO DO}{p_end}

{pstd}Compute unconditional segregation indices using {it:np} method
inference is done by boostrap using default 200 bootstrap iterations,
database used is at individual-level format{p_end}
{phang2}{cmd:. segregsmall id_unit I_minority, format(indiv) method(np)}{p_end}

{pstd}Compute conditional segregation indices using {it:beta} method
with individual-level covariates,
inference is done by boostrap using default 200 bootstrap iterations,
including single units (with one individual) in the analysis,
database used is at individual-level format{p_end}
{phang2}{cmd:. segregsmall id_unit I_minority W, format(indiv) method(beta) conditional(indiv) withsingle}{p_end}

{pstd}Compute conditional segregation indices using {it:np} method
with individual-level covariates,
without inference (only estimation),
excluding cells (unit x type) with only one individual from the analysis,
independence is assumed between {it:K} and {it:p},
database used is at individual-level format{p_end}
{phang2}{cmd:. segregsmall id_unit I_minority W, format(indiv) method(np) conditional(indiv) excludingsinglepertype independencekp}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
The object stored in {cmd:e()} by {cmd:segregsmall} depends on the options,
notably {opt conditional} but not only.
The following presentation gathers the stored objects
according to the type of information they convey
(relative to the data used in the analysis,
relative to the method used,
relative to the results of the estimation and inference made)
and according to the specifided options. 

{p 2 4 2}{hi:Information relative to the data used:}{p_end}

{p 2 4 2}{it:Unconditional analysis}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Scalars}{p_end}
{synopt:{cmd:e(I_withsingle)}}dummy for the use of option {opt withsingle}{p_end}
{synopt:{cmd:e(nb_units_total)}}number of units in the database used, including single units{p_end}
{synopt:{cmd:e(nb_units_single)}}number of units with only one individual in the database used{p_end}
{synopt:{cmd:e(nb_units_studied)}}number of units included in the analysis, depending on option {opt withsingle} (It corresponds to the number of observations){p_end}
{synopt:{cmd:e(nb_individuals)}}number of individuals included in the analysis (i.e. within the units included in the analysis){p_end}
{synopt:{cmd:e(nb_minority_individuals)}}number of minority (or reference) individuals included in the analysis{p_end}
{synopt:{cmd:e(prop_minority_hat)}}estimated proportion of minority (or reference) group (equal to nb_minority_individuals / nb_individuals){p_end}	
{synopt:{cmd:e(K_max)}}maximal size i.e. number of individuals {it:K} in a unit across the units included in the analysis{p_end}
{synopt:{cmd:e(nb_K_with_obs)}}number of distinct unique unit sizes present across the units included in the analysis{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(list_K_with_obs)}}list of the unique unit sizes present across the units included in the analysis{p_end}

{p 2 4 2}{it:Conditional analysis}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Scalars}{p_end}
{synopt:{cmd:e(nb_types)}}number of distinct types considered in the conditional analysis{p_end}
{synopt:{cmd:e(I_unit_level_characteristic)}}dummy equal to 1 if the characteristics considered in the conditional analysis are defined at unit-level{p_end}
{synopt:{cmd:e(I_withsingle)}}dummy for the use of option {opt withsingle}{p_end}
{synopt:{cmd:e(I_excludingsinglepertype)}}dummy for the use of option {opt excludingsinglepertype}{p_end}
{synopt:{cmd:e(nb_units_total)}}number of units in the database used, including single units, across all types{p_end}
{synopt:{cmd:e(nb_units_single)}}number of units with only one individual in the database used, across all types{p_end}
{synopt:{cmd:e(nb_units_studied)}}number of units included in the analysis, depending on option {opt withsingle}, across all types{p_end}
{synopt:{cmd:e(nb_individuals)}}number of individuals included in the analysis (i.e. within the units included in the analysis), across all types{p_end}
{synopt:{cmd:e(nb_minority_individuals)}}number of minority (or reference) individuals included in the analysis, across all types{p_end}
{synopt:{cmd:e(prop_minority_hat)}}estimated proportion of minority (or reference) group, across all types (equal to nb_minority_individuals / nb_individuals){p_end}	
{synopt:{cmd:e(K_max)}}maximal size i.e. number of individuals {it:K} in a unit across the units included in the analysis, across all types{p_end}
{synopt:{cmd:e(nb_K_with_obs)}}number of distinct unique unit sizes present across the units included in the analysis, across all types{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(list_K_with_obs)}}list of the unique unit sizes present across the units included in the analysis, across all types{p_end}
{synopt:{cmd:e(type_frequencies)}}number of studied units of each type (unit-level), 
or number of individuals within studied units of each type (individual-level){p_end}
{synopt:{cmd:e(type_probabilities)}}type_frequencies expressed in proportions{p_end}
{synopt:{cmd:e(summary_info_data_per_type)}}summary information about data 
for each type.
Row {it:i}-th correspond to type {it:i} i.e. modality equal to {it:i} for the type
variable (either {it:Z} or [it:W} depending on the level of the characteristics).
For each type, the columns show different information as regards data used per type.
With unit-level covariates - The columns have the following meanings:
1: I_withsingle;
2: nb_units_total;
3: nb_units_single;
4: nb_units_studied;
5: nb_individuals;
6: nb_minority_individuals;
7: prop_minority_hat;
8: K_max;
9: nb_K_with_obs.
With individual-level covariates - In this case, note that, for each type,
the equivalent of a "unit" in the unconditional analysis is now a "cell"
defined as the intersection of a unit and a given individual type.
For instance, if we have two individuals types (W = 1 or W = 2),
a given unit with, say, 5 individuals of type 1 and 3 individuals of type 2,
is divided into two cells, one containing the 5 individuals of type 1,
the other containing the 3 individuals of type 2.
The matrix gives some information about the number of cells per type.
As for units, a cell is called single whenever it contains only one individual.
The columns have the following meanings:
1: dummy equal to 1 if single cells are included in the analysis (default),
0 otherwise (option {opt excludingsinglepertype};
2: total number of cells in the data used
(equivalent of nb_units_total where units are cells);
3: total number of single cells in the data used 
(equivalent of nb_units_single where units are cells);
4: total number of cells included in the analysis
(equivalent of nb_units_studied where units are cells);
5: nb_individuals (within the included cells);
6: nb_minority_individuals (within the included cells);
7: prop_minority_hat;
8: maximal size i.e. number of individuals {it:K} in a cell across the cells included in the analysis
(equivalent of K_max where units are cells);
9: number of distinct unique cell sizes present across the cells included in the analysis
(equivalent of nb_K_with_obs where units are cells).{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(nb_units_studied_per_type)}}with unit-level covariates only,
column 4th of summary_info_data_per_type: number of units included in the analysis per type{p_end}
{synopt:{cmd:e(nb_cells_studied_per_type)}}with individual-level covariates only,
column 4th of summary_info_data_per_type: number of cells included in the analysis per type{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Scalars}{p_end}
{synopt:{cmd:e(nb_cells_studied_sum_across_type)}}with unit-level covariates only, total number of cells
included in the analysis (sum across types){p_end}
{synopt:{cmd:e(nb_single_cells_sum_across_type)}}with unit-level covariates only, total number of single cells
(sum across types){p_end}

{p 2 4 2}{hi:Information relative to the method used:}{p_end}

{p 2 4 2}{it:Both for unconditional and conditional analysis}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Scalars}{p_end}
{synopt:{cmd:e(I_method_np)}}dummy for the use of method {it:np}{p_end}
{synopt:{cmd:e(I_method_beta)}}dummy for the use of method {it:beta}{p_end}
{synopt:{cmd:e(I_method_ct)}}dummy for the use of method {it:ct}{p_end}
{synopt:{cmd:e(I_conditional)}}dummy for conditional analysis (i.e. use of option {opt conditional}{p_end}
{synopt:{cmd:e(I_unit_level_characteristic)}}dummy for the use of unit-level covariates in case of conditional analysis{p_end}
{synopt:{cmd:e(I_hyp_independenceKp)}}dummy for the assumption of independence between {it:K} and {it:p}{p_end}
{synopt:{cmd:e(I_noci)}}dummy for the use of option {opt noci}{p_end}
{synopt:{cmd:e(nb_bootstrap_repetition)}}number of bootstrap repetitions used for inference and the test of binomial assumption, the case being{p_end}
{synopt:{cmd:e(specified_level)}}if option {opt level} is used, it indicates the specified level{p_end}
{synopt:{cmd:e(I_testbinomial)}}dummy for the use of option {opt testbinomial}{p_end}
{synopt:{cmd:e(nb_ct_repetition)}}number of repetitions used to perform the CT correction in {it:ct} method{p_end}
{synopt:{cmd:e(b_atkinson)}}parameter (often denoted "b" or "beta") used for the Atkinson index{p_end}

{p 2 4 2}{hi:Information relative to the estimation and inference performed:}{p_end}

{p 2 4 2}{it:np method - Unconditional analysis}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Scalars}{p_end}
{synopt:{cmd:e(I_constrained_case)}}dummy equal to 1 in "constrained case",
0 for "unconstrained case".
In constrained case, the distribution of the {it:p} variable is identified
in the data;
in this case, all the segregation indices are point-identified.
Otherwise, they are only partially identified
(except for the Coworker which is point-identified as long as there is 
no single units in the analysis).
We refer to the original paper (D'Haultfoeuille and Rathelot 2017) for 
details.{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(estimates_ci)}}matrix that stores the estimates of the 
identification bounds, plus, in case of inference, confidence intervals
at the usual 90%, 95%, 99% (and also at the specified level through
option {opt level} if any).
The first two columns of the matrix indicate the parameter of interest considered.
The first one, "Index", denotes the index with 1 = Duncan,
2 = Theil,
3.b = Atkinson with parameter "b",
4 = Coworker.
In case of multiple unit sizes and without assuming independence between
{it:K} and {it:p}, the indices are obtained through aggregation of 
the indices computed for each unit size.
The aggregation can be done at the unit-level level meaning that all units
enter the aggregation with the same weight; or at the individual-level 
in which units are weighted according to their size (with larger units
accounting more in the final aggregated index).
Again, we refer to the original paper (D'Haultfoeuille and Rathelot 2017) for 
details.
The second column, "I_indiv_weight", contains dummies equal to 0 for unit-level
aggregation (unweighted) and 1 for individual-level aggregation (weighted).
As the additional computation cost is negligible compared to the estimation,
the command {cmd:segregsmall} computes both indices.
With {it:np} method, due to the partial identification analysis,
there are two ways of constructing
confidence intervals for the segregation indices.
In case of inference, before the confidence intervals, the fifth column
"I_ciboundary" is a dummy equal to 0 if the CIs are constructed using
the "interior" CIs, and equal to 1 if the CIs are constructed using
the "boundary" CIs.
We refer to the original paper (D'Haultfoeuille and Rathelot 2017) for 
details.{p_end}
{synopt:{cmd:e(info_distribution_of_p)}}stores information about the distribution
of {it:p} obtained through the estimation.
We provide a short explanation below and refer to the original paper
(D'Haultfoeuille and Rathelot 2017) for details.
Basically, in the data, for each unit size {it:K}=k,
we can identify the first k moments of the distribution of {it:p}.
Those moments provide information about the distribution of {it:p}
conditional on {it:K}=k.
In general, a distribution is uniquely defined by its entire list of moment
m1, m2, m3, ....
But it happens that when the distribution happens to be on the frontier
of the moment space that a finite number of first moments are sufficient
to uniquely defined the distribution.
Furthermore, in such cases, the distribution have to be a discrete distribution
with specific support points (a.k.a. nodes) and masses.
Without assuming independence between {it:K} and {it:p}, 
for each unit size {it:K}=k
the matrix info_distribution_of_p stores what is known about the distribution
of {it:p} conditional on {it:K}=k.
It should be read three rows by three row.
The first "K_info" introduces the information for the unit size considered,
size k, number of units with that size in the data included in the analysis,
proportion of units wit that size among all included in the analysis.
The fourth column is a dummy equal to 0 when the estimated distribution {it:p}|{it:K}=k
is not uniquely defined, that is its corresponding vector of estimated k first moments
is not in the boundary of the moment space. 
In this case, the two following rows of the matrix are empty (missing values .).
The dummy is equal to 1 when 
the estimated distribution {it:p}|{it:K}=k is uniquely defined
by the estimated k first moments because the vector happens to be in the boundary
of the moment space; it is furthermore a discrete distribution.
In this case, the fifth column "nb_nodes" indicate the number of support points
of the estimated distribution {it:p}|{it:K}=k,
and the two following rows report the support points or nodes and the
associated masses.
When assuming independence between {it:K} and {it:p}
(option {opt independencekp}),
we can identify up to the K_max (the maximal size among the units included
in the analysis) first moment of the unconditional distribution of {it:p}.
The matrix info_distribution_of_p has three rows.
The first one reminds information about the data used: K_max, number of 
units included in the analysis, dummy variable whether the distribution of 
{it:p} is constrained or not. If constrained the number of nodes
and the estimated distribution of {it:p}.
The motivation of info_distribution_of_p is that a segregation measure
is a function of the joint distribution of ({it:K},{it:p}) (or of {it:p}) only.
By providing those estimated distribution in the constrained case, it allows
users to consider other segregations indices.{p_end}
{synopt:{cmd:e(test_binomial_results)}}when the binomial test is done, which 
is equivalent to doing inference with {it:np} method (as the additional
computational cost is negligible compared to the bootstrap procedure),
the result of the test is stored in that matrix.
With {it:np} method, the option {opt testbinomial} only forces the display
of the matrix test_binomial_results in the log, but it is always stored.
It has one row and two columns, the first one stores the value of the test
statistic, the second one the p-value of the test where the null hypothesis
H0 is the binomial assumption: conditional on unit size {it:K} and 
probability than an arbitrary individual in the unit belongs to the 
minority (or reference) group, the number {it:X} of minority (or reference)
individuals in the unit is distributed according to a Binomial({it:K}, {it:p}).{p_end}

{p 2 4 2}{it:np method - Conditional analysis}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Scalars}{p_end}
{synopt:{cmd:e(I_constrained_case)}}dummy equal to 1 when
it is a "constrained case" for each type; and therefore the aggregated
conditional indices are point-identified.{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(I_constrained_case_per_type)}}equivalent to I_constrained_case
in unconditional analysis type by type.{p_end}
{synopt:{cmd:e(estimates_ci_aggregated)}}equivalent of estimates_ci
for the conditional aggregated index.{p_end}
{synopt:{cmd:e(estimates_ci_type_1)}}equivalent to matrix estimates_ci
(that stores the estimates of the identification bounds and confidence intervals
for the indices) for type 1.{p_end}
{synopt:{cmd:e(estimates_ci_type_2)}}idem for type 2.{p_end}
{synopt:{cmd:e(estimates_ci_type_#)}}and so forth for
each type up to the number of types used in the conditional analysis.{p_end}

{p 2 4 2}{it:For the other methods, the stored objects are overall similar. They are reported below stressing only the differences with {it:np} method if necessary}.{p_end}

{p 2 4 2}{it:beta method - Unconditional analysis}:
Both {it:beta} and {it:np} methods estimate exactly the same parameters of interest
but with different estimation techniques:
{it:beta} method assumes a Beta distribution for {it:p}|{it:K} while
{it:np} method leaves it unrestricted (non-parametric approach).
As regards stored objects, the main differences come from the fact that the
issue of partial identification in {it:np} method is absent with {it:beta} method
since the Beta parametrization ensures point-identification.
Therefore, with {it:beta} method, the matrices that store estimation 
and inference results directly report the estimated indices and their
confidence intervals.
Besides, the identification and estimation techniques used in {it:np} requires
conditions that are not satisfied by the Gini index, hence it is not reported
with {it:np} method. It is estimated with {it:beta} method.
{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(estimates_ci)}}analoguous to estimates_ci with {it:np} method,
except that the estimated indices replace the estimates of the identification bounds.{p_end}
{synopt:{cmd:e(info_distribution_of_p)}}{it:beta} method assumes a beta
distribution for {it:p} conditional on {it:K}.
For each unit size k (a row of the matrix), info_distribution_of_p
reports the estimated parameter of the beta distribution:
the first column "K" is the unit size;
the second column "nb_units" reports the number of units of size k
included in the analysis;
the third column "prop_unit" is the proportion of such units among
all units included in the analysis;
the fourth column "nb_comp" is for consistency with the original method
in Rathelot (2012) which uses a mixture of Beta distribution for 
the distribution of {it:p}. Here, the method uses a Beta distribution
only (which can be seen as a mixture of betas with 1 component, hence
the values of that column).
Actually, simulations reveal that the indices obtained via
a Beta distribution or a mixture of Beta with two components or more are 
extremely close in most settings.
The package uses only a Beta distribution to decrease the computational cost.
Besides, the {it:np} and {it:beta} methods estimate exactly the same parameters
of interest (contrary to {it:ct}).
Therefore, the comparison of both estimates enables to check whether the Beta
distribution assumption for {it:p} is sensible;
the fifth and sixth columns reports the estimates of the two shape parameters of
a Beta distribution (often called "alpha" and "beta").
When independence between {it:K} and {it:p} is assumed with option 
{opt independencekp}, the estimation of the parameters of the distribution
of {it:p}, assumed to be a Beta, is performed gathering all units whatever
their size. estimates_ci contains then only one row that reports the
estimates of the two parameters of the Beta distribution and informations
about the number of units used for that estimation.{p_end}
{synopt:{cmd:e(test_binomial_results)}} similar to the {it:np} method
when option {opt testbinomial} is specified.
It is not recommended however since it has the same computational cost
as the estimation and inference with {it:np} method, so it would be 
preferable to obtain the result of the test and the estimation from the
{it:np} method.{p_end}

{p 2 4 2}{it:beta method - Conditional analysis}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(estimates_ci_aggregated)}}same results as estimates_ci aggregated
across types.{p_end}
{synopt:{cmd:e(estimates_ci_type_#)}}same results as estimates_ci for type #.{p_end}

{p 2 4 2}{it:ct method - Unconditional analysis}:
Contrary to methods {it:np} and {it:beta}, {it:ct} estimate a slightly
distinct parameter of interest: the segregation indices considered
with {it:ct} method are defined as functions of the empirical proportions {it:X}/{it:K}.
Given those definitions, the different weights for aggregation (unit or individual)
used in {it:np} and {it:beta} methods are irrelevant.
The estimation is made gathering all units, independent of their sizes.
The method proposed by Carrington and Troske (1997) does not provide
confidence intervals for the segregation indices.
As a consequence, the matrices that stores estimation results differ (cf. below).{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(estimates_ci)}}stores the result of estimation using {it:ct} method.
Each row corresponds to a given index among Duncan, Theil, Atkinson, Coworker,
and Gini.
The first column "index" is simply an encoding of the index (as in the {it:np} method).
The second column "prop_based" reports the index computed
as an inequality index on the list of empirical proportions X/K across units;
it is sometimes called the "naive" index.
The third column "mean_under_ra" (ra standing for random allocation) reports
the mean over {opt recpt(#)} repetitions
of the index in the counterfactual situation of randomness where
individuals are allocated randomly across units, keeping fixed the number of 
units and their sizes.
The fourth column "CT_corr" reports the corrected index; the correction
uses the naive index and the expected value of the index under random allocation.
As a tribute to Cortese, Falk, and Cohen (1976) 
(that appears to be the first to propose the idea of a comparison to randomness
for segregation indices), the sixth column shows the original segregation score proposed by
Cortese and coauthors.
That score is defined as the standardized version of the index: that is
the proportion-based index minus its expected value under random allocation
divided by the standard deviation of the index under random allocation. 
The latter quantity is reported in the fifth columb "sd_under_ra".
The remaining colums show the empirical quantiles of the segregation indices
obtained under random allocation ({opt recpt(#)} draws of indices under random allocation).
By default, the 1%, 5%, 10%, 90%, 95%, and 99% percentiles are reported
(sixth to twelfth columns).
With option {opt level}, the 1-{opt level(#)}-th and {opt level(#)}-th quantiles are also reported
in the thirteenth and fourteenth columns respectively.

{p_end}
{synopt:{cmd:e(test_binomial_results)}} similar to the {it:np} method
when option {opt testbinomial} is specified.
It is not recommended however since it has the same computational cost
as the estimation and inference with {it:np} method, so it would be 
preferable to obtain the result of the test and the estimation from the
{it:np} method.{p_end}

{p 2 4 2}{it:ct method - Conditional analysis}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 32 2: Matrices}{p_end}
{synopt:{cmd:e(estimates_ci_aggregated)}}same results as estimates_ci aggregated
across types.{p_end}
{synopt:{cmd:e(estimates_ci_type_#)}}same results as estimates_ci for type #.{p_end}

{p 2 4 2}{hi:Miscellaneous:}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:segregsmall}{p_end}
{synopt:{cmd:e(cmd_arguments)}}arguments of the command as typed{p_end}


{marker authors}{...}
{title:Authors}

{pstd}Lucas Girard (CREST), Xavier D'Haultfoeuille (CREST), Roland Rathelot (Warwick){p_end}
{pstd}Support: lucas.girard@ensae.fr{p_end}


{marker references}{...}
{title:References}

{pstd}
Carrington, W.J. and Troske, K.R., 1997.
On measuring segregation in samples with small units.
Journal of Business & Economic Statistics, 15(4), pp.402-409.

{pstd}
Cortese, C. F., Falk, R. F. and Cohen, J. K., 1976.
Further considerations on the methodological analysis of segregation indices.
American sociological review, pp.630–637.

{pstd}
D'Haultfoeuille, X. and Rathelot, R., 2017.
Measuring segregation on small units: A partial identification analysis.
Quantitative Economics, 8(1), pp.39-73.

{pstd}
Rathelot, R., 2012. 
Measuring segregation when units are small: a parametric approach. 
Journal of Business & Economic Statistics, 30(4), pp.546-553.
