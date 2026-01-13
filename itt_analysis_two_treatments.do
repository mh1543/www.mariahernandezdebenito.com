/*******************************************************************************
* ITT ANALYSIS FOR RCT WITH TWO TREATMENTS
* Author: Maria Hernandez de Benito
* Date: January 2026
* Purpose: Evaluate Intent-to-Treat effects in an RCT with two treatment arms
*******************************************************************************/

clear all
set more off
set matsize 10000

* Set working directory
cd "YOUR_WORKING_DIRECTORY"

* Install required packages (uncomment if needed)
* ssc install estout, replace
* ssc install ivreg2, replace
* ssc install ranktest, replace
* ssc install coefplot, replace

/*******************************************************************************
* 1. LOAD DATA
*******************************************************************************/

use "YOUR_DATA_FILE.dta", clear

* Define key variables
global outcome "outcome_var"           // Main outcome variable
global treatment1 "treat1"             // Treatment 1 indicator (0/1)
global treatment2 "treat2"             // Treatment 2 indicator (0/1)
global control "control"               // Control group indicator (0/1)

* Baseline covariates for balance checks and control
global covariates "age female education income"

* Additional controls (if needed)
global controls "$covariates"

/*******************************************************************************
* 2. DESCRIPTIVE STATISTICS
*******************************************************************************/

* Summary statistics by treatment group
eststo clear
estpost tabstat $outcome $covariates, by(treatment_group) statistics(mean sd) columns(statistics)
esttab using "tables/summary_stats.tex", replace ///
    cells("mean(fmt(3)) sd(fmt(3))") noobs nonumber ///
    title("Summary Statistics by Treatment Group")

/*******************************************************************************
* 3. BALANCE CHECKS
*******************************************************************************/

* Test randomization balance across treatment groups
eststo clear

foreach var of global covariates {
    * Treatment 1 vs Control
    reg `var' $treatment1, robust
    eststo t1_`var'

    * Treatment 2 vs Control
    reg `var' $treatment2, robust
    eststo t2_`var'

    * Joint test: Treatment 1 vs Treatment 2 vs Control
    reg `var' $treatment1 $treatment2, robust
    test $treatment1 = $treatment2 = 0
    eststo joint_`var'
}

* Export balance table
esttab t1_* using "tables/balance_treatment1.tex", replace ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    title("Balance Test: Treatment 1 vs Control")

esttab t2_* using "tables/balance_treatment2.tex", replace ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    title("Balance Test: Treatment 2 vs Control")

/*******************************************************************************
* 4. MAIN ITT ESTIMATES
*******************************************************************************/

eststo clear

* Model 1: Simple comparison (no controls)
eststo itt1: reg $outcome $treatment1 $treatment2, robust
estadd local controls "No"
estadd local fe "No"

* Model 2: With baseline controls
eststo itt2: reg $outcome $treatment1 $treatment2 $controls, robust
estadd local controls "Yes"
estadd local fe "No"

* Model 3: With fixed effects (if applicable, e.g., strata or cluster FE)
* eststo itt3: areg $outcome $treatment1 $treatment2 $controls, absorb(strata_id) robust
* estadd local controls "Yes"
* estadd local fe "Strata"

* Test equality of treatment effects
test $treatment1 = $treatment2

* Export ITT results
esttab itt* using "tables/itt_main.tex", replace ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    keep($treatment1 $treatment2) ///
    scalars("controls Controls" "fe Fixed Effects" "N Observations" "r2 R-squared") ///
    title("Intent-to-Treat Effects") ///
    mtitles("No Controls" "With Controls" "With FE") ///
    addnote("Robust standard errors in parentheses." ///
            "* p<0.10, ** p<0.05, *** p<0.01")

* Display results
estout itt*, cells(b(star fmt(3)) se(par fmt(3))) ///
    stats(controls fe N r2, fmt(0 0 0 3)) ///
    legend label varlabels(_cons Constant)

/*******************************************************************************
* 5. HETEROGENEITY ANALYSIS
*******************************************************************************/

* Define subgroups for heterogeneity analysis
global het_vars "female high_education"  // Modify as needed

eststo clear
local counter = 1

foreach hetvar of global het_vars {

    * Create interaction terms
    gen treat1_x_`hetvar' = $treatment1 * `hetvar'
    gen treat2_x_`hetvar' = $treatment2 * `hetvar'

    * Regression with interactions
    eststo het`counter': reg $outcome $treatment1 $treatment2 ///
        `hetvar' treat1_x_`hetvar' treat2_x_`hetvar' $controls, robust

    * Test heterogeneity
    test treat1_x_`hetvar' = treat2_x_`hetvar' = 0

    estadd local hetvar "`hetvar'"
    local counter = `counter' + 1
}

* Export heterogeneity results
esttab het* using "tables/itt_heterogeneity.tex", replace ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    keep($treatment1 $treatment2 *_x_*) ///
    scalars("hetvar Heterogeneity Variable" "N Observations" "r2 R-squared") ///
    title("ITT Heterogeneous Effects")

/*******************************************************************************
* 6. SUBGROUP ANALYSIS
*******************************************************************************/

* Analyze by subgroups separately
eststo clear

* By gender (example)
eststo male: reg $outcome $treatment1 $treatment2 $controls if female == 0, robust
estadd local sample "Male"

eststo female: reg $outcome $treatment1 $treatment2 $controls if female == 1, robust
estadd local sample "Female"

* Export subgroup results
esttab male female using "tables/itt_subgroups.tex", replace ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    keep($treatment1 $treatment2) ///
    scalars("sample Sample" "N Observations" "r2 R-squared") ///
    title("ITT Effects by Subgroup")

/*******************************************************************************
* 7. ROBUSTNESS CHECKS
*******************************************************************************/

eststo clear

* Robustness 1: Alternative outcome specification (if applicable)
* eststo robust1: reg log_outcome $treatment1 $treatment2 $controls, robust
* estadd local spec "Log outcome"

* Robustness 2: Clustered standard errors (if applicable)
* eststo robust2: reg $outcome $treatment1 $treatment2 $controls, cluster(cluster_id)
* estadd local spec "Clustered SE"

* Robustness 3: Trimming outliers
summarize $outcome, detail
local p99 = r(p99)
local p1 = r(p1)
eststo robust3: reg $outcome $treatment1 $treatment2 $controls ///
    if $outcome >= `p1' & $outcome <= `p99', robust
estadd local spec "Trimmed outliers"

* Robustness 4: Alternative control specifications
eststo robust4: reg $outcome $treatment1 $treatment2, robust
estadd local spec "No controls"

* Export robustness results
esttab robust* using "tables/itt_robustness.tex", replace ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    keep($treatment1 $treatment2) ///
    scalars("spec Specification" "N Observations" "r2 R-squared") ///
    title("Robustness Checks")

/*******************************************************************************
* 8. MULTIPLE HYPOTHESIS TESTING CORRECTIONS
*******************************************************************************/

* Bonferroni correction for multiple outcomes (if applicable)
* Number of outcomes tested
local num_outcomes = 2  // Two treatments

* Calculate adjusted p-values
matrix pvals = J(`num_outcomes', 2, .)
matrix colnames pvals = "Original_p" "Bonferroni_p"
matrix rownames pvals = "Treatment1" "Treatment2"

* Store original p-values from main specification
quietly reg $outcome $treatment1 $treatment2 $controls, robust
test $treatment1 = 0
matrix pvals[1,1] = r(p)
matrix pvals[1,2] = min(r(p) * `num_outcomes', 1)

test $treatment2 = 0
matrix pvals[2,1] = r(p)
matrix pvals[2,2] = min(r(p) * `num_outcomes', 1)

* Display adjusted p-values
matrix list pvals

/*******************************************************************************
* 9. TREATMENT EFFECT VISUALIZATION
*******************************************************************************/

* Coefficient plot comparing treatment effects
eststo clear
eststo: quietly reg $outcome $treatment1 $treatment2 $controls, robust

coefplot, keep($treatment1 $treatment2) ///
    xline(0, lcolor(red) lpattern(dash)) ///
    title("ITT Treatment Effects with 95% CI") ///
    xtitle("Effect Size") ///
    coeflabels($treatment1 = "Treatment 1" $treatment2 = "Treatment 2") ///
    graphregion(color(white))
graph export "figures/itt_coefplot.png", replace

/*******************************************************************************
* 10. POWER CALCULATIONS (POST-HOC)
*******************************************************************************/

* Calculate minimum detectable effect (MDE) given sample size
* This requires actual sample sizes and variance
quietly summarize $outcome if $control == 1
local sd_control = r(sd)
local n_control = r(N)

quietly summarize $outcome if $treatment1 == 1
local n_treat1 = r(N)

quietly summarize $outcome if $treatment2 == 1
local n_treat2 = r(N)

* Display sample information
display "Sample sizes:"
display "Control: `n_control'"
display "Treatment 1: `n_treat1'"
display "Treatment 2: `n_treat2'"
display "Control group SD: `sd_control'"

/*******************************************************************************
* 11. ATTRITION ANALYSIS (IF APPLICABLE)
*******************************************************************************/

* Check if attrition is balanced across treatment arms
* gen attrition = (missing($outcome))
*
* reg attrition $treatment1 $treatment2 $controls, robust
* eststo attrition_check
*
* esttab attrition_check using "tables/attrition_analysis.tex", replace ///
*     se star(* 0.10 ** 0.05 *** 0.01) ///
*     title("Attrition Analysis")

/*******************************************************************************
* 12. EXPORT SUMMARY TABLE
*******************************************************************************/

* Create final summary table with key results
eststo clear

eststo col1: quietly reg $outcome $treatment1 $treatment2, robust
estadd local controls "No"
quietly test $treatment1 = $treatment2
estadd scalar p_equal = r(p)

eststo col2: quietly reg $outcome $treatment1 $treatment2 $controls, robust
estadd local controls "Yes"
quietly test $treatment1 = $treatment2
estadd scalar p_equal = r(p)

* Calculate effect sizes (Cohen's d)
quietly summarize $outcome if $control == 1
local sd_pooled = r(sd)

quietly reg $outcome $treatment1 $treatment2 $controls, robust
local coef_t1 = _b[$treatment1]
local coef_t2 = _b[$treatment2]
local cohend_t1 = `coef_t1' / `sd_pooled'
local cohend_t2 = `coef_t2' / `sd_pooled'

estadd scalar cohend_t1 = `cohend_t1': col2
estadd scalar cohend_t2 = `cohend_t2': col2

* Final table
esttab col* using "tables/itt_final.tex", replace ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    keep($treatment1 $treatment2) ///
    scalars("controls Controls" "p_equal P-value: T1=T2" ///
            "cohend_t1 Cohen's d (T1)" "cohend_t2 Cohen's d (T2)" ///
            "N Observations" "r2 R-squared") ///
    title("Intent-to-Treat Effects: Final Results") ///
    addnote("Robust standard errors in parentheses." ///
            "* p<0.10, ** p<0.05, *** p<0.01" ///
            "Cohen's d calculated using control group SD")

/*******************************************************************************
* END OF DO FILE
*******************************************************************************/

log close
