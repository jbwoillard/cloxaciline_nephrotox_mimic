# MSSA Bacteremia — ASP vs Cefazolin Nephrotoxicity (MIMIC-IV v3.1)

Causal machine-learning analysis of nephrotoxicity associated with anti-staphylococcal penicillins (ASP: nafcillin/oxacillin) versus cefazolin in methicillin-susceptible *Staphylococcus aureus* (MSSA) bacteremia, using MIMIC-IV v3.1.

---

## Scientific Question

Among adults with MSSA bacteremia, is initial treatment with anti-staphylococcal penicillins associated with a higher short-term risk of acute kidney injury (AKI) than cefazolin, and can treatment-effect heterogeneity be evaluated within a derivation–test individualized treatment-effect framework?

---

## Current Repository Scope

This repository contains only the **final analysis code** and the **final manuscript-oriented report**.

It includes:

- the R scripts required to build the cohort and analysis dataset
- the final unified derivation–test causal ML pipeline
- supplementary robustness and revision analyses
- the final report source (`.Rmd`) and rendered HTML report

It does **not** include:

- MIMIC-IV source data
- derived datasets
- intermediate results
- tables
- figures


All output directories such as `derived/`, `results/`, `tables/`, and `figures/` are generated locally at runtime and are not tracked in the repository.

---

## Final Analysis Framework

The repository reflects the **final unified derivation–test pipeline** used for the main analysis.

### Cohort construction
Scripts `01` to `05` build the analysis dataset:
1. MSSA bacteremia candidate cohort
2. Treatment assignment and time zero
3. AKI outcome derivation
4. Baseline covariate construction
5. Final analysis dataset assembly

### Primary causal ML framework
Scripts `12` to `14` implement the main modeling strategy:
- strict 70/30 derivation–test split
- 5-fold cross-validation in the derivation set only
- candidate model comparison across multiple causal ML architectures
- locked held-out test-set evaluation of the selected model
- post-selection full-cohort refit reported secondarily

### Supplementary analyses
Scripts `15` and `16` provide:
- estimator robustness analysis across alternative causal estimators
- additional manuscript revision analyses, including supplementary sensitivity analyses and diagnostic outputs

---

## Main Findings

### Selected model
Among five candidate models evaluated in the derivation cohort, the selected model was:

**Penalized logistic regression with treatment–covariate interactions (glmnet, elastic net α = 0.5)**

### Primary validated estimate (held-out test set, N = 133)
**DR ATE = 0.153 (95% CI −0.183 to 0.416)**

### Secondary post-selection full-cohort refit (N = 445)
**DR ATE = 0.160 (95% CI 0.069 to 0.267)**

Both estimates were directionally consistent with a higher AKI risk under ASP than under cefazolin.

---

## Cohort Summary

| Step | N |
|---|---:|
| *S. aureus* blood cultures (MIMIC-IV v3.1) | 2,223 |
| After exclusion of oxacillin-resistant isolates | 1,001 |
| With qualifying ASP or cefazolin treatment | 458 |
| After new-user restriction | 446 |
| Analyzable cohort (non-missing primary outcome) | 445 |

Treatment groups in the analyzable cohort:
- ASP: N = 138
- Cefazolin: N = 307

---

## Repository Structure

```text
.
├── README.md
├── R/
│   ├── 01_build_mssa_cohort.R
│   ├── 02_define_treatment_and_time_zero.R
│   ├── 03_build_outcomes_aki.R
│   ├── 04_build_covariates.R
│   ├── 05_analysis_descriptive.R
│   ├── 12_unified_derivation_test_pipeline.R
│   ├── 13_candidate_model_selection_unified.R
│   ├── 14_selected_model_primary_and_ite_evaluation.R
│   ├── 15_ate_robustness_across_estimators.R
│   └── 16_additional_analyses.R
└── report/
    ├── final_report_unified_ml_selection.Rmd
    └── final_report_unified_ml_selection.html
Script Overview
Script	Role	Main local outputs
01_build_mssa_cohort.R	Build MSSA bacteremia candidate cohort	derived/mssa_candidates.rds
02_define_treatment_and_time_zero.R	Define exposure and time zero	derived/cohort_treated.rds
03_build_outcomes_aki.R	Build AKI outcomes	derived/outcomes.rds
04_build_covariates.R	Build baseline covariates	derived/covariates.rds
05_analysis_descriptive.R	Merge cohort, outcomes, and covariates	derived/analysis_dataset.rds
12_unified_derivation_test_pipeline.R	Create derivation/test split and imputation rules	derived/derivation_set_unified.rds, derived/test_set_unified.rds, derived/imputation_rules_unified.rds
13_candidate_model_selection_unified.R	Candidate model comparison in derivation set	local model-selection outputs
14_selected_model_primary_and_ite_evaluation.R	Primary test-set evaluation and full-cohort refit	local selected-model outputs
15_ate_robustness_across_estimators.R	ATE robustness across estimators	local supplementary outputs
16_additional_analyses.R	Additional manuscript and revision analyses	local supplementary outputs
How to Run the Pipeline

Scripts are intended to be run sequentially.

source("R/01_build_mssa_cohort.R")
source("R/02_define_treatment_and_time_zero.R")
source("R/03_build_outcomes_aki.R")
source("R/04_build_covariates.R")
source("R/05_analysis_descriptive.R")

source("R/12_unified_derivation_test_pipeline.R")
source("R/13_candidate_model_selection_unified.R")
source("R/14_selected_model_primary_and_ite_evaluation.R")
source("R/15_ate_robustness_across_estimators.R")
source("R/16_additional_analyses.R")

rmarkdown::render("report/final_report_unified_ml_selection.Rmd")
Paths and Local Data Layout

The scripts are intended to use project-relative paths (for example via here::here()), rather than absolute machine-specific paths.

Expected local layout:

<project-root>/
├── mimic-iv-3.1/
│   ├── hosp/
│   └── icu/
├── derived/
├── results/
├── tables/
├── figures/
├── R/
└── report/

Users should adapt local path configuration if needed.

Data Access

This repository contains code only.

MIMIC-IV data are not included and must not be redistributed through this repository. Access to MIMIC-IV requires credentialed access through PhysioNet and compliance with the corresponding data use agreement.

Expected MIMIC-IV local structure:

mimic-iv-3.1/
├── hosp/
│   ├── microbiologyevents.csv.gz
│   ├── admissions.csv.gz
│   ├── patients.csv.gz
│   ├── emar.csv.gz
│   ├── labevents.csv.gz
│   ├── diagnoses_icd.csv.gz
│   └── prescriptions.csv.gz
└── icu/
    └── icustays.csv.gz
Final Report

The final manuscript-oriented report is available in:

report/final_report_unified_ml_selection.Rmd
report/final_report_unified_ml_selection.html

The HTML file is the rendered version corresponding to the final repository state.

Methodological References
Buell KG, Spicer AB, Casey JD, et al. Individualized treatment effects of oxygen targets in critically ill adults. JAMA. 2024;331(14):1195–1204.
Munroe ES, Spicer A, Castellvi-Font A, et al. A framework for deriving and validating individualized treatment effects in critical care. Lancet Respiratory Medicine. 2025;13:556–568.
VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. Ann Intern Med. 2017;167(4):268–274.
KDIGO Clinical Practice Guideline for Acute Kidney Injury. Kidney Int Suppl. 2012;2:1–138.
