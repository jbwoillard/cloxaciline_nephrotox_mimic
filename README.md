# MSSA Bacteremia — ASP vs Cefazolin Nephrotoxicity (MIMIC-IV v3.1)

Causal machine-learning analysis of nephrotoxicity associated with anti-staphylococcal penicillins (ASP: nafcillin/oxacillin) versus cefazolin in MSSA bacteremia, using MIMIC-IV v3.1.

---

## Scientific Question

Among adults with MSSA bacteremia, is treatment with anti-staphylococcal penicillins associated with a higher short-term risk of acute kidney injury (AKI) than cefazolin, and can treatment-effect heterogeneity be identified and validated using a derivation–test individualized treatment-effect framework?

---

## Analysis Framework

The repository implements two sequential layers:

**Layer 1 — Data preparation and initial causal analysis (scripts 01–10):** Cohort construction, outcome derivation, covariate assembly, a first causal forest estimate (grf), 8 sensitivity analyses, and reviewer-requested supplementary analyses (missing-data audit, MIA robustness, E-value, overlap diagnostics).

**Layer 2 — Unified derivation/test causal ML pipeline (scripts 12–14):** A strict 70/30 derivation–test split, 5-fold CV model selection across five candidate causal ML models in the derivation set only, followed by a single locked evaluation of the selected model in the held-out test set. This layer supersedes the primary ATE estimate from Layer 1 and is the main framework for both average and individualized treatment effect reporting.

---

## Main Findings

### Selected model
Five candidate models were compared in the derivation set (5-fold CV, N=312). The best-performing model across all three selection metrics was:

**Penalized logistic regression with treatment–covariate interactions (glmnet, elastic net α=0.5)**

### Primary validated estimate (held-out test set, N=133)
**DR ATE = 0.153 (95% CI −0.183 to 0.416)**

### Secondary post-selection full-cohort refit (N=445)
**DR ATE = 0.160 (95% CI 0.069 to 0.267)**

Both estimates are directionally consistent with higher AKI risk under ASP than cefazolin. Predicted individualized treatment effects (ITEs) were all positive in the test set, suggesting variation in the magnitude rather than direction of harm.

---

## Cohort

| Step | N |
|---|---:|
| *S. aureus* blood cultures (MIMIC-IV v3.1) | 2,223 |
| After exclusion of oxacillin-resistant isolates (MRSA) | 1,001 |
| With qualifying treatment (ASP or CEF in window) | 458 |
| After new-user design (exclude prior exposure) | 446 |
| Analyzable cohort (non-missing primary outcome) | 445 |

Treatment groups (analyzable cohort): ASP N=138, Cefazolin N=307.

---

## Repository Structure

```
.
├── R/
│   ├── 01_build_mssa_cohort.R
│   ├── 02_define_treatment_and_time_zero.R
│   ├── 03_build_outcomes_aki.R
│   ├── 04_build_covariates.R
│   ├── 05_analysis_descriptive.R
│   ├── 06_primary_causal_grf.R
│   ├── 07_sensitivity_analyses.R
│   ├── 08_generate_tables_figures.R
│   ├── 09_render_report.R
│   ├── 10_reviewer_requested_supplementary_analyses.R
│   ├── 12_unified_derivation_test_pipeline.R
│   ├── 13_candidate_model_selection_unified.R
│   └── 14_selected_model_primary_and_ite_evaluation.R
├── report/
│   └── final_report_unified_ml_selection.Rmd
├── install_r_packages.R
└── README.md
```

Directories `derived/`, `results/`, `tables/`, and `figures/` are generated at runtime and are not tracked in the repository. MIMIC-IV source data is not included (see [Data access](#data-access)).

---

## Scripts

### Layer 1 — Data preparation and initial analysis

| Script | Input | Output |
|---|---|---|
| `01_build_mssa_cohort.R` | `microbiologyevents`, `admissions`, `patients` (MIMIC-IV hosp) | `derived/mssa_candidates.rds` |
| `02_define_treatment_and_time_zero.R` | `mssa_candidates.rds`, `emar` | `derived/cohort_treated.rds` |
| `03_build_outcomes_aki.R` | `cohort_treated.rds`, `labevents` | `derived/outcomes.rds` |
| `04_build_covariates.R` | `cohort_treated.rds`, `labevents`, `diagnoses_icd`, `icustays`, `prescriptions` | `derived/covariates.rds`, `tables/covariate_dictionary.csv` |
| `05_analysis_descriptive.R` | `cohort_treated.rds`, `outcomes.rds`, `covariates.rds` | `derived/analysis_dataset.rds`, `tables/table1.csv`, `tables/raw_outcomes.csv`, `tables/smd_before_weighting.csv`, `figures/raw_aki_rates.png` |
| `06_primary_causal_grf.R` | `analysis_dataset.rds` | `results/model_objects.rds`, `results/ate_results.csv`, `results/cate_predictions.csv`, `results/variable_importance.csv`, `results/best_linear_projection.csv`, `derived/analysis_dataset_with_cate.rds` |
| `07_sensitivity_analyses.R` | `analysis_dataset.rds`, `model_objects.rds`, `emar` | `results/sensitivity_results.csv` |
| `08_generate_tables_figures.R` | `analysis_dataset.rds`, `analysis_dataset_with_cate.rds`, `model_objects.rds`, `sensitivity_results.csv`, other tables | `figures/` (10 PNG), `tables/ate_summary.csv` |
| `09_render_report.R` | all derived/results/tables | `report/final_report.md` |
| `10_reviewer_requested_supplementary_analyses.R` | `analysis_dataset.rds`, `model_objects.rds`, `report/final_report_unified_ml_selection.Rmd` | `tables/reviewer_missingness_summary.csv`, `tables/supplement_mia_vs_primary.csv`, `tables/supplement_evalue.csv`, `tables/supplement_overlap_summary.csv`, `results/mia_model_objects.rds`, `figures/mia_cate_distribution.png`, `figures/supplement_propensity_overlap_reviewer.png` |

### Layer 2 — Unified derivation/test pipeline

| Script | Input | Output |
|---|---|---|
| `12_unified_derivation_test_pipeline.R` | `analysis_dataset.rds` | `derived/derivation_set_unified.rds`, `derived/test_set_unified.rds`, `derived/imputation_rules_unified.rds`, `tables/derivation_test_split_summary.csv` |
| `13_candidate_model_selection_unified.R` | `derivation_set_unified.rds` | `tables/candidate_model_performance_derivation.csv`, `results/derivation_model_selection_metrics.rds` |
| `14_selected_model_primary_and_ite_evaluation.R` | `derivation_set_unified.rds`, `test_set_unified.rds`, `imputation_rules_unified.rds`, `derivation_model_selection_metrics.rds`, `analysis_dataset.rds` | `results/selected_unified_model.rds`, `tables/test_set_selected_model_results.csv`, `tables/test_set_ite_strata_effects.csv`, `tables/full_cohort_refit_selected_model.csv`, `figures/fig1_updated_flowchart.png` – `fig6_*.png`, `figures/supp_*.png` |

### Report

`report/final_report_unified_ml_selection.Rmd` generates the full HTML report from the outputs of Layer 2 (scripts 12–14). It loads `selected_unified_model.rds` and all tables produced by scripts 13–14.

---

## Candidate Models (script 13)

Five models compared by 5-fold CV in the derivation set:

1. **Penalized logistic regression** — outcome model with all covariates, W, and all W×covariate interactions; ITEs = predicted P(AKI|W=1,X) − P(AKI|W=0,X); elastic net α=0.5 (glmnet).
2. **T-learner (XGBoost)** — separate outcome models per arm; ITE = μ₁(X) − μ₀(X).
3. **X-learner (XGBoost)** — two-stage: fits μ₁, μ₀ then pseudo-outcomes D₁, D₀; ITE blended by propensity.
4. **R-learner (manual Robinson decomposition, XGBoost)** — the "Rboost" architecture of Buell et al. JAMA 2024; fits residualised outcome on residualised treatment.
5. **Causal Forest (grf)** — generalised random forest with automatic hyperparameter tuning.

**Selection metrics** (all DR/AIPW-adjusted, out-of-fold):
- Adjusted Qini (primary weight ×3)
- RATE / AUTOC (weight ×2)
- C-for-benefit via PS-matched pairs (weight ×1)

---

## Study Design

- Retrospective new-user, active-comparator cohort study
- Index date (t0): first qualifying ASP or cefazolin dose within a prespecified window around the first positive *S. aureus* blood culture
- New-user design: prior ASP/CEF exposure excluded
- **Primary outcome**: 7-day AKI by KDIGO creatinine criteria (peak ≥ 1.5× baseline, or absolute rise ≥ 0.3 mg/dL)
- **Covariates (N=34)**: demographics, renal function (creatinine, eGFR), severity markers, comorbidities (ICD-coded), ICU status, vasopressors, mechanical ventilation, nephrotoxic co-exposures, missingness indicators for bilirubin and lactate
- **Missing data**: median imputation (continuous) / mode imputation (binary) derived from the derivation set only; missingness indicators for variables with >5% missing

---

## Installation

```r
source("install_r_packages.R")
```

This installs all required packages, including `grf`, `xgboost`, and `glmnet`.

---

## Execution

Scripts must be run in order. Before running, set the `WORK_DIR` variable at the top of each script to your local working directory containing the MIMIC-IV data and derived outputs:

```r
WORK_DIR <- "/path/to/your/working/directory"
```

Then source each script in sequence:

```r
# Layer 1 — data and initial analysis
source("R/01_build_mssa_cohort.R")
source("R/02_define_treatment_and_time_zero.R")
source("R/03_build_outcomes_aki.R")
source("R/04_build_covariates.R")
source("R/05_analysis_descriptive.R")
source("R/06_primary_causal_grf.R")
source("R/07_sensitivity_analyses.R")
source("R/08_generate_tables_figures.R")
source("R/09_render_report.R")
source("R/10_reviewer_requested_supplementary_analyses.R")

# Layer 2 — unified derivation/test pipeline
source("R/12_unified_derivation_test_pipeline.R")
source("R/13_candidate_model_selection_unified.R")
source("R/14_selected_model_primary_and_ite_evaluation.R")

# Final report
rmarkdown::render("report/final_report_unified_ml_selection.Rmd")
```

> **Note**: Scripts 06–10 are retained for reproducibility of the initial analysis and reviewer responses. The primary results are those produced by scripts 12–14. There is no script 11 in this repository.

---

## Data Access

This repository does **not** include MIMIC-IV source data or any derived data files. Access to MIMIC-IV requires completion of a data use agreement through [PhysioNet](https://physionet.org/content/mimiciv/).

Expected local structure under `WORK_DIR`:

```
<WORK_DIR>/
└── mimic-iv-3.1/
    ├── hosp/
    │   ├── microbiologyevents.csv.gz
    │   ├── admissions.csv.gz
    │   ├── patients.csv.gz
    │   ├── emar.csv.gz
    │   ├── labevents.csv.gz
    │   ├── diagnoses_icd.csv.gz
    │   ├── procedures_icd.csv.gz
    │   └── prescriptions.csv.gz
    └── icu/
        └── icustays.csv.gz
```

---

## Methodological References

- Buell KG, Spicer AB, Casey JD, et al. Individualized treatment effects of oxygen targets in critically ill adults. *JAMA* 2024;331(14):1195–1204.
- Munroe ES, Spicer A, Castellvi-Font A, et al. A framework for deriving and validating individualized treatment effects in critical care. *Lancet Respir Med* 2025;13:556–568.
- VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med* 2017;167(4):268–274.
- KDIGO Clinical Practice Guideline for Acute Kidney Injury. *Kidney Int Suppl* 2012;2:1–138.
