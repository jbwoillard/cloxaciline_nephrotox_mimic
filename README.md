# MSSA Bacteremia — ASP vs. Cefazolin Nephrotoxicity Study (MIMIC-IV v3.1)

Target trial emulation reproducing the CloCeBa Lancet signal using MIMIC-IV.

## Scientific Question

Does anti-staphylococcal penicillin (ASP: nafcillin/oxacillin) cause more AKI than cefazolin in MSSA bacteremia, as observed in the CloCeBa trial?

## Key Result

**ATE (overlap-weighted) = +0.181 (95% CI: 0.090–0.271, p=0.0001)** ASP associated with \~18 percentage points higher absolute AKI risk vs. cefazolin.

## How to Run the Pipeline

``` bash
# 1. Set up environment (R 4.4+, packages below)
cd /path/to/cyrielle_mimic_cloxa

# 2. Run each script in order
Rscript R/00_audit_data_and_drug_feasibility.R
Rscript R/01_build_mssa_cohort.R
Rscript R/02_define_treatment_and_time_zero.R
Rscript R/03_build_outcomes_aki.R
Rscript R/04_build_covariates.R
Rscript R/05_analysis_descriptive.R
Rscript R/06_primary_causal_grf.R
Rscript R/07_sensitivity_analyses.R
Rscript R/08_generate_tables_figures.R
Rscript R/09_render_report.R
```

## R Package Dependencies

``` r
install.packages(c(
  "data.table", "duckdb", "DBI", "dplyr", "tidyr", "stringr",
  "lubridate", "ggplot2", "grf", "tableone", "cobalt",
  "rmarkdown", "broom", "purrr", "janitor"
))
```

## Data

-   MIMIC-IV v3.1 (all tables in `mimic-iv-3.1/hosp/` and `mimic-iv-3.1/icu/`, `.csv.gz` format)
-   Reference paper: PIIS0140673625016241.pdf (Lancet CloCeBa trial)

## Results Location

| Output                  | Location                          |
|-------------------------|-----------------------------------|
| Analysis dataset        | `derived/analysis_dataset.rds`    |
| Causal forest model     | `results/model_objects.rds`       |
| ATE estimates           | `results/ate_results.csv`         |
| CATE predictions        | `results/cate_predictions.csv`    |
| Sensitivity analyses    | `results/sensitivity_results.csv` |
| Table 1                 | `tables/table1.csv`               |
| All figures (10 PNGs)   | `figures/`                        |
| Final report (RMarkdown) | `report/final_report.Rmd`          |

## Cohort Summary

| Step                      | N       |
|---------------------------|---------|
| SA blood cultures         | 2,223   |
| MSSA (exclude MRSA)       | 1,001   |
| With qualifying treatment | 458     |
| New-user (final cohort)   | **446** |

-   ASP arm (nafcillin/oxacillin): n=138
-   Cefazolin arm: n=308

## Design

-   New-user, active-comparator cohort study
-   t0 = first ASP or cefazolin administration within [-48h, +5d] of positive blood culture
-   AKI: KDIGO creatinine criteria (peak ≥1.5× baseline OR Δ≥0.3 mg/dL)
-   Causal method: Generalized Random Forest (`grf::causal_forest`)
-   33 pre-t0 covariates, seed=42
