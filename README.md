# MSSA Bacteremia — ASP vs Cefazolin Nephrotoxicity Study (MIMIC-IV v3.1)

Causal machine-learning analysis of nephrotoxicity associated with anti-staphylococcal penicillins (ASP: nafcillin/oxacillin) versus cefazolin in MSSA bacteremia using MIMIC-IV.

## Scientific Question

Among adults with MSSA bacteremia, is treatment with anti-staphylococcal penicillins associated with a higher short-term risk of acute kidney injury (AKI) than cefazolin, and can treatment-effect heterogeneity be identified using a derivation–test individualized treatment-effect framework?

## Current Analysis Framework

This repository now uses a **unified derivation–test causal machine-learning pipeline**.

Key features:
- same clinical cohort and outcome definitions as the initial target trial emulation
- strict **derivation (70%) / held-out test (30%)** split
- multiple candidate causal/effect-based models compared in derivation
- model selection using:
  - adjusted qini
  - RATE
  - C-for-benefit
- only the **selected best model** is evaluated in the held-out test set
- the same selected model is used for both:
  - average treatment effect (ATE)
  - individualized treatment effects (ITE)

## Main Findings

### Primary validated estimate
In the held-out test set, the primary doubly robust average treatment effect (absolute risk difference for ASP vs cefazolin on 7-day AKI) was:

**DR ATE = 0.153 (95% CI −0.183 to 0.416)**

### Secondary post-selection full-cohort refit
After model selection, refitting the selected model on the full analyzable cohort yielded:

**DR ATE = 0.160 (95% CI 0.069 to 0.267)**

### Selected model
The best-performing model in the derivation set was:

**Penalized logistic regression with treatment–covariate interactions**

### Interpretation
The updated unified analysis remains directionally consistent with higher AKI risk under ASP than under cefazolin, but the primary held-out test-set estimate is imprecise. Predicted individualized treatment effects were all positive, suggesting variation in the magnitude of harm rather than a subgroup with predicted renal benefit from ASP.

## Cohort Summary

| Step | N |
|---|---:|
| *S. aureus* blood cultures | 2,223 |
| After exclusion of oxacillin-resistant isolates | 1,001 |
| With qualifying treatment | 458 |
| Eligible descriptive cohort | 446 |
| Analyzable cohort (non-missing primary outcome) | 445 |

Treatment groups:
- ASP (nafcillin/oxacillin): n=138
- Cefazolin: n=307 in analyzable cohort

## Study Design

- Retrospective new-user, active-comparator cohort study
- Index date = first ASP or cefazolin administration within the prespecified treatment window around the positive blood culture
- Primary outcome = 7-day AKI using KDIGO creatinine criteria
- Covariates = pre-treatment demographics, renal function, severity markers, comorbidities, and nephrotoxic co-exposures
- Effect-based model selection performed in derivation only
- Final evaluation performed in a locked held-out test set

## Main Scripts

| Script | Purpose |
|---|---|
| `R/13_candidate_model_selection_unified.R` | Candidate model comparison in derivation |
| `R/14_selected_model_primary_and_ite_evaluation.R` | Held-out test-set evaluation of the selected model |

## Main Outputs

| Output | Location |
|---|---|
| Final unified report (R Markdown) | `report/final_report_unified_ml_selection.Rmd` |
| Final unified report (HTML) | `report/final_report_unified_ml_selection.html` |
| Methodology alignment note | `report/methodology_alignment_note_v2.md` |
| Updated flowchart | `figures/fig1_updated_flowchart.png` |
| Crude AKI rates | `figures/fig2_crude_aki_by_treatment.png` |
| Candidate model selection | `figures/fig3_candidate_model_selection.png` |
| ITE distribution in test set | `figures/fig4_ite_distribution_test.png` |
| ITE tertile effects in test set | `figures/fig5_ite_strata_test.png` |
| Calibration and qini in test set | `figures/fig6_calibration_and_qini_test.png` |

## Data access

This repository does **not** include MIMIC-IV source data or local derived data files. Users must obtain access to MIMIC-IV independently through PhysioNet and recreate the analysis environment locally.

## Notes

This repository reflects the **final unified derivation–test analysis**, which supersedes the earlier full-cohort generalized random forest workflow.
