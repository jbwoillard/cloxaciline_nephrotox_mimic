# Cohort Size Audit
Date: 2026-03-12

## Summary
- `analysis_dataset.rds`: **446 rows** (enrolled cohort, as reported in flowchart)
- Models (primary grf, MIA, sensitivity): **445 rows** (analyzed cohort, complete outcome)
- Discrepancy: **1 patient** excluded due to missing AKI outcome

## Root cause
In `R/06_primary_causal_grf.R`, line:
```r
complete_rows <- !is.na(anal$AKI_7d)
```
This filters out 1 patient with `NA` for `AKI_7d`.

## Discrepant patient
| subject_id | hadm_id | W | AKI_7d | creat_baseline | delta_creat_7d |
|------------|---------|---|--------|----------------|----------------|
| 12942901 | 22378688 | 0 (CEF) | NA | NA | NA |

This patient had **no creatinine measurements** at all (neither baseline nor follow-up), making it impossible to define the AKI outcome. Exclusion is clinically appropriate and methodologically required.

## Verdict
- 446 = enrolled cohort (appears in flowchart, Table 1 header, cohort steps)
- 445 = analyzed cohort (N used in all model fits, sensitivity analyses, and report ATE estimates)
- The report must consistently distinguish these two numbers:
  - "446 patients met all eligibility criteria and were enrolled"
  - "445 patients were included in the primary analysis (1 excluded: missing creatinine data)"

## Files confirming N=445 for models
- `results/model_objects.rds`: `length(W_use) = 445`
- `results/mia_model_objects.rds`: `length(cf_mia_full$W.orig) = 445`
- `results/sensitivity_results.csv`: N = 445
- `derived/analysis_dataset_with_cate.rds`: 445 rows
