# Primary ATE Audit
Date: 2026-03-12

## Summary of findings

There are **three distinct ATE values** in the project, arising from two different
model objects. Only one is the true primary estimate.

| Model | Object | N | Strategy | ATE (overlap) | SE | 95% CI | Primary? |
|-------|--------|---|----------|:---:|:---:|:---:|:---:|
| `cf_train` | `mobj$ate_overlap` → `ate_results.csv` | 355 | Median imputation — **TRAIN SET only** | 0.216 | 0.052 | [0.115, 0.318] | **NO** |
| `cf_full` | `mobj$ate_full` → `sensitivity_results.csv` | 445 | Median imputation — **FULL dataset** | 0.181 | 0.046 | [0.090, 0.271] | **YES** |
| `cf_mia_full` | `mia_obj$mia_full_overlap` | 445 | grf native MIA — **FULL dataset** | 0.181 | 0.046 | [0.091, 0.272] | Supplement only |

## Root cause of the 0.216 value

In `R/06_primary_causal_grf.R`:

```r
# Step 3: fit on training set
cf <- causal_forest(X=X_train, W=W_train, Y=Y_train, ...)  # N=355

# Step 4: compute ATE → stored in ate_results.csv
ate_overlap <- average_treatment_effect(cf, target.sample="overlap")  # 0.216
ate_results <- data.frame(estimand=..., estimate=c(...ate_overlap...))
fwrite(ate_results, "results/ate_results.csv")  # ← SAVES 0.216

# Step 10: refit on full data → TRUE PRIMARY
cf_full <- causal_forest(X=X_use, W=W_use, Y=Y_use, ...)  # N=445
ate_full <- average_treatment_effect(cf_full, target.sample="overlap")  # 0.181
# ate_full is NOT written to ate_results.csv — only saved in model_objects.rds
```

**The bug**: `ate_results.csv` was written before the full-dataset refit in step 10.
It contains the train-set ATE (0.216), not the full-dataset ATE (0.181).

## What each downstream file uses

| File | Value used | Source | Correct? |
|------|-----------|--------|----------|
| `report/final_report.md` | **0.181** | Human writing (from sensitivity script output) | ✓ YES |
| `results/sensitivity_results.csv` primary row | **0.181** | `model_objs$ate_full` in script 07 | ✓ YES |
| `results/ate_results.csv` | 0.216 | `cf_train` ATE | ✗ WRONG — should be 0.181 |
| `tables/supplement_mia_vs_primary.csv` "Original" row | 0.216 | `mobj$ate_overlap` in script 10 | ✗ WRONG — should be 0.181 |

## What the Methods text says

`report/final_report.Rmd` Methods section:
> "Missing data: median imputation + missingness indicators"
> "Causal forest (grf package), 2000 trees, tuned parameters"

This describes **median imputation + indicators** as the primary strategy, i.e. `cf_full`.
The primary estimate is therefore the `cf_full` result = **0.181**.

## Verdict

- **Primary ATE = 0.181** (95% CI: 0.090–0.271, p=0.0001), from `cf_full`, N=445
- **MIA sensitivity ATE = 0.181** (95% CI: 0.091–0.272, p<0.0001), from `cf_mia_full`, N=445
- The two are virtually identical — the MIA sensitivity analysis confirms complete robustness
- **0.216 should not appear as a primary result anywhere** — it is a train-set artifact
