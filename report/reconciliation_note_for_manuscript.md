# Reconciliation Note for Manuscript
Date: 2026-03-12

---

## Q1. Why is the cohort 446 in some places and 445 in others?

**446** = total patients enrolled (met all eligibility criteria). This is the correct number for
the flowchart, Table 1, and any sentence describing cohort construction.

**445** = patients included in the **primary analysis**. One patient (subject_id 12942901,
cefazolin arm) had no creatinine measurements whatsoever (neither baseline nor follow-up),
making it impossible to define the primary outcome (AKI at 7 days). This patient is correctly
excluded from all model fitting by the line `complete_rows <- !is.na(anal$AKI_7d)` in
`R/06_primary_causal_grf.R`.

**Fix**: The manuscript should consistently say:
> "446 patients met inclusion criteria; 445 were included in the primary analysis
> (1 excluded: no creatinine data available to define the AKI outcome)."

---

## Q2. Which ATE is the true primary estimate?

**ATE = 0.181 (95% CI: 0.090–0.271, p=0.0001)** is the true primary estimate.

This comes from `cf_full`, the causal forest refit on the **full analyzed cohort (N=445)**
using **median imputation + missingness indicators** for missing covariates, as described
in the Methods. This value is already correctly stated in `report/final_report.md` and used
in `results/sensitivity_results.csv`.

---

## Q3. Which ATE is the missing-data sensitivity analysis?

**ATE = 0.181 (95% CI: 0.091–0.272, p<0.0001)** from `cf_mia_full` using grf native
MIA (Missing Incorporated in Attributes), N=445. This is the supplementary robustness check.

**The two estimates are virtually identical** (difference < 0.001), confirming complete
robustness of the primary finding to the missing-data strategy.

---

## Q4. What must be updated for consistency?

### Files that are CORRECT and need no changes
| File | Value | Status |
|------|-------|--------|
| `report/final_report.md` | ATE=0.181, primary N=445 | ✓ Correct |
| `results/sensitivity_results.csv` | primary ATE=0.181 | ✓ Correct |
| `results/model_objects.rds` | ate_full=0.181 | ✓ Correct |

### Files that need patching

#### 1. `results/ate_results.csv`
Currently stores train-set ATE (0.216). Must be replaced with full-dataset ATE (0.181):
- ATE_global: 0.219 → replace with cf_full global ATE
- ATT_treated: 0.216 → replace with cf_full treated ATE
- ATE_overlap: 0.216 → **replace with 0.181** (cf_full, ate_full)

#### 2. `tables/supplement_mia_vs_primary.csv`
"Original" row currently shows 0.216 (train-set artifact). Must be replaced with 0.181
(cf_full, N=445). After patch, both rows will show ~0.181, which is the correct conclusion:
the MIA sensitivity result is nearly identical to the primary.

#### 3. `report/final_report.Rmd` — Supplementary section S1
The comparison table in S1 shows "Original: 0.216" — must be corrected to "Original: 0.181".
The interpretation paragraph must be updated accordingly (the difference between strategies
is now negligible, <0.001 rather than ~0.035).

#### 4. Flowchart / cohort description
Anywhere the text says "N=446" for the analysis population, add clarification:
- "446 enrolled, **445 analyzed** (1 excluded: missing creatinine data)"

### Items that do NOT need changes
- Primary ATE in abstract, results, and discussion (already 0.181) ✓
- Sensitivity analyses table (already uses 0.181) ✓
- E-value (computed from MIA model propensity scores, not affected) ✓
- Overlap diagnostics (from cf_full propensities, not affected) ✓
