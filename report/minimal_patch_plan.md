# Minimal Patch Plan
Date: 2026-03-12

## Overview
Two files contain incorrect values (train-set ATE 0.216 instead of full-data ATE 0.181).
The patch is minimal: replace those values, update the supplement comparison text.
The main manuscript text is already correct.

## Patch 1: results/ate_results.csv
**Replace** with cf_full full-dataset estimates (from model_objects.rds$ate_full).
- ATE_overlap: 0.216285 → **0.18071**
- SE: 0.05175888 → **0.04615949**
- CI_lower: 0.1148376 → **0.09023736**
- CI_upper: 0.3177324 → **0.27118264**
- p_value: 2.93e-05 → **9.04e-05**
(ATE_global and ATT_treated from cf_full also need updating — see patch script)

## Patch 2: tables/supplement_mia_vs_primary.csv
**Replace** "Original" row ATE from 0.2163 → **0.1807**
After patch, both rows will read ~0.181, correctly showing that MIA and primary agree.

## Patch 3: report/final_report.Rmd — supplement section S1
Update comparison table and interpretation paragraph:
- Change "Original: 0.216 [0.115, 0.318]" → "Original: 0.181 [0.090, 0.271]"
- Update interpretation: the two strategies give virtually identical results (<0.001 difference)

## Patch 4: Cohort N clarification
Anywhere "N=446 patients" is used for the analytical model population, clarify:
"446 enrolled, 445 analyzed (1 excluded: no creatinine data to define AKI outcome)"
Key locations: Methods section, Results cohort table, any in-text ATE statements.

## What does NOT change
- Primary ATE in abstract, main results, discussion: already 0.181 ✓
- Sensitivity analyses table: already uses 0.181 ✓
- E-value results: derived from MIA propensity scores, correct ✓
- Overlap/positivity tables: based on cf_full propensities, correct ✓
- Figures: CATE distribution, love plot, sensitivity forest — all correct ✓
