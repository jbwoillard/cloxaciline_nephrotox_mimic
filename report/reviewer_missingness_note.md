# Missing Data Audit — Internal Review
Date: 2026-03-12

## Summary
Analysis dataset: N=445 (ASP: 138, CEF: 307)

## Variables with missing data

| Variable | Overall % missing | ASP % missing | CEF % missing | Handling |
|----------|------------------|---------------|---------------|---------|
| lactate | 59.1% | 54.3% | 61.2% | Median imputation + missingness indicator |
| bilirubin | 40.0% | 32.6% | 43.3% | Median imputation + missingness indicator |
| platelets |  2.9% |  1.4% |  3.6% | Median imputation only |
| bun |  2.9% |  2.9% |  2.9% | Median imputation only |
| creat_baseline |  2.5% |  2.2% |  2.6% | Median imputation only |
| egfr_baseline |  2.5% |  2.2% |  2.6% | Median imputation only |
| wbc |  2.5% |  0.7% |  3.3% | Median imputation only |
| glucose |  2.5% |  2.2% |  2.6% | Median imputation only |
| bicarbonate |  2.5% |  2.2% |  2.6% | Median imputation only |
| sodium |  2.5% |  2.2% |  2.6% | Median imputation only |
| potassium |  2.5% |  2.2% |  2.6% | Median imputation only |

## Key observations
- Lactate and bilirubin have the highest missingness (59% and 40% respectively).
- These two variables received both median imputation AND a binary missingness indicator
  in the original causal forest model.
- All other covariates had <5% missing and received median imputation only.
- grf natively supports missing values via MIA (Missing Incorporated in Attributes),
  which can handle missingness more flexibly than simple median imputation.
- A supplementary MIA-based analysis was run as Task B.
