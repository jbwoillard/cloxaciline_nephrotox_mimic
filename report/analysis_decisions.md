# Analysis Decisions Log

---

## Internal review — supplementary analyses (2026-03-12)

### Decision S1: Missing-data robustness (MIA)
- Original strategy: median imputation + binary missingness indicators for vars with >5% missing
- Reviewer request: test grf native MIA handling as robustness check
- MIA result: ATE_overlap = 0.181 (95% CI: 0.091–0.272) vs. 0.216 (0.115–0.318) original
- Decision: primary analysis unchanged. MIA analysis reported as supplementary S1.
  Results are very similar; no reason to replace the original strategy.

### Decision S2: E-value
- EValue R package not installed; E-value computed manually using VanderWeele & Ding (2017) formula
- Overlap-weighted adjusted RR derived from primary causal forest propensity scores
- Adjusted RR = 1.50 (95% CI: 1.17–1.90)
- E-value (point estimate): 2.36; E-value (CI lower): 1.61
- Interpretation: unmeasured confounding would need to be strong to negate the finding.

### Decision S3: Positivity
- PS range: [0.208, 0.484]; no values below 0.10 or above 0.90
- Narrow PS range due to class imbalance (69% cefazolin), NOT positivity violation
- Overlap assumption is met; overlap estimand is appropriate

### Reproducibility note
- All supplementary analyses use seed=42, same as primary analysis
- Script: R/10_reviewer_requested_supplementary_analyses.R
- Existing primary results NOT modified
