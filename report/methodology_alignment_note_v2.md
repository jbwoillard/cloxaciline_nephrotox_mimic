# Methodology Alignment Note v2

**Project:** Target Trial Emulation — ASP vs Cefazolin Nephrotoxicity in MSSA Bacteremia (MIMIC-IV v3.1)
**Date:** 2026-03-27
**Author:** Cyrielle (analysis by Claude Code)

---

## 1. Key Methods from Buell et al. JAMA 2024 (Paper 1)

**Full citation:** Buell KG, Spicer AB, Casey JD, et al. Individualized Treatment Effects of Oxygen Targets in Mechanically Ventilated Critically Ill Adults. *JAMA*. 2024;331(14):1195–1204.

### Study design
- Secondary analysis of two temporally and geographically distinct randomized trials (PILOT, n=1682; ICU-ROX, n=965).
- Model derivation in PILOT; external validation in ICU-ROX.
- Effect-based analysis of heterogeneity of treatment effect (HTE), following the PATH statement framework.

### Candidate model evaluation
- Six machine learning algorithms evaluated in the PILOT derivation cohort using **5-fold cross-validation**.
- Evaluation metric: **mean out-of-sample adjusted qini statistic** (primary discriminability metric).
- The best-performing algorithm (Rboost, an XGBoost implementation of the R-learner) was selected for full-cohort derivation and external validation.
- Five Rboost models were constructed with different seed initializations for stability.

### Model selection metrics
- **Adjusted Qini**: gain in survival from patients being ranked by predicted ITE (ascending benefit order); area above the diagonal in the uplift curve.
- **C-for-benefit**: probability of concordance between predicted and observed benefit in matched pairs; C > 0.5 indicates discrimination better than random.
- Both metrics computed in the **held-out** ICU-ROX validation cohort, not the derivation cohort.

### HTE strata reporting
- Patients ranked by predicted ITE and categorized into **tertiles** (lower, middle, upper thirds) using cutoffs determined in the validation cohort.
- Within-tertile observed absolute risk difference reported with 95% confidence intervals.
- **Likelihood ratio test** for interaction between tertile and randomized group to formally test effect modification.
- Calibration assessed by comparing mean predicted ITE to observed risk difference within each tertile.

### Imputation
- Missing predictor data imputed using **bagged trees derived from the derivation cohort only** (PILOT), then applied to the validation cohort — no leakage.

### Cautious language
- "Predicted individualized treatment effects" (not "true ITEs") throughout.
- Explicitly stated: true individual treatment effects are unmeasurable (counterfactual problem).
- External prospective validation called for before clinical use.

---

## 2. Key Methods from Munroe / Spicer / Castellvi-Font et al. Lancet Respiratory Medicine 2025 (Paper 2)

**Full citation:** Munroe ES, Spicer A, Castellvi-Font A, et al. Evidence-based personalised medicine in critical care: a framework for quantifying and applying individualised treatment effects in patients who are critically ill. *Lancet Respir Med*. 2025;13:556–568.

### Framework overview
- Comprehensive review of approaches to HTE in critical care: subgroup, data-derived subgroups, risk-based, and effect-based modelling.
- Effect-based (ITE) models carry high risk of overfitting and spurious HTE detection — rigorous derivation and validation are essential.
- Prespecification of the approach, included variables, and rationale for expected HTE is recommended.

### Derivation/test separation principle
- Model building (derivation) and model testing (validation) must be conducted in separate datasets.
- In trial settings, derivation and external validation in distinct trials is the gold standard (as in Buell et al.).
- When a single dataset is used, a held-out test set must be locked until after model selection.
- Test set evaluation must not inform any modelling decisions.

### Candidate model comparison
- Multiple effect-based models compared using out-of-sample performance metrics.
- The model with the best performance metrics in derivation is selected; performance is then assessed in the locked test set.
- Metrics used: adjusted qini (discrimination of ITE ordering), RATE/AUTOC (rank-weighted ATE), C-for-benefit (concordance of predicted vs observed benefit).

### HTE strata
- Patients stratified into tertiles or quartiles by predicted ITE.
- Observed treatment effects (absolute risk differences, 95% CIs) compared across strata.
- Formal test of effect modification across strata.
- Calibration plot: mean predicted ITE per tertile vs observed DR ATE per tertile.

### Calibration and discrimination
- Calibration: mean predicted ITE per stratum vs observed stratum-level DR ATE.
- Discrimination: adjusted qini > 0 and C-for-benefit > 0.5 constitute evidence of meaningful discrimination.
- RATE (Rank Average Treatment Effect / AUTOC): positive value indicates patients ranked higher by predicted ITE respond more to treatment.

### Observational language
- HTE analyses identify *associations* between characteristics and differential treatment effects, not causation.
- Predicted ITEs are approximations (CATEs); true individual treatment effects remain unmeasurable.
- Residual confounding must be acknowledged explicitly in observational studies.

---

## 3. What Was Borrowed Directly for the Present Study

| Element | Source | Implementation |
|---|---|---|
| Effect-based ITE modelling | Buell et al., Munroe et al. | Causal ML models (T-, X-, R-learner, grf, penalized logistic) |
| Derivation/test split before any modelling | Both papers | 70/30 stratified split (by treatment x AKI status), seed=42, locked test set |
| 5-fold CV in derivation for model selection | Buell et al. | 5 folds, stratified by treatment x AKI |
| Adjusted qini as primary selection metric | Buell et al. | DR-adjusted qini (doubly robust pseudo-outcomes) |
| RATE/AUTOC as secondary metric | Munroe et al. appendix | Manual AUTOC (cumulative mean DR effect above overall mean) |
| C-for-benefit | Buell et al., Munroe et al. | PS-matched pairs, observational analogue |
| Tertile strata analysis | Buell et al. | Tertiles by predicted ITE in test set; DR ATE per tertile with 95% CI |
| Calibration plot | Buell et al. | Mean predicted ITE vs observed DR ATE per tertile |
| Imputation rules from derivation only | Buell et al. | Median/mode imputation fit on derivation, applied to test without leakage |
| Bagged trees for imputation | Buell et al. | Adapted: median imputation with missingness indicators (smaller N) |
| Cautious observational language | Munroe et al. | "Predicted ITE", "CATE approximation", counterfactual limitation stated |

---

## 4. What Was Adapted for Observational Data (Not RCT)

### No randomization — propensity score required
- In the RCT setting, propensity is known by design (0.5 for 1:1 randomization). In MIMIC-IV, treatment (ASP vs cefazolin) was prescribed by clinicians, creating potential confounding by indication.
- A propensity score model (logistic regression) was estimated from the derivation set covariates.
- All DR pseudo-outcomes use this estimated propensity — the doubly robust (AIPW) formula accounts for confounding by augmenting the outcome difference with an inverse-probability weight.

### DR pseudo-outcomes replace direct outcome differences
- In an RCT, within-stratum risk differences equal the ATE within that stratum. In observational data, DR pseudo-outcomes are used throughout (qini, RATE, strata calibration) to adjust for confounding.

### Overlap weights instead of ATE weights
- The primary causal estimate uses overlap (ATO) weighting rather than IPTW, to target the estimand for the overlap population — patients for whom there is equipoise regarding treatment choice. This is the clinically relevant population.
- This is consistent with the target trial emulation framework and avoids extreme weights.

### No external validation cohort
- Buell et al. used two distinct trials. Here, a single observational cohort (MIMIC-IV) is available. Internal validation using a locked hold-out set is the maximum achievable without external data.
- This is explicitly acknowledged as a limitation.

### Sample size considerations
- PILOT: n=1682; ICU-ROX: n=965. Present study: N=445 (analyzable). This is substantially smaller.
- Regularized models (penalized logistic, shallow XGBoost trees, tuned grf) are preferred over unconstrained models to reduce overfitting.
- Imputation uses simpler median rules (bagged trees would overfit with 34 covariates and N~310 derivation observations).

### Conservative interpretation of CATE heterogeneity
- With N=445, statistical power to detect HTE is low. All HTE findings are described as hypothesis-generating.
- The primary output is the DR ATE in the test set; ITE strata analyses are secondary.

---

## 5. Why the Previous Mixed-Framework Was Replaced by a Unified One

### Problem with the previous framework
The earlier analysis (scripts 06–10) mixed two conceptually distinct frameworks:

1. **Primary ATE**: Estimated using `grf::causal_forest` on the full cohort (no derivation/test split), with overlap weighting. This was ATE = 0.181 [0.090, 0.271].

2. **ITE/CATE analysis**: Estimated separately using the same grf model but with CATE predictions from the same training data — this constitutes in-sample evaluation of predicted ITEs, which is known to produce optimistic discrimination metrics.

3. **Supplementary analyses (script 10)**: Added a separate MIA (missing indicator augmentation) comparison, a separate R-learner comparison, and manually computed E-values — all appended post-hoc without a unified framework.

### Problems identified
- **No model selection**: The grf model was used as the single model without comparison to alternatives using out-of-sample metrics.
- **No locked test set**: ITE predictions were evaluated on the same data used for fitting — circular assessment.
- **Inconsistent ATE**: The `ate_results.csv` reported ATE=0.216 (train-set estimate), while `sensitivity_results.csv` reported 0.181 (full-cohort). These discrepancies arose from the mixed framework.
- **No calibration**: The calibration of predicted ITEs against observed outcomes in held-out data was never assessed.

### Solution: unified framework
The unified framework, following Buell et al. and Munroe et al., implements:
1. A single 70/30 derivation/test split at the outset.
2. Imputation rules derived from the derivation set only.
3. Five-fold CV in the derivation set to select among 5 candidate causal models.
4. Model selection by weighted rank on adjusted qini, RATE, and C-for-benefit.
5. Selected model refit on the full derivation set; evaluated once on the locked test set.
6. Both the primary DR ATE and all ITE/HTE metrics come from the same selected model applied to the same test set.
7. Post-selection, the selected model is refit on all 445 analyzable patients for a secondary full-cohort ATE estimate — clearly labelled as post-validation and not the primary estimate.

This approach eliminates the circular evaluation problem, provides consistent estimates, and aligns with the published methodological framework for ITE analysis in critical care.
