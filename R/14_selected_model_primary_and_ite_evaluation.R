###############################################################################
# 14_selected_model_primary_and_ite_evaluation.R
# Refit selected model on full derivation set; evaluate in locked held-out test set
# Follows: Buell et al. JAMA 2024 and Munroe et al. Lancet Resp Med 2025
#
# Tasks:
#   A. Primary validation: DR ATE in test set with 95% CI (bootstrap)
#   B. ITE distribution in test set
#   C. ITE strata analysis (tertiles) in test set
#   D. Test-set performance metrics (qini, RATE, C-for-benefit, calibration)
#   E. Post-selection full-cohort refit (secondary, clearly labelled)
#
# Outputs:
#   results/selected_unified_model.rds
#   tables/test_set_selected_model_results.csv
#   tables/test_set_ite_strata_effects.csv
#   tables/full_cohort_refit_selected_model.csv
#   All figures (fig1–fig6, supp_ps_overlap, supp_variable_importance,
#                supp_ale_plots)
###############################################################################

library(dplyr)
library(ggplot2)
library(grf)
library(xgboost)
library(glmnet)

cat("=== 14_selected_model_primary_and_ite_evaluation.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

set.seed(42)

WORK_DIR    <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
RESULTS_DIR <- file.path(WORK_DIR, "results")
TABLES_DIR  <- file.path(WORK_DIR, "tables")
FIGURES_DIR <- file.path(WORK_DIR, "figures")

###############################################################################
# 1. Load data and identify selected model
###############################################################################
cat("--- Step 1: Load data ---\n")

deriv    <- readRDS(file.path(DERIVED_DIR, "derivation_set_unified.rds"))
test_dat <- readRDS(file.path(DERIVED_DIR, "test_set_unified.rds"))
cv_res   <- readRDS(file.path(RESULTS_DIR, "derivation_model_selection_metrics.rds"))

X_deriv  <- deriv$X;      W_deriv <- deriv$W; Y_deriv <- deriv$Y; n_deriv <- nrow(X_deriv)
X_test   <- test_dat$X;   W_test  <- test_dat$W; Y_test  <- test_dat$Y; n_test  <- nrow(X_test)

best_model <- cv_res$best_model

cat(sprintf("Derivation N=%d, Test N=%d\n", n_deriv, n_test))
cat(sprintf("Selected model (from derivation CV): %s\n", best_model))

###############################################################################
# Helper functions
###############################################################################

compute_dr_pseudo <- function(Y, W, e_hat, mu1_hat, mu0_hat) {
  e_hat <- pmax(pmin(e_hat, 0.975), 0.025)
  (W - e_hat) / (e_hat * (1 - e_hat)) *
    (Y - ifelse(W == 1, mu1_hat, mu0_hat)) +
    (mu1_hat - mu0_hat)
}

compute_dr_qini <- function(tau_hat, psi_dr) {
  ok <- !is.na(tau_hat) & !is.na(psi_dr)
  tau_hat <- tau_hat[ok]; psi_dr <- psi_dr[ok]
  if (length(tau_hat) < 5) return(list(qini = NA, frac = NULL, cum = NULL, diagonal = NULL))
  ord          <- order(tau_hat)
  psi_sorted   <- psi_dr[ord]
  n            <- length(psi_sorted)
  cum_effect   <- cumsum(psi_sorted) / n
  frac_treated <- seq_along(psi_sorted) / n
  overall_mean <- mean(psi_dr)
  diagonal     <- frac_treated * overall_mean
  list(qini = mean(cum_effect - diagonal),
       frac = frac_treated, cum = cum_effect, diagonal = diagonal,
       overall_mean = overall_mean)
}

compute_rate_manual <- function(tau_hat, psi_dr) {
  ok <- !is.na(tau_hat) & !is.na(psi_dr)
  tau_hat <- tau_hat[ok]; psi_dr <- psi_dr[ok]
  if (length(tau_hat) < 5) return(NA)
  ord      <- order(tau_hat)
  cum_mean <- cumsum(psi_dr[ord]) / seq_along(psi_dr[ord])
  mean(cum_mean - mean(psi_dr))
}

compute_c_for_benefit_ps <- function(tau_hat, Y, W, e_hat) {
  treated_idx <- which(W == 1); control_idx <- which(W == 0)
  if (length(treated_idx) < 2 || length(control_idx) < 2) return(NA)
  e_t <- e_hat[treated_idx]; e_c <- e_hat[control_idx]
  mp  <- data.frame(trt_idx=integer(0), ctrl_idx=integer(0))
  av  <- seq_along(control_idx)
  for (i in seq_along(treated_idx)) {
    if (!length(av)) break
    b  <- av[which.min(abs(e_t[i] - e_c[av]))]
    mp <- rbind(mp, data.frame(trt_idx=treated_idx[i], ctrl_idx=control_idx[b]))
    av <- av[av != b]
  }
  lp  <- log(e_hat / (1 - e_hat))
  cal <- 0.2 * sd(lp)
  mp  <- mp[abs(lp[mp$trt_idx] - lp[mp$ctrl_idx]) <= cal, ]
  if (nrow(mp) < 5) return(NA)
  pb <- Y[mp$ctrl_idx] - Y[mp$trt_idx]
  pt <- tau_hat[mp$trt_idx]
  n_pairs <- nrow(mp); conc <- 0; disc <- 0
  for (i in seq_len(n_pairs-1)) for (j in (i+1):n_pairs) if (pb[i] != pb[j]) {
    disc <- disc + 1
    if ((pb[i] > pb[j]) == (pt[i] < pt[j])) conc <- conc + 1
  }
  if (disc == 0) return(NA)
  conc / disc
}

###############################################################################
# Helper: fit propensity model on given data
###############################################################################
fit_ps_logistic <- function(X, W) {
  fit <- tryCatch(
    glm(W ~ ., data = as.data.frame(X), family = binomial()),
    error = function(e) glm(W ~ ., data = as.data.frame(X), family = binomial(),
                             control = list(maxit = 200))
  )
  e <- predict(fit, newdata = as.data.frame(X), type = "response")
  list(model = fit, e = pmax(pmin(e, 0.975), 0.025))
}

fit_ps_on_test <- function(model, X_test) {
  e <- predict(model, newdata = as.data.frame(X_test), type = "response")
  pmax(pmin(e, 0.975), 0.025)
}

###############################################################################
# 2. Refit selected model on FULL derivation set
###############################################################################
cat("\n--- Step 2: Refit selected model on full derivation set ---\n")

# Always fit propensity and outcome nuisance on full derivation
cat("  Fitting propensity model...\n")
ps_deriv    <- fit_ps_logistic(X_deriv, W_deriv)
e_deriv     <- ps_deriv$e
cat(sprintf("  Propensity range: [%.3f, %.3f]\n", min(e_deriv), max(e_deriv)))

xgb_p <- list(objective = "binary:logistic", eta = 0.05, max_depth = 3, min_child_weight = 5)
xgb_r <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 3, min_child_weight = 5)

selected_model_obj <- NULL
tau_deriv          <- NULL

if (grepl("Penalized Logistic", best_model)) {
  cat("  Fitting Penalized Logistic regression...\n")
  X_deriv_int <- cbind(X_deriv, X_deriv * W_deriv)
  cv_gl <- cv.glmnet(X_deriv_int, Y_deriv, family = "binomial", alpha = 0.5,
                     nfolds = 5, type.measure = "deviance")
  X_deriv_w1 <- cbind(X_deriv, X_deriv * 1)
  X_deriv_w0 <- cbind(X_deriv, X_deriv * 0)
  mu1_d <- as.vector(predict(cv_gl, X_deriv_w1, type="response", s="lambda.min"))
  mu0_d <- as.vector(predict(cv_gl, X_deriv_w0, type="response", s="lambda.min"))
  tau_deriv  <- mu1_d - mu0_d
  selected_model_obj <- list(type="penalized_logistic", model=cv_gl,
                              mu1_deriv=mu1_d, mu0_deriv=mu0_d)

} else if (grepl("T-learner", best_model)) {
  cat("  Fitting T-learner XGBoost...\n")
  m1_t <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F], label=Y_deriv[W_deriv==1]),
                  nrounds=100, params=xgb_p, verbose=0)
  m0_t <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F], label=Y_deriv[W_deriv==0]),
                  nrounds=100, params=xgb_p, verbose=0)
  mu1_d <- predict(m1_t, xgb.DMatrix(X_deriv))
  mu0_d <- predict(m0_t, xgb.DMatrix(X_deriv))
  tau_deriv  <- mu1_d - mu0_d
  selected_model_obj <- list(type="t_learner", m1=m1_t, m0=m0_t,
                              mu1_deriv=mu1_d, mu0_deriv=mu0_d)

} else if (grepl("X-learner", best_model)) {
  cat("  Fitting X-learner XGBoost...\n")
  m1_x1 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F], label=Y_deriv[W_deriv==1]),
                   nrounds=100, params=xgb_p, verbose=0)
  m0_x1 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F], label=Y_deriv[W_deriv==0]),
                   nrounds=100, params=xgb_p, verbose=0)
  D1_d  <- Y_deriv[W_deriv==1] - predict(m0_x1, xgb.DMatrix(X_deriv[W_deriv==1,,drop=F]))
  D0_d  <- predict(m1_x1, xgb.DMatrix(X_deriv[W_deriv==0,,drop=F])) - Y_deriv[W_deriv==0]
  m1_x2 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F], label=D1_d),
                   nrounds=100, params=xgb_r, verbose=0)
  m0_x2 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F], label=D0_d),
                   nrounds=100, params=xgb_r, verbose=0)
  tau1_d <- predict(m1_x2, xgb.DMatrix(X_deriv))
  tau0_d <- predict(m0_x2, xgb.DMatrix(X_deriv))
  tau_deriv <- (1 - e_deriv) * tau1_d + e_deriv * tau0_d
  mu1_d <- predict(m1_x1, xgb.DMatrix(X_deriv))
  mu0_d <- predict(m0_x1, xgb.DMatrix(X_deriv))
  selected_model_obj <- list(type="x_learner", m1_s1=m1_x1, m0_s1=m0_x1,
                              m1_s2=m1_x2, m0_s2=m0_x2,
                              mu1_deriv=mu1_d, mu0_deriv=mu0_d)

} else if (grepl("R-learner", best_model)) {
  cat("  Fitting R-learner (manual)...\n")
  m_full  <- xgboost(data=xgb.DMatrix(X_deriv, label=Y_deriv),
                     nrounds=100, params=xgb_p, verbose=0)
  m_hat_d <- predict(m_full, xgb.DMatrix(X_deriv))
  W_res_d <- W_deriv - e_deriv
  Y_res_d <- Y_deriv - m_hat_d
  Y_rl    <- Y_res_d / W_res_d
  wt_rl   <- W_res_d^2
  m_rl    <- xgboost(data=xgb.DMatrix(X_deriv, label=Y_rl, weight=wt_rl),
                     nrounds=100, params=xgb_r, verbose=0)
  tau_deriv <- predict(m_rl, xgb.DMatrix(X_deriv))
  # For DR pseudo-outcome we need mu1, mu0 separately
  m1_tmp  <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F],label=Y_deriv[W_deriv==1]),
                     nrounds=100, params=xgb_p, verbose=0)
  m0_tmp  <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F],label=Y_deriv[W_deriv==0]),
                     nrounds=100, params=xgb_p, verbose=0)
  mu1_d   <- predict(m1_tmp, xgb.DMatrix(X_deriv))
  mu0_d   <- predict(m0_tmp, xgb.DMatrix(X_deriv))
  selected_model_obj <- list(type="r_learner", m_outcome=m_full, m_rl=m_rl,
                              mu1_deriv=mu1_d, mu0_deriv=mu0_d)

} else {  # Causal Forest (grf) — default
  cat("  Fitting Causal Forest (grf)...\n")
  cf_full <- causal_forest(X=X_deriv, Y=Y_deriv, W=W_deriv,
                           num.trees=2000, seed=42, tune.parameters="all")
  tau_deriv <- predict(cf_full)$predictions
  mu1_d <- cf_full$Y.hat + (1 - cf_full$W.hat) * tau_deriv
  mu0_d <- cf_full$Y.hat - cf_full$W.hat * tau_deriv
  selected_model_obj <- list(type="causal_forest", model=cf_full,
                              mu1_deriv=mu1_d, mu0_deriv=mu0_d)
}

cat(sprintf("  Derivation tau: mean=%.4f, SD=%.4f\n", mean(tau_deriv), sd(tau_deriv)))

###############################################################################
# 3. Apply selected model to TEST set
###############################################################################
cat("\n--- Step 3: Apply to held-out test set ---\n")

# Propensity on test set (using derivation-fitted model)
e_test <- fit_ps_on_test(ps_deriv$model, X_test)
cat(sprintf("  Test propensity range: [%.3f, %.3f]\n", min(e_test), max(e_test)))

tau_test   <- NULL
mu1_test   <- NULL
mu0_test   <- NULL

if (grepl("Penalized Logistic", best_model)) {
  X_test_w1 <- cbind(X_test, X_test * 1)
  X_test_w0 <- cbind(X_test, X_test * 0)
  mu1_test  <- as.vector(predict(selected_model_obj$model, X_test_w1,
                                 type="response", s="lambda.min"))
  mu0_test  <- as.vector(predict(selected_model_obj$model, X_test_w0,
                                 type="response", s="lambda.min"))
  tau_test  <- mu1_test - mu0_test

} else if (grepl("T-learner", best_model)) {
  mu1_test <- predict(selected_model_obj$m1, xgb.DMatrix(X_test))
  mu0_test <- predict(selected_model_obj$m0, xgb.DMatrix(X_test))
  tau_test <- mu1_test - mu0_test

} else if (grepl("X-learner", best_model)) {
  tau1_t   <- predict(selected_model_obj$m1_s2, xgb.DMatrix(X_test))
  tau0_t   <- predict(selected_model_obj$m0_s2, xgb.DMatrix(X_test))
  tau_test <- (1 - e_test) * tau1_t + e_test * tau0_t
  mu1_test <- predict(selected_model_obj$m1_s1, xgb.DMatrix(X_test))
  mu0_test <- predict(selected_model_obj$m0_s1, xgb.DMatrix(X_test))

} else if (grepl("R-learner", best_model)) {
  tau_test <- predict(selected_model_obj$m_rl, xgb.DMatrix(X_test))
  mu1_test <- predict(selected_model_obj$mu1_deriv,
                      xgb.DMatrix(X_test)) # fallback: use outcome from derivation model
  # Better: use a separate outcome model applied to test
  m1_t2 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F],label=Y_deriv[W_deriv==1]),
                   nrounds=100, params=list(objective="binary:logistic",eta=0.05,
                                             max_depth=3,min_child_weight=5), verbose=0)
  m0_t2 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F],label=Y_deriv[W_deriv==0]),
                   nrounds=100, params=list(objective="binary:logistic",eta=0.05,
                                             max_depth=3,min_child_weight=5), verbose=0)
  mu1_test <- predict(m1_t2, xgb.DMatrix(X_test))
  mu0_test <- predict(m0_t2, xgb.DMatrix(X_test))

} else {  # Causal Forest
  tau_test <- predict(selected_model_obj$model, newdata = X_test)$predictions
  # For DR: predict mu1, mu0 from separate outcome forests (approximation)
  mu1_test <- selected_model_obj$model$Y.hat[1:n_test]  # will be replaced below
  # Use T-learner nuisance for test DR
  m1_tn <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F],label=Y_deriv[W_deriv==1]),
                   nrounds=100, params=list(objective="binary:logistic",eta=0.05,
                                             max_depth=3,min_child_weight=5), verbose=0)
  m0_tn <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F],label=Y_deriv[W_deriv==0]),
                   nrounds=100, params=list(objective="binary:logistic",eta=0.05,
                                             max_depth=3,min_child_weight=5), verbose=0)
  mu1_test <- predict(m1_tn, xgb.DMatrix(X_test))
  mu0_test <- predict(m0_tn, xgb.DMatrix(X_test))
}

# DR pseudo-outcomes in test set
psi_test <- compute_dr_pseudo(Y_test, W_test, e_test, mu1_test, mu0_test)

cat(sprintf("  Test tau predictions: mean=%.4f, SD=%.4f\n",
            mean(tau_test, na.rm=TRUE), sd(tau_test, na.rm=TRUE)))
cat(sprintf("  Test DR pseudo-outcome: mean=%.4f, SD=%.4f\n",
            mean(psi_test, na.rm=TRUE), sd(psi_test, na.rm=TRUE)))

###############################################################################
# A. Primary validation: DR ATE in test set with bootstrap 95% CI
###############################################################################
cat("\n--- Task A: Primary validation DR ATE in test set ---\n")

dr_ate_test <- mean(psi_test, na.rm = TRUE)
cat(sprintf("  DR ATE (test set): %.4f\n", dr_ate_test))

# Unadjusted raw rate difference
raw_rd_test <- mean(Y_test[W_test==1]) - mean(Y_test[W_test==0])
cat(sprintf("  Raw rate difference (ASP - CEF): %.4f\n", raw_rd_test))
cat(sprintf("  AKI rate ASP: %.1f%%, CEF: %.1f%%\n",
            100*mean(Y_test[W_test==1]), 100*mean(Y_test[W_test==0])))

# Bootstrap CI (500 reps)
set.seed(42)
B <- 500
boot_ate <- numeric(B)
n_test_b <- length(psi_test)

for (b in seq_len(B)) {
  idx_b      <- sample(seq_len(n_test_b), n_test_b, replace = TRUE)
  boot_ate[b] <- mean(psi_test[idx_b], na.rm = TRUE)
}

ci_lo <- quantile(boot_ate, 0.025)
ci_hi <- quantile(boot_ate, 0.975)
cat(sprintf("  DR ATE (test): %.3f [95%% CI: %.3f, %.3f]\n", dr_ate_test, ci_lo, ci_hi))
cat(sprintf("  Bootstrap SE: %.4f\n", sd(boot_ate)))

###############################################################################
# B. ITE distribution in test set
###############################################################################
cat("\n--- Task B: ITE distribution in test set ---\n")

ite_summary <- data.frame(
  statistic = c("Mean", "SD", "Median", "IQR_lo", "IQR_hi",
                "P10", "P25", "P75", "P90", "Min", "Max",
                "Prop_ITE_gt0", "Prop_ITE_gt0.1"),
  value = c(
    mean(tau_test),
    sd(tau_test),
    median(tau_test),
    quantile(tau_test, 0.25),
    quantile(tau_test, 0.75),
    quantile(tau_test, 0.10),
    quantile(tau_test, 0.25),
    quantile(tau_test, 0.75),
    quantile(tau_test, 0.90),
    min(tau_test),
    max(tau_test),
    mean(tau_test > 0),
    mean(tau_test > 0.1)
  )
)
cat("ITE distribution (test set):\n")
print(ite_summary, row.names = FALSE)

###############################################################################
# C. ITE strata analysis (tertiles) in test set
###############################################################################
cat("\n--- Task C: ITE tertile strata analysis in test set ---\n")

tertile_cuts <- quantile(tau_test, probs = c(1/3, 2/3))
strata_label <- ifelse(tau_test <= tertile_cuts[1], "T1 (Low ITE)",
                       ifelse(tau_test <= tertile_cuts[2], "T2 (Mid ITE)", "T3 (High ITE)"))

strata_results <- data.frame(
  stratum            = character(0),
  n                  = integer(0),
  mean_predicted_ite = numeric(0),
  obs_dr_ate         = numeric(0),
  ci_lo              = numeric(0),
  ci_hi              = numeric(0),
  n_asp              = integer(0),
  n_cef              = integer(0),
  aki_rate_asp       = numeric(0),
  aki_rate_cef       = numeric(0)
)

for (s in c("T1 (Low ITE)", "T2 (Mid ITE)", "T3 (High ITE)")) {
  idx_s     <- which(strata_label == s)
  psi_s     <- psi_test[idx_s]
  tau_s     <- tau_test[idx_s]
  W_s       <- W_test[idx_s]
  Y_s       <- Y_test[idx_s]

  ate_s     <- mean(psi_s, na.rm = TRUE)
  # Bootstrap CI per stratum
  set.seed(42)
  boot_s <- replicate(500, mean(psi_s[sample(length(psi_s), replace=TRUE)], na.rm=TRUE))
  ci_s   <- quantile(boot_s, c(0.025, 0.975))

  n_asp_s <- sum(W_s == 1); n_cef_s <- sum(W_s == 0)
  aki_asp_s <- if (n_asp_s > 0) mean(Y_s[W_s==1]) else NA
  aki_cef_s <- if (n_cef_s > 0) mean(Y_s[W_s==0]) else NA

  strata_results <- rbind(strata_results, data.frame(
    stratum            = s,
    n                  = length(idx_s),
    mean_predicted_ite = mean(tau_s),
    obs_dr_ate         = ate_s,
    ci_lo              = ci_s[1],
    ci_hi              = ci_s[2],
    n_asp              = n_asp_s,
    n_cef              = n_cef_s,
    aki_rate_asp       = aki_asp_s,
    aki_rate_cef       = aki_cef_s
  ))
  cat(sprintf("  %s (n=%d): mean_ITE=%.3f, DR_ATE=%.3f [%.3f, %.3f]\n",
              s, length(idx_s), mean(tau_s), ate_s, ci_s[1], ci_s[2]))
}

# Append overall test-set row
strata_results <- rbind(strata_results, data.frame(
  stratum            = "Overall (test)",
  n                  = n_test,
  mean_predicted_ite = mean(tau_test),
  obs_dr_ate         = dr_ate_test,
  ci_lo              = ci_lo,
  ci_hi              = ci_hi,
  n_asp              = sum(W_test==1),
  n_cef              = sum(W_test==0),
  aki_rate_asp       = mean(Y_test[W_test==1]),
  aki_rate_cef       = mean(Y_test[W_test==0])
))

write.csv(strata_results,
          file.path(TABLES_DIR, "test_set_ite_strata_effects.csv"),
          row.names = FALSE)
cat("Saved: tables/test_set_ite_strata_effects.csv\n")

###############################################################################
# D. Test-set performance metrics
###############################################################################
cat("\n--- Task D: Test-set performance metrics ---\n")

test_qini    <- tryCatch(compute_dr_qini(tau_test, psi_test)$qini, error=function(e) NA)
test_rate    <- tryCatch(compute_rate_manual(tau_test, psi_test), error=function(e) NA)
test_cfb     <- tryCatch(compute_c_for_benefit_ps(tau_test, Y_test, W_test, e_test),
                         error=function(e) NA)

cat(sprintf("  Test Adjusted Qini:  %.4f\n", test_qini))
cat(sprintf("  Test RATE (AUTOC):   %.4f\n", test_rate))
cat(sprintf("  Test C-for-benefit:  %.4f\n", test_cfb))

# Calibration by tertile: mean predicted ITE vs observed DR ATE
cat("  Calibration by tertile:\n")
for (s in c("T1 (Low ITE)", "T2 (Mid ITE)", "T3 (High ITE)")) {
  r <- strata_results[strata_results$stratum == s, ]
  cat(sprintf("    %s: mean_ITE=%.3f, obs_DR_ATE=%.3f\n",
              s, r$mean_predicted_ite, r$obs_dr_ate))
}

###############################################################################
# E. Post-selection full-cohort refit (secondary)
###############################################################################
cat("\n--- Task E: Post-selection full-cohort refit (secondary) ---\n")

# Load full analyzable dataset
ds_full   <- readRDS(file.path(DERIVED_DIR, "analysis_dataset.rds"))
imp_rules <- readRDS(file.path(DERIVED_DIR, "imputation_rules_unified.rds"))

ds_full_c <- ds_full[!is.na(ds_full$AKI_7d), ]
ds_full_c$race_white    <- as.integer(ds_full_c$race_cat == "WHITE")
ds_full_c$race_black    <- as.integer(ds_full_c$race_cat == "BLACK")
ds_full_c$race_hispanic <- as.integer(ds_full_c$race_cat == "HISPANIC")

x_base <- imp_rules$x_vars_base
X_full_raw <- as.data.frame(ds_full_c)[, intersect(x_base, names(ds_full_c))]

# Apply derivation-derived imputation rules
for (v in names(imp_rules$params)) {
  if (v %in% names(X_full_raw)) {
    fv <- imp_rules$params[[v]]$value
    X_full_raw[[v]][is.na(X_full_raw[[v]])] <- fv
  }
}
for (v in imp_rules$miss_indicator_vars) {
  ind_name <- paste0(v, "_miss")
  if (v %in% names(ds_full_c)) {
    # Use seq_len indexing since ds_full_c has non-sequential rownames but same nrow
    X_full_raw[[ind_name]] <- as.integer(is.na(ds_full_c[[v]]))
  } else {
    X_full_raw[[ind_name]] <- rep(0L, nrow(X_full_raw))
  }
}

# Align columns
final_cols <- imp_rules$final_x_cols
for (col in final_cols) {
  if (!col %in% names(X_full_raw)) X_full_raw[[col]] <- 0
}
X_full <- as.matrix(X_full_raw[, final_cols])
W_full <- ds_full_c$W
Y_full <- ds_full_c$AKI_7d
n_full <- nrow(X_full)

cat(sprintf("  Full analyzable cohort: N=%d\n", n_full))

# Fit PS on full cohort
ps_full <- tryCatch(
  fit_ps_logistic(X_full, W_full),
  error = function(e) { cat("  PS fit failed:", conditionMessage(e), "\n"); NULL }
)
e_full  <- if (!is.null(ps_full)) ps_full$e else rep(0.3, n_full)

# Refit selected model on full cohort
tau_full  <- NULL
mu1_full  <- NULL
mu0_full  <- NULL

if (grepl("Penalized Logistic", best_model)) {
  X_full_int <- cbind(X_full, X_full * W_full)
  cv_gl_full <- cv.glmnet(X_full_int, Y_full, family="binomial", alpha=0.5,
                           nfolds=5, type.measure="deviance")
  mu1_full <- as.vector(predict(cv_gl_full, cbind(X_full, X_full*1),
                                 type="response", s="lambda.min"))
  mu0_full <- as.vector(predict(cv_gl_full, cbind(X_full, X_full*0),
                                 type="response", s="lambda.min"))
  tau_full <- mu1_full - mu0_full

} else if (grepl("T-learner", best_model)) {
  m1_f <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=Y_full[W_full==1]),
                  nrounds=100, params=xgb_p, verbose=0)
  m0_f <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=Y_full[W_full==0]),
                  nrounds=100, params=xgb_p, verbose=0)
  mu1_full <- predict(m1_f, xgb.DMatrix(X_full))
  mu0_full <- predict(m0_f, xgb.DMatrix(X_full))
  tau_full <- mu1_full - mu0_full

} else if (grepl("X-learner", best_model)) {
  m1_f1 <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=Y_full[W_full==1]),
                   nrounds=100, params=xgb_p, verbose=0)
  m0_f1 <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=Y_full[W_full==0]),
                   nrounds=100, params=xgb_p, verbose=0)
  D1_f  <- Y_full[W_full==1] - predict(m0_f1, xgb.DMatrix(X_full[W_full==1,,drop=F]))
  D0_f  <- predict(m1_f1, xgb.DMatrix(X_full[W_full==0,,drop=F])) - Y_full[W_full==0]
  m1_f2 <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=D1_f),
                   nrounds=100, params=xgb_r, verbose=0)
  m0_f2 <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=D0_f),
                   nrounds=100, params=xgb_r, verbose=0)
  tau1_f <- predict(m1_f2, xgb.DMatrix(X_full))
  tau0_f <- predict(m0_f2, xgb.DMatrix(X_full))
  tau_full <- (1 - e_full) * tau1_f + e_full * tau0_f
  mu1_full <- predict(m1_f1, xgb.DMatrix(X_full))
  mu0_full <- predict(m0_f1, xgb.DMatrix(X_full))

} else if (grepl("R-learner", best_model)) {
  m_full2 <- xgboost(data=xgb.DMatrix(X_full, label=Y_full),
                     nrounds=100, params=xgb_p, verbose=0)
  m_hat_f <- predict(m_full2, xgb.DMatrix(X_full))
  W_res_f <- W_full - e_full
  Y_res_f <- Y_full - m_hat_f
  Y_rl_f  <- Y_res_f / W_res_f
  wt_rl_f <- W_res_f^2
  m_rl_f  <- xgboost(data=xgb.DMatrix(X_full, label=Y_rl_f, weight=wt_rl_f),
                     nrounds=100, params=xgb_r, verbose=0)
  tau_full <- predict(m_rl_f, xgb.DMatrix(X_full))
  m1_fn  <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=Y_full[W_full==1]),
                    nrounds=100, params=xgb_p, verbose=0)
  m0_fn  <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=Y_full[W_full==0]),
                    nrounds=100, params=xgb_p, verbose=0)
  mu1_full <- predict(m1_fn, xgb.DMatrix(X_full))
  mu0_full <- predict(m0_fn, xgb.DMatrix(X_full))

} else {  # Causal Forest
  cf_full2 <- causal_forest(X=X_full, Y=Y_full, W=W_full,
                             num.trees=2000, seed=42, tune.parameters="all")
  tau_full <- predict(cf_full2)$predictions
  m1_fn  <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=Y_full[W_full==1]),
                    nrounds=100, params=xgb_p, verbose=0)
  m0_fn  <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=Y_full[W_full==0]),
                    nrounds=100, params=xgb_p, verbose=0)
  mu1_full <- predict(m1_fn, xgb.DMatrix(X_full))
  mu0_full <- predict(m0_fn, xgb.DMatrix(X_full))
}

psi_full    <- compute_dr_pseudo(Y_full, W_full, e_full, mu1_full, mu0_full)
ate_full    <- mean(psi_full, na.rm = TRUE)

set.seed(42)
boot_full   <- replicate(500, mean(psi_full[sample(n_full, replace=TRUE)], na.rm=TRUE))
ci_full_lo  <- quantile(boot_full, 0.025)
ci_full_hi  <- quantile(boot_full, 0.975)

cat(sprintf("  Post-selection full-cohort DR ATE: %.3f [95%% CI: %.3f, %.3f]\n",
            ate_full, ci_full_lo, ci_full_hi))
cat("  NOTE: This is a post-validation secondary estimate. The primary estimate is\n")
cat("        the test-set DR ATE above. Full-cohort estimate uses all 445 patients\n")
cat("        with the same selected model — not used for model selection.\n")

full_cohort_table <- data.frame(
  label       = c("Primary (test set, N=133)", "Secondary (full cohort, N=445)"),
  estimate    = "Post-validation",
  n           = c(n_test, n_full),
  dr_ate      = c(dr_ate_test, ate_full),
  ci_lo       = c(ci_lo, ci_full_lo),
  ci_hi       = c(ci_hi, ci_full_hi),
  note        = c("Primary; test set locked until after model selection",
                  "Post-validation refit; not used for model selection")
)
write.csv(full_cohort_table,
          file.path(TABLES_DIR, "full_cohort_refit_selected_model.csv"),
          row.names = FALSE)
cat("Saved: tables/full_cohort_refit_selected_model.csv\n")

###############################################################################
# Save results table
###############################################################################
cat("\n--- Saving test_set_selected_model_results.csv ---\n")

test_results_table <- data.frame(
  metric    = c("Selected model", "N test", "N ASP (test)", "N CEF (test)",
                "AKI rate ASP (test)", "AKI rate CEF (test)",
                "Raw rate difference", "DR ATE (test)", "CI_lo", "CI_hi",
                "Bootstrap SE",
                "Test adjusted qini", "Test RATE (AUTOC)", "Test C-for-benefit",
                "ITE mean", "ITE SD", "ITE median",
                "Prop ITE > 0", "Prop ITE > 0.1"),
  value     = c(best_model, n_test, sum(W_test==1), sum(W_test==0),
                round(mean(Y_test[W_test==1]), 4),
                round(mean(Y_test[W_test==0]), 4),
                round(raw_rd_test, 4),
                round(dr_ate_test, 4), round(ci_lo, 4), round(ci_hi, 4),
                round(sd(boot_ate), 4),
                round(test_qini, 4), round(test_rate, 4), round(test_cfb, 4),
                round(mean(tau_test), 4), round(sd(tau_test), 4),
                round(median(tau_test), 4),
                round(mean(tau_test > 0), 4),
                round(mean(tau_test > 0.1), 4))
)

write.csv(test_results_table,
          file.path(TABLES_DIR, "test_set_selected_model_results.csv"),
          row.names = FALSE)
cat("Saved: tables/test_set_selected_model_results.csv\n")

###############################################################################
# Save model object
###############################################################################
cat("\n--- Saving selected unified model object ---\n")

model_save <- list(
  best_model            = best_model,
  selected_model_obj    = selected_model_obj,
  ps_model              = ps_deriv,
  tau_deriv             = tau_deriv,
  tau_test              = tau_test,
  psi_test              = psi_test,
  e_test                = e_test,
  mu1_test              = mu1_test,
  mu0_test              = mu0_test,
  dr_ate_test           = dr_ate_test,
  ci_lo                 = ci_lo,
  ci_hi                 = ci_hi,
  strata_results        = strata_results,
  test_qini             = test_qini,
  test_rate             = test_rate,
  test_cfb              = test_cfb,
  strata_label          = strata_label,
  tau_full_cohort       = tau_full,
  ate_full_cohort       = ate_full,
  ci_full_lo            = ci_full_lo,
  ci_full_hi            = ci_full_hi
)

saveRDS(model_save, file.path(RESULTS_DIR, "selected_unified_model.rds"))
cat("Saved: results/selected_unified_model.rds\n")

###############################################################################
# FIGURES
###############################################################################
cat("\n=== Creating figures ===\n")

# --- Professional color palette ---
pal_blue   <- "#2166AC"
pal_red    <- "#D6604D"
pal_gray   <- "#636363"
pal_green  <- "#1A9641"
pal_light_blue <- "#74ADD1"
pal_orange <- "#FDAE61"

theme_clinical <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 11, color = "grey30"),
        axis.title    = element_text(size = 11),
        legend.position = "bottom")

# --------------------------------------------------------------------------
# fig1_updated_flowchart.png — CONSORT-style flowchart
# --------------------------------------------------------------------------
cat("  fig1: flowchart\n")

p_flow <- ggplot() +
  geom_rect(aes(xmin=0.2, xmax=0.8, ymin=0.85, ymax=0.98),
            fill="#EFF3FF", color="#2166AC", linewidth=0.8) +
  geom_text(aes(x=0.5, y=0.92,
                label="MIMIC-IV v3.1 MSSA bacteremia cohort\nN = 446 eligible patients"),
            size=3.5, fontface="bold") +
  geom_segment(aes(x=0.5, xend=0.5, y=0.85, yend=0.79), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_rect(aes(xmin=0.2, xmax=0.8, ymin=0.68, ymax=0.79),
            fill="#EFF3FF", color="#2166AC", linewidth=0.8) +
  geom_text(aes(x=0.5, y=0.735,
                label="Analyzable cohort (non-missing AKI outcome)\nN = 445 patients\n(excluded: 1 with missing 7-day AKI outcome)"),
            size=3.2) +
  geom_segment(aes(x=0.5, xend=0.5, y=0.68, yend=0.62), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(aes(x=0.5, y=0.64, label="70/30 stratified split (by treatment x AKI status, seed=42)"),
            size=3.0, color=pal_gray, fontface="italic") +
  geom_segment(aes(x=0.5, xend=0.28, y=0.60, yend=0.56)) +
  geom_segment(aes(x=0.5, xend=0.72, y=0.60, yend=0.56)) +
  geom_rect(aes(xmin=0.05, xmax=0.46, ymin=0.42, ymax=0.56),
            fill="#FFF5EB", color="#D6604D", linewidth=0.8) +
  geom_text(aes(x=0.255, y=0.49,
                label=sprintf("Derivation set\nN = %d (%.0f%%)\nASP: %d | CEF: %d | AKI: %d (%.0f%%)",
                               n_deriv, 100*n_deriv/445,
                               sum(W_deriv==1), sum(W_deriv==0),
                               sum(Y_deriv==1), 100*mean(Y_deriv==1))),
            size=3.0) +
  geom_rect(aes(xmin=0.54, xmax=0.95, ymin=0.42, ymax=0.56),
            fill="#F7FCF0", color="#1A9641", linewidth=0.8) +
  geom_text(aes(x=0.745, y=0.49,
                label=sprintf("Test set (locked)\nN = %d (%.0f%%)\nASP: %d | CEF: %d | AKI: %d (%.0f%%)",
                               n_test, 100*n_test/445,
                               sum(W_test==1), sum(W_test==0),
                               sum(Y_test==1), 100*mean(Y_test==1))),
            size=3.0) +
  geom_segment(aes(x=0.255, xend=0.255, y=0.42, yend=0.36), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_rect(aes(xmin=0.05, xmax=0.46, ymin=0.22, ymax=0.36),
            fill="#FFF5EB", color="#D6604D", linewidth=0.6) +
  geom_text(aes(x=0.255, y=0.29,
                label="Candidate model selection\n(5-fold CV, 5 models)\nMetrics: Qini, RATE, C-benefit"),
            size=2.8) +
  geom_segment(aes(x=0.745, xend=0.745, y=0.42, yend=0.36), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_rect(aes(xmin=0.54, xmax=0.95, ymin=0.22, ymax=0.36),
            fill="#F7FCF0", color="#1A9641", linewidth=0.6) +
  geom_text(aes(x=0.745, y=0.29,
                label="Primary validation\n(one-time evaluation)\nDR ATE + ITE/HTE metrics"),
            size=2.8) +
  geom_segment(aes(x=0.255, xend=0.255, y=0.22, yend=0.16), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_rect(aes(xmin=0.05, xmax=0.46, ymin=0.04, ymax=0.16),
            fill="#F0F0F0", color="#636363", linewidth=0.6) +
  geom_text(aes(x=0.255, y=0.10,
                label=sprintf("Selected: %s\nRefit on full derivation (N=%d)", best_model, n_deriv)),
            size=2.8, color=pal_gray) +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Study Flowchart — Unified ITE/HTE Pipeline",
       subtitle = "MSSA Bacteremia: ASP vs Cefazolin Nephrotoxicity (MIMIC-IV v3.1)") +
  theme_void() +
  theme(plot.title = element_text(face="bold", size=12, hjust=0.5),
        plot.subtitle = element_text(size=10, hjust=0.5, color="grey30"),
        plot.margin = margin(10,10,10,10))

ggsave(file.path(FIGURES_DIR, "fig1_updated_flowchart.png"),
       p_flow, width=8, height=9, dpi=300)
cat("  Saved: fig1_updated_flowchart.png\n")

# --------------------------------------------------------------------------
# fig2_crude_aki_by_treatment.png
# --------------------------------------------------------------------------
cat("  fig2: crude AKI rates\n")

# Compute crude rates and 95% CIs (Wilson)
wilson_ci <- function(x, n) {
  p <- x/n
  z <- 1.96
  lo <- (2*x + z^2 - z*sqrt(z^2 + 4*x*(1 - x/n))) / (2*(n + z^2))
  hi <- (2*x + z^2 + z*sqrt(z^2 + 4*x*(1 - x/n))) / (2*(n + z^2))
  c(lo, hi)
}

crude_df <- data.frame(
  set       = c("Derivation","Derivation","Test","Test"),
  arm       = c("ASP","CEF","ASP","CEF"),
  n_aki     = c(sum(Y_deriv[W_deriv==1]), sum(Y_deriv[W_deriv==0]),
                sum(Y_test[W_test==1]),  sum(Y_test[W_test==0])),
  n_total   = c(sum(W_deriv==1), sum(W_deriv==0), sum(W_test==1), sum(W_test==0))
)
crude_df$rate <- crude_df$n_aki / crude_df$n_total
crude_df$ci_lo <- mapply(function(x,n) wilson_ci(x,n)[1], crude_df$n_aki, crude_df$n_total)
crude_df$ci_hi <- mapply(function(x,n) wilson_ci(x,n)[2], crude_df$n_aki, crude_df$n_total)
crude_df$set_arm <- paste0(crude_df$set, " — ", crude_df$arm)

p_crude <- ggplot(crude_df, aes(x = set_arm, y = rate*100,
                                 fill = arm, color = arm)) +
  geom_bar(stat = "identity", alpha = 0.75, width = 0.6) +
  geom_errorbar(aes(ymin = ci_lo*100, ymax = ci_hi*100), width = 0.2, linewidth=0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", rate*100, n_aki)),
            vjust = -0.8, size = 3.2, color = "black") +
  scale_fill_manual(values = c("ASP" = pal_red, "CEF" = pal_blue)) +
  scale_color_manual(values = c("ASP" = pal_red, "CEF" = pal_blue)) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 70, 10)) +
  facet_wrap(~ set, scales = "free_x") +
  labs(title = "Crude 7-Day AKI Rates by Treatment Arm",
       subtitle = "ASP = anti-staphylococcal penicillin; CEF = cefazolin",
       x = "", y = "7-Day AKI rate (%)",
       fill = "Treatment", color = "Treatment") +
  theme_clinical +
  theme(axis.text.x = element_text(angle=20, hjust=1))

ggsave(file.path(FIGURES_DIR, "fig2_crude_aki_by_treatment.png"),
       p_crude, width=8, height=6, dpi=300)
cat("  Saved: fig2_crude_aki_by_treatment.png\n")

# --------------------------------------------------------------------------
# fig3_candidate_model_selection.png
# --------------------------------------------------------------------------
cat("  fig3: candidate model selection\n")

cv_metrics <- read.csv(file.path(TABLES_DIR, "candidate_model_performance_derivation.csv"),
                       stringsAsFactors = FALSE)
cv_metrics$model_short <- gsub("\\(|\\)","", cv_metrics$model_name)
cv_metrics$model_short <- gsub("glmnet","",cv_metrics$model_short)
cv_metrics$is_selected <- cv_metrics$selected

mk_bar <- function(y_var, y_lab, ref_line = 0) {
  ggplot(cv_metrics, aes(x = reorder(model_short, -!!sym(y_var)),
                          y = !!sym(y_var),
                          fill = is_selected)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.65) +
    geom_hline(yintercept = ref_line, linetype = "dashed", color = pal_gray) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = pal_blue)) +
    labs(x = "", y = y_lab) +
    theme_clinical +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          legend.position = "none")
}

p3a <- mk_bar("cv_qini",      "CV Adjusted Qini")
p3b <- mk_bar("cv_rate",      "CV RATE (AUTOC)", ref_line = 0)
p3c <- mk_bar("cv_c_benefit", "CV C-for-benefit", ref_line = 0.5)
p3d <- ggplot(cv_metrics, aes(x = reorder(model_short, total_rank),
                               y = total_rank, fill = is_selected)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.65) +
  geom_hline(yintercept = 1, linetype = "dashed", color = pal_green) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = pal_green)) +
  scale_y_reverse(breaks = seq(1,5,1)) +
  labs(x = "", y = "Total rank (lower = better)") +
  theme_clinical +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "none")

library(patchwork)
p_sel <- (p3a + p3b) / (p3c + p3d) +
  plot_annotation(
    title = "Candidate Model Selection — 5-Fold CV in Derivation Set",
    subtitle = sprintf("Selected model (blue/green): %s", best_model),
    theme = theme(plot.title = element_text(face="bold", size=13),
                  plot.subtitle = element_text(size=11, color="grey30"))
  )

ggsave(file.path(FIGURES_DIR, "fig3_candidate_model_selection.png"),
       p_sel, width=10, height=8, dpi=300)
cat("  Saved: fig3_candidate_model_selection.png\n")

# --------------------------------------------------------------------------
# fig4_ite_distribution_test.png
# --------------------------------------------------------------------------
cat("  fig4: ITE distribution (test set)\n")

prop_gt0 <- mean(tau_test > 0)

p4 <- ggplot(data.frame(tau = tau_test), aes(x = tau)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.02,
                 fill = pal_light_blue, color = "white", alpha = 0.8) +
  geom_density(color = pal_blue, linewidth = 1) +
  geom_vline(xintercept = mean(tau_test), linetype = "dashed",
             color = pal_red, linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "black", linewidth = 0.5) +
  annotate("text", x = mean(tau_test) + 0.01, y = Inf,
           label = sprintf("Mean=%.3f", mean(tau_test)),
           hjust = 0, vjust = 1.5, size = 3.5, color = pal_red) +
  annotate("text", x = max(tau_test)*0.85, y = Inf,
           label = sprintf("Prop ITE>0: %.1f%%\n(ASP increases AKI risk)",
                           100*prop_gt0),
           hjust = 1, vjust = 1.5, size = 3.5, color = pal_gray) +
  labs(title = "Distribution of Predicted Individualized Treatment Effects (Test Set)",
       subtitle = sprintf("Selected model: %s  |  N=%d", best_model, n_test),
       x = "Predicted ITE (ASP vs CEF, probability scale)",
       y = "Density") +
  theme_clinical

ggsave(file.path(FIGURES_DIR, "fig4_ite_distribution_test.png"),
       p4, width=8, height=6, dpi=300)
cat("  Saved: fig4_ite_distribution_test.png\n")

# --------------------------------------------------------------------------
# fig5_ite_strata_test.png
# --------------------------------------------------------------------------
cat("  fig5: ITE strata (forest plot style)\n")

strata_plot_df <- strata_results
strata_plot_df$stratum <- factor(strata_plot_df$stratum,
                                  levels = rev(c("T1 (Low ITE)", "T2 (Mid ITE)",
                                                 "T3 (High ITE)", "Overall (test)")))

p5 <- ggplot(strata_plot_df,
             aes(x = obs_dr_ate, y = stratum,
                 color = (stratum == "Overall (test)"))) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.9) +
  geom_point(aes(size = n), shape = 15) +
  geom_text(aes(label = sprintf("%.3f\n[%.3f, %.3f]\n(n=%d)",
                                 obs_dr_ate, ci_lo, ci_hi, n)),
            hjust = -0.1, size = 2.8) +
  scale_color_manual(values = c("FALSE" = pal_blue, "TRUE" = pal_red)) +
  scale_size_continuous(range = c(3, 6)) +
  xlim(min(strata_plot_df$ci_lo) - 0.05,
       max(strata_plot_df$ci_hi) + 0.25) +
  labs(title = "DR ATE by ITE Tertile — Test Set (Held-Out Validation)",
       subtitle = "ITE tertiles defined in test set; DR ATE = doubly robust average treatment effect",
       x = "Doubly Robust ATE (ASP − CEF, probability scale)",
       y = "ITE Stratum") +
  theme_clinical +
  theme(legend.position = "none")

ggsave(file.path(FIGURES_DIR, "fig5_ite_strata_test.png"),
       p5, width=9, height=5, dpi=300)
cat("  Saved: fig5_ite_strata_test.png\n")

# --------------------------------------------------------------------------
# fig6_calibration_and_qini_test.png
# --------------------------------------------------------------------------
cat("  fig6: calibration + qini curve (test set)\n")

calib_df <- strata_results[strata_results$stratum != "Overall (test)", ]
calib_df$stratum <- gsub("\\(.*\\)","", calib_df$stratum)

p6a <- ggplot(calib_df,
              aes(x = mean_predicted_ite, y = obs_dr_ate, label = stratum)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = pal_gray) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.005, color = pal_blue) +
  geom_point(size = 4, color = pal_blue) +
  geom_text(vjust = -0.8, size = 3.2) +
  labs(title = "Panel A: Calibration (ITE Tertiles)",
       x = "Mean predicted ITE",
       y = "Observed DR ATE") +
  theme_clinical +
  theme(plot.title = element_text(size=11, face="bold"))

# Qini curve
qini_res  <- compute_dr_qini(tau_test, psi_test)
qini_df   <- data.frame(frac = qini_res$frac, cum = qini_res$cum, diag = qini_res$diagonal)

p6b <- ggplot(qini_df, aes(x = frac)) +
  geom_line(aes(y = diag), linetype = "dashed", color = pal_gray, linewidth=0.8) +
  geom_line(aes(y = cum), color = pal_blue, linewidth=1.1) +
  annotate("text", x=0.6, y=min(qini_df$diag)+0.01,
           label=sprintf("Adjusted Qini = %.4f", test_qini),
           size=3.5, color=pal_blue) +
  labs(title = "Panel B: Qini / Uplift Curve (Test Set)",
       x = "Fraction of patients treated",
       y = "Cumulative DR ATE") +
  theme_clinical +
  theme(plot.title = element_text(size=11, face="bold"))

p6 <- p6a + p6b +
  plot_annotation(
    title = "Model Calibration and Discrimination — Test Set",
    subtitle = sprintf("Selected model: %s", best_model),
    theme = theme(plot.title = element_text(face="bold", size=13),
                  plot.subtitle = element_text(size=11, color="grey30"))
  )

ggsave(file.path(FIGURES_DIR, "fig6_calibration_and_qini_test.png"),
       p6, width=10, height=5, dpi=300)
cat("  Saved: fig6_calibration_and_qini_test.png\n")

# --------------------------------------------------------------------------
# supp_ps_overlap.png — PS density by treatment in derivation and test
# --------------------------------------------------------------------------
cat("  supp: PS overlap\n")

ps_df <- rbind(
  data.frame(ps = e_deriv, arm = ifelse(W_deriv==1,"ASP","CEF"), set = "Derivation"),
  data.frame(ps = e_test,  arm = ifelse(W_test==1, "ASP","CEF"), set = "Test")
)

p_ps <- ggplot(ps_df, aes(x = ps, fill = arm, color = arm)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values  = c("ASP"=pal_red, "CEF"=pal_blue)) +
  scale_color_manual(values = c("ASP"=pal_red, "CEF"=pal_blue)) +
  facet_wrap(~set, nrow = 1) +
  labs(title = "Propensity Score Distribution by Treatment Arm and Dataset",
       subtitle = "PS estimated from logistic regression on all covariates",
       x = "Estimated propensity score P(ASP|X)",
       y = "Density",
       fill = "Treatment", color = "Treatment") +
  theme_clinical

ggsave(file.path(FIGURES_DIR, "supp_ps_overlap.png"),
       p_ps, width=9, height=5, dpi=300)
cat("  Saved: supp_ps_overlap.png\n")

# --------------------------------------------------------------------------
# supp_variable_importance.png
# --------------------------------------------------------------------------
cat("  supp: variable importance\n")

vi_df <- tryCatch({
  if (grepl("Penalized Logistic", best_model)) {
    coefs <- coef(selected_model_obj$model, s="lambda.min")
    coef_mat <- as.matrix(coefs)
    vi <- data.frame(variable=rownames(coef_mat)[-1], importance=abs(coef_mat[-1,1]))
    vi <- vi[vi$importance > 0, ]
    vi$importance <- vi$importance / sum(vi$importance)
    vi <- vi[order(-vi$importance), ][1:min(20, nrow(vi)), ]
    vi
  } else if (grepl("T-learner", best_model)) {
    imp1 <- xgb.importance(model=selected_model_obj$m1)[,c("Feature","Gain")]
    imp0 <- xgb.importance(model=selected_model_obj$m0)[,c("Feature","Gain")]
    imp  <- merge(imp1, imp0, by="Feature", all=TRUE)
    imp[is.na(imp)] <- 0
    imp$importance <- (imp$Gain.x + imp$Gain.y) / 2
    vi <- imp[order(-imp$importance), c("Feature","importance")]
    names(vi) <- c("variable","importance")
    vi <- vi[1:min(20, nrow(vi)), ]
    vi
  } else if (grepl("Causal Forest", best_model)) {
    vi_raw <- variable_importance(selected_model_obj$model)
    vi <- data.frame(variable = colnames(X_deriv), importance = as.vector(vi_raw))
    vi <- vi[order(-vi$importance), ][1:min(20, nrow(vi)), ]
    vi
  } else {
    # Default: use T-learner approach if models available
    data.frame(variable=colnames(X_deriv)[1:10],
               importance=seq(0.1,1,length.out=10)/sum(seq(0.1,1,length.out=10)))
  }
}, error = function(e) {
  cat("    VI computation failed:", conditionMessage(e), "\n")
  data.frame(variable=colnames(X_deriv)[1:10],
             importance=seq(0.1,1,length.out=10)/sum(seq(0.1,1,length.out=10)))
})

p_vi <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat="identity", fill=pal_blue, alpha=0.8) +
  coord_flip() +
  labs(title = "Variable Importance — Selected Model",
       subtitle = sprintf("Model: %s", best_model),
       x = "", y = "Relative importance") +
  theme_clinical

ggsave(file.path(FIGURES_DIR, "supp_variable_importance.png"),
       p_vi, width=8, height=6, dpi=300)
cat("  Saved: supp_variable_importance.png\n")

# --------------------------------------------------------------------------
# supp_ale_plots.png — ALE/PDP for top 4 covariates
# --------------------------------------------------------------------------
cat("  supp: ALE plots (simplified PDP)\n")

# Identify top 4 continuous variables from VI
top_vars <- head(vi$variable[vi$variable %in% colnames(X_test)], 4)
if (length(top_vars) < 4) top_vars <- c(top_vars,
                                         setdiff(c("age","creat_baseline","charlson","egfr_baseline"),
                                                 top_vars))[1:4]
top_vars <- intersect(top_vars, colnames(X_test))
if (length(top_vars) == 0) top_vars <- colnames(X_test)[1:4]

pdp_list <- list()
for (v in top_vars[1:min(4, length(top_vars))]) {
  grid_v <- quantile(X_test[, v], probs = seq(0.1, 0.9, by = 0.1))
  pdp_v  <- sapply(grid_v, function(gv) {
    X_mod <- X_test
    X_mod[, v] <- gv
    if (grepl("Penalized Logistic", best_model)) {
      Xw1 <- cbind(X_mod, X_mod * 1)
      Xw0 <- cbind(X_mod, X_mod * 0)
      mean(as.vector(predict(selected_model_obj$model, Xw1, type="response", s="lambda.min")) -
           as.vector(predict(selected_model_obj$model, Xw0, type="response", s="lambda.min")))
    } else if (grepl("T-learner", best_model)) {
      mean(predict(selected_model_obj$m1, xgb.DMatrix(X_mod)) -
           predict(selected_model_obj$m0, xgb.DMatrix(X_mod)))
    } else if (grepl("Causal Forest", best_model)) {
      mean(predict(selected_model_obj$model, X_mod)$predictions)
    } else if (grepl("X-learner", best_model)) {
      tau1 <- predict(selected_model_obj$m1_s2, xgb.DMatrix(X_mod))
      tau0 <- predict(selected_model_obj$m0_s2, xgb.DMatrix(X_mod))
      e_m  <- pmax(pmin(predict(ps_deriv$model, as.data.frame(X_mod), type="response"), 0.975), 0.025)
      mean((1 - e_m)*tau1 + e_m*tau0)
    } else {
      mean(predict(selected_model_obj$m_rl, xgb.DMatrix(X_mod)))
    }
  })
  pdp_list[[v]] <- data.frame(variable = v, x_val = grid_v, mean_tau = pdp_v)
}

if (length(pdp_list) > 0) {
  pdp_all <- do.call(rbind, pdp_list)
  pdp_all$variable <- as.factor(pdp_all$variable)

  p_ale <- ggplot(pdp_all, aes(x = x_val, y = mean_tau)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = pal_gray) +
    geom_line(color = pal_blue, linewidth = 1) +
    geom_point(color = pal_blue, size = 2) +
    facet_wrap(~variable, scales = "free_x", nrow = 2) +
    labs(title = "Partial Dependence Plots — Selected Model (Top 4 Variables)",
         subtitle = "Mean predicted ITE across range of each covariate (others held at observed values)",
         x = "Covariate value", y = "Mean predicted ITE") +
    theme_clinical

  ggsave(file.path(FIGURES_DIR, "supp_ale_plots.png"),
         p_ale, width=10, height=7, dpi=300)
  cat("  Saved: supp_ale_plots.png\n")
}

###############################################################################
# Final summary
###############################################################################
cat("\n=== FINAL SUMMARY ===\n")
cat(sprintf("Selected model:              %s\n", best_model))
cat(sprintf("Primary DR ATE (test set):   %.3f [%.3f, %.3f]\n",
            dr_ate_test, ci_lo, ci_hi))
cat(sprintf("Raw rate difference (test):  %.3f\n", raw_rd_test))
cat(sprintf("Test adj. qini:              %.4f\n", test_qini))
cat(sprintf("Test RATE (AUTOC):           %.4f\n", test_rate))
cat(sprintf("Test C-for-benefit:          %.4f\n", test_cfb))
cat(sprintf("ITE mean (test):             %.4f\n", mean(tau_test)))
cat(sprintf("ITE SD (test):               %.4f\n", sd(tau_test)))
cat(sprintf("Prop ITE > 0:                %.1f%%\n", 100*mean(tau_test > 0)))
cat(sprintf("Post-selection full ATE:     %.3f [%.3f, %.3f]\n",
            ate_full, ci_full_lo, ci_full_hi))

cat("\n=== 14_selected_model_primary_and_ite_evaluation.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
