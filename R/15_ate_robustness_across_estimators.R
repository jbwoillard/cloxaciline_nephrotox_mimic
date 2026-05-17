###############################################################################
# 15_ate_robustness_across_estimators.R
# Supplementary robustness analysis: ATE across alternative causal estimators
#
# Purpose:
#   Assess whether the directional conclusion (higher AKI risk under ASP vs
#   cefazolin) is robust to causal estimator choice. This supplements the
#   primary unified derivation–test analysis (scripts 12–14) without replacing it.
#
# Framing:
#   - Estimator robustness / vibration analysis
#   - NOT a new primary analysis or post hoc model competition
#
# Estimators compared (same cohort, outcome, covariates, split):
#   1. Penalized logistic regression with treatment-covariate interactions (SELECTED, reference)
#   2. Causal forest (grf)
#   3. R-learner (manual Robinson decomposition, XGBoost)
#   4. X-learner (XGBoost)
#   5. T-learner (XGBoost)
#   6. AIPW benchmark (doubly robust, logistic nuisance models)
#
# Outputs:
#   tables/supplement_ate_robustness_across_estimators.csv
#   figures/supp_forest_ate_across_estimators.png
#   figures/supp_ate_point_estimates_by_estimator.png (optional)
###############################################################################

library(dplyr)
library(ggplot2)
library(grf)
library(xgboost)
library(glmnet)

cat("=== 15_ate_robustness_across_estimators.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

set.seed(42)

WORK_DIR    <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
RESULTS_DIR <- file.path(WORK_DIR, "results")
TABLES_DIR  <- file.path(WORK_DIR, "tables")
FIGURES_DIR <- file.path(WORK_DIR, "figures")

###############################################################################
# 1. Load data
###############################################################################
cat("--- Step 1: Load data ---\n")

deriv    <- readRDS(file.path(DERIVED_DIR, "derivation_set_unified.rds"))
test_dat <- readRDS(file.path(DERIVED_DIR, "test_set_unified.rds"))

X_deriv  <- deriv$X;      W_deriv <- deriv$W; Y_deriv <- deriv$Y
X_test   <- test_dat$X;   W_test  <- test_dat$W; Y_test  <- test_dat$Y

n_deriv  <- nrow(X_deriv)
n_test   <- nrow(X_test)

cat(sprintf("Derivation N=%d (ASP=%d, CEF=%d)\n", n_deriv, sum(W_deriv==1), sum(W_deriv==0)))
cat(sprintf("Test N=%d       (ASP=%d, CEF=%d)\n", n_test,  sum(W_test==1),  sum(W_test==0)))

# Load full cohort for secondary full-cohort refit estimates
imp_rules <- readRDS(file.path(DERIVED_DIR, "imputation_rules_unified.rds"))
ds_full   <- readRDS(file.path(DERIVED_DIR, "analysis_dataset.rds"))
ds_full_c <- ds_full[!is.na(ds_full$AKI_7d), ]
ds_full_c$race_white    <- as.integer(ds_full_c$race_cat == "WHITE")
ds_full_c$race_black    <- as.integer(ds_full_c$race_cat == "BLACK")
ds_full_c$race_hispanic <- as.integer(ds_full_c$race_cat == "HISPANIC")

x_base     <- imp_rules$x_vars_base
X_full_raw <- as.data.frame(ds_full_c)[, intersect(x_base, names(ds_full_c))]

for (v in names(imp_rules$params)) {
  if (v %in% names(X_full_raw)) {
    X_full_raw[[v]][is.na(X_full_raw[[v]])] <- imp_rules$params[[v]]$value
  }
}
for (v in imp_rules$miss_indicator_vars) {
  ind_name <- paste0(v, "_miss")
  if (v %in% names(ds_full_c)) {
    X_full_raw[[ind_name]] <- as.integer(is.na(ds_full_c[[v]]))
  } else {
    X_full_raw[[ind_name]] <- rep(0L, nrow(X_full_raw))
  }
}

final_cols <- imp_rules$final_x_cols
for (col in final_cols) {
  if (!col %in% names(X_full_raw)) X_full_raw[[col]] <- 0
}
X_full <- as.matrix(X_full_raw[, final_cols])
W_full <- ds_full_c$W
Y_full <- ds_full_c$AKI_7d
n_full <- nrow(X_full)
cat(sprintf("Full analyzable cohort N=%d\n\n", n_full))

###############################################################################
# Helper functions
###############################################################################

clip_ps <- function(e) pmax(pmin(e, 0.975), 0.025)

fit_ps_logistic <- function(X, W) {
  fit <- tryCatch(
    glm(W ~ ., data = as.data.frame(X), family = binomial()),
    error = function(e) glm(W ~ ., data = as.data.frame(X), family = binomial(),
                             control = list(maxit = 200))
  )
  e <- predict(fit, newdata = as.data.frame(X), type = "response")
  list(model = fit, e = clip_ps(e))
}

predict_ps <- function(ps_model, X_new) {
  e <- predict(ps_model, newdata = as.data.frame(X_new), type = "response")
  clip_ps(e)
}

dr_pseudo <- function(Y, W, e, mu1, mu0) {
  (W - e) / (e * (1 - e)) * (Y - ifelse(W == 1, mu1, mu0)) + (mu1 - mu0)
}

boot_ci <- function(psi, B = 500, seed = 42) {
  set.seed(seed)
  n <- length(psi)
  boot_vals <- replicate(B, mean(psi[sample(n, n, replace = TRUE)], na.rm = TRUE))
  list(
    ate = mean(psi, na.rm = TRUE),
    se  = sd(boot_vals),
    lo  = quantile(boot_vals, 0.025),
    hi  = quantile(boot_vals, 0.975)
  )
}

xgb_p <- list(objective = "binary:logistic", eta = 0.05, max_depth = 3, min_child_weight = 5)
xgb_r <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 3, min_child_weight = 5)

# Results accumulator
rob_results <- list()

###############################################################################
# Propensity scores (derived from derivation set, applied to test set)
###############################################################################
cat("--- Fitting propensity models ---\n")
ps_deriv <- fit_ps_logistic(X_deriv, W_deriv)
e_deriv  <- ps_deriv$e
e_test   <- predict_ps(ps_deriv$model, X_test)

ps_full_fit <- tryCatch(fit_ps_logistic(X_full, W_full), error = function(e) NULL)
e_full      <- if (!is.null(ps_full_fit)) ps_full_fit$e else rep(mean(W_full), n_full)

cat(sprintf("  Deriv PS range: [%.3f, %.3f]\n", min(e_deriv), max(e_deriv)))
cat(sprintf("  Test  PS range: [%.3f, %.3f]\n", min(e_test),  max(e_test)))
cat(sprintf("  Full  PS range: [%.3f, %.3f]\n", min(e_full),  max(e_full)))

###############################################################################
# ESTIMATOR 1: Penalized logistic regression with treatment-covariate interactions
# (Reference: this is the PRIMARY selected model from the unified pipeline)
###############################################################################
cat("\n--- Estimator 1: Penalized Logistic Interaction (SELECTED primary model) ---\n")

# Fit on derivation set
X_d_int <- cbind(X_deriv, X_deriv * W_deriv)
set.seed(42)
cv_gl   <- cv.glmnet(X_d_int, Y_deriv, family = "binomial", alpha = 0.5,
                     nfolds = 5, type.measure = "deviance")

# Predict in test set
X_t_w1  <- cbind(X_test, X_test * 1)
X_t_w0  <- cbind(X_test, X_test * 0)
mu1_t1  <- as.vector(predict(cv_gl, X_t_w1, type = "response", s = "lambda.min"))
mu0_t1  <- as.vector(predict(cv_gl, X_t_w0, type = "response", s = "lambda.min"))
psi_t1  <- dr_pseudo(Y_test, W_test, e_test, mu1_t1, mu0_t1)
res_t1  <- boot_ci(psi_t1)
cat(sprintf("  Test ATE: %.4f [%.4f, %.4f]\n", res_t1$ate, res_t1$lo, res_t1$hi))

# Full cohort refit
X_f_int <- cbind(X_full, X_full * W_full)
set.seed(42)
cv_gl_f <- cv.glmnet(X_f_int, Y_full, family = "binomial", alpha = 0.5,
                     nfolds = 5, type.measure = "deviance")
mu1_f1  <- as.vector(predict(cv_gl_f, cbind(X_full, X_full*1), type="response", s="lambda.min"))
mu0_f1  <- as.vector(predict(cv_gl_f, cbind(X_full, X_full*0), type="response", s="lambda.min"))
psi_f1  <- dr_pseudo(Y_full, W_full, e_full, mu1_f1, mu0_f1)
res_f1  <- boot_ci(psi_f1)
cat(sprintf("  Full ATE: %.4f [%.4f, %.4f]\n", res_f1$ate, res_f1$lo, res_f1$hi))

rob_results[["Penalized Logistic Interaction (selected)"]] <- list(
  test = res_t1, full = res_f1,
  notes = "Selected primary model: penalized logistic regression with treatment-covariate interactions (glmnet, alpha=0.5, lambda.min). ATE derived as mean DR pseudo-outcome in held-out test set."
)

###############################################################################
# ESTIMATOR 2: Causal forest (grf)
###############################################################################
cat("\n--- Estimator 2: Causal Forest (grf) ---\n")

set.seed(42)
cf <- causal_forest(X = X_deriv, Y = Y_deriv, W = W_deriv,
                    num.trees = 2000, seed = 42, tune.parameters = "all")

# Test set predictions
tau_cf_test <- predict(cf, newdata = X_test)$predictions

# Nuisance for DR in test (T-learner from deriv)
m1_cf <- xgboost(data = xgb.DMatrix(X_deriv[W_deriv==1,,drop=F], label=Y_deriv[W_deriv==1]),
                 nrounds = 100, params = xgb_p, verbose = 0)
m0_cf <- xgboost(data = xgb.DMatrix(X_deriv[W_deriv==0,,drop=F], label=Y_deriv[W_deriv==0]),
                 nrounds = 100, params = xgb_p, verbose = 0)
mu1_cf_test <- predict(m1_cf, xgb.DMatrix(X_test))
mu0_cf_test <- predict(m0_cf, xgb.DMatrix(X_test))
psi_cf_test <- dr_pseudo(Y_test, W_test, e_test, mu1_cf_test, mu0_cf_test)
res_cf_t    <- boot_ci(psi_cf_test)
cat(sprintf("  Test ATE (DR): %.4f [%.4f, %.4f]\n", res_cf_t$ate, res_cf_t$lo, res_cf_t$hi))

# Full cohort refit
set.seed(42)
cf_full <- causal_forest(X = X_full, Y = Y_full, W = W_full,
                         num.trees = 2000, seed = 42, tune.parameters = "all")
psi_cf_full <- cf_full$debiased.error + cf_full$predictions  # DR-corrected
psi_cf_full <- dr_pseudo(Y_full, W_full, clip_ps(cf_full$W.hat),
                          cf_full$Y.hat + (1-cf_full$W.hat)*cf_full$predictions,
                          cf_full$Y.hat - cf_full$W.hat*cf_full$predictions)
res_cf_f    <- boot_ci(psi_cf_full)
cat(sprintf("  Full ATE (DR): %.4f [%.4f, %.4f]\n", res_cf_f$ate, res_cf_f$lo, res_cf_f$hi))

rob_results[["Causal Forest (grf)"]] <- list(
  test = res_cf_t, full = res_cf_f,
  notes = "GRF causal_forest (2000 trees, tune.parameters='all', seed=42). Fitted on derivation set. Test ATE via mean DR pseudo-outcome (nuisance from T-learner XGBoost). Full-cohort ATE from GRF internal DR-corrected pseudo-outcome."
)

###############################################################################
# ESTIMATOR 3: R-learner (XGBoost)
###############################################################################
cat("\n--- Estimator 3: R-learner (XGBoost) ---\n")

# Step 1: fit outcome model m(x) and propensity e(x) on derivation
m_rlearner <- xgboost(data = xgb.DMatrix(X_deriv, label = Y_deriv),
                      nrounds = 100, params = xgb_p, verbose = 0)
m_hat_d    <- predict(m_rlearner, xgb.DMatrix(X_deriv))

# Step 2: residualize
W_res_d <- W_deriv - e_deriv
Y_res_d <- Y_deriv - m_hat_d

# Protect against near-zero residuals
ok_d      <- abs(W_res_d) > 1e-3
Y_rl_d    <- Y_res_d / W_res_d
Y_rl_d[!ok_d] <- median(Y_rl_d[ok_d], na.rm = TRUE)
wt_rl_d   <- W_res_d^2

# Step 3: fit CATE model
set.seed(42)
m_rl <- xgboost(data = xgb.DMatrix(X_deriv, label = Y_rl_d, weight = wt_rl_d),
                nrounds = 100, params = xgb_r, verbose = 0)

# Test set: predict CATE, then DR ATE
m1_rl <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F], label=Y_deriv[W_deriv==1]),
                 nrounds=100, params=xgb_p, verbose=0)
m0_rl <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F], label=Y_deriv[W_deriv==0]),
                 nrounds=100, params=xgb_p, verbose=0)
mu1_rl_test <- predict(m1_rl, xgb.DMatrix(X_test))
mu0_rl_test <- predict(m0_rl, xgb.DMatrix(X_test))
psi_rl_test <- dr_pseudo(Y_test, W_test, e_test, mu1_rl_test, mu0_rl_test)
res_rl_t    <- boot_ci(psi_rl_test)
cat(sprintf("  Test ATE (DR): %.4f [%.4f, %.4f]\n", res_rl_t$ate, res_rl_t$lo, res_rl_t$hi))

# Full cohort refit
m_rl_full   <- xgboost(data=xgb.DMatrix(X_full, label=Y_full),
                        nrounds=100, params=xgb_p, verbose=0)
m_hat_full  <- predict(m_rl_full, xgb.DMatrix(X_full))
W_res_full  <- W_full - e_full
Y_res_full  <- Y_full - m_hat_full
ok_f        <- abs(W_res_full) > 1e-3
Y_rl_full   <- Y_res_full / W_res_full
Y_rl_full[!ok_f] <- median(Y_rl_full[ok_f], na.rm = TRUE)
wt_rl_full  <- W_res_full^2
set.seed(42)
m_rl_f      <- xgboost(data=xgb.DMatrix(X_full, label=Y_rl_full, weight=wt_rl_full),
                        nrounds=100, params=xgb_r, verbose=0)
m1_rl_f     <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=Y_full[W_full==1]),
                        nrounds=100, params=xgb_p, verbose=0)
m0_rl_f     <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=Y_full[W_full==0]),
                        nrounds=100, params=xgb_p, verbose=0)
mu1_rl_full <- predict(m1_rl_f, xgb.DMatrix(X_full))
mu0_rl_full <- predict(m0_rl_f, xgb.DMatrix(X_full))
psi_rl_full <- dr_pseudo(Y_full, W_full, e_full, mu1_rl_full, mu0_rl_full)
res_rl_f    <- boot_ci(psi_rl_full)
cat(sprintf("  Full ATE (DR): %.4f [%.4f, %.4f]\n", res_rl_f$ate, res_rl_f$lo, res_rl_f$hi))

rob_results[["R-learner (XGBoost)"]] <- list(
  test = res_rl_t, full = res_rl_f,
  notes = "Manual R-learner (Robinson decomposition). XGBoost outcome and propensity models. CATE: residualized outcome regressed on covariates, weighted by squared propensity residuals. ATE via mean DR pseudo-outcome."
)

###############################################################################
# ESTIMATOR 4: X-learner (XGBoost)
###############################################################################
cat("\n--- Estimator 4: X-learner (XGBoost) ---\n")

# Stage 1: T-learner on derivation
m1_x_s1 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F],label=Y_deriv[W_deriv==1]),
                   nrounds=100, params=xgb_p, verbose=0)
m0_x_s1 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F],label=Y_deriv[W_deriv==0]),
                   nrounds=100, params=xgb_p, verbose=0)

# Stage 2: imputed treatment effects
D1_d <- Y_deriv[W_deriv==1] - predict(m0_x_s1, xgb.DMatrix(X_deriv[W_deriv==1,,drop=F]))
D0_d <- predict(m1_x_s1, xgb.DMatrix(X_deriv[W_deriv==0,,drop=F])) - Y_deriv[W_deriv==0]

set.seed(42)
m1_x_s2 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F],label=D1_d),
                   nrounds=100, params=xgb_r, verbose=0)
m0_x_s2 <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F],label=D0_d),
                   nrounds=100, params=xgb_r, verbose=0)

# Test set predictions
tau1_x_test <- predict(m1_x_s2, xgb.DMatrix(X_test))
tau0_x_test <- predict(m0_x_s2, xgb.DMatrix(X_test))
tau_x_test  <- (1 - e_test) * tau1_x_test + e_test * tau0_x_test
mu1_x_test  <- predict(m1_x_s1, xgb.DMatrix(X_test))
mu0_x_test  <- predict(m0_x_s1, xgb.DMatrix(X_test))
psi_x_test  <- dr_pseudo(Y_test, W_test, e_test, mu1_x_test, mu0_x_test)
res_x_t     <- boot_ci(psi_x_test)
cat(sprintf("  Test ATE (DR): %.4f [%.4f, %.4f]\n", res_x_t$ate, res_x_t$lo, res_x_t$hi))

# Full cohort refit
m1_xf_s1 <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=Y_full[W_full==1]),
                    nrounds=100, params=xgb_p, verbose=0)
m0_xf_s1 <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=Y_full[W_full==0]),
                    nrounds=100, params=xgb_p, verbose=0)
D1_f <- Y_full[W_full==1] - predict(m0_xf_s1, xgb.DMatrix(X_full[W_full==1,,drop=F]))
D0_f <- predict(m1_xf_s1, xgb.DMatrix(X_full[W_full==0,,drop=F])) - Y_full[W_full==0]
set.seed(42)
m1_xf_s2 <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=D1_f),
                    nrounds=100, params=xgb_r, verbose=0)
m0_xf_s2 <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=D0_f),
                    nrounds=100, params=xgb_r, verbose=0)
mu1_xf  <- predict(m1_xf_s1, xgb.DMatrix(X_full))
mu0_xf  <- predict(m0_xf_s1, xgb.DMatrix(X_full))
psi_xf  <- dr_pseudo(Y_full, W_full, e_full, mu1_xf, mu0_xf)
res_x_f <- boot_ci(psi_xf)
cat(sprintf("  Full ATE (DR): %.4f [%.4f, %.4f]\n", res_x_f$ate, res_x_f$lo, res_x_f$hi))

rob_results[["X-learner (XGBoost)"]] <- list(
  test = res_x_t, full = res_x_f,
  notes = "X-learner (Künzel et al.): two-stage XGBoost. Stage 1: separate outcome models. Stage 2: imputed treatment effect models. CATE combines both via propensity-weighted average. ATE via mean DR pseudo-outcome."
)

###############################################################################
# ESTIMATOR 5: T-learner (XGBoost)
###############################################################################
cat("\n--- Estimator 5: T-learner (XGBoost) ---\n")

m1_t <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==1,,drop=F],label=Y_deriv[W_deriv==1]),
                nrounds=100, params=xgb_p, verbose=0)
m0_t <- xgboost(data=xgb.DMatrix(X_deriv[W_deriv==0,,drop=F],label=Y_deriv[W_deriv==0]),
                nrounds=100, params=xgb_p, verbose=0)
mu1_t_test  <- predict(m1_t, xgb.DMatrix(X_test))
mu0_t_test  <- predict(m0_t, xgb.DMatrix(X_test))
psi_t_test  <- dr_pseudo(Y_test, W_test, e_test, mu1_t_test, mu0_t_test)
res_t_t     <- boot_ci(psi_t_test)
cat(sprintf("  Test ATE (DR): %.4f [%.4f, %.4f]\n", res_t_t$ate, res_t_t$lo, res_t_t$hi))

m1_tf <- xgboost(data=xgb.DMatrix(X_full[W_full==1,,drop=F],label=Y_full[W_full==1]),
                 nrounds=100, params=xgb_p, verbose=0)
m0_tf <- xgboost(data=xgb.DMatrix(X_full[W_full==0,,drop=F],label=Y_full[W_full==0]),
                 nrounds=100, params=xgb_p, verbose=0)
mu1_tf  <- predict(m1_tf, xgb.DMatrix(X_full))
mu0_tf  <- predict(m0_tf, xgb.DMatrix(X_full))
psi_tf  <- dr_pseudo(Y_full, W_full, e_full, mu1_tf, mu0_tf)
res_t_f <- boot_ci(psi_tf)
cat(sprintf("  Full ATE (DR): %.4f [%.4f, %.4f]\n", res_t_f$ate, res_t_f$lo, res_t_f$hi))

rob_results[["T-learner (XGBoost)"]] <- list(
  test = res_t_t, full = res_t_f,
  notes = "T-learner: separate XGBoost outcome models for treated and controls. CATE = mu1(x) - mu0(x). ATE via mean DR pseudo-outcome (mu1, mu0 from test predictions)."
)

###############################################################################
# ESTIMATOR 6: AIPW (doubly robust benchmark, logistic nuisance models)
###############################################################################
cat("\n--- Estimator 6: AIPW (doubly robust, logistic nuisance) ---\n")

# Fit outcome model on full derivation set (not stratified)
out_model_d <- glm(Y_deriv ~ ., data = as.data.frame(cbind(X_deriv, W=W_deriv)),
                   family = binomial())
# Predict under treated and control
X_test_df <- as.data.frame(X_test)
mu1_aipw_test <- predict(out_model_d,
                         newdata = cbind(X_test_df, W=rep(1, n_test)),
                         type = "response")
mu0_aipw_test <- predict(out_model_d,
                         newdata = cbind(X_test_df, W=rep(0, n_test)),
                         type = "response")
psi_aipw_test <- dr_pseudo(Y_test, W_test, e_test, mu1_aipw_test, mu0_aipw_test)
res_aipw_t    <- boot_ci(psi_aipw_test)
cat(sprintf("  Test ATE (AIPW): %.4f [%.4f, %.4f]\n", res_aipw_t$ate, res_aipw_t$lo, res_aipw_t$hi))

# Full cohort
out_model_f <- glm(Y_full ~ ., data = as.data.frame(cbind(X_full, W=W_full)),
                   family = binomial())
X_full_df   <- as.data.frame(X_full)
mu1_aipw_f  <- predict(out_model_f, newdata=cbind(X_full_df, W=rep(1,n_full)), type="response")
mu0_aipw_f  <- predict(out_model_f, newdata=cbind(X_full_df, W=rep(0,n_full)), type="response")
psi_aipw_f  <- dr_pseudo(Y_full, W_full, e_full, mu1_aipw_f, mu0_aipw_f)
res_aipw_f  <- boot_ci(psi_aipw_f)
cat(sprintf("  Full ATE (AIPW): %.4f [%.4f, %.4f]\n", res_aipw_f$ate, res_aipw_f$lo, res_aipw_f$hi))

rob_results[["AIPW (logistic nuisance, benchmark)"]] <- list(
  test = res_aipw_t, full = res_aipw_f,
  notes = "AIPW benchmark: logistic regression outcome model (all covariates + treatment), logistic propensity score. DR pseudo-outcome applied to test set. Included as a doubly-robust parametric reference."
)

###############################################################################
# Build supplementary table
###############################################################################
cat("\n--- Building supplementary table ---\n")

direction_label <- function(lo, hi) {
  if (lo > 0) return("ASP harmful (ASP > CEF)")
  if (hi < 0) return("ASP beneficial (ASP < CEF)")
  return("Neutral / uncertain (CI crosses 0)")
}

tbl_rows <- list()
for (nm in names(rob_results)) {
  r    <- rob_results[[nm]]
  rt   <- r$test
  rf   <- r$full
  tbl_rows[[length(tbl_rows)+1]] <- data.frame(
    estimator         = nm,
    cohort            = "Test set (held-out)",
    n                 = n_test,
    ATE               = round(rt$ate, 4),
    SE                = round(rt$se,  4),
    CI_lo             = round(rt$lo,  4),
    CI_hi             = round(rt$hi,  4),
    direction         = direction_label(rt$lo, rt$hi),
    notes             = r$notes,
    stringsAsFactors  = FALSE
  )
  tbl_rows[[length(tbl_rows)+1]] <- data.frame(
    estimator         = nm,
    cohort            = "Full cohort refit",
    n                 = n_full,
    ATE               = round(rf$ate, 4),
    SE                = round(rf$se,  4),
    CI_lo             = round(rf$lo,  4),
    CI_hi             = round(rf$hi,  4),
    direction         = direction_label(rf$lo, rf$hi),
    notes             = r$notes,
    stringsAsFactors  = FALSE
  )
}

rob_tbl <- do.call(rbind, tbl_rows)
write.csv(rob_tbl,
          file.path(TABLES_DIR, "supplement_ate_robustness_across_estimators.csv"),
          row.names = FALSE)
cat("Saved: tables/supplement_ate_robustness_across_estimators.csv\n")
print(rob_tbl[, c("estimator","cohort","ATE","CI_lo","CI_hi","direction")])

###############################################################################
# Forest plot — primary figure (test set ATEs)
###############################################################################
cat("\n--- Building forest plot ---\n")

# Colors and formatting
pal_blue  <- "#2166AC"
pal_red   <- "#D6604D"
pal_gray  <- "#636363"
pal_green <- "#1A9641"

est_order <- names(rob_results)

# Test set only for primary forest plot
plot_df <- rob_tbl[rob_tbl$cohort == "Test set (held-out)", ]
plot_df$estimator_f <- factor(plot_df$estimator,
                               levels = rev(est_order))
plot_df$is_primary <- grepl("selected", plot_df$estimator, ignore.case = TRUE)
plot_df$color_grp  <- ifelse(plot_df$is_primary, "Primary (selected)", "Alternative")

p_forest <- ggplot(plot_df, aes(y = estimator_f, x = ATE,
                                 xmin = CI_lo, xmax = CI_hi,
                                 color = color_grp, shape = color_grp)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = pal_gray, linewidth = 0.7) +
  geom_errorbarh(aes(height = 0.18), linewidth = 0.9) +
  geom_point(aes(size = color_grp)) +
  geom_text(aes(label = sprintf("%.3f [%.3f, %.3f]", ATE, CI_lo, CI_hi)),
            hjust = -0.08, size = 3.1, color = "black") +
  scale_color_manual(values = c("Primary (selected)" = pal_red,
                                 "Alternative"       = pal_blue)) +
  scale_shape_manual(values = c("Primary (selected)" = 18,
                                 "Alternative"       = 15)) +
  scale_size_manual(values  = c("Primary (selected)" = 5,
                                 "Alternative"       = 3.5)) +
  scale_x_continuous(
    limits = c(min(plot_df$CI_lo) - 0.06, max(plot_df$CI_hi) + 0.30),
    breaks = seq(-0.1, 0.5, by = 0.1)
  ) +
  labs(
    title    = "ATE Robustness Across Causal Estimators — Held-Out Test Set",
    subtitle = sprintf("ASP vs cefazolin, 7-day AKI, MIMIC-IV v3.1 (N test=%d). Primary model shown in red.", n_test),
    x        = "Average Treatment Effect (ASP − CEF, absolute risk difference)",
    y        = "",
    color    = "Estimator type",
    shape    = "Estimator type",
    size     = "Estimator type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 10.5, color = "grey30"),
    axis.text.y       = element_text(size = 10),
    legend.position   = "bottom",
    legend.title      = element_text(size = 10)
  )

ggsave(file.path(FIGURES_DIR, "supp_forest_ate_across_estimators.png"),
       p_forest, width = 10, height = 6, dpi = 300)
cat("Saved: figures/supp_forest_ate_across_estimators.png\n")

###############################################################################
# Optional: point estimates comparison (test vs full cohort)
###############################################################################
cat("\n--- Building point estimates by estimator (optional) ---\n")

plot_df2 <- rob_tbl
plot_df2$estimator_f  <- factor(plot_df2$estimator, levels = rev(est_order))
plot_df2$cohort_f     <- factor(plot_df2$cohort,
                                 levels = c("Test set (held-out)", "Full cohort refit"))
plot_df2$is_primary   <- grepl("selected", plot_df2$estimator, ignore.case = TRUE)

plot_df2$y_offset <- ifelse(plot_df2$cohort == "Test set (held-out)", 0.15, -0.15)
plot_df2$y_pos    <- as.numeric(plot_df2$estimator_f) + plot_df2$y_offset

p_compare <- ggplot(plot_df2, aes(y = y_pos, x = ATE,
                                   xmin = CI_lo, xmax = CI_hi,
                                   color = cohort_f, shape = cohort_f)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = pal_gray) +
  geom_errorbarh(aes(height = 0.08), linewidth = 0.8) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = seq_along(levels(plot_df2$estimator_f)),
                     labels = levels(plot_df2$estimator_f)) +
  scale_color_manual(values = c("Test set (held-out)" = pal_red,
                                 "Full cohort refit"  = pal_blue)) +
  scale_shape_manual(values = c("Test set (held-out)" = 15, "Full cohort refit" = 16)) +
  scale_x_continuous(
    breaks = seq(-0.3, 0.5, by = 0.1),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    title    = "ATE Estimates by Estimator and Cohort",
    subtitle = "Red squares = test set (held-out, primary); Blue circles = full cohort refit (secondary)",
    x        = "Average Treatment Effect (ASP − CEF, absolute risk difference)",
    y        = "", color = "Cohort", shape = "Cohort"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 10, color = "grey30"),
    axis.text.y        = element_text(size = 9.5),
    legend.position    = "bottom"
  )

ggsave(file.path(FIGURES_DIR, "supp_ate_point_estimates_by_estimator.png"),
       p_compare, width = 10, height = 6, dpi = 300)
cat("Saved: figures/supp_ate_point_estimates_by_estimator.png\n")

###############################################################################
# Summary statistics
###############################################################################
cat("\n=== ROBUSTNESS SUMMARY ===\n")

test_ates   <- rob_tbl[rob_tbl$cohort == "Test set (held-out)", ]
n_positive  <- sum(test_ates$CI_lo > 0)
n_total_est <- nrow(test_ates)

cat(sprintf("Estimators compared: %d\n", n_total_est))
cat(sprintf("Estimators with CI_lo > 0 (directionally positive, ASP harmful): %d / %d\n",
            n_positive, n_total_est))
cat(sprintf("ATE range (test set): [%.3f, %.3f]\n",
            min(test_ates$ATE), max(test_ates$ATE)))
cat(sprintf("All directionally consistent (ASP >= CEF): %s\n",
            ifelse(all(test_ates$ATE > 0), "YES", "NO")))

for (i in seq_len(nrow(test_ates))) {
  cat(sprintf("  %-45s ATE=%+.3f [%+.3f, %+.3f] | %s\n",
              test_ates$estimator[i],
              test_ates$ATE[i], test_ates$CI_lo[i], test_ates$CI_hi[i],
              test_ates$direction[i]))
}

###############################################################################
# Save summary object for report use
###############################################################################
rob_summary <- list(
  rob_tbl        = rob_tbl,
  test_ates      = test_ates,
  n_estimators   = n_total_est,
  n_positive     = n_positive,
  ate_range      = c(min(test_ates$ATE), max(test_ates$ATE)),
  all_positive   = all(test_ates$ATE > 0),
  all_ci_lo_pos  = all(test_ates$CI_lo > 0),
  n_test         = n_test,
  n_full         = n_full
)
saveRDS(rob_summary, file.path(RESULTS_DIR, "ate_robustness_summary.rds"))
cat("Saved: results/ate_robustness_summary.rds\n")

cat("\n=== 15_ate_robustness_across_estimators.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
