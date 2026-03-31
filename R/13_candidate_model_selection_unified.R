###############################################################################
# 13_candidate_model_selection_unified.R
# Candidate causal ML model selection via 5-fold CV in derivation set
# Follows: Buell et al. JAMA 2024 (Rboost / R-learner selection framework)
#          Munroe et al. Lancet Resp Med 2025 (effect-based modelling framework)
#
# Candidate models:
#   1. Penalized logistic regression with treatment-covariate interactions (glmnet)
#   2. T-learner with XGBoost
#   3. X-learner with XGBoost
#   4. R-learner (manual Robinson decomposition with XGBoost)
#   5. Causal Forest (grf)
#
# Selection metrics (all out-of-fold, derivation only):
#   - Adjusted Qini (DR doubly-robust version)
#   - RATE / AUTOC (rank-weighted DR ATE)
#   - C-for-benefit (observational analogue via PS matching)
#
# Outputs:
#   tables/candidate_model_performance_derivation.csv  (overwrite)
#   results/derivation_model_selection_metrics.rds
###############################################################################

library(dplyr)
library(grf)
library(xgboost)
library(glmnet)

cat("=== 13_candidate_model_selection_unified.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

set.seed(42)

WORK_DIR    <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
RESULTS_DIR <- file.path(WORK_DIR, "results")
TABLES_DIR  <- file.path(WORK_DIR, "tables")

###############################################################################
# 1. Load derivation set
###############################################################################
cat("--- Step 1: Load derivation set ---\n")
deriv <- readRDS(file.path(DERIVED_DIR, "derivation_set_unified.rds"))

X <- deriv$X
W <- deriv$W
Y <- deriv$Y
n <- nrow(X)

cat(sprintf("Derivation N=%d, p=%d covariates\n", n, ncol(X)))
cat(sprintf("ASP: %d (%.1f%%), CEF: %d (%.1f%%)\n",
            sum(W==1), 100*mean(W==1), sum(W==0), 100*mean(W==0)))
cat(sprintf("AKI: %d (%.1f%%)\n", sum(Y==1), 100*mean(Y==1)))

###############################################################################
# Helper functions
###############################################################################

# DR pseudo-outcome (AIPW augmented)
compute_dr_pseudo <- function(Y, W, e_hat, mu1_hat, mu0_hat) {
  e_hat <- pmax(pmin(e_hat, 0.975), 0.025)
  psi   <- (W - e_hat) / (e_hat * (1 - e_hat)) *
            (Y - ifelse(W == 1, mu1_hat, mu0_hat)) +
            (mu1_hat - mu0_hat)
  return(psi)
}

# DR-adjusted qini coefficient
# Patients ranked ascending by tau_hat (most negative = most ASP-beneficial)
# Qini = mean(cumulative DR ATE - diagonal)
compute_dr_qini <- function(tau_hat, psi_dr) {
  ok  <- !is.na(tau_hat) & !is.na(psi_dr)
  tau_hat <- tau_hat[ok]; psi_dr <- psi_dr[ok]
  if (length(tau_hat) < 5) return(list(qini = NA))
  ord          <- order(tau_hat)
  psi_sorted   <- psi_dr[ord]
  n_obs        <- length(psi_sorted)
  cum_effect   <- cumsum(psi_sorted) / n_obs
  frac_treated <- seq_along(psi_sorted) / n_obs
  overall_mean <- mean(psi_dr)
  diagonal     <- frac_treated * overall_mean
  qini_auc     <- mean(cum_effect - diagonal)
  list(qini = qini_auc, frac = frac_treated, cum = cum_effect,
       diagonal = diagonal, overall_mean = overall_mean)
}

# RATE / AUTOC (manual rank-weighted DR ATE)
compute_rate_manual <- function(tau_hat, psi_dr) {
  ok  <- !is.na(tau_hat) & !is.na(psi_dr)
  tau_hat <- tau_hat[ok]; psi_dr <- psi_dr[ok]
  if (length(tau_hat) < 5) return(NA)
  ord          <- order(tau_hat)
  psi_sorted   <- psi_dr[ord]
  cum_mean     <- cumsum(psi_sorted) / seq_along(psi_sorted)
  overall_mean <- mean(psi_dr)
  autoc        <- mean(cum_mean - overall_mean)
  return(autoc)
}

# C-for-benefit (PS-matched pairs, observational)
compute_c_for_benefit_ps <- function(tau_hat, Y, W, e_hat) {
  treated_idx <- which(W == 1)
  control_idx <- which(W == 0)
  if (length(treated_idx) < 2 || length(control_idx) < 2) return(NA)
  e_treated   <- e_hat[treated_idx]
  e_control   <- e_hat[control_idx]

  # Greedy 1:1 PS matching without replacement
  matched_pairs <- data.frame(trt_idx = integer(0), ctrl_idx = integer(0))
  avail_ctrl    <- seq_along(control_idx)

  for (i in seq_along(treated_idx)) {
    if (length(avail_ctrl) == 0) break
    dist_ps <- abs(e_treated[i] - e_control[avail_ctrl])
    best    <- avail_ctrl[which.min(dist_ps)]
    matched_pairs <- rbind(matched_pairs,
                           data.frame(trt_idx  = treated_idx[i],
                                      ctrl_idx = control_idx[best]))
    avail_ctrl <- avail_ctrl[avail_ctrl != best]
  }

  # Caliper: 0.2 * SD of logit(PS)
  logit_ps <- log(e_hat / (1 - e_hat))
  caliper  <- 0.2 * sd(logit_ps)
  logit_diff <- abs(logit_ps[matched_pairs$trt_idx] - logit_ps[matched_pairs$ctrl_idx])
  matched_pairs <- matched_pairs[logit_diff <= caliper, ]

  if (nrow(matched_pairs) < 5) {
    cat("    C-for-benefit: fewer than 5 matched pairs after caliper — using DR rank fallback\n")
    # Fallback: rank tau_hat against psi (DR pseudo-outcome as proxy for observed benefit)
    ok <- !is.na(tau_hat) & !is.na(e_hat)
    if (sum(ok) < 5) return(NA)
    # Spearman-like: negative correlation between tau_hat and psi means those with
    # lower (more beneficial) tau have higher DR benefit
    cor_val <- cor(-tau_hat[ok], Y[ok], method = "kendall")
    return(0.5 + cor_val / 2)  # map [-1,1] to [0,1]
  }

  # Observed benefit per pair: Y_control - Y_treated (1 if treatment reduced AKI)
  pair_benefit <- Y[matched_pairs$ctrl_idx] - Y[matched_pairs$trt_idx]
  pair_tau     <- tau_hat[matched_pairs$trt_idx]

  n_pairs    <- nrow(matched_pairs)
  concordant <- 0
  total_disc <- 0

  for (i in seq_len(n_pairs - 1)) {
    for (j in (i + 1):n_pairs) {
      if (pair_benefit[i] != pair_benefit[j]) {
        total_disc <- total_disc + 1
        if ((pair_benefit[i] > pair_benefit[j]) == (pair_tau[i] < pair_tau[j])) {
          concordant <- concordant + 1
        }
      }
    }
  }
  if (total_disc == 0) return(NA)
  return(concordant / total_disc)
}

# Nuisance model (outcome + propensity, fit on training fold)
fit_nuisance_xgb <- function(X_tr, W_tr, Y_tr, X_val) {
  X_tr_w1 <- X_tr[W_tr == 1, , drop = FALSE]
  X_tr_w0 <- X_tr[W_tr == 0, , drop = FALSE]
  Y_tr_w1 <- Y_tr[W_tr == 1]
  Y_tr_w0 <- Y_tr[W_tr == 0]

  xgb_params <- list(objective = "binary:logistic", eta = 0.05,
                     max_depth = 3, min_child_weight = 5)

  fit_m1 <- xgboost(data = xgb.DMatrix(X_tr_w1, label = Y_tr_w1),
                    nrounds = 100, params = xgb_params, verbose = 0)
  fit_m0 <- xgboost(data = xgb.DMatrix(X_tr_w0, label = Y_tr_w0),
                    nrounds = 100, params = xgb_params, verbose = 0)

  mu1_val <- predict(fit_m1, xgb.DMatrix(X_val))
  mu0_val <- predict(fit_m0, xgb.DMatrix(X_val))

  # Propensity: logistic regression (more stable than XGBoost at small N)
  fit_ps <- tryCatch(
    glm(W_tr ~ ., data = as.data.frame(X_tr), family = binomial()),
    error = function(e) {
      glm(W_tr ~ ., data = as.data.frame(X_tr), family = binomial(),
          control = list(maxit = 100))
    }
  )
  e_val <- predict(fit_ps, newdata = as.data.frame(X_val), type = "response")

  list(mu1 = mu1_val, mu0 = mu0_val, e = e_val,
       mu1_model = fit_m1, mu0_model = fit_m0, ps_model = fit_ps)
}

###############################################################################
# 2. Setup stratified 5-fold CV
###############################################################################
cat("\n--- Step 2: Setup 5-fold CV (stratified by W x Y) ---\n")

K           <- 5
stratum_cv  <- paste0("W", W, "_Y", Y)
folds       <- vector("list", K)

set.seed(42)
for (s in unique(stratum_cv)) {
  idx_s         <- which(stratum_cv == s)
  idx_shuffled  <- sample(idx_s)
  chunk_sizes   <- diff(round(seq(0, length(idx_shuffled), length.out = K + 1)))
  start <- 1
  for (k in seq_len(K)) {
    end_k      <- start + chunk_sizes[k] - 1
    folds[[k]] <- c(folds[[k]], idx_shuffled[start:end_k])
    start      <- end_k + 1
  }
}

cat(sprintf("Fold sizes: %s\n", paste(sapply(folds, length), collapse = ", ")))

###############################################################################
# 3. Candidate models definition
###############################################################################
cat("\n--- Step 3: Candidate models ---\n")

model_names <- c(
  "Penalized Logistic (glmnet)",
  "T-learner XGBoost",
  "X-learner XGBoost",
  "R-learner (manual)",
  "Causal Forest (grf)"
)

for (i in seq_along(model_names)) cat(sprintf("  %d. %s\n", i, model_names[i]))

###############################################################################
# 4. 5-fold CV loop
###############################################################################
cat("\n--- Step 4: 5-fold CV (all models) ---\n")

oof_tau     <- matrix(NA, nrow = n, ncol = length(model_names))
colnames(oof_tau) <- model_names

oof_e_hat   <- rep(NA, n)
oof_mu1_hat <- rep(NA, n)
oof_mu0_hat <- rep(NA, n)

for (k in seq_len(K)) {
  cat(sprintf("\n  === Fold %d/%d ===\n", k, K))
  val_idx <- folds[[k]]
  tr_idx  <- setdiff(seq_len(n), val_idx)

  X_tr <- X[tr_idx,  , drop = FALSE]
  X_vl <- X[val_idx, , drop = FALSE]
  W_tr <- W[tr_idx];  W_vl <- W[val_idx]
  Y_tr <- Y[tr_idx];  Y_vl <- Y[val_idx]

  # Fit nuisance models
  nuis <- tryCatch(
    fit_nuisance_xgb(X_tr, W_tr, Y_tr, X_vl),
    error = function(e) {
      cat(sprintf("    Nuisance failed: %s\n", conditionMessage(e))); NULL
    }
  )
  if (is.null(nuis)) {
    cat("    Skipping fold — nuisance model failure\n"); next
  }
  oof_e_hat[val_idx]   <- nuis$e
  oof_mu1_hat[val_idx] <- nuis$mu1
  oof_mu0_hat[val_idx] <- nuis$mu0

  # --------------------------------------------------------------------------
  # Model 1: Penalized logistic regression with treatment-covariate interactions
  # --------------------------------------------------------------------------
  tryCatch({
    # Include W and all X*W interaction features
    X_tr_int    <- cbind(X_tr, X_tr * W_tr)
    X_vl_w1     <- cbind(X_vl, X_vl * 1)
    X_vl_w0     <- cbind(X_vl, X_vl * 0)

    cv_gl       <- cv.glmnet(X_tr_int, Y_tr, family = "binomial", alpha = 0.5,
                              nfolds = 3, type.measure = "deviance")
    mu1_gl      <- as.vector(predict(cv_gl, X_vl_w1, type = "response", s = "lambda.min"))
    mu0_gl      <- as.vector(predict(cv_gl, X_vl_w0, type = "response", s = "lambda.min"))
    oof_tau[val_idx, 1] <- mu1_gl - mu0_gl
    cat("    1. Penalized Logistic: OK\n")
  }, error = function(e) cat(sprintf("    1. Penalized Logistic FAILED: %s\n", conditionMessage(e))))

  # --------------------------------------------------------------------------
  # Model 2: T-learner (separate XGBoost per arm)
  # --------------------------------------------------------------------------
  tryCatch({
    xgb_p <- list(objective = "binary:logistic", eta = 0.05, max_depth = 3,
                  min_child_weight = 5)
    m1_t  <- xgboost(data = xgb.DMatrix(X_tr[W_tr==1,,drop=F], label=Y_tr[W_tr==1]),
                     nrounds=100, params=xgb_p, verbose=0)
    m0_t  <- xgboost(data = xgb.DMatrix(X_tr[W_tr==0,,drop=F], label=Y_tr[W_tr==0]),
                     nrounds=100, params=xgb_p, verbose=0)
    oof_tau[val_idx, 2] <- predict(m1_t, xgb.DMatrix(X_vl)) -
                           predict(m0_t, xgb.DMatrix(X_vl))
    cat("    2. T-learner XGBoost: OK\n")
  }, error = function(e) cat(sprintf("    2. T-learner FAILED: %s\n", conditionMessage(e))))

  # --------------------------------------------------------------------------
  # Model 3: X-learner (two-stage XGBoost)
  # --------------------------------------------------------------------------
  tryCatch({
    xgb_p  <- list(objective = "binary:logistic", eta = 0.05, max_depth = 3,
                   min_child_weight = 5)
    xgb_r  <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 3,
                   min_child_weight = 5)
    # Stage 1: mu1, mu0
    m1_x1  <- xgboost(data = xgb.DMatrix(X_tr[W_tr==1,,drop=F], label=Y_tr[W_tr==1]),
                      nrounds=100, params=xgb_p, verbose=0)
    m0_x1  <- xgboost(data = xgb.DMatrix(X_tr[W_tr==0,,drop=F], label=Y_tr[W_tr==0]),
                      nrounds=100, params=xgb_p, verbose=0)
    # Stage 2 pseudo-outcomes
    D1     <- Y_tr[W_tr==1] - predict(m0_x1, xgb.DMatrix(X_tr[W_tr==1,,drop=F]))
    D0     <- predict(m1_x1, xgb.DMatrix(X_tr[W_tr==0,,drop=F])) - Y_tr[W_tr==0]
    m1_x2  <- xgboost(data = xgb.DMatrix(X_tr[W_tr==1,,drop=F], label=D1),
                      nrounds=100, params=xgb_r, verbose=0)
    m0_x2  <- xgboost(data = xgb.DMatrix(X_tr[W_tr==0,,drop=F], label=D0),
                      nrounds=100, params=xgb_r, verbose=0)
    tau1_v <- predict(m1_x2, xgb.DMatrix(X_vl))
    tau0_v <- predict(m0_x2, xgb.DMatrix(X_vl))
    e_clip  <- pmax(pmin(nuis$e, 0.975), 0.025)
    oof_tau[val_idx, 3] <- (1 - e_clip) * tau1_v + e_clip * tau0_v
    cat("    3. X-learner XGBoost: OK\n")
  }, error = function(e) cat(sprintf("    3. X-learner FAILED: %s\n", conditionMessage(e))))

  # --------------------------------------------------------------------------
  # Model 4: R-learner (manual Robinson decomposition)
  # --------------------------------------------------------------------------
  tryCatch({
    xgb_p   <- list(objective = "binary:logistic", eta = 0.05, max_depth = 3,
                    min_child_weight = 5)
    xgb_r   <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 3,
                    min_child_weight = 5)
    # m(x) = E[Y|X] (ignoring W) — use XGBoost
    m_full  <- xgboost(data = xgb.DMatrix(X_tr, label=Y_tr),
                       nrounds=100, params=xgb_p, verbose=0)
    m_hat   <- predict(m_full, xgb.DMatrix(X_tr))
    # e(x) = P(W=1|X) — logistic regression
    e_hat_tr <- tryCatch(
      predict(glm(W_tr ~ ., data=as.data.frame(X_tr), family=binomial()),
              newdata=as.data.frame(X_tr), type="response"),
      error = function(e2) pmax(pmin(nuis$e, 0.975), 0.025)
    )
    e_hat_tr <- pmax(pmin(e_hat_tr, 0.975), 0.025)
    # R-learner pseudo-outcome
    Y_res    <- Y_tr - m_hat
    W_res    <- W_tr - e_hat_tr
    Y_rl     <- Y_res / W_res
    wt_rl    <- W_res^2
    m_rl     <- xgboost(data = xgb.DMatrix(X_tr, label=Y_rl, weight=wt_rl),
                        nrounds=100, params=xgb_r, verbose=0)
    oof_tau[val_idx, 4] <- predict(m_rl, xgb.DMatrix(X_vl))
    cat("    4. R-learner (manual): OK\n")
  }, error = function(e) cat(sprintf("    4. R-learner FAILED: %s\n", conditionMessage(e))))

  # --------------------------------------------------------------------------
  # Model 5: Causal Forest (grf)
  # --------------------------------------------------------------------------
  tryCatch({
    cf_k <- causal_forest(X = X_tr, Y = Y_tr, W = W_tr,
                          num.trees = 1000, seed = 42 + k,
                          tune.parameters = "all")
    oof_tau[val_idx, 5] <- predict(cf_k, newdata = X_vl)$predictions
    cat("    5. Causal Forest (grf): OK\n")
  }, error = function(e) cat(sprintf("    5. Causal Forest FAILED: %s\n", conditionMessage(e))))

  cat(sprintf("  Fold %d complete\n", k))
}

###############################################################################
# 5. Compute DR pseudo-outcomes for full derivation set
###############################################################################
cat("\n--- Step 5: Compute DR pseudo-outcomes ---\n")

cat(sprintf("Missing nuisance — e_hat: %d, mu1_hat: %d, mu0_hat: %d\n",
            sum(is.na(oof_e_hat)), sum(is.na(oof_mu1_hat)), sum(is.na(oof_mu0_hat))))

psi_dr  <- compute_dr_pseudo(Y, W, oof_e_hat, oof_mu1_hat, oof_mu0_hat)
dr_ate  <- mean(psi_dr, na.rm = TRUE)
cat(sprintf("Derivation DR ATE (cross-fitted, all models): %.4f\n", dr_ate))
cat(sprintf("DR pseudo-outcome: mean=%.4f, SD=%.4f\n",
            mean(psi_dr, na.rm=TRUE), sd(psi_dr, na.rm=TRUE)))

###############################################################################
# 6. Evaluate each model
###############################################################################
cat("\n--- Step 6: Evaluate candidate models ---\n")

metrics <- data.frame(
  model_name    = model_names,
  cv_qini       = NA_real_,
  cv_rate       = NA_real_,
  cv_c_benefit  = NA_real_,
  n_complete    = NA_integer_,
  stringsAsFactors = FALSE
)

for (m in seq_along(model_names)) {
  tau_m       <- oof_tau[, m]
  n_ok        <- sum(!is.na(tau_m))
  metrics$n_complete[m] <- n_ok

  if (n_ok < 20) {
    cat(sprintf("  %s: too few predictions (%d), skipping\n", model_names[m], n_ok))
    next
  }

  ok      <- !is.na(tau_m) & !is.na(psi_dr) & !is.na(oof_e_hat)
  tau_ok  <- tau_m[ok]
  psi_ok  <- psi_dr[ok]
  e_ok    <- oof_e_hat[ok]
  Y_ok    <- Y[ok]
  W_ok    <- W[ok]

  # Adjusted Qini
  q_res <- tryCatch(compute_dr_qini(tau_ok, psi_ok), error = function(e) list(qini = NA))
  metrics$cv_qini[m] <- q_res$qini

  # RATE / AUTOC
  metrics$cv_rate[m] <- tryCatch(compute_rate_manual(tau_ok, psi_ok), error = function(e) NA)

  # C-for-benefit
  metrics$cv_c_benefit[m] <- tryCatch(
    compute_c_for_benefit_ps(tau_ok, Y_ok, W_ok, e_ok),
    error = function(e) NA
  )

  cat(sprintf("  %s: qini=%.4f, RATE=%.4f, C-benefit=%.4f (n=%d)\n",
              model_names[m],
              ifelse(is.na(metrics$cv_qini[m]),      NA, metrics$cv_qini[m]),
              ifelse(is.na(metrics$cv_rate[m]),       NA, metrics$cv_rate[m]),
              ifelse(is.na(metrics$cv_c_benefit[m]),  NA, metrics$cv_c_benefit[m]),
              n_ok))
}

###############################################################################
# 7. Rank and select best model
###############################################################################
cat("\n--- Step 7: Rank and select best model ---\n")

# Higher is better for all three metrics
metrics$rank_qini   <- rank(-metrics$cv_qini,      na.last = "keep", ties.method = "min")
metrics$rank_rate   <- rank(-metrics$cv_rate,       na.last = "keep", ties.method = "min")
metrics$rank_c_benefit <- rank(-metrics$cv_c_benefit, na.last = "keep", ties.method = "min")

# Tiebreaker: qini (weight 3) > RATE (weight 2) > C-benefit (weight 1)
metrics$total_rank_score <- 3 * metrics$rank_qini +
                             2 * metrics$rank_rate +
                             1 * metrics$rank_c_benefit
metrics$total_rank  <- rank(metrics$total_rank_score, na.last = "keep", ties.method = "min")

cat("\nModel performance table (sorted by total rank):\n")
print(metrics[order(metrics$total_rank), ], row.names = FALSE)

best_idx   <- which.min(metrics$total_rank)
best_model <- model_names[best_idx]
cat(sprintf("\n*** SELECTED MODEL: %s (total_rank=%d) ***\n",
            best_model, metrics$total_rank[best_idx]))
cat(sprintf("  CV qini=%.4f, RATE=%.4f, C-benefit=%.4f\n",
            metrics$cv_qini[best_idx],
            metrics$cv_rate[best_idx],
            metrics$cv_c_benefit[best_idx]))

metrics$selected <- (metrics$model_name == best_model)

###############################################################################
# 8. Save outputs
###############################################################################
cat("\n--- Step 8: Save outputs ---\n")

write.csv(metrics,
          file.path(TABLES_DIR, "candidate_model_performance_derivation.csv"),
          row.names = FALSE)
cat("Saved: tables/candidate_model_performance_derivation.csv\n")

results_obj <- list(
  oof_tau         = oof_tau,
  oof_e_hat       = oof_e_hat,
  oof_mu1_hat     = oof_mu1_hat,
  oof_mu0_hat     = oof_mu0_hat,
  psi_dr          = psi_dr,
  dr_ate          = dr_ate,
  metrics         = metrics,
  best_model      = best_model,
  best_model_idx  = best_idx,
  folds           = folds,
  K               = K
)

saveRDS(results_obj,
        file.path(RESULTS_DIR, "derivation_model_selection_metrics.rds"))
cat("Saved: results/derivation_model_selection_metrics.rds\n")

cat("\n=== FINAL SUMMARY ===\n")
cat(sprintf("Best model:        %s\n", best_model))
cat(sprintf("CV Adjusted Qini:  %.4f\n", metrics$cv_qini[best_idx]))
cat(sprintf("CV RATE (AUTOC):   %.4f\n", metrics$cv_rate[best_idx]))
cat(sprintf("CV C-for-benefit:  %.4f\n", metrics$cv_c_benefit[best_idx]))
cat(sprintf("Derivation DR ATE: %.4f\n", dr_ate))

cat("\n=== 13_candidate_model_selection_unified.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
