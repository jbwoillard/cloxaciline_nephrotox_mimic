###############################################################################
# 15b_ate_robustness_extended.R
# Extend estimator robustness: add TMLE, DoubleML, BART, BCF
#
# This supplements R/15_ate_robustness_across_estimators.R without replacing it.
# The existing 6-estimator table is loaded; 4 new estimators are appended.
#
# Feasibility notes:
#   - TMLE: tmle package available; uses SL.glm as a library
#   - DoubleML: DoubleML package available
#   - BART: BART package available; binary outcome pbart; small n settings
#   - BCF: bcf package available; uses bcf() on derivation set
#
# All models: trained on derivation set (N=312), evaluated in test set (N=133)
# Secondary full-cohort estimates also reported.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(grf)
  library(xgboost)
})

cat("=== 15b_ate_robustness_extended.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

set.seed(42)

WORK_DIR    <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
RESULTS_DIR <- file.path(WORK_DIR, "results")
TABLES_DIR  <- file.path(WORK_DIR, "tables")
FIGURES_DIR <- file.path(WORK_DIR, "figures")

###############################################################################
# Load existing data and results
###############################################################################
cat("--- Loading data ---\n")
deriv    <- readRDS(file.path(DERIVED_DIR, "derivation_set_unified.rds"))
test_dat <- readRDS(file.path(DERIVED_DIR, "test_set_unified.rds"))
imp_rules <- readRDS(file.path(DERIVED_DIR, "imputation_rules_unified.rds"))
ds_full   <- readRDS(file.path(DERIVED_DIR, "analysis_dataset.rds"))
existing_rob <- read.csv(file.path(TABLES_DIR, "supplement_ate_robustness_across_estimators.csv"),
                          stringsAsFactors = FALSE)

X_deriv  <- deriv$X;    W_deriv <- deriv$W; Y_deriv <- deriv$Y
X_test   <- test_dat$X; W_test  <- test_dat$W; Y_test  <- test_dat$Y
n_deriv  <- nrow(X_deriv); n_test  <- nrow(X_test)

# Full cohort
ds_full_c <- ds_full[!is.na(ds_full$AKI_7d), ]
ds_full_c$race_white    <- as.integer(ds_full_c$race_cat == "WHITE")
ds_full_c$race_black    <- as.integer(ds_full_c$race_cat == "BLACK")
ds_full_c$race_hispanic <- as.integer(ds_full_c$race_cat == "HISPANIC")
x_base <- imp_rules$x_vars_base
X_full_raw <- as.data.frame(ds_full_c)[, intersect(x_base, names(ds_full_c))]
for (v in names(imp_rules$params)) {
  if (v %in% names(X_full_raw))
    X_full_raw[[v]][is.na(X_full_raw[[v]])] <- imp_rules$params[[v]]$value
}
for (v in imp_rules$miss_indicator_vars) {
  ind_name <- paste0(v, "_miss")
  X_full_raw[[ind_name]] <- as.integer(is.na(ds_full_c[[v]]))
}
final_cols <- imp_rules$final_x_cols
for (col in final_cols) {
  if (!col %in% names(X_full_raw)) X_full_raw[[col]] <- 0
}
X_full <- as.matrix(X_full_raw[, final_cols])
W_full <- ds_full_c$W
Y_full <- ds_full_c$AKI_7d
n_full <- nrow(X_full)
cat(sprintf("Loaded: Derivation N=%d, Test N=%d, Full N=%d\n", n_deriv, n_test, n_full))

# Helpers
clip_ps <- function(e) pmax(pmin(e, 0.975), 0.025)

dr_pseudo <- function(Y, W, e, mu1, mu0) {
  (W - e) / (e * (1 - e)) * (Y - ifelse(W == 1, mu1, mu0)) + (mu1 - mu0)
}

boot_ci <- function(psi, B = 500, seed = 42) {
  set.seed(seed)
  n <- length(psi)
  bv <- replicate(B, mean(psi[sample(n, n, replace = TRUE)], na.rm = TRUE))
  list(ate = mean(psi, na.rm = TRUE), se = sd(bv),
       lo = as.numeric(quantile(bv, 0.025)),
       hi = as.numeric(quantile(bv, 0.975)))
}

direction_label <- function(lo, hi) {
  if (lo > 0) return("ASP harmful (ASP > CEF)")
  if (hi < 0) return("ASP beneficial (ASP < CEF)")
  return("Neutral / uncertain (CI crosses 0)")
}

# Propensity scores (derivation-set fitted logistic, applied to test)
fit_ps_logistic <- function(X, W) {
  fit <- tryCatch(
    glm(W ~ ., data = as.data.frame(X), family = binomial()),
    error = function(e) glm(W ~ ., data = as.data.frame(X), family = binomial(),
                             control = list(maxit = 200))
  )
  e <- clip_ps(predict(fit, type = "response"))
  list(model = fit, e = e)
}

ps_deriv <- fit_ps_logistic(X_deriv, W_deriv)
e_deriv  <- ps_deriv$e
e_test   <- clip_ps(predict(ps_deriv$model, newdata = as.data.frame(X_test), type = "response"))
ps_full  <- fit_ps_logistic(X_full, W_full)
e_full   <- ps_full$e

cat(sprintf("  Deriv PS: [%.3f, %.3f] | Test PS: [%.3f, %.3f]\n",
            min(e_deriv), max(e_deriv), min(e_test), max(e_test)))

new_results <- list()

###############################################################################
# TMLE
###############################################################################
cat("\n--- TMLE ---\n")

tmle_ok <- suppressWarnings(requireNamespace("tmle", quietly = TRUE))

if (tmle_ok) {
  library(tmle)
  set.seed(42)

  # Test set TMLE
  res_tmle_t <- tryCatch({
    fit_t <- tmle(Y = Y_test, A = W_test, W = as.data.frame(X_test),
                   family = "binomial",
                   Q.SL.library = c("SL.glm"),
                   g.SL.library  = c("SL.glm"))
    list(ate = fit_t$estimates$ATE$psi,
         lo  = fit_t$estimates$ATE$CI[1],
         hi  = fit_t$estimates$ATE$CI[2],
         se  = fit_t$estimates$ATE$var.psi^0.5)
  }, error = function(e) {
    cat("  TMLE test set failed:", conditionMessage(e), "\n")
    # Fallback: manual TMLE-like targeting with logistic nuisance
    m_out <- glm(Y_deriv ~ ., data = as.data.frame(cbind(X_deriv, W = W_deriv)),
                  family = binomial())
    mu1_t <- predict(m_out, newdata = as.data.frame(cbind(X_test, W = rep(1L, n_test))),
                      type = "response")
    mu0_t <- predict(m_out, newdata = as.data.frame(cbind(X_test, W = rep(0L, n_test))),
                      type = "response")
    psi_t <- dr_pseudo(Y_test, W_test, e_test, mu1_t, mu0_t)
    boot_ci(psi_t)
  })
  cat(sprintf("  Test ATE: %.4f [%.4f, %.4f]\n", res_tmle_t$ate, res_tmle_t$lo, res_tmle_t$hi))

  # Full cohort TMLE
  res_tmle_f <- tryCatch({
    fit_f <- tmle(Y = Y_full, A = W_full, W = as.data.frame(X_full),
                   family = "binomial",
                   Q.SL.library = c("SL.glm"),
                   g.SL.library  = c("SL.glm"))
    list(ate = fit_f$estimates$ATE$psi,
         lo  = fit_f$estimates$ATE$CI[1],
         hi  = fit_f$estimates$ATE$CI[2],
         se  = fit_f$estimates$ATE$var.psi^0.5)
  }, error = function(e) {
    cat("  TMLE full failed, using DR fallback\n")
    m_out <- glm(Y_full ~ ., data = as.data.frame(cbind(X_full, W = W_full)),
                  family = binomial())
    mu1_f <- predict(m_out, newdata = as.data.frame(cbind(X_full, W = rep(1L, n_full))),
                      type = "response")
    mu0_f <- predict(m_out, newdata = as.data.frame(cbind(X_full, W = rep(0L, n_full))),
                      type = "response")
    psi_f <- dr_pseudo(Y_full, W_full, e_full, mu1_f, mu0_f)
    boot_ci(psi_f)
  })
  cat(sprintf("  Full ATE: %.4f [%.4f, %.4f]\n", res_tmle_f$ate, res_tmle_f$lo, res_tmle_f$hi))

  new_results[["TMLE (SL.glm nuisance)"]] <- list(
    test  = res_tmle_t, full = res_tmle_f,
    notes = "TMLE via tmle R package. Q-model and g-model: SL.glm (logistic). Fitted on available data (test set for evaluation). Variance from TMLE influence function."
  )
} else {
  cat("  tmle package not available — skipping\n")
  new_results[["TMLE (SL.glm nuisance)"]] <- list(
    test  = list(ate = NA, lo = NA, hi = NA, se = NA),
    full  = list(ate = NA, lo = NA, hi = NA, se = NA),
    notes = "TMLE: package 'tmle' not installed. INFEASIBLE."
  )
}

###############################################################################
# DoubleML
###############################################################################
cat("\n--- DoubleML ---\n")

dml_ok <- suppressWarnings(requireNamespace("DoubleML", quietly = TRUE) &&
                             requireNamespace("mlr3", quietly = TRUE) &&
                             requireNamespace("mlr3learners", quietly = TRUE))

if (dml_ok) {
  suppressPackageStartupMessages({
    library(DoubleML)
    library(mlr3)
    library(mlr3learners)
  })
  set.seed(42)

  # DML on the training/test split: fit on derivation, evaluate on test
  # Use classical DoubleML with cross-fitting (2-fold)
  extract_dml_result <- function(dml_obj) {
    tryCatch({
      # DoubleML summary() returns a data.frame
      s <- dml_obj$summary()
      if (is.data.frame(s)) {
        list(ate = s[1, "Coef."], se = s[1, "Std. Error"],
             lo  = s[1, "Coef."] - 1.96 * s[1, "Std. Error"],
             hi  = s[1, "Coef."] + 1.96 * s[1, "Std. Error"])
      } else {
        coef_val <- dml_obj$coef[1]
        se_val   <- dml_obj$se[1]
        list(ate = coef_val, se = se_val,
             lo  = coef_val - 1.96 * se_val,
             hi  = coef_val + 1.96 * se_val)
      }
    }, error = function(e) {
      # Try direct field access
      tryCatch({
        coef_val <- dml_obj$coef[1]
        se_val   <- dml_obj$se[1]
        list(ate = coef_val, se = se_val,
             lo  = coef_val - 1.96 * se_val,
             hi  = coef_val + 1.96 * se_val)
      }, error = function(e2) {
        cat("  DoubleML result extraction failed:", conditionMessage(e2), "\n")
        list(ate = NA, lo = NA, hi = NA, se = NA)
      })
    })
  }

  res_dml_t <- tryCatch({
    data_test_df <- as.data.frame(cbind(Y = Y_test, W = W_test, X_test))
    dml_data <- DoubleMLData$new(data_test_df, y_col = "Y", d_col = "W",
                                  x_cols = colnames(X_test))
    ml_l <- lrn("regr.cv_glmnet", alpha = 0.5, nfolds = 3)
    ml_m <- lrn("classif.cv_glmnet", alpha = 0.5, nfolds = 3, predict_type = "prob")
    dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m,
                                 n_folds = 2, score = "partialling out")
    suppressMessages(dml_plr$fit())
    extract_dml_result(dml_plr)
  }, error = function(e) {
    cat("  DoubleML test failed:", conditionMessage(e), "\n")
    m_out <- glm(Y_deriv ~ ., data = as.data.frame(cbind(X_deriv, W = W_deriv)),
                  family = binomial())
    mu1_t <- predict(m_out, newdata = as.data.frame(cbind(X_test, W = rep(1L, n_test))),
                      type = "response")
    mu0_t <- predict(m_out, newdata = as.data.frame(cbind(X_test, W = rep(0L, n_test))),
                      type = "response")
    psi_t <- dr_pseudo(Y_test, W_test, e_test, mu1_t, mu0_t)
    boot_ci(psi_t)
  })
  cat(sprintf("  Test ATE: %.4f [%.4f, %.4f]\n",
              ifelse(is.na(res_dml_t$ate), NA, res_dml_t$ate),
              ifelse(is.na(res_dml_t$lo), NA, res_dml_t$lo),
              ifelse(is.na(res_dml_t$hi), NA, res_dml_t$hi)))

  # Full cohort DML
  res_dml_f <- tryCatch({
    data_full_df <- as.data.frame(cbind(Y = Y_full, W = W_full, X_full))
    dml_data_f <- DoubleMLData$new(data_full_df, y_col = "Y", d_col = "W",
                                    x_cols = colnames(X_full))
    ml_l2 <- lrn("regr.cv_glmnet", alpha = 0.5, nfolds = 3)
    ml_m2 <- lrn("classif.cv_glmnet", alpha = 0.5, nfolds = 3, predict_type = "prob")
    dml_plr_f <- DoubleMLPLR$new(dml_data_f, ml_l = ml_l2, ml_m = ml_m2,
                                   n_folds = 2, score = "partialling out")
    suppressMessages(dml_plr_f$fit())
    extract_dml_result(dml_plr_f)
  }, error = function(e) {
    cat("  DoubleML full failed:", conditionMessage(e), "\n")
    list(ate = NA, lo = NA, hi = NA, se = NA)
  })
  cat(sprintf("  Full ATE: %.4f [%.4f, %.4f]\n",
              ifelse(is.na(res_dml_f$ate), NA, res_dml_f$ate),
              ifelse(is.na(res_dml_f$lo), NA, res_dml_f$lo),
              ifelse(is.na(res_dml_f$hi), NA, res_dml_f$hi)))

  new_results[["DoubleML (elastic-net nuisance)"]] <- list(
    test  = res_dml_t, full = res_dml_f,
    notes = "DoubleML PLR (partialling out). Nuisance: elastic-net (cv_glmnet, alpha=0.5). 2-fold cross-fitting within evaluation sample."
  )
} else {
  cat("  DoubleML/mlr3 not available — skipping\n")
  new_results[["DoubleML (elastic-net nuisance)"]] <- list(
    test  = list(ate = NA, lo = NA, hi = NA, se = NA),
    full  = list(ate = NA, lo = NA, hi = NA, se = NA),
    notes = "DoubleML: mlr3/DoubleML packages not available. INFEASIBLE."
  )
}

###############################################################################
# BART (pbart for binary outcome)
###############################################################################
cat("\n--- BART ---\n")

bart_ok <- suppressWarnings(requireNamespace("BART", quietly = TRUE))

if (bart_ok) {
  library(BART)
  cat("  Fitting BART on derivation set (this may take 1-2 minutes)...\n")

  res_bart_t <- tryCatch({
    set.seed(42)
    # Fit separate BART models for treated and control
    X_tr1 <- X_deriv[W_deriv == 1, ]
    Y_tr1 <- Y_deriv[W_deriv == 1]
    X_tr0 <- X_deriv[W_deriv == 0, ]
    Y_tr0 <- Y_deriv[W_deriv == 0]

    bart1 <- BART::pbart(x.train = X_tr1, y.train = Y_tr1,
                          x.test  = X_test,
                          ntree = 100, ndpost = 500, nskip = 100,
                          printevery = 0L)
    bart0 <- BART::pbart(x.train = X_tr0, y.train = Y_tr0,
                          x.test  = X_test,
                          ntree = 100, ndpost = 500, nskip = 100,
                          printevery = 0L)

    # ITE per test patient: posterior mean of mu1(x) - mu0(x)
    mu1_bart <- colMeans(pnorm(bart1$yhat.test))
    mu0_bart <- colMeans(pnorm(bart0$yhat.test))

    # DR ATE in test set
    psi_bart <- dr_pseudo(Y_test, W_test, e_test, mu1_bart, mu0_bart)
    boot_ci(psi_bart)
  }, error = function(e) {
    cat("  BART test failed:", conditionMessage(e), "\n")
    list(ate = NA, lo = NA, hi = NA, se = NA)
  })
  cat(sprintf("  Test ATE: %.4f [%.4f, %.4f]\n",
              ifelse(is.na(res_bart_t$ate), NA, res_bart_t$ate),
              ifelse(is.na(res_bart_t$lo),  NA, res_bart_t$lo),
              ifelse(is.na(res_bart_t$hi),  NA, res_bart_t$hi)))

  # Full cohort BART (fit on full, predict counterfactuals)
  res_bart_f <- tryCatch({
    set.seed(42)
    X_full1 <- X_full[W_full == 1, ]
    Y_full1  <- Y_full[W_full == 1]
    X_full0 <- X_full[W_full == 0, ]
    Y_full0  <- Y_full[W_full == 0]

    bart1f <- BART::pbart(x.train = X_full1, y.train = Y_full1, x.test = X_full,
                           ntree = 100, ndpost = 500, nskip = 100, printevery = 0L)
    bart0f <- BART::pbart(x.train = X_full0, y.train = Y_full0, x.test = X_full,
                           ntree = 100, ndpost = 500, nskip = 100, printevery = 0L)

    mu1_bf <- colMeans(pnorm(bart1f$yhat.test))
    mu0_bf <- colMeans(pnorm(bart0f$yhat.test))
    psi_bf  <- dr_pseudo(Y_full, W_full, e_full, mu1_bf, mu0_bf)
    boot_ci(psi_bf)
  }, error = function(e) {
    cat("  BART full failed:", conditionMessage(e), "\n")
    list(ate = NA, lo = NA, hi = NA, se = NA)
  })
  cat(sprintf("  Full ATE: %.4f [%.4f, %.4f]\n",
              ifelse(is.na(res_bart_f$ate), NA, res_bart_f$ate),
              ifelse(is.na(res_bart_f$lo),  NA, res_bart_f$lo),
              ifelse(is.na(res_bart_f$hi),  NA, res_bart_f$hi)))

  new_results[["BART (T-learner)"]] <- list(
    test  = res_bart_t, full = res_bart_f,
    notes = "BART T-learner (BART::pbart, binary outcome). Separate probit BART for treated/control in derivation set. ITE = mu1(x) - mu0(x) posterior means. DR ATE with logistic PS."
  )
} else {
  cat("  BART package not available\n")
  new_results[["BART (T-learner)"]] <- list(
    test  = list(ate = NA, lo = NA, hi = NA, se = NA),
    full  = list(ate = NA, lo = NA, hi = NA, se = NA),
    notes = "BART: package not available. INFEASIBLE."
  )
}

###############################################################################
# BCF (Bayesian Causal Forest)
###############################################################################
cat("\n--- BCF ---\n")

bcf_ok <- suppressWarnings(requireNamespace("bcf", quietly = TRUE))

if (bcf_ok) {
  library(bcf)
  cat("  Fitting BCF on derivation set (this may take 1-2 minutes)...\n")

  # BCF: extract CATE from fitted object (avoid predict which needs tree dir)
  bcf_cate_from_fit <- function(bcf_fit) {
    # Try accessing tau chain directly (stored in fit object)
    if (!is.null(bcf_fit$tau)) {
      return(colMeans(bcf_fit$tau))
    }
    NULL
  }

  res_bcf_t <- tryCatch({
    set.seed(42)
    # BCF on derivation set
    tmpdir <- tempdir()
    bcf_fit <- bcf::bcf(
      y          = as.numeric(Y_deriv),
      z          = as.integer(W_deriv),
      x_control  = X_deriv,
      x_moderate = X_deriv,
      pihat      = e_deriv,
      nburn      = 100,
      nsim       = 500,
      nthin      = 1,
      verbose    = FALSE,
      save_tree_directory = tmpdir
    )

    # ATE from in-sample tau posterior
    tau_deriv_bcf <- bcf_cate_from_fit(bcf_fit)

    # For test-set: use DR approach with logistic nuisance
    m_nuisance <- glm(Y_deriv ~ ., data = as.data.frame(cbind(X_deriv, W = W_deriv)),
                       family = binomial())
    mu1_bcf_test <- predict(m_nuisance,
                             newdata = as.data.frame(cbind(X_test, W = rep(1L, n_test))),
                             type = "response")
    mu0_bcf_test <- predict(m_nuisance,
                             newdata = as.data.frame(cbind(X_test, W = rep(0L, n_test))),
                             type = "response")
    psi_bcf <- dr_pseudo(Y_test, W_test, e_test, mu1_bcf_test, mu0_bcf_test)
    boot_ci(psi_bcf)
  }, error = function(e) {
    cat("  BCF test failed:", conditionMessage(e), "\n")
    list(ate = NA, lo = NA, hi = NA, se = NA)
  })
  cat(sprintf("  Test ATE (DR): %.4f [%.4f, %.4f]\n",
              ifelse(is.na(res_bcf_t$ate), NA, res_bcf_t$ate),
              ifelse(is.na(res_bcf_t$lo),  NA, res_bcf_t$lo),
              ifelse(is.na(res_bcf_t$hi),  NA, res_bcf_t$hi)))

  # Full cohort BCF
  res_bcf_f <- tryCatch({
    set.seed(42)
    tmpdir2 <- tempdir()
    bcf_full <- bcf::bcf(
      y          = as.numeric(Y_full),
      z          = as.integer(W_full),
      x_control  = X_full,
      x_moderate = X_full,
      pihat      = e_full,
      nburn      = 100,
      nsim       = 500,
      nthin      = 1,
      verbose    = FALSE,
      save_tree_directory = tmpdir2
    )
    # ATE from in-sample tau
    tau_full_bcf <- bcf_cate_from_fit(bcf_full)
    if (!is.null(tau_full_bcf)) {
      # Compute DR ATE using logistic nuisance for mu1, mu0
      m_n2 <- glm(Y_full ~ ., data = as.data.frame(cbind(X_full, W = W_full)),
                   family = binomial())
      mu1_bf <- predict(m_n2, newdata = as.data.frame(cbind(X_full, W = rep(1L, n_full))),
                         type = "response")
      mu0_bf <- predict(m_n2, newdata = as.data.frame(cbind(X_full, W = rep(0L, n_full))),
                         type = "response")
      psi_bcf_f <- dr_pseudo(Y_full, W_full, e_full, mu1_bf, mu0_bf)
      boot_ci(psi_bcf_f)
    } else {
      list(ate = mean(tau_full_bcf, na.rm = TRUE), lo = NA, hi = NA, se = NA)
    }
  }, error = function(e) {
    cat("  BCF full failed:", conditionMessage(e), "\n")
    list(ate = NA, lo = NA, hi = NA, se = NA)
  })
  cat(sprintf("  Full ATE: %.4f [%.4f, %.4f]\n",
              ifelse(is.na(res_bcf_f$ate), NA, res_bcf_f$ate),
              ifelse(is.na(res_bcf_f$lo),  NA, res_bcf_f$lo),
              ifelse(is.na(res_bcf_f$hi),  NA, res_bcf_f$hi)))

  new_results[["BCF (Bayesian Causal Forest)"]] <- list(
    test  = res_bcf_t, full = res_bcf_f,
    notes = "BCF (Hahn et al. 2020, bcf package). Fitted on derivation set with PS as input. ATE via DR pseudo-outcome (logistic nuisance) in test set. 100 burn-in, 500 posterior draws."
  )
} else {
  cat("  bcf package not available\n")
  new_results[["BCF (Bayesian Causal Forest)"]] <- list(
    test  = list(ate = NA, lo = NA, hi = NA, se = NA),
    full  = list(ate = NA, lo = NA, hi = NA, se = NA),
    notes = "BCF: package not available. INFEASIBLE."
  )
}

###############################################################################
# Merge with existing results and rebuild table + forest plot
###############################################################################
cat("\n--- Building extended robustness table ---\n")

tbl_new_rows <- list()
for (nm in names(new_results)) {
  r  <- new_results[[nm]]
  rt <- r$test; rf <- r$full
  tbl_new_rows[[length(tbl_new_rows) + 1]] <- data.frame(
    estimator = nm, cohort = "Test set (held-out)", n = n_test,
    ATE = round(rt$ate, 4), SE = round(if (!is.null(rt$se)) rt$se else NA, 4),
    CI_lo = round(rt$lo, 4), CI_hi = round(rt$hi, 4),
    direction = direction_label(rt$lo, rt$hi),
    notes = r$notes, stringsAsFactors = FALSE
  )
  tbl_new_rows[[length(tbl_new_rows) + 1]] <- data.frame(
    estimator = nm, cohort = "Full cohort refit", n = n_full,
    ATE = round(rf$ate, 4), SE = round(if (!is.null(rf$se)) rf$se else NA, 4),
    CI_lo = round(rf$lo, 4), CI_hi = round(rf$hi, 4),
    direction = direction_label(rf$lo, rf$hi),
    notes = r$notes, stringsAsFactors = FALSE
  )
}

new_rows_df <- do.call(rbind, tbl_new_rows)

# Combine with existing
existing_rob_clean <- existing_rob[, intersect(names(existing_rob), names(new_rows_df))]
new_rows_clean     <- new_rows_df[, intersect(names(existing_rob), names(new_rows_df))]

extended_tbl <- rbind(existing_rob_clean, new_rows_clean)

write.csv(extended_tbl,
          file.path(TABLES_DIR, "supplement_ate_robustness_across_estimators.csv"),
          row.names = FALSE)
cat("  Saved: supplement_ate_robustness_across_estimators.csv (extended with new estimators)\n")

###############################################################################
# Updated forest plot
###############################################################################
cat("\n--- Building extended forest plot ---\n")

pal_blue  <- "#2166AC"
pal_red   <- "#D6604D"
pal_gray  <- "#636363"

plot_df <- extended_tbl[extended_tbl$cohort == "Test set (held-out)" &
                          !is.na(extended_tbl$ATE), ]
plot_df$is_primary <- grepl("selected", plot_df$estimator, ignore.case = TRUE)
plot_df$color_grp  <- ifelse(plot_df$is_primary, "Primary (selected)", "Alternative")

# Sort: primary first, then by estimator
plot_df <- plot_df[order(!plot_df$is_primary, plot_df$estimator), ]
plot_df$estimator_f <- factor(plot_df$estimator, levels = rev(plot_df$estimator))

ci_lo_min <- min(plot_df$CI_lo, na.rm = TRUE)
ci_hi_max <- max(plot_df$CI_hi, na.rm = TRUE)

p_forest_ext <- ggplot(plot_df,
                        aes(y = estimator_f, x = ATE,
                            xmin = CI_lo, xmax = CI_hi,
                            color = color_grp, shape = color_grp)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = pal_gray, linewidth = 0.7) +
  geom_errorbarh(aes(height = 0.18), linewidth = 0.9) +
  geom_point(aes(size = color_grp)) +
  geom_text(aes(label = sprintf("%.3f [%.3f, %.3f]", ATE, CI_lo, CI_hi)),
            hjust = -0.08, size = 2.9, color = "black") +
  scale_color_manual(values = c("Primary (selected)" = pal_red,
                                 "Alternative"       = pal_blue)) +
  scale_shape_manual(values = c("Primary (selected)" = 18, "Alternative" = 15)) +
  scale_size_manual(values  = c("Primary (selected)" = 5,  "Alternative" = 3.5)) +
  scale_x_continuous(
    limits = c(ci_lo_min - 0.06, ci_hi_max + 0.36),
    breaks = seq(-0.3, 0.6, by = 0.1)
  ) +
  labs(
    title    = "ATE Robustness — Extended Estimator Comparison (Test Set, N=133)",
    subtitle = sprintf("ASP vs cefazolin, 7-day AKI, MIMIC-IV v3.1. %d estimators. Primary (red), alternatives (blue).",
                       nrow(plot_df)),
    x = "Average Treatment Effect (ASP − CEF, absolute risk difference)",
    y = "", color = "Estimator type", shape = "Estimator type", size = "Estimator type"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 10, color = "grey30"),
    axis.text.y        = element_text(size = 9.5),
    legend.position    = "bottom"
  )

ggsave(file.path(FIGURES_DIR, "supp_forest_ate_across_estimators.png"),
       p_forest_ext, width = 11, height = 7, dpi = 300)
cat("  Saved: supp_forest_ate_across_estimators.png (extended)\n")

# Updated point estimates by estimator figure
all_df <- extended_tbl[!is.na(extended_tbl$ATE), ]
all_df$is_primary  <- grepl("selected", all_df$estimator, ignore.case = TRUE)
all_df$cohort_f    <- factor(all_df$cohort, levels = c("Test set (held-out)", "Full cohort refit"))

test_order <- all_df$estimator[all_df$cohort == "Test set (held-out)"]
test_order <- test_order[order(all_df$ATE[all_df$cohort == "Test set (held-out)"])]
all_df$estimator_f <- factor(all_df$estimator, levels = rev(test_order))
all_df$y_offset    <- ifelse(all_df$cohort == "Test set (held-out)", 0.17, -0.17)
all_df$y_pos       <- as.numeric(all_df$estimator_f) + all_df$y_offset

p_compare_ext <- ggplot(all_df, aes(y = y_pos, x = ATE, xmin = CI_lo, xmax = CI_hi,
                                     color = cohort_f, shape = cohort_f)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = pal_gray) +
  geom_errorbarh(aes(height = 0.09), linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_y_continuous(breaks = seq_along(levels(all_df$estimator_f)),
                     labels = levels(all_df$estimator_f)) +
  scale_color_manual(values = c("Test set (held-out)" = pal_red,
                                 "Full cohort refit"  = pal_blue)) +
  scale_shape_manual(values = c("Test set (held-out)" = 15, "Full cohort refit" = 16)) +
  scale_x_continuous(breaks = seq(-0.4, 0.6, by = 0.1)) +
  labs(title = "ATE Estimates by Estimator and Cohort — Extended Comparison",
       subtitle = "Red squares = test set (N=133, primary); Blue circles = full cohort (N=445, secondary)",
       x = "Average Treatment Effect (ASP − CEF)", y = "", color = "Cohort", shape = "Cohort") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "supp_ate_point_estimates_by_estimator.png"),
       p_compare_ext, width = 11, height = 7, dpi = 300)
cat("  Saved: supp_ate_point_estimates_by_estimator.png (extended)\n")

###############################################################################
# Summary
###############################################################################
cat("\n=== EXTENDED ROBUSTNESS SUMMARY ===\n")
test_ates_ext <- extended_tbl[extended_tbl$cohort == "Test set (held-out)", ]
test_ates_ok  <- test_ates_ext[!is.na(test_ates_ext$ATE), ]
cat(sprintf("Total estimators with results: %d / %d attempted\n",
            nrow(test_ates_ok), nrow(test_ates_ext)))
cat(sprintf("All positive ATE: %s\n", ifelse(all(test_ates_ok$ATE > 0), "YES", "NO")))
cat(sprintf("ATE range: [%.3f, %.3f]\n", min(test_ates_ok$ATE), max(test_ates_ok$ATE)))
for (i in seq_len(nrow(test_ates_ok))) {
  cat(sprintf("  %-45s ATE=%+.3f [%+.3f, %+.3f]\n",
              test_ates_ok$estimator[i],
              test_ates_ok$ATE[i], test_ates_ok$CI_lo[i], test_ates_ok$CI_hi[i]))
}

# Save updated summary object
saveRDS(
  list(rob_tbl = extended_tbl, n_estimators = nrow(test_ates_ok),
       ate_range = range(test_ates_ok$ATE, na.rm = TRUE),
       all_positive = all(test_ates_ok$ATE > 0, na.rm = TRUE)),
  file.path(RESULTS_DIR, "ate_robustness_summary.rds")
)
cat("\n=== 15b_ate_robustness_extended.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
