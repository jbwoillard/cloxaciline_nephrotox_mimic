###############################################################################
# 06_primary_causal_grf.R
# Primary causal analysis using Generalized Random Forest (causal forest)
# Output: results/model_objects.rds, results/ate_results.csv,
#         results/cate_predictions.csv, results/best_linear_projection.csv
###############################################################################

library(grf)
library(dplyr)
library(data.table)
library(ggplot2)

cat("=== 06_primary_causal_grf.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

set.seed(42)

WORK_DIR   <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
RESULTS_DIR <- file.path(WORK_DIR, "results")
FIGURES_DIR <- file.path(WORK_DIR, "figures")

# Load analysis dataset
anal <- readRDS(file.path(DERIVED_DIR, "analysis_dataset.rds"))
cat(sprintf("Analysis dataset: %d rows\n", nrow(anal)))

###############################################################################
# 1. Prepare X matrix (covariates)
###############################################################################
cat("\n--- Step 1: Prepare X matrix ---\n")

# Define covariate columns for causal forest
# All pre-t0 variables
x_vars <- c(
  "age", "female",
  # Labs
  "creat_baseline", "egfr_baseline", "wbc", "platelets",
  "bilirubin", "lactate", "glucose", "bicarbonate", "sodium",
  "potassium", "bun",
  # Comorbidities
  "ckd_any", "diabetes", "hypertension", "heart_failure",
  "liver_disease", "cancer", "copd",
  "charlson",
  # Clinical state
  "icu_at_t0", "vasopressor", "mech_vent",
  # Co-exposures
  "concomitant_vancomycin", "concomitant_aminoglycoside",
  "concomitant_loop_diuretic",
  # Time
  "hours_bc_to_t0",
  # Race (binary indicators)
  "aki_at_baseline"
)

x_vars_exist <- intersect(x_vars, names(anal))
cat(sprintf("Using %d covariates\n", length(x_vars_exist)))

# Handle race as dummy variable
anal[, race_white := as.integer(race_cat == "WHITE")]
anal[, race_black := as.integer(race_cat == "BLACK")]
anal[, race_hispanic := as.integer(race_cat == "HISPANIC")]
x_vars_exist <- c(x_vars_exist, "race_white", "race_black", "race_hispanic")

# Build X matrix with median imputation
X_full <- as.matrix(anal[, ..x_vars_exist])
cat(sprintf("X matrix dimensions: %d x %d\n", nrow(X_full), ncol(X_full)))

# Count missingness
n_miss <- sum(is.na(X_full))
cat(sprintf("Total missing values in X: %d (%.1f%%)\n",
            n_miss, 100*n_miss/length(X_full)))

# Median imputation + missingness indicators for high-miss variables
X_imp <- X_full
miss_indicators <- matrix(0, nrow=nrow(X_full), ncol=0)
miss_cols <- character(0)

for(j in seq_len(ncol(X_full))) {
  col_name <- colnames(X_full)[j]
  pct_miss <- mean(is.na(X_full[,j]))
  if(pct_miss > 0) {
    med_val <- median(X_full[,j], na.rm=TRUE)
    X_imp[is.na(X_imp[,j]), j] <- med_val
    if(pct_miss > 0.05) {  # Add indicator for >5% missing
      miss_ind <- as.numeric(is.na(X_full[,j]))
      miss_indicators <- cbind(miss_indicators, miss_ind)
      miss_cols <- c(miss_cols, paste0(col_name, "_miss"))
    }
  }
}

if(ncol(miss_indicators) > 0) {
  colnames(miss_indicators) <- miss_cols
  X_imp <- cbind(X_imp, miss_indicators)
  cat(sprintf("Added %d missingness indicators\n", length(miss_cols)))
}

cat(sprintf("Final X matrix: %d x %d\n", nrow(X_imp), ncol(X_imp)))

# W and Y vectors
# Only use rows with complete outcome
complete_rows <- !is.na(anal$AKI_7d)
cat(sprintf("Complete outcomes: %d/%d\n", sum(complete_rows), nrow(anal)))

X_use <- X_imp[complete_rows, ]
W_use <- anal$W[complete_rows]
Y_use <- anal$AKI_7d[complete_rows]
Y_cont_use <- anal$delta_creat_7d[complete_rows]  # continuous Y

cat(sprintf("Final N for analysis: %d (ASP: %d, CEF: %d)\n",
            length(W_use), sum(W_use==1), sum(W_use==0)))

###############################################################################
# 2. Train/test split (stratified by W)
###############################################################################
cat("\n--- Step 2: Train/test split (80/20) ---\n")

set.seed(42)
n_total <- length(W_use)

# Stratified split
idx_asp <- which(W_use == 1)
idx_cef <- which(W_use == 0)

train_asp <- sample(idx_asp, floor(0.8 * length(idx_asp)))
train_cef <- sample(idx_cef, floor(0.8 * length(idx_cef)))
train_idx <- sort(c(train_asp, train_cef))
test_idx  <- setdiff(seq_len(n_total), train_idx)

cat(sprintf("Train N: %d (ASP: %d, CEF: %d)\n",
            length(train_idx),
            sum(W_use[train_idx]==1),
            sum(W_use[train_idx]==0)))
cat(sprintf("Test N: %d (ASP: %d, CEF: %d)\n",
            length(test_idx),
            sum(W_use[test_idx]==1),
            sum(W_use[test_idx]==0)))

X_train <- X_use[train_idx, ]
W_train <- W_use[train_idx]
Y_train <- Y_use[train_idx]

X_test <- X_use[test_idx, ]
W_test <- W_use[test_idx]
Y_test <- Y_use[test_idx]

###############################################################################
# 3. Fit causal forest on training set
###############################################################################
cat("\n--- Step 3: Fit causal forest ---\n")

set.seed(42)
cf <- causal_forest(
  X = X_train,
  Y = Y_train,
  W = W_train,
  num.trees = 2000,
  seed = 42,
  tune.parameters = "all"
)

cat("Causal forest fitted\n")
cat(sprintf("  Num trees: %d\n", cf$num.trees))

# Overall OOB predictions
cat("\nOOB predictions summary (training set):\n")
tau_oob <- predict(cf)$predictions
cat(sprintf("  Mean CATE (OOB): %.4f\n", mean(tau_oob)))
cat(sprintf("  SD CATE (OOB): %.4f\n", sd(tau_oob)))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(tau_oob), max(tau_oob)))

###############################################################################
# 4. ATE estimation
###############################################################################
cat("\n--- Step 4: ATE estimation ---\n")

# Global ATE
ate_global <- average_treatment_effect(cf, target.sample="all")
cat(sprintf("\nATE (global): %.4f (SE: %.4f, 95%% CI: [%.4f, %.4f])\n",
            ate_global["estimate"],
            ate_global["std.err"],
            ate_global["estimate"] - 1.96*ate_global["std.err"],
            ate_global["estimate"] + 1.96*ate_global["std.err"]))

# ATE in treated (ATT)
ate_treated <- average_treatment_effect(cf, target.sample="treated")
cat(sprintf("ATT (treated): %.4f (SE: %.4f, 95%% CI: [%.4f, %.4f])\n",
            ate_treated["estimate"],
            ate_treated["std.err"],
            ate_treated["estimate"] - 1.96*ate_treated["std.err"],
            ate_treated["estimate"] + 1.96*ate_treated["std.err"]))

# ATE with overlap weighting
ate_overlap <- average_treatment_effect(cf, target.sample="overlap")
cat(sprintf("ATE (overlap): %.4f (SE: %.4f, 95%% CI: [%.4f, %.4f])\n",
            ate_overlap["estimate"],
            ate_overlap["std.err"],
            ate_overlap["estimate"] - 1.96*ate_overlap["std.err"],
            ate_overlap["estimate"] + 1.96*ate_overlap["std.err"]))

# Save ATE results
ate_results <- data.frame(
  estimand = c("ATE_global", "ATT_treated", "ATE_overlap"),
  estimate = c(ate_global["estimate"], ate_treated["estimate"], ate_overlap["estimate"]),
  std_err  = c(ate_global["std.err"],  ate_treated["std.err"],  ate_overlap["std.err"]),
  ci_lower = c(ate_global["estimate"] - 1.96*ate_global["std.err"],
               ate_treated["estimate"] - 1.96*ate_treated["std.err"],
               ate_overlap["estimate"] - 1.96*ate_overlap["std.err"]),
  ci_upper = c(ate_global["estimate"] + 1.96*ate_global["std.err"],
               ate_treated["estimate"] + 1.96*ate_treated["std.err"],
               ate_overlap["estimate"] + 1.96*ate_overlap["std.err"]),
  p_value  = c(2*pnorm(-abs(ate_global["estimate"]/ate_global["std.err"])),
               2*pnorm(-abs(ate_treated["estimate"]/ate_treated["std.err"])),
               2*pnorm(-abs(ate_overlap["estimate"]/ate_overlap["std.err"])))
)
rownames(ate_results) <- NULL
cat("\nATE results table:\n")
print(ate_results)

fwrite(ate_results, file.path(RESULTS_DIR, "ate_results.csv"))
cat("Saved: results/ate_results.csv\n")

###############################################################################
# 5. CATE predictions on test set
###############################################################################
cat("\n--- Step 5: CATE on test set ---\n")

cate_test <- predict(cf, newdata=X_test, estimate.variance=TRUE)
tau_test <- cate_test$predictions
var_test <- cate_test$variance.estimates

cat(sprintf("CATE predictions on test set:\n"))
cat(sprintf("  Mean: %.4f\n", mean(tau_test)))
cat(sprintf("  SD: %.4f\n", sd(tau_test)))
cat(sprintf("  IQR: [%.4f, %.4f]\n", quantile(tau_test, 0.25), quantile(tau_test, 0.75)))

# Build CATE predictions dataframe
# Get original row indices for test set
anal_complete <- anal[complete_rows]
anal_test <- anal_complete[test_idx]

cate_df <- data.frame(
  subject_id = anal_test$subject_id,
  hadm_id    = anal_test$hadm_id,
  W          = W_test,
  Y_observed = Y_test,
  CATE       = tau_test,
  CATE_se    = sqrt(pmax(var_test, 0)),
  CATE_lower = tau_test - 1.96*sqrt(pmax(var_test, 0)),
  CATE_upper = tau_test + 1.96*sqrt(pmax(var_test, 0))
)

fwrite(cate_df, file.path(RESULTS_DIR, "cate_predictions.csv"))
cat(sprintf("Saved: results/cate_predictions.csv (%d rows)\n", nrow(cate_df)))

###############################################################################
# 6. Also get CATE on ALL complete observations
###############################################################################
cat("\n--- Step 6: CATE on all complete observations ---\n")

cate_all <- predict(cf, newdata=X_use, estimate.variance=TRUE)
tau_all <- cate_all$predictions

anal_complete_dt <- anal[complete_rows]
anal_complete_dt[, CATE := tau_all]
anal_complete_dt[, CATE_se := sqrt(pmax(cate_all$variance.estimates, 0))]

cat(sprintf("CATE on all: mean=%.4f, SD=%.4f\n", mean(tau_all), sd(tau_all)))
cat(sprintf("Proportion with CATE > 0 (ASP harmful): %.1f%%\n",
            100*mean(tau_all > 0)))

###############################################################################
# 7. Variable importance
###############################################################################
cat("\n--- Step 7: Variable importance ---\n")

vi <- variable_importance(cf)
vi_df <- data.frame(
  variable = colnames(X_train),
  importance = as.numeric(vi)
)
vi_df <- vi_df[order(-vi_df$importance), ]

cat("Top 15 most important variables:\n")
print(head(vi_df, 15))

fwrite(vi_df, file.path(RESULTS_DIR, "variable_importance.csv"))
cat("Saved: results/variable_importance.csv\n")

###############################################################################
# 8. Best linear projection (heterogeneity)
###############################################################################
cat("\n--- Step 8: Best linear projection ---\n")

# Select key heterogeneity variables for BLP
blp_vars <- c("age", "creat_baseline", "egfr_baseline", "ckd_any",
              "icu_at_t0", "vasopressor", "charlson", "wbc",
              "concomitant_vancomycin", "aki_at_baseline")
blp_vars_exist <- intersect(blp_vars, colnames(X_train))
A_train <- X_train[, blp_vars_exist, drop=FALSE]

# Handle any remaining NA in A
A_train_imp <- A_train
for(j in seq_len(ncol(A_train_imp))) {
  A_train_imp[is.na(A_train_imp[,j]), j] <- median(A_train_imp[,j], na.rm=TRUE)
}

blp_fit <- best_linear_projection(cf, A=A_train_imp)
cat("Best Linear Projection:\n")
print(blp_fit)

# Save BLP - extract from printed output via capture
blp_summary <- capture.output(print(blp_fit))
# Use the raw blp_fit object attributes
blp_df <- tryCatch({
  # blp_fit is an lmtest object or similar
  sm <- summary.default(blp_fit)
  cf_mat <- sm[["coefficients"]]
  df_out <- as.data.frame(cf_mat)
  df_out$variable <- rownames(df_out)
  df_out
}, error = function(e) {
  # Fallback: just save as text
  data.frame(variable="blp_output", info=paste(blp_summary, collapse="\n"))
})
fwrite(blp_df, file.path(RESULTS_DIR, "best_linear_projection.csv"))
cat("Saved: results/best_linear_projection.csv\n")

###############################################################################
# 9. Propensity score (from regression forest)
###############################################################################
cat("\n--- Step 9: Propensity scores ---\n")

# Get propensity scores from the forest's W.hat
W_hat_train <- cf$W.hat
cat(sprintf("Propensity score (W.hat) range: [%.3f, %.3f]\n",
            min(W_hat_train), max(W_hat_train)))
cat(sprintf("Mean W.hat (should be ~%.2f): %.3f\n",
            mean(W_train), mean(W_hat_train)))

###############################################################################
# 10. Fit causal forest on full data for final estimates
###############################################################################
cat("\n--- Step 10: Refit on full dataset ---\n")

set.seed(42)
cf_full <- causal_forest(
  X = X_use,
  Y = Y_use,
  W = W_use,
  num.trees = 2000,
  seed = 42,
  tune.parameters = "all"
)

ate_full <- average_treatment_effect(cf_full, target.sample="overlap")
cat(sprintf("ATE (overlap, full N=%d): %.4f (SE: %.4f, p=%.4f)\n",
            length(Y_use),
            ate_full["estimate"],
            ate_full["std.err"],
            2*pnorm(-abs(ate_full["estimate"]/ate_full["std.err"]))))

# Calibration test
cal_test <- test_calibration(cf_full)
cat("\nCalibration test:\n")
print(cal_test)

###############################################################################
# 11. Save model objects
###############################################################################
cat("\n--- Step 11: Save model objects ---\n")

model_objects <- list(
  cf_train = cf,
  cf_full  = cf_full,
  X_use = X_use,
  W_use = W_use,
  Y_use = Y_use,
  train_idx = train_idx,
  test_idx  = test_idx,
  ate_global  = ate_global,
  ate_treated = ate_treated,
  ate_overlap = ate_overlap,
  ate_full    = ate_full,
  cate_all_predictions = tau_all,
  variable_importance = vi_df,
  blp_fit = blp_fit,
  covariate_names = colnames(X_use)
)

saveRDS(model_objects, file.path(RESULTS_DIR, "model_objects.rds"))
cat("Saved: results/model_objects.rds\n")

# Also add CATE to full analysis dataset
anal_complete_dt[, W_hat := cf_full$W.hat]
anal_complete_dt[, Y_hat := cf_full$Y.hat]
saveRDS(anal_complete_dt, file.path(DERIVED_DIR, "analysis_dataset_with_cate.rds"))
cat("Saved: derived/analysis_dataset_with_cate.rds\n")

cat("\n=== 06_primary_causal_grf.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
