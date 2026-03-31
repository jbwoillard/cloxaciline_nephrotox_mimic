###############################################################################
# 12_unified_derivation_test_pipeline.R
# UNIFIED derivation/test split for the consolidated ITE/HTE pipeline
# Follows: Buell et al. JAMA 2024 and Munroe et al. Lancet Resp Med 2025
#
# Outputs:
#   derived/derivation_set_unified.rds
#   derived/test_set_unified.rds
#   derived/imputation_rules_unified.rds
#   tables/derivation_test_split_summary.csv
###############################################################################

library(dplyr)

cat("=== 12_unified_derivation_test_pipeline.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

set.seed(42)

WORK_DIR    <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
TABLES_DIR  <- file.path(WORK_DIR, "tables")

###############################################################################
# 1. Load analysis dataset
###############################################################################
cat("--- Step 1: Load dataset ---\n")
ds <- readRDS(file.path(DERIVED_DIR, "analysis_dataset.rds"))
cat(sprintf("Full dataset: %d rows, %d columns\n", nrow(ds), ncol(ds)))
cat(sprintf("Eligible N: %d (all rows in analysis_dataset.rds)\n", nrow(ds)))

# Primary outcome: AKI_7d (1 NA exists)
# W: binary treatment (1=ASP, 0=CEF)
ds_complete <- ds[!is.na(ds$AKI_7d), ]
cat(sprintf("Analyzable N (non-missing 7-day AKI outcome): %d\n", nrow(ds_complete)))
cat(sprintf("  Excluded (missing AKI_7d): %d\n", nrow(ds) - nrow(ds_complete)))

cat("Treatment distribution:\n")
print(table(ds_complete$W, useNA = "always"))
cat("Outcome (AKI_7d) distribution:\n")
print(table(ds_complete$AKI_7d, useNA = "always"))

###############################################################################
# 2. Define covariate set (all pre-time-zero baseline variables)
###############################################################################
cat("\n--- Step 2: Define covariate set ---\n")

x_vars <- c(
  "age", "female",
  "creat_baseline", "egfr_baseline",
  "wbc", "platelets", "bilirubin", "lactate",
  "glucose", "bicarbonate", "sodium", "potassium", "bun",
  "ckd_any", "diabetes", "hypertension", "heart_failure",
  "liver_disease", "cancer", "copd", "charlson",
  "icu_at_t0", "vasopressor", "mech_vent",
  "concomitant_vancomycin", "concomitant_aminoglycoside",
  "concomitant_loop_diuretic",
  "aki_at_baseline",
  "hours_bc_to_t0.x"
)

# Race dummy encoding
ds_complete$race_white    <- as.integer(ds_complete$race_cat == "WHITE")
ds_complete$race_black    <- as.integer(ds_complete$race_cat == "BLACK")
ds_complete$race_hispanic <- as.integer(ds_complete$race_cat == "HISPANIC")
x_vars <- c(x_vars, "race_white", "race_black", "race_hispanic")

x_vars_exist <- intersect(x_vars, names(ds_complete))
cat(sprintf("Covariates available: %d\n", length(x_vars_exist)))

###############################################################################
# 3. Stratified 70/30 derivation/test split
#    Stratify by treatment arm (W) x AKI status to ensure balance
###############################################################################
cat("\n--- Step 3: Stratified 70/30 split (seed=42) ---\n")

ds_complete$stratum <- paste0("W", ds_complete$W, "_Y", ds_complete$AKI_7d)
cat("Stratum sizes before split:\n")
print(table(ds_complete$stratum))

set.seed(42)
strata_levels <- unique(ds_complete$stratum)
deriv_idx <- c()
for (s in strata_levels) {
  idx_s     <- which(ds_complete$stratum == s)
  n_deriv_s <- round(0.70 * length(idx_s))
  deriv_idx <- c(deriv_idx, sample(idx_s, n_deriv_s))
}
deriv_idx <- sort(deriv_idx)
test_idx  <- setdiff(seq_len(nrow(ds_complete)), deriv_idx)

derivation_set <- ds_complete[deriv_idx, ]
test_set       <- ds_complete[test_idx,  ]

cat(sprintf("Derivation set: N=%d (ASP=%d, CEF=%d, AKI=%d, AKI%%=%.1f%%)\n",
            nrow(derivation_set),
            sum(derivation_set$W == 1),
            sum(derivation_set$W == 0),
            sum(derivation_set$AKI_7d == 1),
            100 * mean(derivation_set$AKI_7d == 1)))
cat(sprintf("Test set:       N=%d (ASP=%d, CEF=%d, AKI=%d, AKI%%=%.1f%%)\n",
            nrow(test_set),
            sum(test_set$W == 1),
            sum(test_set$W == 0),
            sum(test_set$AKI_7d == 1),
            100 * mean(test_set$AKI_7d == 1)))

###############################################################################
# 4. Imputation rules derived from derivation set ONLY (no leakage)
#    Strategy: median imputation for continuous, mode for binary
#    Missingness indicators for variables with >5% missing in derivation
###############################################################################
cat("\n--- Step 4: Derive imputation rules from derivation set only ---\n")

X_deriv_raw <- derivation_set[, x_vars_exist]
X_test_raw  <- test_set[,       x_vars_exist]

miss_pct_deriv <- sapply(X_deriv_raw, function(x) mean(is.na(x)))
cat("Missingness in derivation set (variables with any missing):\n")
miss_present <- sort(miss_pct_deriv[miss_pct_deriv > 0], decreasing = TRUE)
print(round(miss_present, 3))

# Mode helper for binary variables
mode_val <- function(x) {
  x <- x[!is.na(x)]
  as.numeric(names(sort(table(x), decreasing=TRUE))[1])
}

imputation_params <- list()
miss_indicator_vars <- c()

X_deriv_imp <- X_deriv_raw
X_test_imp  <- X_test_raw

for (v in x_vars_exist) {
  pct_miss_v <- mean(is.na(X_deriv_raw[[v]]))
  if (pct_miss_v > 0) {
    # Determine if binary (0/1 only in non-missing values)
    vals_nonmiss <- X_deriv_raw[[v]][!is.na(X_deriv_raw[[v]])]
    is_binary <- all(vals_nonmiss %in% c(0, 1))

    if (is_binary) {
      fill_val <- mode_val(X_deriv_raw[[v]])
      fill_type <- "mode"
    } else {
      fill_val <- median(X_deriv_raw[[v]], na.rm = TRUE)
      fill_type <- "median"
    }

    imputation_params[[v]] <- list(
      value    = fill_val,
      type     = fill_type,
      pct_miss = pct_miss_v
    )

    # Impute derivation set
    X_deriv_imp[[v]][is.na(X_deriv_imp[[v]])] <- fill_val
    # Impute test set using derivation-derived value (no leakage)
    X_test_imp[[v]][is.na(X_test_imp[[v]])]   <- fill_val

    # Missingness indicator for variables with >5% missing in derivation
    if (pct_miss_v > 0.05) {
      miss_indicator_vars <- c(miss_indicator_vars, v)
    }
  }
}

# Append missingness indicators
for (v in miss_indicator_vars) {
  ind_name <- paste0(v, "_miss")
  X_deriv_imp[[ind_name]] <- as.integer(is.na(X_deriv_raw[[v]]))
  X_test_imp[[ind_name]]  <- as.integer(is.na(X_test_raw[[v]]))
  cat(sprintf("  Added missingness indicator: %s\n", ind_name))
}

cat(sprintf("Imputed derivation matrix: %d x %d\n", nrow(X_deriv_imp), ncol(X_deriv_imp)))
cat(sprintf("Imputed test matrix:       %d x %d\n", nrow(X_test_imp),  ncol(X_test_imp)))
cat(sprintf("Remaining NAs - derivation: %d, test: %d\n",
            sum(is.na(X_deriv_imp)), sum(is.na(X_test_imp))))

###############################################################################
# 5. Build final objects
###############################################################################
cat("\n--- Step 5: Build and save objects ---\n")

final_x_cols <- names(X_deriv_imp)
cat(sprintf("Final covariate columns: %d\n", length(final_x_cols)))

derivation_final <- list(
  data       = derivation_set,
  X          = as.matrix(X_deriv_imp),
  W          = derivation_set$W,
  Y          = derivation_set$AKI_7d,
  x_cols     = final_x_cols,
  subject_id = derivation_set$subject_id,
  hadm_id    = derivation_set$hadm_id,
  stratum    = derivation_set$stratum,
  n          = nrow(derivation_set),
  n_asp      = sum(derivation_set$W == 1),
  n_cef      = sum(derivation_set$W == 0),
  n_aki      = sum(derivation_set$AKI_7d == 1),
  pct_aki    = mean(derivation_set$AKI_7d == 1)
)

test_final <- list(
  data       = test_set,
  X          = as.matrix(X_test_imp),
  W          = test_set$W,
  Y          = test_set$AKI_7d,
  x_cols     = final_x_cols,
  subject_id = test_set$subject_id,
  hadm_id    = test_set$hadm_id,
  stratum    = test_set$stratum,
  n          = nrow(test_set),
  n_asp      = sum(test_set$W == 1),
  n_cef      = sum(test_set$W == 0),
  n_aki      = sum(test_set$AKI_7d == 1),
  pct_aki    = mean(test_set$AKI_7d == 1)
)

imputation_rules <- list(
  params              = imputation_params,
  miss_indicator_vars = miss_indicator_vars,
  x_vars_base         = x_vars_exist,
  final_x_cols        = final_x_cols
)

saveRDS(derivation_final, file.path(DERIVED_DIR, "derivation_set_unified.rds"))
cat("Saved: derived/derivation_set_unified.rds\n")

saveRDS(test_final, file.path(DERIVED_DIR, "test_set_unified.rds"))
cat("Saved: derived/test_set_unified.rds\n")

saveRDS(imputation_rules, file.path(DERIVED_DIR, "imputation_rules_unified.rds"))
cat("Saved: derived/imputation_rules_unified.rds\n")

###############################################################################
# 6. Create split summary table
###############################################################################
cat("\n--- Step 6: Create split summary table ---\n")

split_summary <- data.frame(
  set     = c("Derivation", "Test", "Total"),
  n_total = c(nrow(derivation_set), nrow(test_set), nrow(ds_complete)),
  n_asp   = c(sum(derivation_set$W==1), sum(test_set$W==1), sum(ds_complete$W==1)),
  n_cef   = c(sum(derivation_set$W==0), sum(test_set$W==0), sum(ds_complete$W==0)),
  n_aki   = c(sum(derivation_set$AKI_7d==1), sum(test_set$AKI_7d==1), sum(ds_complete$AKI_7d==1)),
  pct_aki = c(100*mean(derivation_set$AKI_7d==1),
              100*mean(test_set$AKI_7d==1),
              100*mean(ds_complete$AKI_7d==1))
)
split_summary$pct_aki <- round(split_summary$pct_aki, 1)

cat("Split summary:\n")
print(split_summary)

write.csv(split_summary,
          file.path(TABLES_DIR, "derivation_test_split_summary.csv"),
          row.names = FALSE)
cat("Saved: tables/derivation_test_split_summary.csv (overwritten with unified version)\n")

cat("\n=== SUMMARY ===\n")
cat(sprintf("Eligible N:         %d\n", nrow(ds)))
cat(sprintf("Analyzable N:       %d (excluded %d with missing AKI outcome)\n",
            nrow(ds_complete), nrow(ds) - nrow(ds_complete)))
cat(sprintf("Derivation set:     N=%d (%.0f%%)\n", nrow(derivation_set), 100*nrow(derivation_set)/nrow(ds_complete)))
cat(sprintf("Test set:           N=%d (%.0f%%)\n", nrow(test_set), 100*nrow(test_set)/nrow(ds_complete)))
cat(sprintf("Covariates (final): %d\n", length(final_x_cols)))

cat("\n=== 12_unified_derivation_test_pipeline.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
