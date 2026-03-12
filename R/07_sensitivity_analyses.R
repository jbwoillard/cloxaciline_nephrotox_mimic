###############################################################################
# 07_sensitivity_analyses.R
# Sensitivity analyses for robustness
# Output: results/sensitivity_results.csv
###############################################################################

library(grf)
library(dplyr)
library(data.table)
library(lubridate)

cat("=== 07_sensitivity_analyses.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

set.seed(42)

WORK_DIR    <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
RESULTS_DIR <- file.path(WORK_DIR, "results")
HOSP_DIR    <- file.path(WORK_DIR, "mimic-iv-3.1/hosp")

# Load analysis dataset
anal <- readRDS(file.path(DERIVED_DIR, "analysis_dataset.rds"))
model_objs <- readRDS(file.path(RESULTS_DIR, "model_objects.rds"))

cat(sprintf("Analysis dataset: %d rows\n", nrow(anal)))
cat(sprintf("ASP: %d, CEF: %d\n", sum(anal$W==1), sum(anal$W==0)))

# Reference: main analysis ATE
ate_main <- model_objs$ate_full
cat(sprintf("\nMain ATE (overlap): %.4f (SE: %.4f)\n",
            ate_main["estimate"], ate_main["std.err"]))

###############################################################################
# Helper: run causal forest and extract ATE
###############################################################################
run_causal_forest <- function(data, Y_col, label, additional_info="") {
  cat(sprintf("\n  Running: %s (N=%d)\n", label, nrow(data)))

  # Build X matrix same as main analysis
  x_vars <- c(
    "age", "female",
    "creat_baseline", "egfr_baseline", "wbc", "platelets",
    "bilirubin", "lactate", "glucose", "bicarbonate", "sodium",
    "potassium", "bun",
    "ckd_any", "diabetes", "hypertension", "heart_failure",
    "liver_disease", "cancer", "copd",
    "charlson",
    "icu_at_t0", "vasopressor", "mech_vent",
    "concomitant_vancomycin", "concomitant_aminoglycoside",
    "concomitant_loop_diuretic",
    "hours_bc_to_t0", "aki_at_baseline",
    "race_white", "race_black", "race_hispanic"
  )

  # Add race dummies if missing
  if(!"race_white" %in% names(data)) {
    data[, race_white := as.integer(race_cat == "WHITE")]
    data[, race_black := as.integer(race_cat == "BLACK")]
    data[, race_hispanic := as.integer(race_cat == "HISPANIC")]
  }

  x_vars_exist <- intersect(x_vars, names(data))
  X <- as.matrix(data[, ..x_vars_exist])

  # Median imputation
  for(j in seq_len(ncol(X))) {
    X[is.na(X[,j]), j] <- median(X[,j], na.rm=TRUE)
  }

  W <- data$W
  Y <- data[[Y_col]]

  # Keep complete
  comp <- !is.na(Y) & !is.na(W)
  X <- X[comp, ]
  W <- W[comp]
  Y <- Y[comp]

  n <- length(Y)
  cat(sprintf("    Complete N: %d (W=1: %d, W=0: %d)\n",
              n, sum(W==1), sum(W==0)))

  if(sum(W==1) < 10 || sum(W==0) < 10 || n < 50) {
    cat("    Insufficient data - skipping\n")
    return(data.frame(
      analysis = label,
      outcome = Y_col,
      N = n,
      N_ASP = sum(W==1),
      N_CEF = sum(W==0),
      ATE_estimate = NA,
      ATE_se = NA,
      ATE_lower = NA,
      ATE_upper = NA,
      ATE_pvalue = NA,
      note = "Insufficient data"
    ))
  }

  set.seed(42)
  cf <- tryCatch(
    causal_forest(X=X, Y=Y, W=W, num.trees=1000, seed=42,
                  tune.parameters="all"),
    error = function(e) {
      cat(sprintf("    Error fitting forest: %s\n", e$message))
      return(NULL)
    }
  )

  if(is.null(cf)) {
    return(data.frame(
      analysis=label, outcome=Y_col, N=n,
      N_ASP=sum(W==1), N_CEF=sum(W==0),
      ATE_estimate=NA, ATE_se=NA, ATE_lower=NA, ATE_upper=NA, ATE_pvalue=NA,
      note="Forest error"
    ))
  }

  ate <- average_treatment_effect(cf, target.sample="overlap")

  cat(sprintf("    ATE: %.4f (SE: %.4f, p=%.4f)\n",
              ate["estimate"], ate["std.err"],
              2*pnorm(-abs(ate["estimate"]/ate["std.err"]))))

  data.frame(
    analysis = label,
    outcome = Y_col,
    N = n,
    N_ASP = sum(W==1),
    N_CEF = sum(W==0),
    ATE_estimate = ate["estimate"],
    ATE_se       = ate["std.err"],
    ATE_lower    = ate["estimate"] - 1.96*ate["std.err"],
    ATE_upper    = ate["estimate"] + 1.96*ate["std.err"],
    ATE_pvalue   = 2*pnorm(-abs(ate["estimate"]/ate["std.err"])),
    note = additional_info
  )
}

###############################################################################
# Initialize results list
###############################################################################
sens_results <- list()

# Add main analysis as reference
sens_results[["main"]] <- data.frame(
  analysis = "Primary analysis",
  outcome = "AKI_7d",
  N = sum(!is.na(anal$AKI_7d)),
  N_ASP = sum(anal$W==1 & !is.na(anal$AKI_7d)),
  N_CEF = sum(anal$W==0 & !is.na(anal$AKI_7d)),
  ATE_estimate = model_objs$ate_full["estimate"],
  ATE_se       = model_objs$ate_full["std.err"],
  ATE_lower    = model_objs$ate_full["estimate"] - 1.96*model_objs$ate_full["std.err"],
  ATE_upper    = model_objs$ate_full["estimate"] + 1.96*model_objs$ate_full["std.err"],
  ATE_pvalue   = 2*pnorm(-abs(model_objs$ate_full["estimate"]/model_objs$ate_full["std.err"])),
  note = "Main analysis (overlap weighting)"
)

###############################################################################
# A. As-treated: require >= 48h on initial treatment (no early switching)
###############################################################################
cat("\n=== A. As-treated: exclude early switchers (<48h) ===\n")

# Switchers: received both drugs in window
# For as-treated: exclude anyone who switched within 48h
cohort_full <- readRDS(file.path(DERIVED_DIR, "cohort_treated.rds"))

# Define early switchers from emar data
library(duckdb); library(DBI)
emar_path <- file.path(HOSP_DIR, "emar.csv.gz")
con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
dbExecute(con, "SET threads=4; SET memory_limit='8GB';")

subj_ids <- unique(anal$subject_id)
sid_list <- paste(subj_ids, collapse=",")

emar_sw <- dbGetQuery(con, paste0("
  SELECT subject_id, hadm_id, medication, event_txt, charttime
  FROM read_csv_auto('", emar_path, "', compression='gzip')
  WHERE subject_id IN (", sid_list, ")
    AND (UPPER(medication) LIKE '%NAFCILLIN%'
      OR UPPER(medication) LIKE '%OXACILLIN%'
      OR UPPER(medication) LIKE '%CEFAZOLIN%')
    AND event_txt IN ('Administered', 'Applied', 'Started', 'Restarted', 'New Bag')
"))
dbDisconnect(con, shutdown=TRUE)

emar_sw <- as.data.table(emar_sw)
emar_sw[, charttime := as.POSIXct(charttime)]
emar_sw[, drug_group := fifelse(
  grepl("NAFCILLIN|OXACILLIN", toupper(medication)) &
    !grepl("DICLOXACILLIN", toupper(medication)),
  "ASP", "CEF"
)]

# Merge with cohort to get t0 and W
emar_sw2 <- merge(emar_sw, cohort_full[, .(subject_id, hadm_id, t0, W)],
                  by=c("subject_id","hadm_id"))
emar_sw2[, hours_from_t0 := as.numeric(difftime(charttime, t0, units="hours"))]

# Find those who switched within 48h
switched_48h <- emar_sw2[
  hours_from_t0 >= 0 & hours_from_t0 <= 48 &
  ((W==1 & drug_group=="CEF") | (W==0 & drug_group=="ASP")),
  .(subject_id, hadm_id, early_switch=TRUE)
]
switched_48h <- unique(switched_48h)
cat(sprintf("Early switchers (<48h): %d\n", nrow(switched_48h)))

anal_as_treated <- merge(anal, switched_48h, by=c("subject_id","hadm_id"), all.x=TRUE)
anal_as_treated[is.na(early_switch), early_switch := FALSE]
anal_as_treated <- anal_as_treated[early_switch==FALSE]
cat(sprintf("As-treated cohort: %d (removed %d)\n", nrow(anal_as_treated), nrow(anal)-nrow(anal_as_treated)))

sens_results[["A_as_treated"]] <- run_causal_forest(
  anal_as_treated, "AKI_7d", "A. As-treated (exclude early switchers)",
  "Excluded switchers within 48h"
)

###############################################################################
# B. Exclude AKI at baseline
###############################################################################
cat("\n=== B. Exclude baseline AKI ===\n")

anal_no_aki0 <- anal[is.na(aki_at_baseline) | aki_at_baseline==0]
cat(sprintf("Excluded baseline AKI: %d → %d\n", nrow(anal), nrow(anal_no_aki0)))

sens_results[["B_no_baseline_aki"]] <- run_causal_forest(
  anal_no_aki0, "AKI_7d", "B. Exclude AKI at baseline",
  "AKI at baseline defined as creatinine rise >= 0.3 or >= 1.5x pre-admission"
)

###############################################################################
# C. Exclude early switchers (<24h) - stricter
###############################################################################
cat("\n=== C. Exclude early switchers (<24h) ===\n")

switched_24h <- emar_sw2[
  hours_from_t0 >= 0 & hours_from_t0 <= 24 &
  ((W==1 & drug_group=="CEF") | (W==0 & drug_group=="ASP")),
  .(subject_id, hadm_id, early_switch_24h=TRUE)
]
switched_24h <- unique(switched_24h)
cat(sprintf("Early switchers (<24h): %d\n", nrow(switched_24h)))

anal_no_early_switch <- merge(anal, switched_24h, by=c("subject_id","hadm_id"), all.x=TRUE)
anal_no_early_switch[is.na(early_switch_24h), early_switch_24h := FALSE]
anal_no_early_switch <- anal_no_early_switch[early_switch_24h==FALSE]

sens_results[["C_no_switch_24h"]] <- run_causal_forest(
  anal_no_early_switch, "AKI_7d", "C. Exclude early switchers (<24h)",
  "Excluded those switching therapy within 24h of t0"
)

###############################################################################
# D. Continuous outcome: delta_creat_7d
###############################################################################
cat("\n=== D. Continuous outcome: delta creatinine at 7d ===\n")

sens_results[["D_delta_creat"]] <- run_causal_forest(
  anal, "delta_creat_7d", "D. Continuous outcome (delta Cr 7d)",
  "Continuous: peak minus baseline creatinine (mg/dL) at 7 days"
)

###############################################################################
# E. Composite outcome: AKI or death within 7d
###############################################################################
cat("\n=== E. Composite outcome: AKI or death within 7d ===\n")

sens_results[["E_composite"]] <- run_causal_forest(
  anal, "composite_aki_death", "E. Composite outcome (AKI or 7d death)",
  "Composite of AKI at 7d or death within 7d"
)

###############################################################################
# F. MSSA confirmed only (exclude unknown susceptibility)
###############################################################################
cat("\n=== F. MSSA confirmed only ===\n")

anal_mssa_confirmed <- anal[mssa_confirmed==TRUE]
cat(sprintf("MSSA confirmed only: %d (from %d)\n", nrow(anal_mssa_confirmed), nrow(anal)))

sens_results[["F_mssa_confirmed"]] <- run_causal_forest(
  anal_mssa_confirmed, "AKI_7d", "F. MSSA confirmed only",
  "Restricted to patients with confirmed oxacillin-susceptible S. aureus"
)

###############################################################################
# G. Exclude switchers completely (neither drug in both)
###############################################################################
cat("\n=== G. Exclude ALL switchers (received both drugs) ===\n")

anal_no_switch <- anal[received_both==FALSE]
cat(sprintf("No switchers cohort: %d\n", nrow(anal_no_switch)))

sens_results[["G_no_switch"]] <- run_causal_forest(
  anal_no_switch, "AKI_7d", "G. Exclude all switchers (received both drugs)",
  "Active comparator: monotherapy only"
)

###############################################################################
# H. Restricted to ICU patients
###############################################################################
cat("\n=== H. ICU patients only ===\n")

anal_icu <- anal[icu_at_t0==1]
cat(sprintf("ICU at t0: %d\n", nrow(anal_icu)))

sens_results[["H_icu_only"]] <- run_causal_forest(
  anal_icu, "AKI_7d", "H. ICU patients only",
  "Restricted to patients in ICU at time zero"
)

###############################################################################
# Combine all results
###############################################################################
cat("\n=== Combining all sensitivity results ===\n")

all_sens <- rbindlist(lapply(sens_results, as.data.table), fill=TRUE)
rownames(all_sens) <- NULL

# Reorder columns
cols_order <- c("analysis", "outcome", "N", "N_ASP", "N_CEF",
                "ATE_estimate", "ATE_se", "ATE_lower", "ATE_upper", "ATE_pvalue", "note")
cols_exist <- intersect(cols_order, names(all_sens))
all_sens <- all_sens[, ..cols_exist]

cat("\nAll sensitivity analysis results:\n")
print(all_sens[, .(analysis, N, N_ASP, N_CEF,
                    ATE=round(ATE_estimate,3),
                    SE=round(ATE_se,3),
                    CI_lower=round(ATE_lower,3),
                    CI_upper=round(ATE_upper,3),
                    pvalue=round(ATE_pvalue,4))])

fwrite(all_sens, file.path(RESULTS_DIR, "sensitivity_results.csv"))
cat("\nSaved: results/sensitivity_results.csv\n")

cat("\n=== 07_sensitivity_analyses.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
