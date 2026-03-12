###############################################################################
# 04_build_covariates.R
# Build all pre-t0 covariates for the analysis
# Output: derived/covariates.rds, tables/covariate_dictionary.csv
###############################################################################

library(duckdb)
library(DBI)
library(dplyr)
library(lubridate)
library(data.table)
library(stringr)

cat("=== 04_build_covariates.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

WORK_DIR <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
HOSP_DIR <- file.path(WORK_DIR, "mimic-iv-3.1/hosp")
ICU_DIR  <- file.path(WORK_DIR, "mimic-iv-3.1/icu")
DERIVED_DIR <- file.path(WORK_DIR, "derived")
TABLES_DIR  <- file.path(WORK_DIR, "tables")

# Load cohort
cohort <- readRDS(file.path(DERIVED_DIR, "cohort_treated.rds"))
cat(sprintf("Cohort loaded: %d patients\n", nrow(cohort)))

subj_ids <- unique(cohort$subject_id)
hadm_ids <- unique(cohort$hadm_id)
sid_list  <- paste(subj_ids, collapse=",")
hadm_list <- paste(hadm_ids, collapse=",")

open_duck <- function() {
  con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
  dbExecute(con, "SET threads=4; SET memory_limit='8GB';")
  con
}

###############################################################################
# 1. Demographics (from cohort)
###############################################################################
cat("\n--- Step 1: Demographics ---\n")

covars <- cohort[, .(
  subject_id, hadm_id, t0,
  age, gender, race,
  hours_bc_to_t0,
  admission_type
)]

# Gender binary
covars[, female := as.integer(gender == "F")]

# Race categories
covars[, race_cat := fcase(
  grepl("WHITE", toupper(race)), "WHITE",
  grepl("BLACK|AFRICAN", toupper(race)), "BLACK",
  grepl("HISPANIC|LATINO", toupper(race)), "HISPANIC",
  grepl("ASIAN", toupper(race)), "ASIAN",
  default = "OTHER/UNKNOWN"
)]

cat("Race categories:\n")
print(table(covars$race_cat))

###############################################################################
# 2. Baseline labs from labevents
###############################################################################
cat("\n--- Step 2: Baseline labs ---\n")

# Lab itemids:
# 50912 = creatinine (already have from outcomes)
# 51301 = WBC
# 51755 = WBC (alternate)
# 51265 = platelets
# 50885 = bilirubin total
# 50813 = lactate
# 50931 = glucose
# 50882 = bicarbonate
# 50971 = potassium
# 50983 = sodium
# 51006 = BUN
# 50893 = calcium

lab_itemids <- c(50912, 51301, 51755, 51265, 50885, 50813,
                 50931, 50882, 50971, 50983, 51006, 50893)
lab_ids_str <- paste(lab_itemids, collapse=",")

con <- open_duck()
lab_path <- file.path(HOSP_DIR, "labevents.csv.gz")

lab_query <- paste0("
  SELECT subject_id, hadm_id, itemid, charttime, valuenum
  FROM read_csv_auto('", lab_path, "', compression='gzip')
  WHERE itemid IN (", lab_ids_str, ")
    AND subject_id IN (", sid_list, ")
    AND valuenum IS NOT NULL
    AND valuenum > 0
")

labs_raw <- dbGetQuery(con, lab_query)
dbDisconnect(con, shutdown=TRUE)

cat(sprintf("Lab rows loaded: %d\n", nrow(labs_raw)))

labs_raw <- as.data.table(labs_raw)
labs_raw[, charttime := as.POSIXct(charttime)]

# Merge with cohort to get t0
labs_t0 <- merge(labs_raw, cohort[, .(subject_id, hadm_id, t0)],
                 by=c("subject_id","hadm_id"))
labs_t0[, hours_from_t0 := as.numeric(difftime(charttime, t0, units="hours"))]

# Baseline window: [-48h, 0h] before t0
labs_baseline <- labs_t0[hours_from_t0 >= -48 & hours_from_t0 <= 0]

cat(sprintf("Lab rows in baseline window: %d\n", nrow(labs_baseline)))

# For each patient/hadm/itemid, take median value in baseline window
labs_summary <- labs_baseline[, .(
  value = median(valuenum, na.rm=TRUE)
), by=.(subject_id, hadm_id, itemid)]

# Pivot wide
labs_wide <- dcast(labs_summary, subject_id + hadm_id ~ itemid, value.var="value")

# Rename columns
itemid_names <- c(
  "50912" = "creat_baseline",
  "51301" = "wbc",
  "51755" = "wbc_alt",
  "51265" = "platelets",
  "50885" = "bilirubin",
  "50813" = "lactate",
  "50931" = "glucose",
  "50882" = "bicarbonate",
  "50971" = "potassium",
  "50983" = "sodium",
  "51006" = "bun",
  "50893" = "calcium"
)

# Rename existing columns
existing_cols <- intersect(names(itemid_names), names(labs_wide))
setnames(labs_wide, old=existing_cols, new=itemid_names[existing_cols])

# Combine WBC
if("wbc" %in% names(labs_wide) && "wbc_alt" %in% names(labs_wide)) {
  labs_wide[is.na(wbc) & !is.na(wbc_alt), wbc := wbc_alt]
  labs_wide[, wbc_alt := NULL]
} else if("wbc_alt" %in% names(labs_wide)) {
  setnames(labs_wide, "wbc_alt", "wbc")
}

cat("Lab variables available:\n")
print(names(labs_wide))
cat("\nLab missingness:\n")
for(v in names(labs_wide)[-(1:2)]) {
  cat(sprintf("  %s: missing %d/%d (%.1f%%)\n",
              v, sum(is.na(labs_wide[[v]])), nrow(labs_wide),
              100*mean(is.na(labs_wide[[v]]))))
}

# Merge into covars
covars <- merge(covars, labs_wide, by=c("subject_id","hadm_id"), all.x=TRUE)

###############################################################################
# 3. Compute eGFR (CKD-EPI 2021)
###############################################################################
cat("\n--- Step 3: Compute eGFR ---\n")

# CKD-EPI 2021 (race-free)
# eGFR = 142 * min(Scr/kappa, 1)^alpha * max(Scr/kappa, 1)^(-1.200) * 0.9938^Age * (1.012 if female)
# Female: kappa=0.7, alpha=-0.241; Male: kappa=0.9, alpha=-0.302

covars[, egfr_baseline := {
  kappa <- fifelse(female==1, 0.7, 0.9)
  alpha <- fifelse(female==1, -0.241, -0.302)
  cr <- creat_baseline / kappa
  egfr_val <- 142 * pmin(cr, 1)^alpha * pmax(cr, 1)^(-1.200) * 0.9938^age * (1 + 0.012*female)
  round(egfr_val, 1)
}]

cat("eGFR baseline summary:\n")
print(summary(covars$egfr_baseline))

# CKD flag from eGFR
covars[, ckd_egfr := as.integer(!is.na(egfr_baseline) & egfr_baseline < 60)]

###############################################################################
# 4. Comorbidities from diagnoses_icd
###############################################################################
cat("\n--- Step 4: Comorbidities from ICD codes ---\n")

con <- open_duck()
diag_path <- file.path(HOSP_DIR, "diagnoses_icd.csv.gz")

diag_query <- paste0("
  SELECT subject_id, hadm_id, icd_code, icd_version
  FROM read_csv_auto('", diag_path, "', compression='gzip')
  WHERE hadm_id IN (", hadm_list, ")
")

diag <- dbGetQuery(con, diag_query)
dbDisconnect(con, shutdown=TRUE)

cat(sprintf("Diagnoses loaded: %d rows\n", nrow(diag)))

diag <- as.data.table(diag)
diag[, icd9 := icd_version == 9]
diag[, icd10 := icd_version == 10]

# ICD code matching function
has_icd <- function(dt, subject_id_vec, hadm_id_vec, icd10_patterns=NULL, icd9_patterns=NULL) {
  # dt has subject_id, hadm_id, icd_code, icd9, icd10
  matches <- rep(FALSE, length(subject_id_vec))

  key_dt <- data.table(subject_id=subject_id_vec, hadm_id=hadm_id_vec, idx=seq_along(subject_id_vec))
  merged_dt <- merge(key_dt, dt, by=c("subject_id","hadm_id"))

  if(!is.null(icd10_patterns) && nrow(merged_dt[icd10==TRUE]) > 0) {
    icd10_match <- merged_dt[icd10==TRUE][, any(grepl(paste(icd10_patterns, collapse="|"), icd_code)), by=idx]
    matches[icd10_match[V1==TRUE, idx]] <- TRUE
  }
  if(!is.null(icd9_patterns) && nrow(merged_dt[icd9==TRUE]) > 0) {
    icd9_match <- merged_dt[icd9==TRUE][, any(grepl(paste(icd9_patterns, collapse="|"), icd_code)), by=idx]
    matches[icd9_match[V1==TRUE, idx]] <- TRUE
  }
  as.integer(matches)
}

# CKD: N18 (ICD10), 585 (ICD9)
covars[, ckd_icd := has_icd(diag, subject_id, hadm_id,
                              icd10_patterns="^N18",
                              icd9_patterns="^585")]
# Diabetes: E10-E14 (ICD10), 250 (ICD9)
covars[, diabetes := has_icd(diag, subject_id, hadm_id,
                               icd10_patterns="^E1[0-4]",
                               icd9_patterns="^250")]
# Hypertension: I10-I15 (ICD10), 401-405 (ICD9)
covars[, hypertension := has_icd(diag, subject_id, hadm_id,
                                   icd10_patterns="^I1[0-5]",
                                   icd9_patterns="^40[1-5]")]
# Heart failure: I50 (ICD10), 428 (ICD9)
covars[, heart_failure := has_icd(diag, subject_id, hadm_id,
                                    icd10_patterns="^I50",
                                    icd9_patterns="^428")]
# Liver disease: K70-K77 (ICD10), 571 (ICD9)
covars[, liver_disease := has_icd(diag, subject_id, hadm_id,
                                    icd10_patterns="^K7[0-7]",
                                    icd9_patterns="^571")]
# Cancer: C00-C99 (ICD10), 140-209 (ICD9)
covars[, cancer := has_icd(diag, subject_id, hadm_id,
                             icd10_patterns="^C[0-9]",
                             icd9_patterns="^1[4-9][0-9]|^20[0-9]")]
# COPD: J44 (ICD10), 496 (ICD9)
covars[, copd := has_icd(diag, subject_id, hadm_id,
                          icd10_patterns="^J44",
                          icd9_patterns="^496")]
# HIV/AIDS: B20-B24, Z21 (ICD10), 042-044 (ICD9)
covars[, hiv := has_icd(diag, subject_id, hadm_id,
                         icd10_patterns="^B2[0-4]|^Z21",
                         icd9_patterns="^04[2-4]")]

cat("Comorbidity prevalence:\n")
for(v in c("ckd_icd","diabetes","hypertension","heart_failure","liver_disease","cancer","copd","hiv")) {
  cat(sprintf("  %s: %d (%.1f%%)\n", v, sum(covars[[v]]==1, na.rm=TRUE),
              100*mean(covars[[v]]==1, na.rm=TRUE)))
}

# Charlson comorbidity index (simplified)
covars[, charlson_raw := (
  ckd_icd * 2 +
  diabetes * 1 +
  heart_failure * 1 +
  liver_disease * 3 +  # severe liver disease
  cancer * 2 +
  copd * 1 +
  hiv * 6
)]

###############################################################################
# 5. ICU status
###############################################################################
cat("\n--- Step 5: ICU status ---\n")

con <- open_duck()
icu_path <- file.path(ICU_DIR, "icustays.csv.gz")

icu <- dbGetQuery(con, paste0("
  SELECT subject_id, hadm_id, stay_id,
         intime, outtime, first_careunit, last_careunit, los
  FROM read_csv_auto('", icu_path, "', compression='gzip')
  WHERE hadm_id IN (", hadm_list, ")
"))
dbDisconnect(con, shutdown=TRUE)

icu <- as.data.table(icu)
icu[, intime := as.POSIXct(intime)]
icu[, outtime := as.POSIXct(outtime)]

cat(sprintf("ICU stays loaded: %d\n", nrow(icu)))

# ICU at t0: in ICU when t0 occurs
icu_status <- merge(icu, cohort[, .(subject_id, hadm_id, t0)],
                    by=c("subject_id","hadm_id"))
icu_status[, in_icu_at_t0 := (intime <= t0 & (is.na(outtime) | outtime >= t0))]
icu_status[, icu_before_t0 := intime <= t0]

icu_per_hadm <- icu_status[, .(
  icu_at_t0 = as.integer(any(in_icu_at_t0)),
  icu_before_t0 = as.integer(any(icu_before_t0)),
  icu_los_before_t0 = sum(pmin(as.numeric(difftime(pmin(outtime, t0[1]), intime, units="hours")),
                                as.numeric(difftime(t0[1], intime, units="hours")), na.rm=TRUE), na.rm=TRUE)
), by=.(subject_id, hadm_id)]

covars <- merge(covars, icu_per_hadm, by=c("subject_id","hadm_id"), all.x=TRUE)
covars[is.na(icu_at_t0), icu_at_t0 := 0L]
covars[is.na(icu_before_t0), icu_before_t0 := 0L]

cat(sprintf("ICU at t0: %d (%.1f%%)\n",
            sum(covars$icu_at_t0==1), 100*mean(covars$icu_at_t0==1)))

###############################################################################
# 6. Vasopressors before/at t0
###############################################################################
cat("\n--- Step 6: Vasopressors ---\n")

# Use prescriptions for simplicity (emar would be more accurate but slow)
con <- open_duck()
rx_path <- file.path(HOSP_DIR, "prescriptions.csv.gz")

vaso_query <- paste0("
  SELECT subject_id, hadm_id, starttime, stoptime, drug
  FROM read_csv_auto('", rx_path, "', compression='gzip')
  WHERE hadm_id IN (", hadm_list, ")
    AND (
      UPPER(drug) LIKE '%NOREPINEPHRINE%'
      OR UPPER(drug) LIKE '%VASOPRESSIN%'
      OR UPPER(drug) LIKE '%EPINEPHRINE%'
      OR UPPER(drug) LIKE '%PHENYLEPHRINE%'
      OR UPPER(drug) LIKE '%DOPAMINE%'
    )
")

vaso_rx <- dbGetQuery(con, vaso_query)
dbDisconnect(con, shutdown=TRUE)

vaso_rx <- as.data.table(vaso_rx)
vaso_rx[, starttime := as.POSIXct(starttime)]
vaso_rx[, stoptime := as.POSIXct(stoptime)]

cat(sprintf("Vasopressor rows: %d\n", nrow(vaso_rx)))

# Vasopressor before or at t0
vaso_with_t0 <- merge(vaso_rx, cohort[, .(subject_id, hadm_id, t0)],
                      by=c("subject_id","hadm_id"))
vaso_with_t0[, vaso_before_t0 := starttime <= t0]

vaso_per_hadm <- vaso_with_t0[vaso_before_t0==TRUE, .(
  vasopressor = 1L
), by=.(subject_id, hadm_id)]
vaso_per_hadm <- unique(vaso_per_hadm)

covars <- merge(covars, vaso_per_hadm, by=c("subject_id","hadm_id"), all.x=TRUE)
covars[is.na(vasopressor), vasopressor := 0L]

cat(sprintf("Vasopressor before t0: %d (%.1f%%)\n",
            sum(covars$vasopressor==1), 100*mean(covars$vasopressor==1)))

###############################################################################
# 7. Mechanical ventilation
###############################################################################
cat("\n--- Step 7: Mechanical ventilation ---\n")

con <- open_duck()
proc_path <- file.path(HOSP_DIR, "procedures_icd.csv.gz")

# ICD procedures for mechanical ventilation
# ICD10: 5A19XXX series, or ICD9: 96.7x
vent_query <- paste0("
  SELECT subject_id, hadm_id, icd_code, icd_version
  FROM read_csv_auto('", proc_path, "', compression='gzip')
  WHERE hadm_id IN (", hadm_list, ")
    AND (
      (icd_version = 10 AND icd_code LIKE '5A19%')
      OR (icd_version = 9 AND icd_code LIKE '967%')
      OR (icd_version = 9 AND icd_code LIKE '9670%')
      OR (icd_version = 9 AND icd_code LIKE '9671%')
      OR (icd_version = 9 AND icd_code LIKE '9672%')
    )
")

vent_proc <- dbGetQuery(con, vent_query)
dbDisconnect(con, shutdown=TRUE)

cat(sprintf("Ventilation procedure rows: %d\n", nrow(vent_proc)))

if(nrow(vent_proc) > 0) {
  vent_proc <- as.data.table(vent_proc)
  vent_per_hadm <- vent_proc[, .(mech_vent = 1L), by=.(subject_id, hadm_id)]
  vent_per_hadm <- unique(vent_per_hadm)
  covars <- merge(covars, vent_per_hadm, by=c("subject_id","hadm_id"), all.x=TRUE)
}
if(!"mech_vent" %in% names(covars)) covars[, mech_vent := 0L]
covars[is.na(mech_vent), mech_vent := 0L]

cat(sprintf("Mechanical ventilation in admission: %d (%.1f%%)\n",
            sum(covars$mech_vent==1), 100*mean(covars$mech_vent==1)))

###############################################################################
# 8. Nephrotoxic co-exposures before/at t0
###############################################################################
cat("\n--- Step 8: Nephrotoxic co-exposures ---\n")

con <- open_duck()

nephrotox_query <- paste0("
  SELECT subject_id, hadm_id, starttime, stoptime, drug
  FROM read_csv_auto('", rx_path, "', compression='gzip')
  WHERE hadm_id IN (", hadm_list, ")
    AND (
      UPPER(drug) LIKE '%VANCOMYCIN%'
      OR UPPER(drug) LIKE '%GENTAMICIN%'
      OR UPPER(drug) LIKE '%TOBRAMYCIN%'
      OR UPPER(drug) LIKE '%AMIKACIN%'
      OR UPPER(drug) LIKE '%FUROSEMIDE%'
      OR UPPER(drug) LIKE '%TORSEMIDE%'
      OR UPPER(drug) LIKE '%BUMETANIDE%'
      OR UPPER(drug) LIKE '%COLISTIN%'
      OR UPPER(drug) LIKE '%POLYMYXIN%'
      OR UPPER(drug) LIKE '%AMPHOTERICIN%'
      OR UPPER(drug) LIKE '%CISPLATIN%'
      OR UPPER(drug) LIKE '%CYCLOSPORINE%'
      OR UPPER(drug) LIKE '%TACROLIMUS%'
    )
")

nephrotox_rx <- dbGetQuery(con, nephrotox_query)
dbDisconnect(con, shutdown=TRUE)

nephrotox_rx <- as.data.table(nephrotox_rx)
nephrotox_rx[, starttime := as.POSIXct(starttime)]

cat(sprintf("Nephrotoxic drug rows: %d\n", nrow(nephrotox_rx)))
cat("Nephrotoxic drug distribution:\n")
print(sort(table(toupper(nephrotox_rx$drug)), decreasing=TRUE)[1:10])

# Classify
nephrotox_rx[, drug_class := fcase(
  grepl("VANCOMYCIN", toupper(drug)), "vancomycin",
  grepl("GENTAMICIN|TOBRAMYCIN|AMIKACIN", toupper(drug)), "aminoglycoside",
  grepl("FUROSEMIDE|TORSEMIDE|BUMETANIDE", toupper(drug)), "loop_diuretic",
  grepl("COLISTIN|POLYMYXIN", toupper(drug)), "colistin_polymyxin",
  grepl("AMPHOTERICIN", toupper(drug)), "amphotericin",
  default = "other_nephrotox"
)]

# Before/at t0
nephrotox_with_t0 <- merge(nephrotox_rx, cohort[, .(subject_id, hadm_id, t0)],
                            by=c("subject_id","hadm_id"))
nephrotox_before <- nephrotox_with_t0[starttime <= t0]

# Per-drug-class indicator
nephrotox_per_hadm <- nephrotox_before[, .(
  concomitant_vancomycin = as.integer(any(drug_class=="vancomycin")),
  concomitant_aminoglycoside = as.integer(any(drug_class=="aminoglycoside")),
  concomitant_loop_diuretic = as.integer(any(drug_class=="loop_diuretic")),
  concomitant_nephrotox_any = 1L
), by=.(subject_id, hadm_id)]
nephrotox_per_hadm <- unique(nephrotox_per_hadm)

covars <- merge(covars, nephrotox_per_hadm, by=c("subject_id","hadm_id"), all.x=TRUE)
for(v in c("concomitant_vancomycin","concomitant_aminoglycoside",
           "concomitant_loop_diuretic","concomitant_nephrotox_any")) {
  covars[is.na(get(v)), (v) := 0L]
}

cat("Nephrotoxic co-exposures before t0:\n")
for(v in c("concomitant_vancomycin","concomitant_aminoglycoside",
           "concomitant_loop_diuretic","concomitant_nephrotox_any")) {
  cat(sprintf("  %s: %d (%.1f%%)\n", v,
              sum(covars[[v]]==1, na.rm=TRUE),
              100*mean(covars[[v]]==1, na.rm=TRUE)))
}

###############################################################################
# 9. Final cleaning
###############################################################################
cat("\n--- Step 9: Final covariate dataset ---\n")

# Combined CKD (ICD or eGFR)
covars[, ckd_any := as.integer(ckd_icd==1 | (!is.na(ckd_egfr) & ckd_egfr==1))]
covars[is.na(ckd_any), ckd_any := 0L]

# Severity: Charlson updated
covars[, charlson := charlson_raw + (ckd_icd * 2)]

# Select final columns
cov_cols <- c(
  "subject_id", "hadm_id", "t0",
  # Demographics
  "age", "female", "race_cat",
  # Labs
  "creat_baseline", "egfr_baseline", "wbc", "platelets",
  "bilirubin", "lactate", "glucose", "bicarbonate", "sodium", "potassium", "bun",
  # Comorbidities
  "ckd_icd", "ckd_egfr", "ckd_any",
  "diabetes", "hypertension", "heart_failure", "liver_disease",
  "cancer", "copd", "hiv", "charlson",
  # Clinical state
  "icu_at_t0", "icu_before_t0",
  "vasopressor", "mech_vent",
  # Co-exposures
  "concomitant_vancomycin", "concomitant_aminoglycoside",
  "concomitant_loop_diuretic", "concomitant_nephrotox_any",
  # Time variables
  "hours_bc_to_t0", "admission_type"
)

# Keep only existing columns
cov_cols_exist <- intersect(cov_cols, names(covars))
covariates_out <- covars[, ..cov_cols_exist]

cat(sprintf("Covariates: %d rows, %d variables\n", nrow(covariates_out), ncol(covariates_out)))

# Missingness report
cat("\nMissingness report:\n")
miss_df <- data.frame(
  variable = names(covariates_out),
  n_missing = sapply(covariates_out, function(x) sum(is.na(x))),
  pct_missing = round(sapply(covariates_out, function(x) 100*mean(is.na(x))), 1)
)
print(miss_df[miss_df$n_missing > 0, ])

###############################################################################
# 10. Save
###############################################################################
saveRDS(covariates_out, file.path(DERIVED_DIR, "covariates.rds"))
cat(sprintf("\nSaved: derived/covariates.rds\n"))

# Covariate dictionary
dict <- data.frame(
  variable = names(covariates_out),
  description = c(
    "Subject ID", "Admission ID", "Time zero",
    "Age at admission (years)", "Female sex (0/1)", "Race category",
    "Baseline creatinine (mg/dL)", "Baseline eGFR (mL/min/1.73m2)",
    "WBC (K/uL)", "Platelets (K/uL)",
    "Total bilirubin (mg/dL)", "Lactate (mmol/L)", "Glucose (mg/dL)",
    "Bicarbonate (mEq/L)", "Sodium (mEq/L)", "Potassium (mEq/L)", "BUN (mg/dL)",
    "CKD by ICD code", "CKD by eGFR<60", "CKD any definition",
    "Diabetes mellitus", "Hypertension", "Heart failure",
    "Liver disease", "Cancer", "COPD", "HIV/AIDS",
    "Charlson comorbidity index (simplified)",
    "In ICU at t0", "ICU admission before t0",
    "Vasopressor before t0", "Mechanical ventilation in admission",
    "Concomitant vancomycin", "Concomitant aminoglycoside",
    "Concomitant loop diuretic", "Any nephrotoxic co-exposure",
    "Hours from blood culture to t0", "Admission type"
  )[seq_along(names(covariates_out))]
)

fwrite(dict, file.path(TABLES_DIR, "covariate_dictionary.csv"))
cat("Saved: tables/covariate_dictionary.csv\n")

cat("\n=== 04_build_covariates.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
