###############################################################################
# 03_build_outcomes_aki.R
# Build AKI outcomes using KDIGO creatinine criteria
# Output: derived/outcomes.rds
###############################################################################

library(duckdb)
library(DBI)
library(dplyr)
library(lubridate)
library(data.table)

cat("=== 03_build_outcomes_aki.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

WORK_DIR <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
HOSP_DIR <- file.path(WORK_DIR, "mimic-iv-3.1/hosp")
DERIVED_DIR <- file.path(WORK_DIR, "derived")

# Load cohort
cohort <- readRDS(file.path(DERIVED_DIR, "cohort_treated.rds"))
cat(sprintf("Cohort loaded: %d patients\n", nrow(cohort)))

subj_ids <- unique(cohort$subject_id)
sid_list <- paste(subj_ids, collapse=",")

###############################################################################
# 1. Load creatinine from labevents (itemid=50912)
###############################################################################
cat("\n--- Step 1: Load creatinine from labevents ---\n")

con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
dbExecute(con, "SET threads=4; SET memory_limit='8GB';")

lab_path <- file.path(HOSP_DIR, "labevents.csv.gz")

# Peek columns first
lab_cols <- dbGetQuery(con, paste0(
  "SELECT * FROM read_csv_auto('", lab_path, "', compression='gzip') LIMIT 0"
))
cat("Labevents columns:", paste(names(lab_cols), collapse=", "), "\n")

# Load creatinine only (itemid 50912) for our subjects
creat_query <- paste0("
  SELECT subject_id, hadm_id, itemid, charttime, valuenum, valueuom, flag
  FROM read_csv_auto('", lab_path, "', compression='gzip')
  WHERE itemid = 50912
    AND subject_id IN (", sid_list, ")
    AND valuenum IS NOT NULL
    AND valuenum > 0
    AND valuenum < 50
")

creat_raw <- dbGetQuery(con, creat_query)
dbDisconnect(con, shutdown=TRUE)

cat(sprintf("Creatinine rows loaded: %d\n", nrow(creat_raw)))
cat(sprintf("Unique subjects: %d\n", length(unique(creat_raw$subject_id))))

creat_raw <- as.data.table(creat_raw)
creat_raw[, charttime := as.POSIXct(charttime)]

###############################################################################
# 2. For each patient, compute baseline and peak creatinine
###############################################################################
cat("\n--- Step 2: Compute baseline and peak creatinine ---\n")

# Merge creatinine with cohort to get t0
creat_with_t0 <- merge(
  creat_raw,
  cohort[, .(subject_id, hadm_id, t0, dischtime, hospital_expire_flag, dod)],
  by=c("subject_id","hadm_id"),
  all.x=FALSE
)

cat(sprintf("Creatinine rows linked to cohort: %d\n", nrow(creat_with_t0)))

# Compute hours from t0
creat_with_t0[, hours_from_t0 := as.numeric(difftime(charttime, t0, units="hours"))]

# BASELINE: creatinine in [-48h, +6h] relative to t0
# Prefer measurements closest to and before t0
baseline_window <- creat_with_t0[hours_from_t0 >= -48 & hours_from_t0 <= 6]

# Best baseline = measurement closest to t0 (prefer pre-t0 if available)
# Strategy: take lowest value in window (KDIGO reference)
baseline <- baseline_window[, .(
  baseline_creat = min(valuenum, na.rm=TRUE),
  baseline_creat_median = median(valuenum, na.rm=TRUE),
  baseline_creat_time = charttime[which.min(valuenum)],
  n_baseline_measures = .N
), by=.(subject_id, hadm_id)]

cat(sprintf("Patients with baseline creatinine: %d\n", nrow(baseline)))

# PEAK: maximum creatinine in [t0, t0+7days]
peak_window <- creat_with_t0[hours_from_t0 >= 0 & hours_from_t0 <= 168]

peak_creat <- peak_window[, .(
  peak_creat = max(valuenum, na.rm=TRUE),
  peak_creat_time = charttime[which.max(valuenum)],
  n_peak_measures = .N
), by=.(subject_id, hadm_id)]

cat(sprintf("Patients with peak creatinine in 7d: %d\n", nrow(peak_creat)))

# Also get creatinine at day 2 (48h) for trajectory
day2_window <- creat_with_t0[hours_from_t0 >= 24 & hours_from_t0 <= 72]
creat_day2 <- day2_window[, .(
  creat_day2 = median(valuenum, na.rm=TRUE),
  n_day2_measures = .N
), by=.(subject_id, hadm_id)]

###############################################################################
# 3. Compute AKI_7d using KDIGO criteria
###############################################################################
cat("\n--- Step 3: Compute AKI outcomes ---\n")

# Combine baseline and peak
aki_data <- merge(baseline, peak_creat, by=c("subject_id","hadm_id"), all=TRUE)
aki_data <- merge(aki_data, creat_day2, by=c("subject_id","hadm_id"), all=TRUE)

# KDIGO AKI criteria:
# 1. Increase >= 0.3 mg/dL within 48h (same-day data)
# 2. Increase to >= 1.5x baseline within 7 days
aki_data[, delta_creat_7d := peak_creat - baseline_creat]
aki_data[, creat_ratio := peak_creat / baseline_creat]

aki_data[, AKI_7d := fifelse(
  !is.na(baseline_creat) & !is.na(peak_creat) &
    (delta_creat_7d >= 0.3 | creat_ratio >= 1.5),
  1L, 0L
)]

# Stage 1: >=0.3 or >=1.5x
# Stage 2: >=2x
# Stage 3: >=3x or >=4.0 mg/dL with acute rise, or RRT
aki_data[, AKI_stage := fifelse(
  is.na(baseline_creat) | is.na(peak_creat), NA_integer_,
  fifelse(creat_ratio >= 3 | (creat_ratio >= 2.5 & peak_creat >= 4.0), 3L,
  fifelse(creat_ratio >= 2, 2L,
  fifelse(delta_creat_7d >= 0.3 | creat_ratio >= 1.5, 1L, 0L)))
)]

cat("AKI_7d distribution:\n")
print(table(aki_data$AKI_7d, useNA="ifany"))
cat("\nAKI stage distribution:\n")
print(table(aki_data$AKI_stage, useNA="ifany"))

###############################################################################
# 4. Mortality outcomes
###############################################################################
cat("\n--- Step 4: Mortality outcomes ---\n")

# Join with cohort for mortality
outcomes_full <- merge(
  aki_data,
  cohort[, .(subject_id, hadm_id, t0, dischtime, hospital_expire_flag, dod, W)],
  by=c("subject_id","hadm_id"),
  all.y=TRUE
)

# In-hospital mortality
outcomes_full[, died_inhosp := as.integer(hospital_expire_flag == 1)]

# 7-day mortality
outcomes_full[, died_7d := as.integer(
  (!is.na(dod) & as.numeric(difftime(dod, t0, units="days")) <= 7) |
  (hospital_expire_flag == 1 & !is.na(t0) & !is.na(dischtime) &
     as.numeric(difftime(dischtime, t0, units="days")) <= 7)
)]
outcomes_full[is.na(died_7d), died_7d := 0L]

# Composite: AKI_7d OR died_7d
outcomes_full[, composite_aki_death := as.integer(
  (!is.na(AKI_7d) & AKI_7d == 1) | died_7d == 1
)]

cat("Died in hospital:\n")
print(table(outcomes_full$died_inhosp, useNA="ifany"))
cat("\nDied within 7d:\n")
print(table(outcomes_full$died_7d, useNA="ifany"))
cat("\nComposite AKI/death:\n")
print(table(outcomes_full$composite_aki_death, useNA="ifany"))

###############################################################################
# 5. Baseline AKI (at t0)
###############################################################################
cat("\n--- Step 5: Baseline AKI ---\n")

# AKI at baseline: using pre-admission creatinine
# Look for creatinine in [-7 days, -48h] before t0 as "pre-admission"
preadm_window <- creat_with_t0[hours_from_t0 >= -168 & hours_from_t0 < -48]
preadm <- preadm_window[, .(
  preadm_creat = min(valuenum, na.rm=TRUE),
  n_preadm = .N
), by=.(subject_id, hadm_id)]

outcomes_full <- merge(outcomes_full, preadm, by=c("subject_id","hadm_id"), all.x=TRUE)

# AKI at baseline (relative to pre-admission)
outcomes_full[, aki_at_baseline := fifelse(
  !is.na(preadm_creat) & !is.na(baseline_creat) &
    (baseline_creat / preadm_creat >= 1.5 | baseline_creat - preadm_creat >= 0.3),
  1L, 0L
)]

cat("AKI at baseline distribution:\n")
print(table(outcomes_full$aki_at_baseline, useNA="ifany"))

###############################################################################
# 6. Summary
###############################################################################
cat("\n--- Outcomes summary by treatment group ---\n")

# By W group
summary_tbl <- outcomes_full[!is.na(W), .(
  N = .N,
  AKI_7d_n = sum(AKI_7d==1, na.rm=TRUE),
  AKI_7d_pct = round(mean(AKI_7d==1, na.rm=TRUE)*100, 1),
  Died_7d_n = sum(died_7d==1, na.rm=TRUE),
  Composite_n = sum(composite_aki_death==1, na.rm=TRUE),
  Missing_creat = sum(is.na(AKI_7d))
), by=W]

print(summary_tbl)

###############################################################################
# 7. Save
###############################################################################
outcomes_out <- outcomes_full[, .(
  subject_id, hadm_id, W, t0,
  baseline_creat, peak_creat, creat_day2,
  delta_creat_7d, creat_ratio,
  AKI_7d, AKI_stage,
  n_baseline_measures, n_peak_measures,
  aki_at_baseline, preadm_creat,
  died_inhosp, died_7d, composite_aki_death
)]

saveRDS(outcomes_out, file.path(DERIVED_DIR, "outcomes.rds"))
cat(sprintf("\nSaved: derived/outcomes.rds (%d rows)\n", nrow(outcomes_out)))

cat("\n=== 03_build_outcomes_aki.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
