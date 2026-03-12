###############################################################################
# 01_build_mssa_cohort.R
# Build the MSSA bacteremia candidate cohort
# Output: derived/mssa_candidates.rds
###############################################################################

library(duckdb)
library(DBI)
library(dplyr)
library(lubridate)
library(data.table)
library(stringr)

cat("=== 01_build_mssa_cohort.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

WORK_DIR <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
HOSP_DIR <- file.path(WORK_DIR, "mimic-iv-3.1/hosp")
DERIVED_DIR <- file.path(WORK_DIR, "derived")

# Helper: open DuckDB connection
open_duck <- function() {
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  dbExecute(con, "SET threads=4;")
  dbExecute(con, "SET memory_limit='8GB';")
  con
}

###############################################################################
# 1. Load microbiologyevents: SA blood cultures
###############################################################################
cat("--- Step 1: Load microbiologyevents ---\n")
con <- open_duck()

micro_path <- file.path(HOSP_DIR, "microbiologyevents.csv.gz")

# First peek at columns
cols_query <- sprintf(
  "SELECT * FROM read_csv_auto('%s', compression='gzip') LIMIT 0",
  micro_path
)
cols_df <- dbGetQuery(con, cols_query)
cat("Microbiologyevents columns:", paste(names(cols_df), collapse=", "), "\n")

# Get SA blood cultures
sa_query <- paste0("
  SELECT subject_id, hadm_id, micro_specimen_id,
         chartdate, charttime,
         spec_type_desc, org_name,
         ab_name, interpretation,
         test_name
  FROM read_csv_auto('", micro_path, "', compression='gzip')
  WHERE spec_type_desc = 'BLOOD CULTURE'
    AND (UPPER(org_name) LIKE '%STAPH AUREUS%' OR UPPER(org_name) LIKE '%STAPHYLOCOCCUS AUREUS%')
")

sa_blood <- dbGetQuery(con, sa_query)
cat(sprintf("SA blood culture rows: %d\n", nrow(sa_blood)))
cat(sprintf("Unique patients: %d\n", length(unique(sa_blood$subject_id))))
cat(sprintf("Unique admissions: %d\n", length(unique(sa_blood$hadm_id))))

dbDisconnect(con, shutdown=TRUE)

# Convert to data.table for processing
sa_blood <- as.data.table(sa_blood)
sa_blood[, charttime := as.POSIXct(charttime)]
sa_blood[, chartdate := as.Date(chartdate)]

###############################################################################
# 2. Identify MSSA based on oxacillin susceptibility
###############################################################################
cat("\n--- Step 2: Determine MSSA status ---\n")

# Get rows with oxacillin susceptibility
ox_susc <- sa_blood[!is.na(ab_name) & toupper(ab_name) == "OXACILLIN",
                    .(subject_id, hadm_id, micro_specimen_id, interpretation)]

cat(sprintf("Rows with oxacillin susceptibility: %d\n", nrow(ox_susc)))
cat("Interpretation distribution:\n")
print(table(ox_susc$interpretation, useNA="ifany"))

# Per patient/hadm: if ANY specimen has R → MRSA; if all S → MSSA
ox_per_hadm <- ox_susc[, .(
  has_S = any(interpretation == "S", na.rm=TRUE),
  has_R = any(interpretation == "R", na.rm=TRUE)
), by=.(subject_id, hadm_id)]

ox_per_hadm[, mssa_confirmed := fifelse(has_R, FALSE, fifelse(has_S, TRUE, NA))]

cat(sprintf("\nHadm with MSSA confirmed: %d\n", sum(ox_per_hadm$mssa_confirmed == TRUE, na.rm=TRUE)))
cat(sprintf("Hadm with MRSA confirmed: %d\n", sum(ox_per_hadm$mssa_confirmed == FALSE, na.rm=TRUE)))

###############################################################################
# 3. Get first blood culture time per patient/hadm
###############################################################################
cat("\n--- Step 3: First blood culture time ---\n")

# Get unique specimen events (not susceptibility rows - those duplicate)
# Use distinct specimen per patient/hadm
unique_cultures <- sa_blood[, .(
  first_culture_time = min(charttime, na.rm=TRUE),
  first_culture_date = min(chartdate, na.rm=TRUE),
  n_cultures = uniqueN(micro_specimen_id)
), by=.(subject_id, hadm_id)]

cat(sprintf("Unique patient/hadm pairs with SA blood culture: %d\n", nrow(unique_cultures)))

# Some charttime may be NA - fill with chartdate at noon if needed
unique_cultures[is.na(first_culture_time) & !is.na(first_culture_date),
                first_culture_time := as.POSIXct(paste(first_culture_date, "12:00:00"))]

###############################################################################
# 4. Load admissions and patients
###############################################################################
cat("\n--- Step 4: Load admissions and patients ---\n")
con <- open_duck()

adm_path <- file.path(HOSP_DIR, "admissions.csv.gz")
pat_path <- file.path(HOSP_DIR, "patients.csv.gz")

adm <- dbGetQuery(con, sprintf(
  "SELECT subject_id, hadm_id, admittime, dischtime, race, admission_type, hospital_expire_flag
   FROM read_csv_auto('%s', compression='gzip')", adm_path))

pat <- dbGetQuery(con, sprintf(
  "SELECT subject_id, gender, anchor_age, anchor_year, dod
   FROM read_csv_auto('%s', compression='gzip')", pat_path))

dbDisconnect(con, shutdown=TRUE)

adm <- as.data.table(adm)
pat <- as.data.table(pat)

adm[, admittime := as.POSIXct(admittime)]
adm[, dischtime := as.POSIXct(dischtime)]
pat[, dod := as.Date(dod)]

cat(sprintf("Admissions loaded: %d\n", nrow(adm)))
cat(sprintf("Patients loaded: %d\n", nrow(pat)))

###############################################################################
# 5. Build candidate table
###############################################################################
cat("\n--- Step 5: Build candidate table ---\n")

# Merge: unique cultures + admissions + patients
cands <- merge(unique_cultures, adm[, .(subject_id, hadm_id, admittime, dischtime,
                                         race, admission_type, hospital_expire_flag)],
               by=c("subject_id","hadm_id"), all.x=TRUE)

cands <- merge(cands, pat[, .(subject_id, gender, anchor_age, anchor_year, dod)],
               by="subject_id", all.x=TRUE)

# For admissions without hadm_id match (hadm_id is NA), we keep but flag
cat(sprintf("Candidates before eligibility filtering: %d\n", nrow(cands)))
cat(sprintf("With missing hadm_id: %d\n", sum(is.na(cands$hadm_id))))

# Compute approximate age at admission
cands[!is.na(admittime) & !is.na(anchor_year),
      age_at_adm := anchor_age + (year(admittime) - anchor_year)]
# If admittime missing, use anchor_age
cands[is.na(age_at_adm), age_at_adm := anchor_age]

###############################################################################
# 6. Add MSSA classification
###############################################################################
cands <- merge(cands, ox_per_hadm[, .(subject_id, hadm_id, mssa_confirmed,
                                       has_S, has_R)],
               by=c("subject_id","hadm_id"), all.x=TRUE)

# Candidates without oxacillin result
cands[is.na(mssa_confirmed), mssa_confirmed := NA]  # unknown

cat(sprintf("\nMSSA confirmed (oxacillin S): %d\n",
            sum(cands$mssa_confirmed == TRUE, na.rm=TRUE)))
cat(sprintf("MRSA confirmed (oxacillin R): %d\n",
            sum(cands$mssa_confirmed == FALSE, na.rm=TRUE)))
cat(sprintf("Unknown susceptibility: %d\n",
            sum(is.na(cands$mssa_confirmed))))

###############################################################################
# 7. Eligibility filtering
###############################################################################
cat("\n--- Step 6: Eligibility filtering ---\n")

n0 <- nrow(cands)
cat(sprintf("Starting N: %d\n", n0))

# Age >= 18
cands_elig <- cands[!is.na(age_at_adm) & age_at_adm >= 18]
cat(sprintf("After age >= 18: %d (removed %d)\n", nrow(cands_elig), n0 - nrow(cands_elig)))

# Must have hadm_id
n1 <- nrow(cands_elig)
cands_elig <- cands_elig[!is.na(hadm_id)]
cat(sprintf("After requiring hadm_id: %d (removed %d)\n", nrow(cands_elig), n1 - nrow(cands_elig)))

# Exclude MRSA confirmed (we keep MSSA + unknown)
n2 <- nrow(cands_elig)
cands_elig <- cands_elig[is.na(mssa_confirmed) | mssa_confirmed == TRUE]
cat(sprintf("After excluding MRSA: %d (removed %d)\n", nrow(cands_elig), n2 - nrow(cands_elig)))

cat(sprintf("\nFinal candidate pool: %d patients, %d admissions\n",
            length(unique(cands_elig$subject_id)),
            nrow(cands_elig)))

cat(sprintf("  MSSA confirmed: %d\n", sum(cands_elig$mssa_confirmed == TRUE, na.rm=TRUE)))
cat(sprintf("  Unknown susceptibility: %d\n", sum(is.na(cands_elig$mssa_confirmed))))

###############################################################################
# 8. Clean up and save
###############################################################################
cat("\n--- Step 7: Save output ---\n")

# Standardize columns
mssa_candidates <- cands_elig[, .(
  subject_id, hadm_id,
  mssa_confirmed,
  first_blood_culture_time = first_culture_time,
  first_blood_culture_date = first_culture_date,
  n_cultures,
  admittime, dischtime,
  age = age_at_adm,
  gender,
  race,
  admission_type,
  hospital_expire_flag,
  dod,
  anchor_age, anchor_year
)]

# Sort
setorder(mssa_candidates, subject_id, hadm_id)

cat("Summary of mssa_candidates:\n")
print(str(mssa_candidates))
cat("\nAge distribution:\n")
print(summary(mssa_candidates$age))
cat("\nGender distribution:\n")
print(table(mssa_candidates$gender, useNA="ifany"))
cat("\nRace distribution (top 10):\n")
print(head(sort(table(mssa_candidates$race, useNA="ifany"), decreasing=TRUE), 10))

# Save
saveRDS(mssa_candidates, file.path(DERIVED_DIR, "mssa_candidates.rds"))
cat(sprintf("\nSaved: %s/mssa_candidates.rds\n", DERIVED_DIR))
cat(sprintf("Rows: %d\n", nrow(mssa_candidates)))

cat("\n=== 01_build_mssa_cohort.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
