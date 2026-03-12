###############################################################################
# 02_define_treatment_and_time_zero.R
# Define treatment group and time zero from emar table
# Output: derived/cohort_treated.rds
###############################################################################

library(duckdb)
library(DBI)
library(dplyr)
library(lubridate)
library(data.table)
library(stringr)

cat("=== 02_define_treatment_and_time_zero.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

WORK_DIR <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
HOSP_DIR <- file.path(WORK_DIR, "mimic-iv-3.1/hosp")
DERIVED_DIR <- file.path(WORK_DIR, "derived")

# Load candidates
mssa_cands <- readRDS(file.path(DERIVED_DIR, "mssa_candidates.rds"))
cat(sprintf("MSSA candidates loaded: %d\n", nrow(mssa_cands)))

subj_ids <- unique(mssa_cands$subject_id)
cat(sprintf("Unique subject_ids: %d\n", length(subj_ids)))

###############################################################################
# 1. Load EMAR for ASP and cefazolin administrations
###############################################################################
cat("\n--- Step 1: Load EMAR table ---\n")

con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
dbExecute(con, "SET threads=4; SET memory_limit='8GB';")

emar_path <- file.path(HOSP_DIR, "emar.csv.gz")

# Peek at columns
emar_cols <- dbGetQuery(con, paste0(
  "SELECT * FROM read_csv_auto('", emar_path, "', compression='gzip') LIMIT 0"
))
cat("EMAR columns:", paste(names(emar_cols), collapse=", "), "\n")

# Build subject_id list for SQL IN clause
sid_list <- paste(subj_ids, collapse=",")

# Query ASP and cefazolin from emar filtered to our subjects
emar_query <- paste0("
  SELECT subject_id, hadm_id, emar_id, emar_seq,
         medication, event_txt, charttime,
         scheduletime
  FROM read_csv_auto('", emar_path, "', compression='gzip')
  WHERE subject_id IN (", sid_list, ")
    AND (
      UPPER(medication) LIKE '%NAFCILLIN%'
      OR UPPER(medication) LIKE '%OXACILLIN%'
      OR UPPER(medication) LIKE '%CEFAZOLIN%'
    )
    AND event_txt IN ('Administered', 'Applied', 'Started', 'Restarted', 'New Bag')
")

emar_drugs <- dbGetQuery(con, emar_query)
dbDisconnect(con, shutdown=TRUE)

cat(sprintf("EMAR rows (ASP+CEF for our subjects): %d\n", nrow(emar_drugs)))
cat("Medication distribution:\n")
print(sort(table(emar_drugs$medication, useNA="ifany"), decreasing=TRUE)[1:20])

# Convert to data.table
emar_drugs <- as.data.table(emar_drugs)
emar_drugs[, charttime := as.POSIXct(charttime)]
emar_drugs[, scheduletime := as.POSIXct(scheduletime)]

# Classify drug group
emar_drugs[, drug_group := fifelse(
  grepl("NAFCILLIN|OXACILLIN", toupper(medication)) &
    !grepl("DICLOXACILLIN|CLOXACILLIN", toupper(medication)),
  "ASP",
  fifelse(grepl("CEFAZOLIN", toupper(medication)), "CEF", NA_character_)
)]

# Exclude dicloxacillin/cloxacillin (oral ASPs)
emar_drugs <- emar_drugs[!is.na(drug_group)]

cat(sprintf("\nAfter classification - ASP: %d, CEF: %d\n",
            sum(emar_drugs$drug_group=="ASP"),
            sum(emar_drugs$drug_group=="CEF")))
cat("Unique subjects with ASP:", length(unique(emar_drugs[drug_group=="ASP", subject_id])), "\n")
cat("Unique subjects with CEF:", length(unique(emar_drugs[drug_group=="CEF", subject_id])), "\n")

###############################################################################
# 2. Link treatment to MSSA episode
###############################################################################
cat("\n--- Step 2: Link treatment to bacteremia episode ---\n")

# Merge emar with candidates to get first_blood_culture_time
emar_with_bc <- merge(
  emar_drugs,
  mssa_cands[, .(subject_id, hadm_id, first_blood_culture_time, admittime, dischtime)],
  by=c("subject_id","hadm_id"),
  all.x=FALSE  # only keep rows matching our candidates
)

cat(sprintf("After linking to candidate hadm_ids: %d rows\n", nrow(emar_with_bc)))

# Compute time difference from first blood culture
emar_with_bc[, hours_from_bc := as.numeric(difftime(charttime, first_blood_culture_time, units="hours"))]

# Treatment window: [-48h, +5 days (=120h)] of first blood culture
emar_in_window <- emar_with_bc[hours_from_bc >= -48 & hours_from_bc <= 120]
cat(sprintf("EMAR rows within [-48h, +120h] window: %d\n", nrow(emar_in_window)))

###############################################################################
# 3. New-user design: exclude patients with prior ASP/CEF exposure
###############################################################################
cat("\n--- Step 3: New-user design - check for prior exposure ---\n")

# Prior exposure = drug given BEFORE the window (more than 48h before blood culture)
prior_exposure <- emar_with_bc[hours_from_bc < -48, .(
  prior_ASP = any(drug_group == "ASP"),
  prior_CEF = any(drug_group == "CEF")
), by=.(subject_id, hadm_id)]

cat(sprintf("Patients with prior ASP exposure: %d\n",
            sum(prior_exposure$prior_ASP, na.rm=TRUE)))
cat(sprintf("Patients with prior CEF exposure: %d\n",
            sum(prior_exposure$prior_CEF, na.rm=TRUE)))

###############################################################################
# 4. For each candidate, get first qualifying treatment
###############################################################################
cat("\n--- Step 4: Identify first treatment ---\n")

# For each subject/hadm, find first drug administration within window
first_treat <- emar_in_window[
  !is.na(charttime),
  .(
    first_time = min(charttime, na.rm=TRUE),
    drug_group_first = drug_group[which.min(charttime)],
    drug_name_first = medication[which.min(charttime)]
  ),
  by=.(subject_id, hadm_id)
]

cat(sprintf("Candidates with any treatment in window: %d\n", nrow(first_treat)))
cat("First treatment group:\n")
print(table(first_treat$drug_group_first))

###############################################################################
# 5. Classify W and t0
###############################################################################
cat("\n--- Step 5: Assign W and t0 ---\n")

# Merge back
cohort <- merge(
  mssa_cands,
  first_treat[, .(subject_id, hadm_id,
                  t0 = first_time,
                  W = fifelse(drug_group_first=="ASP", 1L, 0L),
                  drug_group = drug_group_first,
                  drug_name = drug_name_first)],
  by=c("subject_id","hadm_id"),
  all.x=FALSE
)

cat(sprintf("Cohort with treatment: %d\n", nrow(cohort)))
cat(sprintf("  ASP (W=1): %d\n", sum(cohort$W==1)))
cat(sprintf("  CEF (W=0): %d\n", sum(cohort$W==0)))

# Remove those with prior exposure
cohort <- merge(
  cohort,
  prior_exposure,
  by=c("subject_id","hadm_id"),
  all.x=TRUE
)
cohort[is.na(prior_ASP), prior_ASP := FALSE]
cohort[is.na(prior_CEF), prior_CEF := FALSE]

n_before_prior <- nrow(cohort)
# Remove if prior exposure to SAME drug group
cohort <- cohort[!(W==1 & prior_ASP) & !(W==0 & prior_CEF)]
cat(sprintf("\nAfter removing prior exposed: %d (removed %d)\n",
            nrow(cohort), n_before_prior - nrow(cohort)))

###############################################################################
# 6. Check for co-administration (received both) and handle
###############################################################################
cat("\n--- Step 6: Check for switching/co-administration ---\n")

# For each patient/hadm in window, check if they received BOTH
both_drugs <- emar_in_window[, .(
  had_ASP = any(drug_group=="ASP"),
  had_CEF = any(drug_group=="CEF")
), by=.(subject_id, hadm_id)]

both_drugs[, received_both := had_ASP & had_CEF]
cat(sprintf("Patients receiving both drugs in window: %d\n",
            sum(both_drugs$received_both)))

# Mark switchers - those who switch within 24h of t0
# (Keep in cohort but flag)
cohort <- merge(cohort,
                both_drugs[, .(subject_id, hadm_id, received_both)],
                by=c("subject_id","hadm_id"), all.x=TRUE)
cohort[is.na(received_both), received_both := FALSE]
cat(sprintf("Switchers in final cohort: %d\n", sum(cohort$received_both)))

# For primary analysis, "intention to treat" from first drug:
# patients are already classified by first drug

###############################################################################
# 7. Compute time from blood culture to treatment
###############################################################################
cohort[, hours_bc_to_t0 := as.numeric(difftime(t0, first_blood_culture_time, units="hours"))]
cat("\nHours from blood culture to t0:\n")
print(summary(cohort$hours_bc_to_t0))

###############################################################################
# 8. Final cohort summary
###############################################################################
cat("\n--- Final cohort summary ---\n")
cat(sprintf("Total: %d (ASP: %d, CEF: %d)\n",
            nrow(cohort), sum(cohort$W==1), sum(cohort$W==0)))
cat(sprintf("MSSA confirmed: %d, Unknown: %d\n",
            sum(cohort$mssa_confirmed==TRUE, na.rm=TRUE),
            sum(is.na(cohort$mssa_confirmed))))
cat(sprintf("Switchers (received both): %d\n", sum(cohort$received_both)))

# Select output columns
cohort_out <- cohort[, .(
  subject_id, hadm_id,
  W, t0, drug_name, drug_group,
  received_both,
  mssa_confirmed,
  first_blood_culture_time,
  hours_bc_to_t0,
  admittime, dischtime,
  age, gender, race,
  admission_type, hospital_expire_flag,
  dod, prior_ASP, prior_CEF
)]

setorder(cohort_out, subject_id, hadm_id)

saveRDS(cohort_out, file.path(DERIVED_DIR, "cohort_treated.rds"))
cat(sprintf("\nSaved: derived/cohort_treated.rds (%d rows)\n", nrow(cohort_out)))

cat("\n=== 02_define_treatment_and_time_zero.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
