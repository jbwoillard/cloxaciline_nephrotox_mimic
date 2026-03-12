## ============================================================
## 00_audit_data_and_drug_feasibility.R
## Audit MIMIC-IV tables and drug/MSSA feasibility
## ============================================================

library(data.table)
library(duckdb)
library(DBI)
library(stringr)

# ── paths ─────────────────────────────────────────────────────
DATA  <- "mimic-iv-3.1"
HOSP  <- file.path(DATA, "hosp")
ICU   <- file.path(DATA, "icu")
TABS  <- "tables"
RPT   <- "report"
LOGS  <- "logs"
dir.create(TABS, showWarnings = FALSE)
dir.create(RPT,  showWarnings = FALSE)
dir.create(LOGS, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(LOGS, "00_audit.log"), append = TRUE)
}

log_msg("=== STEP 0: AUDIT STARTED ===")

# ── 1. Available tables ───────────────────────────────────────
hosp_tables <- list.files(HOSP, pattern = "\\.csv(\\.gz)?$")
icu_tables  <- list.files(ICU,  pattern = "\\.csv(\\.gz)?$")
log_msg("HOSP tables: ", paste(hosp_tables, collapse = ", "))
log_msg("ICU  tables: ", paste(icu_tables,  collapse = ", "))

# ── 2. DuckDB connection for efficient queries ────────────────
con <- dbConnect(duckdb())

# Helper: query gzipped CSV
gz_query <- function(con, file_path, sql_template) {
  # DuckDB can read .csv.gz natively
  sql <- gsub("__FILE__", shQuote(file_path), sql_template)
  dbGetQuery(con, sql)
}

# ── 3. Drug feasibility audit ─────────────────────────────────
log_msg("Starting drug feasibility audit ...")

# Search terms (case-insensitive)
drug_terms <- c("cloxacillin", "cloxacillin", "oxacillin", "nafcillin",
                "cefazolin", "cefazoline")
pattern_asp   <- "(?i)(cloxacillin|oxacillin|nafcillin)"
pattern_cef   <- "(?i)(cefazolin)"

## 3a. prescriptions table
log_msg("  Querying prescriptions ...")
rx_file <- file.path(HOSP, "prescriptions.csv.gz")

rx_asp <- dbGetQuery(con, sprintf(
  "SELECT drug, COUNT(*) AS n_rx,
          COUNT(DISTINCT subject_id) AS n_patients,
          COUNT(DISTINCT hadm_id)    AS n_hadm
   FROM read_csv_auto('%s')
   WHERE regexp_matches(drug, '(?i)(cloxacillin|oxacillin|nafcillin)')
   GROUP BY drug ORDER BY n_rx DESC", rx_file))

rx_cef <- dbGetQuery(con, sprintf(
  "SELECT drug, COUNT(*) AS n_rx,
          COUNT(DISTINCT subject_id) AS n_patients,
          COUNT(DISTINCT hadm_id)    AS n_hadm
   FROM read_csv_auto('%s')
   WHERE regexp_matches(drug, '(?i)(cefazolin)')
   GROUP BY drug ORDER BY n_rx DESC", rx_file))

log_msg("    ASP prescriptions rows: ", nrow(rx_asp))
log_msg("    CEF prescriptions rows: ", nrow(rx_cef))
print(rx_asp)
print(rx_cef)

## 3b. pharmacy table
log_msg("  Querying pharmacy ...")
ph_file <- file.path(HOSP, "pharmacy.csv.gz")

ph_asp <- dbGetQuery(con, sprintf(
  "SELECT medication, COUNT(*) AS n_rx,
          COUNT(DISTINCT subject_id) AS n_patients,
          COUNT(DISTINCT hadm_id)    AS n_hadm
   FROM read_csv_auto('%s')
   WHERE regexp_matches(medication, '(?i)(cloxacillin|oxacillin|nafcillin)')
   GROUP BY medication ORDER BY n_rx DESC", ph_file))

ph_cef <- dbGetQuery(con, sprintf(
  "SELECT medication, COUNT(*) AS n_rx,
          COUNT(DISTINCT subject_id) AS n_patients,
          COUNT(DISTINCT hadm_id)    AS n_hadm
   FROM read_csv_auto('%s')
   WHERE regexp_matches(medication, '(?i)(cefazolin)')
   GROUP BY medication ORDER BY n_rx DESC", ph_file))

log_msg("    ASP pharmacy rows: ", nrow(ph_asp))
log_msg("    CEF pharmacy rows: ", nrow(ph_cef))
print(ph_asp)
print(ph_cef)

## 3c. emar table
log_msg("  Querying emar ...")
emar_file <- file.path(HOSP, "emar.csv.gz")

emar_asp <- dbGetQuery(con, sprintf(
  "SELECT medication, COUNT(*) AS n_admin,
          COUNT(DISTINCT subject_id) AS n_patients,
          COUNT(DISTINCT hadm_id)    AS n_hadm
   FROM read_csv_auto('%s')
   WHERE regexp_matches(medication, '(?i)(cloxacillin|oxacillin|nafcillin)')
   GROUP BY medication ORDER BY n_admin DESC", emar_file))

emar_cef <- dbGetQuery(con, sprintf(
  "SELECT medication, COUNT(*) AS n_admin,
          COUNT(DISTINCT subject_id) AS n_patients,
          COUNT(DISTINCT hadm_id)    AS n_hadm
   FROM read_csv_auto('%s')
   WHERE regexp_matches(medication, '(?i)(cefazolin)')
   GROUP BY medication ORDER BY n_admin DESC", emar_file))

log_msg("    ASP emar rows: ", nrow(emar_asp))
log_msg("    CEF emar rows: ", nrow(emar_cef))
print(emar_asp)
print(emar_cef)

## 3d. inputevents (ICU)
log_msg("  Querying inputevents (ICU) + d_items for drug labels ...")
items_file <- file.path(ICU, "d_items.csv.gz")
inp_file   <- file.path(ICU, "inputevents.csv.gz")

items_asp <- dbGetQuery(con, sprintf(
  "SELECT itemid, label, category
   FROM read_csv_auto('%s')
   WHERE regexp_matches(label, '(?i)(cloxacillin|oxacillin|nafcillin|cefazolin)')",
  items_file))

log_msg("    ICU drug items found: ", nrow(items_asp))
print(items_asp)

if (nrow(items_asp) > 0) {
  itemid_list <- paste(items_asp$itemid, collapse = ", ")
  inp_drugs <- dbGetQuery(con, sprintf(
    "SELECT i.itemid, d.label,
            COUNT(*) AS n_admin,
            COUNT(DISTINCT i.subject_id) AS n_patients,
            COUNT(DISTINCT i.stay_id)    AS n_stays
     FROM read_csv_auto('%s') i
     JOIN read_csv_auto('%s') d USING (itemid)
     WHERE i.itemid IN (%s)
     GROUP BY i.itemid, d.label ORDER BY n_admin DESC",
    inp_file, items_file, itemid_list))
  log_msg("    ICU inputevents drug admins: ", nrow(inp_drugs))
  print(inp_drugs)
} else {
  inp_drugs <- data.frame()
  log_msg("    No ICU drug items matched.")
}

# ── 4. Build drug_feasibility.csv ─────────────────────────────
build_feasibility_row <- function(source, drug_group, df_asp, df_cef) {
  make_row <- function(grp, df) {
    if (nrow(df) == 0) {
      data.frame(source=source, drug_group=grp,
                 n_rows=0, n_patients=0, n_hadm=0, drugs_found="none",
                 stringsAsFactors=FALSE)
    } else {
      ncol_rx  <- if ("n_rx"    %in% names(df)) "n_rx"
                  else if ("n_admin" %in% names(df)) "n_admin" else NA
      drug_col <- if ("drug" %in% names(df)) "drug"
                  else if ("medication" %in% names(df)) "medication"
                  else if ("label" %in% names(df)) "label" else NA
      data.frame(
        source     = source,
        drug_group = grp,
        n_rows     = if (!is.na(ncol_rx)) sum(df[[ncol_rx]]) else NA,
        n_patients = sum(df$n_patients),
        n_hadm     = if ("n_hadm" %in% names(df)) sum(df$n_hadm) else NA,
        drugs_found= if (!is.na(drug_col)) paste(df[[drug_col]], collapse="; ") else "",
        stringsAsFactors=FALSE
      )
    }
  }
  rbind(make_row("ASP (oxacillin/nafcillin)", df_asp),
        make_row("Cefazolin",                  df_cef))
}

feas <- rbind(
  build_feasibility_row("prescriptions", "ASP/CEF", rx_asp,   rx_cef),
  build_feasibility_row("pharmacy",      "ASP/CEF", ph_asp,   ph_cef),
  build_feasibility_row("emar",          "ASP/CEF", emar_asp, emar_cef)
)
feas <- feas[seq(1, nrow(feas), 2), ]  # one row per source/group combo - restructure
# Rebuild properly
feas_asp <- rbind(
  data.frame(source="prescriptions", drug_group="ASP",
             n_rows=ifelse(nrow(rx_asp)>0, sum(rx_asp$n_rx),0),
             n_patients=ifelse(nrow(rx_asp)>0,sum(rx_asp$n_patients),0),
             n_hadm=ifelse(nrow(rx_asp)>0,sum(rx_asp$n_hadm),0),
             drugs_found=ifelse(nrow(rx_asp)>0,paste(rx_asp$drug,collapse="; "),"none")),
  data.frame(source="pharmacy",      drug_group="ASP",
             n_rows=ifelse(nrow(ph_asp)>0, sum(ph_asp$n_rx),0),
             n_patients=ifelse(nrow(ph_asp)>0,sum(ph_asp$n_patients),0),
             n_hadm=ifelse(nrow(ph_asp)>0,sum(ph_asp$n_hadm),0),
             drugs_found=ifelse(nrow(ph_asp)>0,paste(ph_asp$medication,collapse="; "),"none")),
  data.frame(source="emar",          drug_group="ASP",
             n_rows=ifelse(nrow(emar_asp)>0, sum(emar_asp$n_admin),0),
             n_patients=ifelse(nrow(emar_asp)>0,sum(emar_asp$n_patients),0),
             n_hadm=ifelse(nrow(emar_asp)>0,sum(emar_asp$n_hadm),0),
             drugs_found=ifelse(nrow(emar_asp)>0,paste(emar_asp$medication,collapse="; "),"none"))
)
feas_cef <- rbind(
  data.frame(source="prescriptions", drug_group="Cefazolin",
             n_rows=ifelse(nrow(rx_cef)>0, sum(rx_cef$n_rx),0),
             n_patients=ifelse(nrow(rx_cef)>0,sum(rx_cef$n_patients),0),
             n_hadm=ifelse(nrow(rx_cef)>0,sum(rx_cef$n_hadm),0),
             drugs_found=ifelse(nrow(rx_cef)>0,paste(rx_cef$drug,collapse="; "),"none")),
  data.frame(source="pharmacy",      drug_group="Cefazolin",
             n_rows=ifelse(nrow(ph_cef)>0, sum(ph_cef$n_rx),0),
             n_patients=ifelse(nrow(ph_cef)>0,sum(ph_cef$n_patients),0),
             n_hadm=ifelse(nrow(ph_cef)>0,sum(ph_cef$n_hadm),0),
             drugs_found=ifelse(nrow(ph_cef)>0,paste(ph_cef$medication,collapse="; "),"none")),
  data.frame(source="emar",          drug_group="Cefazolin",
             n_rows=ifelse(nrow(emar_cef)>0, sum(emar_cef$n_admin),0),
             n_patients=ifelse(nrow(emar_cef)>0,sum(emar_cef$n_patients),0),
             n_hadm=ifelse(nrow(emar_cef)>0,sum(emar_cef$n_hadm),0),
             drugs_found=ifelse(nrow(emar_cef)>0,paste(emar_cef$medication,collapse="; "),"none"))
)
drug_feas <- rbind(feas_asp, feas_cef)
write.csv(drug_feas, file.path(TABS, "drug_feasibility.csv"), row.names = FALSE)
log_msg("  Saved: tables/drug_feasibility.csv")
print(drug_feas)

# ── 5. MSSA feasibility ───────────────────────────────────────
log_msg("Starting MSSA feasibility audit ...")
micro_file <- file.path(HOSP, "microbiologyevents.csv.gz")

# Get column names first
micro_cols <- dbGetQuery(con, sprintf(
  "DESCRIBE SELECT * FROM read_csv_auto('%s') LIMIT 1", micro_file))
log_msg("  microbiologyevents columns: ", paste(micro_cols$column_name, collapse=", "))

# SA bacteremia: blood cultures with Staphylococcus aureus
sa_blood <- dbGetQuery(con, sprintf(
  "SELECT spec_type_desc, org_name,
          COUNT(*) AS n_events,
          COUNT(DISTINCT subject_id) AS n_patients,
          COUNT(DISTINCT hadm_id)    AS n_hadm
   FROM read_csv_auto('%s')
   WHERE regexp_matches(org_name, '(?i)(staphylococcus aureus)')
     AND regexp_matches(spec_type_desc, '(?i)(blood)')
   GROUP BY spec_type_desc, org_name
   ORDER BY n_events DESC", micro_file))

log_msg("  SA blood cultures: ", nrow(sa_blood), " rows")
print(sa_blood)

# Check for susceptibility data (methicillin/oxacillin)
sa_suscept <- dbGetQuery(con, sprintf(
  "SELECT ab_name, interpretation,
          COUNT(*) AS n,
          COUNT(DISTINCT subject_id) AS n_patients
   FROM read_csv_auto('%s')
   WHERE regexp_matches(org_name, '(?i)(staphylococcus aureus)')
     AND regexp_matches(spec_type_desc, '(?i)(blood)')
     AND ab_name IS NOT NULL
     AND regexp_matches(ab_name, '(?i)(oxacillin|methicillin|cloxacillin|nafcillin)')
   GROUP BY ab_name, interpretation
   ORDER BY ab_name, n DESC", micro_file))

log_msg("  SA susceptibility (anti-staph): ", nrow(sa_suscept), " rows")
print(sa_suscept)

# All SA blood: total unique patients
sa_total <- dbGetQuery(con, sprintf(
  "SELECT COUNT(DISTINCT subject_id) AS n_patients_sa_blood,
          COUNT(DISTINCT hadm_id)    AS n_hadm_sa_blood
   FROM read_csv_auto('%s')
   WHERE regexp_matches(org_name, '(?i)(staphylococcus aureus)')
     AND regexp_matches(spec_type_desc, '(?i)(blood)')", micro_file))

log_msg("  Total SA blood culture patients: ", sa_total$n_patients_sa_blood)
log_msg("  Total SA blood culture hadm: ",     sa_total$n_hadm_sa_blood)

# MSSA specifically (susceptible to oxacillin = MSSA)
mssa_count <- dbGetQuery(con, sprintf(
  "SELECT COUNT(DISTINCT subject_id) AS n_patients_mssa,
          COUNT(DISTINCT hadm_id)    AS n_hadm_mssa
   FROM read_csv_auto('%s')
   WHERE regexp_matches(org_name, '(?i)(staphylococcus aureus)')
     AND regexp_matches(spec_type_desc, '(?i)(blood)')
     AND ab_name = 'OXACILLIN'
     AND interpretation = 'S'", micro_file))

log_msg("  Confirmed MSSA patients (oxacillin S): ", mssa_count$n_patients_mssa)

# Save microbiology feasibility
micro_feas <- rbind(
  data.frame(criterion="SA_blood_all",
             n_patients=sa_total$n_patients_sa_blood,
             n_hadm=sa_total$n_hadm_sa_blood,
             note="Blood cultures positive for Staphylococcus aureus"),
  data.frame(criterion="MSSA_confirmed_oxacillin_S",
             n_patients=mssa_count$n_patients_mssa,
             n_hadm=mssa_count$n_hadm_mssa,
             note="Oxacillin susceptible (S) = MSSA")
)
write.csv(micro_feas, file.path(TABS, "microbiology_feasibility.csv"), row.names=FALSE)
log_msg("  Saved: tables/microbiology_feasibility.csv")

dbDisconnect(con, shutdown=TRUE)

# ── 6. Feasibility summary report ────────────────────────────
log_msg("Writing feasibility summary ...")

summary_txt <- c(
  "# Feasibility Summary",
  paste("Date:", Sys.time()),
  "",
  "## 1. Data Availability",
  paste("MIMIC-IV version: 3.1"),
  paste("HOSP tables available:", paste(hosp_tables, collapse=", ")),
  paste("ICU tables available:",  paste(icu_tables,  collapse=", ")),
  "",
  "## 2. Drug Feasibility",
  "",
  "### ASP (Anti-Staphylococcal Penicillins: oxacillin/nafcillin)",
  paste("- prescriptions: n_patients =", feas_asp$n_patients[1], "| n_rows =", feas_asp$n_rows[1]),
  paste("- pharmacy:      n_patients =", feas_asp$n_patients[2], "| n_rows =", feas_asp$n_rows[2]),
  paste("- emar:          n_patients =", feas_asp$n_patients[3], "| n_rows =", feas_asp$n_rows[3]),
  paste("- Drugs found:", paste(unique(unlist(strsplit(feas_asp$drugs_found, "; "))), collapse="; ")),
  "",
  "### Cefazolin",
  paste("- prescriptions: n_patients =", feas_cef$n_patients[1], "| n_rows =", feas_cef$n_rows[1]),
  paste("- pharmacy:      n_patients =", feas_cef$n_patients[2], "| n_rows =", feas_cef$n_rows[2]),
  paste("- emar:          n_patients =", feas_cef$n_patients[3], "| n_rows =", feas_cef$n_rows[3]),
  paste("- Drugs found:", paste(unique(unlist(strsplit(feas_cef$drugs_found, "; "))), collapse="; ")),
  "",
  "## 3. MSSA Feasibility",
  paste("- SA blood culture patients (total):", sa_total$n_patients_sa_blood),
  paste("- Confirmed MSSA (oxacillin S):", mssa_count$n_patients_mssa),
  "",
  "## 4. Decision",
  ifelse(nrow(rx_asp) > 0 & nrow(rx_cef) > 0,
    paste0("**FEASIBLE**: Both ASP (oxacillin/nafcillin) and cefazolin are present. ",
           "Cloxacillin is NOT present in MIMIC-IV (US database). ",
           "Will use **ASP (oxacillin + nafcillin) vs cefazolin** design."),
    "**REVIEW REQUIRED**: drug availability insufficient for planned analysis."),
  "",
  "Note: Cloxacillin is not used in the US and is absent from MIMIC-IV. ",
  "The study will use oxacillin and/or nafcillin as anti-staphylococcal penicillin (ASP) comparator,",
  "consistent with the conceptual replication approach described in the prompt.",
  "",
  "## 5. Tables produced",
  "- tables/drug_feasibility.csv",
  "- tables/microbiology_feasibility.csv"
)

writeLines(summary_txt, file.path(RPT, "feasibility_summary.md"))
log_msg("  Saved: report/feasibility_summary.md")
log_msg("=== STEP 0: AUDIT COMPLETE ===")
