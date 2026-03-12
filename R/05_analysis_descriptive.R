###############################################################################
# 05_analysis_descriptive.R
# Build analysis dataset and generate Table 1
# Output: derived/analysis_dataset.rds, tables/table1.csv, tables/raw_outcomes.csv
###############################################################################

library(dplyr)
library(data.table)
library(tableone)
library(ggplot2)
library(lubridate)

cat("=== 05_analysis_descriptive.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

WORK_DIR <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
TABLES_DIR  <- file.path(WORK_DIR, "tables")
FIGURES_DIR <- file.path(WORK_DIR, "figures")

# Load all components
cohort   <- readRDS(file.path(DERIVED_DIR, "cohort_treated.rds"))
outcomes <- readRDS(file.path(DERIVED_DIR, "outcomes.rds"))
covars   <- readRDS(file.path(DERIVED_DIR, "covariates.rds"))

cat(sprintf("Cohort: %d\nOutcomes: %d\nCovariates: %d\n",
            nrow(cohort), nrow(outcomes), nrow(covars)))

###############################################################################
# 1. Build analysis dataset
###############################################################################
cat("\n--- Step 1: Build analysis dataset ---\n")

# Join all
anal <- merge(
  cohort[, .(subject_id, hadm_id, W, t0, drug_name, drug_group,
             received_both, mssa_confirmed, hours_bc_to_t0,
             admittime, dischtime, hospital_expire_flag, dod)],
  outcomes[, .(subject_id, hadm_id,
               baseline_creat, peak_creat, creat_day2,
               delta_creat_7d, creat_ratio,
               AKI_7d, AKI_stage, aki_at_baseline, preadm_creat,
               died_inhosp, died_7d, composite_aki_death,
               n_baseline_measures, n_peak_measures)],
  by=c("subject_id","hadm_id"),
  all.x=TRUE
)

anal <- merge(
  anal,
  covars[, !c("t0"), with=FALSE],
  by=c("subject_id","hadm_id"),
  all.x=TRUE
)

cat(sprintf("Analysis dataset: %d rows, %d cols\n", nrow(anal), ncol(anal)))
cat(sprintf("ASP (W=1): %d, CEF (W=0): %d\n", sum(anal$W==1), sum(anal$W==0)))

# Treatment label
anal[, treatment := fifelse(W==1, "ASP (nafcillin/oxacillin)", "Cefazolin")]

# Hospital LOS
anal[, hosp_los_days := as.numeric(difftime(dischtime, admittime, units="days"))]

# Save
saveRDS(anal, file.path(DERIVED_DIR, "analysis_dataset.rds"))
cat(sprintf("Saved: derived/analysis_dataset.rds\n"))

###############################################################################
# 2. Table 1: baseline characteristics by treatment group
###############################################################################
cat("\n--- Step 2: Table 1 ---\n")

# Recode categorical variables
anal[, drug_group_label := fifelse(W==1, "ASP", "CEF")]
anal[, icu_at_t0_label := fifelse(icu_at_t0==1, "Yes", "No")]
anal[, mssa_confirmed_label := fifelse(is.na(mssa_confirmed), "Unknown",
                                        fifelse(mssa_confirmed==TRUE, "Confirmed", "Probable"))]

# Variables for Table 1
tab1_vars <- c(
  # Demographics
  "age", "female", "race_cat",
  # Clinical
  "admission_type", "icu_at_t0", "mssa_confirmed_label",
  # Comorbidities
  "ckd_any", "egfr_baseline", "creat_baseline",
  "diabetes", "hypertension", "heart_failure",
  "liver_disease", "cancer", "copd",
  "charlson",
  # Severity
  "vasopressor", "mech_vent",
  # Labs
  "wbc", "platelets", "lactate", "bilirubin",
  "sodium", "potassium", "bun", "bicarbonate",
  # Co-exposures
  "concomitant_vancomycin", "concomitant_aminoglycoside",
  "concomitant_loop_diuretic",
  # Time
  "hours_bc_to_t0"
)

# Keep vars that exist
tab1_vars_exist <- intersect(tab1_vars, names(anal))

# Categorical variables
cat_vars <- c("female", "race_cat", "admission_type", "icu_at_t0",
              "mssa_confirmed_label",
              "ckd_any", "diabetes", "hypertension", "heart_failure",
              "liver_disease", "cancer", "copd",
              "vasopressor", "mech_vent",
              "concomitant_vancomycin", "concomitant_aminoglycoside",
              "concomitant_loop_diuretic")
cat_vars_exist <- intersect(cat_vars, tab1_vars_exist)

# Create Table 1 using tableone
tab1_obj <- CreateTableOne(
  vars = tab1_vars_exist,
  strata = "drug_group_label",
  data = as.data.frame(anal),
  factorVars = cat_vars_exist,
  test = TRUE,
  smd = TRUE
)

tab1_print <- print(tab1_obj, smd=TRUE, quote=FALSE, noSpaces=TRUE,
                    printToggle=FALSE, showAllLevels=FALSE)
cat("Table 1 preview:\n")
print(tab1_print)

# Save
tab1_df <- as.data.frame(tab1_print)
tab1_df$variable <- rownames(tab1_df)
fwrite(tab1_df, file.path(TABLES_DIR, "table1.csv"))
cat("\nSaved: tables/table1.csv\n")

###############################################################################
# 3. SMD analysis
###############################################################################
cat("\n--- Step 3: SMD analysis ---\n")

# Extract SMD from tableone
smd_vals <- ExtractSmd(tab1_obj)
cat("SMD for key variables:\n")
print(round(smd_vals, 3))

# Save SMD
smd_mat <- as.matrix(smd_vals)
smd_df <- data.frame(
  variable = rownames(smd_mat),
  SMD = as.numeric(smd_mat[,1])
)
fwrite(smd_df, file.path(TABLES_DIR, "smd_before_weighting.csv"))
cat("Saved: tables/smd_before_weighting.csv\n")

# Count imbalanced variables (SMD > 0.1)
n_imbalanced <- sum(abs(smd_vals) > 0.1, na.rm=TRUE)
cat(sprintf("\nVariables with SMD > 0.1: %d/%d\n", n_imbalanced, length(smd_vals)))

###############################################################################
# 4. Raw outcome rates
###############################################################################
cat("\n--- Step 4: Raw outcome rates ---\n")

raw_outcomes <- anal[, .(
  N = .N,
  # AKI
  AKI_7d_n = sum(AKI_7d==1, na.rm=TRUE),
  AKI_7d_pct = round(100*mean(AKI_7d==1, na.rm=TRUE), 1),
  AKI_stage1_n = sum(AKI_stage==1, na.rm=TRUE),
  AKI_stage2_n = sum(AKI_stage==2, na.rm=TRUE),
  AKI_stage3_n = sum(AKI_stage==3, na.rm=TRUE),
  # Mortality
  died_7d_n = sum(died_7d==1, na.rm=TRUE),
  died_7d_pct = round(100*mean(died_7d==1, na.rm=TRUE), 1),
  died_inhosp_n = sum(died_inhosp==1, na.rm=TRUE),
  died_inhosp_pct = round(100*mean(died_inhosp==1, na.rm=TRUE), 1),
  # Composite
  composite_n = sum(composite_aki_death==1, na.rm=TRUE),
  composite_pct = round(100*mean(composite_aki_death==1, na.rm=TRUE), 1),
  # Creatinine
  mean_baseline_creat = round(mean(baseline_creat, na.rm=TRUE), 2),
  mean_peak_creat = round(mean(peak_creat, na.rm=TRUE), 2),
  mean_delta_creat = round(mean(delta_creat_7d, na.rm=TRUE), 2)
), by=.(W, treatment)]

print(raw_outcomes)
fwrite(raw_outcomes, file.path(TABLES_DIR, "raw_outcomes.csv"))
cat("Saved: tables/raw_outcomes.csv\n")

# Chi-square test for AKI
aki_tab <- table(anal$W, anal$AKI_7d)
cat("\nAKI chi-square test:\n")
print(chisq.test(aki_tab))

# Crude OR
if(all(dim(aki_tab) == c(2,2))) {
  or_crude <- (aki_tab[2,2]/aki_tab[2,1]) / (aki_tab[1,2]/aki_tab[1,1])
  cat(sprintf("Crude OR (ASP vs CEF): %.2f\n", or_crude))
}

###############################################################################
# 5. Creatinine trajectory summary
###############################################################################
cat("\n--- Step 5: Creatinine trajectory ---\n")

creat_traj <- anal[, .(
  mean_creat = mean(c(baseline_creat, creat_day2, peak_creat), na.rm=TRUE),
  baseline = mean(baseline_creat, na.rm=TRUE),
  day2 = mean(creat_day2, na.rm=TRUE),
  peak = mean(peak_creat, na.rm=TRUE)
), by=treatment]

cat("Creatinine trajectory (mean mg/dL):\n")
print(creat_traj)

###############################################################################
# 6. Quick summary figure: AKI rates
###############################################################################
cat("\n--- Step 6: Save summary figures ---\n")

# AKI bar plot
aki_summary <- anal[!is.na(AKI_7d), .(
  pct_AKI = 100*mean(AKI_7d==1)
), by=treatment]

p_aki <- ggplot(aki_summary, aes(x=treatment, y=pct_AKI, fill=treatment)) +
  geom_col(width=0.5, alpha=0.8) +
  geom_text(aes(label=sprintf("%.1f%%", pct_AKI)), vjust=-0.5, size=4) +
  scale_fill_manual(values=c("ASP (nafcillin/oxacillin)"="#e74c3c",
                              "Cefazolin"="#3498db")) +
  labs(title="AKI incidence at 7 days by treatment group",
       subtitle="Unadjusted rates",
       x="Treatment", y="AKI incidence (%)") +
  theme_bw(base_size=12) +
  theme(legend.position="none") +
  ylim(0, 60)

ggsave(file.path(FIGURES_DIR, "raw_aki_rates.png"), p_aki,
       width=6, height=5, dpi=150)
cat("Saved: figures/raw_aki_rates.png\n")

###############################################################################
# 7. Overall summary stats
###############################################################################
cat("\n--- Overall analysis dataset summary ---\n")
cat(sprintf("Total patients: %d\n", nrow(anal)))
cat(sprintf("ASP (W=1): %d (%.1f%%)\n", sum(anal$W==1), 100*mean(anal$W==1)))
cat(sprintf("CEF (W=0): %d (%.1f%%)\n", sum(anal$W==0), 100*mean(anal$W==0)))
cat(sprintf("Overall AKI_7d: %d/%d (%.1f%%)\n",
            sum(anal$AKI_7d==1, na.rm=TRUE),
            sum(!is.na(anal$AKI_7d)),
            100*mean(anal$AKI_7d==1, na.rm=TRUE)))
cat(sprintf("Overall in-hospital mortality: %d/%d (%.1f%%)\n",
            sum(anal$died_inhosp==1, na.rm=TRUE),
            nrow(anal),
            100*mean(anal$died_inhosp==1, na.rm=TRUE)))

cat("\n=== 05_analysis_descriptive.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
