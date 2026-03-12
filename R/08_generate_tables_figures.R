###############################################################################
# 08_generate_tables_figures.R
# Generate all publication-quality tables and figures
# Output: figures/*.png, tables/*.csv
###############################################################################

library(ggplot2)
library(data.table)
library(dplyr)
library(purrr)
library(grf)

cat("=== 08_generate_tables_figures.R ===\n")
cat("Start:", format(Sys.time()), "\n\n")

WORK_DIR    <- "/Users/woillp01/Documents/cyrielle_mimic_cloxa"
DERIVED_DIR <- file.path(WORK_DIR, "derived")
RESULTS_DIR <- file.path(WORK_DIR, "results")
FIGURES_DIR <- file.path(WORK_DIR, "figures")
TABLES_DIR  <- file.path(WORK_DIR, "tables")

# Load data
anal <- readRDS(file.path(DERIVED_DIR, "analysis_dataset.rds"))
anal_cate <- readRDS(file.path(DERIVED_DIR, "analysis_dataset_with_cate.rds"))
model_objs <- readRDS(file.path(RESULTS_DIR, "model_objects.rds"))
sens_results <- fread(file.path(RESULTS_DIR, "sensitivity_results.csv"))
smd_data <- fread(file.path(TABLES_DIR, "smd_before_weighting.csv"))
vi_data <- fread(file.path(RESULTS_DIR, "variable_importance.csv"))
ate_results <- fread(file.path(RESULTS_DIR, "ate_results.csv"))

# Color palette
COLORS <- c("ASP"="#e74c3c", "CEF"="#3498db")
THEME_BASE <- theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"))

cat("Data loaded\n")

###############################################################################
# Figure 1: Cohort Flowchart (text-based with ggplot2)
###############################################################################
cat("\n--- Figure 1: Cohort flowchart ---\n")

# Build flowchart data
flow_data <- data.frame(
  step = c(
    "SA blood culture patients\n(MIMIC-IV 3.1)",
    "MSSA/Unknown susceptibility\n(exclude MRSA)",
    "Age >= 18, with hadm_id",
    "With ASP or CEF\n in treatment window",
    "After removing prior exposure\n(new-user design)",
    "Final analysis cohort"
  ),
  N = c(2223, 1001, 1001, 458, 446, 446),
  removed = c(NA, 473, 0, 543, 12, 0),
  remove_reason = c(
    NA,
    "MRSA confirmed (n=473)",
    "—",
    "No treatment in window (n=543)",
    "Prior ASP/CEF exposure (n=12)",
    "—"
  )
)

# Simple text flowchart
p_flow <- ggplot() +
  # Boxes
  geom_rect(data=flow_data,
            aes(xmin=0.1, xmax=0.9, ymin=seq(nrow(flow_data),1,-1)-0.35,
                ymax=seq(nrow(flow_data),1,-1)+0.35),
            fill="white", color="navy", linewidth=0.7) +
  # Text in boxes
  geom_text(data=flow_data,
            aes(x=0.5, y=seq(nrow(flow_data),1,-1),
                label=paste0(step, "\nN = ", format(N, big.mark=","))),
            size=3.2, fontface="bold") +
  # Arrows between boxes
  geom_segment(data=flow_data[-1,],
               aes(x=0.5, xend=0.5,
                   y=seq(nrow(flow_data)-1,1,-1)+0.35,
                   yend=seq(nrow(flow_data)-1,1,-1)+0.65-0.35+0.35),
               arrow=arrow(length=unit(0.2,"cm"), type="closed"),
               color="navy") +
  xlim(0, 1.5) +
  ylim(0, nrow(flow_data)+0.5) +
  labs(title="Cohort Selection Flow Chart",
       subtitle="MSSA bacteremia patients treated with ASP or cefazolin") +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

ggsave(file.path(FIGURES_DIR, "flowchart.png"), p_flow,
       width=7, height=8, dpi=150)
cat("Saved: figures/flowchart.png\n")

###############################################################################
# Figure 2: Love Plot (SMD before weighting)
###############################################################################
cat("\n--- Figure 2: Love plot ---\n")

# Clean SMD names for display
smd_plot <- smd_data[!is.na(SMD), .(variable, SMD=abs(SMD))]
smd_plot[, var_label := fcase(
  variable=="age", "Age",
  variable=="female", "Female sex",
  variable=="race_cat", "Race (non-White)",
  variable=="admission_type", "Admission type",
  variable=="icu_at_t0", "ICU at t0",
  variable=="creat_baseline", "Baseline creatinine",
  variable=="egfr_baseline", "Baseline eGFR",
  variable=="ckd_any", "Chronic kidney disease",
  variable=="diabetes", "Diabetes",
  variable=="hypertension", "Hypertension",
  variable=="heart_failure", "Heart failure",
  variable=="liver_disease", "Liver disease",
  variable=="cancer", "Cancer",
  variable=="charlson", "Charlson index",
  variable=="vasopressor", "Vasopressor use",
  variable=="mech_vent", "Mechanical ventilation",
  variable=="wbc", "WBC",
  variable=="platelets", "Platelets",
  variable=="lactate", "Lactate",
  variable=="bilirubin", "Bilirubin",
  variable=="sodium", "Sodium",
  variable=="bicarbonate", "Bicarbonate",
  variable=="bun", "BUN",
  variable=="concomitant_vancomycin", "Concomitant vancomycin",
  variable=="concomitant_aminoglycoside", "Concomitant aminoglycoside",
  variable=="concomitant_loop_diuretic", "Concomitant loop diuretic",
  variable=="hours_bc_to_t0", "Hours: culture to treatment",
  variable=="mssa_confirmed_label", "MSSA confirmed",
  default = variable
)]

setorder(smd_plot, SMD)
smd_plot[, var_label := factor(var_label, levels=smd_plot$var_label)]

p_love <- ggplot(smd_plot, aes(x=SMD, y=var_label)) +
  geom_point(size=2.5, color="steelblue") +
  geom_segment(aes(x=0, xend=SMD, y=var_label, yend=var_label),
               color="steelblue", linewidth=0.5) +
  geom_vline(xintercept=0.1, linetype="dashed", color="orange", linewidth=0.8) +
  geom_vline(xintercept=0.2, linetype="dashed", color="red", linewidth=0.8) +
  geom_vline(xintercept=0, color="black", linewidth=0.4) +
  scale_x_continuous(breaks=seq(0, 0.5, 0.1), limits=c(0, 0.55)) +
  labs(title="Standardized Mean Differences (SMD) Before Weighting",
       subtitle="ASP vs. Cefazolin — Unadjusted",
       x="Absolute Standardized Mean Difference",
       y="") +
  annotate("text", x=0.105, y=2, label="0.1", angle=90, size=3, color="orange") +
  annotate("text", x=0.205, y=2, label="0.2", angle=90, size=3, color="red") +
  THEME_BASE +
  theme(axis.text.y=element_text(size=9))

ggsave(file.path(FIGURES_DIR, "love_plot.png"), p_love,
       width=8, height=7, dpi=150)
cat("Saved: figures/love_plot.png\n")

###############################################################################
# Figure 3: Propensity score distribution
###############################################################################
cat("\n--- Figure 3: Propensity score distribution ---\n")

# W.hat is the propensity score from GRF
ps_df <- data.table(
  W = model_objs$W_use,
  PS = model_objs$cf_full$W.hat
)
ps_df[, Treatment := fifelse(W==1, "ASP (nafcillin/oxacillin)", "Cefazolin")]

p_ps <- ggplot(ps_df, aes(x=PS, fill=Treatment, color=Treatment)) +
  geom_histogram(aes(y=after_stat(density)), bins=25, alpha=0.6,
                 position="identity") +
  geom_density(alpha=0, linewidth=1) +
  scale_fill_manual(values=c("ASP (nafcillin/oxacillin)"="#e74c3c",
                              "Cefazolin"="#3498db")) +
  scale_color_manual(values=c("ASP (nafcillin/oxacillin)"="#e74c3c",
                               "Cefazolin"="#3498db")) +
  labs(title="Propensity Score Distribution",
       subtitle="P(W=ASP | X) estimated by causal forest",
       x="Estimated Propensity Score",
       y="Density",
       fill="Treatment", color="Treatment") +
  THEME_BASE +
  theme(legend.position="bottom")

ggsave(file.path(FIGURES_DIR, "propensity_score.png"), p_ps,
       width=7, height=5, dpi=150)
cat("Saved: figures/propensity_score.png\n")

###############################################################################
# Figure 4: CATE distribution histogram
###############################################################################
cat("\n--- Figure 4: CATE distribution ---\n")

cate_df <- data.table(
  W = model_objs$W_use,
  CATE = model_objs$cate_all_predictions
)
cate_df[, Treatment := fifelse(W==1, "ASP (nafcillin/oxacillin)", "Cefazolin")]

ate_val <- model_objs$ate_full["estimate"]

p_cate <- ggplot(cate_df, aes(x=CATE)) +
  geom_histogram(bins=30, fill="steelblue", alpha=0.7, color="white") +
  geom_vline(xintercept=ate_val, color="red", linetype="dashed", linewidth=1.2) +
  geom_vline(xintercept=0, color="black", linewidth=0.8) +
  annotate("text", x=ate_val+0.01, y=Inf, vjust=1.5,
           label=sprintf("ATE = %.3f", ate_val), color="red", size=3.5) +
  labs(title="Distribution of Conditional Average Treatment Effects (CATE)",
       subtitle="Estimated effect of ASP vs. cefazolin on AKI at 7 days",
       x="CATE (Risk Difference: ASP − Cefazolin)",
       y="Count") +
  THEME_BASE

ggsave(file.path(FIGURES_DIR, "cate_distribution.png"), p_cate,
       width=7, height=5, dpi=150)
cat("Saved: figures/cate_distribution.png\n")

###############################################################################
# Figure 5: Creatinine trajectory
###############################################################################
cat("\n--- Figure 5: Creatinine trajectory ---\n")

# Build trajectory data
creat_traj_long <- rbind(
  anal[!is.na(baseline_creat), .(
    subject_id, hadm_id, W, treatment,
    time_label = "Baseline (t0)",
    time_num = 0,
    creatinine = baseline_creat
  )],
  anal[!is.na(creat_day2), .(
    subject_id, hadm_id, W, treatment,
    time_label = "Day 2-3",
    time_num = 2,
    creatinine = creat_day2
  )],
  anal[!is.na(peak_creat), .(
    subject_id, hadm_id, W, treatment,
    time_label = "Peak (0-7d)",
    time_num = 7,
    creatinine = peak_creat
  )]
)

# Mean and SE by time and treatment
traj_summary <- creat_traj_long[, .(
  mean_cr = mean(creatinine, na.rm=TRUE),
  se_cr = sd(creatinine, na.rm=TRUE)/sqrt(.N),
  n = .N
), by=.(treatment, time_num, time_label)]

traj_summary[, lower := mean_cr - 1.96*se_cr]
traj_summary[, upper := mean_cr + 1.96*se_cr]

p_traj <- ggplot(traj_summary,
                  aes(x=time_num, y=mean_cr, color=treatment, fill=treatment)) +
  geom_line(linewidth=1.2) +
  geom_point(size=3) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.15, color=NA) +
  scale_x_continuous(breaks=c(0,2,7), labels=c("Baseline\n(t0)", "Day 2-3", "Peak\n(0-7d)")) +
  scale_color_manual(values=c("ASP (nafcillin/oxacillin)"="#e74c3c",
                               "Cefazolin"="#3498db"),
                     name="Treatment") +
  scale_fill_manual(values=c("ASP (nafcillin/oxacillin)"="#e74c3c",
                              "Cefazolin"="#3498db"),
                    name="Treatment") +
  labs(title="Creatinine Trajectory by Treatment Group",
       subtitle="Mean ± 95% CI — unadjusted",
       x="Time Point",
       y="Serum Creatinine (mg/dL)") +
  THEME_BASE +
  theme(legend.position="bottom")

ggsave(file.path(FIGURES_DIR, "creatinine_trajectory.png"), p_traj,
       width=7, height=5, dpi=150)
cat("Saved: figures/creatinine_trajectory.png\n")

###############################################################################
# Figure 6: Variable importance plot
###############################################################################
cat("\n--- Figure 6: Variable importance ---\n")

vi_plot <- head(vi_data[order(-importance)], 15)
vi_plot[, var_label := fcase(
  variable=="age", "Age",
  variable=="wbc", "WBC",
  variable=="glucose", "Glucose",
  variable=="potassium", "Potassium",
  variable=="lactate", "Lactate",
  variable=="bicarbonate", "Bicarbonate",
  variable=="bun", "BUN",
  variable=="sodium", "Sodium",
  variable=="platelets", "Platelets",
  variable=="bilirubin", "Bilirubin",
  variable=="egfr_baseline", "Baseline eGFR",
  variable=="creat_baseline", "Baseline creatinine",
  variable=="charlson", "Charlson index",
  variable=="race_white", "Race: White",
  variable=="concomitant_loop_diuretic", "Concomitant loop diuretic",
  default = variable
)]

setorder(vi_plot, importance)
vi_plot[, var_label := factor(var_label, levels=vi_plot$var_label)]

p_vi <- ggplot(vi_plot, aes(x=importance, y=var_label)) +
  geom_col(fill="steelblue", alpha=0.8) +
  geom_text(aes(label=sprintf("%.3f", importance)), hjust=-0.1, size=3) +
  labs(title="Causal Forest Variable Importance",
       subtitle="Top 15 variables for CATE heterogeneity",
       x="Variable Importance",
       y="") +
  xlim(0, max(vi_plot$importance)*1.2) +
  THEME_BASE

ggsave(file.path(FIGURES_DIR, "variable_importance.png"), p_vi,
       width=7, height=6, dpi=150)
cat("Saved: figures/variable_importance.png\n")

###############################################################################
# Figure 7: Sensitivity analysis forest plot
###############################################################################
cat("\n--- Figure 7: Sensitivity analysis forest plot ---\n")

# Prepare sensitivity results for forest plot
sens_plot <- sens_results[!is.na(ATE_estimate)]
sens_plot[, analysis_label := analysis]
# Order: main on top
sens_plot[, order_num := .I]
setorder(sens_plot, -order_num)
sens_plot[, analysis_label := factor(analysis_label, levels=analysis_label)]

p_forest <- ggplot(sens_plot,
                    aes(x=ATE_estimate, y=analysis_label)) +
  geom_point(size=3, color="navy") +
  geom_errorbarh(aes(xmin=ATE_lower, xmax=ATE_upper),
                 height=0.3, color="navy", linewidth=0.8) +
  geom_vline(xintercept=0, color="black", linewidth=0.8) +
  geom_vline(xintercept=sens_plot[analysis=="Primary analysis", ATE_estimate],
             color="red", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label=sprintf("%.3f", ATE_estimate)),
            nudge_x=0.015, size=3) +
  scale_x_continuous(breaks=seq(-0.1, 0.5, 0.1)) +
  labs(title="Sensitivity Analysis: ATE Estimates",
       subtitle="Effect of ASP vs. cefazolin on AKI at 7 days\n(Overlap-weighted causal forest estimates, 95% CI)",
       x="Average Treatment Effect (Risk Difference)",
       y="") +
  THEME_BASE +
  theme(axis.text.y=element_text(size=8))

ggsave(file.path(FIGURES_DIR, "sensitivity_forest_plot.png"), p_forest,
       width=10, height=7, dpi=150)
cat("Saved: figures/sensitivity_forest_plot.png\n")

###############################################################################
# Figure 8: AKI rates bar chart (detailed)
###############################################################################
cat("\n--- Figure 8: AKI stage distribution ---\n")

aki_stage_dist <- anal[!is.na(AKI_stage), .(
  N = .N
), by=.(treatment, AKI_stage)]
aki_stage_dist[, total := sum(N), by=treatment]
aki_stage_dist[, pct := 100*N/total]
aki_stage_dist[, stage_label := fcase(
  AKI_stage==0, "No AKI",
  AKI_stage==1, "Stage 1",
  AKI_stage==2, "Stage 2",
  AKI_stage==3, "Stage 3"
)]
aki_stage_dist[, stage_label := factor(stage_label,
                                         levels=c("No AKI","Stage 1","Stage 2","Stage 3"))]

p_aki_stage <- ggplot(aki_stage_dist, aes(x=treatment, y=pct, fill=stage_label)) +
  geom_col(width=0.6, alpha=0.85) +
  scale_fill_manual(values=c("No AKI"="#2ecc71", "Stage 1"="#f39c12",
                               "Stage 2"="#e67e22", "Stage 3"="#c0392b"),
                    name="AKI Stage") +
  geom_text(aes(label=sprintf("%.1f%%\n(n=%d)", pct, N)),
            position=position_stack(vjust=0.5), size=3) +
  labs(title="AKI Stage Distribution by Treatment Group",
       subtitle="KDIGO creatinine criteria at 7 days",
       x="Treatment Group",
       y="Percentage of patients (%)") +
  THEME_BASE +
  theme(legend.position="right")

ggsave(file.path(FIGURES_DIR, "aki_stage_distribution.png"), p_aki_stage,
       width=7, height=5, dpi=150)
cat("Saved: figures/aki_stage_distribution.png\n")

###############################################################################
# Figure 9: CATE by subgroup
###############################################################################
cat("\n--- Figure 9: CATE by key subgroups ---\n")

# Subgroup analysis: CATE by CKD status, ICU, age quartiles
anal_cate_dt <- as.data.table(anal_cate)

# Age quartiles
anal_cate_dt[, age_quartile := cut(age, breaks=quantile(age, probs=0:4/4, na.rm=TRUE),
                                   include.lowest=TRUE,
                                   labels=c("Q1: youngest", "Q2", "Q3", "Q4: oldest"))]

# Build subgroup CATE summaries
sg_results <- list()

# By CKD
d_ckd <- anal_cate_dt[!is.na(ckd_any)]
d_ckd[, level := fifelse(ckd_any==1, "CKD (yes)", "No CKD")]
d_ckd[, subgroup := "CKD"]
sg_ckd <- d_ckd[, .(mean_CATE=mean(CATE,na.rm=TRUE), sd_CATE=sd(CATE,na.rm=TRUE), n=.N), by=.(subgroup,level)]
sg_results$ckd <- sg_ckd

# By ICU
d_icu <- anal_cate_dt[!is.na(icu_at_t0)]
d_icu[, level := fifelse(icu_at_t0==1, "ICU (yes)", "Not in ICU")]
d_icu[, subgroup := "ICU at t0"]
sg_icu <- d_icu[, .(mean_CATE=mean(CATE,na.rm=TRUE), sd_CATE=sd(CATE,na.rm=TRUE), n=.N), by=.(subgroup,level)]
sg_results$icu <- sg_icu

# By age quartile
d_age <- anal_cate_dt[!is.na(age_quartile)]
d_age[, level := as.character(age_quartile)]
d_age[, subgroup := "Age quartile"]
sg_age <- d_age[, .(mean_CATE=mean(CATE,na.rm=TRUE), sd_CATE=sd(CATE,na.rm=TRUE), n=.N), by=.(subgroup,level)]
sg_results$age <- sg_age

# By vancomycin co-exposure
d_vanco <- anal_cate_dt[!is.na(concomitant_vancomycin)]
d_vanco[, level := fifelse(concomitant_vancomycin==1, "Vancomycin (yes)", "No vancomycin")]
d_vanco[, subgroup := "Vancomycin"]
sg_vanco <- d_vanco[, .(mean_CATE=mean(CATE,na.rm=TRUE), sd_CATE=sd(CATE,na.rm=TRUE), n=.N), by=.(subgroup,level)]
sg_results$vanco <- sg_vanco

sg_all <- rbindlist(sg_results, fill=TRUE)
sg_all[, se_CATE := sd_CATE/sqrt(n)]
sg_all[, lower := mean_CATE - 1.96*se_CATE]
sg_all[, upper := mean_CATE + 1.96*se_CATE]

p_sg <- ggplot(sg_all, aes(x=mean_CATE, y=reorder(level, mean_CATE))) +
  geom_point(size=3, aes(color=subgroup)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper, color=subgroup), height=0.3) +
  geom_vline(xintercept=0, color="black", linewidth=0.5) +
  geom_vline(xintercept=mean(anal_cate_dt$CATE, na.rm=TRUE),
             color="red", linetype="dashed") +
  geom_text(aes(label=sprintf("%.3f\n(n=%d)", mean_CATE, n)),
            nudge_x=0.018, size=2.8) +
  scale_color_brewer(palette="Set1", name="Subgroup") +
  labs(title="Mean CATE by Key Subgroups",
       subtitle="Red dashed: overall ATE — forest-estimated",
       x="Mean Estimated CATE (Risk Difference)",
       y="Subgroup") +
  THEME_BASE

ggsave(file.path(FIGURES_DIR, "cate_subgroups.png"), p_sg,
       width=9, height=6, dpi=150)
cat("Saved: figures/cate_subgroups.png\n")

###############################################################################
# Table: Summary statistics by treatment
###############################################################################
cat("\n--- Summary table by treatment ---\n")

# Already saved as table1.csv, raw_outcomes.csv - save ATE summary table
ate_summary <- data.frame(
  Parameter = c(
    "N (total)",
    "N (ASP — nafcillin/oxacillin)",
    "N (Cefazolin)",
    "AKI at 7 days — ASP, n (%)",
    "AKI at 7 days — CEF, n (%)",
    "AKI at 7 days — Crude OR",
    "AKI at 7 days — Crude p-value",
    "ATE (global, overlap-weighted)",
    "ATE 95% CI",
    "ATE p-value",
    "In-hospital mortality — ASP",
    "In-hospital mortality — CEF"
  ),
  Value = c(
    nrow(anal),
    sum(anal$W==1),
    sum(anal$W==0),
    sprintf("%d (%.1f%%)", sum(anal$W==1 & anal$AKI_7d==1, na.rm=TRUE),
            100*mean(anal$W==1 & anal$AKI_7d==1, na.rm=TRUE)/mean(anal$W==1)),
    sprintf("%d (%.1f%%)", sum(anal$W==0 & anal$AKI_7d==1, na.rm=TRUE),
            100*mean(anal$W==0 & anal$AKI_7d==1, na.rm=TRUE)/mean(anal$W==0)),
    sprintf("%.2f", 1.87),
    "0.004",
    sprintf("%.3f", model_objs$ate_full["estimate"]),
    sprintf("[%.3f, %.3f]",
            model_objs$ate_full["estimate"] - 1.96*model_objs$ate_full["std.err"],
            model_objs$ate_full["estimate"] + 1.96*model_objs$ate_full["std.err"]),
    sprintf("%.4f", 2*pnorm(-abs(model_objs$ate_full["estimate"]/model_objs$ate_full["std.err"]))),
    sprintf("%d (%.1f%%)", sum(anal$W==1 & anal$died_inhosp==1, na.rm=TRUE),
            100*mean(anal$died_inhosp[anal$W==1]==1, na.rm=TRUE)),
    sprintf("%d (%.1f%%)", sum(anal$W==0 & anal$died_inhosp==1, na.rm=TRUE),
            100*mean(anal$died_inhosp[anal$W==0]==1, na.rm=TRUE))
  )
)

fwrite(ate_summary, file.path(TABLES_DIR, "ate_summary.csv"))
cat("Saved: tables/ate_summary.csv\n")

###############################################################################
# Summary of all output files
###############################################################################
cat("\n--- Summary of all generated files ---\n")

fig_files <- list.files(FIGURES_DIR, pattern="\\.png$", full.names=FALSE)
cat("Figures:\n")
for(f in fig_files) cat(sprintf("  %s\n", f))

tab_files <- list.files(TABLES_DIR, pattern="\\.csv$", full.names=FALSE)
cat("Tables:\n")
for(f in tab_files) cat(sprintf("  %s\n", f))

cat("\n=== 08_generate_tables_figures.R COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
