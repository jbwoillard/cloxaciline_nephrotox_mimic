# Feasibility Summary

Date: 2026-03-10 16:18:43.223691

## 1. Data Availability

MIMIC-IV version: 3.1 HOSP tables available: admissions.csv.gz, d_hcpcs.csv.gz, d_icd_diagnoses.csv.gz, d_icd_procedures.csv.gz, d_labitems.csv.gz, diagnoses_icd.csv.gz, drgcodes.csv.gz, emar_detail.csv.gz, emar.csv.gz, hcpcsevents.csv.gz, labevents.csv.gz, microbiologyevents.csv.gz, omr.csv.gz, patients.csv.gz, pharmacy.csv.gz, poe_detail.csv.gz, poe.csv.gz, prescriptions.csv.gz, procedures_icd.csv.gz, provider.csv.gz, services.csv.gz, transfers.csv.gz ICU tables available: caregiver.csv.gz, chartevents.csv.gz, d_items.csv.gz, datetimeevents.csv.gz, icustays.csv.gz, ingredientevents.csv.gz, inputevents.csv.gz, outputevents.csv.gz, procedureevents.csv.gz

## 2. Drug Feasibility

### ASP (Anti-Staphylococcal Penicillins: oxacillin/nafcillin)

-   prescriptions: n_patients = 2535 \| n_rows = 4247
-   pharmacy: n_patients = 2539 \| n_rows = 4239
-   emar: n_patients = 849 \| n_rows = 27769
-   Drugs found: Nafcillin; DiCLOXacillin; Nafcillin Desensitization; Oxacillin; Nafcillin Graded Challenge; Oxacillin Desensitization; nafcillin

### Cefazolin

-   prescriptions: n_patients = 46358 \| n_rows = 69154
-   pharmacy: n_patients = 46339 \| n_rows = 68649
-   emar: n_patients = 29737 \| n_rows = 153866
-   Drugs found: CefazoLIN; CeFAZolin; ceFAZolin; CeFAZolin Desensitization; CefazoLIN Desensitization; Cefazolin; CeFAZolin Duple 2g/50mL 50mL BAG; Cefazolin 2 g; CeFAZolin Graded Challenge; CeFAZolin-Heparin Lock (For HD/Pheresis Catheter); CefazoLIN-Heparin Lock (For HD/Pheresis Catheter); Cefazolin 2 gr; Cefazolin 2 gr or Placebo

## 3. MSSA Feasibility

-   SA blood culture patients (total): 0
-   Confirmed MSSA (oxacillin S): 0

## 4. Decision

**FEASIBLE**: Both ASP (oxacillin/nafcillin) and cefazolin are present. Cloxacillin is NOT present in MIMIC-IV (US database). Will use **ASP (oxacillin + nafcillin) vs cefazolin** design.

Note: Cloxacillin is not used in the US and is absent from MIMIC-IV. The study will use oxacillin and/or nafcillin as anti-staphylococcal penicillin (ASP) comparator, consistent with the conceptual replication approach described in the prompt.

## 5. Tables produced

-   tables/drug_feasibility.csv
-   tables/microbiology_feasibility.csv
