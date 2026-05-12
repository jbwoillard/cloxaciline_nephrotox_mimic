repos <- "https://cloud.r-project.org"

pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "lubridate",
  "readr",
  "vroom",
  "arrow",
  "duckdb",
  "DBI",
  "dbplyr",
  "janitor",
  "purrr",
  "ggplot2",
  "patchwork",
  "gt",
  "gtsummary",
  "tableone",
  "cobalt",
  "survival",
  "broom",
  "broom.helpers",
  "rmarkdown",
  "knitr",
  "yaml",
  "here",
  "qs",
  "grf",
  "xgboost",
  "glmnet"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  install.packages(to_install, repos = repos, Ncpus = max(1, parallel::detectCores() - 1))
}

cat("\nPackages installés. Vérification:\n")
library(grf)
library(xgboost)
library(glmnet)
packageVersion("grf")
packageVersion("xgboost")
packageVersion("glmnet")
sessionInfo()
