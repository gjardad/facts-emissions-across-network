###############################################################################
# preprocess/build_deployment_panel.R
#
# PURPOSE
#   Build a minimal firm-year panel for deployment (non-ETS) firms covering
#   the firm-year observations present in the B2B data (i.e., those that
#   receive a proxy value from run_subsampled_en.R). Contains only the
#   variables needed to run diagnostics and statistics locally without
#   requiring the full annual accounts or B2B data.
#
#   Copy this file to local 1 once; use it alongside any proxy/emissions
#   panel generated on RMD.
#
# INPUT
#   {PROC_DATA}/deployment_proxy_avg.RData   — defines deployment firm-years
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   {PROC_DATA}/deployment_panel.RData
#     deployment_panel: (vat, year, nace5d, revenue)
#
# RUNS ON: RMD
###############################################################################

# ── Paths ────────────────────────────────────────────────────────────────────
if (tolower(Sys.info()[["user"]]) == "jardang") {
  REPO_DIR <- "C:/Users/jardang/Documents/facts-emissions-across-network"
} else if (tolower(Sys.info()[["user"]]) == "jota_") {
  REPO_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
                        error = function(e) normalizePath(getwd(), winslash = "/"))
  while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
} else {
  stop("Define REPO_DIR for this user.")
}
source(file.path(REPO_DIR, "paths.R"))

library(dplyr)

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading deployment proxy...\n")
load(file.path(PROC_DATA, "deployment_proxy_avg.RData"))

cat("Loading annual accounts...\n")
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))

# ── Build panel ───────────────────────────────────────────────────────────────
deployment_panel <- proxy_avg %>%
  select(vat, year) %>%
  left_join(
    df_annual_accounts_selected_sample_key_variables %>%
      select(vat, year, nace5d, revenue),
    by = c("vat", "year")
  )

n_missing_nace    <- sum(is.na(deployment_panel$nace5d))
n_missing_revenue <- sum(is.na(deployment_panel$revenue))

cat("Deployment panel:", nrow(deployment_panel), "firm-years,",
    n_distinct(deployment_panel$vat), "firms\n")
cat("  Missing nace5d: ", n_missing_nace,
    sprintf("(%.1f%%)\n", 100 * n_missing_nace / nrow(deployment_panel)))
cat("  Missing revenue:", n_missing_revenue,
    sprintf("(%.1f%%)\n", 100 * n_missing_revenue / nrow(deployment_panel)))

# ── Save ─────────────────────────────────────────────────────────────────────
OUT_PATH <- file.path(PROC_DATA, "deployment_panel.RData")
save(deployment_panel, file = OUT_PATH)
cat("\nSaved:", OUT_PATH, "\n")
