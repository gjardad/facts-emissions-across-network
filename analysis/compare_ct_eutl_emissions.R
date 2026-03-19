###############################################################################
# analysis/compare_ct_eutl_emissions.R
#
# PURPOSE
#   Compare Climate TRACE emission estimates against EUTL verified emissions
#   at the installation-year level for matched facilities.
#
# INPUTS
#   {PROC_DATA}/ct_eutl_match.RData         (from match_ct_to_eutl.R)
#   {RAW_DATA}/EUTL/Oct_2024_version/compliance.csv
#
# RUNS ON: local 1
###############################################################################

# ── Paths ────────────────────────────────────────────────────────────────────
REPO_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) normalizePath(getwd(), winslash = "/")
)
while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
source(file.path(REPO_DIR, "paths.R"))

library(dplyr)

# ── Load data ────────────────────────────────────────────────────────────────
load(file.path(PROC_DATA, "ct_eutl_match.RData"))  # ct, ct_eutl_match, be_inst

compliance <- read.csv(
  file.path(RAW_DATA, "EUTL", "Oct_2024_version", "compliance.csv"),
  stringsAsFactors = FALSE
)

# EUTL verified emissions at installation-year level
eutl <- compliance %>%
  filter(!is.na(verified)) %>%
  select(installation_id, year, eutl_emissions = verified)


# ── Merge CT and EUTL at installation-year level ─────────────────────────────

# CT emissions by matched installation-year
ct_inst <- ct %>%
  inner_join(ct_eutl_match %>% select(source_id, installation_id),
             by = "source_id") %>%
  group_by(installation_id, year) %>%
  summarise(ct_emissions = sum(ct_emissions, na.rm = TRUE), .groups = "drop")

# Join
compared <- ct_inst %>%
  inner_join(eutl, by = c("installation_id", "year")) %>%
  filter(ct_emissions > 0, eutl_emissions > 0)

cat("Matched installation-years with both CT > 0 and EUTL > 0:", nrow(compared), "\n")
cat("Unique installations:", n_distinct(compared$installation_id), "\n")
cat("Year range:", range(compared$year), "\n\n")


# ── Comparison metrics ───────────────────────────────────────────────────────

y    <- compared$eutl_emissions
yhat <- compared$ct_emissions

# RMSE and nRMSE
rmse <- sqrt(mean((y - yhat)^2))
nrmse_sd <- rmse / sd(y)

# Correlation
pearson_r   <- cor(y, yhat, use = "complete.obs")
spearman_rho <- cor(y, yhat, method = "spearman", use = "complete.obs")

# Log-scale correlation
pearson_log <- cor(log(y), log(yhat), use = "complete.obs")

# Ratio CT / EUTL
ratio <- yhat / y
median_ratio <- median(ratio)
mean_ratio   <- mean(ratio)

# Median absolute percentage difference
median_apd <- median(abs(yhat - y) / y)

cat("-- Comparison metrics --\n")
cat(sprintf("N observations:          %d\n", nrow(compared)))
cat(sprintf("RMSE:                    %.0f\n", rmse))
cat(sprintf("nRMSE (sd):              %.3f\n", nrmse_sd))
cat(sprintf("Pearson r:               %.3f\n", pearson_r))
cat(sprintf("Pearson r (log scale):   %.3f\n", pearson_log))
cat(sprintf("Spearman rho:            %.3f\n", spearman_rho))
cat(sprintf("Median CT/EUTL ratio:    %.3f\n", median_ratio))
cat(sprintf("Mean CT/EUTL ratio:      %.3f\n", mean_ratio))
cat(sprintf("Median abs pct diff:     %.3f\n", median_apd))

cat("\n-- EUTL emissions summary --\n")
print(summary(y))

cat("\n-- CT emissions summary --\n")
print(summary(yhat))

cat("\n-- CT/EUTL ratio distribution --\n")
print(summary(ratio))

# By-year breakdown
cat("\n-- By year --\n")
by_year <- compared %>%
  group_by(year) %>%
  summarise(
    n = n(),
    spearman_rho = cor(eutl_emissions, ct_emissions, method = "spearman"),
    median_ratio = median(ct_emissions / eutl_emissions),
    median_apd = median(abs(ct_emissions - eutl_emissions) / eutl_emissions),
    .groups = "drop"
  )
print(as.data.frame(by_year))
