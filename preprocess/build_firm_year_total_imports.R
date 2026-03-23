###############################################################################
# preprocess/build_firm_year_total_imports.R
#
# PURPOSE
#   Compute total import value per firm-year from customs data.
#   Used as the import-cost component of the Leontief cost denominator in
#   compute_b_loop.R. Total imports (not separated by end-use category) are
#   used as a proxy for intermediate import costs. Capital goods imports are
#   included, which slightly overstates the denominator — a conservative
#   approximation that keeps A row sums further below 1.
#
# INPUT
#   {PROC_DATA}/df_trade.RData
#
# OUTPUT
#   {PROC_DATA}/firm_year_total_imports.RData
#     firm_year_total_imports : data frame (vat, year, total_imports)
#
# RUNS ON: RMD
###############################################################################

# ── Paths ─────────────────────────────────────────────────────────────────────
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

load(file.path(PROC_DATA, "df_trade.RData"))

firm_year_total_imports <- df_trade %>%
  filter(flow == "I") %>%
  group_by(vat_ano, year) %>%
  summarise(total_imports = sum(cn_value, na.rm = TRUE), .groups = "drop") %>%
  rename(vat = vat_ano)

save(firm_year_total_imports,
     file = file.path(PROC_DATA, "firm_year_total_imports.RData"))

cat("Saved firm_year_total_imports:", nrow(firm_year_total_imports), "rows\n")
cat("  Years:", min(firm_year_total_imports$year), "-",
    max(firm_year_total_imports$year), "\n")
cat("  Firms:", n_distinct(firm_year_total_imports$vat), "\n")
