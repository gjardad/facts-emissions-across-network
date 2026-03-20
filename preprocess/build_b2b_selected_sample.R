###############################################################################
# preprocess/build_b2b_selected_sample.R
#
# PURPOSE
#   Restrict the full B2B network to transactions where both buyer and seller
#   appear in the annual accounts selected sample. Produces a clean
#   buyer-seller-year panel used by all downstream proxy scripts.
#
# INPUT
#   {PROC_DATA}/df_b2b.RData
#   {PROC_DATA}/annual_accounts_selected_sample.RData
#
# OUTPUT
#   {PROC_DATA}/b2b_selected_sample.RData
#     Contains: df_b2b_selected_sample (data.frame)
#     Columns: vat_i_ano (seller), vat_j_ano (buyer), year, corr_sales_ij
#
# RUNS ON: RMD (requires full B2B data)
# COPIED FROM: inferring_emissions/preprocess/build_b2b_selected_sample.R
###############################################################################

# ── Paths ────────────────────────────────────────────────────────────────────
if (tolower(Sys.info()[["user"]]) == "jardang") {
  REPO_DIR <- "C:/Users/jardang/Documents/facts-emissions-across-network"
} else if (tolower(Sys.info()[["user"]]) == "jota_"){
  REPO_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
                        error = function(e) normalizePath(getwd(), winslash = "/"))
  while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
} else {
  stop("Define REPO_DIR for this user.")
}
source(file.path(REPO_DIR, "paths.R"))

library(tidyverse)
library(dplyr)

# ── Load data ────────────────────────────────────────────────────────────────
load(paste0(PROC_DATA, "/df_b2b.RData"))
load(paste0(PROC_DATA, "/annual_accounts_selected_sample.RData"))

# ── Create B2B with only transactions between firms in selected sample ────────

# Drop original (uncorrected) sales; keep only corrected number
df_b2b <- df_b2b %>%
  select(-sales_ij)

unique_vat_by_year <- df_annual_accounts_selected_sample %>%
  group_by(year) %>%
  summarise(unique_vat = list(unique(vat_ano))) %>%
  pull(unique_vat)

# starts at i=3 because df_annual_accounts_selected_sample covers 2000, 2001
# as well (first two elements of unique_vat_by_year are 2000, 2001)
i <- 2
df_b2b_selected_sample <- c()

for (y in 2002:2022) {
  i <- i + 1

  temp_data <- df_b2b %>%
    filter(year == y) %>%
    filter(vat_i_ano %in% unique_vat_by_year[[i]],
           vat_j_ano %in% unique_vat_by_year[[i]])

  df_b2b_selected_sample <- bind_rows(df_b2b_selected_sample, temp_data)
}

# ── Save ─────────────────────────────────────────────────────────────────────
save(df_b2b_selected_sample, file = paste0(PROC_DATA, "/b2b_selected_sample.RData"))
cat("Saved: b2b_selected_sample.RData\n")
cat("  Rows:", nrow(df_b2b_selected_sample), "\n")
cat("  Years: 2002–2022\n")
