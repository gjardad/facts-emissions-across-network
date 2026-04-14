###############################################################################
# analysis/imjv_topx_comparison.R
#
# PURPOSE
#   For the 56 IMJV TP firm-years, compare three classification rules at
#   different top-X% cutoffs:
#     (a) Top X% of proxy_mean within CRF-year
#     (b) Top X% of revenue within CRF-year
#     (c) Top X% of proxy_mean within CRF-year-SIZE BIN
#         (size bins = revenue deciles within CRF-year)
#
# INPUT
#   {RAW_DATA}/IMJV/imjv_co2_lucht.tsv
#   {RAW_DATA}/IMJV/crosswalk/imjv_cbb_vat_ano.csv
#   {RAW_DATA}/IMJV/crosswalk/non_ets_cbbs.csv
#   {PROC_DATA}/allocation_glo_balanced/alloc_YYYY.RData
#   {PROC_DATA}/extensive_margin_calibration.RData
#   {PROC_DATA}/deployment_panel.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData (FULL)
#
# RUNS ON: RMD (needs full annual accounts and deployment panel)
###############################################################################

# -- Paths -------------------------------------------------------------------
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

cat("===================================================================\n")
cat("  IMJV top-X% classification comparison\n")
cat("===================================================================\n\n")


# =============================================================================
# 1. Load IMJV → non-ETS firm-years
# =============================================================================
cat("-- Loading IMJV data ------------------------------------------------\n")

imjv <- read.delim(file.path(RAW_DATA, "IMJV", "imjv_co2_lucht.tsv"),
                    stringsAsFactors = FALSE,
                    colClasses = c(cbb_number = "character", cbb_current = "character"))
xwalk <- read.csv(file.path(RAW_DATA, "IMJV", "crosswalk", "imjv_cbb_vat_ano.csv"),
                   stringsAsFactors = FALSE,
                   colClasses = c(cbb = "character", vat_ano = "character"))
imjv <- imjv %>%
  left_join(xwalk %>% select(cbb, vat_ano), by = c("cbb_number" = "cbb")) %>%
  filter(!is.na(vat_ano))
imjv_vy <- imjv %>%
  group_by(vat_ano, year) %>%
  summarise(imjv_co2_t = sum(emission_kg, na.rm = TRUE) / 1000, .groups = "drop")
vat_names <- imjv %>%
  group_by(vat_ano) %>%
  summarise(firm_name = first(firm_name), .groups = "drop")

non_ets <- read.csv(file.path(RAW_DATA, "IMJV", "crosswalk", "non_ets_cbbs.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(cbb = "character", vat_ano = "character"))
non_ets_vats <- unique(non_ets$vat_ano)


# =============================================================================
# 2. Load allocation → proxy values for all imputed firms
# =============================================================================
cat("-- Loading allocation files -----------------------------------------\n")

alloc_list <- list()
for (yr in 2005:2021) {
  load(file.path(PROC_DATA, "allocation_glo_balanced", sprintf("alloc_%d.RData", yr)))
  alloc_list[[as.character(yr)]] <- year_firms
  rm(year_firms)
}
alloc <- bind_rows(alloc_list)
rm(alloc_list)

alloc_imp <- alloc %>%
  filter(source == "imputed") %>%
  select(vat, year, crf_group, scope1, p_i, proxy_mean_i)

cat("  Imputed firm-years:", nrow(alloc_imp), "\n")


# =============================================================================
# 3. Load revenue from annual accounts (FULL on RMD)
# =============================================================================
cat("-- Loading annual accounts ------------------------------------------\n")

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  select(vat, year, revenue)
rm(df_annual_accounts_selected_sample_key_variables)
cat("  Accounts rows:", nrow(accounts), "\n")

# Join revenue to imputed firms
alloc_imp <- alloc_imp %>%
  left_join(accounts, by = c("vat", "year"))

n_with_rev <- sum(!is.na(alloc_imp$revenue))
cat("  Imputed firm-years with revenue:", n_with_rev,
    "(", round(100 * n_with_rev / nrow(alloc_imp), 1), "%)\n\n")


# =============================================================================
# 4. Compute percentiles within CRF-year
# =============================================================================
cat("-- Computing percentiles --------------------------------------------\n")

# (a) Proxy percentile within CRF-year
alloc_imp <- alloc_imp %>%
  group_by(crf_group, year) %>%
  mutate(
    n_crf_yr = n(),
    proxy_pctile_crf = rank(proxy_mean_i, ties.method = "average") / n()
  ) %>%
  ungroup()

# (b) Revenue percentile within CRF-year (among firms with non-NA revenue)
alloc_imp <- alloc_imp %>%
  group_by(crf_group, year) %>%
  mutate(
    rev_pctile_crf = {
      r <- rep(NA_real_, n())
      ok <- !is.na(revenue)
      if (sum(ok) > 1) {
        r[ok] <- rank(revenue[ok], ties.method = "average") / sum(ok)
      }
      r
    }
  ) %>%
  ungroup()

# (c) Size bins within CRF-year, then proxy percentile within CRF-year-size bin
alloc_imp <- alloc_imp %>%
  group_by(crf_group, year) %>%
  mutate(
    size_bin = {
      sb <- rep(NA_integer_, n())
      ok <- !is.na(revenue)
      if (sum(ok) >= 10) {
        # Deciles of revenue within CRF-year
        sb[ok] <- as.integer(cut(revenue[ok],
                                  breaks = quantile(revenue[ok],
                                                     probs = seq(0, 1, by = 0.1),
                                                     na.rm = TRUE),
                                  include.lowest = TRUE, labels = FALSE))
      } else if (sum(ok) >= 2) {
        # Too few for deciles; use 2 bins (above/below median)
        sb[ok] <- as.integer(revenue[ok] >= median(revenue[ok], na.rm = TRUE)) + 1L
      }
      sb
    }
  ) %>%
  ungroup()

# Proxy percentile within CRF-year-size bin
alloc_imp <- alloc_imp %>%
  group_by(crf_group, year, size_bin) %>%
  mutate(
    n_bin = n(),
    proxy_pctile_bin = if (n() >= 2 & !all(is.na(size_bin))) {
      rank(proxy_mean_i, ties.method = "average") / n()
    } else {
      rep(NA_real_, n())
    }
  ) %>%
  ungroup()

cat("  Done.\n\n")


# =============================================================================
# 5. Identify IMJV TP firm-years
# =============================================================================

tp <- imjv_vy %>%
  filter(vat_ano %in% non_ets_vats, year >= 2005, year <= 2021, imjv_co2_t > 0) %>%
  inner_join(alloc_imp, by = c("vat_ano" = "vat", "year")) %>%
  filter(scope1 > 0) %>%
  left_join(vat_names, by = "vat_ano")

cat("TP firm-years:", nrow(tp), "\n")
cat("With revenue:", sum(!is.na(tp$revenue)), "\n")
cat("With size_bin:", sum(!is.na(tp$size_bin)), "\n\n")


# =============================================================================
# 6. Top-X% comparison table
# =============================================================================

cutoffs <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)

# For a fair comparison, restrict to firm-years that have all three measures
tp_complete <- tp %>% filter(!is.na(revenue), !is.na(size_bin), !is.na(proxy_pctile_bin))
n_complete <- nrow(tp_complete)

cat("===================================================================\n")
cat("  TOP-X% CLASSIFICATION COMPARISON\n")
cat("  (restricted to firm-years with revenue and size bin: n =", n_complete, ")\n")
cat("===================================================================\n\n")

cat(sprintf("%8s %18s %18s %18s\n",
            "Top X%",
            "Proxy/CRF-yr",
            "Revenue/CRF-yr",
            "Proxy/CRF-yr-size"))
cat(strrep("-", 66), "\n")

for (x in cutoffs) {
  threshold <- 1 - x  # percentile threshold (e.g., top 10% = pctile >= 0.90)
  n_proxy_crf  <- sum(tp_complete$proxy_pctile_crf >= threshold)
  n_rev_crf    <- sum(tp_complete$rev_pctile_crf >= threshold, na.rm = TRUE)
  n_proxy_bin  <- sum(tp_complete$proxy_pctile_bin >= threshold, na.rm = TRUE)
  cat(sprintf("%7.0f%% %8d / %-4d %8d / %-4d %8d / %-4d\n",
              x * 100,
              n_proxy_crf, n_complete,
              n_rev_crf, n_complete,
              n_proxy_bin, n_complete))
}

cat("\n")

# Also report on the full 56 (not restricted to complete cases)
cat("===================================================================\n")
cat("  SAME TABLE, UNRESTRICTED (all 56 TP firm-years)\n")
cat("  Revenue/CRF-yr uses only firm-years with revenue (n =",
    sum(!is.na(tp$rev_pctile_crf)), ")\n")
cat("  Proxy/CRF-yr-size uses only firm-years with size bin (n =",
    sum(!is.na(tp$proxy_pctile_bin)), ")\n")
cat("===================================================================\n\n")

cat(sprintf("%8s %18s %18s %18s\n",
            "Top X%",
            "Proxy/CRF-yr",
            "Revenue/CRF-yr",
            "Proxy/CRF-yr-size"))
cat(strrep("-", 66), "\n")

for (x in cutoffs) {
  threshold <- 1 - x
  n_proxy_crf  <- sum(tp$proxy_pctile_crf >= threshold)
  n_rev_crf    <- sum(tp$rev_pctile_crf >= threshold, na.rm = TRUE)
  n_proxy_bin  <- sum(tp$proxy_pctile_bin >= threshold, na.rm = TRUE)
  denom_rev    <- sum(!is.na(tp$rev_pctile_crf))
  denom_bin    <- sum(!is.na(tp$proxy_pctile_bin))
  cat(sprintf("%7.0f%% %8d / %-4d %8d / %-4d %8d / %-4d\n",
              x * 100,
              n_proxy_crf, nrow(tp),
              n_rev_crf, denom_rev,
              n_proxy_bin, denom_bin))
}

cat("\n")

# Detail for proxy/CRF-yr-size: show which size bins the IMJV firms land in
cat("===================================================================\n")
cat("  SIZE BIN DETAIL FOR IMJV TP FIRM-YEARS\n")
cat("===================================================================\n\n")

tp_with_bin <- tp %>% filter(!is.na(size_bin))
cat(sprintf("%-28s %5s %8s %4s %6s %8s %8s\n",
            "Firm", "Year", "CRF", "Bin", "n_bin", "prx_pct", "rev_pct"))
cat(strrep("-", 75), "\n")
tp_sorted <- tp_with_bin %>% arrange(crf_group, year, size_bin)
for (i in seq_len(nrow(tp_sorted))) {
  r <- tp_sorted[i, ]
  cat(sprintf("%-28s %5d %8s %4d %6d %8.3f %8.3f\n",
              substr(r$firm_name, 1, 28), r$year, r$crf_group,
              r$size_bin, r$n_bin,
              r$proxy_pctile_bin,
              ifelse(is.na(r$rev_pctile_crf), NA, r$rev_pctile_crf)))
}

cat("\n===================================================================\n")
cat("  Done.\n")
cat("===================================================================\n")
