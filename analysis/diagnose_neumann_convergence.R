###############################################################################
# analysis/diagnose_neumann_convergence.R
#
# PURPOSE
#   Compare Neumann series at K=50 vs K=500 on full data to quantify how much
#   upstream emission mass is missed by stopping at 50 iterations.
#   Uses a single draw (b=1) and a single year (worst-case) for efficiency.
#
# OUTPUT
#   Prints diagnostics to console. No files saved.
#
# RUNS ON: RMD (full data) or local (downsampled)
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
library(Matrix)

CARBON_PRICE <- c(
  "2005" = 25.29, "2006" = 21.53, "2007" =  0.86, "2008" = 25.74,
  "2009" = 18.41, "2010" = 18.98, "2011" = 18.08, "2012" =  9.49,
  "2013" =  5.94, "2014" =  7.89, "2015" =  8.52, "2016" =  5.92,
  "2017" =  6.63, "2018" = 18.55, "2019" = 27.84, "2020" = 27.94,
  "2021" = 62.25, "2022" = 85.51
)

# ── Load data ────────────────────────────────────────────────────────────────
cat("── Loading data ─────────────────────────────────────────────\n")

load(file.path(PROC_DATA, "b2b_selected_sample.RData"))
b2b_all <- df_b2b_selected_sample
rm(df_b2b_selected_sample)

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl_all <- firm_year_belgian_euets %>% select(vat, year, emissions)
rm(firm_year_belgian_euets)

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts_all <- df_annual_accounts_selected_sample_key_variables %>%
  select(vat, year, wage_bill) %>%
  mutate(wage_bill = pmax(coalesce(wage_bill, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

load(file.path(PROC_DATA, "firm_year_total_imports.RData"))
if ("vat_ano" %in% names(firm_year_total_imports))
  firm_year_total_imports <- rename(firm_year_total_imports, vat = vat_ano)
imports_all <- firm_year_total_imports %>% select(vat, year, total_imports)
rm(firm_year_total_imports)

cat("  Data loaded.\n\n")

# ── Run diagnostic for multiple years ────────────────────────────────────────
DIAG_YEARS <- c(2005, 2010, 2015, 2022)

results <- list()

for (t in DIAG_YEARS) {
  cat(sprintf("── Year %d ──────────────────────────────────────────────────\n", t))

  b2b_t      <- b2b_all %>% filter(year == t)
  eutl_t     <- eutl_all %>% filter(year == t)
  accounts_t <- accounts_all %>% filter(year == t)
  imports_t  <- imports_all %>% filter(year == t)

  all_vats <- sort(unique(c(b2b_t$vat_i_ano, b2b_t$vat_j_ano)))
  N <- length(all_vats)
  vat_idx <- setNames(seq_len(N), all_vats)

  # Aggregate B2B
  b2b_agg <- b2b_t %>%
    group_by(vat_i_ano, vat_j_ano) %>%
    summarise(sales = sum(corr_sales_ij, na.rm = TRUE), .groups = "drop") %>%
    filter(sales > 0)

  # Cost vector
  cost_vec <- rep(1e-6, N)

  b2b_rowsums <- b2b_agg %>%
    group_by(vat_j_ano) %>%
    summarise(domestic_inputs = sum(sales), .groups = "drop")
  buyer_idx <- vat_idx[b2b_rowsums$vat_j_ano]
  ok_buyer <- !is.na(buyer_idx)
  cost_vec[buyer_idx[ok_buyer]] <- cost_vec[buyer_idx[ok_buyer]] +
    b2b_rowsums$domestic_inputs[ok_buyer]

  acc_idx <- match(accounts_t$vat, all_vats)
  ok_acc <- !is.na(acc_idx) & accounts_t$wage_bill > 0
  cost_vec[acc_idx[ok_acc]] <- cost_vec[acc_idx[ok_acc]] + accounts_t$wage_bill[ok_acc]

  imp_idx <- match(imports_t$vat, all_vats)
  ok_imp <- !is.na(imp_idx) & !is.na(imports_t$total_imports) & imports_t$total_imports > 0
  cost_vec[imp_idx[ok_imp]] <- cost_vec[imp_idx[ok_imp]] + imports_t$total_imports[ok_imp]

  ets_cost_idx <- match(eutl_t$vat, all_vats)
  ok_ets_c <- !is.na(ets_cost_idx) & !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  cost_vec[ets_cost_idx[ok_ets_c]] <- cost_vec[ets_cost_idx[ok_ets_c]] +
    eutl_t$emissions[ok_ets_c] * CARBON_PRICE[as.character(t)]

  # Build A
  row_i <- vat_idx[b2b_agg$vat_j_ano]
  col_j <- vat_idx[b2b_agg$vat_i_ano]
  ok_ij <- !is.na(row_i) & !is.na(col_j)
  A <- sparseMatrix(
    i = row_i[ok_ij], j = col_j[ok_ij],
    x = b2b_agg$sales[ok_ij] / cost_vec[row_i[ok_ij]],
    dims = c(N, N)
  )

  rs <- rowSums(A)
  cat(sprintf("  N firms: %d\n", N))
  cat(sprintf("  Max row sum: %.7f\n", max(rs)))
  cat(sprintf("  Firms with row sum > 0.99: %d (%.1f%%)\n",
              sum(rs > 0.99), 100 * mean(rs > 0.99)))
  cat(sprintf("  Firms with row sum > 0.999: %d\n", sum(rs > 0.999)))

  # ETS emission intensities (no imputation — just ETS for diagnostic)
  eps <- rep(0, N)
  ets_idx <- match(eutl_t$vat, all_vats)
  ok_ets <- !is.na(ets_idx) & !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  eps[ets_idx[ok_ets]] <- eutl_t$emissions[ok_ets] / cost_vec[ets_idx[ok_ets]]

  # Run Neumann to K=50
  m_50 <- as.numeric(eps)
  term <- as.numeric(eps)
  for (k in 1:50) {
    term <- as.numeric(A %*% term)
    m_50 <- m_50 + term
  }

  # Continue to K=500
  m_500 <- m_50
  term_500 <- term
  for (k in 51:500) {
    term_500 <- as.numeric(A %*% term_500)
    m_500 <- m_500 + term_500
  }

  # Compare upstream emissions
  upstream_50  <- pmax(cost_vec * (m_50 - eps), 0)
  upstream_500 <- pmax(cost_vec * (m_500 - eps), 0)

  diff <- upstream_500 - upstream_50
  has_up <- upstream_500 > 0
  rel_diff <- ifelse(has_up, diff / upstream_500, 0)

  total_50  <- sum(upstream_50)
  total_500 <- sum(upstream_500)

  cat(sprintf("\n  --- K=50 vs K=500 ---\n"))
  cat(sprintf("  Total upstream (K=50):  %.0f tonnes\n", total_50))
  cat(sprintf("  Total upstream (K=500): %.0f tonnes\n", total_500))
  cat(sprintf("  Aggregate gap: %.4f%%\n", 100 * (total_500 - total_50) / total_500))
  cat(sprintf("  Firm-level relative gap:\n"))
  cat(sprintf("    Median: %.4f%%\n", 100 * median(rel_diff[has_up])))
  cat(sprintf("    90th pctile: %.4f%%\n", 100 * quantile(rel_diff[has_up], 0.9)))
  cat(sprintf("    99th pctile: %.4f%%\n", 100 * quantile(rel_diff[has_up], 0.99)))
  cat(sprintf("    Max: %.2f%% (%.0f tonnes)\n", 100 * max(rel_diff), max(diff)))
  cat(sprintf("  Firms with >1%% gap: %d\n\n", sum(rel_diff > 0.01)))

  results[[as.character(t)]] <- data.frame(
    year = t, N = N, max_rowsum = max(rs),
    n_rs_gt99 = sum(rs > 0.99), n_rs_gt999 = sum(rs > 0.999),
    total_up_50 = total_50, total_up_500 = total_500,
    agg_gap_pct = 100 * (total_500 - total_50) / total_500,
    median_gap_pct = 100 * median(rel_diff[has_up]),
    p99_gap_pct = 100 * quantile(rel_diff[has_up], 0.99),
    max_gap_pct = 100 * max(rel_diff),
    n_gt1pct = sum(rel_diff > 0.01)
  )
}

cat("\n══════════════════════════════════════════════════════════════\n")
cat("Summary across years:\n")
print(bind_rows(results))
cat("══════════════════════════════════════════════════════════════\n")
