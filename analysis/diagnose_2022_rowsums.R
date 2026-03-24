###############################################################################
# analysis/diagnose_2022_rowsums.R
#
# PURPOSE
#   1. Run K=50 vs K=500 Neumann diagnostic for 2020 and 2021
#   2. Investigate why 2022 has 51% of firms with row sum > 0.99
#      and at least one firm with row sum = 1.0
#
# RUNS ON: RMD
###############################################################################

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

# ── Helper: build A and cost_vec for a given year ────────────────────────────
build_year <- function(t, b2b_all, accounts_all, imports_all, eutl_all) {
  b2b_t      <- b2b_all %>% filter(year == t)
  eutl_t     <- eutl_all %>% filter(year == t)
  accounts_t <- accounts_all %>% filter(year == t)
  imports_t  <- imports_all %>% filter(year == t)

  all_vats <- sort(unique(c(b2b_t$vat_i_ano, b2b_t$vat_j_ano)))
  N <- length(all_vats)
  vat_idx <- setNames(seq_len(N), all_vats)

  b2b_agg <- b2b_t %>%
    group_by(vat_i_ano, vat_j_ano) %>%
    summarise(sales = sum(corr_sales_ij, na.rm = TRUE), .groups = "drop") %>%
    filter(sales > 0)

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

  row_i <- vat_idx[b2b_agg$vat_j_ano]
  col_j <- vat_idx[b2b_agg$vat_i_ano]
  ok_ij <- !is.na(row_i) & !is.na(col_j)
  A <- sparseMatrix(
    i = row_i[ok_ij], j = col_j[ok_ij],
    x = b2b_agg$sales[ok_ij] / cost_vec[row_i[ok_ij]],
    dims = c(N, N)
  )

  eps <- rep(0, N)
  ets_idx <- match(eutl_t$vat, all_vats)
  ok_ets <- !is.na(ets_idx) & !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  eps[ets_idx[ok_ets]] <- eutl_t$emissions[ok_ets] / cost_vec[ets_idx[ok_ets]]

  list(A = A, cost_vec = cost_vec, eps = eps, all_vats = all_vats, N = N,
       b2b_rowsums = b2b_rowsums, accounts_t = accounts_t, imports_t = imports_t,
       eutl_t = eutl_t, vat_idx = vat_idx)
}

# ── Helper: run K=50 vs K=500 diagnostic ────────────────────────────────────
run_neumann_diag <- function(yr) {
  cat(sprintf("\n── Year %d: K=50 vs K=500 ─────────────────────────────────\n", yr))
  env <- build_year(yr, b2b_all, accounts_all, imports_all, eutl_all)
  rs <- rowSums(env$A)
  cat(sprintf("  N firms: %d | Max row sum: %.7f\n", env$N, max(rs)))

  m_50 <- as.numeric(env$eps); term <- as.numeric(env$eps)
  for (k in 1:50) { term <- as.numeric(env$A %*% term); m_50 <- m_50 + term }
  m_500 <- m_50; term_500 <- term
  for (k in 51:500) { term_500 <- as.numeric(env$A %*% term_500); m_500 <- m_500 + term_500 }

  up_50  <- pmax(env$cost_vec * (m_50 - env$eps), 0)
  up_500 <- pmax(env$cost_vec * (m_500 - env$eps), 0)
  has_up <- up_500 > 0
  rel_diff <- ifelse(has_up, (up_500 - up_50) / up_500, 0)

  cat(sprintf("  Aggregate gap: %.4f%%\n", 100 * (sum(up_500) - sum(up_50)) / sum(up_500)))
  cat(sprintf("  Median firm gap: %.4f%% | 99th: %.4f%% | Max: %.2f%%\n",
              100 * median(rel_diff[has_up]),
              100 * quantile(rel_diff[has_up], 0.99),
              100 * max(rel_diff)))
  cat(sprintf("  Firms with >1%% gap: %d\n", sum(rel_diff > 0.01)))
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: K=50 vs K=500 for 2020 and 2021
# ══════════════════════════════════════════════════════════════════════════════
cat("══════════════════════════════════════════════════════════════\n")
cat("PART 1: Neumann convergence diagnostic for 2020 and 2021\n")
cat("══════════════════════════════════════════════════════════════\n")

run_neumann_diag(2020)
run_neumann_diag(2021)

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: Investigate 2022 row sums
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n══════════════════════════════════════════════════════════════\n")
cat("PART 2: Investigating 2022 row sums\n")
cat("══════════════════════════════════════════════════════════════\n")

# ── 2a. Coverage comparison: 2020 vs 2021 vs 2022 ───────────────────────────
cat("\n── 2a. Data coverage across years ─────────────────────────\n")

for (t in 2020:2022) {
  b2b_t <- b2b_all %>% filter(year == t)
  b2b_vats <- unique(c(b2b_t$vat_i_ano, b2b_t$vat_j_ano))
  acc_t <- accounts_all %>% filter(year == t)
  imp_t <- imports_all %>% filter(year == t)

  in_acc <- b2b_vats %in% acc_t$vat
  acc_matched <- acc_t %>% filter(vat %in% b2b_vats)
  has_wb <- acc_matched$wage_bill > 0
  in_imp <- b2b_vats %in% imp_t$vat[imp_t$total_imports > 0]

  cat(sprintf("\n  Year %d:\n", t))
  cat(sprintf("    B2B firms: %d\n", length(b2b_vats)))
  cat(sprintf("    In annual accounts: %d (%.1f%%)\n",
              sum(in_acc), 100 * mean(in_acc)))
  cat(sprintf("    With wage_bill > 0: %d (%.1f%% of matched)\n",
              sum(has_wb), 100 * mean(has_wb)))
  cat(sprintf("    With imports > 0: %d (%.1f%%)\n",
              sum(in_imp), 100 * mean(in_imp)))
  cat(sprintf("    Accounts firm-years total: %d\n", nrow(acc_t)))
}

# ── 2b. Build 2022, inspect firms with row sum >= 1 and > 0.99 ─────────────
cat("\n── 2b. Firms with row sum = 1.0 in 2022 ───────────────────\n")

env22 <- build_year(2022, b2b_all, accounts_all, imports_all, eutl_all)
rs22 <- rowSums(env22$A)

# Firms with row sum >= 0.9999999 (effectively 1.0)
exact1 <- which(rs22 >= 1.0 - 1e-7)
cat(sprintf("  Firms with row sum >= 0.9999999: %d\n\n", length(exact1)))

if (length(exact1) > 0) {
  for (i in exact1[1:min(10, length(exact1))]) {
    v <- env22$all_vats[i]
    dom_b2b <- 0
    idx_b <- which(env22$b2b_rowsums$vat_j_ano == v)
    if (length(idx_b)) dom_b2b <- env22$b2b_rowsums$domestic_inputs[idx_b]

    acc_row <- accounts_all %>% filter(year == 2022, vat == v)
    wb <- if (nrow(acc_row)) acc_row$wage_bill[1] else NA

    imp_row <- imports_all %>% filter(year == 2022, vat == v)
    imp_val <- if (nrow(imp_row) && !is.na(imp_row$total_imports[1])) imp_row$total_imports[1] else 0

    eutl_row <- eutl_all %>% filter(year == 2022, vat == v)
    permit <- if (nrow(eutl_row) && !is.na(eutl_row$emissions[1]) && eutl_row$emissions[1] > 0)
      eutl_row$emissions[1] * CARBON_PRICE["2022"] else 0

    cat(sprintf("  VAT: %s\n", substr(v, 1, 16)))
    cat(sprintf("    Row sum: %.10f\n", rs22[i]))
    cat(sprintf("    Domestic B2B: %15.0f\n", dom_b2b))
    cat(sprintf("    Wage bill:    %15.0f  (in accounts: %s)\n",
                ifelse(is.na(wb), 0, wb), ifelse(nrow(acc_row) > 0, "YES", "NO")))
    cat(sprintf("    Imports:      %15.0f\n", imp_val))
    cat(sprintf("    Permit cost:  %15.0f\n", permit))
    cat(sprintf("    Cost vec:     %15.0f\n", env22$cost_vec[i]))
    cat(sprintf("    B2B / cost:   %.10f\n\n", dom_b2b / env22$cost_vec[i]))
  }
}

# ── 2c. Profile of firms with row sum > 0.99 in 2022 vs 2021 ────────────────
cat("── 2c. Row sum > 0.99: 2022 vs 2021 comparison ───────────\n")

env21 <- build_year(2021, b2b_all, accounts_all, imports_all, eutl_all)
rs21 <- rowSums(env21$A)

cat(sprintf("\n  2021: %d / %d firms with rs > 0.99 (%.1f%%)\n",
            sum(rs21 > 0.99), env21$N, 100 * mean(rs21 > 0.99)))
cat(sprintf("  2022: %d / %d firms with rs > 0.99 (%.1f%%)\n",
            sum(rs22 > 0.99), env22$N, 100 * mean(rs22 > 0.99)))

# Among 2022 firms with rs > 0.99, what fraction have wage_bill = 0 or missing?
high22_vats <- env22$all_vats[rs22 > 0.99]
acc22 <- accounts_all %>% filter(year == 2022)
high22_in_acc <- high22_vats %in% acc22$vat
high22_acc <- acc22 %>% filter(vat %in% high22_vats)
high22_wb0 <- sum(high22_acc$wage_bill == 0 | is.na(high22_acc$wage_bill))

cat(sprintf("\n  2022 firms with rs > 0.99 (n = %d):\n", length(high22_vats)))
cat(sprintf("    In accounts: %d (%.1f%%)\n",
            sum(high22_in_acc), 100 * mean(high22_in_acc)))
cat(sprintf("    wage_bill = 0 or NA: %d (%.1f%% of matched)\n",
            high22_wb0, 100 * high22_wb0 / sum(high22_in_acc)))
cat(sprintf("    With imports > 0: %d (%.1f%%)\n",
            sum(high22_vats %in% imports_all$vat[imports_all$year == 2022 & imports_all$total_imports > 0]),
            100 * mean(high22_vats %in% imports_all$vat[imports_all$year == 2022 & imports_all$total_imports > 0])))

# Wage bill distribution for rs > 0.99 firms
cat("\n  Wage bill distribution for rs > 0.99 firms (2022):\n")
wb_high22 <- high22_acc$wage_bill
cat(sprintf("    Min: %.0f | Median: %.0f | Mean: %.0f | Max: %.0f\n",
            min(wb_high22, na.rm=T), median(wb_high22, na.rm=T),
            mean(wb_high22, na.rm=T), max(wb_high22, na.rm=T)))

# Domestic B2B distribution for rs > 0.99 firms
dom_b2b_high22 <- sapply(high22_vats, function(v) {
  idx <- which(env22$b2b_rowsums$vat_j_ano == v)
  if (length(idx)) env22$b2b_rowsums$domestic_inputs[idx] else 0
})
cat(sprintf("  Domestic B2B distribution for rs > 0.99 firms (2022):\n"))
cat(sprintf("    Min: %.0f | Median: %.0f | Mean: %.0f | Max: %.0f\n",
            min(dom_b2b_high22), median(dom_b2b_high22),
            mean(dom_b2b_high22), max(dom_b2b_high22)))

# Ratio wage_bill / domestic_b2b for these firms
ratio_wb_b2b <- wb_high22 / dom_b2b_high22[match(high22_acc$vat, high22_vats)]
cat(sprintf("  wage_bill / domestic_B2B ratio for rs > 0.99 (2022):\n"))
cat(sprintf("    Min: %.6f | Median: %.6f | Mean: %.6f | Max: %.6f\n",
            min(ratio_wb_b2b, na.rm=T), median(ratio_wb_b2b, na.rm=T),
            mean(ratio_wb_b2b, na.rm=T), max(ratio_wb_b2b, na.rm=T)))

cat("\n══════════════════════════════════════════════════════════════\n")
cat("Done.\n")
