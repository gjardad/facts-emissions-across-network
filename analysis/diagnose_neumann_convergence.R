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
DIAG_YEARS   <- c(2005, 2010, 2015, 2021)
CHECKPOINTS  <- c(50L, 100L, 200L, 300L, 500L)
K_MAX        <- max(CHECKPOINTS)

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
  cat(sprintf("  N firms: %d | Max row sum: %.7f\n", N, max(rs)))

  # ETS emission intensities (no imputation — just ETS for diagnostic)
  eps <- rep(0, N)
  ets_idx <- match(eutl_t$vat, all_vats)
  ok_ets <- !is.na(ets_idx) & !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  eps[ets_idx[ok_ets]] <- eutl_t$emissions[ok_ets] / cost_vec[ets_idx[ok_ets]]

  # Run Neumann with snapshots at each checkpoint
  m    <- as.numeric(eps)
  term <- as.numeric(eps)
  snapshots <- list()

  for (k in seq_len(K_MAX)) {
    term <- as.numeric(A %*% term)
    m    <- m + term
    if (k %in% CHECKPOINTS) {
      upstream_k <- pmax(cost_vec * (m - eps), 0)
      snapshots[[as.character(k)]] <- sum(upstream_k)
      cat(sprintf("    K=%3d done\n", k))
    }
  }

  # Reference = K_MAX
  total_ref <- snapshots[[as.character(K_MAX)]]

  cat(sprintf("\n  --- Cumulative upstream at each K (reference = K=%d) ---\n", K_MAX))
  cat(sprintf("  %-6s  %15s  %10s  %10s\n", "K", "total_upstream", "pct_of_ref", "gap_pct"))

  for (ck in CHECKPOINTS) {
    total_k <- snapshots[[as.character(ck)]]
    pct     <- 100 * total_k / total_ref
    gap     <- 100 * (total_ref - total_k) / total_ref
    cat(sprintf("  K=%-4d  %15.0f  %9.4f%%  %9.4f%%\n", ck, total_k, pct, gap))

    results[[length(results) + 1L]] <- data.frame(
      year = t, K = ck, total_upstream = total_k,
      pct_of_ref = pct, gap_pct = gap
    )
  }
  cat("\n")
}

cat("══════════════════════════════════════════════════════════════\n")
cat("Summary:\n")
print(bind_rows(results))
cat("══════════════════════════════════════════════════════════════\n")
