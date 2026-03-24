###############################################################################
# analysis/compute_b_loop.R
#
# PURPOSE
#   For each of B = 50 subsample draws (from run_subsampled_en.R):
#     1. NIR calibration: classify deployment emitters, distribute E_sy_deploy
#     2. Build emission intensity vector for all firms (ETS + imputed)
#     3. Compute upstream embodied emissions via Leontief / Neumann series
#     4. Compute RQ1–RQ5 statistics
#   Average statistics across B draws (point estimates) + cross-B s.d.
#
#   EFFICIENCY: For each year, A is built once and reused across all B draws.
#   The Neumann series m ≈ (I-A)^{-1} ε = ε + Aε + A²ε + ... is run B times
#   per year with different ε vectors (B cheap sparse mat-vec iterations, not
#   B full inversions).
#
#   For Belgian data, the series converges in K ≈ 3–5 iterations because the
#   domestic row sums of A are small (Belgium is highly import-intensive, so
#   most intermediate inputs are imported and absent from B2B).
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData
#   {PROC_DATA}/b2b_selected_sample.RData
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/deployment_panel.RData              (vat, year, nace5d)
#   {PROC_DATA}/nir_calibration_targets.RData       (E_NIR_kt by CRF group × year)
#   {PROC_DATA}/firm_year_total_imports.RData        (vat, year, total_imports)
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv  (NACE 2d → crf_group)
#
# OUTPUT  — copy b_loop_statistics.RData to local 1
#   {PROC_DATA}/b_loop_statistics.RData
#     stats2d_summary : point estimates + cross-B s.d. by NACE 2-digit × year
#     stats5d_summary : same at NACE 5-digit × year
#     all_stats2d     : full B × sector-year draw-level data (keep on RMD)
#     all_stats5d     : same at NACE 5-digit
#     conv_summary    : Neumann convergence diagnostics by year
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
library(Matrix)
library(parallel)
library(doParallel)
library(foreach)

# ── Parameters ────────────────────────────────────────────────────────────────
YEARS         <- 2005:2021  # 2022 excluded: wage_bill data quality issue (see A_ITERATION.md)
NEUMANN_MAXIT <- 50L    # max terms in Neumann series per (year, draw)
NEUMANN_TOL   <- 1e-8   # convergence: max(|term_k|) / (max(|m|) + eps) < tol
MIN_N_STATS   <- 3L     # min firms per sector-year to compute a statistic

# EU ETS annual average carbon price (EUR / tonne CO2), 2005–2022.
# Source: ICAP Allowance Price Explorer, secondary market daily prices.
#   Raw data: RAW_DATA/raw/icap_systems.json (downloaded 2026-03-24)
#   Annual averages computed from daily "Secondary Market" observations.
# Phase 1 (2005–2007) allowances expired worthless → 2007 price ~0.
CARBON_PRICE <- c(
  "2005" = 25.29,
  "2006" = 21.53,
  "2007" =  0.86,   # Phase 1 collapse
  "2008" = 25.74,
  "2009" = 18.41,
  "2010" = 18.98,
  "2011" = 18.08,
  "2012" =  9.49,
  "2013" =  5.94,
  "2014" =  7.89,
  "2015" =  8.52,
  "2016" =  5.92,
  "2017" =  6.63,
  "2018" = 18.55,
  "2019" = 27.84,
  "2020" = 27.94,
  "2021" = 62.25,
  "2022" = 85.51
)
N_CORES_SET <- if (tolower(Sys.info()[["user"]]) == "jardang") {
  40L    # RMD: set explicitly (detectCores may undercount)
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 2L)
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("  B-LOOP: NIR + LEONTIEF/NEUMANN + RQ1-RQ5\n")
cat("  Years:", min(YEARS), "–", max(YEARS),
    "| Neumann tol:", NEUMANN_TOL, "| max iter:", NEUMANN_MAXIT, "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")


# =============================================================================
# SECTION 1: Load data
# =============================================================================
cat("── Loading data ─────────────────────────────────────────────\n")

load(file.path(PROC_DATA, "deployment_proxy_list.RData"))
B       <- length(proxy_list)
N_CORES <- min(B, N_CORES_SET)
cat("  proxy_list: B =", B, "draws\n")
cat("  Cores for draw parallelism:", N_CORES, "\n")

load(file.path(PROC_DATA, "b2b_selected_sample.RData"))
b2b <- df_b2b_selected_sample %>% filter(year %in% YEARS)
rm(df_b2b_selected_sample)
cat("  B2B rows (years >= 2005):", nrow(b2b), "\n")

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)
cat("  EUTL firm-years:", nrow(eutl), "\n")

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, wage_bill) %>%
  mutate(nace2d   = substr(nace5d, 1, 2),
         wage_bill = pmax(coalesce(wage_bill, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)
cat("  Accounts firm-years:", nrow(accounts), "\n")

load(file.path(PROC_DATA, "firm_year_total_imports.RData"))
if ("vat_ano" %in% names(firm_year_total_imports))
  firm_year_total_imports <- rename(firm_year_total_imports, vat = vat_ano)
imports <- firm_year_total_imports %>%
  filter(year %in% YEARS) %>%
  select(vat, year, total_imports)
rm(firm_year_total_imports)
cat("  Import firm-years:", nrow(imports), "\n")

load(file.path(PROC_DATA, "deployment_panel.RData"))

nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(nace2d = "character")
) %>%
  select(nace2d, crf_group)

deploy_nace <- deployment_panel %>%
  mutate(nace2d = substr(nace5d, 1, 2)) %>%
  select(vat, nace2d, nace5d) %>%
  distinct(vat, .keep_all = TRUE) %>%
  left_join(nace_crf, by = "nace2d")
rm(deployment_panel)

cat("\n")


# =============================================================================
# SECTION 2: NIR calibration targets (fixed across draws)
# =============================================================================
cat("── NIR calibration targets ──────────────────────────────────\n")

load(file.path(PROC_DATA, "nir_calibration_targets.RData"))  # nir_targets: (crf_group, year, E_NIR_kt)

# E_ETS aggregated to CRF group × year (using the same NACE → crf_group crosswalk)
E_ETS_group <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group), year %in% YEARS) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions, na.rm = TRUE), .groups = "drop")

E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%         # kt → t CO2
  left_join(E_ETS_group, by = c("crf_group", "year")) %>%
  mutate(E_ETS    = coalesce(E_ETS, 0),
         E_deploy = pmax(E_NIR - E_ETS, 0))

n_floored <- sum(E_deploy_panel$E_NIR < E_deploy_panel$E_ETS, na.rm = TRUE)
if (n_floored > 0)
  cat("  WARNING:", n_floored, "group-years where E_ETS > E_NIR → floored to 0\n")

cat("  E_deploy built:", nrow(E_deploy_panel), "group-years\n\n")


# =============================================================================
# SECTION 3: Helper functions
# =============================================================================

gini <- function(x) {
  x <- x[!is.na(x) & x > 0]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  x <- sort(x)
  2 * sum(x * seq_len(n)) / (n * sum(x)) - (n + 1L) / n
}

pct_ratio <- function(x, p_hi, p_lo) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2L) return(NA_real_)
  q <- quantile(x, c(p_lo, p_hi), names = FALSE)
  if (q[1L] == 0) return(NA_real_)
  q[2L] / q[1L]
}

safe_spearman <- function(x, y, min_n = MIN_N_STATS) {
  keep <- !is.na(x) & !is.na(y) & x > 0 & y > 0
  if (sum(keep) < min_n) return(NA_real_)
  cor(x[keep], y[keep], method = "spearman")
}

# Neumann series: m = (I - A)^{-1} ε ≈ ε + Aε + A²ε + ...
# Returns m (total embodied emissions per unit output) and convergence info.
neumann_series <- function(A, epsilon) {
  m    <- as.numeric(epsilon)
  term <- as.numeric(epsilon)
  k_final  <- 0L
  rel_final <- Inf

  for (k in seq_len(NEUMANN_MAXIT)) {
    term      <- as.numeric(A %*% term)
    m         <- m + term
    rel_final <- max(abs(term)) / (max(abs(m)) + 1e-15)
    k_final   <- k
    if (rel_final < NEUMANN_TOL) break
  }

  list(m = m, k = k_final, converged = rel_final < NEUMANN_TOL,
       rel_err = rel_final)
}

# Statistics for one (draw, year): sector-level summaries at NACE 2d and 5d
compute_stats_bt <- function(firms_bt) {
  # firms_bt: (vat, nace2d, nace5d, scope1, upstream, euets)

  do_stats <- function(df) {
    data.frame(
      n_firms      = nrow(df),
      s1_gini      = gini(df$scope1),
      s1_p90p10    = pct_ratio(df$scope1, 0.9, 0.1),
      s1_p75p25    = pct_ratio(df$scope1, 0.75, 0.25),
      s1_var_log   = var(log(df$scope1[df$scope1 > 0] + 1)),
      up_gini      = gini(df$upstream),
      up_p90p10    = pct_ratio(df$upstream, 0.9, 0.1),
      up_p75p25    = pct_ratio(df$upstream, 0.75, 0.25),
      up_var_log   = var(log(df$upstream[df$upstream > 0] + 1)),
      gini_ratio   = gini(df$upstream) / gini(df$scope1),
      spearman_all = safe_spearman(df$scope1, df$upstream),
      spearman_ets = safe_spearman(df$scope1[df$euets == 1L],
                                    df$upstream[df$euets == 1L])
    )
  }

  stats2d <- firms_bt %>%
    filter(scope1 > 0 | upstream > 0) %>%
    group_by(nace2d) %>%
    filter(n() >= MIN_N_STATS) %>%
    group_modify(~do_stats(.x)) %>%
    ungroup()

  stats5d <- firms_bt %>%
    filter(scope1 > 0 | upstream > 0) %>%
    group_by(nace5d) %>%
    filter(n() >= MIN_N_STATS) %>%
    group_modify(~do_stats(.x)) %>%
    ungroup()

  list(stats2d = stats2d, stats5d = stats5d)
}


# =============================================================================
# SECTION 4: Main loop — year outer, draw inner
# =============================================================================
cat("── Main loop: years × draws ─────────────────────────────────\n\n")

stats2d_all <- list()
stats5d_all <- list()
conv_rows   <- list()
t0_total    <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()
  cat(sprintf("Year %d ", t))

  # ── Slice data for year t ─────────────────────────────────────────────────
  b2b_t      <- b2b[b2b$year == t, ]
  eutl_t     <- eutl[eutl$year == t, ]
  accounts_t <- accounts[accounts$year == t, ]
  imports_t  <- imports[imports$year == t, ]
  E_dep_t    <- E_deploy_panel[E_deploy_panel$year == t,
                                c("crf_group", "E_deploy")]

  # Proxy slices for year t: list of B vectors (vat → proxy value)
  proxy_t <- lapply(proxy_list, function(px) {
    sub <- px[px$year == t, c("vat", "proxy")]
    setNames(sub$proxy, sub$vat)
  })

  # ── Build A matrix for year t (shared across all B draws) ─────────────────
  all_vats <- sort(unique(c(b2b_t$vat_i_ano, b2b_t$vat_j_ano)))
  N        <- length(all_vats)
  vat_idx  <- setNames(seq_len(N), all_vats)

  # Aggregate B2B to (seller, buyer) — needed for both cost_vec and A
  b2b_agg <- b2b_t %>%
    group_by(vat_i_ano, vat_j_ano) %>%
    summarise(sales = sum(corr_sales_ij, na.rm = TRUE), .groups = "drop") %>%
    filter(sales > 0)

  # ── Cost denominator: wage_bill + domestic B2B inputs + total imports ──────
  # Under zero-profit pricing this equals revenue. Using it guarantees A row
  # sums < 1 by construction: domestic B2B is a strict subset of total cost.
  cost_vec <- rep(1e-6, N)

  # Domestic B2B input cost (row sums of purchase matrix; vat_j_ano = buyer)
  b2b_rowsums <- b2b_agg %>%
    group_by(vat_j_ano) %>%
    summarise(domestic_inputs = sum(sales), .groups = "drop")
  buyer_idx <- vat_idx[b2b_rowsums$vat_j_ano]
  ok_buyer  <- !is.na(buyer_idx)
  cost_vec[buyer_idx[ok_buyer]] <- cost_vec[buyer_idx[ok_buyer]] +
    b2b_rowsums$domestic_inputs[ok_buyer]

  # Wage bill
  acc_idx_c <- match(accounts_t$vat, all_vats)
  ok_acc    <- !is.na(acc_idx_c) & accounts_t$wage_bill > 0
  cost_vec[acc_idx_c[ok_acc]] <- cost_vec[acc_idx_c[ok_acc]] +
    accounts_t$wage_bill[ok_acc]

  # Total imports
  imp_cost_idx <- match(imports_t$vat, all_vats)
  ok_imp_cost  <- !is.na(imp_cost_idx) & !is.na(imports_t$total_imports) &
                  imports_t$total_imports > 0
  cost_vec[imp_cost_idx[ok_imp_cost]] <- cost_vec[imp_cost_idx[ok_imp_cost]] +
    imports_t$total_imports[ok_imp_cost]

  # ETS permit costs: emissions_i × carbon_price_t (regulated firms only)
  ets_cost_idx <- match(eutl_t$vat, all_vats)
  ok_ets_c     <- !is.na(ets_cost_idx) & !is.na(eutl_t$emissions) &
                  eutl_t$emissions > 0
  cost_vec[ets_cost_idx[ok_ets_c]] <- cost_vec[ets_cost_idx[ok_ets_c]] +
    eutl_t$emissions[ok_ets_c] * CARBON_PRICE[as.character(t)]

  row_i  <- vat_idx[b2b_agg$vat_j_ano]   # buyer  → row
  col_j  <- vat_idx[b2b_agg$vat_i_ano]   # seller → column
  ok_ij  <- !is.na(row_i) & !is.na(col_j)

  # A_{ij} = purchases_{ij} / cost_i; row sums < 1 guaranteed by construction
  A <- sparseMatrix(
    i    = row_i[ok_ij],
    j    = col_j[ok_ij],
    x    = b2b_agg$sales[ok_ij] / cost_vec[row_i[ok_ij]],
    dims = c(N, N)
  )

  max_rowsum <- max(rowSums(A))
  if (max_rowsum >= 1)
    cat(sprintf("\n  WARNING year %d: max row sum of A = %.4f >= 1\n", t, max_rowsum))

  # ── ETS emission intensities: emissions per unit cost (fixed across draws) ──
  eps_ets <- rep(0, N)
  ets_idx <- match(eutl_t$vat, all_vats)
  ok_ets  <- !is.na(ets_idx) & !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  eps_ets[ets_idx[ok_ets]] <- eutl_t$emissions[ok_ets] / cost_vec[ets_idx[ok_ets]]
  ets_vats_t <- eutl_t$vat[ok_ets & !is.na(ets_idx)]

  # ── Draw loop (parallelised over B draws via PSOCK + foreach) ─────────────
  # A, cost_vec, eps_ets, etc. are exported to each worker once per year.
  # Each draw gets its own ε vector and Neumann computation.
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr); library(Matrix) })
  clusterExport(cl, c("A", "cost_vec", "eps_ets", "all_vats", "ets_vats_t",
                       "proxy_t", "deploy_nace", "E_dep_t", "accounts_t",
                       "neumann_series", "compute_stats_bt",
                       "gini", "pct_ratio", "safe_spearman",
                       "NEUMANN_MAXIT", "NEUMANN_TOL", "MIN_N_STATS", "t"),
                 envir = environment())

  draw_results <- foreach(b = seq_len(B),
                           .packages = c("dplyr", "Matrix")) %dopar% {

    proxy_b <- proxy_t[[b]]   # named vector vat → proxy_avg

    deploy_b <- data.frame(
      vat       = names(proxy_b),
      proxy_avg = as.numeric(proxy_b),
      stringsAsFactors = FALSE
    ) %>%
      filter(proxy_avg > 0) %>%
      left_join(deploy_nace, by = "vat") %>%
      filter(!is.na(crf_group)) %>%
      left_join(E_dep_t, by = "crf_group") %>%
      mutate(E_deploy = coalesce(E_deploy, 0)) %>%
      group_by(crf_group) %>%
      mutate(
        sinh_sum    = sum(sinh(proxy_avg)),
        w           = if_else(sinh_sum > 0, sinh(proxy_avg) / sinh_sum, 0),
        emissions_b = E_deploy * w
      ) %>%
      ungroup() %>%
      filter(emissions_b > 0)

    # Full emission intensity vector: ETS (fixed) + imputed deployment
    eps_b   <- eps_ets
    imp_idx <- match(deploy_b$vat, all_vats)
    ok_imp  <- !is.na(imp_idx)
    eps_b[imp_idx[ok_imp]] <- eps_b[imp_idx[ok_imp]] +
      deploy_b$emissions_b[ok_imp] / cost_vec[imp_idx[ok_imp]]

    # Neumann series: m ≈ (I - A)^{-1} ε
    ns <- neumann_series(A, eps_b)

    # upstream_i = cost_i * (m_i - eps_i), floored at 0
    # scope1_i   = cost_i * eps_i = emissions_i (by construction)
    upstream_b <- pmax(cost_vec * (ns$m - eps_b), 0)
    scope1_b   <- cost_vec * eps_b

    firms_bt <- data.frame(
      vat      = all_vats,
      scope1   = scope1_b,
      upstream = upstream_b,
      stringsAsFactors = FALSE
    ) %>%
      left_join(accounts_t %>% select(vat, nace2d, nace5d), by = "vat") %>%
      filter(!is.na(nace5d), scope1 > 0 | upstream > 0) %>%
      mutate(euets = as.integer(vat %in% ets_vats_t))

    st <- compute_stats_bt(firms_bt)

    list(
      stats2d = st$stats2d %>% mutate(year = t, draw = b),
      stats5d = st$stats5d %>% mutate(year = t, draw = b),
      k       = ns$k,
      rel_err = ns$rel_err
    )
  }

  stopCluster(cl)

  # Collect results from all draws
  k_vec   <- integer(B)
  rel_vec <- numeric(B)
  for (b in seq_len(B)) {
    res <- draw_results[[b]]
    stats2d_all[[length(stats2d_all) + 1L]] <- res$stats2d
    stats5d_all[[length(stats5d_all) + 1L]] <- res$stats5d
    k_vec[b]   <- res$k
    rel_vec[b] <- res$rel_err
  }

  conv_rows[[length(conv_rows) + 1L]] <- data.frame(
    year       = t,
    max_rowsum = max_rowsum,
    med_k      = median(k_vec),
    max_k      = max(k_vec),
    pct_conv   = mean(k_vec < NEUMANN_MAXIT) * 100,
    stringsAsFactors = FALSE
  )

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  cat(sprintf("(%s s)\n", elapsed))
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 5: Convergence summary
# =============================================================================
conv_summary <- bind_rows(conv_rows)
cat("── Neumann convergence ──────────────────────────────────────\n")
print(conv_summary)
cat("\n")


# =============================================================================
# SECTION 6: Aggregate across B draws
# =============================================================================
cat("── Aggregating across draws ─────────────────────────────────\n")

all_stats2d <- bind_rows(stats2d_all)
all_stats5d <- bind_rows(stats5d_all)

stat_cols <- c("s1_gini", "s1_p90p10", "s1_p75p25", "s1_var_log",
               "up_gini", "up_p90p10", "up_p75p25", "up_var_log",
               "gini_ratio", "spearman_all", "spearman_ets")

summarise_draws <- function(df, grp) {
  df %>%
    group_by(across(all_of(grp))) %>%
    summarise(
      across(all_of(stat_cols),
             list(mean = ~mean(.x, na.rm = TRUE),
                  sd   = ~sd(.x,   na.rm = TRUE)),
             .names = "{.col}__{.fn}"),
      n_draws = n(),
      .groups = "drop"
    )
}

stats2d_summary <- summarise_draws(all_stats2d, c("nace2d", "year"))
stats5d_summary <- summarise_draws(all_stats5d, c("nace5d", "year"))

cat("  NACE 2-digit:", nrow(stats2d_summary), "sector-years\n")
cat("  NACE 5-digit:", nrow(stats5d_summary), "sector-years\n\n")


# =============================================================================
# SECTION 7: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "b_loop_statistics.RData")

save(
  stats2d_summary,   # ← copy to local 1
  stats5d_summary,   # ← copy to local 1
  all_stats2d,       # keep on RMD (full draw-level data)
  all_stats5d,
  conv_summary,
  file = OUT_PATH
)

cat("══════════════════════════════════════════════\n")
cat("Saved:", OUT_PATH, "\n")
cat("  stats2d_summary:", nrow(stats2d_summary), "rows\n")
cat("  stats5d_summary:", nrow(stats5d_summary), "rows\n")
cat("  Total time:", total_time, "min\n")
cat("══════════════════════════════════════════════\n")
