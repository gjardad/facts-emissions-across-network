###############################################################################
# analysis/build_upstream_emissions.R
#
# PURPOSE
#   Build firm-level scope 1 and upstream embodied emissions for each year
#   and each of B = 50 subsample draws. This is the expensive step (A matrix
#   construction + Neumann series). Downstream RQ scripts load the output
#   and compute statistics cheaply.
#
#   For each year:
#     1. Build A matrix from B2B (fixed across draws)
#     2. Build cost_vec (fixed across draws)
#     3. For each draw b = 1..B:
#        a. NIR calibration → deployment ε
#        b. Neumann series: m = (I - A)^{-1} ε
#        c. scope1_i = cost_i × ε_i;  upstream_i = cost_i × (m_i - ε_i)
#        d. Decompose upstream into across-sector and within-sector components
#           at both NACE 2-digit and 5-digit granularity:
#             upstream_across_Xd_i = cost_i × Σ_j A_{ij} × m̄_{s_Xd(j)}
#             upstream_within_Xd_i = upstream_i - upstream_across_Xd_i
#           where m̄_s is the unweighted sector-average Neumann multiplier.
#     4. Save firm-level (vat, scope1, upstream, upstream_across_2d,
#        upstream_within_2d, upstream_across_5d, upstream_within_5d,
#        euets) per draw
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData
#   {PROC_DATA}/b2b_selected_sample.RData
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/deployment_panel.RData
#   {PROC_DATA}/nir_calibration_targets.RData
#   {PROC_DATA}/firm_year_total_imports.RData
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# OUTPUT
#   {PROC_DATA}/upstream_emissions/firms_YYYY.RData  (one file per year)
#     firms_by_draw : list of B data.frames, each with columns
#                     (vat, scope1, upstream, upstream_across_2d,
#                      upstream_within_2d, upstream_across_5d,
#                      upstream_within_5d, euets)
#                     filtered to scope1 > 0 | upstream > 0
#     conv_info     : data.frame with (draw, k, rel_err, converged)
#     max_rowsum    : scalar, max row sum of A for this year
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
YEARS         <- 2005:2021  # 2022 excluded: wage_bill data quality issue
NEUMANN_MAXIT <- 50L
NEUMANN_TOL   <- 1e-8
MIN_N_STATS   <- 3L

# EU ETS annual average carbon price (EUR / tonne CO2), 2005–2022.
# Source: ICAP Allowance Price Explorer, secondary market daily prices.
CARBON_PRICE <- c(
  "2005" = 25.29, "2006" = 21.53, "2007" =  0.86, "2008" = 25.74,
  "2009" = 18.41, "2010" = 18.98, "2011" = 18.08, "2012" =  9.49,
  "2013" =  5.94, "2014" =  7.89, "2015" =  8.52, "2016" =  5.92,
  "2017" =  6.63, "2018" = 18.55, "2019" = 27.84, "2020" = 27.94,
  "2021" = 62.25, "2022" = 85.51
)

N_CORES_SET <- if (tolower(Sys.info()[["user"]]) == "jardang") {
  40L
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 2L)
}

OUT_DIR <- file.path(PROC_DATA, "upstream_emissions")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  BUILD UPSTREAM EMISSIONS: A matrix + Neumann per year × draw\n")
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

load(file.path(PROC_DATA, "nir_calibration_targets.RData"))

E_ETS_group <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group), year %in% YEARS) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions, na.rm = TRUE), .groups = "drop")

E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
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


# =============================================================================
# SECTION 4: Main loop — year outer, draw inner
# =============================================================================
cat("── Main loop: years × draws ─────────────────────────────────\n\n")

conv_all <- list()
t0_total <- Sys.time()

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

  proxy_t <- lapply(proxy_list, function(px) {
    sub <- px[px$year == t, c("vat", "proxy")]
    setNames(sub$proxy, sub$vat)
  })

  # ── Build A matrix for year t (shared across all B draws) ─────────────────
  all_vats <- sort(unique(c(b2b_t$vat_i_ano, b2b_t$vat_j_ano)))
  N        <- length(all_vats)
  vat_idx  <- setNames(seq_len(N), all_vats)

  b2b_agg <- b2b_t %>%
    group_by(vat_i_ano, vat_j_ano) %>%
    summarise(sales = sum(corr_sales_ij, na.rm = TRUE), .groups = "drop") %>%
    filter(sales > 0)

  # ── Cost denominator: wage_bill + domestic B2B inputs + total imports ──────
  cost_vec <- rep(1e-6, N)

  b2b_rowsums <- b2b_agg %>%
    group_by(vat_j_ano) %>%
    summarise(domestic_inputs = sum(sales), .groups = "drop")
  buyer_idx <- vat_idx[b2b_rowsums$vat_j_ano]
  ok_buyer  <- !is.na(buyer_idx)
  cost_vec[buyer_idx[ok_buyer]] <- cost_vec[buyer_idx[ok_buyer]] +
    b2b_rowsums$domestic_inputs[ok_buyer]

  acc_idx_c <- match(accounts_t$vat, all_vats)
  ok_acc    <- !is.na(acc_idx_c) & accounts_t$wage_bill > 0
  cost_vec[acc_idx_c[ok_acc]] <- cost_vec[acc_idx_c[ok_acc]] +
    accounts_t$wage_bill[ok_acc]

  imp_cost_idx <- match(imports_t$vat, all_vats)
  ok_imp_cost  <- !is.na(imp_cost_idx) & !is.na(imports_t$total_imports) &
                  imports_t$total_imports > 0
  cost_vec[imp_cost_idx[ok_imp_cost]] <- cost_vec[imp_cost_idx[ok_imp_cost]] +
    imports_t$total_imports[ok_imp_cost]

  ets_cost_idx <- match(eutl_t$vat, all_vats)
  ok_ets_c     <- !is.na(ets_cost_idx) & !is.na(eutl_t$emissions) &
                  eutl_t$emissions > 0
  cost_vec[ets_cost_idx[ok_ets_c]] <- cost_vec[ets_cost_idx[ok_ets_c]] +
    eutl_t$emissions[ok_ets_c] * CARBON_PRICE[as.character(t)]

  row_i  <- vat_idx[b2b_agg$vat_j_ano]
  col_j  <- vat_idx[b2b_agg$vat_i_ano]
  ok_ij  <- !is.na(row_i) & !is.na(col_j)

  A <- sparseMatrix(
    i    = row_i[ok_ij],
    j    = col_j[ok_ij],
    x    = b2b_agg$sales[ok_ij] / cost_vec[row_i[ok_ij]],
    dims = c(N, N)
  )

  max_rowsum <- max(rowSums(A))
  if (max_rowsum >= 1)
    cat(sprintf("\n  WARNING year %d: max row sum of A = %.4f >= 1\n", t, max_rowsum))

  # ── NACE sectors for each firm in all_vats (for upstream decomposition) ────
  acc_match  <- match(all_vats, accounts_t$vat)
  nace2d_vec <- accounts_t$nace2d[acc_match]
  nace5d_vec <- accounts_t$nace5d[acc_match]
  # nace{2,5}d_vec[i] = NA if firm i is not in accounts

  # ── ETS emission intensities (fixed across draws) ─────────────────────────
  eps_ets <- rep(0, N)
  ets_idx <- match(eutl_t$vat, all_vats)
  ok_ets  <- !is.na(ets_idx) & !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  eps_ets[ets_idx[ok_ets]] <- eutl_t$emissions[ok_ets] / cost_vec[ets_idx[ok_ets]]
  ets_vats_t <- eutl_t$vat[ok_ets & !is.na(ets_idx)]

  # ── Draw loop (parallelised) ──────────────────────────────────────────────
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr); library(Matrix) })
  clusterExport(cl, c("A", "cost_vec", "eps_ets", "all_vats", "ets_vats_t",
                       "nace2d_vec", "nace5d_vec",
                       "proxy_t", "deploy_nace", "E_dep_t",
                       "neumann_series",
                       "NEUMANN_MAXIT", "NEUMANN_TOL"),
                envir = environment())

  draw_results <- foreach(b = seq_len(B),
                           .packages = c("dplyr", "Matrix")) %dopar% {

    proxy_b <- proxy_t[[b]]

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

    # scope1 and upstream in levels (tonnes CO2)
    upstream_b <- pmax(cost_vec * (ns$m - eps_b), 0)
    scope1_b   <- cost_vec * eps_b

    # ── Decompose upstream: across-sector vs within-sector ──────────────────
    # For each granularity (NACE 2d, 5d):
    #   m_bar[j] = mean(m[firms in same sector as j])
    #   upstream_across = cost × A %*% m_bar
    #   upstream_within = upstream - upstream_across
    # For firms without a known sector, m_bar[j] = m[j] (no decomposition).
    m_vec <- ns$m

    # NACE 2-digit decomposition
    m_bar_2d    <- m_vec
    has_nace_2d <- !is.na(nace2d_vec)
    if (any(has_nace_2d)) {
      sm_2d <- tapply(m_vec[has_nace_2d], nace2d_vec[has_nace_2d], mean)
      m_bar_2d[has_nace_2d] <- sm_2d[nace2d_vec[has_nace_2d]]
    }
    upstream_across_2d_b <- pmax(cost_vec * as.numeric(A %*% m_bar_2d), 0)
    upstream_within_2d_b <- upstream_b - upstream_across_2d_b

    # NACE 5-digit decomposition
    m_bar_5d    <- m_vec
    has_nace_5d <- !is.na(nace5d_vec)
    if (any(has_nace_5d)) {
      sm_5d <- tapply(m_vec[has_nace_5d], nace5d_vec[has_nace_5d], mean)
      m_bar_5d[has_nace_5d] <- sm_5d[nace5d_vec[has_nace_5d]]
    }
    upstream_across_5d_b <- pmax(cost_vec * as.numeric(A %*% m_bar_5d), 0)
    upstream_within_5d_b <- upstream_b - upstream_across_5d_b

    # Keep only firms with positive emissions
    firms_b <- data.frame(
      vat                = all_vats,
      scope1             = scope1_b,
      upstream           = upstream_b,
      upstream_across_2d = upstream_across_2d_b,
      upstream_within_2d = upstream_within_2d_b,
      upstream_across_5d = upstream_across_5d_b,
      upstream_within_5d = upstream_within_5d_b,
      stringsAsFactors   = FALSE
    )
    firms_b$euets <- as.integer(firms_b$vat %in% ets_vats_t)
    firms_b <- firms_b[firms_b$scope1 > 0 | firms_b$upstream > 0, ]

    list(
      firms   = firms_b,
      k       = ns$k,
      rel_err = ns$rel_err,
      converged = ns$converged
    )
  }

  stopCluster(cl)

  # Collect results
  firms_by_draw <- vector("list", B)
  conv_info     <- data.frame(draw = integer(B), k = integer(B),
                              rel_err = numeric(B), converged = logical(B))
  for (b in seq_len(B)) {
    firms_by_draw[[b]] <- draw_results[[b]]$firms
    conv_info$draw[b]      <- b
    conv_info$k[b]         <- draw_results[[b]]$k
    conv_info$rel_err[b]   <- draw_results[[b]]$rel_err
    conv_info$converged[b] <- draw_results[[b]]$converged
  }

  # Save year file
  out_path <- file.path(OUT_DIR, sprintf("firms_%d.RData", t))
  save(firms_by_draw, conv_info, max_rowsum, file = out_path)

  conv_all[[length(conv_all) + 1L]] <- cbind(year = t, max_rowsum = max_rowsum,
                                             conv_info)

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  n_firms_avg <- round(mean(sapply(firms_by_draw, nrow)))
  cat(sprintf("(%s s | ~%d firms/draw)\n", elapsed, n_firms_avg))
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 5: Convergence summary
# =============================================================================
conv_summary <- bind_rows(conv_all) %>%
  group_by(year) %>%
  summarise(
    max_rowsum = first(max_rowsum),
    med_k      = median(k),
    max_k      = max(k),
    pct_conv   = mean(converged) * 100,
    .groups    = "drop"
  )

cat("── Neumann convergence ──────────────────────────────────────\n")
print(conv_summary)

# Save convergence summary alongside year files
save(conv_summary, file = file.path(OUT_DIR, "conv_summary.RData"))

cat(sprintf("\n══════════════════════════════════════════════════════════════\n"))
cat("Saved:", OUT_DIR, "\n")
cat("  ", length(YEARS), "year files (firms_YYYY.RData)\n")
cat("  conv_summary.RData\n")
cat("  Total time:", total_time, "min\n")
cat("══════════════════════════════════════════════════════════════\n")
