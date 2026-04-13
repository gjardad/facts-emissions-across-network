###############################################################################
# analysis/build_upstream_emissions_glo.R
#
# PURPOSE
#   Build firm-level upstream embodied emissions for each year and each of
#   B perturbation draws (from the uncertainty pipeline) plus a deterministic
#   point-estimate draw (from the GLO allocation).
#
#   For each year:
#     1. Build A matrix from B2B (fixed across draws)
#     2. Build cost_vec (fixed across draws)
#     3. Factorize (I - A) once via sparse LU
#     4. For each draw b = 0..B (b=0 = deterministic, b=1..B = perturbation):
#        a. Build emission intensity vector eps_b from that draw's scope 1
#        b. Direct solve: m = (I - A)^{-1} eps  (back-solve using LU)
#        c. scope1_i = cost_i * eps_i; upstream_i = cost_i * (m_i - eps_i)
#        d. Decompose upstream into across-sector and within-sector components
#           at both NACE 2-digit and 5-digit granularity
#     5. Save firm-level results per year
#
# INPUT
#   {PROC_DATA}/uncertainty_draws_<scheme>/draw_BBBB.RData (B perturbation draws)
#   {PROC_DATA}/allocation_glo_<scheme>/alloc_YYYY.RData   (deterministic alloc)
#   {PROC_DATA}/b2b_selected_sample.RData
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/firm_year_total_imports.RData
#
# OUTPUT
#   {PROC_DATA}/upstream_emissions_glo_<scheme>/firms_YYYY.RData (one per year)
#     firms_by_draw : list of (B+1) data.frames (index 1 = deterministic,
#                     indices 2..(B+1) = perturbation draws)
#                     Columns: vat, scope1, upstream, upstream_across_2d,
#                     upstream_within_2d, upstream_across_5d,
#                     upstream_within_5d, upstream_across_crf,
#                     upstream_within_crf, source
#     max_rowsum    : max row sum of A for this year
#
# RUNS ON: RMD (40 cores) or local 1 (downsampled, fewer cores)
###############################################################################

# -- Paths --------------------------------------------------------------------
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
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)
library(Matrix)
library(parallel)
library(doParallel)
library(foreach)

# -- Parameters ---------------------------------------------------------------
WEIGHT_SCHEME <- "balanced"
YEARS         <- 2005:2007   # TEST: set to 2005:2021 for full run
SOLVER <- "direct"   # direct sparse LU solve of (I - A) x = eps

# EU ETS annual average carbon price (EUR / tonne CO2), 2005-2021.
CARBON_PRICE <- c(
  "2005" = 25.29, "2006" = 21.53, "2007" =  0.86, "2008" = 25.74,
  "2009" = 18.41, "2010" = 18.98, "2011" = 18.08, "2012" =  9.49,
  "2013" =  5.94, "2014" =  7.89, "2015" =  8.52, "2016" =  5.92,
  "2017" =  6.63, "2018" = 18.55, "2019" = 27.84, "2020" = 27.94,
  "2021" = 62.25
)

N_CORES_SET <- if (tolower(Sys.info()[["user"]]) == "jardang") {
  40L
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 2L)
}

DRAW_DIR  <- file.path(PROC_DATA, sprintf("uncertainty_draws_%s", WEIGHT_SCHEME))
ALLOC_DIR <- file.path(PROC_DATA, sprintf("allocation_glo_%s", WEIGHT_SCHEME))
OUT_DIR   <- file.path(PROC_DATA, sprintf("upstream_emissions_glo_%s", WEIGHT_SCHEME))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("===================================================================\n")
cat("  BUILD UPSTREAM EMISSIONS (GLO + uncertainty draws)\n")
cat("  WEIGHT_SCHEME =", WEIGHT_SCHEME, "\n")
cat("  Years:", min(YEARS), "--", max(YEARS), "| Solver:", SOLVER, "\n")
cat("  Cores:", N_CORES_SET, "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load shared data (B2B, accounts, imports)
# =============================================================================
cat("-- Loading shared data ----------------------------------------------\n")

load(file.path(PROC_DATA, "b2b_selected_sample.RData"))
b2b <- df_b2b_selected_sample %>% filter(year %in% YEARS)
rm(df_b2b_selected_sample)
cat("  B2B rows:", nrow(b2b), "\n")

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
  mutate(nace2d   = make_nace2d(nace5d),
         wage_bill = pmax(coalesce(wage_bill, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

# CRF crosswalk: nace2d -> crf_group
nace_crf <- read.csv(file.path(REPO_DIR, "preprocess", "crosswalks",
                                "nace_crf_crosswalk.csv"),
                      colClasses = c(nace2d = "character"))[, c("nace2d", "crf_group")]
nace_crf <- rbind(nace_crf, data.frame(nace2d = "17/18", crf_group = "paper"))
accounts <- merge(accounts, nace_crf, by = "nace2d", all.x = TRUE)
cat("  Accounts firm-years:", nrow(accounts), "\n")

load(file.path(PROC_DATA, "firm_year_total_imports.RData"))
if ("vat_ano" %in% names(firm_year_total_imports))
  firm_year_total_imports <- rename(firm_year_total_imports, vat = vat_ano)
imports <- firm_year_total_imports %>%
  filter(year %in% YEARS) %>%
  select(vat, year, total_imports)
rm(firm_year_total_imports)
cat("  Import firm-years:", nrow(imports), "\n\n")


# =============================================================================
# SECTION 2: Preprocess — reshape draws into per-year lists (option c)
# =============================================================================
cat("-- Preprocessing: reshape draws into per-year lists -----------------\n")

draw_files <- sort(list.files(DRAW_DIR, pattern = "^draw_\\d+\\.RData$",
                               full.names = TRUE))
B <- length(draw_files)
cat("  Found", B, "draw files\n")

# Initialize per-year structure
draws_by_year <- setNames(
  lapply(YEARS, function(t) vector("list", B)),
  as.character(YEARS)
)

t0_pre <- Sys.time()

for (b in seq_len(B)) {
  load(draw_files[b])  # loads draw_firms
  for (t in YEARS) {
    yr_slice <- draw_firms[draw_firms$year == t,
                            c("vat", "scope1_b", "source")]
    names(yr_slice)[names(yr_slice) == "scope1_b"] <- "scope1"
    draws_by_year[[as.character(t)]][[b]] <- yr_slice
  }
  if (b %% 25 == 0)
    cat(sprintf("  Loaded draw %d / %d\n", b, B))
}
rm(draw_firms)

# Load deterministic allocation as "draw 0" for each year
for (t in YEARS) {
  alloc_path <- file.path(ALLOC_DIR, sprintf("alloc_%d.RData", t))
  if (file.exists(alloc_path)) {
    load(alloc_path)  # loads year_firms
    det_slice <- year_firms[, c("vat", "scope1", "source")]
    draws_by_year[[as.character(t)]] <- c(list(det_slice),
                                           draws_by_year[[as.character(t)]])
    rm(year_firms)
  }
}

elapsed_pre <- round(as.numeric(difftime(Sys.time(), t0_pre, units = "secs")), 1)
cat(sprintf("  Preprocessing done in %.1fs\n", elapsed_pre))
cat(sprintf("  Draws per year: %d (1 deterministic + %d perturbation)\n\n",
            length(draws_by_year[[1]]), B))


# =============================================================================
# SECTION 3: Main loop — year outer, draw inner (parallelized)
#   Factorize (I - A) once per year via sparse LU, then back-solve per draw.
# =============================================================================
cat("-- Main loop: years x draws ----------------------------------------\n\n")

N_DRAWS_TOTAL <- B + 1L  # 1 deterministic + B perturbation
N_CORES <- min(N_DRAWS_TOTAL, N_CORES_SET)

diag_all <- list()
t0_total <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()
  cat(sprintf("Year %d ", t))

  year_draws <- draws_by_year[[as.character(t)]]
  if (is.null(year_draws) || length(year_draws) == 0) {
    cat("-- SKIPPED (no draws)\n")
    next
  }

  # -- Slice data for year t --------------------------------------------------
  b2b_t      <- b2b[b2b$year == t, ]
  eutl_t     <- eutl[eutl$year == t, ]
  accounts_t <- accounts[accounts$year == t, ]
  imports_t  <- imports[imports$year == t, ]

  # -- Build A matrix for year t (shared across all draws) --------------------
  all_vats <- sort(unique(c(b2b_t$vat_i_ano, b2b_t$vat_j_ano)))
  N        <- length(all_vats)
  vat_idx  <- setNames(seq_len(N), all_vats)

  b2b_agg <- b2b_t %>%
    group_by(vat_i_ano, vat_j_ano) %>%
    summarise(sales = sum(corr_sales_ij, na.rm = TRUE), .groups = "drop") %>%
    filter(sales > 0)

  # -- Cost denominator: wage_bill + domestic B2B inputs + total imports -------
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

  # -- Factorize (I - A) once for this year -----------------------------------
  IminusA <- Diagonal(N) - A
  LU      <- lu(IminusA)        # sparse LU factorization (reused across draws)

  # -- Sector vectors for upstream decomposition ------------------------------
  acc_match  <- match(all_vats, accounts_t$vat)
  nace2d_vec <- accounts_t$nace2d[acc_match]
  nace5d_vec <- accounts_t$nace5d[acc_match]
  crf_vec    <- accounts_t$crf_group[acc_match]

  # -- Draw loop (parallelized) -----------------------------------------------
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr); library(Matrix) })
  clusterExport(cl, c("A", "LU", "cost_vec", "all_vats",
                       "nace2d_vec", "nace5d_vec", "crf_vec",
                       "year_draws"),
                envir = environment())

  n_draws <- length(year_draws)  # B + 1

  draw_results <- foreach(b = seq_len(n_draws),
                           .packages = c("dplyr", "Matrix")) %dopar% {

    alloc_b <- year_draws[[b]]
    eps_b   <- rep(0, length(all_vats))
    alloc_idx <- match(alloc_b$vat, all_vats)
    ok_alloc  <- !is.na(alloc_idx)
    eps_b[alloc_idx[ok_alloc]] <- alloc_b$scope1[ok_alloc] /
                                   cost_vec[alloc_idx[ok_alloc]]

    # Source vector
    source_vec <- rep("none", length(all_vats))
    source_vec[alloc_idx[ok_alloc]] <- alloc_b$source[ok_alloc]

    # Direct solve: m = (I - A)^{-1} eps
    m_vec <- as.numeric(solve(LU, eps_b))

    upstream_b <- pmax(cost_vec * (m_vec - eps_b), 0)
    scope1_b   <- cost_vec * eps_b

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

    # CRF group decomposition
    m_bar_crf    <- m_vec
    has_crf      <- !is.na(crf_vec)
    if (any(has_crf)) {
      sm_crf <- tapply(m_vec[has_crf], crf_vec[has_crf], mean)
      m_bar_crf[has_crf] <- sm_crf[crf_vec[has_crf]]
    }
    upstream_across_crf_b <- pmax(cost_vec * as.numeric(A %*% m_bar_crf), 0)
    upstream_within_crf_b <- upstream_b - upstream_across_crf_b

    firms_b <- data.frame(
      vat                = all_vats,
      scope1             = scope1_b,
      upstream           = upstream_b,
      upstream_across_2d = upstream_across_2d_b,
      upstream_within_2d = upstream_within_2d_b,
      upstream_across_5d  = upstream_across_5d_b,
      upstream_within_5d  = upstream_within_5d_b,
      upstream_across_crf = upstream_across_crf_b,
      upstream_within_crf = upstream_within_crf_b,
      source              = source_vec,
      stringsAsFactors   = FALSE
    )
    firms_b <- firms_b[firms_b$scope1 > 0 | firms_b$upstream > 0, ]

    list(firms = firms_b)
  }

  stopCluster(cl)

  # Collect results
  firms_by_draw <- lapply(draw_results, `[[`, "firms")

  out_path <- file.path(OUT_DIR, sprintf("firms_%d.RData", t))
  save(firms_by_draw, max_rowsum, file = out_path)

  diag_all[[length(diag_all) + 1L]] <- data.frame(
    year = t, max_rowsum = max_rowsum, N = N, n_draws = n_draws
  )

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  n_firms_avg <- round(mean(sapply(firms_by_draw, nrow)))
  cat(sprintf("(%s s | ~%d firms/draw | %d draws)\n", elapsed, n_firms_avg, n_draws))

  rm(year_draws, LU, IminusA)
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 4: Summary
# =============================================================================
diag_summary <- bind_rows(diag_all)

cat("-- Diagnostics ---------------------------------------------------------\n")
print(diag_summary)

save(diag_summary, file = file.path(OUT_DIR, "diag_summary.RData"))

cat(sprintf("\n===================================================================\n"))
cat("Saved:", OUT_DIR, "\n")
cat(sprintf("  %d year files (firms_YYYY.RData)\n", length(YEARS)))
cat("  diag_summary.RData\n")
cat(sprintf("  Total time: %.1f min\n", total_time))
cat("===================================================================\n")
