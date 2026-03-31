###############################################################################
# analysis/b_loop_scope1_dispersion.R
#
# PURPOSE
#   RQ1: How dispersed is carbon productivity across firms within sectors?
#
#   For each of B subsample draws (from run_subsampled_en.R):
#     1. NIR calibration → distribute sector-year deployment emissions
#     2. Assign scope 1 emissions: ETS observed + deployment imputed
#     3. Compute within-sector carbon productivity dispersion statistics
#   Average statistics across B draws (point estimates) + cross-B s.d.
#
#   Scope 1 = raw emissions (ETS observed or NIR-imputed). No A matrix or
#   Neumann series needed: the eps/cost round-trip in compute_b_loop.R
#   cancels out (scope1_i = cost_i × emissions_i / cost_i = emissions_i).
#
#   Carbon productivity = revenue / scope1 (EUR per tonne CO2).
#   Benchmarks: Lyubich et al (2018), De Lyon & Dechezlepretre (2025).
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData     (B proxy draws)
#   {PROC_DATA}/firm_year_belgian_euets.RData    (ETS observed emissions)
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/deployment_panel.RData           (vat, year, nace5d)
#   {PROC_DATA}/nir_calibration_targets.RData    (E_NIR_kt by CRF group × year)
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# OUTPUT
#   {PROC_DATA}/b_loop_scope1_dispersion.RData
#     stats2d_summary : point estimates + cross-B s.d. by NACE 2-digit × year
#     stats5d_summary : same at NACE 5-digit × year
#     all_stats2d     : full B × sector-year draw-level data
#     all_stats5d     : same at NACE 5-digit
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
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# ── Parameters ────────────────────────────────────────────────────────────────
YEARS       <- 2005:2021  # 2022 excluded: wage_bill data quality issue
MIN_N_STATS <- 3L         # min firms per sector-year to compute a statistic

N_CORES_SET <- if (tolower(Sys.info()[["user"]]) == "jardang") {
  40L
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 2L)
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("  RQ1: SCOPE 1 / CARBON PRODUCTIVITY DISPERSION\n")
cat("  Years:", min(YEARS), "–", max(YEARS), "\n")
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

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)
cat("  EUTL firm-years:", nrow(eutl), "\n")

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)
cat("  Accounts firm-years:", nrow(accounts), "\n")

load(file.path(PROC_DATA, "deployment_panel.RData"))

nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(nace2d = "character")
) %>%
  select(nace2d, crf_group)

deploy_nace <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
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

compute_stats_bt <- function(firms_bt) {
  # firms_bt: (vat, nace2d, nace5d, scope1, revenue, euets)

  do_stats <- function(df) {
    # Carbon productivity: revenue / scope1 (firms with both > 0)
    ok_cp   <- df$scope1 > 0 & df$revenue > 0
    cp      <- df$revenue[ok_cp] / df$scope1[ok_cp]
    n_cp    <- length(cp)
    log_cp  <- if (n_cp >= 2L) log(cp) else NA_real_

    data.frame(
      n_firms      = nrow(df),
      n_cp         = n_cp,
      # Scope 1 emission levels
      s1_gini      = gini(df$scope1),
      s1_p90p10    = pct_ratio(df$scope1, 0.9, 0.1),
      s1_p75p25    = pct_ratio(df$scope1, 0.75, 0.25),
      s1_var_log   = var(log(df$scope1[df$scope1 > 0] + 1)),
      # Carbon productivity (revenue / scope1)
      cp_p90p10    = if (n_cp >= 2L) pct_ratio(cp, 0.9, 0.1) else NA_real_,
      cp_p75p25    = if (n_cp >= 2L) pct_ratio(cp, 0.75, 0.25) else NA_real_,
      cp_var_log   = if (n_cp >= 2L) var(log_cp) else NA_real_,
      cp_p9010_log = if (n_cp >= 2L) diff(quantile(log_cp, c(0.1, 0.9))) else NA_real_
    )
  }

  stats2d <- firms_bt %>%
    filter(scope1 > 0) %>%
    group_by(nace2d) %>%
    filter(n() >= MIN_N_STATS) %>%
    group_modify(~do_stats(.x)) %>%
    ungroup()

  stats5d <- firms_bt %>%
    filter(scope1 > 0) %>%
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
t0_total    <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()
  cat(sprintf("Year %d ", t))

  eutl_t     <- eutl[eutl$year == t, ]
  accounts_t <- accounts[accounts$year == t, ]
  E_dep_t    <- E_deploy_panel[E_deploy_panel$year == t,
                                c("crf_group", "E_deploy")]

  proxy_t <- lapply(proxy_list, function(px) {
    sub <- px[px$year == t, c("vat", "proxy")]
    setNames(sub$proxy, sub$vat)
  })

  # ETS firms: vat → emissions (fixed across draws)
  ok_ets <- !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  ets_vats_t      <- eutl_t$vat[ok_ets]
  ets_emissions_t <- setNames(eutl_t$emissions[ok_ets], eutl_t$vat[ok_ets])

  # ── Draw loop ──────────────────────────────────────────────────────────────
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr) })
  clusterExport(cl, c("ets_vats_t", "ets_emissions_t",
                       "proxy_t", "deploy_nace", "E_dep_t", "accounts_t",
                       "compute_stats_bt", "gini", "pct_ratio",
                       "MIN_N_STATS", "t"),
                envir = environment())

  draw_results <- foreach(b = seq_len(B), .packages = "dplyr") %dopar% {

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

    # Scope 1 = raw emissions
    ets_df <- data.frame(
      vat = ets_vats_t, scope1 = ets_emissions_t[ets_vats_t],
      euets = 1L, stringsAsFactors = FALSE
    )
    deploy_df <- data.frame(
      vat = deploy_b$vat, scope1 = deploy_b$emissions_b,
      euets = 0L, stringsAsFactors = FALSE
    )

    firms_bt <- bind_rows(ets_df, deploy_df) %>%
      group_by(vat) %>%
      summarise(scope1 = sum(scope1), euets = max(euets), .groups = "drop") %>%
      left_join(accounts_t %>% select(vat, nace2d, nace5d, revenue), by = "vat") %>%
      filter(!is.na(nace5d), scope1 > 0)

    st <- compute_stats_bt(firms_bt)

    list(
      stats2d = st$stats2d %>% mutate(year = t, draw = b),
      stats5d = st$stats5d %>% mutate(year = t, draw = b)
    )
  }

  stopCluster(cl)

  for (b in seq_len(B)) {
    stats2d_all[[length(stats2d_all) + 1L]] <- draw_results[[b]]$stats2d
    stats5d_all[[length(stats5d_all) + 1L]] <- draw_results[[b]]$stats5d
  }

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  cat(sprintf("(%s s)\n", elapsed))
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 5: Aggregate across B draws
# =============================================================================
cat("── Aggregating across draws ─────────────────────────────────\n")

all_stats2d <- bind_rows(stats2d_all)
all_stats5d <- bind_rows(stats5d_all)

stat_cols <- c("s1_gini", "s1_p90p10", "s1_p75p25", "s1_var_log",
               "cp_p90p10", "cp_p75p25", "cp_var_log", "cp_p9010_log")

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
# SECTION 6: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "b_loop_scope1_dispersion.RData")

save(
  stats2d_summary,
  stats5d_summary,
  all_stats2d,
  all_stats5d,
  file = OUT_PATH
)

cat("══════════════════════════════════════════════════════════════\n")
cat("Saved:", OUT_PATH, "\n")
cat("  stats2d_summary:", nrow(stats2d_summary), "rows\n")
cat("  stats5d_summary:", nrow(stats5d_summary), "rows\n")
cat("  Total time:", total_time, "min\n")
cat("══════════════════════════════════════════════════════════════\n")
