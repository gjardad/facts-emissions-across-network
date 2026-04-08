###############################################################################
# analysis/b_loop_scope1_dispersion_pi.R
#
# PURPOSE
#   Compute scope 1 dispersion measures with prediction intervals that account
#   for two sources of uncertainty:
#     Source 1: estimation uncertainty (B subsample draws of EN coefficients)
#     Source 2: prediction error (residual noise in proxy, calibrated from CV)
#
#   For each (year, draw b, replication r), perturbs deployment proxies by
#   empirical OOS residuals, re-runs Pareto allocation, and computes dispersion.
#   Produces B x R realizations per sector-year, enabling variance decomposition.
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData   (B proxy vectors)
#   {PROC_DATA}/deployment_draw_params.RData  (GPA params + threshold per draw)
#   {PROC_DATA}/oos_residuals.RData           (empirical prediction errors)
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/deployment_panel.RData
#   {PROC_DATA}/nir_calibration_targets.RData
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# OUTPUT
#   {PROC_DATA}/b_loop_scope1_dispersion_pi.RData
#     pi_stats2d : full B x R draw-level dispersion at NACE 2-digit x year
#     pi_summary2d : point estimates + PIs + variance decomposition
#
# RUNS ON: local 1 (lightweight — no EN fitting, no B2B)
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
source(file.path(REPO_DIR, "utils", "calibration_pareto.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# -- Parameters ---------------------------------------------------------------
YEARS       <- 2005:2021
R_REPS      <- 50L          # inner replications per draw (Source 2)
MIN_N_STATS <- 3L           # min firms for dispersion stats
DEPLOY_CAP  <- 30000        # physics-based upper bound (30 kt CO2/year)
SEED_OFFSET <- 999000L      # seeds for perturbation RNG

N_CORES_SET <- if (tolower(Sys.info()[["user"]]) == "jardang") {
  40L
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 2L)
}

cat("===================================================================\n")
cat("  SCOPE 1 DISPERSION WITH PREDICTION INTERVALS\n")
cat("  Years:", min(YEARS), "--", max(YEARS), "\n")
cat("  R (inner replications):", R_REPS, "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load data
# =============================================================================
cat("-- Loading data -----------------------------------------------------\n")

# Proxy list (B draws)
load(file.path(PROC_DATA, "deployment_proxy_list.RData"))
B <- length(proxy_list)
cat("  proxy_list: B =", B, "draws\n")

# Draw-specific GPA + threshold
draw_params_path <- file.path(PROC_DATA, "deployment_draw_params.RData")
if (file.exists(draw_params_path)) {
  load(draw_params_path)
  USE_DRAW_PARAMS <- TRUE
  cat("  draw_params:", length(draw_params), "draws\n")
} else {
  USE_DRAW_PARAMS <- FALSE
  cat("  draw_params not found -- computing global GPA fallback\n")
  FALLBACK_THRESHOLD <- 0.62
}

# OOS residuals (empirical prediction errors)
oos_path <- file.path(PROC_DATA, "oos_residuals.RData")
if (file.exists(oos_path)) {
  load(oos_path)
  epsilon_pool <- oos_residuals$epsilon
  cat("  OOS residuals:", length(epsilon_pool), "observations\n")
  cat("    mean =", round(mean(epsilon_pool), 4),
      " sd =", round(sd(epsilon_pool), 4), "\n")
} else {
  cat("  OOS residuals not found -- using parametric approximation (sd=2.5)\n")
  epsilon_pool <- rnorm(10000, mean = 0, sd = 2.5)
}

# ETS firms
load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)
cat("  EUTL firm-years:", nrow(eutl), "\n")

# Accounts (revenue + NACE)
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)
cat("  Accounts firm-years:", nrow(accounts), "\n")

# Deployment panel (for NACE + CRF mapping)
load(file.path(PROC_DATA, "deployment_panel.RData"))
nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(nace2d = "character")
) %>% select(nace2d, crf_group)

deploy_nace <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  select(vat, nace2d, nace5d) %>%
  distinct(vat, .keep_all = TRUE) %>%
  left_join(nace_crf, by = "nace2d")
rm(deployment_panel)

# NIR calibration targets
load(file.path(PROC_DATA, "nir_calibration_targets.RData"))

E_ETS_group <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group), year %in% YEARS) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions, na.rm = TRUE), .groups = "drop")

# Pre-ETS backcast (reuse from allocation script)
# Load the allocation file for year 2005 to get pre_ets structure
# For simplicity, recompute pre-ETS inline (same logic as b_allocation_pareto.R)
eutl_pos <- eutl %>% filter(!is.na(emissions), emissions > 0)
ets_entry <- eutl_pos %>%
  group_by(vat) %>%
  summarise(first_ets_year = min(year), .groups = "drop")
late_entrants <- ets_entry %>% filter(first_ets_year > min(YEARS))

anchor <- eutl_pos %>%
  inner_join(late_entrants, by = "vat") %>%
  group_by(vat) %>%
  arrange(year) %>%
  slice_head(n = 2) %>%
  summarise(e_anchor = mean(emissions), anchor_year = min(year), .groups = "drop") %>%
  left_join(
    accounts %>% group_by(vat) %>%
      summarise(nace2d = names(which.max(table(nace2d))), .groups = "drop"),
    by = "vat"
  ) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group))

sector_year_mu <- eutl_pos %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d)) %>%
  group_by(nace2d, year) %>%
  summarise(mu = mean(log(emissions)), .groups = "drop")

pre_ets_list <- list()
for (i in seq_len(nrow(anchor))) {
  v <- anchor$vat[i]; s <- anchor$nace2d[i]; T_entry <- anchor$anchor_year[i]
  e_anch <- anchor$e_anchor[i]; crf <- anchor$crf_group[i]
  mu_T <- sector_year_mu$mu[sector_year_mu$nace2d == s &
                             sector_year_mu$year == T_entry]
  if (length(mu_T) == 0) next
  pre_years <- accounts$year[accounts$vat == v & accounts$year < T_entry &
                              accounts$year %in% YEARS]
  pre_years <- unique(pre_years)
  if (length(pre_years) == 0) next
  for (yr in pre_years) {
    mu_t <- sector_year_mu$mu[sector_year_mu$nace2d == s &
                               sector_year_mu$year == yr]
    e_hat <- if (length(mu_t) > 0) e_anch * exp(mu_t - mu_T) else e_anch
    pre_ets_list[[length(pre_ets_list) + 1]] <- data.frame(
      vat = v, year = yr, nace2d = s, crf_group = crf,
      emissions = e_hat, stringsAsFactors = FALSE
    )
  }
}
pre_ets <- if (length(pre_ets_list) > 0) bind_rows(pre_ets_list) else
  data.frame(vat = character(0), year = integer(0), nace2d = character(0),
             crf_group = character(0), emissions = numeric(0))

E_pre_ets <- pre_ets %>%
  group_by(crf_group, year) %>%
  summarise(E_pre_ets = sum(emissions), n_pre_ets = n(), .groups = "drop")

# E_deploy panel
E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(E_ETS_group, by = c("crf_group", "year")) %>%
  left_join(E_pre_ets, by = c("crf_group", "year")) %>%
  mutate(E_ETS     = coalesce(E_ETS, 0),
         E_pre_ets = coalesce(E_pre_ets, 0),
         E_deploy  = pmax(E_NIR - E_ETS - E_pre_ets, 0))

cat("  E_deploy panel:", nrow(E_deploy_panel), "group-years\n")
cat("  Pre-ETS backcast:", nrow(pre_ets), "firm-years\n")

# GPA fallback if draw_params not available
if (!USE_DRAW_PARAMS) {
  eutl_with_nace <- eutl %>%
    filter(emissions > 0) %>%
    left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
    filter(!is.na(nace2d))
  fallback_ref <- build_reference_dist(eutl_with_nace$emissions,
                                        eutl_with_nace$year,
                                        eutl_with_nace$nace2d)
  fallback_gpa <- fit_gpa(fallback_ref)
  if (is.null(fallback_gpa)) stop("GPA fitting failed")
  cat("  GPA fallback: xi =", round(fallback_gpa["xi"], 4),
      " alpha =", round(fallback_gpa["alpha"], 4),
      " k =", round(fallback_gpa["k"], 4), "\n")
  draw_params <- replicate(B, list(gpa_params = fallback_gpa,
                                     threshold_pct = FALLBACK_THRESHOLD),
                            simplify = FALSE)
  rm(eutl_with_nace, fallback_ref)
}
cat("\n")


# =============================================================================
# SECTION 2: Helper functions
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

compute_stats_sector <- function(df) {
  ok_cp  <- df$scope1 > 0 & df$revenue > 0
  cp     <- df$revenue[ok_cp] / df$scope1[ok_cp]
  n_cp   <- length(cp)
  log_cp <- if (n_cp >= 2L) log(cp) else NA_real_

  data.frame(
    n_firms      = nrow(df),
    n_cp         = n_cp,
    s1_gini      = gini(df$scope1),
    s1_p90p10    = pct_ratio(df$scope1, 0.9, 0.1),
    s1_p75p25    = pct_ratio(df$scope1, 0.75, 0.25),
    s1_var_log   = var(log(df$scope1[df$scope1 > 0] + 1)),
    cp_p90p10    = if (n_cp >= 2L) pct_ratio(cp, 0.9, 0.1) else NA_real_,
    cp_p75p25    = if (n_cp >= 2L) pct_ratio(cp, 0.75, 0.25) else NA_real_,
    cp_var_log   = if (n_cp >= 2L) var(log_cp) else NA_real_,
    cp_p9010_log = if (n_cp >= 2L) diff(quantile(log_cp, c(0.1, 0.9))) else NA_real_
  )
}

# Allocate emissions for a single CRF group given perturbed proxy values.
# Returns data.frame (vat, scope1) for deployment firms, or NULL if no allocation.
allocate_crf_group <- function(deploy_vats, proxy_vals, n_observed,
                                E_deploy, gpa_params, threshold, deploy_cap) {
  if (length(proxy_vals) == 0 || E_deploy <= 0) return(NULL)

  proxy_levels_full <- pmax(sinh(proxy_vals), 0)
  thr_current       <- threshold
  MAX_ITER          <- 20L

  for (iter in seq_len(MAX_ITER)) {
    proxy_levels <- proxy_levels_full

    # Apply threshold
    pos_mask <- proxy_levels > 0
    if (sum(pos_mask) >= 2 && thr_current > 0) {
      tau <- quantile(proxy_levels[pos_mask], probs = thr_current, names = FALSE)
      proxy_levels[proxy_levels <= tau] <- 0
    }

    emit_mask     <- proxy_levels > 0
    n_deploy_emit <- sum(emit_mask)

    if (n_deploy_emit == 0) {
      if (thr_current > 1e-4) { thr_current <- thr_current / 2; next }
      return(NULL)
    }

    # Combined ranking
    n_total    <- n_observed + n_deploy_emit
    emit_idx   <- which(emit_mask)
    emit_proxy <- proxy_levels[emit_idx]
    ranks      <- rank(emit_proxy, ties.method = "average")

    # GPA weights
    pw <- pareto_weights(ranks, n_total, gpa_params)

    if (is.null(pw)) {
      total_sinh <- sum(sinh(proxy_vals[emit_idx]))
      if (total_sinh > 0) {
        emissions_out <- E_deploy * sinh(proxy_vals[emit_idx]) / total_sinh
      } else {
        emissions_out <- rep(E_deploy / n_deploy_emit, n_deploy_emit)
      }
    } else {
      emissions_out <- E_deploy * pw / sum(pw)
    }

    # Physics-based cap
    if (max(emissions_out) >= deploy_cap) {
      if (thr_current > 1e-4) { thr_current <- thr_current / 2; next }
      cap  <- deploy_cap * (1 - 1e-6)
      over <- which(emissions_out > cap)
      freed <- sum(emissions_out[over] - cap)
      emissions_out[over] <- cap
      under <- which(emissions_out <= cap)
      if (length(under) > 0) {
        emissions_out[under] <- emissions_out[under] +
          freed * emissions_out[under] / sum(emissions_out[under])
      }
    }
    break
  }

  data.frame(
    vat    = deploy_vats[emit_idx],
    scope1 = as.numeric(emissions_out),
    euets  = 0L,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# SECTION 3: Main double loop — year outer, (draw, rep) inner
# =============================================================================
cat("-- Main loop: years x draws x reps ---------------------------------\n\n")

stats2d_all <- list()
t0_total    <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()

  # -- Year-specific fixed data ------------------------------------------------
  eutl_t  <- eutl[eutl$year == t, ]
  E_dep_t <- E_deploy_panel[E_deploy_panel$year == t, c("crf_group", "E_deploy")]

  accounts_t <- accounts[accounts$year == t, ]

  # ETS firms with positive emissions
  ok_ets <- !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  ets_vats_t      <- eutl_t$vat[ok_ets]
  ets_emissions_t <- setNames(eutl_t$emissions[ok_ets], eutl_t$vat[ok_ets])

  ets_crf_t <- eutl_t[ok_ets, ] %>%
    left_join(accounts_t %>% distinct(vat, nace2d), by = "vat") %>%
    left_join(nace_crf, by = "nace2d") %>%
    filter(!is.na(crf_group)) %>%
    select(vat, emissions, crf_group)

  # Pre-ETS firms
  pre_ets_t <- pre_ets[pre_ets$year == t, ]
  pre_ets_vats_t      <- pre_ets_t$vat
  pre_ets_emissions_t <- setNames(pre_ets_t$emissions, pre_ets_t$vat)
  pre_ets_crf_t       <- pre_ets_t %>% select(vat, emissions, crf_group)

  # Proxy vectors per draw (named vectors: vat -> proxy)
  proxy_t <- lapply(proxy_list, function(px) {
    sub <- px[px$year == t, c("vat", "proxy")]
    setNames(sub$proxy, sub$vat)
  })

  # -- Parallel loop over (draw, rep) pairs ------------------------------------
  N_CORES <- min(B, N_CORES_SET)
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr); library(lmom) })

  # Build task grid: each task is (draw_b, rep_r) where r=0 is unperturbed
  tasks <- expand.grid(b = seq_len(B), r = 0:R_REPS)

  # Export year-specific data
  draw_gpa_list <- lapply(draw_params, `[[`, "gpa_params")
  draw_thr_list <- sapply(draw_params, `[[`, "threshold_pct")

  clusterExport(cl, c(
    "ets_vats_t", "ets_emissions_t", "ets_crf_t",
    "pre_ets_vats_t", "pre_ets_emissions_t", "pre_ets_crf_t",
    "proxy_t", "deploy_nace", "E_dep_t", "nace_crf",
    "draw_gpa_list", "draw_thr_list",
    "epsilon_pool", "DEPLOY_CAP", "SEED_OFFSET", "t",
    "accounts_t", "MIN_N_STATS",
    "pareto_weights", "allocate_crf_group",
    "gini", "pct_ratio", "compute_stats_sector"
  ), envir = environment())

  year_results <- foreach(
    task_idx = seq_len(nrow(tasks)),
    .packages = c("dplyr", "lmom")
  ) %dopar% {

    b <- tasks$b[task_idx]
    r <- tasks$r[task_idx]

    proxy_b     <- proxy_t[[b]]
    gpa_b       <- draw_gpa_list[[b]]
    threshold_b <- draw_thr_list[b]
    if (is.na(threshold_b)) threshold_b <- 0

    # Perturb proxy if r > 0
    if (r > 0) {
      set.seed(SEED_OFFSET + t * 1000L + b * 100L + r)
      eps <- sample(epsilon_pool, length(proxy_b), replace = TRUE)
      proxy_b <- pmax(proxy_b + eps, 0)
    }

    # Build deployment firms with CRF group
    deploy_b <- data.frame(
      vat       = names(proxy_b),
      proxy_val = as.numeric(proxy_b),
      stringsAsFactors = FALSE
    ) %>%
      filter(proxy_val > 0, !(vat %in% pre_ets_vats_t)) %>%
      left_join(deploy_nace, by = "vat") %>%
      filter(!is.na(crf_group)) %>%
      left_join(E_dep_t, by = "crf_group") %>%
      mutate(E_deploy = coalesce(E_deploy, 0))

    all_crf_groups <- union(
      unique(deploy_b$crf_group),
      union(unique(ets_crf_t$crf_group), unique(pre_ets_crf_t$crf_group))
    )

    # Allocate within each CRF group
    deploy_imputed <- list()
    for (cg in all_crf_groups) {
      cg_deploy  <- deploy_b[deploy_b$crf_group == cg, ]
      cg_ets     <- ets_crf_t[ets_crf_t$crf_group == cg, ]
      cg_pre_ets <- pre_ets_crf_t[pre_ets_crf_t$crf_group == cg, ]
      n_observed <- nrow(cg_ets) + nrow(cg_pre_ets)
      E_cg       <- if (nrow(cg_deploy) > 0) cg_deploy$E_deploy[1] else {
        ed <- E_dep_t$E_deploy[E_dep_t$crf_group == cg]
        if (length(ed) > 0) ed[1] else 0
      }

      if (nrow(cg_deploy) == 0 || E_cg <= 0) next

      res <- allocate_crf_group(
        deploy_vats = cg_deploy$vat,
        proxy_vals  = cg_deploy$proxy_val,
        n_observed  = n_observed,
        E_deploy    = E_cg,
        gpa_params  = gpa_b,
        threshold   = threshold_b,
        deploy_cap  = DEPLOY_CAP
      )
      if (!is.null(res)) deploy_imputed[[length(deploy_imputed) + 1]] <- res
    }

    # Combine ETS + pre_ets + deployment
    ets_df <- data.frame(
      vat = ets_vats_t, scope1 = ets_emissions_t[ets_vats_t],
      euets = 1L, stringsAsFactors = FALSE
    )
    pre_ets_df <- if (length(pre_ets_vats_t) > 0) {
      data.frame(vat = pre_ets_vats_t,
                 scope1 = pre_ets_emissions_t[pre_ets_vats_t],
                 euets = 0L, stringsAsFactors = FALSE)
    } else {
      data.frame(vat = character(0), scope1 = numeric(0), euets = integer(0))
    }
    deploy_df <- if (length(deploy_imputed) > 0) bind_rows(deploy_imputed) else
      data.frame(vat = character(0), scope1 = numeric(0), euets = integer(0))

    firms_btr <- bind_rows(ets_df, pre_ets_df, deploy_df) %>%
      group_by(vat) %>%
      summarise(scope1 = sum(scope1), euets = max(euets), .groups = "drop") %>%
      filter(scope1 > 0) %>%
      left_join(accounts_t %>% select(vat, nace2d, nace5d, revenue), by = "vat") %>%
      filter(!is.na(nace5d))

    # Compute dispersion stats (NACE 2-digit only for PI)
    stats2d <- firms_btr %>%
      filter(scope1 > 0) %>%
      group_by(nace2d) %>%
      filter(n() >= MIN_N_STATS) %>%
      group_modify(~compute_stats_sector(.x)) %>%
      ungroup() %>%
      mutate(year = t, draw = b, rep = r)

    stats2d
  }

  stopCluster(cl)

  # Collect results for this year
  stats2d_year <- bind_rows(year_results)
  stats2d_all[[length(stats2d_all) + 1L]] <- stats2d_year

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  n_rows  <- nrow(stats2d_year)
  cat(sprintf("Year %d (%s s | %d stat rows)\n", t, elapsed, n_rows))

  rm(year_results, stats2d_year, proxy_t)
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 4: Aggregate — variance decomposition
# =============================================================================
cat("-- Aggregating: variance decomposition ------------------------------\n")

pi_stats2d <- bind_rows(stats2d_all)
rm(stats2d_all)

stat_cols <- c("s1_gini", "s1_p90p10", "s1_p75p25", "s1_var_log",
               "cp_p90p10", "cp_p75p25", "cp_var_log", "cp_p9010_log")

# For each (nace2d, year), compute:
#   - Point estimate: mean across all B x R realizations
#   - 90% PI: [p5, p95]
#   - Variance decomposition: Source 1 (across-B) + Source 2 (within-B)
pi_summary2d <- pi_stats2d %>%
  group_by(nace2d, year) %>%
  summarise(
    across(all_of(stat_cols), list(
      # Point estimate (grand mean)
      mean = ~mean(.x, na.rm = TRUE),
      # Total SD
      sd   = ~sd(.x, na.rm = TRUE),
      # 90% PI bounds
      p05  = ~quantile(.x, 0.05, na.rm = TRUE),
      p95  = ~quantile(.x, 0.95, na.rm = TRUE),
      # 50% PI bounds
      p25  = ~quantile(.x, 0.25, na.rm = TRUE),
      p75  = ~quantile(.x, 0.75, na.rm = TRUE)
    ), .names = "{.col}__{.fn}"),
    n_total = n(),
    .groups = "drop"
  )

# Variance decomposition: need draw-level aggregation first
# Source 1 = Var_b[E_r[stat | b]]  (between-draw variance of draw means)
# Source 2 = E_b[Var_r[stat | b]]  (average within-draw variance)
draw_means <- pi_stats2d %>%
  group_by(nace2d, year, draw) %>%
  summarise(
    across(all_of(stat_cols), ~mean(.x, na.rm = TRUE), .names = "{.col}__draw_mean"),
    .groups = "drop"
  )

draw_vars <- pi_stats2d %>%
  group_by(nace2d, year, draw) %>%
  summarise(
    across(all_of(stat_cols), ~var(.x, na.rm = TRUE), .names = "{.col}__draw_var"),
    .groups = "drop"
  )

# Source 1: variance of draw means
source1 <- draw_means %>%
  group_by(nace2d, year) %>%
  summarise(
    across(ends_with("__draw_mean"),
           ~var(.x, na.rm = TRUE),
           .names = "{.col}"),
    .groups = "drop"
  )
names(source1) <- sub("__draw_mean$", "__var_s1", names(source1))

# Source 2: mean of within-draw variances
source2 <- draw_vars %>%
  group_by(nace2d, year) %>%
  summarise(
    across(ends_with("__draw_var"),
           ~mean(.x, na.rm = TRUE),
           .names = "{.col}"),
    .groups = "drop"
  )
names(source2) <- sub("__draw_var$", "__var_s2", names(source2))

# Merge decomposition into summary
pi_summary2d <- pi_summary2d %>%
  left_join(source1, by = c("nace2d", "year")) %>%
  left_join(source2, by = c("nace2d", "year"))

# Add percentage shares
for (sc in stat_cols) {
  v1 <- pi_summary2d[[paste0(sc, "__var_s1")]]
  v2 <- pi_summary2d[[paste0(sc, "__var_s2")]]
  vtot <- v1 + v2
  pi_summary2d[[paste0(sc, "__pct_s1")]] <- ifelse(vtot > 0, 100 * v1 / vtot, NA_real_)
  pi_summary2d[[paste0(sc, "__pct_s2")]] <- ifelse(vtot > 0, 100 * v2 / vtot, NA_real_)
}

cat("  pi_stats2d:", nrow(pi_stats2d), "rows (full B x R draw-level)\n")
cat("  pi_summary2d:", nrow(pi_summary2d), "sector-years\n\n")

# Quick summary of variance shares
cat("-- Variance decomposition summary (median across sector-years) -----\n")
for (sc in c("s1_gini", "s1_p90p10", "cp_var_log")) {
  pct1 <- median(pi_summary2d[[paste0(sc, "__pct_s1")]], na.rm = TRUE)
  pct2 <- median(pi_summary2d[[paste0(sc, "__pct_s2")]], na.rm = TRUE)
  cat(sprintf("  %-15s: Source 1 (estimation) = %5.1f%%  Source 2 (prediction) = %5.1f%%\n",
              sc, pct1, pct2))
}


# =============================================================================
# SECTION 5: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "b_loop_scope1_dispersion_pi.RData")

save(
  pi_stats2d,
  pi_summary2d,
  file = OUT_PATH
)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat("  pi_stats2d:   ", nrow(pi_stats2d), "rows\n")
cat("  pi_summary2d: ", nrow(pi_summary2d), "rows\n")
cat("  Total time:   ", total_time, "min\n")
cat("  B =", B, "| R =", R_REPS, "| B x R =", B * (R_REPS + 1), "\n")
cat("===================================================================\n")
