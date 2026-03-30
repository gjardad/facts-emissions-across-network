###############################################################################
# analysis/b_loop_scope1_dispersion_pareto.R
#
# PURPOSE
#   RQ1: Carbon productivity dispersion using Pareto (GPA) redistribution.
#
#   Same as b_loop_scope1_dispersion.R but replaces sinh-weighted allocation
#   with rank-based GPA quantile mapping within each CRF group × year cell.
#
#   The GPA shape (within-group emission dispersion) is estimated from pooled
#   within-NACE-2d-year log-emission deviations of ETS emitters. This shape
#   is applied to CRF group × year cells at deployment — same approximation
#   as assuming the pooled GPA applies to any grouping of emitters.
#
#   Pipeline:
#     0. Pre-ETS backcast: for firms entering ETS after 2005, estimate
#        pre-entry emissions using sector-year FE and first-two-year anchor.
#        Subtract E_pre_ets from E_deploy; treat pre_ets firms as observed.
#     1. NIR calibration → E_deploy = E_NIR - E_ETS - E_pre_ets per CRF × year
#     2. Within each CRF group × year:
#        a. Rank deployment firms by proxy (sinh-transformed to levels)
#        b. Apply CV threshold: zero out firms below threshold percentile
#        c. Map remaining ranks to GPA quantiles → exp(w) shares
#        d. Calibrate to E_deploy
#        e. Apply 30 kt physics-based cap (replaces min-ETS constraint)
#     3. Combine ETS + pre_ets (observed) + deployment (Pareto-imputed)
#     4. Compute within-NACE-sector carbon productivity dispersion
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/deployment_panel.RData
#   {PROC_DATA}/nir_calibration_targets.RData
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# OUTPUT
#   {PROC_DATA}/b_loop_scope1_dispersion_pareto.RData
#
# RUNS ON: RMD (or local with downsampled data for testing)
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
source(file.path(REPO_DIR, "utils", "calibration_pareto.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# ── Parameters ────────────────────────────────────────────────────────────────
YEARS         <- 2005:2021
MIN_N_STATS   <- 3L
FALLBACK_THRESHOLD <- 0.62  # used only if draw_params not available
FALLBACK_GPA      <- TRUE   # set to FALSE once draw_params is generated on RMD
DEPLOY_CAP    <- 30000      # physics-based upper bound (30 kt CO2/year)

N_CORES_SET <- if (tolower(Sys.info()[["user"]]) == "jardang") {
  40L
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 2L)
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("  RQ1: CARBON PRODUCTIVITY DISPERSION (PARETO REDISTRIBUTION)\n")
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
# SECTION 2: Pre-ETS backcast
# =============================================================================
# Firms entering ETS after 2005 have unobserved emissions in earlier years.
# We backcast their emissions using their first two observed ETS years as
# anchor, adjusted by sector-year fixed effects in log emissions.
# These "pre_ets" emissions are subtracted from E_deploy and the firms are
# treated as observed emitters alongside ETS firms.
cat("── Pre-ETS backcast ─────────────────────────────────────────\n")

eutl_pos <- eutl %>% filter(!is.na(emissions), emissions > 0)

# First ETS year per firm
ets_entry <- eutl_pos %>%
  group_by(vat) %>%
  summarise(first_ets_year = min(year), .groups = "drop")

late_entrants <- ets_entry %>% filter(first_ets_year > min(YEARS))

# Anchor: average of first two observed emissions
anchor <- eutl_pos %>%
  inner_join(late_entrants, by = "vat") %>%
  group_by(vat) %>%
  arrange(year) %>%
  slice_head(n = 2) %>%
  summarise(
    e_anchor    = mean(emissions),
    anchor_year = min(year),
    .groups = "drop"
  ) %>%
  left_join(
    accounts %>%
      group_by(vat) %>%
      summarise(nace2d = names(which.max(table(nace2d))), .groups = "drop"),
    by = "vat"
  ) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group))

# Sector-year means of log emissions (for FE adjustment)
sector_year_mu <- eutl_pos %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d)) %>%
  group_by(nace2d, year) %>%
  summarise(mu = mean(log(emissions)), .groups = "drop")

# Build backcast table
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

# Aggregate E_pre_ets by CRF group × year
E_pre_ets <- pre_ets %>%
  group_by(crf_group, year) %>%
  summarise(E_pre_ets = sum(emissions), n_pre_ets = n(), .groups = "drop")

cat("  Late entrants:", nrow(late_entrants), "firms\n")
cat("  Backcasted firm-years:", nrow(pre_ets), "\n")
cat("  Total E_pre_ets:", round(sum(pre_ets$emissions) / 1e6, 1), "Mt\n\n")


# =============================================================================
# SECTION 3: NIR calibration targets (fixed across draws)
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
  left_join(E_pre_ets, by = c("crf_group", "year")) %>%
  mutate(E_ETS      = coalesce(E_ETS, 0),
         E_pre_ets  = coalesce(E_pre_ets, 0),
         E_deploy   = pmax(E_NIR - E_ETS - E_pre_ets, 0))

n_floored <- sum(E_deploy_panel$E_NIR < (E_deploy_panel$E_ETS + E_deploy_panel$E_pre_ets),
                  na.rm = TRUE)
if (n_floored > 0)
  cat("  WARNING:", n_floored, "group-years where E_ETS + E_pre_ets > E_NIR → floored to 0\n")

cat("  E_deploy built:", nrow(E_deploy_panel), "group-years\n\n")


# =============================================================================
# SECTION 3: GPA + threshold — draw-specific or global fallback
# =============================================================================
cat("── GPA + threshold parameters ──────────────────────────────\n")

draw_params_path <- file.path(PROC_DATA, "deployment_draw_params.RData")
if (file.exists(draw_params_path)) {
  load(draw_params_path)  # draw_params: list of B elements
  USE_DRAW_PARAMS <- TRUE
  thresholds <- sapply(draw_params, `[[`, "threshold_pct")
  cat("  Draw-specific params loaded:", length(draw_params), "draws\n")
  cat("  Threshold: mean =", round(mean(thresholds, na.rm = TRUE), 3),
      " sd =", round(sd(thresholds, na.rm = TRUE), 3), "\n")
} else {
  USE_DRAW_PARAMS <- FALSE
  cat("  draw_params not found — using global fallback\n")
  cat("  Fallback threshold:", FALLBACK_THRESHOLD, "\n")

  # Compute global GPA from all ETS emitters
  eutl_with_nace <- eutl %>%
    filter(emissions > 0) %>%
    left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
    filter(!is.na(nace2d))
  ref_dist <- build_reference_dist(eutl_with_nace$emissions,
                                    eutl_with_nace$year,
                                    eutl_with_nace$nace2d)
  fallback_gpa <- fit_gpa(ref_dist)
  if (is.null(fallback_gpa))
    stop("GPA fitting failed — cannot proceed")
  cat("  GPA (global): xi =", round(fallback_gpa["xi"], 4),
      " alpha =", round(fallback_gpa["alpha"], 4),
      " k =", round(fallback_gpa["k"], 4), "\n")
}
cat("\n")


# =============================================================================
# SECTION 4: Helper functions
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
  do_stats <- function(df) {
    ok_cp   <- df$scope1 > 0 & df$revenue > 0
    cp      <- df$revenue[ok_cp] / df$scope1[ok_cp]
    n_cp    <- length(cp)
    log_cp  <- if (n_cp >= 2L) log(cp) else NA_real_

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
# SECTION 5: Main loop — year outer, draw inner
# =============================================================================
cat("── Main loop: years × draws ─────────────────────────────────\n\n")

stats2d_all <- list()
stats5d_all <- list()
flags_all   <- list()
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

  # ETS firms with positive emissions
  ok_ets <- !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  ets_vats_t      <- eutl_t$vat[ok_ets]
  ets_emissions_t <- setNames(eutl_t$emissions[ok_ets], eutl_t$vat[ok_ets])

  # ETS emitters with CRF group assignment (firm-level, for combined ranking)
  ets_crf_t <- eutl_t[ok_ets, ] %>%
    left_join(accounts_t %>% distinct(vat, nace2d), by = "vat") %>%
    left_join(nace_crf, by = "nace2d") %>%
    filter(!is.na(crf_group)) %>%
    select(vat, emissions, crf_group)

  # Pre-ETS firms for this year (backcasted emissions, treated as observed)
  pre_ets_t <- pre_ets[pre_ets$year == t, ]
  if (nrow(pre_ets_t) > 0) {
    pre_ets_vats_t      <- pre_ets_t$vat
    pre_ets_emissions_t <- setNames(pre_ets_t$emissions, pre_ets_t$vat)
    pre_ets_crf_t       <- pre_ets_t %>% select(vat, emissions, crf_group)
  } else {
    pre_ets_vats_t      <- character(0)
    pre_ets_emissions_t <- setNames(numeric(0), character(0))
    pre_ets_crf_t       <- data.frame(vat = character(0), emissions = numeric(0),
                                       crf_group = character(0))
  }

  # ── Draw loop ──────────────────────────────────────────────────────────────
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr); library(lmom) })
  # Draw-specific or fallback params
  if (USE_DRAW_PARAMS) {
    draw_gpa_list <- lapply(draw_params, `[[`, "gpa_params")
    draw_thr_list <- sapply(draw_params, `[[`, "threshold_pct")
  } else {
    draw_gpa_list <- replicate(B, fallback_gpa, simplify = FALSE)
    draw_thr_list <- rep(FALLBACK_THRESHOLD, B)
  }

  clusterExport(cl, c("ets_vats_t", "ets_emissions_t", "ets_crf_t",
                       "pre_ets_vats_t", "pre_ets_emissions_t", "pre_ets_crf_t",
                       "proxy_t", "deploy_nace", "E_dep_t", "accounts_t",
                       "draw_gpa_list", "draw_thr_list",
                       "compute_stats_bt", "gini", "pct_ratio",
                       "pareto_weights",
                       "MIN_N_STATS", "DEPLOY_CAP", "t"),
                envir = environment())

  draw_results <- foreach(b = seq_len(B), .packages = c("dplyr", "lmom")) %dopar% {

    proxy_b     <- proxy_t[[b]]
    gpa_b       <- draw_gpa_list[[b]]
    threshold_b <- draw_thr_list[b]
    if (is.na(threshold_b)) threshold_b <- 0

    # Build deployment firms with CRF group info
    deploy_b <- data.frame(
      vat       = names(proxy_b),
      proxy_val = as.numeric(proxy_b),
      stringsAsFactors = FALSE
    ) %>%
      filter(proxy_val > 0) %>%
      left_join(deploy_nace, by = "vat") %>%
      filter(!is.na(crf_group)) %>%
      left_join(E_dep_t, by = "crf_group") %>%
      mutate(E_deploy = coalesce(E_deploy, 0))

    # All CRF groups with deployment, ETS, or pre-ETS firms
    all_crf_groups <- union(unique(deploy_b$crf_group),
                             union(unique(ets_crf_t$crf_group),
                                   unique(pre_ets_crf_t$crf_group)))

    # ── Pareto redistribution within each CRF group ────────────────────────
    deploy_imputed <- list()
    alloc_flags    <- list()

    for (cg in all_crf_groups) {
      cg_deploy  <- deploy_b[deploy_b$crf_group == cg, ]
      cg_ets     <- ets_crf_t[ets_crf_t$crf_group == cg, ]
      cg_pre_ets <- pre_ets_crf_t[pre_ets_crf_t$crf_group == cg, ]
      n_ets      <- nrow(cg_ets)
      n_pre_ets  <- nrow(cg_pre_ets)
      n_observed <- n_ets + n_pre_ets
      E_cg       <- if (nrow(cg_deploy) > 0) cg_deploy$E_deploy[1] else {
        ed <- E_dep_t$E_deploy[E_dep_t$crf_group == cg]
        if (length(ed) > 0) ed[1] else 0
      }

      # ── Pure observed: no deployment firms ─────────────────────────────
      if (nrow(cg_deploy) == 0) {
        if (n_observed > 0) {
          alloc_flags[[length(alloc_flags) + 1]] <- data.frame(
            crf_group = cg, alloc_flag = "pure_ets", stringsAsFactors = FALSE)
        }
        next
      }

      if (E_cg <= 0) next

      # ── Determine sector type ──────────────────────────────────────────
      is_mixed <- n_observed > 0

      # ── Iterative threshold + allocation ───────────────────────────────
      proxy_levels_full <- pmax(sinh(cg_deploy$proxy_val), 0)
      thr_current       <- threshold_b
      flag              <- if (is_mixed) "mixed_ok" else "pure_deploy"
      MAX_ITER          <- 20L

      for (iter in seq_len(MAX_ITER)) {
        proxy_levels <- proxy_levels_full  # reset each iteration

        # Apply threshold: zero out firms below percentile cutoff
        pos_mask <- proxy_levels > 0
        if (sum(pos_mask) >= 2 && thr_current > 0) {
          tau <- quantile(proxy_levels[pos_mask], probs = thr_current,
                          names = FALSE)
          proxy_levels[proxy_levels <= tau] <- 0
        }

        emit_mask      <- proxy_levels > 0
        n_deploy_emit  <- sum(emit_mask)

        if (n_deploy_emit == 0) {
          # Threshold wiped everyone — halve and retry
          if (thr_current > 1e-4) {
            thr_current <- thr_current / 2
            if (is_mixed) flag <- "mixed_lowered"
            next
          }
          # threshold = 0 and still no emitters → skip
          flag <- if (is_mixed) "mixed_no_deploy" else "pure_deploy"
          break
        }

        # Combined ranking: deployment ranks 1..n_deploy_emit,
        # observed (ETS + pre_ets) implicitly at n_deploy_emit+1..n_total
        n_total    <- n_observed + n_deploy_emit
        emit_idx   <- which(emit_mask)
        emit_proxy <- proxy_levels[emit_idx]
        ranks      <- rank(emit_proxy, ties.method = "average")

        # GPA weights
        pw <- pareto_weights(ranks, n_total, gpa_b)

        if (is.null(pw)) {
          # GPA failed → sinh fallback
          flag <- "gpa_fallback_sinh"
          total_sinh <- sum(sinh(cg_deploy$proxy_val[emit_idx]))
          if (total_sinh > 0) {
            emissions_b <- E_cg * sinh(cg_deploy$proxy_val[emit_idx]) / total_sinh
          } else {
            emissions_b <- rep(E_cg / n_deploy_emit, n_deploy_emit)
          }
        } else {
          emissions_b <- E_cg * pw / sum(pw)
        }

        # Check physics-based upper-bound constraint
        if (max(emissions_b) >= DEPLOY_CAP) {
          if (thr_current > 1e-4) {
            # Lower threshold to include more firms
            thr_current <- thr_current / 2
            flag <- if (is_mixed) "mixed_lowered" else "pure_deploy_lowered"
            next
          }
          # threshold ≈ 0, constraint still violated → cap
          flag <- if (is_mixed) "mixed_capped" else "pure_deploy_capped"
          cap  <- DEPLOY_CAP * (1 - 1e-6)
          over <- which(emissions_b > cap)
          freed <- sum(emissions_b[over] - cap)
          emissions_b[over] <- cap
          under <- which(emissions_b <= cap)
          if (length(under) > 0) {
            emissions_b[under] <- emissions_b[under] +
              freed * emissions_b[under] / sum(emissions_b[under])
          }
        }
        break  # constraint satisfied or capped
      }

      # Record flag
      alloc_flags[[length(alloc_flags) + 1]] <- data.frame(
        crf_group = cg, alloc_flag = flag, stringsAsFactors = FALSE)

      # Store imputed emissions (if any emitters survived)
      if (n_deploy_emit > 0) {
        deploy_imputed[[length(deploy_imputed) + 1]] <- data.frame(
          vat    = cg_deploy$vat[emit_idx],
          scope1 = as.numeric(emissions_b),
          euets  = 0L,
          stringsAsFactors = FALSE
        )
      }
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

    firms_bt <- bind_rows(ets_df, pre_ets_df, deploy_df) %>%
      group_by(vat) %>%
      summarise(scope1 = sum(scope1), euets = max(euets), .groups = "drop") %>%
      left_join(accounts_t %>% select(vat, nace2d, nace5d, revenue), by = "vat") %>%
      filter(!is.na(nace5d), scope1 > 0)

    st <- compute_stats_bt(firms_bt)

    flags_df <- if (length(alloc_flags) > 0) bind_rows(alloc_flags) else
      data.frame(crf_group = character(0), alloc_flag = character(0))

    list(
      stats2d    = st$stats2d %>% mutate(year = t, draw = b),
      stats5d    = st$stats5d %>% mutate(year = t, draw = b),
      alloc_flags = flags_df %>% mutate(year = t, draw = b)
    )
  }

  stopCluster(cl)

  for (b in seq_len(B)) {
    stats2d_all[[length(stats2d_all) + 1L]] <- draw_results[[b]]$stats2d
    stats5d_all[[length(stats5d_all) + 1L]] <- draw_results[[b]]$stats5d
    flags_all[[length(flags_all) + 1L]]     <- draw_results[[b]]$alloc_flags
  }

  # Save intermediate checkpoint for this year
  year_stats2d <- bind_rows(lapply(draw_results, `[[`, "stats2d"))
  year_stats5d <- bind_rows(lapply(draw_results, `[[`, "stats5d"))
  year_flags   <- bind_rows(lapply(draw_results, `[[`, "alloc_flags"))
  save(year_stats2d, year_stats5d, year_flags,
       file = file.path(PROC_DATA,
                        sprintf("b_loop_pareto_year_%d.RData", t)))

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  cat(sprintf("(%s s) — saved checkpoint\n", elapsed))
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 6: Aggregate across B draws
# =============================================================================
cat("── Aggregating across draws ─────────────────────────────────\n")

all_stats2d    <- bind_rows(stats2d_all)
all_stats5d    <- bind_rows(stats5d_all)
all_alloc_flags <- bind_rows(flags_all)

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
cat("  NACE 5-digit:", nrow(stats5d_summary), "sector-years\n")

# Allocation flag summary
flag_summary <- all_alloc_flags %>%
  count(alloc_flag, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))
cat("\n  Allocation flags:\n")
for (i in seq_len(nrow(flag_summary))) {
  cat(sprintf("    %-25s %6d (%5.1f%%)\n",
              flag_summary$alloc_flag[i], flag_summary$n[i], flag_summary$pct[i]))
}
cat("\n")


# =============================================================================
# SECTION 7: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "b_loop_scope1_dispersion_pareto.RData")

save(
  stats2d_summary,
  stats5d_summary,
  all_stats2d,
  all_stats5d,
  all_alloc_flags,
  file = OUT_PATH
)

cat("══════════════════════════════════════════════════════════════\n")
cat("Saved:", OUT_PATH, "\n")
cat("  stats2d_summary:", nrow(stats2d_summary), "rows\n")
cat("  stats5d_summary:", nrow(stats5d_summary), "rows\n")
cat("  Total time:", total_time, "min\n")
cat("══════════════════════════════════════════════════════════════\n")
