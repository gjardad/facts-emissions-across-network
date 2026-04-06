###############################################################################
# analysis/b_allocation_proportional.R
#
# PURPOSE
#   Assign scope 1 emissions to all Belgian firms for each year and each of
#   B subsample draws, using PROPORTIONAL redistribution (shares proportional
#   to sinh-transformed EN proxy).
#
#   This is a companion to b_allocation_pareto.R. The only difference is
#   the redistribution step: instead of GPA quantile weights, emission shares
#   within each CRF group are proportional to sinh(proxy).
#
#   Pipeline per year × draw:
#     0. Pre-ETS backcast: identical to b_allocation_pareto.R
#     1. NIR calibration: E_deploy = E_NIR - E_ETS - E_pre_ets per CRF × year
#     2. Within each CRF group × year:
#        a. Rank deployment firms by proxy (sinh-transformed to levels)
#        b. Apply CV threshold: zero out firms below threshold percentile
#        c. Proportional shares: share_i = sinh(proxy_i) / sum(sinh(proxy_j))
#        d. Calibrate to E_deploy
#        e. Apply 30 kt physics-based cap
#     3. Combine ETS + pre_ets (observed) + deployment (proportional-imputed)
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/deployment_panel.RData
#   {PROC_DATA}/nir_calibration_targets.RData
#   {PROC_DATA}/deployment_draw_params.RData  (optional; only threshold used)
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# OUTPUT
#   {PROC_DATA}/allocation_proportional/alloc_YYYY.RData  (one file per year)
#     year_firms_by_draw : list of B data.frames, each with columns
#                          (vat, scope1, euets)
#     year_flags         : data.frame with allocation diagnostics
#                          (crf_group, alloc_flag, n_deploy_total,
#                           n_deploy_emit, thr_initial, thr_effective,
#                           year, draw)
#
# RUNS ON: RMD (or local with downsampled data for testing)
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
library(parallel)
library(doParallel)
library(foreach)

# -- Parameters ---------------------------------------------------------------
YEARS              <- 2005:2021
FALLBACK_THRESHOLD <- 0.62
DEPLOY_CAP         <- 30000  # physics-based upper bound (30 kt CO2/year)

N_CORES_SET <- if (tolower(Sys.info()[["user"]]) == "jardang") {
  40L
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 2L)
}

OUT_DIR <- file.path(PROC_DATA, "allocation_proportional")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("===================================================================\n")
cat("  PROPORTIONAL EMISSION ALLOCATION\n")
cat("  Years:", min(YEARS), "--", max(YEARS), "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load data
# =============================================================================
cat("-- Loading data -----------------------------------------------------\n")

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
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))
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
# SECTION 2: Pre-ETS backcast (identical to b_allocation_pareto.R)
# =============================================================================
cat("-- Pre-ETS backcast -------------------------------------------------\n")

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

# Aggregate E_pre_ets by CRF group x year
E_pre_ets <- pre_ets %>%
  group_by(crf_group, year) %>%
  summarise(E_pre_ets = sum(emissions), n_pre_ets = n(), .groups = "drop")

cat("  Late entrants:", nrow(late_entrants), "firms\n")
cat("  Backcasted firm-years:", nrow(pre_ets), "\n")
cat("  Total E_pre_ets:", round(sum(pre_ets$emissions) / 1e6, 1), "Mt\n\n")


# =============================================================================
# SECTION 3: NIR calibration targets (fixed across draws)
# =============================================================================
cat("-- NIR calibration targets ------------------------------------------\n")

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
  cat("  WARNING:", n_floored, "group-years where E_ETS + E_pre_ets > E_NIR (floored to 0)\n")

cat("  E_deploy built:", nrow(E_deploy_panel), "group-years\n\n")


# =============================================================================
# SECTION 4: Threshold parameters -- draw-specific or global fallback
# =============================================================================
cat("-- Threshold parameters ---------------------------------------------\n")

draw_params_path <- file.path(PROC_DATA, "deployment_draw_params.RData")
if (file.exists(draw_params_path)) {
  load(draw_params_path)
  USE_DRAW_PARAMS <- TRUE
  thresholds <- sapply(draw_params, `[[`, "threshold_pct")
  cat("  Draw-specific params loaded:", length(draw_params), "draws\n")
  cat("  Threshold: mean =", round(mean(thresholds, na.rm = TRUE), 3),
      " sd =", round(sd(thresholds, na.rm = TRUE), 3), "\n")
} else {
  USE_DRAW_PARAMS <- FALSE
  cat("  draw_params not found -- using global fallback\n")
  cat("  Fallback threshold:", FALLBACK_THRESHOLD, "\n")
}
cat("\n")


# =============================================================================
# SECTION 5: Main loop -- year outer, draw inner (parallelised)
# =============================================================================
cat("-- Main loop: years x draws ----------------------------------------\n\n")

t0_total <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()
  cat(sprintf("Year %d ", t))

  eutl_t     <- eutl[eutl$year == t, ]
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

  # ETS emitters with CRF group assignment
  ets_crf_t <- eutl_t[ok_ets, ] %>%
    left_join(accounts[accounts$year == t, ] %>% distinct(vat, nace2d), by = "vat") %>%
    left_join(nace_crf, by = "nace2d") %>%
    filter(!is.na(crf_group)) %>%
    select(vat, emissions, crf_group)

  # Pre-ETS firms for this year
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

  # -- Draw loop (parallelised) ------------------------------------------------
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr) })

  if (USE_DRAW_PARAMS) {
    draw_thr_list <- sapply(draw_params, `[[`, "threshold_pct")
  } else {
    draw_thr_list <- rep(FALLBACK_THRESHOLD, B)
  }

  clusterExport(cl, c("ets_vats_t", "ets_emissions_t", "ets_crf_t",
                       "pre_ets_vats_t", "pre_ets_emissions_t", "pre_ets_crf_t",
                       "proxy_t", "deploy_nace", "E_dep_t",
                       "draw_thr_list",
                       "DEPLOY_CAP", "t"),
                envir = environment())

  draw_results <- foreach(b = seq_len(B), .packages = "dplyr") %dopar% {

    proxy_b     <- proxy_t[[b]]
    threshold_b <- draw_thr_list[b]
    if (is.na(threshold_b)) threshold_b <- 0

    # Build deployment firms with CRF group info
    # Exclude pre_ets firms to avoid double-counting
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

    all_crf_groups <- union(unique(deploy_b$crf_group),
                             union(unique(ets_crf_t$crf_group),
                                   unique(pre_ets_crf_t$crf_group)))

    # -- Proportional redistribution within each CRF group --------------------
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

      # -- Pure observed: no deployment firms ----------------------------------
      if (nrow(cg_deploy) == 0) {
        if (n_observed > 0) {
          alloc_flags[[length(alloc_flags) + 1]] <- data.frame(
            crf_group      = cg,
            alloc_flag     = "pure_ets",
            n_deploy_total = 0L,
            n_deploy_emit  = 0L,
            thr_initial    = NA_real_,
            thr_effective  = NA_real_,
            stringsAsFactors = FALSE)
        }
        next
      }

      if (E_cg <= 0) next

      # -- Iterative threshold + allocation ------------------------------------
      proxy_levels_full <- pmax(sinh(cg_deploy$proxy_val), 0)
      thr_current       <- threshold_b
      flag              <- if (n_observed > 0) "mixed_ok" else "pure_deploy"
      MAX_ITER          <- 20L

      for (iter in seq_len(MAX_ITER)) {
        proxy_levels <- proxy_levels_full

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
          if (thr_current > 1e-4) {
            thr_current <- thr_current / 2
            if (n_observed > 0) flag <- "mixed_lowered"
            next
          }
          flag <- if (n_observed > 0) "mixed_no_deploy" else "pure_deploy"
          break
        }

        emit_idx   <- which(emit_mask)

        # *** PROPORTIONAL SHARES (key difference from Pareto) ***
        total_sinh <- sum(sinh(cg_deploy$proxy_val[emit_idx]))
        if (total_sinh > 0) {
          emissions_b <- E_cg * sinh(cg_deploy$proxy_val[emit_idx]) / total_sinh
        } else {
          emissions_b <- rep(E_cg / n_deploy_emit, n_deploy_emit)
        }

        # Physics-based upper-bound constraint
        if (max(emissions_b) >= DEPLOY_CAP) {
          if (thr_current > 1e-4) {
            thr_current <- thr_current / 2
            flag <- if (n_observed > 0) "mixed_lowered" else "pure_deploy_lowered"
            next
          }
          flag <- if (n_observed > 0) "mixed_capped" else "pure_deploy_capped"
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
        break
      }

      alloc_flags[[length(alloc_flags) + 1]] <- data.frame(
        crf_group      = cg,
        alloc_flag     = flag,
        n_deploy_total = nrow(cg_deploy),
        n_deploy_emit  = n_deploy_emit,
        thr_initial    = threshold_b,
        thr_effective  = thr_current,
        stringsAsFactors = FALSE)

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

    firms_b <- bind_rows(ets_df, pre_ets_df, deploy_df) %>%
      group_by(vat) %>%
      summarise(scope1 = sum(scope1), euets = max(euets), .groups = "drop") %>%
      filter(scope1 > 0)

    flags_df <- if (length(alloc_flags) > 0) bind_rows(alloc_flags) else
      data.frame(crf_group = character(0), alloc_flag = character(0),
                 n_deploy_total = integer(0), n_deploy_emit = integer(0),
                 thr_initial = numeric(0), thr_effective = numeric(0))

    list(firms = firms_b, flags = flags_df %>% mutate(year = t, draw = b))
  }

  stopCluster(cl)

  # Collect results
  year_firms_by_draw <- vector("list", B)
  year_flags_list    <- vector("list", B)
  for (b in seq_len(B)) {
    year_firms_by_draw[[b]] <- draw_results[[b]]$firms
    year_flags_list[[b]]    <- draw_results[[b]]$flags
  }
  year_flags <- bind_rows(year_flags_list)

  # Save checkpoint
  out_path <- file.path(OUT_DIR, sprintf("alloc_%d.RData", t))
  save(year_firms_by_draw, year_flags, file = out_path)

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  n_firms_avg <- round(mean(sapply(year_firms_by_draw, nrow)))
  cat(sprintf("(%s s | ~%d firms/draw) -- saved\n", elapsed, n_firms_avg))
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 6: Allocation flag summary (across all years)
# =============================================================================
cat("-- Allocation flag summary ------------------------------------------\n")

all_flags <- list()
for (t in YEARS) {
  f <- file.path(OUT_DIR, sprintf("alloc_%d.RData", t))
  if (file.exists(f)) {
    load(f)
    all_flags[[length(all_flags) + 1]] <- year_flags
  }
}
all_flags_df <- bind_rows(all_flags)

flag_summary <- all_flags_df %>%
  count(alloc_flag, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))

for (i in seq_len(nrow(flag_summary))) {
  cat(sprintf("  %-25s %6d (%5.1f%%)\n",
              flag_summary$alloc_flag[i], flag_summary$n[i], flag_summary$pct[i]))
}

cat(sprintf("\n===================================================================\n"))
cat("Saved to:", OUT_DIR, "\n")
cat("  ", length(YEARS), "year files (alloc_YYYY.RData)\n")
cat("  Total time:", total_time, "min\n")
cat("===================================================================\n")
