###############################################################################
# analysis/b_allocation_glo.R
#
# PURPOSE
#   Assign scope 1 emissions to all Belgian firms for each year, using a
#   SINGLE GLO allocation (not B draws). The 200 elastic-net repeats are
#   collapsed to two firm-year-level statistics:
#
#     p_i           = share of repeats with positive proxy
#     proxy_mean_i  = mean proxy across repeats
#
#   For each (CRF group x year) cell:
#     1. Pre-ETS backcast (lifted from b_allocation_pareto.R).
#     2. NIR calibration: E_deploy = E_NIR - E_ETS - E_pre_ets.
#     3. Baseline emitter classification: Dhat_i = 1{p_i >= tau_star}.
#     4. Joint within-cell ranking (rank 1 = highest):
#          a. true ETS firms ranked by observed emissions (desc)
#          b. pre-ETS firms ranked by backcasted emissions (desc)
#          c. imputed firms (Dhat_i = 1) ranked by proxy_mean_i (desc)
#     5. Compute exp_w_k = glo_weights(rank_k, n_cell, glo_par) for the joint
#        cell using a single fitted GLO.
#     6. For imputed firms only:
#            share_i = exp_w_i / sum_{j in imputed} exp_w_j
#            emis_i  = share_i * E_deploy
#        ETS and pre_ets firms keep their observed/backcasted emissions.
#
#   No physics-based cap is applied -- the diagnostic script counts how often
#   imputed levels exceed reasonable bounds.
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData
#   {PROC_DATA}/firm_year_belgian_euets.RData   (provides nace5d for ETS firms)
#   {PROC_DATA}/deployment_panel.RData          (provides nace5d for non-ETS)
#   {PROC_DATA}/nir_calibration_targets.RData
#   {PROC_DATA}/glo_reference_params.RData
#   {PROC_DATA}/pooled_youden_tau.RData
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# NACE source for the CRF mapping:
#   - ETS firms (true ETS + pre-ETS backcast): nace5d on
#     firm_year_belgian_euets (which itself is sourced from annual accounts).
#     Coverage is structurally ~94% of distinct ETS vats; firms with NA
#     nace2d are dropped.
#   - Non-ETS firms (deployment): nace5d on deployment_panel (annual accounts).
#   This avoids any dependency on annual_accounts_selected_sample_key_variables
#   on local 1, which is downsampled.
#
# OUTPUT
#   {PROC_DATA}/allocation_glo/alloc_YYYY.RData  (one file per year)
#     year_firms : data.frame with columns
#       vat, year, crf_group, scope1, source ("ets"/"pre_ets"/"imputed"),
#       rank_in_cell, n_cell, n_imputed_in_cell, p_i, proxy_mean_i
#
# RUNS ON: RMD or local 1 (downsampled).
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
source(file.path(REPO_DIR, "utils", "glo_fit.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

suppressPackageStartupMessages({
  library(dplyr)
})

# -- Parameters ---------------------------------------------------------------
YEARS <- 2005:2021

OUT_DIR <- file.path(PROC_DATA, "allocation_glo")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("===================================================================\n")
cat("  GLO EMISSION ALLOCATION (single GLO, single allocation)\n")
cat("  Years:", min(YEARS), "--", max(YEARS), "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load data
# =============================================================================
cat("-- Loading data -----------------------------------------------------\n")

load(file.path(PROC_DATA, "deployment_proxy_list.RData"))
B <- length(proxy_list)
cat("  proxy_list: B =", B, "draws\n")

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
# NACE for ETS firms is sourced from firm_year_belgian_euets$nace5d
# (which itself was joined from annual accounts upstream). This avoids a
# dependency on the downsampled annual_accounts file on local 1, and matches
# the inferring_emissions convention of using annual-accounts NACE for ETS
# firms (NOT nace_id_from_eutl, which is the installation-level NACE).
# Coverage is structurally ~94% of distinct ETS vats; firms with NA nace5d
# are dropped from cell ranking and from E_ETS with a warning below.
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  select(vat, year, emissions, nace2d)
rm(firm_year_belgian_euets)

n_ets_na_nace <- sum(is.na(eutl$nace2d))
cat("  EUTL firm-years:", nrow(eutl),
    "  (with NA nace2d:", n_ets_na_nace, ")\n")

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

load(file.path(PROC_DATA, "glo_reference_params.RData"))
cat("  GLO params: xi =", round(glo_par["xi"], 4),
    " alpha =", round(glo_par["alpha"], 4),
    " k =", round(glo_par["k"], 4), "\n")

load(file.path(PROC_DATA, "pooled_youden_tau.RData"))
cat("  tau_star (pooled-mixed Youden) =", round(tau_star, 4), "\n\n")


# =============================================================================
# SECTION 2: Collapse 200 repeats to firm-year statistics
# =============================================================================
cat("-- Collapsing repeats to (p_i, proxy_mean_i) ------------------------\n")

# proxy_list[[b]] is data.frame(vat, year, proxy)
proxy_long <- bind_rows(proxy_list, .id = "draw")
cat("  Stacked proxy rows:", nrow(proxy_long), "\n")

proxy_summary <- proxy_long %>%
  group_by(vat, year) %>%
  summarise(
    p_i          = mean(proxy > 0, na.rm = TRUE),
    proxy_mean_i = mean(proxy,     na.rm = TRUE),
    .groups = "drop"
  )
rm(proxy_long)
cat("  Distinct (vat, year):", nrow(proxy_summary), "\n\n")


# =============================================================================
# SECTION 3: Pre-ETS backcast (lifted from b_allocation_pareto.R)
# =============================================================================
cat("-- Pre-ETS backcast -------------------------------------------------\n")

eutl_pos <- eutl %>% filter(!is.na(emissions), emissions > 0)

ets_entry <- eutl_pos %>%
  group_by(vat) %>%
  summarise(first_ets_year = min(year), .groups = "drop")

late_entrants <- ets_entry %>% filter(first_ets_year > min(YEARS))

# nace2d for the anchor comes directly from eutl (which itself sourced it
# from annual accounts upstream). Take the modal non-NA nace2d per vat.
vat_nace2d <- eutl %>%
  filter(!is.na(nace2d)) %>%
  group_by(vat) %>%
  summarise(nace2d = names(which.max(table(nace2d))), .groups = "drop")

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
  left_join(vat_nace2d, by = "vat") %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group))

# Sector-year mean of log emissions for the FE adjustment in the backcast.
# Built from ETS firms with non-NA nace2d (from EUTL).
sector_year_mu <- eutl_pos %>%
  filter(!is.na(nace2d)) %>%
  group_by(nace2d, year) %>%
  summarise(mu = mean(log(emissions)), .groups = "drop")

# Years to backcast for each late entrant: all YEARS strictly before its
# ETS entry. (Original used annual_accounts firm-year coverage to limit this
# to years the firm appeared in accounts; that source is downsampled on
# local 1 and structurally incomplete for ETS firms anyway. We assume the
# firm exists for all YEARS prior to its ETS entry.)
pre_ets_list <- list()
for (i in seq_len(nrow(anchor))) {
  v <- anchor$vat[i]; s <- anchor$nace2d[i]; T_entry <- anchor$anchor_year[i]
  e_anch <- anchor$e_anchor[i]; crf <- anchor$crf_group[i]
  mu_T <- sector_year_mu$mu[sector_year_mu$nace2d == s &
                             sector_year_mu$year == T_entry]
  if (length(mu_T) == 0) next
  pre_years <- seq(min(YEARS), T_entry - 1L)
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

cat("  Late entrants:", nrow(late_entrants), "firms\n")
cat("  Backcasted firm-years:", nrow(pre_ets), "\n")
cat("  Total E_pre_ets:", round(sum(pre_ets$emissions) / 1e6, 1), "Mt\n\n")


# =============================================================================
# SECTION 4: NIR calibration targets
# =============================================================================
cat("-- NIR calibration targets ------------------------------------------\n")

load(file.path(PROC_DATA, "nir_calibration_targets.RData"))

E_ETS_group <- eutl %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group), year %in% YEARS) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions, na.rm = TRUE), .groups = "drop")

E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(E_ETS_group, by = c("crf_group", "year")) %>%
  left_join(E_pre_ets, by = c("crf_group", "year")) %>%
  mutate(E_ETS     = coalesce(E_ETS, 0),
         E_pre_ets = coalesce(E_pre_ets, 0),
         E_deploy  = pmax(E_NIR - E_ETS - E_pre_ets, 0))

n_floored <- sum(E_deploy_panel$E_NIR <
                   (E_deploy_panel$E_ETS + E_deploy_panel$E_pre_ets), na.rm = TRUE)
if (n_floored > 0)
  cat("  WARNING:", n_floored,
      "group-years where E_ETS + E_pre_ets > E_NIR (floored to 0)\n")

cat("  E_deploy built:", nrow(E_deploy_panel), "group-years\n\n")


# =============================================================================
# SECTION 5: Per-year allocation
# =============================================================================
cat("-- Allocation -------------------------------------------------------\n")

# Helper: allocate one (crf_group x year) cell.
# Returns a data.frame with one row per firm in the cell (ETS, pre_ets, imputed).
allocate_cell <- function(cg, yr,
                          ets_cg, preets_cg, imputed_cg, E_dep) {
  # ETS: rank by observed emissions, descending
  ets_sorted <- ets_cg[order(-ets_cg$emissions), , drop = FALSE]
  # pre-ETS: rank by backcasted emissions, descending
  pre_sorted <- preets_cg[order(-preets_cg$emissions), , drop = FALSE]
  # imputed (Dhat = 1): rank by proxy_mean_i, descending
  imp_sorted <- imputed_cg[order(-imputed_cg$proxy_mean_i), , drop = FALSE]

  n_ets     <- nrow(ets_sorted)
  n_preets  <- nrow(pre_sorted)
  n_imputed <- nrow(imp_sorted)
  n_cell    <- n_ets + n_preets + n_imputed

  if (n_cell == 0) return(NULL)

  # Joint ranks: ETS first (1..n_ets), then pre-ETS, then imputed.
  ets_ranks <- if (n_ets > 0)     seq_len(n_ets)                                    else integer(0)
  pre_ranks <- if (n_preets > 0)  n_ets + seq_len(n_preets)                         else integer(0)
  imp_ranks <- if (n_imputed > 0) n_ets + n_preets + seq_len(n_imputed)             else integer(0)

  # Imputed shares: GLO weights at joint ranks, normalized within imputed only.
  if (n_imputed > 0 && E_dep > 0) {
    exp_w_imp <- glo_weights(imp_ranks, n_cell, glo_par)
    if (any(is.na(exp_w_imp)) || sum(exp_w_imp) == 0) {
      # Fallback: equal split
      imp_share <- rep(1 / n_imputed, n_imputed)
    } else {
      imp_share <- exp_w_imp / sum(exp_w_imp)
    }
    imp_emis <- E_dep * imp_share
  } else {
    imp_emis <- numeric(n_imputed)
  }

  rows <- list()
  if (n_ets > 0) {
    rows[[length(rows) + 1]] <- data.frame(
      vat               = ets_sorted$vat,
      year              = yr,
      crf_group         = cg,
      scope1            = ets_sorted$emissions,
      source            = "ets",
      rank_in_cell      = ets_ranks,
      n_cell            = n_cell,
      n_imputed_in_cell = n_imputed,
      p_i               = NA_real_,
      proxy_mean_i      = NA_real_,
      stringsAsFactors  = FALSE
    )
  }
  if (n_preets > 0) {
    rows[[length(rows) + 1]] <- data.frame(
      vat               = pre_sorted$vat,
      year              = yr,
      crf_group         = cg,
      scope1            = pre_sorted$emissions,
      source            = "pre_ets",
      rank_in_cell      = pre_ranks,
      n_cell            = n_cell,
      n_imputed_in_cell = n_imputed,
      p_i               = NA_real_,
      proxy_mean_i      = NA_real_,
      stringsAsFactors  = FALSE
    )
  }
  if (n_imputed > 0) {
    rows[[length(rows) + 1]] <- data.frame(
      vat               = imp_sorted$vat,
      year              = yr,
      crf_group         = cg,
      scope1            = imp_emis,
      source            = "imputed",
      rank_in_cell      = imp_ranks,
      n_cell            = n_cell,
      n_imputed_in_cell = n_imputed,
      p_i               = imp_sorted$p_i,
      proxy_mean_i      = imp_sorted$proxy_mean_i,
      stringsAsFactors  = FALSE
    )
  }
  bind_rows(rows)
}

t0_total <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()
  cat(sprintf("Year %d ", t))

  # ---- ETS firms with positive emissions, mapped to CRF group ----
  eutl_t <- eutl[eutl$year == t & !is.na(eutl$emissions) & eutl$emissions > 0, ]
  ets_crf_t <- eutl_t %>%
    left_join(nace_crf, by = "nace2d") %>%
    filter(!is.na(crf_group)) %>%
    select(vat, emissions, crf_group)

  # ---- Pre-ETS firms for this year ----
  pre_ets_t <- pre_ets[pre_ets$year == t, c("vat", "emissions", "crf_group")]

  # ---- Imputed candidates: deployment firms with Dhat_i = 1 ----
  proxy_t <- proxy_summary[proxy_summary$year == t, ]
  imputed_t <- proxy_t %>%
    filter(p_i >= tau_star, !(vat %in% pre_ets_t$vat)) %>%
    left_join(deploy_nace %>% select(vat, crf_group), by = "vat") %>%
    filter(!is.na(crf_group)) %>%
    select(vat, crf_group, p_i, proxy_mean_i)

  # ---- E_deploy by CRF group for this year ----
  E_dep_t <- E_deploy_panel[E_deploy_panel$year == t,
                             c("crf_group", "E_deploy")]

  all_crf_groups <- sort(unique(c(ets_crf_t$crf_group,
                                   pre_ets_t$crf_group,
                                   imputed_t$crf_group)))

  cell_rows <- list()
  for (cg in all_crf_groups) {
    ets_cg     <- ets_crf_t[ets_crf_t$crf_group == cg, ]
    pre_cg     <- pre_ets_t[pre_ets_t$crf_group == cg, ]
    imp_cg     <- imputed_t[imputed_t$crf_group == cg, ]
    E_dep_cg   <- E_dep_t$E_deploy[E_dep_t$crf_group == cg]
    if (length(E_dep_cg) == 0) E_dep_cg <- 0

    cr <- allocate_cell(cg, t, ets_cg, pre_cg, imp_cg, E_dep_cg)
    if (!is.null(cr)) cell_rows[[length(cell_rows) + 1]] <- cr
  }

  year_firms <- if (length(cell_rows) > 0) bind_rows(cell_rows) else
    data.frame(vat = character(0), year = integer(0),
               crf_group = character(0), scope1 = numeric(0),
               source = character(0), rank_in_cell = integer(0),
               n_cell = integer(0), n_imputed_in_cell = integer(0),
               p_i = numeric(0), proxy_mean_i = numeric(0))

  out_path <- file.path(OUT_DIR, sprintf("alloc_%d.RData", t))
  save(year_firms, file = out_path)

  dt <- round(as.numeric(difftime(Sys.time(), t0_year, units = "secs")), 1)
  cat(sprintf("[%4d firms, %3d cells] (%.1fs)\n",
              nrow(year_firms),
              length(unique(year_firms$crf_group)), dt))
}

cat("\nTotal time:",
    round(as.numeric(difftime(Sys.time(), t0_total, units = "mins")), 1), "min\n")
cat("Output dir:", OUT_DIR, "\n")
