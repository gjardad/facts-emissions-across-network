###############################################################################
# analysis/uncertainty_propagation.R
#
# PURPOSE
#   Run B perturbation draws of the deployment imputation pipeline. For each
#   draw, perturb the extensive margin (Bernoulli flips around the deployed
#   threshold) and the intensive margin (latent-normal draws of the within-
#   cell rank), then re-allocate sector-year totals using the joint-cell GLO
#   rule consistent with b_allocation_glo.R.
#
#   See deployment_prediction_interval_step_by_step.md for the framework.
#
# INPUT
#   {PROC_DATA}/deployment_proxy_list.RData          (200 EN repeats)
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/deployment_panel.RData
#   {PROC_DATA}/nir_calibration_targets.RData
#   {PROC_DATA}/glo_reference_params.RData
#   {PROC_DATA}/threshold_calibration.RData
#   {PROC_DATA}/extensive_margin_calibration.RData
#   {PROC_DATA}/intensive_margin_calibration.RData
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# OUTPUT
#   {PROC_DATA}/uncertainty_draws/draw_BBBB.RData    (one file per draw)
#     draw_firms : data.frame with columns
#       vat, year, crf_group, scope1_b, source ("ets"/"pre_ets"/"imputed"),
#       Dtilde, draw
#
# RUNS ON: local 1
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
B_DRAWS       <- 200L
YEARS         <- 2005:2021
BASE_SEED     <- 2026L
WEIGHT_SCHEME <- "balanced"

OUT_DIR <- file.path(PROC_DATA, sprintf("uncertainty_draws_%s", WEIGHT_SCHEME))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("===================================================================\n")
cat("  UNCERTAINTY PROPAGATION (deployment perturbation)\n")
cat("  B =", B_DRAWS, "draws | years:", min(YEARS), "--", max(YEARS), "\n")
cat("  WEIGHT_SCHEME =", WEIGHT_SCHEME, "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load all inputs
# =============================================================================

cat("-- Loading data -----------------------------------------------------\n")

load(file.path(PROC_DATA, "deployment_proxy_list.RData"))
M_REPEATS <- length(proxy_list)
cat("  proxy_list: M =", M_REPEATS, "EN repeats\n")

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))
rm(df_annual_accounts_selected_sample_key_variables)

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

load(file.path(PROC_DATA, "glo_reference_params.RData"))
cat("  GLO params: xi =", round(glo_par["xi"], 4),
    " alpha =", round(glo_par["alpha"], 4),
    " k =", round(glo_par["k"], 4), "\n")

load(file.path(PROC_DATA, "extensive_margin_calibration.RData"))
stopifnot(WEIGHT_SCHEME %in% SCHEMES)
sch_coefs       <- ext_coefs[[WEIGHT_SCHEME]]
sch_q_function  <- ext_q_function[[WEIGHT_SCHEME]]
sch_q_star      <- q_star_pooled[[WEIGHT_SCHEME]]
cat("  ext logistic [", WEIGHT_SCHEME, "]: alpha =", round(sch_coefs["alpha"], 3),
    " beta_p =", round(sch_coefs["beta_p"], 3),
    " beta_ap =", round(sch_coefs["beta_ap"], 3), "\n")
cat("  q_star_pooled =", round(sch_q_star, 4), "\n")

load(file.path(PROC_DATA,
  sprintf("intensive_margin_calibration_%s.RData", WEIGHT_SCHEME)))
cat("  int model: a =", round(mean_coefs["a"], 3),
    " b =", round(mean_coefs["b"], 3),
    " sigma =", round(sigma_const, 3), "\n")
cat("  branch B: ", length(branch_b_z), " empirical z values\n\n")


# =============================================================================
# SECTION 2: Collapse 200 EN repeats to (p_i, proxy_mean_i)
# =============================================================================
cat("-- Collapsing repeats to (p_i, proxy_mean_i) ------------------------\n")

proxy_long <- bind_rows(proxy_list, .id = "draw")
rm(proxy_list); gc()

proxy_summary <- proxy_long %>%
  group_by(vat, year) %>%
  summarise(
    p_i          = mean(proxy > 0, na.rm = TRUE),
    proxy_mean_i = mean(proxy,     na.rm = TRUE),
    .groups = "drop"
  )
rm(proxy_long); gc()

cat("  Distinct (vat, year):", nrow(proxy_summary), "\n\n")


# =============================================================================
# SECTION 3: Pre-ETS backcast (lifted from b_allocation_glo.R)
# =============================================================================
cat("-- Pre-ETS backcast -------------------------------------------------\n")

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
  summarise(E_pre_ets = sum(emissions), .groups = "drop")

cat("  Backcasted pre-ETS firm-years:", nrow(pre_ets), "\n\n")


# =============================================================================
# SECTION 4: NIR calibration targets and E_deploy panel
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
  mutate(E_ETS     = coalesce(E_ETS, 0),
         E_pre_ets = coalesce(E_pre_ets, 0),
         E_deploy  = pmax(E_NIR - E_ETS - E_pre_ets, 0))

cat("  E_deploy built:", nrow(E_deploy_panel), "group-years\n\n")


# =============================================================================
# SECTION 5: Build the deployment imputed panel
# =============================================================================
cat("-- Building imputed panel -------------------------------------------\n")

# All firm-years that are NOT ETS firms in that year and NOT pre-ETS in that
# year. These are the "imputed" candidates. Among them, Dhat = 1 iff
# q_i >= q_star_pooled, where q_i is from the logistic on (p_i, log(1+proxy_mean)).

ets_vy <- eutl_pos %>% distinct(vat, year)
pre_ets_vy <- pre_ets %>% distinct(vat, year)

imputed_panel <- proxy_summary %>%
  anti_join(ets_vy,     by = c("vat", "year")) %>%
  anti_join(pre_ets_vy, by = c("vat", "year")) %>%
  left_join(deploy_nace %>% select(vat, crf_group), by = "vat") %>%
  filter(!is.na(crf_group)) %>%
  mutate(q_i  = sch_q_function(p_i, proxy_mean_i),
         Dhat = as.integer(q_i >= sch_q_star))

cat("  Imputed firm-years (all):", nrow(imputed_panel), "\n")
cat("  Baseline Dhat = 1:       ", sum(imputed_panel$Dhat == 1), "\n")
cat("  Baseline Dhat = 0:       ", sum(imputed_panel$Dhat == 0), "\n\n")


# =============================================================================
# SECTION 6: Precompute baseline uhat_i / zhat_i for Dhat = 1 imputed firms
# =============================================================================
cat("-- Computing baseline uhat_i / zhat_i for Dhat = 1 imputed firms ----\n")

# uhat_i = within-imputed percentile of proxy_mean_i (rank 1 = highest)
compute_pct <- function(x, n) {
  r <- rank(-x, ties.method = "average")
  (r - 0.5) / n
}

baseline_dhat1 <- imputed_panel %>%
  filter(Dhat == 1) %>%
  group_by(crf_group, year) %>%
  mutate(
    uhat_i = compute_pct(proxy_mean_i, n()),
    zhat_i = qnorm(uhat_i)
  ) %>%
  ungroup()

# Predicted mean and constant sigma for branch A draws
baseline_dhat1$mu_i    <- mean_coefs["a"] + mean_coefs["b"] * baseline_dhat1$zhat_i
baseline_dhat1$sigma_i <- sigma_const

cat("  Baseline Dhat = 1 with uhat_i computed:", nrow(baseline_dhat1), "\n\n")


# =============================================================================
# SECTION 7: Per-cell allocator for one draw
# =============================================================================
# Inputs to allocate_cell_draw:
#   ets_cg     : data.frame(vat, emissions) for ETS firms in this CRF-year
#   pre_cg     : data.frame(vat, emissions) for pre-ETS firms
#   imp_em     : data.frame(vat, ztilde) for SIMULATED imputed emitters
#                                                (only those with Dtilde = 1)
#   E_dep      : E_NIR - E_ETS - E_pre_ets, residual to allocate to imputed
#
# Output: data.frame with vat, scope1_b, source, rank_in_cell, n_cell

allocate_cell_draw <- function(cg, yr, ets_cg, pre_cg, imp_em, E_dep) {
  ets_sorted <- ets_cg[order(-ets_cg$emissions), , drop = FALSE]
  pre_sorted <- pre_cg[order(-pre_cg$emissions), , drop = FALSE]
  imp_sorted <- imp_em[order(-imp_em$ztilde), , drop = FALSE]

  n_ets     <- nrow(ets_sorted)
  n_preets  <- nrow(pre_sorted)
  n_imputed <- nrow(imp_sorted)
  n_cell    <- n_ets + n_preets + n_imputed

  if (n_cell == 0) return(NULL)

  ets_ranks <- if (n_ets > 0)     seq_len(n_ets)                            else integer(0)
  pre_ranks <- if (n_preets > 0)  n_ets + seq_len(n_preets)                 else integer(0)
  imp_ranks <- if (n_imputed > 0) n_ets + n_preets + seq_len(n_imputed)     else integer(0)

  if (n_imputed > 0 && E_dep > 0) {
    exp_w_imp <- glo_weights(imp_ranks, n_cell, glo_par)
    if (any(is.na(exp_w_imp)) || sum(exp_w_imp) == 0) {
      imp_share <- rep(1 / n_imputed, n_imputed)
    } else {
      imp_share <- exp_w_imp / sum(exp_w_imp)
    }
    imp_emis <- E_dep * imp_share
  } else {
    imp_emis <- numeric(n_imputed)
  }

  rows <- list()
  if (n_ets > 0) rows[[length(rows) + 1]] <- data.frame(
    vat = ets_sorted$vat, year = yr, crf_group = cg,
    scope1_b = ets_sorted$emissions, source = "ets",
    rank_in_cell = ets_ranks, n_cell = n_cell, stringsAsFactors = FALSE)
  if (n_preets > 0) rows[[length(rows) + 1]] <- data.frame(
    vat = pre_sorted$vat, year = yr, crf_group = cg,
    scope1_b = pre_sorted$emissions, source = "pre_ets",
    rank_in_cell = pre_ranks, n_cell = n_cell, stringsAsFactors = FALSE)
  if (n_imputed > 0) rows[[length(rows) + 1]] <- data.frame(
    vat = imp_sorted$vat, year = yr, crf_group = cg,
    scope1_b = imp_emis, source = "imputed",
    rank_in_cell = imp_ranks, n_cell = n_cell, stringsAsFactors = FALSE)

  bind_rows(rows)
}


# =============================================================================
# SECTION 8: Main perturbation loop
# =============================================================================
cat("-- Running B =", B_DRAWS, "draws ------------------------------------\n\n")

t0_total <- Sys.time()

for (b in seq_len(B_DRAWS)) {
  set.seed(BASE_SEED + b)
  t0_draw <- Sys.time()

  # ---- 8.1 Extensive margin: draw Dtilde for all imputed firms ----
  imputed_panel$Dtilde <- as.integer(runif(nrow(imputed_panel)) < imputed_panel$q_i)

  # ---- 8.2 Intensive margin: draw ztilde for simulated emitters ----
  # Two cases:
  #   (a) Dhat = 1, Dtilde = 1: use branch A latent-normal model.
  #   (b) Dhat = 0, Dtilde = 1: sample from branch B empirical distribution.

  sim_emit_a_idx <- which(imputed_panel$Dhat == 1 & imputed_panel$Dtilde == 1)
  sim_emit_b_idx <- which(imputed_panel$Dhat == 0 & imputed_panel$Dtilde == 1)

  ztilde_vec <- rep(NA_real_, nrow(imputed_panel))

  if (length(sim_emit_a_idx) > 0) {
    a_match <- match(
      paste(imputed_panel$vat[sim_emit_a_idx], imputed_panel$year[sim_emit_a_idx]),
      paste(baseline_dhat1$vat,             baseline_dhat1$year)
    )
    mu_a    <- baseline_dhat1$mu_i[a_match]
    sigma_a <- baseline_dhat1$sigma_i[a_match]
    ztilde_vec[sim_emit_a_idx] <- rnorm(length(sim_emit_a_idx),
                                         mean = mu_a, sd = sigma_a)
  }

  if (length(sim_emit_b_idx) > 0) {
    ztilde_vec[sim_emit_b_idx] <- sample(branch_b_z,
                                          size = length(sim_emit_b_idx),
                                          replace = TRUE)
  }

  imputed_panel$ztilde <- ztilde_vec

  # ---- 8.3 Per-year, per-CRF allocation ----
  draw_cell_rows <- list()

  for (t in YEARS) {
    eutl_t <- eutl[eutl$year == t & !is.na(eutl$emissions) & eutl$emissions > 0, ]
    ets_crf_t <- eutl_t %>%
      left_join(accounts[accounts$year == t, ] %>% distinct(vat, nace2d),
                by = "vat") %>%
      left_join(nace_crf, by = "nace2d") %>%
      filter(!is.na(crf_group)) %>%
      select(vat, emissions, crf_group)

    pre_ets_t <- pre_ets[pre_ets$year == t, c("vat", "emissions", "crf_group")]

    imp_t_simem <- imputed_panel[imputed_panel$year == t &
                                  imputed_panel$Dtilde == 1,
                                  c("vat", "crf_group", "ztilde")]

    E_dep_t <- E_deploy_panel[E_deploy_panel$year == t,
                               c("crf_group", "E_deploy")]

    all_cgs <- sort(unique(c(ets_crf_t$crf_group,
                             pre_ets_t$crf_group,
                             imp_t_simem$crf_group)))

    for (cg in all_cgs) {
      ets_cg <- ets_crf_t[ets_crf_t$crf_group == cg, ]
      pre_cg <- pre_ets_t[pre_ets_t$crf_group == cg, ]
      imp_cg <- imp_t_simem[imp_t_simem$crf_group == cg, ]
      E_dep_cg <- E_dep_t$E_deploy[E_dep_t$crf_group == cg]
      if (length(E_dep_cg) == 0) E_dep_cg <- 0

      cr <- allocate_cell_draw(cg, t, ets_cg, pre_cg, imp_cg, E_dep_cg)
      if (!is.null(cr)) draw_cell_rows[[length(draw_cell_rows) + 1]] <- cr
    }
  }

  draw_firms <- if (length(draw_cell_rows) > 0) bind_rows(draw_cell_rows) else
    data.frame(vat = character(0), year = integer(0),
               crf_group = character(0), scope1_b = numeric(0),
               source = character(0), rank_in_cell = integer(0),
               n_cell = integer(0))
  draw_firms$draw <- b

  out_path <- file.path(OUT_DIR, sprintf("draw_%04d.RData", b))
  save(draw_firms, file = out_path)

  dt <- round(as.numeric(difftime(Sys.time(), t0_draw, units = "secs")), 1)
  cat(sprintf("  Draw %3d / %d: %d firm-years (%.1fs)\n",
              b, B_DRAWS, nrow(draw_firms), dt))
}

cat("\nTotal time:",
    round(as.numeric(difftime(Sys.time(), t0_total, units = "mins")), 1), "min\n")
cat("Output dir:", OUT_DIR, "\n")
