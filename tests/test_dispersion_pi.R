###############################################################################
# tests/test_dispersion_pi.R
#
# Quick smoke test for b_loop_scope1_dispersion_pi.R
# Runs on 1 year, 3 draws, 5 reps to verify the pipeline works end-to-end.
###############################################################################

# -- Paths --------------------------------------------------------------------
REPO_DIR <- normalizePath("c:/Users/jota_/Documents/facts-emissions-across-network",
                           winslash = "/")
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "calibration_pareto.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)

# -- Override parameters for quick test ---------------------------------------
YEARS       <- 2015L          # single year
B_TEST      <- 3L             # only 3 draws
R_REPS      <- 5L             # only 5 reps per draw
MIN_N_STATS <- 3L
DEPLOY_CAP  <- 30000
SEED_OFFSET <- 999000L

cat("===================================================================\n")
cat("  TEST: Dispersion PI pipeline (1 year, B=3, R=5)\n")
cat("===================================================================\n\n")

# -- Load data ----------------------------------------------------------------
cat("-- Loading data -----------------------------------------------------\n")

load(file.path(PROC_DATA, "deployment_proxy_list.RData"))
B_full <- length(proxy_list)
proxy_list <- proxy_list[1:B_TEST]  # subset to 3 draws
B <- B_TEST
cat("  proxy_list: using", B, "of", B_full, "draws\n")

# GPA fallback (no draw_params available)
FALLBACK_THRESHOLD <- 0.62

# Synthetic OOS residuals
set.seed(42)
epsilon_pool <- rnorm(5000, mean = 0, sd = 2.5)
cat("  Synthetic residuals: n=5000, sd=2.5\n")

# ETS
load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)
cat("  EUTL firm-years:", nrow(eutl), "\n")

# Accounts
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)
cat("  Accounts firm-years:", nrow(accounts), "\n")

# Deployment panel + CRF
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

# NIR targets
load(file.path(PROC_DATA, "nir_calibration_targets.RData"))

E_ETS_group <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group)) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions, na.rm = TRUE), .groups = "drop")

# Pre-ETS backcast (simplified: recompute)
eutl_pos <- eutl %>% filter(!is.na(emissions), emissions > 0)
ets_entry <- eutl_pos %>%
  group_by(vat) %>%
  summarise(first_ets_year = min(year), .groups = "drop")
late_entrants <- ets_entry %>% filter(first_ets_year > min(YEARS))

anchor <- eutl_pos %>%
  inner_join(late_entrants, by = "vat") %>%
  group_by(vat) %>% arrange(year) %>% slice_head(n = 2) %>%
  summarise(e_anchor = mean(emissions), anchor_year = min(year), .groups = "drop") %>%
  left_join(accounts %>% group_by(vat) %>%
              summarise(nace2d = names(which.max(table(nace2d))), .groups = "drop"),
            by = "vat") %>%
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
  mu_T <- sector_year_mu$mu[sector_year_mu$nace2d == s & sector_year_mu$year == T_entry]
  if (length(mu_T) == 0) next
  pre_years <- accounts$year[accounts$vat == v & accounts$year < T_entry &
                              accounts$year %in% YEARS]
  if (length(pre_years) == 0) next
  for (yr in unique(pre_years)) {
    mu_t <- sector_year_mu$mu[sector_year_mu$nace2d == s & sector_year_mu$year == yr]
    e_hat <- if (length(mu_t) > 0) e_anch * exp(mu_t - mu_T) else e_anch
    pre_ets_list[[length(pre_ets_list) + 1]] <- data.frame(
      vat = v, year = yr, nace2d = s, crf_group = crf,
      emissions = e_hat, stringsAsFactors = FALSE)
  }
}
pre_ets <- if (length(pre_ets_list) > 0) bind_rows(pre_ets_list) else
  data.frame(vat = character(0), year = integer(0), nace2d = character(0),
             crf_group = character(0), emissions = numeric(0))
E_pre_ets <- pre_ets %>%
  group_by(crf_group, year) %>%
  summarise(E_pre_ets = sum(emissions), n_pre_ets = n(), .groups = "drop")

E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(E_ETS_group, by = c("crf_group", "year")) %>%
  left_join(E_pre_ets, by = c("crf_group", "year")) %>%
  mutate(E_ETS     = coalesce(E_ETS, 0),
         E_pre_ets = coalesce(E_pre_ets, 0),
         E_deploy  = pmax(E_NIR - E_ETS - E_pre_ets, 0))

# GPA fallback
eutl_with_nace <- eutl %>%
  filter(emissions > 0) %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d))
fallback_ref <- build_reference_dist(eutl_with_nace$emissions,
                                      eutl_with_nace$year,
                                      eutl_with_nace$nace2d)
fallback_gpa <- fit_gpa(fallback_ref)
stopifnot(!is.null(fallback_gpa))
cat("  GPA: xi =", round(fallback_gpa["xi"], 4), "\n")

draw_params <- replicate(B, list(gpa_params = fallback_gpa,
                                   threshold_pct = FALLBACK_THRESHOLD),
                          simplify = FALSE)
cat("  Data loading complete\n\n")


# -- Helper functions ---------------------------------------------------------
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

allocate_crf_group <- function(deploy_vats, proxy_vals, n_observed,
                                E_deploy, gpa_params, threshold, deploy_cap) {
  if (length(proxy_vals) == 0 || E_deploy <= 0) return(NULL)
  proxy_levels_full <- pmax(sinh(proxy_vals), 0)
  thr_current <- threshold
  MAX_ITER <- 20L
  for (iter in seq_len(MAX_ITER)) {
    proxy_levels <- proxy_levels_full
    pos_mask <- proxy_levels > 0
    if (sum(pos_mask) >= 2 && thr_current > 0) {
      tau <- quantile(proxy_levels[pos_mask], probs = thr_current, names = FALSE)
      proxy_levels[proxy_levels <= tau] <- 0
    }
    emit_mask <- proxy_levels > 0
    n_deploy_emit <- sum(emit_mask)
    if (n_deploy_emit == 0) {
      if (thr_current > 1e-4) { thr_current <- thr_current / 2; next }
      return(NULL)
    }
    n_total <- n_observed + n_deploy_emit
    emit_idx <- which(emit_mask)
    emit_proxy <- proxy_levels[emit_idx]
    ranks <- rank(emit_proxy, ties.method = "average")
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
    if (max(emissions_out) >= deploy_cap) {
      if (thr_current > 1e-4) { thr_current <- thr_current / 2; next }
      cap <- deploy_cap * (1 - 1e-6)
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
  data.frame(vat = deploy_vats[emit_idx], scope1 = as.numeric(emissions_out),
             euets = 0L, stringsAsFactors = FALSE)
}


# -- Run on single year, sequentially -----------------------------------------
cat("-- Running double loop (sequential) --------------------------------\n")
t <- YEARS[1]
t0 <- Sys.time()

eutl_t  <- eutl[eutl$year == t, ]
E_dep_t <- E_deploy_panel[E_deploy_panel$year == t, c("crf_group", "E_deploy")]
accounts_t <- accounts[accounts$year == t, ]

ok_ets <- !is.na(eutl_t$emissions) & eutl_t$emissions > 0
ets_vats_t      <- eutl_t$vat[ok_ets]
ets_emissions_t <- setNames(eutl_t$emissions[ok_ets], eutl_t$vat[ok_ets])

ets_crf_t <- eutl_t[ok_ets, ] %>%
  left_join(accounts_t %>% distinct(vat, nace2d), by = "vat") %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group)) %>%
  select(vat, emissions, crf_group)

pre_ets_t <- pre_ets[pre_ets$year == t, ]
pre_ets_vats_t      <- pre_ets_t$vat
pre_ets_emissions_t <- setNames(pre_ets_t$emissions, pre_ets_t$vat)
pre_ets_crf_t       <- pre_ets_t %>% select(vat, emissions, crf_group)

proxy_t <- lapply(proxy_list, function(px) {
  sub <- px[px$year == t, c("vat", "proxy")]
  setNames(sub$proxy, sub$vat)
})

draw_gpa_list <- lapply(draw_params, `[[`, "gpa_params")
draw_thr_list <- sapply(draw_params, `[[`, "threshold_pct")

stats_all <- list()

for (b in seq_len(B)) {
  for (r in 0:R_REPS) {
    proxy_b <- proxy_t[[b]]
    gpa_b <- draw_gpa_list[[b]]
    threshold_b <- draw_thr_list[b]

    if (r > 0) {
      set.seed(SEED_OFFSET + t * 1000L + b * 100L + r)
      eps <- sample(epsilon_pool, length(proxy_b), replace = TRUE)
      proxy_b <- pmax(proxy_b + eps, 0)
    }

    deploy_b <- data.frame(
      vat = names(proxy_b), proxy_val = as.numeric(proxy_b),
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

    deploy_imputed <- list()
    for (cg in all_crf_groups) {
      cg_deploy  <- deploy_b[deploy_b$crf_group == cg, ]
      cg_ets     <- ets_crf_t[ets_crf_t$crf_group == cg, ]
      cg_pre_ets <- pre_ets_crf_t[pre_ets_crf_t$crf_group == cg, ]
      n_observed <- nrow(cg_ets) + nrow(cg_pre_ets)
      E_cg <- if (nrow(cg_deploy) > 0) cg_deploy$E_deploy[1] else {
        ed <- E_dep_t$E_deploy[E_dep_t$crf_group == cg]
        if (length(ed) > 0) ed[1] else 0
      }
      if (nrow(cg_deploy) == 0 || E_cg <= 0) next
      res <- allocate_crf_group(cg_deploy$vat, cg_deploy$proxy_val,
                                 n_observed, E_cg, gpa_b, threshold_b, DEPLOY_CAP)
      if (!is.null(res)) deploy_imputed[[length(deploy_imputed) + 1]] <- res
    }

    ets_df <- data.frame(vat = ets_vats_t, scope1 = ets_emissions_t[ets_vats_t],
                          euets = 1L, stringsAsFactors = FALSE)
    pre_ets_df <- if (length(pre_ets_vats_t) > 0) {
      data.frame(vat = pre_ets_vats_t, scope1 = pre_ets_emissions_t[pre_ets_vats_t],
                 euets = 0L, stringsAsFactors = FALSE)
    } else data.frame(vat = character(0), scope1 = numeric(0), euets = integer(0))
    deploy_df <- if (length(deploy_imputed) > 0) bind_rows(deploy_imputed) else
      data.frame(vat = character(0), scope1 = numeric(0), euets = integer(0))

    firms_btr <- bind_rows(ets_df, pre_ets_df, deploy_df) %>%
      group_by(vat) %>%
      summarise(scope1 = sum(scope1), euets = max(euets), .groups = "drop") %>%
      filter(scope1 > 0) %>%
      left_join(accounts_t %>% select(vat, nace2d, nace5d, revenue), by = "vat") %>%
      filter(!is.na(nace5d))

    stats2d <- firms_btr %>%
      filter(scope1 > 0) %>%
      group_by(nace2d) %>%
      filter(n() >= MIN_N_STATS) %>%
      group_modify(~compute_stats_sector(.x)) %>%
      ungroup() %>%
      mutate(year = t, draw = b, rep = r)

    stats_all[[length(stats_all) + 1]] <- stats2d
  }
  cat(sprintf("  Draw %d complete (%d reps)\n", b, R_REPS + 1))
}

elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
pi_stats2d <- bind_rows(stats_all)
cat(sprintf("\nDone in %s s | %d stat rows\n\n", elapsed, nrow(pi_stats2d)))


# -- Validate output structure ------------------------------------------------
cat("-- Validation -------------------------------------------------------\n")

n_draws <- n_distinct(pi_stats2d$draw)
n_reps  <- n_distinct(pi_stats2d$rep)
n_sectors <- n_distinct(pi_stats2d$nace2d)

cat("  Draws:", n_draws, "(expected", B, ")\n")
cat("  Reps per draw:", n_reps, "(expected", R_REPS + 1, "including r=0)\n")
cat("  Sectors:", n_sectors, "\n")
cat("  Total rows:", nrow(pi_stats2d), "\n")

stopifnot(n_draws == B)
stopifnot(n_reps == R_REPS + 1)
cat("  PASS: draw and rep counts match\n")

# Check that perturbed stats differ from unperturbed
r0_stats <- pi_stats2d %>% filter(rep == 0)
rN_stats <- pi_stats2d %>% filter(rep > 0)

cat("\n  Unperturbed (r=0) Gini by sector:\n")
for (i in seq_len(min(5, nrow(r0_stats)))) {
  cat(sprintf("    NACE %s, draw %d: Gini = %.4f\n",
              r0_stats$nace2d[i], r0_stats$draw[i], r0_stats$s1_gini[i]))
}

cat("\n  Perturbed vs unperturbed Gini (sector", r0_stats$nace2d[1], "):\n")
sec1 <- r0_stats$nace2d[1]
gini_r0 <- r0_stats$s1_gini[r0_stats$nace2d == sec1 & r0_stats$draw == 1]
gini_rN <- rN_stats$s1_gini[rN_stats$nace2d == sec1 & rN_stats$draw == 1]
cat(sprintf("    r=0: %.4f | r>0 mean: %.4f, sd: %.4f\n",
            gini_r0, mean(gini_rN, na.rm = TRUE), sd(gini_rN, na.rm = TRUE)))

# Variance decomposition for one sector
cat("\n  Quick variance decomposition (sector", sec1, "):\n")
sec_data <- pi_stats2d %>% filter(nace2d == sec1)
total_var <- var(sec_data$s1_gini, na.rm = TRUE)
draw_means <- sec_data %>%
  group_by(draw) %>%
  summarise(m = mean(s1_gini, na.rm = TRUE), .groups = "drop")
source1 <- var(draw_means$m)
draw_vars <- sec_data %>%
  group_by(draw) %>%
  summarise(v = var(s1_gini, na.rm = TRUE), .groups = "drop")
source2 <- mean(draw_vars$v)
cat(sprintf("    Total var: %.6f\n", total_var))
cat(sprintf("    Source 1 (estimation): %.6f (%.1f%%)\n",
            source1, 100 * source1 / (source1 + source2)))
cat(sprintf("    Source 2 (prediction): %.6f (%.1f%%)\n",
            source2, 100 * source2 / (source1 + source2)))

cat("\n===================================================================\n")
cat("  TEST PASSED\n")
cat("===================================================================\n")
