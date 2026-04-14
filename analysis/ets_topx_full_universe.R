###############################################################################
# analysis/ets_topx_full_universe.R
#
# PURPOSE
#   Merge ETS training firms with deployment firms to create the full universe
#   of Belgian firms. Rank all firms within NACE2d-year and CRF-year by q_i
#   and by revenue. Report how many ETS emitters fall in the top X%.
#
#   This is the real deployment scenario: can the proxy identify ETS emitters
#   among all Belgian firms in their sector?
#
# INPUT
#   {PROC_DATA}/repeated_cv_proxy_crf_asinh.RData  (200 repeats for training)
#   {PROC_DATA}/training_sample.RData              (ETS + zero-sector firms)
#   {PROC_DATA}/deployment_proxy_list.RData         (deployment proxy draws)
#   {PROC_DATA}/deployment_panel.RData              (nace5d for deployment)
#   {PROC_DATA}/extensive_margin_calibration.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData (FULL)
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# RUNS ON: RMD
###############################################################################

# -- Paths -------------------------------------------------------------------
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

cat("===================================================================\n")
cat("  Top-X% classification: full universe (ETS + deployment)\n")
cat("===================================================================\n\n")


# =============================================================================
# 1. Load logistic calibration
# =============================================================================
load(file.path(PROC_DATA, "extensive_margin_calibration.RData"))
q_fn   <- ext_q_function[["balanced"]]
q_star <- q_star_pooled[["balanced"]]
cat("q_star =", round(q_star, 4), "\n\n")


# =============================================================================
# 2. ETS / training firms: proxy from 200 repeated cross-fits
# =============================================================================
cat("-- Loading training proxy -------------------------------------------\n")

load(file.path(PROC_DATA, "training_sample.RData"))
ts <- training_sample; rm(training_sample)

load(file.path(PROC_DATA, "repeated_cv_proxy_crf_asinh.RData"))
train_proxy <- data.frame(
  vat        = firmyear_index$vat,
  year       = firmyear_index$year,
  p_i        = rowMeans(proxy_matrix > 0, na.rm = TRUE),
  proxy_mean = rowMeans(proxy_matrix, na.rm = TRUE),
  stringsAsFactors = FALSE
)
rm(proxy_matrix, firmyear_index, repeat_diagnostics, repeat_timing,
   repeated_cv_proxy_panel, syt, crf_group_map)

train_proxy$q_i <- q_fn(train_proxy$p_i, train_proxy$proxy_mean)

# Merge with training sample for emissions, nace, revenue
train_panel <- ts %>%
  inner_join(train_proxy, by = c("vat", "year")) %>%
  mutate(
    emitter = as.integer(emissions > 0),
    source  = ifelse(euets == 1, "ets", "zero_sector")
  ) %>%
  select(vat, year, nace2d, revenue, emitter, source, q_i, p_i, proxy_mean)

rm(ts, train_proxy)
cat("  Training firm-years:", nrow(train_panel), "\n")
cat("  ETS emitters:", sum(train_panel$source == "ets" & train_panel$emitter == 1), "\n\n")


# =============================================================================
# 3. Deployment firms: proxy from deployment_proxy_list
# =============================================================================
cat("-- Loading deployment proxy -----------------------------------------\n")

load(file.path(PROC_DATA, "deployment_proxy_list.RData"))
B <- length(proxy_list)
cat("  Draws:", B, "\n")

# Collapse to p_i, proxy_mean per (vat, year)
proxy_long <- bind_rows(proxy_list, .id = "draw")
rm(proxy_list)

deploy_proxy <- proxy_long %>%
  group_by(vat, year) %>%
  summarise(
    p_i        = mean(proxy > 0, na.rm = TRUE),
    proxy_mean = mean(proxy, na.rm = TRUE),
    .groups    = "drop"
  )
rm(proxy_long)

deploy_proxy$q_i <- q_fn(deploy_proxy$p_i, deploy_proxy$proxy_mean)
cat("  Deployment firm-years:", nrow(deploy_proxy), "\n")


# =============================================================================
# 4. Get NACE and revenue for deployment firms
# =============================================================================
cat("-- Loading deployment panel + annual accounts -----------------------\n")

load(file.path(PROC_DATA, "deployment_panel.RData"))
deploy_nace <- deployment_panel %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = substr(nace5d, 1, 2))
rm(deployment_panel)

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  select(vat, year, revenue)
rm(df_annual_accounts_selected_sample_key_variables)

deploy_panel <- deploy_proxy %>%
  left_join(deploy_nace, by = c("vat", "year")) %>%
  left_join(accounts, by = c("vat", "year")) %>%
  mutate(
    emitter = 0L,         # deployment firms have no observed emissions
    source  = "deploy"
  ) %>%
  select(vat, year, nace2d, revenue, emitter, source, q_i, p_i, proxy_mean)

rm(deploy_proxy, deploy_nace, accounts)
cat("  Deploy panel with NACE:", sum(!is.na(deploy_panel$nace2d)), "/",
    nrow(deploy_panel), "\n\n")


# =============================================================================
# 5. Stack into full universe
# =============================================================================
cat("-- Building full universe -------------------------------------------\n")

# Remove deployment firm-years that overlap with training (late entrants)
deploy_panel <- deploy_panel %>%
  anti_join(train_panel, by = c("vat", "year"))

universe <- bind_rows(train_panel, deploy_panel)
cat("  Total firm-years:", nrow(universe), "\n")
cat("  ETS emitters:", sum(universe$source == "ets" & universe$emitter == 1), "\n")
cat("  Zero-sector:", sum(universe$source == "zero_sector"), "\n")
cat("  Deployment:", sum(universe$source == "deploy"), "\n")
cat("  With NACE2d:", sum(!is.na(universe$nace2d)), "\n")
cat("  With revenue:", sum(!is.na(universe$revenue)), "\n\n")

# NACE-CRF crosswalk
nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(nace2d = "character")
) %>% select(nace2d, crf_group)

universe <- universe %>%
  left_join(nace_crf, by = "nace2d")


# =============================================================================
# 6. Compute percentile ranks
# =============================================================================
cat("-- Computing percentile ranks ---------------------------------------\n")

# Within CRF-year
universe <- universe %>%
  filter(!is.na(crf_group)) %>%
  group_by(crf_group, year) %>%
  mutate(
    n_crf_yr       = n(),
    q_pctile_crf   = rank(q_i, ties.method = "average") / n(),
    rev_pctile_crf = {
      r <- rep(NA_real_, n())
      ok <- !is.na(revenue)
      if (sum(ok) > 1) r[ok] <- rank(revenue[ok], ties.method = "average") / sum(ok)
      r
    }
  ) %>%
  ungroup()

# Within NACE2d-year
universe <- universe %>%
  filter(!is.na(nace2d)) %>%
  group_by(nace2d, year) %>%
  mutate(
    n_nace2d_yr       = n(),
    q_pctile_nace2d   = rank(q_i, ties.method = "average") / n(),
    rev_pctile_nace2d = {
      r <- rep(NA_real_, n())
      ok <- !is.na(revenue)
      if (sum(ok) > 1) r[ok] <- rank(revenue[ok], ties.method = "average") / sum(ok)
      r
    }
  ) %>%
  ungroup()

cat("  Done.\n\n")


# =============================================================================
# 7. Top-X% comparison
# =============================================================================

ets_emitters <- universe %>% filter(source == "ets", emitter == 1)
n_ets <- nrow(ets_emitters)

cat("===================================================================\n")
cat("  RESULTS: How many ETS emitters fall in top X%?\n")
cat("  ETS emitter firm-years:", n_ets, "\n")
cat("===================================================================\n\n")

cutoffs <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)

cat("--- Within CRF-year ---\n\n")
cat(sprintf("%8s %15s %15s\n", "Top X%", "q_i", "Revenue"))
cat(strrep("-", 42), "\n")
for (x in cutoffs) {
  threshold <- 1 - x
  n_q   <- sum(ets_emitters$q_pctile_crf >= threshold, na.rm = TRUE)
  n_rev <- sum(ets_emitters$rev_pctile_crf >= threshold, na.rm = TRUE)
  cat(sprintf("%7.0f%% %7d / %-5d %7d / %-5d\n",
              x * 100, n_q, n_ets, n_rev, n_ets))
}

cat("\n--- Within NACE2d-year ---\n\n")
cat(sprintf("%8s %15s %15s\n", "Top X%", "q_i", "Revenue"))
cat(strrep("-", 42), "\n")
for (x in cutoffs) {
  threshold <- 1 - x
  n_q   <- sum(ets_emitters$q_pctile_nace2d >= threshold, na.rm = TRUE)
  n_rev <- sum(ets_emitters$rev_pctile_nace2d >= threshold, na.rm = TRUE)
  cat(sprintf("%7.0f%% %7d / %-5d %7d / %-5d\n",
              x * 100, n_q, n_ets, n_rev, n_ets))
}

# Cell size summary
cat("\n--- Cell sizes for ETS emitters ---\n\n")
cat("CRF-year:\n")
cat("  Min:   ", min(ets_emitters$n_crf_yr), "\n")
cat("  Median:", median(ets_emitters$n_crf_yr), "\n")
cat("  Max:   ", max(ets_emitters$n_crf_yr), "\n")

cat("\nNACE2d-year:\n")
cat("  Min:   ", min(ets_emitters$n_nace2d_yr), "\n")
cat("  Median:", median(ets_emitters$n_nace2d_yr), "\n")
cat("  Max:   ", max(ets_emitters$n_nace2d_yr), "\n")

# Breakdown by sector type
cat("\n--- Breakdown: mixed vs non-mixed sectors ---\n\n")
mixed_nace <- c("17", "18", "19", "24", "25")
ets_mixed    <- ets_emitters %>% filter(nace2d %in% mixed_nace)
ets_nonmixed <- ets_emitters %>% filter(!(nace2d %in% mixed_nace))

cat("Mixed sectors (", nrow(ets_mixed), "ETS emitters):\n")
cat(sprintf("%8s %15s %15s\n", "Top X%", "q_i/NACE2d", "Rev/NACE2d"))
cat(strrep("-", 42), "\n")
for (x in cutoffs) {
  threshold <- 1 - x
  n_q   <- sum(ets_mixed$q_pctile_nace2d >= threshold, na.rm = TRUE)
  n_rev <- sum(ets_mixed$rev_pctile_nace2d >= threshold, na.rm = TRUE)
  cat(sprintf("%7.0f%% %7d / %-5d %7d / %-5d\n",
              x * 100, n_q, nrow(ets_mixed), n_rev, nrow(ets_mixed)))
}

cat("\nNon-mixed sectors (", nrow(ets_nonmixed), "ETS emitters):\n")
cat(sprintf("%8s %15s %15s\n", "Top X%", "q_i/NACE2d", "Rev/NACE2d"))
cat(strrep("-", 42), "\n")
for (x in cutoffs) {
  threshold <- 1 - x
  n_q   <- sum(ets_nonmixed$q_pctile_nace2d >= threshold, na.rm = TRUE)
  n_rev <- sum(ets_nonmixed$rev_pctile_nace2d >= threshold, na.rm = TRUE)
  cat(sprintf("%7.0f%% %7d / %-5d %7d / %-5d\n",
              x * 100, n_q, nrow(ets_nonmixed), n_rev, nrow(ets_nonmixed)))
}

cat("\n===================================================================\n")
cat("  Done.\n")
cat("===================================================================\n")
