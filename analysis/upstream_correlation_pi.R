###############################################################################
# analysis/upstream_correlation_pi.R
#
# PURPOSE
#   RQ3/4: within-sector correlation and rank correlation between scope 1 and
#   upstream embodied emissions, with prediction intervals from the B
#   perturbation draws plus the deterministic point estimate.
#
#   For each sector x year (at three granularities: NACE 2-digit, NACE 5-digit,
#   and CRF category):
#     - Point estimate: statistics from draw 0 (deterministic GLO allocation
#       propagated through the Leontief model)
#     - Prediction interval: 2.5th and 97.5th percentiles across B draws
#
# INPUT
#   {PROC_DATA}/upstream_emissions_glo_<scheme>/firms_YYYY.RData
#
# OUTPUT
#   {PROC_DATA}/upstream_correlation_pi_<scheme>.RData
#     corr_point_*  : point-estimate correlations by sector x year (from draw 0)
#     corr_draws_*  : all B draw-level correlation stats
#     corr_summary_*: point estimate + 90%/95% PIs per sector-year
#     (* = 2d, 5d, crf for the three granularities)
#
# RUNS ON: local 1 or RMD (reads upstream_emissions output, no B2B needed)
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

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# -- Parameters ---------------------------------------------------------------
WEIGHT_SCHEME <- "balanced"
YEARS         <- 2005:2021
MIN_N_CORR    <- 10L

UPSTREAM_DIR <- file.path(PROC_DATA,
  sprintf("upstream_emissions_glo_%s", WEIGHT_SCHEME))

cat("===================================================================\n")
cat("  UPSTREAM CORRELATION WITH PREDICTION INTERVALS\n")
cat("  WEIGHT_SCHEME =", WEIGHT_SCHEME, "\n")
cat("  Years:", min(YEARS), "--", max(YEARS), "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load accounts (NACE codes)
# =============================================================================

cat("-- Loading accounts -------------------------------------------------\n")

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))
rm(df_annual_accounts_selected_sample_key_variables)
setDT(accounts)

# CRF crosswalk
nace_crf <- fread(file.path(REPO_DIR, "preprocess", "crosswalks",
                             "nace_crf_crosswalk.csv"),
                   select = c("nace2d", "crf_group"),
                   colClasses = c(nace2d = "character"))
nace_crf <- rbind(nace_crf,
                   data.table(nace2d = "17/18", crf_group = "paper"))
accounts <- merge(accounts, nace_crf, by = "nace2d", all.x = TRUE)
cat("  Accounts firm-years:", nrow(accounts), "\n\n")


# =============================================================================
# SECTION 2: Correlation helper
# =============================================================================

compute_correlations <- function(firms_dt, group_col) {
  firms_dt <- firms_dt[!is.na(get(group_col))]

  firms_dt[scope1 > 0 & upstream > 0, {
    n <- .N
    if (n >= MIN_N_CORR) {
      .(n_firms        = n,
        spearman_s1_up = cor(scope1, upstream, method = "spearman"),
        pearson_s1_up  = cor(log(scope1), log(upstream), method = "pearson"))
    } else {
      .(n_firms        = n,
        spearman_s1_up = NA_real_,
        pearson_s1_up  = NA_real_)
    }
  }, by = group_col]
}


# =============================================================================
# SECTION 3: Process upstream files
# =============================================================================

cat("-- Processing upstream files ----------------------------------------\n")

GROUP_LEVELS <- c("nace2d", "nace5d", "crf_group")

point_lists <- setNames(lapply(GROUP_LEVELS, function(x) list()), GROUP_LEVELS)
draw_lists  <- setNames(lapply(GROUP_LEVELS, function(x) list()), GROUP_LEVELS)
t0 <- Sys.time()

for (t in YEARS) {
  up_path <- file.path(UPSTREAM_DIR, sprintf("firms_%d.RData", t))
  if (!file.exists(up_path)) {
    cat(sprintf("  Year %d: SKIPPED (no upstream file)\n", t))
    next
  }
  load(up_path)  # loads firms_by_draw, max_rowsum

  accounts_t <- accounts[year == t]
  n_draws <- length(firms_by_draw)

  for (b in seq_len(n_draws)) {
    firms_b <- as.data.table(firms_by_draw[[b]])
    firms_b <- merge(firms_b,
                      accounts_t[, .(vat, nace2d, nace5d, crf_group)],
                      by = "vat", all.x = TRUE)

    for (g in GROUP_LEVELS) {
      st <- compute_correlations(firms_b, g)
      if (nrow(st) == 0) next

      st[, year := t]

      if (b == 1L) {
        point_lists[[g]][[length(point_lists[[g]]) + 1]] <- st
      } else {
        st[, draw := b - 1L]
        draw_lists[[g]][[length(draw_lists[[g]]) + 1]] <- st
      }
    }
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("  Year %d: %d draws (%.1fs elapsed)\n", t, n_draws, elapsed))

  rm(firms_by_draw)
  gc()
}

corr_point_2d  <- rbindlist(point_lists[["nace2d"]])
corr_point_5d  <- rbindlist(point_lists[["nace5d"]])
corr_point_crf <- rbindlist(point_lists[["crf_group"]])
corr_draws_2d  <- rbindlist(draw_lists[["nace2d"]])
corr_draws_5d  <- rbindlist(draw_lists[["nace5d"]])
corr_draws_crf <- rbindlist(draw_lists[["crf_group"]])

cat(sprintf("\n  Point -- 2d: %d | 5d: %d | crf: %d sector-years\n",
            nrow(corr_point_2d), nrow(corr_point_5d), nrow(corr_point_crf)))
cat(sprintf("  Draws -- 2d: %d | 5d: %d | crf: %d\n\n",
            nrow(corr_draws_2d), nrow(corr_draws_5d), nrow(corr_draws_crf)))


# =============================================================================
# SECTION 4: Summarize draws into PIs
# =============================================================================

cat("-- Summarizing draws into prediction intervals ----------------------\n")

stat_cols <- c("spearman_s1_up", "pearson_s1_up")

build_pi <- function(draws_dt, group_col) {
  draws_dt[, {
    out <- list(n_draws = .N)
    for (col in stat_cols) {
      v <- get(col)
      v <- v[!is.na(v)]
      if (length(v) >= 5L) {
        out[[paste0(col, "_median")]] <- median(v)
        out[[paste0(col, "_lo95")]]   <- quantile(v, 0.025, names = FALSE)
        out[[paste0(col, "_hi95")]]   <- quantile(v, 0.975, names = FALSE)
        out[[paste0(col, "_lo90")]]   <- quantile(v, 0.05,  names = FALSE)
        out[[paste0(col, "_hi90")]]   <- quantile(v, 0.95,  names = FALSE)
      } else {
        out[[paste0(col, "_median")]] <- NA_real_
        out[[paste0(col, "_lo95")]]   <- NA_real_
        out[[paste0(col, "_hi95")]]   <- NA_real_
        out[[paste0(col, "_lo90")]]   <- NA_real_
        out[[paste0(col, "_hi90")]]   <- NA_real_
      }
    }
    out
  }, by = c(group_col, "year")]
}

corr_summary_2d  <- merge(build_pi(corr_draws_2d,  "nace2d"),    corr_point_2d,
                           by = c("nace2d", "year"),    all.x = TRUE, suffixes = c("", "_point"))
corr_summary_5d  <- merge(build_pi(corr_draws_5d,  "nace5d"),    corr_point_5d,
                           by = c("nace5d", "year"),    all.x = TRUE, suffixes = c("", "_point"))
corr_summary_crf <- merge(build_pi(corr_draws_crf, "crf_group"), corr_point_crf,
                           by = c("crf_group", "year"), all.x = TRUE, suffixes = c("", "_point"))

cat(sprintf("  Summary -- 2d: %d | 5d: %d | crf: %d sector-years with PIs\n\n",
            nrow(corr_summary_2d), nrow(corr_summary_5d), nrow(corr_summary_crf)))


# =============================================================================
# SECTION 5: Diagnostics
# =============================================================================

cat("-- Diagnostics -------------------------------------------------------\n\n")

cat("Spearman(scope1, upstream) point estimate by NACE 2d (averaged across years):\n")
pt_avg <- corr_point_2d[, .(mean_spearman = round(mean(spearman_s1_up, na.rm = TRUE), 3),
                              mean_pearson  = round(mean(pearson_s1_up, na.rm = TRUE), 3),
                              n_years       = .N),
                          by = nace2d][order(nace2d)]
print(pt_avg)

cat("\nSpearman(scope1, upstream) point estimate by CRF group (averaged across years):\n")
pt_avg_crf <- corr_point_crf[, .(mean_spearman = round(mean(spearman_s1_up, na.rm = TRUE), 3),
                                   mean_pearson  = round(mean(pearson_s1_up, na.rm = TRUE), 3),
                                   n_years       = .N),
                               by = crf_group][order(crf_group)]
print(pt_avg_crf)

cat("\n95% PI width for Spearman (NACE 2d):\n")
corr_summary_2d[, spearman_width := spearman_s1_up_hi95 - spearman_s1_up_lo95]
cat(sprintf("  Median: %.4f\n", median(corr_summary_2d$spearman_width, na.rm = TRUE)))


# =============================================================================
# SECTION 6: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA,
  sprintf("upstream_correlation_pi_%s.RData", WEIGHT_SCHEME))
save(
  corr_point_2d,  corr_draws_2d,  corr_summary_2d,
  corr_point_5d,  corr_draws_5d,  corr_summary_5d,
  corr_point_crf, corr_draws_crf, corr_summary_crf,
  WEIGHT_SCHEME, MIN_N_CORR, stat_cols,
  file = OUT_PATH
)

cat(sprintf("\n===================================================================\n"))
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  2d  -- point: %d | draws: %d | summary: %d\n",
            nrow(corr_point_2d),  nrow(corr_draws_2d),  nrow(corr_summary_2d)))
cat(sprintf("  5d  -- point: %d | draws: %d | summary: %d\n",
            nrow(corr_point_5d),  nrow(corr_draws_5d),  nrow(corr_summary_5d)))
cat(sprintf("  crf -- point: %d | draws: %d | summary: %d\n",
            nrow(corr_point_crf), nrow(corr_draws_crf), nrow(corr_summary_crf)))
cat("===================================================================\n")
