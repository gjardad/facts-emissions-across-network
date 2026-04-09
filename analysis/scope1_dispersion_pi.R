###############################################################################
# analysis/scope1_dispersion_pi.R
#
# PURPOSE
#   RQ1: within-sector dispersion of scope 1 emissions and carbon productivity,
#   with prediction intervals from the B perturbation draws.
#
#   For each sector × year (at three granularities: NACE 2-digit, NACE 5-digit,
#   and CRF category):
#     - Point estimate: dispersion statistics from the deterministic GLO
#       allocation (allocation_glo_<scheme>/).
#     - Prediction interval: 2.5th and 97.5th percentiles of the same
#       statistics across B perturbation draws (uncertainty_draws_<scheme>/).
#
#   Dispersion statistics:
#     s1_gini, s1_p90p10, s1_p75p25, s1_var_log (scope 1 emissions)
#     cp_p90p10, cp_p75p25, cp_var_log, cp_p9010_log (carbon productivity)
#
# INPUT
#   {PROC_DATA}/allocation_glo_<scheme>/alloc_YYYY.RData  (point estimate)
#   {PROC_DATA}/uncertainty_draws_<scheme>/draw_BBBB.RData (B draws)
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   {PROC_DATA}/scope1_dispersion_pi_<scheme>.RData
#     disp_point_*  : point-estimate dispersion by sector x year
#     disp_draws_*  : all B draw-level dispersion stats
#     disp_summary_*: point estimate + 90%/95% PIs per sector-year
#     (* = 2d, 5d, crf for the three granularities)
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
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# -- Parameters ---------------------------------------------------------------
WEIGHT_SCHEME <- "balanced"
YEARS         <- 2005:2021
MIN_N_STATS   <- 3L

ALLOC_DIR <- file.path(PROC_DATA, sprintf("allocation_glo_%s", WEIGHT_SCHEME))
DRAW_DIR  <- file.path(PROC_DATA, sprintf("uncertainty_draws_%s", WEIGHT_SCHEME))

cat("===================================================================\n")
cat("  SCOPE 1 DISPERSION WITH PREDICTION INTERVALS\n")
cat("  WEIGHT_SCHEME =", WEIGHT_SCHEME, "\n")
cat("  Years:", min(YEARS), "--", max(YEARS), "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load accounts (revenue + NACE)
# =============================================================================

cat("-- Loading accounts -------------------------------------------------\n")

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)
setDT(accounts)

# CRF crosswalk: nace2d -> crf_group
nace_crf <- fread(file.path(REPO_DIR, "preprocess", "crosswalks",
                             "nace_crf_crosswalk.csv"),
                   select = c("nace2d", "crf_group"),
                   colClasses = c(nace2d = "character"))
# Add combined "17/18" row (both map to "paper") to match make_nace2d() convention
nace_crf <- rbind(nace_crf,
                   data.table(nace2d = "17/18", crf_group = "paper"))
accounts <- merge(accounts, nace_crf, by = "nace2d", all.x = TRUE)
cat("  Accounts firm-years:", nrow(accounts), "\n\n")


# =============================================================================
# SECTION 2: Dispersion statistics helper functions
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

compute_sector_stats <- function(firms_dt, group_col) {
  # firms_dt must have: scope1, revenue, and the column named by group_col
  # Returns one row per group with dispersion stats
  firms_dt <- firms_dt[scope1 > 0 & !is.na(get(group_col))]

  firms_dt[, .(
    n_firms     = .N,
    s1_gini     = gini(scope1),
    s1_p90p10   = pct_ratio(scope1, 0.9, 0.1),
    s1_p75p25   = pct_ratio(scope1, 0.75, 0.25),
    s1_var_log  = if (.N >= 2L) var(log(scope1 + 1)) else NA_real_,
    cp_p90p10   = {
      ok <- scope1 > 0 & revenue > 0
      if (sum(ok) >= 2L) pct_ratio(revenue[ok] / scope1[ok], 0.9, 0.1) else NA_real_
    },
    cp_p75p25   = {
      ok <- scope1 > 0 & revenue > 0
      if (sum(ok) >= 2L) pct_ratio(revenue[ok] / scope1[ok], 0.75, 0.25) else NA_real_
    },
    cp_var_log  = {
      ok <- scope1 > 0 & revenue > 0
      cp <- revenue[ok] / scope1[ok]
      if (length(cp) >= 2L) var(log(cp)) else NA_real_
    },
    cp_p9010_log = {
      ok <- scope1 > 0 & revenue > 0
      cp <- revenue[ok] / scope1[ok]
      if (length(cp) >= 2L) diff(quantile(log(cp), c(0.1, 0.9))) else NA_real_
    }
  ), by = group_col][n_firms >= MIN_N_STATS]
}


# =============================================================================
# SECTION 3: Point estimate from deterministic allocation
# =============================================================================

cat("-- Point estimate from deterministic allocation ----------------------\n")

GROUP_LEVELS <- c("nace2d", "nace5d", "crf_group")

point_lists <- setNames(lapply(GROUP_LEVELS, function(x) list()), GROUP_LEVELS)

for (t in YEARS) {
  alloc_path <- file.path(ALLOC_DIR, sprintf("alloc_%d.RData", t))
  if (!file.exists(alloc_path)) {
    cat(sprintf("  Year %d: SKIPPED (no allocation file)\n", t))
    next
  }
  load(alloc_path)  # loads year_firms

  setDT(year_firms)
  firms_t <- merge(year_firms[, .(vat, scope1)],
                    accounts[year == t, .(vat, nace2d, nace5d, crf_group, revenue)],
                    by = "vat", all.x = TRUE)

  for (g in GROUP_LEVELS) {
    st <- compute_sector_stats(firms_t, g)
    if (nrow(st) > 0) {
      st[, year := t]
      point_lists[[g]][[length(point_lists[[g]]) + 1]] <- st
    }
  }
}

disp_point_2d  <- rbindlist(point_lists[["nace2d"]])
disp_point_5d  <- rbindlist(point_lists[["nace5d"]])
disp_point_crf <- rbindlist(point_lists[["crf_group"]])
cat(sprintf("  Point estimate — 2d: %d | 5d: %d | crf: %d sector-years\n\n",
            nrow(disp_point_2d), nrow(disp_point_5d), nrow(disp_point_crf)))


# =============================================================================
# SECTION 4: Dispersion across B perturbation draws
# =============================================================================

cat("-- Dispersion across B perturbation draws ---------------------------\n")

draw_files <- sort(list.files(DRAW_DIR, pattern = "^draw_\\d+\\.RData$",
                               full.names = TRUE))
B <- length(draw_files)
cat(sprintf("  Found %d draw files\n", B))

draw_lists <- setNames(lapply(GROUP_LEVELS, function(x) list()), GROUP_LEVELS)
t0 <- Sys.time()

for (b in seq_len(B)) {
  load(draw_files[b])  # loads draw_firms
  setDT(draw_firms)

  # draw_firms has: vat, year, crf_group, scope1_b, source, rank_in_cell, ...
  setnames(draw_firms, "scope1_b", "scope1")

  for (t in YEARS) {
    firms_bt <- draw_firms[year == t]
    if (nrow(firms_bt) == 0) next

    firms_bt <- merge(firms_bt[, .(vat, scope1)],
                       accounts[year == t, .(vat, nace2d, nace5d, crf_group, revenue)],
                       by = "vat", all.x = TRUE)

    for (g in GROUP_LEVELS) {
      st <- compute_sector_stats(firms_bt, g)
      if (nrow(st) > 0) {
        st[, `:=`(year = t, draw = b)]
        draw_lists[[g]][[length(draw_lists[[g]]) + 1]] <- st
      }
    }
  }

  if (b %% 25 == 0) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("  Draw %d / %d (%.1fs elapsed)\n", b, B, elapsed))
  }
}

disp_draws_2d  <- rbindlist(draw_lists[["nace2d"]])
disp_draws_5d  <- rbindlist(draw_lists[["nace5d"]])
disp_draws_crf <- rbindlist(draw_lists[["crf_group"]])
total_time <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat(sprintf("  All draws done in %.1f min\n", total_time))
cat(sprintf("  Draw-level — 2d: %d | 5d: %d | crf: %d\n\n",
            nrow(disp_draws_2d), nrow(disp_draws_5d), nrow(disp_draws_crf)))


# =============================================================================
# SECTION 5: Summarize draws into PIs
# =============================================================================

cat("-- Summarizing draws into prediction intervals ----------------------\n")

stat_cols <- c("s1_gini", "s1_p90p10", "s1_p75p25", "s1_var_log",
               "cp_p90p10", "cp_p75p25", "cp_var_log", "cp_p9010_log")

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

disp_summary_2d  <- merge(build_pi(disp_draws_2d,  "nace2d"),    disp_point_2d,
                           by = c("nace2d", "year"),    all.x = TRUE, suffixes = c("", "_point"))
disp_summary_5d  <- merge(build_pi(disp_draws_5d,  "nace5d"),    disp_point_5d,
                           by = c("nace5d", "year"),    all.x = TRUE, suffixes = c("", "_point"))
disp_summary_crf <- merge(build_pi(disp_draws_crf, "crf_group"), disp_point_crf,
                           by = c("crf_group", "year"), all.x = TRUE, suffixes = c("", "_point"))

cat(sprintf("  Summary — 2d: %d | 5d: %d | crf: %d sector-years with PIs\n\n",
            nrow(disp_summary_2d), nrow(disp_summary_5d), nrow(disp_summary_crf)))


# =============================================================================
# SECTION 6: Diagnostics
# =============================================================================

cat("-- Diagnostics -------------------------------------------------------\n\n")

cat("Point-estimate s1_gini by NACE 2d (averaged across years):\n")
pt_avg <- disp_point_2d[, .(mean_gini = round(mean(s1_gini, na.rm = TRUE), 3),
                              n_years   = .N),
                          by = nace2d][order(nace2d)]
print(pt_avg)

cat("\nPoint-estimate s1_gini by CRF group (averaged across years):\n")
pt_avg_crf <- disp_point_crf[, .(mean_gini = round(mean(s1_gini, na.rm = TRUE), 3),
                                   n_years   = .N),
                               by = crf_group][order(crf_group)]
print(pt_avg_crf)

cat("\n95% PI width for s1_gini (NACE 2d):\n")
disp_summary_2d[, gini_width := s1_gini_hi95 - s1_gini_lo95]
cat(sprintf("  Median: %.4f\n", median(disp_summary_2d$gini_width, na.rm = TRUE)))


# =============================================================================
# SECTION 7: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA,
  sprintf("scope1_dispersion_pi_%s.RData", WEIGHT_SCHEME))
save(
  disp_point_2d,  disp_draws_2d,  disp_summary_2d,
  disp_point_5d,  disp_draws_5d,  disp_summary_5d,
  disp_point_crf, disp_draws_crf, disp_summary_crf,
  WEIGHT_SCHEME, MIN_N_STATS, stat_cols,
  file = OUT_PATH
)

cat(sprintf("\n===================================================================\n"))
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  2d  — point: %d | draws: %d | summary: %d\n",
            nrow(disp_point_2d),  nrow(disp_draws_2d),  nrow(disp_summary_2d)))
cat(sprintf("  5d  — point: %d | draws: %d | summary: %d\n",
            nrow(disp_point_5d),  nrow(disp_draws_5d),  nrow(disp_summary_5d)))
cat(sprintf("  crf — point: %d | draws: %d | summary: %d\n",
            nrow(disp_point_crf), nrow(disp_draws_crf), nrow(disp_summary_crf)))
cat("===================================================================\n")
