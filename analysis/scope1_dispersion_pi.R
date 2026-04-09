###############################################################################
# analysis/scope1_dispersion_pi.R
#
# PURPOSE
#   RQ1: within-sector dispersion of scope 1 emissions and carbon productivity,
#   with prediction intervals from the B perturbation draws.
#
#   For each NACE 2-digit sector × year:
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
#     disp_point    : data.frame with point-estimate dispersion by nace2d x year
#     disp_draws    : data.frame with all B draw-level dispersion stats
#     disp_summary  : data.frame with point estimate + 90%/95% PIs per cell
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

compute_sector_stats <- function(firms_dt) {
  # firms_dt must have: nace2d, scope1, revenue (all with year already fixed)
  # Returns one row per nace2d with dispersion stats
  firms_dt <- firms_dt[scope1 > 0 & !is.na(nace2d)]

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
  ), by = nace2d][n_firms >= MIN_N_STATS]
}


# =============================================================================
# SECTION 3: Point estimate from deterministic allocation
# =============================================================================

cat("-- Point estimate from deterministic allocation ----------------------\n")

disp_point_list <- list()

for (t in YEARS) {
  alloc_path <- file.path(ALLOC_DIR, sprintf("alloc_%d.RData", t))
  if (!file.exists(alloc_path)) {
    cat(sprintf("  Year %d: SKIPPED (no allocation file)\n", t))
    next
  }
  load(alloc_path)  # loads year_firms

  setDT(year_firms)
  firms_t <- merge(year_firms[, .(vat, scope1)],
                    accounts[year == t, .(vat, nace2d, revenue)],
                    by = "vat", all.x = TRUE)

  st <- compute_sector_stats(firms_t)
  if (nrow(st) > 0) {
    st[, year := t]
    disp_point_list[[length(disp_point_list) + 1]] <- st
  }
}

disp_point <- rbindlist(disp_point_list)
cat(sprintf("  Point estimate: %d sector-years\n\n", nrow(disp_point)))


# =============================================================================
# SECTION 4: Dispersion across B perturbation draws
# =============================================================================

cat("-- Dispersion across B perturbation draws ---------------------------\n")

draw_files <- sort(list.files(DRAW_DIR, pattern = "^draw_\\d+\\.RData$",
                               full.names = TRUE))
B <- length(draw_files)
cat(sprintf("  Found %d draw files\n", B))

disp_draw_list <- list()
t0 <- Sys.time()

for (b in seq_len(B)) {
  load(draw_files[b])  # loads draw_firms
  setDT(draw_firms)

  # Join accounts for revenue and nace2d
  # draw_firms has: vat, year, crf_group, scope1_b, source, rank_in_cell, ...
  # Rename scope1_b -> scope1 for the stats function
  setnames(draw_firms, "scope1_b", "scope1")

  for (t in YEARS) {
    firms_bt <- draw_firms[year == t]
    if (nrow(firms_bt) == 0) next

    firms_bt <- merge(firms_bt[, .(vat, scope1)],
                       accounts[year == t, .(vat, nace2d, revenue)],
                       by = "vat", all.x = TRUE)

    st <- compute_sector_stats(firms_bt)
    if (nrow(st) > 0) {
      st[, `:=`(year = t, draw = b)]
      disp_draw_list[[length(disp_draw_list) + 1]] <- st
    }
  }

  if (b %% 25 == 0) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("  Draw %d / %d (%.1fs elapsed)\n", b, B, elapsed))
  }
}

disp_draws <- rbindlist(disp_draw_list)
total_time <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat(sprintf("  All draws done in %.1f min\n", total_time))
cat(sprintf("  Total draw-level sector-years: %d\n\n", nrow(disp_draws)))


# =============================================================================
# SECTION 5: Summarize draws into PIs
# =============================================================================

cat("-- Summarizing draws into prediction intervals ----------------------\n")

stat_cols <- c("s1_gini", "s1_p90p10", "s1_p75p25", "s1_var_log",
               "cp_p90p10", "cp_p75p25", "cp_var_log", "cp_p9010_log")

# Per (nace2d, year): median, 2.5th, 97.5th, 5th, 95th across draws
disp_pi <- disp_draws[, {
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
}, by = .(nace2d, year)]

# Merge point estimate into the PI summary
disp_summary <- merge(disp_pi, disp_point,
                       by = c("nace2d", "year"), all.x = TRUE,
                       suffixes = c("", "_point"))

cat(sprintf("  Summary: %d sector-years with PIs\n\n", nrow(disp_summary)))


# =============================================================================
# SECTION 6: Diagnostics
# =============================================================================

cat("-- Diagnostics -------------------------------------------------------\n\n")

cat("Point-estimate s1_gini by sector (averaged across years):\n")
pt_avg <- disp_point[, .(mean_gini = round(mean(s1_gini, na.rm = TRUE), 3),
                          n_years   = .N),
                      by = nace2d][order(nace2d)]
print(pt_avg)

cat("\nPI width (95%) for s1_gini, averaged across sector-years:\n")
disp_summary[, gini_width := s1_gini_hi95 - s1_gini_lo95]
cat(sprintf("  Median PI width: %.4f\n", median(disp_summary$gini_width, na.rm = TRUE)))
cat(sprintf("  Mean PI width:   %.4f\n", mean(disp_summary$gini_width, na.rm = TRUE)))

cat("\nExample sector-years (first 10 rows of disp_summary):\n")
print(head(disp_summary[order(nace2d, year),
  .(nace2d, year, n_draws,
    s1_gini, s1_gini_median, s1_gini_lo95, s1_gini_hi95)], 10))


# =============================================================================
# SECTION 7: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA,
  sprintf("scope1_dispersion_pi_%s.RData", WEIGHT_SCHEME))
save(disp_point, disp_draws, disp_summary,
     WEIGHT_SCHEME, MIN_N_STATS, stat_cols,
     file = OUT_PATH)

cat(sprintf("\n===================================================================\n"))
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  disp_point:   %d sector-years (deterministic)\n", nrow(disp_point)))
cat(sprintf("  disp_draws:   %d draw x sector-years\n", nrow(disp_draws)))
cat(sprintf("  disp_summary: %d sector-years with PIs\n", nrow(disp_summary)))
cat("===================================================================\n")
