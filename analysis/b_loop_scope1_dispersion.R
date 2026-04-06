###############################################################################
# analysis/b_loop_scope1_dispersion_v2.R
#
# PURPOSE
#   RQ1: Carbon productivity dispersion using Pareto-allocated emissions.
#
#   Loads firm-level scope 1 allocation from b_allocation_pareto.R output,
#   joins revenue and NACE codes, and computes within-sector dispersion
#   statistics across B subsample draws.
#
# INPUT
#   {PROC_DATA}/allocation_pareto/alloc_YYYY.RData  (from b_allocation_pareto.R)
#     year_firms_by_draw : list of B data.frames (vat, scope1, euets)
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   {PROC_DATA}/b_loop_scope1_dispersion_pareto.RData
#     stats2d_summary : point estimates + cross-B s.d. by NACE 2-digit x year
#     stats5d_summary : same at NACE 5-digit x year
#     all_stats2d     : full B x sector-year draw-level data
#     all_stats5d     : same at NACE 5-digit
#
# RUNS ON: RMD or local (lightweight — no allocation, no B2B)
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

# -- Parameters ---------------------------------------------------------------
YEARS        <- 2005:2021
MIN_N_STATS  <- 3L
ALLOC_METHOD <- "proportional"  # "pareto" or "proportional"

ALLOC_DIR <- file.path(PROC_DATA, paste0("allocation_", ALLOC_METHOD))

cat("===================================================================\n")
cat("  RQ1: SCOPE 1 CARBON PRODUCTIVITY DISPERSION\n")
cat("  Allocation method:", ALLOC_METHOD, "\n")
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
cat("  Accounts firm-years:", nrow(accounts), "\n\n")


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
# SECTION 3: Main loop -- year x draw
# =============================================================================
cat("-- Main loop: years x draws ----------------------------------------\n\n")

stats2d_all <- list()
stats5d_all <- list()
t0_total    <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()

  alloc_path <- file.path(ALLOC_DIR, sprintf("alloc_%d.RData", t))
  if (!file.exists(alloc_path)) {
    cat(sprintf("Year %d -- SKIPPED (allocation file not found)\n", t))
    next
  }
  load(alloc_path)  # year_firms_by_draw, year_flags
  B <- length(year_firms_by_draw)

  accounts_t <- accounts[accounts$year == t, ]

  for (b in seq_len(B)) {
    firms_bt <- year_firms_by_draw[[b]] %>%
      left_join(accounts_t %>% select(vat, nace2d, nace5d, revenue), by = "vat") %>%
      filter(!is.na(nace5d), scope1 > 0)

    st <- compute_stats_bt(firms_bt)

    stats2d_all[[length(stats2d_all) + 1L]] <- st$stats2d %>%
      mutate(year = t, draw = b)
    stats5d_all[[length(stats5d_all) + 1L]] <- st$stats5d %>%
      mutate(year = t, draw = b)
  }

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  cat(sprintf("Year %d (%s s)\n", t, elapsed))
  rm(year_firms_by_draw, year_flags)
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 4: Aggregate across B draws
# =============================================================================
cat("-- Aggregating across draws -----------------------------------------\n")

all_stats2d <- bind_rows(stats2d_all)
all_stats5d <- bind_rows(stats5d_all)

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
cat("  NACE 5-digit:", nrow(stats5d_summary), "sector-years\n\n")


# =============================================================================
# SECTION 5: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, paste0("b_loop_scope1_dispersion_", ALLOC_METHOD, ".RData"))

save(
  stats2d_summary,
  stats5d_summary,
  all_stats2d,
  all_stats5d,
  file = OUT_PATH
)

cat("===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat("  stats2d_summary:", nrow(stats2d_summary), "rows\n")
cat("  stats5d_summary:", nrow(stats5d_summary), "rows\n")
cat("  Total time:", total_time, "min\n")
cat("===================================================================\n")
