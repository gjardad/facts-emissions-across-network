###############################################################################
# analysis/b_loop_upstream_dispersion.R
#
# PURPOSE
#   RQ2: How dispersed are upstream network-adjusted emissions across firms
#   within sectors? How does this compare with scope 1 dispersion?
#
#   For each year and each of B subsample draws:
#     1. Load firm-level (scope1, upstream, decomposition) from
#        build_upstream_emissions.R output
#     2. Join revenue and NACE codes from annual accounts
#     3. Compute within-sector dispersion of upstream carbon productivity
#        (revenue / upstream), comparable to RQ1's scope 1 carbon productivity
#     4. Decompose variance of upstream emission intensity into:
#        - across-sector component (input mix × sector-avg supplier intensity)
#        - within-sector component (supplier selection within sectors)
#        at both NACE 2-digit and 5-digit granularity
#
#   All statistics are revenue-based (per EUR of revenue) for comparability
#   with RQ1, even though the underlying Leontief is cost-based.
#
# INPUT
#   {PROC_DATA}/upstream_emissions/firms_YYYY.RData  (from build_upstream_emissions.R)
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   {PROC_DATA}/b_loop_upstream_dispersion.RData
#     stats2d_summary : point estimates + cross-B s.d. by NACE 2-digit × year
#     stats5d_summary : same at NACE 5-digit × year
#     all_stats2d     : full B × sector-year draw-level data
#     all_stats5d     : same at NACE 5-digit
#
# RUNS ON: RMD (or local with downsampled data)
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

library(dplyr)

# ── Parameters ────────────────────────────────────────────────────────────────
YEARS       <- 2005:2021
MIN_N_STATS <- 3L

UP_DIR <- file.path(PROC_DATA, "upstream_emissions")

cat("═══════════════════════════════════════════════════════════════\n")
cat("  RQ2: UPSTREAM CARBON PRODUCTIVITY DISPERSION + DECOMPOSITION\n")
cat("  Years:", min(YEARS), "–", max(YEARS), "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")


# =============================================================================
# SECTION 1: Load accounts (revenue + NACE)
# =============================================================================
cat("── Loading accounts ─────────────────────────────────────────\n")

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = substr(nace5d, 1, 2),
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

safe_spearman <- function(x, y, min_n = MIN_N_STATS) {
  keep <- !is.na(x) & !is.na(y) & x > 0 & y > 0
  if (sum(keep) < min_n) return(NA_real_)
  cor(x[keep], y[keep], method = "spearman")
}

# Variance decomposition of upstream emission intensity into across-sector
# and within-sector components. Returns shares of total variance.
var_decomp <- function(upstream, across, within, revenue, min_n = MIN_N_STATS) {
  ok <- !is.na(upstream) & !is.na(revenue) & upstream > 0 & revenue > 0
  if (sum(ok) < min_n) return(list(var_total = NA_real_, share_across = NA_real_,
                                    share_within = NA_real_, share_cov = NA_real_))
  u <- upstream[ok] / revenue[ok]
  a <- across[ok]  / revenue[ok]
  w <- within[ok]  / revenue[ok]

  var_total  <- var(u)
  var_across <- var(a)
  var_within <- var(w)
  cov_aw     <- cov(a, w)

  # var(u) = var(a) + var(w) + 2*cov(a,w) by construction
  list(
    var_total    = var_total,
    share_across = var_across / var_total,
    share_within = var_within / var_total,
    share_cov    = 2 * cov_aw / var_total
  )
}

compute_stats_bt <- function(firms_bt) {
  # firms_bt: (vat, nace2d, nace5d, scope1, upstream, upstream_across_2d,
  #            upstream_within_2d, upstream_across_5d, upstream_within_5d,
  #            revenue, euets)

  do_stats <- function(df, decomp_suffix) {
    # Which decomposition columns to use (2d or 5d)
    across_col <- df[[paste0("upstream_across_", decomp_suffix)]]
    within_col <- df[[paste0("upstream_within_", decomp_suffix)]]

    # Upstream carbon productivity: revenue / upstream (firms with both > 0)
    ok_up   <- df$upstream > 0 & df$revenue > 0
    up_cp   <- df$revenue[ok_up] / df$upstream[ok_up]
    n_up_cp <- length(up_cp)
    log_up  <- if (n_up_cp >= 2L) log(up_cp) else NA_real_

    # Scope 1 carbon productivity (for direct comparison)
    ok_s1    <- df$scope1 > 0 & df$revenue > 0
    s1_cp    <- df$revenue[ok_s1] / df$scope1[ok_s1]
    n_s1_cp  <- length(s1_cp)
    log_s1   <- if (n_s1_cp >= 2L) log(s1_cp) else NA_real_

    # Variance decomposition
    vd <- var_decomp(df$upstream, across_col, within_col, df$revenue)

    data.frame(
      n_firms       = nrow(df),
      # Upstream carbon productivity dispersion
      n_up_cp       = n_up_cp,
      up_cp_p90p10    = if (n_up_cp >= 2L) pct_ratio(up_cp, 0.9, 0.1) else NA_real_,
      up_cp_p75p25    = if (n_up_cp >= 2L) pct_ratio(up_cp, 0.75, 0.25) else NA_real_,
      up_cp_var_log   = if (n_up_cp >= 2L) var(log_up) else NA_real_,
      up_cp_p9010_log = if (n_up_cp >= 2L) diff(quantile(log_up, c(0.1, 0.9))) else NA_real_,
      up_cp_gini      = if (n_up_cp >= 2L) gini(up_cp) else NA_real_,
      # Scope 1 carbon productivity (for same-table comparison)
      n_s1_cp       = n_s1_cp,
      s1_cp_p90p10    = if (n_s1_cp >= 2L) pct_ratio(s1_cp, 0.9, 0.1) else NA_real_,
      s1_cp_p75p25    = if (n_s1_cp >= 2L) pct_ratio(s1_cp, 0.75, 0.25) else NA_real_,
      s1_cp_var_log   = if (n_s1_cp >= 2L) var(log_s1) else NA_real_,
      s1_cp_p9010_log = if (n_s1_cp >= 2L) diff(quantile(log_s1, c(0.1, 0.9))) else NA_real_,
      s1_cp_gini      = if (n_s1_cp >= 2L) gini(s1_cp) else NA_real_,
      # Rank correlation: scope 1 vs upstream
      spearman_s1_up  = safe_spearman(df$scope1, df$upstream),
      spearman_s1_up_ets = safe_spearman(df$scope1[df$euets == 1L],
                                          df$upstream[df$euets == 1L]),
      # Variance decomposition
      vd_var_total    = vd$var_total,
      vd_share_across = vd$share_across,
      vd_share_within = vd$share_within,
      vd_share_cov    = vd$share_cov
    )
  }

  # NACE 2d stats: use 2d decomposition
  stats2d <- firms_bt %>%
    filter(scope1 > 0 | upstream > 0) %>%
    group_by(nace2d) %>%
    filter(n() >= MIN_N_STATS) %>%
    group_modify(~do_stats(.x, "2d")) %>%
    ungroup()

  # NACE 5d stats: use 5d decomposition
  stats5d <- firms_bt %>%
    filter(scope1 > 0 | upstream > 0) %>%
    group_by(nace5d) %>%
    filter(n() >= MIN_N_STATS) %>%
    group_modify(~do_stats(.x, "5d")) %>%
    ungroup()

  list(stats2d = stats2d, stats5d = stats5d)
}


# =============================================================================
# SECTION 3: Main loop — year × draw
# =============================================================================
cat("── Main loop: years × draws ─────────────────────────────────\n\n")

stats2d_all <- list()
stats5d_all <- list()
t0_total    <- Sys.time()

for (t in YEARS) {
  t0_year <- Sys.time()

  # Load pre-computed upstream emissions for year t
  up_path <- file.path(UP_DIR, sprintf("firms_%d.RData", t))
  if (!file.exists(up_path)) {
    cat(sprintf("Year %d — SKIPPED (file not found)\n", t))
    next
  }
  load(up_path)  # firms_by_draw, conv_info, max_rowsum
  B <- length(firms_by_draw)

  accounts_t <- accounts[accounts$year == t, ]

  for (b in seq_len(B)) {
    firms_bt <- firms_by_draw[[b]] %>%
      left_join(accounts_t %>% select(vat, nace2d, nace5d, revenue), by = "vat") %>%
      filter(!is.na(nace5d))

    st <- compute_stats_bt(firms_bt)

    stats2d_all[[length(stats2d_all) + 1L]] <- st$stats2d %>%
      mutate(year = t, draw = b)
    stats5d_all[[length(stats5d_all) + 1L]] <- st$stats5d %>%
      mutate(year = t, draw = b)
  }

  elapsed <- round(difftime(Sys.time(), t0_year, units = "secs"), 1)
  cat(sprintf("Year %d (%s s)\n", t, elapsed))
  rm(firms_by_draw, conv_info, max_rowsum)
  gc()
}

total_time <- round(difftime(Sys.time(), t0_total, units = "mins"), 1)
cat(sprintf("\nAll years complete in %.1f min\n\n", total_time))


# =============================================================================
# SECTION 4: Aggregate across B draws
# =============================================================================
cat("── Aggregating across draws ─────────────────────────────────\n")

all_stats2d <- bind_rows(stats2d_all)
all_stats5d <- bind_rows(stats5d_all)

stat_cols <- c(
  # Upstream carbon productivity
  "up_cp_p90p10", "up_cp_p75p25", "up_cp_var_log", "up_cp_p9010_log", "up_cp_gini",
  # Scope 1 carbon productivity (for comparison)
  "s1_cp_p90p10", "s1_cp_p75p25", "s1_cp_var_log", "s1_cp_p9010_log", "s1_cp_gini",
  # Rank correlations
  "spearman_s1_up", "spearman_s1_up_ets",
  # Variance decomposition
  "vd_var_total", "vd_share_across", "vd_share_within", "vd_share_cov"
)

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
OUT_PATH <- file.path(PROC_DATA, "b_loop_upstream_dispersion.RData")

save(
  stats2d_summary,
  stats5d_summary,
  all_stats2d,
  all_stats5d,
  file = OUT_PATH
)

cat("══════════════════════════════════════════════════════════════\n")
cat("Saved:", OUT_PATH, "\n")
cat("  stats2d_summary:", nrow(stats2d_summary), "rows\n")
cat("  stats5d_summary:", nrow(stats5d_summary), "rows\n")
cat("  Total time:", total_time, "min\n")
cat("══════════════════════════════════════════════════════════════\n")
