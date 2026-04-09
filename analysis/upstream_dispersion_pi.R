###############################################################################
# analysis/upstream_dispersion_pi.R
#
# PURPOSE
#   RQ2: within-sector dispersion of upstream embodied emissions and upstream
#   carbon productivity, with prediction intervals from the B perturbation
#   draws plus the deterministic point estimate.
#
#   For each NACE 2-digit sector × year:
#     - Point estimate: statistics from draw 0 (deterministic GLO allocation
#       propagated through the Leontief model)
#     - Prediction interval: 2.5th and 97.5th percentiles across B draws
#
#   Also computes: scope 1 vs upstream correlation (Spearman), upstream
#   variance decomposition (across-sector vs within-sector components).
#
# INPUT
#   {PROC_DATA}/upstream_emissions_glo_<scheme>/firms_YYYY.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   {PROC_DATA}/upstream_dispersion_pi_<scheme>.RData
#     disp_point    : point-estimate dispersion by nace2d x year (from draw 0)
#     disp_draws    : all B draw-level dispersion stats
#     disp_summary  : point estimate + 90%/95% PIs per sector-year
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
MIN_N_STATS   <- 3L

UPSTREAM_DIR <- file.path(PROC_DATA,
  sprintf("upstream_emissions_glo_%s", WEIGHT_SCHEME))

cat("===================================================================\n")
cat("  UPSTREAM DISPERSION WITH PREDICTION INTERVALS\n")
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

compute_upstream_stats <- function(firms_dt) {
  # firms_dt must have: nace2d, scope1, upstream, upstream_within_2d, revenue
  firms_dt <- firms_dt[!is.na(nace2d)]

  firms_dt[upstream > 0 | scope1 > 0, .(
    n_firms            = .N,
    # Upstream dispersion
    up_gini            = gini(upstream),
    up_p90p10          = pct_ratio(upstream, 0.9, 0.1),
    up_p75p25          = pct_ratio(upstream, 0.75, 0.25),
    up_var_log         = if (sum(upstream > 0) >= 2L) var(log(upstream[upstream > 0] + 1)) else NA_real_,
    # Upstream carbon productivity (revenue / upstream)
    ucp_p90p10         = {
      ok <- upstream > 0 & revenue > 0
      if (sum(ok) >= 2L) pct_ratio(revenue[ok] / upstream[ok], 0.9, 0.1) else NA_real_
    },
    ucp_p75p25         = {
      ok <- upstream > 0 & revenue > 0
      if (sum(ok) >= 2L) pct_ratio(revenue[ok] / upstream[ok], 0.75, 0.25) else NA_real_
    },
    ucp_var_log        = {
      ok <- upstream > 0 & revenue > 0
      cp <- revenue[ok] / upstream[ok]
      if (length(cp) >= 2L) var(log(cp)) else NA_real_
    },
    ucp_p9010_log      = {
      ok <- upstream > 0 & revenue > 0
      cp <- revenue[ok] / upstream[ok]
      if (length(cp) >= 2L) diff(quantile(log(cp), c(0.1, 0.9))) else NA_real_
    },
    # Scope 1 vs upstream rank correlation
    spearman_s1_up     = {
      ok <- scope1 > 0 & upstream > 0
      if (sum(ok) >= 3L) cor(scope1[ok], upstream[ok], method = "spearman") else NA_real_
    },
    # Within-sector share of upstream (upstream_within_2d / upstream total)
    share_within_2d    = {
      ok <- upstream > 0
      if (sum(ok) >= 2L) mean(upstream_within_2d[ok] / upstream[ok]) else NA_real_
    }
  ), by = nace2d][n_firms >= MIN_N_STATS]
}


# =============================================================================
# SECTION 3: Process upstream files — point estimate (draw 0) + perturbation
# =============================================================================

cat("-- Processing upstream files ----------------------------------------\n")

disp_point_list <- list()
disp_draw_list  <- list()
t0 <- Sys.time()

for (t in YEARS) {
  up_path <- file.path(UPSTREAM_DIR, sprintf("firms_%d.RData", t))
  if (!file.exists(up_path)) {
    cat(sprintf("  Year %d: SKIPPED (no upstream file)\n", t))
    next
  }
  load(up_path)  # loads firms_by_draw, conv_info, max_rowsum

  accounts_t <- accounts[year == t]
  n_draws <- length(firms_by_draw)

  for (b in seq_len(n_draws)) {
    firms_b <- as.data.table(firms_by_draw[[b]])
    firms_b <- merge(firms_b,
                      accounts_t[, .(vat, nace2d, revenue)],
                      by = "vat", all.x = TRUE)

    st <- compute_upstream_stats(firms_b)
    if (nrow(st) == 0) next

    st[, year := t]

    if (b == 1L) {
      # Draw 0 = deterministic point estimate
      disp_point_list[[length(disp_point_list) + 1]] <- st
    } else {
      st[, draw := b - 1L]
      disp_draw_list[[length(disp_draw_list) + 1]] <- st
    }
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("  Year %d: %d draws (%.1fs elapsed)\n", t, n_draws, elapsed))

  rm(firms_by_draw, conv_info)
  gc()
}

disp_point <- rbindlist(disp_point_list)
disp_draws <- rbindlist(disp_draw_list)
cat(sprintf("\n  Point estimate: %d sector-years\n", nrow(disp_point)))
cat(sprintf("  Draw-level: %d draw x sector-years\n\n", nrow(disp_draws)))


# =============================================================================
# SECTION 4: Summarize draws into PIs
# =============================================================================

cat("-- Summarizing draws into prediction intervals ----------------------\n")

stat_cols <- c("up_gini", "up_p90p10", "up_p75p25", "up_var_log",
               "ucp_p90p10", "ucp_p75p25", "ucp_var_log", "ucp_p9010_log",
               "spearman_s1_up", "share_within_2d")

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

disp_summary <- merge(disp_pi, disp_point,
                       by = c("nace2d", "year"), all.x = TRUE,
                       suffixes = c("", "_point"))

cat(sprintf("  Summary: %d sector-years with PIs\n\n", nrow(disp_summary)))


# =============================================================================
# SECTION 5: Diagnostics
# =============================================================================

cat("-- Diagnostics -------------------------------------------------------\n\n")

cat("Point-estimate upstream Gini by sector (averaged across years):\n")
pt_avg <- disp_point[, .(mean_gini = round(mean(up_gini, na.rm = TRUE), 3),
                          n_years   = .N),
                      by = nace2d][order(nace2d)]
print(pt_avg)

cat("\nUpstream carbon productivity p90/p10 (point estimate):\n")
cat(sprintf("  ucp_p90p10 median (across sector-years): %.1f\n",
            median(disp_point$ucp_p90p10, na.rm = TRUE)))
cat(sprintf("  ucp_p9010_log median: %.3f\n",
            median(disp_point$ucp_p9010_log, na.rm = TRUE)))

cat("\nPI width (95%) for up_gini:\n")
disp_summary[, gini_width := up_gini_hi95 - up_gini_lo95]
cat(sprintf("  Median: %.4f\n", median(disp_summary$gini_width, na.rm = TRUE)))

cat("\nSpearman(scope1, upstream) point estimate:\n")
cat(sprintf("  Median: %.3f\n",
            median(disp_point$spearman_s1_up, na.rm = TRUE)))

cat("\nWithin-sector share of upstream (point estimate):\n")
cat(sprintf("  Median: %.3f\n",
            median(disp_point$share_within_2d, na.rm = TRUE)))


# =============================================================================
# SECTION 6: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA,
  sprintf("upstream_dispersion_pi_%s.RData", WEIGHT_SCHEME))
save(disp_point, disp_draws, disp_summary,
     WEIGHT_SCHEME, MIN_N_STATS, stat_cols,
     file = OUT_PATH)

cat(sprintf("\n===================================================================\n"))
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  disp_point:   %d sector-years\n", nrow(disp_point)))
cat(sprintf("  disp_draws:   %d draw x sector-years\n", nrow(disp_draws)))
cat(sprintf("  disp_summary: %d sector-years with PIs\n", nrow(disp_summary)))
cat("===================================================================\n")
