###############################################################################
# analysis/report_uncertainty.R
#
# PURPOSE
#   Aggregate the B perturbation draws produced by uncertainty_propagation.R
#   into per firm-year point estimates and 95% prediction intervals.
#
#   For each firm-year (vat, year), we have B simulated emissions values
#   ytilde_i1, ..., ytilde_iB. We report:
#     - point_mean   : mean across draws
#     - point_median : median across draws
#     - pi_lo_95     : 2.5th percentile across draws
#     - pi_hi_95     : 97.5th percentile across draws
#     - pi_lo_90     : 5th percentile across draws
#     - pi_hi_90     : 95th percentile across draws
#     - p_zero       : fraction of draws with ytilde = 0
#     - source       : "ets" / "pre_ets" / "imputed" (modal across draws)
#
# INPUT
#   {PROC_DATA}/uncertainty_draws/draw_BBBB.RData  (one file per draw)
#
# OUTPUT
#   {PROC_DATA}/firm_year_pi.RData
#     firm_year_pi : data.frame with one row per (vat, year)
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

# -- Parameter ----------------------------------------------------------------
WEIGHT_SCHEME <- "balanced"

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

cat("===================================================================\n")
cat("  REPORT UNCERTAINTY: build firm-year prediction intervals\n")
cat("  WEIGHT_SCHEME =", WEIGHT_SCHEME, "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: List draw files
# =============================================================================

DRAW_DIR <- file.path(PROC_DATA, sprintf("uncertainty_draws_%s", WEIGHT_SCHEME))
draw_files <- sort(list.files(DRAW_DIR, pattern = "^draw_\\d+\\.RData$",
                               full.names = TRUE))
B <- length(draw_files)
cat("Found", B, "draw files in", DRAW_DIR, "\n\n")

if (B == 0) stop("No draw files found.")


# =============================================================================
# SECTION 2: Stack draws into a long data.frame
# =============================================================================

cat("-- Loading draws ----------------------------------------------------\n")

t0 <- Sys.time()
all_draws <- vector("list", B)
for (b in seq_len(B)) {
  load(draw_files[b])  # loads draw_firms
  all_draws[[b]] <- draw_firms[, c("vat", "year", "crf_group", "scope1_b",
                                    "source", "draw")]
  if (b %% 25 == 0) cat(sprintf("  Loaded draw %d / %d\n", b, B))
}
all_long <- bind_rows(all_draws)
rm(all_draws); gc()

cat(sprintf("\nStacked %d rows in %.1fs\n\n",
            nrow(all_long),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))


# =============================================================================
# SECTION 3: Per firm-year aggregation (data.table for speed)
# =============================================================================

cat("-- Aggregating to firm-year level ----------------------------------\n")

# A firm-year may not appear in every draw (drawn as non-emitter → no row).
# Missing rows count as scope1_b = 0. We exploit this by NOT padding zeros
# explicitly and instead computing the summaries from the n_appearances and
# the sorted non-zero values.
#
# Conventions used below:
#   B          : total number of draws
#   n_apps     : number of draws this firm-year appears in (with positive
#                allocation, i.e. drawn as an emitter)
#   B - n_apps : number of "missing" draws → scope1_b = 0 in those draws
#
# For ETS / pre-ETS firms: n_apps == B and scope1_b is constant across draws.
# For imputed firms: n_apps in [0, B] and scope1_b varies.

setDT(all_long)

# Quantile of a vector v with extra (B - n) zeros prepended.
# Equivalent to quantile(c(rep(0, B-n), sort(v)), q, type = 7).
quantile_with_zeros <- function(v_sorted, n, B, q) {
  # Position in [1, B], type 7: h = (B - 1) * q + 1
  h <- (B - 1) * q + 1
  # Number of leading zeros
  n_zero <- B - n
  if (h <= n_zero + 1e-12) return(0)
  # Index into v_sorted (1-based)
  j <- h - n_zero
  j_lo <- floor(j); j_hi <- ceiling(j)
  if (j_lo < 1) j_lo <- 1
  if (j_hi > n) j_hi <- n
  if (j_lo == j_hi) return(v_sorted[j_lo])
  v_sorted[j_lo] + (j - j_lo) * (v_sorted[j_hi] - v_sorted[j_lo])
}

cat("  Computing summaries with data.table...\n")
t1 <- Sys.time()

# Group by (vat, year, crf_group) and compute summaries
# Note: source modal is computed via first() since within a (vat, year) pair
# the source is constant across draws (a firm-year is always either ets,
# pre_ets, or imputed — never mixed across draws).

firm_year_dt <- all_long[, {
  v_sorted <- sort(scope1_b)
  n        <- .N
  list(
    crf_group     = crf_group[1],
    source        = source[1],
    n_appearances = n,
    point_mean    = sum(scope1_b) / B,
    point_median  = quantile_with_zeros(v_sorted, n, B, 0.5),
    pi_lo_95      = quantile_with_zeros(v_sorted, n, B, 0.025),
    pi_hi_95      = quantile_with_zeros(v_sorted, n, B, 0.975),
    pi_lo_90      = quantile_with_zeros(v_sorted, n, B, 0.05),
    pi_hi_90      = quantile_with_zeros(v_sorted, n, B, 0.95),
    p_zero        = (B - n) / B
  )
}, by = .(vat, year)]

cat(sprintf("  Done in %.1fs\n",
            as.numeric(difftime(Sys.time(), t1, units = "secs"))))

firm_year_pi <- as.data.frame(firm_year_dt)


# =============================================================================
# SECTION 4: Diagnostics
# =============================================================================

cat("\n-- Diagnostics ------------------------------------------------------\n\n")

cat("Source breakdown (modal across draws):\n")
print(table(firm_year_pi$source))

cat("\nETS firms (should have p_zero = 0, PI width = 0):\n")
ets_pi <- firm_year_pi[firm_year_pi$source == "ets", ]
cat(sprintf("  n = %d, mean p_zero = %.4f\n", nrow(ets_pi), mean(ets_pi$p_zero)))
cat(sprintf("  mean PI width (95%%) = %.2e (should be 0)\n",
            mean(ets_pi$pi_hi_95 - ets_pi$pi_lo_95)))

cat("\nImputed firms PI summary:\n")
imp_pi <- firm_year_pi[firm_year_pi$source == "imputed", ]
cat(sprintf("  n = %d\n", nrow(imp_pi)))
cat(sprintf("  point_mean: median = %.0f, max = %.2e\n",
            median(imp_pi$point_mean), max(imp_pi$point_mean)))
cat(sprintf("  p_zero distribution:\n"))
print(summary(imp_pi$p_zero))
cat(sprintf("  PI width (95%%): median = %.0f\n",
            median(imp_pi$pi_hi_95 - imp_pi$pi_lo_95)))


# =============================================================================
# SECTION 5: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA,
  sprintf("firm_year_pi_%s.RData", WEIGHT_SCHEME))
save(firm_year_pi, B, WEIGHT_SCHEME, file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  %d firm-years, B = %d draws\n", nrow(firm_year_pi), B))
cat("===================================================================\n")
