###############################################################################
# analysis/upstream_variance_decomp_pi.R
#
# PURPOSE
#   Decompose within-sector variance of upstream emission intensity
#   (upstream / revenue) into:
#     - across-sector component: heterogeneity in average supplier-sector
#       emission intensity (input mix × sector-avg multiplier)
#     - within-sector component: firm-level supplier selection within sectors
#     - covariance term: 2 × cov(across, within) / var(total)
#
#   Identity: var(u) = var(a) + var(w) + 2·cov(a,w)
#             1      = share_across + share_within + share_cov
#
#   Point estimates from draw 0 (deterministic GLO allocation); prediction
#   intervals from quantiles across B perturbation draws.
#
#   Decomposition at NACE 2-digit, 5-digit, and CRF category granularity.
#
# INPUT
#   {PROC_DATA}/upstream_emissions_glo_<scheme>/firms_YYYY.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   {PROC_DATA}/upstream_variance_decomp_pi_<scheme>.RData
#     vd_point_*  : point-estimate decomposition by sector x year
#     vd_draws_*  : all B draw-level decomposition stats
#     vd_summary_*: point estimate + 90%/95% PIs per sector-year
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
MIN_N_STATS   <- 3L

UPSTREAM_DIR <- file.path(PROC_DATA,
  sprintf("upstream_emissions_glo_%s", WEIGHT_SCHEME))

cat("===================================================================\n")
cat("  UPSTREAM VARIANCE DECOMPOSITION WITH PREDICTION INTERVALS\n")
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
nace_crf <- rbind(nace_crf,
                   data.table(nace2d = "17/18", crf_group = "paper"))
accounts <- merge(accounts, nace_crf, by = "nace2d", all.x = TRUE)
cat("  Accounts firm-years:", nrow(accounts), "\n\n")


# =============================================================================
# SECTION 2: Variance decomposition helper
# =============================================================================

# Decompose var(upstream/revenue) into across-sector and within-sector shares.
# upstream = across + within by construction (from build_upstream_emissions_glo.R),
# so var(u) = var(a) + var(w) + 2*cov(a,w).
var_decomp <- function(upstream, across, within, revenue, min_n = MIN_N_STATS) {
  ok <- !is.na(upstream) & !is.na(revenue) & upstream > 0 & revenue > 0
  n <- sum(ok)
  if (n < min_n) {
    return(data.table(
      n_firms      = n,
      var_total    = NA_real_,
      share_across = NA_real_,
      share_within = NA_real_,
      share_cov    = NA_real_
    ))
  }
  u <- upstream[ok] / revenue[ok]
  a <- across[ok]   / revenue[ok]
  w <- within[ok]   / revenue[ok]

  vt <- var(u)
  va <- var(a)
  vw <- var(w)
  ca <- cov(a, w)

  data.table(
    n_firms      = n,
    var_total    = vt,
    share_across = va / vt,
    share_within = vw / vt,
    share_cov    = 2 * ca / vt
  )
}

compute_vd <- function(firms_dt) {
  # firms_dt must have: nace2d, nace5d, crf_group, upstream,
  #   upstream_across_2d/5d/crf, upstream_within_2d/5d/crf, revenue
  firms_dt <- firms_dt[!is.na(nace2d)]

  # NACE 2-digit decomposition
  vd2d <- firms_dt[upstream > 0 | scope1 > 0,
    var_decomp(upstream, upstream_across_2d, upstream_within_2d, revenue),
    by = nace2d
  ][n_firms >= MIN_N_STATS]

  # NACE 5-digit decomposition
  vd5d <- firms_dt[upstream > 0 | scope1 > 0,
    var_decomp(upstream, upstream_across_5d, upstream_within_5d, revenue),
    by = nace5d
  ][n_firms >= MIN_N_STATS]

  # CRF group decomposition
  vdcrf <- firms_dt[upstream > 0 | scope1 > 0,
    var_decomp(upstream, upstream_across_crf, upstream_within_crf, revenue),
    by = crf_group
  ][n_firms >= MIN_N_STATS]

  list(vd2d = vd2d, vd5d = vd5d, vdcrf = vdcrf)
}


# =============================================================================
# SECTION 3: Process upstream files — point estimate (draw 0) + perturbation
# =============================================================================

cat("-- Processing upstream files ----------------------------------------\n")

GRANULARITIES <- c("vd2d", "vd5d", "vdcrf")
GROUP_COLS    <- c(vd2d = "nace2d", vd5d = "nace5d", vdcrf = "crf_group")

point_lists <- setNames(lapply(GRANULARITIES, function(x) list()), GRANULARITIES)
draw_lists  <- setNames(lapply(GRANULARITIES, function(x) list()), GRANULARITIES)
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
                      accounts_t[, .(vat, nace2d, nace5d, crf_group, revenue)],
                      by = "vat", all.x = TRUE)

    vd <- compute_vd(firms_b)

    for (g in GRANULARITIES) {
      vd_g <- vd[[g]]
      if (nrow(vd_g) == 0) next

      if (b == 1L) {
        vd_g[, year := t]
        point_lists[[g]][[length(point_lists[[g]]) + 1]] <- vd_g
      } else {
        vd_g[, `:=`(year = t, draw = b - 1L)]
        draw_lists[[g]][[length(draw_lists[[g]]) + 1]] <- vd_g
      }
    }
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("  Year %d: %d draws (%.1fs elapsed)\n", t, n_draws, elapsed))

  rm(firms_by_draw, conv_info)
  gc()
}

vd_point_2d  <- rbindlist(point_lists[["vd2d"]])
vd_point_5d  <- rbindlist(point_lists[["vd5d"]])
vd_point_crf <- rbindlist(point_lists[["vdcrf"]])
vd_draws_2d  <- rbindlist(draw_lists[["vd2d"]])
vd_draws_5d  <- rbindlist(draw_lists[["vd5d"]])
vd_draws_crf <- rbindlist(draw_lists[["vdcrf"]])

cat(sprintf("\n  Point — 2d: %d | 5d: %d | crf: %d sector-years\n",
            nrow(vd_point_2d), nrow(vd_point_5d), nrow(vd_point_crf)))
cat(sprintf("  Draws — 2d: %d | 5d: %d | crf: %d\n\n",
            nrow(vd_draws_2d), nrow(vd_draws_5d), nrow(vd_draws_crf)))


# =============================================================================
# SECTION 4: Summarize draws into prediction intervals
# =============================================================================

cat("-- Summarizing draws into prediction intervals ----------------------\n")

stat_cols <- c("var_total", "share_across", "share_within", "share_cov")

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

vd_summary_2d  <- merge(build_pi(vd_draws_2d,  "nace2d"),    vd_point_2d,
                         by = c("nace2d", "year"),    all.x = TRUE, suffixes = c("", "_point"))
vd_summary_5d  <- merge(build_pi(vd_draws_5d,  "nace5d"),    vd_point_5d,
                         by = c("nace5d", "year"),    all.x = TRUE, suffixes = c("", "_point"))
vd_summary_crf <- merge(build_pi(vd_draws_crf, "crf_group"), vd_point_crf,
                         by = c("crf_group", "year"), all.x = TRUE, suffixes = c("", "_point"))

cat(sprintf("  Summary — 2d: %d | 5d: %d | crf: %d sector-years with PIs\n\n",
            nrow(vd_summary_2d), nrow(vd_summary_5d), nrow(vd_summary_crf)))


# =============================================================================
# SECTION 5: Diagnostics
# =============================================================================

cat("-- Diagnostics -------------------------------------------------------\n\n")

cat("Variance decomposition shares (NACE 2d, point estimate, averaged across years):\n")
vd_avg <- vd_point_2d[, .(
  share_across = round(mean(share_across, na.rm = TRUE), 3),
  share_within = round(mean(share_within, na.rm = TRUE), 3),
  share_cov    = round(mean(share_cov,    na.rm = TRUE), 3),
  n_years      = .N
), by = nace2d][order(nace2d)]
print(vd_avg)

cat("\nVariance decomposition shares (CRF, point estimate, averaged across years):\n")
vd_avg_crf <- vd_point_crf[, .(
  share_across = round(mean(share_across, na.rm = TRUE), 3),
  share_within = round(mean(share_within, na.rm = TRUE), 3),
  share_cov    = round(mean(share_cov,    na.rm = TRUE), 3),
  n_years      = .N
), by = crf_group][order(crf_group)]
print(vd_avg_crf)

cat("\n95% PI width for share_within (NACE 2d):\n")
vd_summary_2d[, sw_width := share_within_hi95 - share_within_lo95]
cat(sprintf("  Median: %.4f\n", median(vd_summary_2d$sw_width, na.rm = TRUE)))


# =============================================================================
# SECTION 6: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA,
  sprintf("upstream_variance_decomp_pi_%s.RData", WEIGHT_SCHEME))

save(
  vd_point_2d,  vd_point_5d,  vd_point_crf,
  vd_draws_2d,  vd_draws_5d,  vd_draws_crf,
  vd_summary_2d, vd_summary_5d, vd_summary_crf,
  WEIGHT_SCHEME, MIN_N_STATS, stat_cols,
  file = OUT_PATH
)

cat(sprintf("\n===================================================================\n"))
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  2d  — point: %d | draws: %d | summary: %d\n",
            nrow(vd_point_2d),  nrow(vd_draws_2d),  nrow(vd_summary_2d)))
cat(sprintf("  5d  — point: %d | draws: %d | summary: %d\n",
            nrow(vd_point_5d),  nrow(vd_draws_5d),  nrow(vd_summary_5d)))
cat(sprintf("  crf — point: %d | draws: %d | summary: %d\n",
            nrow(vd_point_crf), nrow(vd_draws_crf), nrow(vd_summary_crf)))
cat("===================================================================\n")
