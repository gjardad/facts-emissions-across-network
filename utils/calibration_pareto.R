###############################################################################
# utils/calibration_pareto.R
#
# GPA (Generalized Pareto) redistribution helpers for deployment pipeline.
# Ported from inferring_emissions/analysis/active/diagnostic_redistribution_comparison.R
#
# Functions:
#   build_reference_dist() — estimate within-sector log-emission deviations
#   fit_gpa()              — fit GPA via L-moments
#   pareto_weights()       — rank → GPA quantile → exp(w) shares
###############################################################################

library(lmom)

# ── Build reference distribution from training emitters ──────────────────────
# Year-demean, then sector-demean log(emissions). Pool across all sectors.
# Returns a sorted numeric vector of deviations d_j.
build_reference_dist <- function(emissions, year, nace2d) {
  # Only emitters
  ok <- !is.na(emissions) & emissions > 0
  if (sum(ok) < 20) return(NULL)

  log_y  <- log(emissions[ok])
  yr     <- year[ok]
  sec    <- nace2d[ok]

  # Year-demean
  mu_t   <- tapply(log_y, yr, mean)
  tilde  <- log_y - mu_t[as.character(yr)]

  # Sector-demean
  mu_s   <- tapply(tilde, sec, mean)
  d      <- tilde - mu_s[sec]

  sort(d)
}

# ── Fit GPA via L-moments ────────────────────────────────────────────────────
fit_gpa <- function(ref_dist) {
  if (is.null(ref_dist) || length(ref_dist) < 20) return(NULL)
  lmoms <- samlmu(ref_dist, nmom = 3)
  tryCatch(pelgpa(lmoms), error = function(e) NULL)
}

# ── Pareto weights from ranks ────────────────────────────────────────────────
# Given n_total positions in the sector-year distribution and the ranks of
# the firms to impute, return exp(w_i) weights (unnormalized).
#
# ranks:    integer vector of ranks for firms to impute (1 = lowest)
# n_total:  total number of emitters in the sector-year (ETS + non-ETS)
# gpa_params: output of fit_gpa() / pelgpa()
#
# Returns: numeric vector of exp(w_i) values, or NULL on failure.
pareto_weights <- function(ranks, n_total, gpa_params) {
  if (is.null(gpa_params) || length(ranks) == 0) return(NULL)

  # Hazen plotting positions
  p_i <- (ranks - 0.5) / n_total

  # Clamp to avoid boundary issues
  p_i <- pmin(pmax(p_i, 1e-6), 1 - 1e-6)

  w_i <- tryCatch(
    quagpa(p_i, gpa_params),
    error = function(e) NULL
  )
  if (is.null(w_i) || any(!is.finite(w_i))) return(NULL)

  # Numerical stability: subtract max before exp
  exp(w_i - max(w_i))
}
