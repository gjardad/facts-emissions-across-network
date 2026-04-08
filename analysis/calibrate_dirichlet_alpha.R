###############################################################################
# analysis/calibrate_dirichlet_alpha.R
#
# PURPOSE
#   Build a lookup table alpha(N) for the Dirichlet share perturbation used
#   in prediction intervals.
#
#   For each group size N, find the concentration multiplier alpha such that:
#     w* ~ Dirichlet(alpha * w_ref)
#     Spearman(rank(w_ref), rank(w*)) ≈ RHO_TARGET  (averaged over M sims)
#
#   where w_ref are shares derived from empirical sinh(proxy) quantiles and
#   RHO_TARGET = 0.329 (average OOS rank correlation from CV).
#
# INPUT
#   {PROC_DATA}/repeated_cv_proxy_crf_asinh.RData — for empirical proxy dist
#
# OUTPUT
#   {PROC_DATA}/dirichlet_alpha_lookup.RData
#     alpha_lookup : named numeric vector, alpha per N
#     N_grid       : integer vector of group sizes
#     RHO_TARGET   : scalar
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

cat("===================================================================\n")
cat("  DIRICHLET ALPHA CALIBRATION\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load empirical sinh(proxy) distribution
# =============================================================================

cat("Loading empirical proxy distribution...\n")
load(file.path(PROC_DATA, "repeated_cv_proxy_crf_asinh.RData"))

pos_proxy <- repeated_cv_proxy_panel$proxy_mean[repeated_cv_proxy_panel$proxy_mean > 0]
sinh_proxy_empirical <- sinh(pos_proxy)
cat("  Positive-proxy firms:", length(pos_proxy), "\n")
cat("  sinh(proxy) range:", sprintf("[%.2e, %.2e]\n\n",
    min(sinh_proxy_empirical), max(sinh_proxy_empirical)))

rm(proxy_matrix, repeated_cv_proxy_panel, crf_group_map,
   firmyear_index, repeat_diagnostics, repeat_timing, syt)


# =============================================================================
# SECTION 2: Parameters
# =============================================================================

N_grid     <- c(2L, 3L, 4L, 5L, 7L, 10L, 15L, 20L, 30L, 50L, 75L,
                100L, 150L, 200L, 300L, 500L, 750L, 1000L)
RHO_TARGET <- 0.329
M_SIM      <- 2000L   # Dirichlet draws per evaluation
BISECT_TOL <- 0.005
BISECT_MAX <- 40L


# =============================================================================
# SECTION 3: Helper functions
# =============================================================================

# Reference shares for N firms from empirical sinh(proxy) quantiles
make_ref_shares <- function(N, sinh_empirical) {
  p_i    <- (seq_len(N) - 0.5) / N
  shares <- quantile(sinh_empirical, probs = p_i, names = FALSE)
  shares / sum(shares)
}

# Average Spearman between reference ranks and M Dirichlet-perturbed ranks
avg_spearman <- function(w_ref, alpha_mult, M) {
  N <- length(w_ref)
  if (N < 2) return(1.0)

  alpha_vec <- alpha_mult * w_ref

  # Vectorized Dirichlet: M x N gamma draws, normalize rows
  gamma_mat <- matrix(rgamma(M * N, shape = rep(alpha_vec, each = M), rate = 1),
                       nrow = M, ncol = N)
  dir_mat <- gamma_mat / rowSums(gamma_mat)

  ref_ranks <- rank(w_ref, ties.method = "average")
  rho_sum <- 0
  for (m in seq_len(M)) {
    rho_sum <- rho_sum + cor(ref_ranks,
                              rank(dir_mat[m, ], ties.method = "average"),
                              method = "spearman")
  }
  rho_sum / M
}

# Bisection to find alpha_mult giving avg Spearman ≈ rho_target
find_alpha <- function(w_ref, rho_target, M, tol, max_iter) {
  N <- length(w_ref)
  if (N < 2) return(NA_real_)

  lo <- 0.01
  hi <- 5000

  rho_lo <- avg_spearman(w_ref, lo, M)
  rho_hi <- avg_spearman(w_ref, hi, M)

  # Expand bracket if needed
  while (rho_lo > rho_target && lo > 1e-6) {
    lo <- lo / 10
    rho_lo <- avg_spearman(w_ref, lo, M)
  }
  while (rho_hi < rho_target && hi < 1e8) {
    hi <- hi * 10
    rho_hi <- avg_spearman(w_ref, hi, M)
  }

  for (iter in seq_len(max_iter)) {
    mid <- exp((log(lo) + log(hi)) / 2)
    rho_mid <- avg_spearman(w_ref, mid, M)
    if (abs(rho_mid - rho_target) < tol) return(mid)
    if (rho_mid < rho_target) lo <- mid else hi <- mid
  }

  exp((log(lo) + log(hi)) / 2)
}


# =============================================================================
# SECTION 4: Main calibration
# =============================================================================

cat("Calibrating alpha for", length(N_grid), "group sizes, rho =", RHO_TARGET, "\n")
cat("Simulations per evaluation:", M_SIM, "\n\n")

alpha_lookup <- setNames(numeric(length(N_grid)), paste0("N=", N_grid))

t0 <- Sys.time()

for (i in seq_along(N_grid)) {
  N <- N_grid[i]
  w_ref <- make_ref_shares(N, sinh_proxy_empirical)

  t0_n <- Sys.time()
  alpha_lookup[i] <- find_alpha(w_ref, RHO_TARGET, M_SIM, BISECT_TOL, BISECT_MAX)
  elapsed <- round(difftime(Sys.time(), t0_n, units = "secs"), 1)
  cat(sprintf("  N = %4d: alpha = %8.2f (%s s)\n", N, alpha_lookup[i], elapsed))
}

total_time <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat(sprintf("\nCalibration complete in %.1f s\n\n", total_time))


# =============================================================================
# SECTION 5: Verification
# =============================================================================

cat("-- Verification (10,000 draws) --------------------------------------\n")

for (check_N in c(10, 50, 200)) {
  idx <- which.min(abs(N_grid - check_N))
  w_check <- make_ref_shares(N_grid[idx], sinh_proxy_empirical)
  set.seed(42)
  rho_verify <- avg_spearman(w_check, alpha_lookup[idx], 10000)
  cat(sprintf("  N=%d: alpha=%.2f, verified Spearman=%.4f (target=%.3f)\n",
              N_grid[idx], alpha_lookup[idx], rho_verify, RHO_TARGET))
}


# =============================================================================
# SECTION 6: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA, "dirichlet_alpha_lookup.RData")

save(alpha_lookup, N_grid, RHO_TARGET, M_SIM, file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat("  alpha_lookup:", length(alpha_lookup), "values\n")
cat("  N range:", paste(range(N_grid), collapse = "--"), "\n")
cat("  RHO_TARGET:", RHO_TARGET, "\n")
cat("  Total time:", total_time, "s\n")
cat("===================================================================\n")
