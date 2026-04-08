###############################################################################
# analysis/check_pi_coverage.R
#
# PURPOSE
#   Check the calibration of the branch A latent-normal model on training
#   data. For each branch A firm-year:
#
#     mu_i    = a + b * zhat_i
#     sigma_i = sqrt(exp(d + h * zhat_i))
#
#   The model claims z_i ~ N(mu_i, sigma_i^2). Coverage check: what fraction
#   of z_i fall inside the predicted 95% / 90% / 50% intervals?
#
#   Also reports z-space residual diagnostics and a binned coverage table
#   stratified by zhat_i.
#
# INPUT
#   {PROC_DATA}/training_summary_objects.RData
#   {PROC_DATA}/threshold_calibration.RData
#   {PROC_DATA}/intensive_margin_calibration.RData
#
# OUTPUT
#   Console diagnostics + saved coverage table at
#   {PROC_DATA}/coverage_check.RData
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
cat("  COVERAGE CHECK: branch A latent-normal model in z-space\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load
# =============================================================================

load(file.path(PROC_DATA, "training_summary_objects.RData"))
load(file.path(PROC_DATA, "threshold_calibration.RData"))
load(file.path(PROC_DATA, "intensive_margin_calibration.RData"))

cat("Loaded calibration objects.\n")
cat(sprintf("  mean: a = %.4f, b = %.4f\n", mean_coefs["a"], mean_coefs["b"]))
cat(sprintf("  sigma (constant): %.4f\n\n", sigma_const))


# =============================================================================
# SECTION 2: Rebuild branch A sample (same as in calibrate_intensive_margin.R)
# =============================================================================

ts <- training_summary
ts$tau <- ifelse(ts$primary_crf_group %in% MIXED_CRFS,
                  unname(tau_loso[as.character(ts$primary_crf_group)]),
                  tau_pooled_mixed)
ts$Dhat <- as.integer(ts$p_i >= ts$tau)

emit <- ts[ts$D_i == 1, ]
br_a <- emit[emit$Dhat == 1, ]

cell_n <- ave(br_a$D_i, br_a$primary_crf_group, br_a$year, FUN = length)
br_a <- br_a[cell_n >= 3, ]

compute_pct <- function(x, n) {
  r <- rank(-x, ties.method = "average")
  (r - 0.5) / n
}

br_a$u_i    <- ave(br_a$y_i,          br_a$primary_crf_group, br_a$year,
                    FUN = function(x) compute_pct(x, length(x)))
br_a$uhat_i <- ave(br_a$proxy_mean_i, br_a$primary_crf_group, br_a$year,
                    FUN = function(x) compute_pct(x, length(x)))
br_a$z_i    <- qnorm(br_a$u_i)
br_a$zhat_i <- qnorm(br_a$uhat_i)

cat("Branch A sample:", nrow(br_a), "firm-years\n\n")


# =============================================================================
# SECTION 3: Predicted mean and sigma per firm-year
# =============================================================================

br_a$mu_i    <- mean_coefs["a"] + mean_coefs["b"] * br_a$zhat_i
br_a$sigma_i <- sigma_const

cat("Sigma summary:\n")
print(summary(br_a$sigma_i))


# =============================================================================
# SECTION 4: Coverage at standard PI widths
# =============================================================================

cat("\n-- Coverage of predicted N(mu, sigma^2) intervals --------------------\n\n")

levels <- c(0.50, 0.80, 0.90, 0.95)

cov_results <- data.frame(level = levels, nominal = levels, observed = NA_real_)
for (i in seq_along(levels)) {
  alpha <- 1 - levels[i]
  z_alpha <- qnorm(1 - alpha / 2)
  lo <- br_a$mu_i - z_alpha * br_a$sigma_i
  hi <- br_a$mu_i + z_alpha * br_a$sigma_i
  cov_results$observed[i] <- mean(br_a$z_i >= lo & br_a$z_i <= hi)
}

print(cov_results, row.names = FALSE)

cat("\n  Observed should be close to nominal if the latent-normal model is\n")
cat("  well-calibrated. Higher than nominal = intervals too wide; lower =\n")
cat("  intervals too narrow.\n")


# =============================================================================
# SECTION 5: Standardized residual distribution
# =============================================================================

cat("\n-- Standardized residuals: (z_i - mu_i) / sigma_i --------------------\n\n")

br_a$std_resid <- (br_a$z_i - br_a$mu_i) / br_a$sigma_i

cat("Summary:\n")
print(summary(br_a$std_resid))
cat(sprintf("\n  Mean:     %.4f  (should be ~0)\n", mean(br_a$std_resid)))
cat(sprintf("  SD:       %.4f  (should be ~1)\n", sd(br_a$std_resid)))
cat(sprintf("  |resid| > 1.96: %.4f  (should be ~0.05)\n",
            mean(abs(br_a$std_resid) > 1.96)))
cat(sprintf("  |resid| > 1.64: %.4f  (should be ~0.10)\n",
            mean(abs(br_a$std_resid) > 1.64)))


# =============================================================================
# SECTION 6: Binned coverage by zhat_i
# =============================================================================

cat("\n-- Binned coverage by zhat_i decile -------------------------------\n\n")

br_a$bin <- cut(br_a$zhat_i,
                breaks = quantile(br_a$zhat_i, probs = seq(0, 1, 0.1)),
                include.lowest = TRUE)

binned <- aggregate(br_a$z_i,
                     by = list(bin = br_a$bin),
                     FUN = function(z) {
                       idx <- match(z, br_a$z_i)
                       mu <- br_a$mu_i[idx]
                       sg <- br_a$sigma_i[idx]
                       cov95 <- mean(z >= mu - 1.96 * sg & z <= mu + 1.96 * sg)
                       cov90 <- mean(z >= mu - 1.64 * sg & z <= mu + 1.64 * sg)
                       c(n = length(z), cov95 = cov95, cov90 = cov90)
                     })

binned_df <- data.frame(
  bin   = binned$bin,
  n     = binned$x[, "n"],
  cov95 = binned$x[, "cov95"],
  cov90 = binned$x[, "cov90"]
)
print(binned_df, row.names = FALSE)


# =============================================================================
# SECTION 7: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA, "coverage_check.RData")
save(cov_results, binned_df, br_a, file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat("===================================================================\n")
