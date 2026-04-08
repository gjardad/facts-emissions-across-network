###############################################################################
# analysis/calibrate_extensive_margin.R
#
# PURPOSE
#   Fit a logistic regression of true emitter status on the selection
#   frequency p_i, using the mixed-CRF training subset.
#
#     logit Pr(D_i = 1) = alpha + beta * p_i
#
#   The fitted model is used in deployment to draw Bernoulli emitter
#   indicators around the deterministic baseline (Dhat = 1{p_i >= tau}).
#
# INPUT
#   {PROC_DATA}/training_summary_objects.RData
#
# OUTPUT
#   {PROC_DATA}/extensive_margin_calibration.RData
#     ext_fit         : fitted glm object
#     ext_coefs       : named numeric (alpha, beta)
#     ext_q_function  : function(p) returning Pr(D=1 | p)
#     ext_diagnostics : list with sample size, deviance, calibration table
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
cat("  EXTENSIVE MARGIN CALIBRATION (logit Pr(D=1) = a + b * p_i)\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load
# =============================================================================

load(file.path(PROC_DATA, "training_summary_objects.RData"))

MIXED_CRFS <- c("paper", "refining", "metals")
mixed <- training_summary[training_summary$primary_crf_group %in% MIXED_CRFS, ]

cat("Mixed-CRF firm-years:", nrow(mixed), "\n")
cat("  True emitters:", sum(mixed$D_i == 1), "\n")
cat("  Non-emitters: ", sum(mixed$D_i == 0), "\n\n")


# =============================================================================
# SECTION 2: Fit logistic regression
# =============================================================================

cat("Fitting logistic regression...\n")

ext_fit <- glm(D_i ~ p_i, data = mixed, family = binomial(link = "logit"))

ext_coefs <- coef(ext_fit)
names(ext_coefs) <- c("alpha", "beta")

cat("\nCoefficients:\n")
print(round(ext_coefs, 4))

cat("\nSummary:\n")
print(summary(ext_fit)$coefficients)

cat(sprintf("\nNull deviance:     %.1f on %d df\n",
            ext_fit$null.deviance, ext_fit$df.null))
cat(sprintf("Residual deviance: %.1f on %d df\n",
            ext_fit$deviance, ext_fit$df.residual))


# =============================================================================
# SECTION 3: Build q(p) function
# =============================================================================

ext_q_function <- function(p) {
  eta <- ext_coefs["alpha"] + ext_coefs["beta"] * p
  unname(plogis(eta))
}


# =============================================================================
# SECTION 4: Calibration diagnostics
# =============================================================================

cat("\n-- Calibration: predicted vs observed emitter rate by p_i bin -------\n")

mixed$q_pred <- ext_q_function(mixed$p_i)

# Bin by p_i
breaks <- c(-0.001, 0.05, 0.20, 0.40, 0.60, 0.80, 0.95, 1.001)
mixed$bin <- cut(mixed$p_i, breaks = breaks, include.lowest = TRUE)

calib_tab <- aggregate(cbind(D_i, q_pred) ~ bin, data = mixed,
                        FUN = function(x) c(n = length(x), rate = mean(x)))
calib_df <- data.frame(
  bin       = calib_tab$bin,
  n         = calib_tab$D_i[, "n"],
  obs_rate  = calib_tab$D_i[, "rate"],
  pred_rate = calib_tab$q_pred[, "rate"]
)
print(calib_df, row.names = FALSE)

cat("\n-- Implied aggregate metrics ----------------------------------------\n")
implied_TPR <- mean(mixed$q_pred[mixed$D_i == 1])
implied_FPR <- mean(mixed$q_pred[mixed$D_i == 0])
cat(sprintf("  E[q | D=1] = %.4f\n", implied_TPR))
cat(sprintf("  E[q | D=0] = %.4f\n", implied_FPR))


# =============================================================================
# SECTION 5: Save
# =============================================================================

ext_diagnostics <- list(
  n           = nrow(mixed),
  n_pos       = sum(mixed$D_i == 1),
  n_neg       = sum(mixed$D_i == 0),
  null_dev    = ext_fit$null.deviance,
  resid_dev   = ext_fit$deviance,
  calib_table = calib_df,
  implied_TPR = implied_TPR,
  implied_FPR = implied_FPR
)

OUT_PATH <- file.path(PROC_DATA, "extensive_margin_calibration.RData")
save(ext_fit, ext_coefs, ext_q_function, ext_diagnostics, file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  alpha = %.4f, beta = %.4f\n", ext_coefs["alpha"], ext_coefs["beta"]))
cat("===================================================================\n")
