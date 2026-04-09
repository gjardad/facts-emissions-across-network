###############################################################################
# analysis/calibrate_extensive_margin.R
#
# PURPOSE
#   Fit a logistic regression of true emitter status on (p_i, log(1+proxy_mean))
#   on the mixed-CRF training subset:
#
#     logit Pr(D_i = 1) = alpha + beta_1 * p_i + beta_2 * log(1 + proxy_mean_i)
#
#   Produces:
#     - Pooled-mixed fit (used in deployment for classification and Bernoulli)
#     - LOSO fits per held-out mixed CRF (used to assign Dhat to mixed CRF
#       training observations honestly out-of-sample, for branch A/B sorting)
#     - Youden-optimal q* thresholds (pooled and per-LOSO-fold)
#
#   Note: this script supersedes the old calibrate_threshold.R (which fit a
#   univariate threshold on p_i). The two roles are now combined.
#
# INPUT
#   {PROC_DATA}/training_summary_objects.RData
#
# OUTPUT
#   {PROC_DATA}/extensive_margin_calibration.RData
#     ext_fit_pooled  : glm fit on all mixed CRFs
#     ext_coefs       : named numeric (alpha, beta_p, beta_lp)
#     ext_q_function  : function(p_i, proxy_mean_i) returning q
#     q_star_pooled   : Youden-optimal q* from the pooled fit
#     loso_fits       : list of glm fits per held-out mixed CRF
#     loso_q_star     : named numeric of Youden q* per LOSO fit (for sorting
#                       training mixed-CRF observations into branches)
#     loso_diagnostics: per-fold held-out TPR/FPR/PPV/J table
#     ext_diagnostics : pooled implied metrics + held-out summary
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
cat("  EXTENSIVE MARGIN CALIBRATION\n")
cat("  logit Pr(D=1) = alpha + b_p * p_i + b_lp * log(1+proxy_mean)\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load
# =============================================================================

load(file.path(PROC_DATA, "training_summary_objects.RData"))

MIXED_CRFS <- c("paper", "refining", "metals")
mixed <- training_summary[training_summary$primary_crf_group %in% MIXED_CRFS, ]
mixed$lp <- log1p(pmax(mixed$proxy_mean_i, 0))

cat("Mixed-CRF firm-years:", nrow(mixed), "\n")
cat("  True emitters:", sum(mixed$D_i == 1), "\n")
cat("  Non-emitters: ", sum(mixed$D_i == 0), "\n\n")


# =============================================================================
# SECTION 2: Helper — find Youden-optimal q* on a fitted probability vector
# =============================================================================

youden_max_q <- function(q, emit) {
  grid <- sort(unique(round(q, 6)))
  best_t <- NA_real_; best_j <- -Inf
  best_tpr <- NA_real_; best_fpr <- NA_real_
  np <- sum(emit == 1); nn <- sum(emit == 0)
  for (t in grid) {
    pred <- q >= t
    tpr <- sum(pred & emit == 1) / np
    fpr <- sum(pred & emit == 0) / nn
    j <- tpr - fpr
    if (j > best_j) {
      best_j <- j; best_t <- t; best_tpr <- tpr; best_fpr <- fpr
    }
  }
  list(q_star = best_t, j = best_j, tpr = best_tpr, fpr = best_fpr)
}


# =============================================================================
# SECTION 3: Pooled fit (deployment classifier)
# =============================================================================

cat("-- Pooled fit on all mixed CRFs -------------------------------------\n")

ext_fit_pooled <- glm(D_i ~ p_i + lp, data = mixed,
                       family = binomial(link = "logit"))

ext_coefs <- coef(ext_fit_pooled)
names(ext_coefs) <- c("alpha", "beta_p", "beta_lp")

cat("\nCoefficients:\n")
print(round(ext_coefs, 4))

q_pooled <- predict(ext_fit_pooled, type = "response")
fit_q <- youden_max_q(q_pooled, mixed$D_i)
q_star_pooled <- fit_q$q_star

cat(sprintf("\nYouden-optimal q* (pooled, in-sample): %.4f\n", q_star_pooled))
cat(sprintf("  In-sample J = %.4f, TPR = %.4f, FPR = %.4f\n",
            fit_q$j, fit_q$tpr, fit_q$fpr))

ext_q_function <- function(p, proxy_mean, coefs = ext_coefs) {
  lp <- log1p(pmax(proxy_mean, 0))
  eta <- coefs["alpha"] + coefs["beta_p"] * p + coefs["beta_lp"] * lp
  unname(plogis(eta))
}


# =============================================================================
# SECTION 4: LOSO fits per mixed CRF
# =============================================================================

cat("\n-- LOSO fits (for honest training-side branch A/B sorting) ----------\n\n")

loso_fits   <- list()
loso_q_star <- setNames(numeric(length(MIXED_CRFS)), MIXED_CRFS)
loso_diagnostics <- data.frame(
  held_out  = character(0),
  q_star    = numeric(0),
  train_j   = numeric(0),
  train_tpr = numeric(0),
  train_fpr = numeric(0),
  test_tpr  = numeric(0),
  test_fpr  = numeric(0),
  test_ppv  = numeric(0),
  test_j    = numeric(0),
  stringsAsFactors = FALSE
)

for (held_out in MIXED_CRFS) {
  train <- mixed[mixed$primary_crf_group != held_out, ]
  test  <- mixed[mixed$primary_crf_group == held_out, ]

  fit <- glm(D_i ~ p_i + lp, data = train, family = binomial(link = "logit"))
  loso_fits[[held_out]] <- fit

  q_train <- predict(fit, newdata = train, type = "response")
  fit_qt  <- youden_max_q(q_train, train$D_i)
  loso_q_star[held_out] <- fit_qt$q_star

  q_test  <- predict(fit, newdata = test, type = "response")
  pred_test <- q_test >= fit_qt$q_star
  test_tpr <- sum(pred_test & test$D_i == 1) / sum(test$D_i == 1)
  test_fpr <- sum(pred_test & test$D_i == 0) / sum(test$D_i == 0)
  test_tp  <- sum(pred_test & test$D_i == 1)
  test_fp  <- sum(pred_test & test$D_i == 0)
  test_ppv <- if (test_tp + test_fp > 0) test_tp / (test_tp + test_fp) else NA_real_

  loso_diagnostics <- rbind(loso_diagnostics, data.frame(
    held_out  = held_out,
    q_star    = fit_qt$q_star,
    train_j   = fit_qt$j, train_tpr = fit_qt$tpr, train_fpr = fit_qt$fpr,
    test_tpr  = test_tpr, test_fpr = test_fpr, test_ppv = test_ppv,
    test_j    = test_tpr - test_fpr,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  Held out %-10s q* = %.4f | train J = %.3f | test: TPR = %.3f, FPR = %.3f, PPV = %.3f, J = %.3f\n",
              held_out, fit_qt$q_star, fit_qt$j,
              test_tpr, test_fpr, test_ppv, test_tpr - test_fpr))
}

cat(sprintf("\n  Average held-out: TPR = %.4f, FPR = %.4f, PPV = %.4f, J = %.4f\n",
            mean(loso_diagnostics$test_tpr),
            mean(loso_diagnostics$test_fpr),
            mean(loso_diagnostics$test_ppv),
            mean(loso_diagnostics$test_j)))


# =============================================================================
# SECTION 5: TPR check on non-mixed CRFs (using pooled q*)
# =============================================================================

cat("\n-- TPR of pooled q* on non-mixed CRFs (true emitters only) ----------\n\n")

nonmixed <- training_summary[!(training_summary$primary_crf_group %in% MIXED_CRFS) &
                              training_summary$D_i == 1, ]
nonmixed$q <- ext_q_function(nonmixed$p_i, nonmixed$proxy_mean_i)

cat(sprintf("  Non-mixed true emitters: %d\n", nrow(nonmixed)))
cat(sprintf("  TPR at q* = %.4f: %.4f\n\n",
            q_star_pooled, mean(nonmixed$q >= q_star_pooled)))

cat("  By CRF:\n")
for (crf in sort(unique(nonmixed$primary_crf_group))) {
  sub <- nonmixed[nonmixed$primary_crf_group == crf, ]
  cat(sprintf("    %-12s: %4d emitters, TPR = %.3f\n",
              crf, nrow(sub), mean(sub$q >= q_star_pooled)))
}


# =============================================================================
# SECTION 6: Save
# =============================================================================

ext_diagnostics <- list(
  n           = nrow(mixed),
  n_pos       = sum(mixed$D_i == 1),
  n_neg       = sum(mixed$D_i == 0),
  in_sample_J = fit_q$j,
  loso_avg_tpr = mean(loso_diagnostics$test_tpr),
  loso_avg_fpr = mean(loso_diagnostics$test_fpr),
  loso_avg_ppv = mean(loso_diagnostics$test_ppv),
  loso_avg_j   = mean(loso_diagnostics$test_j)
)

OUT_PATH <- file.path(PROC_DATA, "extensive_margin_calibration.RData")
save(ext_fit_pooled, ext_coefs, ext_q_function,
     q_star_pooled, loso_fits, loso_q_star, loso_diagnostics,
     ext_diagnostics, MIXED_CRFS,
     file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  alpha = %.3f, beta_p = %.3f, beta_lp = %.3f\n",
            ext_coefs["alpha"], ext_coefs["beta_p"], ext_coefs["beta_lp"]))
cat(sprintf("  q_star_pooled = %.4f\n", q_star_pooled))
cat("  loso_q_star:", paste(sprintf("%s=%.4f", names(loso_q_star), loso_q_star),
                            collapse = ", "), "\n")
cat("===================================================================\n")
