###############################################################################
# analysis/calibrate_threshold.R
#
# PURPOSE
#   Compute the deployed extensive-margin thresholds on p_i:
#
#     - tau_loso[paper], tau_loso[refining], tau_loso[metals]:
#         each obtained by holding out one mixed CRF, pooling the other two,
#         and maximizing Youden's J. Used in training to assign Dhat to
#         mixed-CRF observations (honest out-of-sample for branch A/B sorting).
#
#     - tau_pooled_mixed:
#         obtained by pooling all three mixed CRFs and maximizing Youden's J.
#         Used in training to assign Dhat to non-mixed-CRF observations, and
#         in deployment for ALL firms.
#
# INPUT
#   {PROC_DATA}/training_summary_objects.RData
#
# OUTPUT
#   {PROC_DATA}/threshold_calibration.RData
#     tau_loso         : named numeric vector (paper, refining, metals)
#     tau_pooled_mixed : scalar
#     loso_diagnostics : data.frame with per-fold tau, train J, test J/TPR/FPR
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
cat("  THRESHOLD CALIBRATION\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load
# =============================================================================

load(file.path(PROC_DATA, "training_summary_objects.RData"))

MIXED_CRFS <- c("paper", "refining", "metals")
mixed <- training_summary[training_summary$primary_crf_group %in% MIXED_CRFS, ]

cat("Mixed-CRF firm-years:", nrow(mixed), "\n")
cat("  by CRF:\n")
print(table(mixed$primary_crf_group, mixed$D_i))
cat("\n")


# =============================================================================
# SECTION 2: Helper — Youden's J max over a (p, emit) sample
# =============================================================================

youden_max <- function(p, emit) {
  grid <- sort(unique(p))
  best_tau <- NA_real_; best_j <- -Inf
  best_tpr <- NA_real_; best_fpr <- NA_real_
  n_pos <- sum(emit == 1)
  n_neg <- sum(emit == 0)

  for (tau in grid) {
    pred <- p >= tau
    tpr <- sum(pred & emit == 1) / n_pos
    fpr <- sum(pred & emit == 0) / n_neg
    j   <- tpr - fpr
    if (j > best_j) {
      best_j <- j; best_tau <- tau; best_tpr <- tpr; best_fpr <- fpr
    }
  }
  list(tau = best_tau, j = best_j, tpr = best_tpr, fpr = best_fpr)
}


# =============================================================================
# SECTION 3: LOSO thresholds for mixed CRFs
# =============================================================================

cat("-- LOSO thresholds (pool the other two mixed CRFs) ------------------\n\n")

loso_diagnostics <- data.frame(
  held_out  = character(0),
  tau       = numeric(0),
  train_j   = numeric(0),
  train_tpr = numeric(0),
  train_fpr = numeric(0),
  test_j    = numeric(0),
  test_tpr  = numeric(0),
  test_fpr  = numeric(0),
  stringsAsFactors = FALSE
)

tau_loso <- setNames(numeric(length(MIXED_CRFS)), MIXED_CRFS)

for (held_out in MIXED_CRFS) {
  train <- mixed[mixed$primary_crf_group != held_out, ]
  test  <- mixed[mixed$primary_crf_group == held_out, ]

  fit <- youden_max(train$p_i, train$D_i)

  pred_test <- test$p_i >= fit$tau
  test_tpr  <- sum(pred_test & test$D_i == 1) / sum(test$D_i == 1)
  test_fpr  <- sum(pred_test & test$D_i == 0) / sum(test$D_i == 0)
  test_j    <- test_tpr - test_fpr

  tau_loso[held_out] <- fit$tau

  loso_diagnostics <- rbind(loso_diagnostics, data.frame(
    held_out  = held_out,
    tau       = fit$tau,
    train_j   = fit$j, train_tpr = fit$tpr, train_fpr = fit$fpr,
    test_j    = test_j, test_tpr  = test_tpr,  test_fpr  = test_fpr,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  Held out %-10s tau = %.4f | train J = %.3f | test: J = %.3f, TPR = %.3f, FPR = %.3f\n",
              held_out, fit$tau, fit$j, test_j, test_tpr, test_fpr))
}

cat(sprintf("\n  Average held-out J: %.4f\n", mean(loso_diagnostics$test_j)))
cat(sprintf("  Average held-out TPR: %.4f\n", mean(loso_diagnostics$test_tpr)))
cat(sprintf("  Average held-out FPR: %.4f\n\n", mean(loso_diagnostics$test_fpr)))


# =============================================================================
# SECTION 4: Pooled-mixed threshold
# =============================================================================

cat("-- Pooled-mixed threshold (all three CRFs together) -----------------\n\n")

fit_pooled <- youden_max(mixed$p_i, mixed$D_i)
tau_pooled_mixed <- fit_pooled$tau

cat(sprintf("  tau_pooled_mixed = %.4f\n", tau_pooled_mixed))
cat(sprintf("  In-sample J = %.4f, TPR = %.4f, FPR = %.4f\n\n",
            fit_pooled$j, fit_pooled$tpr, fit_pooled$fpr))


# =============================================================================
# SECTION 5: TPR check on non-mixed CRFs (for documentation)
# =============================================================================

cat("-- TPR of tau_pooled_mixed on non-mixed CRFs (true emitters only) ---\n\n")

nonmixed <- training_summary[!(training_summary$primary_crf_group %in% MIXED_CRFS) &
                              training_summary$D_i == 1, ]
cat(sprintf("  Non-mixed true emitters: %d\n", nrow(nonmixed)))
cat(sprintf("  TPR at tau = %.4f: %.4f\n\n",
            tau_pooled_mixed, mean(nonmixed$p_i >= tau_pooled_mixed)))

cat("  By CRF:\n")
for (crf in sort(unique(nonmixed$primary_crf_group))) {
  sub <- nonmixed[nonmixed$primary_crf_group == crf, ]
  cat(sprintf("    %-12s: %4d emitters, TPR = %.3f\n",
              crf, nrow(sub), mean(sub$p_i >= tau_pooled_mixed)))
}


# =============================================================================
# SECTION 6: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA, "threshold_calibration.RData")
save(tau_loso, tau_pooled_mixed, loso_diagnostics, MIXED_CRFS, file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat("  tau_loso:        ", paste(sprintf("%s=%.4f", names(tau_loso), tau_loso), collapse = ", "), "\n")
cat("  tau_pooled_mixed:", round(tau_pooled_mixed, 4), "\n")
cat("===================================================================\n")
