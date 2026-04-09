###############################################################################
# analysis/calibrate_extensive_margin.R
#
# PURPOSE
#   Fit logistic regressions of true emitter status on (p_i, log(1+proxy_mean))
#   on the mixed-CRF training subset:
#
#     logit Pr(D_i = 1) = alpha + beta_p * p_i + beta_lp * log(1 + proxy_mean_i)
#
#   Two weighting schemes are fit and saved side by side:
#
#     - "balanced":   each mixed CRF is weighted equally (down-weights paper,
#                     up-weights refining and metals). Primary classifier
#                     reported in the main results.
#     - "unweighted": vanilla logistic. Reported as a robustness check.
#
#   Each scheme produces:
#     - Pooled fit on all mixed CRFs (used for deployment classification)
#     - LOSO fits per held-out mixed CRF (used to assign Dhat to mixed-CRF
#       training observations honestly out-of-sample)
#     - Youden-optimal q* thresholds (pooled and per-LOSO-fold)
#     - q(p_i, proxy_mean_i) prediction function
#
# INPUT
#   {PROC_DATA}/training_summary_objects.RData
#
# OUTPUT
#   {PROC_DATA}/extensive_margin_calibration.RData
#     Per scheme s in {"balanced", "unweighted"}, the saved object contains:
#       ext_fit_pooled[[s]]    : glm fit on all mixed CRFs
#       ext_coefs[[s]]         : named numeric (alpha, beta_p, beta_lp)
#       ext_q_function[[s]]    : function(p, proxy_mean) returning q
#       q_star_pooled[[s]]     : Youden-optimal q* from the pooled fit
#       loso_fits[[s]]         : list of glm fits per held-out mixed CRF
#       loso_q_star[[s]]       : named numeric of Youden q* per LOSO fit
#       loso_diagnostics[[s]]  : per-fold held-out TPR/FPR/PPV/J table
#       ext_diagnostics[[s]]   : pooled implied metrics + held-out summary
#     Plus: MIXED_CRFS, SCHEMES
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
cat("  Schemes: balanced (primary), unweighted (robustness)\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load
# =============================================================================

load(file.path(PROC_DATA, "training_summary_objects.RData"))

MIXED_CRFS <- c("paper", "refining", "metals")
SCHEMES    <- c("balanced", "unweighted")

mixed <- training_summary[training_summary$primary_crf_group %in% MIXED_CRFS, ]
mixed$lp <- log1p(pmax(mixed$proxy_mean_i, 0))

# CRF-balanced weights: each CRF gets equal total weight, normalized so the
# sum of weights equals N (so coefficient SEs are on the same scale as the
# unweighted fit).
crf_n <- table(mixed$primary_crf_group)
mixed$w_balanced <- 1 / unname(crf_n[as.character(mixed$primary_crf_group)])
mixed$w_balanced <- mixed$w_balanced * nrow(mixed) / sum(mixed$w_balanced)

cat("Mixed-CRF firm-years:", nrow(mixed), "\n")
cat("  True emitters:", sum(mixed$D_i == 1), "\n")
cat("  Non-emitters: ", sum(mixed$D_i == 0), "\n\n")
cat("CRF sizes (raw):\n")
print(crf_n)
cat("\n")


# =============================================================================
# SECTION 2: Helper — Youden-optimal q* on a probability vector
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
# SECTION 3: Per-scheme calibration function
# =============================================================================

calibrate_scheme <- function(scheme) {
  cat(sprintf("\n###################################################\n"))
  cat(sprintf("  Scheme: %s\n", scheme))
  cat(sprintf("###################################################\n\n"))

  if (scheme == "balanced") {
    weights_vec <- mixed$w_balanced
  } else if (scheme == "unweighted") {
    weights_vec <- rep(1, nrow(mixed))
  } else {
    stop("Unknown scheme: ", scheme)
  }

  # ---- Pooled fit ----
  cat("-- Pooled fit on all mixed CRFs -------------------------------------\n")
  fit_pooled <- suppressWarnings(
    glm(D_i ~ p_i + lp, data = mixed, weights = weights_vec,
        family = binomial(link = "logit"))
  )
  coefs <- coef(fit_pooled)
  names(coefs) <- c("alpha", "beta_p", "beta_lp")
  cat("Coefficients:\n"); print(round(coefs, 4))

  q_pooled <- predict(fit_pooled, type = "response")
  fit_q <- youden_max_q(q_pooled, mixed$D_i)
  q_star_pooled <- fit_q$q_star
  cat(sprintf("\nYouden-optimal q* (pooled, in-sample): %.4f\n", q_star_pooled))
  cat(sprintf("  In-sample J = %.4f, TPR = %.4f, FPR = %.4f\n",
              fit_q$j, fit_q$tpr, fit_q$fpr))

  # Build a closure capturing this scheme's coefficients
  q_function <- local({
    .coefs <- coefs
    function(p, proxy_mean) {
      lp <- log1p(pmax(proxy_mean, 0))
      eta <- .coefs["alpha"] + .coefs["beta_p"] * p + .coefs["beta_lp"] * lp
      unname(plogis(eta))
    }
  })

  # ---- LOSO fits ----
  cat("\n-- LOSO fits (for honest training-side branch A/B sorting) ----------\n\n")

  loso_fits <- list()
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
    train_mask <- mixed$primary_crf_group != held_out
    train <- mixed[train_mask, ]
    test  <- mixed[!train_mask, ]

    if (scheme == "balanced") {
      # Recompute balanced weights on the LOSO training subset (two CRFs only)
      crf_n_train <- table(train$primary_crf_group)
      tw <- 1 / unname(crf_n_train[as.character(train$primary_crf_group)])
      tw <- tw * nrow(train) / sum(tw)
    } else {
      tw <- rep(1, nrow(train))
    }

    fit <- suppressWarnings(
      glm(D_i ~ p_i + lp, data = train, weights = tw,
          family = binomial(link = "logit"))
    )
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

  # ---- TPR check on non-mixed CRFs ----
  cat("\n-- TPR of pooled q* on non-mixed CRFs (true emitters only) ----------\n\n")
  nonmixed <- training_summary[!(training_summary$primary_crf_group %in% MIXED_CRFS) &
                                training_summary$D_i == 1, ]
  nonmixed$q <- q_function(nonmixed$p_i, nonmixed$proxy_mean_i)
  cat(sprintf("  Non-mixed true emitters: %d\n", nrow(nonmixed)))
  cat(sprintf("  TPR at q* = %.4f: %.4f\n\n",
              q_star_pooled, mean(nonmixed$q >= q_star_pooled)))
  for (crf in sort(unique(nonmixed$primary_crf_group))) {
    sub <- nonmixed[nonmixed$primary_crf_group == crf, ]
    cat(sprintf("    %-12s: %4d emitters, TPR = %.3f\n",
                crf, nrow(sub), mean(sub$q >= q_star_pooled)))
  }

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

  list(
    fit_pooled       = fit_pooled,
    coefs            = coefs,
    q_function       = q_function,
    q_star_pooled    = q_star_pooled,
    loso_fits        = loso_fits,
    loso_q_star      = loso_q_star,
    loso_diagnostics = loso_diagnostics,
    ext_diagnostics  = ext_diagnostics
  )
}


# =============================================================================
# SECTION 4: Run both schemes
# =============================================================================

results <- lapply(setNames(SCHEMES, SCHEMES), calibrate_scheme)


# =============================================================================
# SECTION 5: Save
# =============================================================================

ext_fit_pooled    <- lapply(results, `[[`, "fit_pooled")
ext_coefs         <- lapply(results, `[[`, "coefs")
ext_q_function    <- lapply(results, `[[`, "q_function")
q_star_pooled     <- lapply(results, `[[`, "q_star_pooled")
loso_fits         <- lapply(results, `[[`, "loso_fits")
loso_q_star       <- lapply(results, `[[`, "loso_q_star")
loso_diagnostics  <- lapply(results, `[[`, "loso_diagnostics")
ext_diagnostics   <- lapply(results, `[[`, "ext_diagnostics")

OUT_PATH <- file.path(PROC_DATA, "extensive_margin_calibration.RData")
save(ext_fit_pooled, ext_coefs, ext_q_function,
     q_star_pooled, loso_fits, loso_q_star, loso_diagnostics,
     ext_diagnostics, MIXED_CRFS, SCHEMES,
     file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
for (s in SCHEMES) {
  cat(sprintf("  [%s] alpha = %.3f, beta_p = %.3f, beta_lp = %.3f, q* = %.4f\n",
              s, ext_coefs[[s]]["alpha"], ext_coefs[[s]]["beta_p"],
              ext_coefs[[s]]["beta_lp"], q_star_pooled[[s]]))
}
cat("===================================================================\n")
