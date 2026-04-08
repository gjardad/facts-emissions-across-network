load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel
panel$p_i <- rowMeans(proxy_matrix > 0)

mixed <- panel[panel$primary_crf_group %in% c("paper", "refining", "metals"), ]
sectors <- c("paper", "refining", "metals")

# Year-demean p_i within a given sample
year_demean <- function(df) {
  df$p_resid <- ave(df$p_i, df$year, FUN = function(x) x - mean(x))
  df
}

# In-sample Youden on a vector of p and emit
youden_max <- function(p, emit) {
  grid <- sort(unique(p))
  best_tau <- NA; best_j <- -Inf
  best_tpr <- NA; best_fpr <- NA
  for (tau in grid) {
    pred <- p >= tau
    tpr <- sum(pred & emit == 1) / sum(emit == 1)
    fpr <- sum(pred & emit == 0) / sum(emit == 0)
    j   <- tpr - fpr
    if (j > best_j) {
      best_j <- j; best_tau <- tau; best_tpr <- tpr; best_fpr <- fpr
    }
  }
  list(tau = best_tau, j = best_j, tpr = best_tpr, fpr = best_fpr)
}

cat("=== LOSO with year-demeaned pooled training ===\n\n")

results <- list()
for (held_out in sectors) {
  train <- mixed[mixed$primary_crf_group != held_out, ]
  test  <- mixed[mixed$primary_crf_group == held_out, ]

  # Year-demean within train (using train's own per-year means)
  train <- year_demean(train)

  # Find threshold on train's residualized scale
  fit <- youden_max(train$p_resid, train$emit)

  # Year-demean test using test's own per-year means (deployment-style:
  # the held-out sector uses its own year structure)
  test <- year_demean(test)

  # Apply threshold
  pred_test <- test$p_resid >= fit$tau
  tpr_test  <- sum(pred_test & test$emit == 1) / sum(test$emit == 1)
  fpr_test  <- sum(pred_test & test$emit == 0) / sum(test$emit == 0)
  j_test    <- tpr_test - fpr_test

  cat(sprintf("Held out: %-10s | tau* = %7.4f (resid scale)\n", held_out, fit$tau))
  cat(sprintf("  Train (pooled, residualized): J = %.4f, TPR = %.4f, FPR = %.4f\n",
              fit$j, fit$tpr, fit$fpr))
  cat(sprintf("  Test (held-out, residualized): J = %.4f, TPR = %.4f, FPR = %.4f\n\n",
              j_test, tpr_test, fpr_test))

  results[[held_out]] <- list(tau = fit$tau, j_test = j_test,
                                tpr_test = tpr_test, fpr_test = fpr_test)
}

cat("Average held-out J:", mean(sapply(results, `[[`, "j_test")), "\n")
cat("Average held-out TPR:", mean(sapply(results, `[[`, "tpr_test")), "\n")
cat("Average held-out FPR:", mean(sapply(results, `[[`, "fpr_test")), "\n")
cat("Average tau (resid scale):", mean(sapply(results, `[[`, "tau")), "\n")
