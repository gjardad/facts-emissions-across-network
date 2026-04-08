load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel
panel$p_i <- rowMeans(proxy_matrix > 0)

mixed <- panel[panel$primary_crf_group %in% c("paper", "refining", "metals"), ]
sectors <- c("paper", "refining", "metals")

youden_max <- function(p, emit) {
  grid <- sort(unique(p))
  best_tau <- NA; best_j <- -Inf; best_tpr <- NA; best_fpr <- NA
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

cat("=== LOSO: pool the other two sectors (raw p_i, no demeaning) ===\n\n")

results <- list()
for (held_out in sectors) {
  train <- mixed[mixed$primary_crf_group != held_out, ]
  test  <- mixed[mixed$primary_crf_group == held_out, ]

  fit <- youden_max(train$p_i, train$emit)

  pred_test <- test$p_i >= fit$tau
  tpr_test  <- sum(pred_test & test$emit == 1) / sum(test$emit == 1)
  fpr_test  <- sum(pred_test & test$emit == 0) / sum(test$emit == 0)
  j_test    <- tpr_test - fpr_test

  cat(sprintf("Held out: %-10s | tau* = %.4f\n", held_out, fit$tau))
  cat(sprintf("  Train (pooled): J = %.4f, TPR = %.4f, FPR = %.4f\n",
              fit$j, fit$tpr, fit$fpr))
  cat(sprintf("  Test (held-out): J = %.4f, TPR = %.4f, FPR = %.4f\n\n",
              j_test, tpr_test, fpr_test))

  results[[held_out]] <- list(tau = fit$tau, j = j_test, tpr = tpr_test, fpr = fpr_test)
}

cat("Average tau:", mean(sapply(results, `[[`, "tau")), "\n")
cat("Average held-out J:", mean(sapply(results, `[[`, "j")), "\n")
cat("Average held-out TPR:", mean(sapply(results, `[[`, "tpr")), "\n")
cat("Average held-out FPR:", mean(sapply(results, `[[`, "fpr")), "\n")
