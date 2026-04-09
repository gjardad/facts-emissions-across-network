suppressPackageStartupMessages(library(data.table))
load("C:/Users/jota_/Documents/NBB_data/processed/training_summary_objects.RData")
setDT(training_summary)

mixed <- training_summary[primary_crf_group %in% c("paper","refining","metals")]
mixed[, lp := log1p(pmax(proxy_mean_i, 0))]
sectors <- c("paper","refining","metals")

# Helper: find threshold q* on a vector of predicted probabilities
youden_max_q <- function(q, emit) {
  grid <- sort(unique(round(q, 5)))
  best_t <- NA; best_j <- -Inf; best_tpr <- NA; best_fpr <- NA
  np <- sum(emit == 1); nn <- sum(emit == 0)
  for (t in grid) {
    pred <- q >= t
    tpr <- sum(pred & emit == 1) / np
    fpr <- sum(pred & emit == 0) / nn
    j <- tpr - fpr
    if (j > best_j) { best_j <- j; best_t <- t; best_tpr <- tpr; best_fpr <- fpr }
  }
  list(t = best_t, j = best_j, tpr = best_tpr, fpr = best_fpr)
}

eval_classifier <- function(formula, label) {
  cat(sprintf("\n=== %s ===\n", label))
  results <- list()
  for (held_out in sectors) {
    train <- mixed[primary_crf_group != held_out]
    test  <- mixed[primary_crf_group == held_out]

    fit <- glm(formula, data = train, family = binomial(link = "logit"))

    # Tune q* on training (in-train Youden) — honest LOSO since training
    # doesn't see the held-out sector
    q_train <- predict(fit, newdata = train, type = "response")
    cut_fit <- youden_max_q(q_train, train$D_i)

    # Apply to held-out sector
    q_test <- predict(fit, newdata = test, type = "response")
    pred_t <- q_test >= cut_fit$t
    tpr_t  <- sum(pred_t & test$D_i == 1) / sum(test$D_i == 1)
    fpr_t  <- sum(pred_t & test$D_i == 0) / sum(test$D_i == 0)
    j_t    <- tpr_t - fpr_t

    cat(sprintf("  Held out %-9s | q* = %.4f | train J = %.3f | test: TPR = %.3f, FPR = %.3f, J = %.3f\n",
                held_out, cut_fit$t, cut_fit$j, tpr_t, fpr_t, j_t))

    results[[held_out]] <- list(j = j_t, tpr = tpr_t, fpr = fpr_t, q_star = cut_fit$t)
  }
  cat(sprintf("  Average held-out: J = %.4f, TPR = %.4f, FPR = %.4f\n",
              mean(sapply(results, `[[`, "j")),
              mean(sapply(results, `[[`, "tpr")),
              mean(sapply(results, `[[`, "fpr"))))
  invisible(results)
}

# Reference: LOSO threshold on raw p_i (matches the .md numbers)
cat("=== Reference: LOSO threshold on raw p_i ===\n")
youden_max_p <- function(p, emit) {
  grid <- sort(unique(p))
  best_t <- NA; best_j <- -Inf; best_tpr <- NA; best_fpr <- NA
  np <- sum(emit == 1); nn <- sum(emit == 0)
  for (t in grid) {
    pred <- p >= t
    tpr <- sum(pred & emit == 1) / np
    fpr <- sum(pred & emit == 0) / nn
    j <- tpr - fpr
    if (j > best_j) { best_j <- j; best_t <- t; best_tpr <- tpr; best_fpr <- fpr }
  }
  list(t = best_t, j = best_j, tpr = best_tpr, fpr = best_fpr)
}
ref_results <- list()
for (held_out in sectors) {
  train <- mixed[primary_crf_group != held_out]
  test  <- mixed[primary_crf_group == held_out]
  ft <- youden_max_p(train$p_i, train$D_i)
  pt <- test$p_i >= ft$t
  tt <- sum(pt & test$D_i == 1) / sum(test$D_i == 1)
  ff <- sum(pt & test$D_i == 0) / sum(test$D_i == 0)
  cat(sprintf("  Held out %-9s | tau = %.4f | train J = %.3f | test: TPR = %.3f, FPR = %.3f, J = %.3f\n",
              held_out, ft$t, ft$j, tt, ff, tt - ff))
  ref_results[[held_out]] <- list(j = tt - ff, tpr = tt, fpr = ff)
}
cat(sprintf("  Average held-out: J = %.4f, TPR = %.4f, FPR = %.4f\n",
            mean(sapply(ref_results, `[[`, "j")),
            mean(sapply(ref_results, `[[`, "tpr")),
            mean(sapply(ref_results, `[[`, "fpr"))))

# Logistic specifications
eval_classifier(D_i ~ p_i, "Logistic on p_i alone")
eval_classifier(D_i ~ lp,  "Logistic on log(1+proxy_mean) alone")
eval_classifier(D_i ~ p_i + lp, "Logistic on p_i + log(1+proxy_mean)")
eval_classifier(D_i ~ p_i * lp, "Logistic on p_i * log(1+proxy_mean) (with interaction)")
