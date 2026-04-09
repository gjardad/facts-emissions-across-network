suppressPackageStartupMessages(library(data.table))
load("C:/Users/jota_/Documents/NBB_data/processed/training_summary_objects.RData")
setDT(training_summary)
mixed <- training_summary[primary_crf_group %in% c("paper","refining","metals")]
mixed[, lp := log1p(pmax(proxy_mean_i, 0))]
sectors <- c("paper","refining","metals")

youden_max_q <- function(q, emit) {
  grid <- sort(unique(round(q, 5)))
  best_t <- NA; best_j <- -Inf
  np <- sum(emit == 1); nn <- sum(emit == 0)
  for (t in grid) {
    pred <- q >= t
    tpr <- sum(pred & emit == 1) / np
    fpr <- sum(pred & emit == 0) / nn
    j <- tpr - fpr
    if (j > best_j) { best_j <- j; best_t <- t }
  }
  best_t
}

# Pool LOSO predictions across the three folds, then compute pooled PPV
pool_loso <- function(formula, label) {
  cat(sprintf("\n=== %s ===\n", label))
  pred_pool <- list()
  for (held_out in sectors) {
    train <- mixed[primary_crf_group != held_out]
    test  <- mixed[primary_crf_group == held_out]
    fit <- glm(formula, data = train, family = binomial(link = "logit"))
    q_train <- predict(fit, newdata = train, type = "response")
    q_star  <- youden_max_q(q_train, train$D_i)
    q_test  <- predict(fit, newdata = test, type = "response")
    pred_pool[[held_out]] <- data.table(
      sec   = held_out,
      D_i   = test$D_i,
      pred  = as.integer(q_test >= q_star)
    )
  }
  pp <- rbindlist(pred_pool)
  TP <- sum(pp$pred == 1 & pp$D_i == 1)
  FP <- sum(pp$pred == 1 & pp$D_i == 0)
  FN <- sum(pp$pred == 0 & pp$D_i == 1)
  TN <- sum(pp$pred == 0 & pp$D_i == 0)
  cat(sprintf("  Pooled: TP=%d FP=%d FN=%d TN=%d\n", TP, FP, FN, TN))
  cat(sprintf("  TPR = %.4f  FPR = %.4f  PPV = TP/(TP+FP) = %.4f\n",
              TP/(TP+FN), FP/(FP+TN), TP/(TP+FP)))
  cat("  PPV by held-out sector:\n")
  for (s in sectors) {
    sub <- pp[sec == s]
    tp <- sum(sub$pred == 1 & sub$D_i == 1)
    fp <- sum(sub$pred == 1 & sub$D_i == 0)
    cat(sprintf("    %-9s TP=%d FP=%d  PPV=%.4f\n",
                s, tp, fp, if (tp+fp > 0) tp/(tp+fp) else NA))
  }
}

# Reference: raw p_i threshold (deployed baseline)
cat("=== Reference: raw p_i threshold (LOSO Youden, pooled) ===\n")
youden_max_p <- function(p, emit) {
  grid <- sort(unique(p))
  best_t <- NA; best_j <- -Inf
  np <- sum(emit == 1); nn <- sum(emit == 0)
  for (t in grid) {
    pred <- p >= t
    tpr <- sum(pred & emit == 1) / np
    fpr <- sum(pred & emit == 0) / nn
    j <- tpr - fpr
    if (j > best_j) { best_j <- j; best_t <- t }
  }
  best_t
}
pp <- list()
for (held_out in sectors) {
  train <- mixed[primary_crf_group != held_out]
  test  <- mixed[primary_crf_group == held_out]
  tau <- youden_max_p(train$p_i, train$D_i)
  pp[[held_out]] <- data.table(
    sec  = held_out,
    D_i  = test$D_i,
    pred = as.integer(test$p_i >= tau)
  )
}
pp <- rbindlist(pp)
TP <- sum(pp$pred == 1 & pp$D_i == 1)
FP <- sum(pp$pred == 1 & pp$D_i == 0)
FN <- sum(pp$pred == 0 & pp$D_i == 1)
TN <- sum(pp$pred == 0 & pp$D_i == 0)
cat(sprintf("  Pooled: TP=%d FP=%d FN=%d TN=%d\n", TP, FP, FN, TN))
cat(sprintf("  TPR = %.4f  FPR = %.4f  PPV = %.4f\n",
            TP/(TP+FN), FP/(FP+TN), TP/(TP+FP)))
cat("  PPV by held-out sector:\n")
for (s in sectors) {
  sub <- pp[sec == s]
  tp <- sum(sub$pred == 1 & sub$D_i == 1)
  fp <- sum(sub$pred == 1 & sub$D_i == 0)
  cat(sprintf("    %-9s TP=%d FP=%d  PPV=%.4f\n",
              s, tp, fp, if (tp+fp > 0) tp/(tp+fp) else NA))
}

pool_loso(D_i ~ lp,        "Logistic on log(1+proxy_mean) alone")
pool_loso(D_i ~ p_i + lp,  "Logistic on p_i + log(1+proxy_mean)")
pool_loso(D_i ~ p_i * lp,  "Logistic on p_i * log(1+proxy_mean) (interaction)")
