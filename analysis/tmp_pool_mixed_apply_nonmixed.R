load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel
panel$p_i <- rowMeans(proxy_matrix > 0)

mixed_crfs <- c("paper", "refining", "metals")
mixed <- panel[panel$primary_crf_group %in% mixed_crfs, ]

# In-sample Youden on pooled mixed sectors
youden_max <- function(p, emit) {
  grid <- sort(unique(p))
  best_tau <- NA; best_j <- -Inf
  for (tau in grid) {
    pred <- p >= tau
    tpr <- sum(pred & emit == 1) / sum(emit == 1)
    fpr <- sum(pred & emit == 0) / sum(emit == 0)
    j   <- tpr - fpr
    if (j > best_j) { best_j <- j; best_tau <- tau }
  }
  list(tau = best_tau, j = best_j)
}

fit <- youden_max(mixed$p_i, mixed$emit)
cat(sprintf("Pooled mixed sectors: tau* = %.4f, J = %.4f\n\n", fit$tau, fit$j))

# Apply tau to mixed sectors (in-sample sanity check)
m_tpr <- mean(mixed$p_i[mixed$emit == 1] >= fit$tau)
m_fpr <- mean(mixed$p_i[mixed$emit == 0] >= fit$tau)
cat(sprintf("Mixed sectors (in-sample): TPR = %.4f, FPR = %.4f\n\n", m_tpr, m_fpr))

# Apply to non-mixed sectors (true emitters only — there are no non-emitters)
nonmixed <- panel[!(panel$primary_crf_group %in% mixed_crfs) & panel$emit == 1, ]
cat(sprintf("Non-mixed sectors: %d true emitters\n\n", nrow(nonmixed)))

tpr_nonmixed <- mean(nonmixed$p_i >= fit$tau)
cat(sprintf("Applying tau = %.4f to non-mixed: TPR = %.4f (%d / %d)\n\n",
            fit$tau, tpr_nonmixed, sum(nonmixed$p_i >= fit$tau), nrow(nonmixed)))

cat("By CRF group:\n")
for (crf in sort(unique(nonmixed$primary_crf_group))) {
  sub <- nonmixed[nonmixed$primary_crf_group == crf, ]
  cat(sprintf("  %-12s: %4d emitters, TPR = %.3f\n",
              crf, nrow(sub), mean(sub$p_i >= fit$tau)))
}
