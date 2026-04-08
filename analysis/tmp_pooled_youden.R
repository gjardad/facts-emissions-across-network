load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel
panel$p_i <- rowMeans(proxy_matrix > 0)

mixed <- panel[panel$primary_crf_group %in% c("paper", "refining", "metals"), ]

# Pooled Youden's J on raw p_i
emit <- mixed$emit
p    <- mixed$p_i
grid <- sort(unique(p))

best_tau <- NA; best_j <- -Inf
for (tau in grid) {
  pred <- p >= tau
  tpr <- sum(pred & emit == 1) / sum(emit == 1)
  fpr <- sum(pred & emit == 0) / sum(emit == 0)
  j   <- tpr - fpr
  if (j > best_j) { best_j <- j; best_tau <- tau }
}
cat(sprintf("Pooled raw p_i: tau* = %.4f, J = %.4f\n", best_tau, best_j))

# Now year-CRF residualized p_i
mixed$p_resid <- ave(mixed$p_i, mixed$year, FUN = function(x) x - mean(x))
p2 <- mixed$p_resid
grid2 <- sort(unique(p2))
best_tau2 <- NA; best_j2 <- -Inf
for (tau in grid2) {
  pred <- p2 >= tau
  tpr <- sum(pred & emit == 1) / sum(emit == 1)
  fpr <- sum(pred & emit == 0) / sum(emit == 0)
  j   <- tpr - fpr
  if (j > best_j2) { best_j2 <- j; best_tau2 <- tau }
}
cat(sprintf("Pooled CRF-year-residualized p_i: tau* = %.4f, J = %.4f\n", best_tau2, best_j2))
