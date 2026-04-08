load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel

# Compute p_i = mean(proxy > 0) across 200 repeats
p_i <- rowMeans(proxy_matrix > 0)
panel$p_i <- p_i

# ETS-only sectors = CRF groups other than mixed (paper, refining, metals)
mixed <- c("paper", "refining", "metals")
ets_only <- panel[!(panel$primary_crf_group %in% mixed) & panel$emit == 1, ]
cat("True emitters in ETS-only sectors:", nrow(ets_only), "\n\n")

cat("By CRF group:\n")
print(table(ets_only$primary_crf_group))

cat("\np_i distribution among ETS-only true emitters:\n")
print(summary(ets_only$p_i))

cat("\nTPR at various thresholds:\n")
for (tau in c(0.95, 0.9725, 0.98, 0.9825, 0.99)) {
  tpr <- mean(ets_only$p_i >= tau)
  cat(sprintf("  tau = %.4f: TPR = %.4f (%d / %d)\n",
              tau, tpr, sum(ets_only$p_i >= tau), nrow(ets_only)))
}

cat("\nTPR by CRF group at tau = 0.98:\n")
for (crf in sort(unique(ets_only$primary_crf_group))) {
  sub <- ets_only[ets_only$primary_crf_group == crf, ]
  tpr <- mean(sub$p_i >= 0.98)
  cat(sprintf("  %-12s: %d emitters, TPR = %.3f\n", crf, nrow(sub), tpr))
}
