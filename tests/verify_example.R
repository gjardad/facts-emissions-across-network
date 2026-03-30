## Verify the numerical examples in DISTRIBUTION_FITTING.md section 5.3 and 5.5
library(lmom)

xi <- -1.05; alpha <- 1.62; k <- -0.10
gpa_params <- c(xi = xi, alpha = alpha, k = k)

cat("=== CV Example (section 5.3): 6 firms, N=6 ===\n")
ranks_cv <- 1:6
N_cv <- 6
p_cv <- (ranks_cv - 0.5) / N_cv
w_cv <- quagpa(p_cv, gpa_params)
omega_cv <- exp(w_cv - max(w_cv))
shares_cv <- omega_cv / sum(omega_cv)
E_target_cv <- 1000
y_cv <- E_target_cv * shares_cv

cat("\nPlotting positions:\n")
cat(sprintf("  p = %s\n", paste(round(p_cv, 4), collapse = ", ")))
cat("\nGPA quantiles w:\n")
cat(sprintf("  w = %s\n", paste(round(w_cv, 3), collapse = ", ")))
cat("\nUnnorm weights:\n")
cat(sprintf("  omega = %s\n", paste(round(omega_cv, 4), collapse = ", ")))
cat(sprintf("  sum(omega) = %.4f\n", sum(omega_cv)))
cat("\nShares and emissions:\n")
for (i in 1:6) {
  cat(sprintf("  Firm %d: share=%.1f%%, y=%.1f\n", i, 100*shares_cv[i], y_cv[i]))
}
cat(sprintf("  Total: %.1f\n", sum(y_cv)))

cat("\n=== Deployment Example (section 5.5): 5 deploy + 3 ETS, N=8 ===\n")
ranks_dep <- 1:5
N_dep <- 8
p_dep <- (ranks_dep - 0.5) / N_dep
w_dep <- quagpa(p_dep, gpa_params)
omega_dep <- exp(w_dep - max(w_dep))
shares_dep <- omega_dep / sum(omega_dep)
E_deploy <- 200
y_dep <- E_deploy * shares_dep

cat("\nPlotting positions:\n")
cat(sprintf("  p = %s\n", paste(round(p_dep, 4), collapse = ", ")))
cat("\nGPA quantiles w:\n")
cat(sprintf("  w = %s\n", paste(round(w_dep, 3), collapse = ", ")))
cat("\nUnnorm weights:\n")
cat(sprintf("  omega = %s\n", paste(round(omega_dep, 4), collapse = ", ")))
cat(sprintf("  sum(omega) = %.4f\n", sum(omega_dep)))
cat("\nShares and emissions:\n")
for (i in 1:5) {
  cat(sprintf("  d%d: share=%.1f%%, y=%.1f\n", i, 100*shares_dep[i], y_dep[i]))
}
cat(sprintf("  Total: %.1f\n", sum(y_dep)))
cat(sprintf("  Max deploy emission: %.1f vs min_ets_emit=50 → %s\n",
            max(y_dep), if(max(y_dep) >= 50) "VIOLATED" else "OK"))
