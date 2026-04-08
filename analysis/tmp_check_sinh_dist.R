load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
pm <- repeated_cv_proxy_panel$proxy_mean
pos <- pm[pm > 0]
cat("Positive proxies:", length(pos), "\n")
print(summary(pos))
s <- sinh(pos)
cat("\nsinh(proxy):\n")
print(summary(s))
cat("\nlog sinh(proxy):\n")
print(summary(log(s)))
cat("\nSD of log sinh(proxy):", sd(log(s)), "\n")
cat("\nQuantiles of sinh(proxy):\n")
print(quantile(s, probs = seq(0, 1, 0.1)))
