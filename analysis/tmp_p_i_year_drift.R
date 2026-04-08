load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel
panel$p_i <- rowMeans(proxy_matrix > 0)

cat("=== Mean p_i by year (all firms) ===\n")
print(round(tapply(panel$p_i, panel$year, mean), 3))

cat("\n=== Mean p_i by year, true non-emitters only ===\n")
ne <- panel[panel$emit == 0, ]
print(round(tapply(ne$p_i, ne$year, mean), 3))

cat("\n=== Mean p_i by year, true emitters only ===\n")
em <- panel[panel$emit == 1, ]
print(round(tapply(em$p_i, em$year, mean), 3))

cat("\n=== Share with p_i > 0 by year, all firms ===\n")
print(round(tapply(panel$p_i > 0, panel$year, mean), 3))

cat("\n=== Share with p_i > 0 by year, true non-emitters only ===\n")
print(round(tapply(ne$p_i > 0, ne$year, mean), 3))

cat("\n=== Share with p_i > 0 by year, true emitters only ===\n")
print(round(tapply(em$p_i > 0, em$year, mean), 3))

cat("\n=== Same, restricted to mixed CRFs ===\n")
mixed <- panel[panel$primary_crf_group %in% c("paper", "refining", "metals"), ]
mne <- mixed[mixed$emit == 0, ]
mem <- mixed[mixed$emit == 1, ]
cat("Mixed non-emitters, mean p_i by year:\n")
print(round(tapply(mne$p_i, mne$year, mean), 3))
cat("\nMixed emitters, mean p_i by year:\n")
print(round(tapply(mem$p_i, mem$year, mean), 3))
