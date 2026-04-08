load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel
mixed <- panel[panel$primary_crf_group %in% c("paper", "refining", "metals"), ]
emit  <- mixed[mixed$emit == 1, ]

cat("True emitters in mixed CRFs:", nrow(emit), "\n\n")

cat("By CRF x year (true emitters only):\n")
tab <- table(emit$primary_crf_group, emit$year)
print(tab)

cat("\nDistribution of CRF-year cell sizes (true emitters):\n")
sizes <- as.vector(tab)
sizes <- sizes[sizes > 0]
print(summary(sizes))
cat("Cells with 1 emitter:", sum(sizes == 1), "/", length(sizes), "\n")
cat("Cells with <=3 emitters:", sum(sizes <= 3), "/", length(sizes), "\n")
cat("Cells with <=5 emitters:", sum(sizes <= 5), "/", length(sizes), "\n")
