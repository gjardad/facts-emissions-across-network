load("C:/Users/jota_/Documents/NBB_data/processed/repeated_cv_proxy_crf_asinh.RData")
panel <- repeated_cv_proxy_panel
panel$p_i <- rowMeans(proxy_matrix > 0)

# Build Dhat using LOSO for mixed CRFs and pooled-mixed for non-mixed
mixed <- c("paper", "refining", "metals")

# LOSO thresholds (from previous run)
tau_loso <- list(paper = 0.980, refining = 0.965, metals = 0.965)
tau_pooled_mixed <- 0.965

panel$tau <- ifelse(panel$primary_crf_group %in% mixed,
                     unlist(tau_loso[as.character(panel$primary_crf_group)]),
                     tau_pooled_mixed)
panel$Dhat <- as.integer(panel$p_i >= panel$tau)

# Branch A: D=1 and Dhat=1
brA <- panel[panel$emit == 1 & panel$Dhat == 1, ]
cat("Branch A total firm-years:", nrow(brA), "\n")
cat("By CRF:\n")
print(table(brA$primary_crf_group))

# Cell-size distribution: cell = (CRF, year) within branch A
cells <- table(brA$primary_crf_group, brA$year)
sizes <- as.vector(cells)
sizes <- sizes[sizes > 0]
cat("\nCell-size (CRF x year, branch A) distribution:\n")
print(summary(sizes))
cat("Total non-empty cells:", length(sizes), "\n")
cat("Cells with 1 emitter:", sum(sizes == 1), "\n")
cat("Cells with <=3:", sum(sizes <= 3), "\n")
cat("Cells with <=5:", sum(sizes <= 5), "\n")

# Show the smallest cells
cat("\nSmallest cells (size <= 3):\n")
small <- which(cells <= 3 & cells > 0, arr.ind = TRUE)
for (k in seq_len(min(20, nrow(small)))) {
  i <- small[k, 1]; j <- small[k, 2]
  cat(sprintf("  %s x %s: %d emitters\n",
              rownames(cells)[i], colnames(cells)[j], cells[i, j]))
}

# Branch B: D=1 and Dhat=0
brB <- panel[panel$emit == 1 & panel$Dhat == 0, ]
cat("\nBranch B (false negatives) total firm-years:", nrow(brB), "\n")
print(table(brB$primary_crf_group))
