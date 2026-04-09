###############################################################################
# analysis/diag_p_zero_distribution.R
#
# PURPOSE
#   Diagnose the distribution of p_zero among imputed firm-years in
#   deployment. Produces a histogram and counts firms in confidence bands.
#
# INPUT
#   {PROC_DATA}/firm_year_pi.RData
#
# OUTPUT
#   {OUTPUT_DIR}/figures/p_zero_histogram.png
#   Console summary
###############################################################################

if (tolower(Sys.info()[["user"]]) == "jardang") {
  REPO_DIR <- "C:/Users/jardang/Documents/facts-emissions-across-network"
} else if (tolower(Sys.info()[["user"]]) == "jota_") {
  REPO_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
                        error = function(e) normalizePath(getwd(), winslash = "/"))
  while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
} else {
  stop("Define REPO_DIR for this user.")
}
source(file.path(REPO_DIR, "paths.R"))

cat("===================================================================\n")
cat("  P_ZERO DISTRIBUTION FOR IMPUTED FIRM-YEARS\n")
cat("===================================================================\n\n")

load(file.path(PROC_DATA, "firm_year_pi.RData"))

imp <- firm_year_pi[firm_year_pi$source == "imputed", ]
cat("Imputed firm-years:", nrow(imp), "\n\n")

cat("p_zero summary:\n")
print(summary(imp$p_zero))

# Confidence bands
bands <- list(
  "p_zero <= 0.05 (almost certainly emitter)" = imp$p_zero <= 0.05,
  "p_zero in (0.05, 0.20] (very likely emitter)" = imp$p_zero > 0.05 & imp$p_zero <= 0.20,
  "p_zero in (0.20, 0.40] (likely emitter)"     = imp$p_zero > 0.20 & imp$p_zero <= 0.40,
  "p_zero in (0.40, 0.60] (uncertain)"          = imp$p_zero > 0.40 & imp$p_zero <= 0.60,
  "p_zero in (0.60, 0.80] (likely non-emitter)" = imp$p_zero > 0.60 & imp$p_zero <= 0.80,
  "p_zero in (0.80, 0.95] (very likely non-emitter)" = imp$p_zero > 0.80 & imp$p_zero <= 0.95,
  "p_zero > 0.95 (almost certainly non-emitter)" = imp$p_zero > 0.95
)

cat("\nConfidence bands:\n")
for (lbl in names(bands)) {
  n <- sum(bands[[lbl]])
  cat(sprintf("  %-55s %7d (%5.1f%%)\n", lbl, n, 100 * n / nrow(imp)))
}

n_confident_emitter <- sum(imp$p_zero <= 0.05)
cat(sprintf("\n=> Confident emitters (p_zero <= 0.05): %d firm-years\n",
            n_confident_emitter))
cat(sprintf("=> Confident emitters (p_zero <= 0.20): %d firm-years\n",
            sum(imp$p_zero <= 0.20)))


# =============================================================================
# Histogram
# =============================================================================

FIG_DIR <- file.path(OUTPUT_DIR, "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

OUT_PNG <- file.path(FIG_DIR, "p_zero_histogram.png")

png(OUT_PNG, width = 1000, height = 700, res = 120)
par(mar = c(5, 5, 4, 2))
hist(imp$p_zero,
     breaks = seq(0, 1, by = 0.025),
     col    = "#4477aa",
     border = "white",
     main   = "Distribution of p_zero across imputed firm-years",
     xlab   = "p_zero (fraction of B=200 draws with zero emissions)",
     ylab   = "Firm-year count",
     cex.axis = 1.0, cex.lab = 1.1, cex.main = 1.2)
abline(v = c(0.05, 0.95), col = "red", lty = 2)
legend("top", legend = c("p_zero = 0.05 (confident emitter)",
                          "p_zero = 0.95 (confident non-emitter)"),
       col = "red", lty = 2, bty = "n", cex = 0.9)
dev.off()

cat("\nHistogram saved to:", OUT_PNG, "\n")
