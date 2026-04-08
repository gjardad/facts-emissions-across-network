###############################################################################
# analysis/characterize_residuals.R
#
# PURPOSE
#   Diagnostic exploration of OOS prediction residuals from run_subsampled_en.R.
#   Characterizes the empirical distribution of epsilon = asinh(y_true) - proxy
#   for held-out ETS emitters across B draws.
#
#   Key questions:
#     1. Shape of the residual distribution (symmetric? heavy tails?)
#     2. Heteroscedasticity: does Var(epsilon) depend on proxy level?
#     3. Sector dependence: do residuals differ across NACE 2-digit sectors?
#     4. Draw dependence: are some draws noisier than others?
#
# INPUT
#   {PROC_DATA}/oos_residuals.RData
#
# RUNS ON: local 1 (after transfer from RMD)
###############################################################################

# -- Paths --------------------------------------------------------------------
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

library(dplyr)


# =============================================================================
# SECTION 1: Load residuals
# =============================================================================
cat("===================================================================\n")
cat("  OOS RESIDUAL CHARACTERIZATION\n")
cat("===================================================================\n\n")

load(file.path(PROC_DATA, "oos_residuals.RData"))
cat("Loaded:", nrow(oos_residuals), "OOS residuals across",
    n_distinct(oos_residuals$draw), "draws\n")
cat("  Unique firms:", n_distinct(oos_residuals$vat), "\n")
cat("  Years:", paste(range(oos_residuals$year), collapse = "--"), "\n\n")


# =============================================================================
# SECTION 2: Overall distribution
# =============================================================================
cat("-- Overall distribution ---------------------------------------------\n")
cat("  Mean epsilon:   ", round(mean(oos_residuals$epsilon), 4), "\n")
cat("  Median epsilon: ", round(median(oos_residuals$epsilon), 4), "\n")
cat("  SD epsilon:     ", round(sd(oos_residuals$epsilon), 4), "\n")
cat("  Skewness:       ", round(moments::skewness(oos_residuals$epsilon), 4), "\n")
cat("  Kurtosis:       ", round(moments::kurtosis(oos_residuals$epsilon), 4), "\n")
cat("  Range:          [", round(min(oos_residuals$epsilon), 2), ",",
    round(max(oos_residuals$epsilon), 2), "]\n")

# Percentiles
pcts <- quantile(oos_residuals$epsilon, probs = c(0.01, 0.05, 0.10, 0.25, 0.50,
                                                     0.75, 0.90, 0.95, 0.99))
cat("  Percentiles:\n")
for (i in seq_along(pcts)) {
  cat(sprintf("    %5s: %7.3f\n", names(pcts)[i], pcts[i]))
}

# Spearman correlation between proxy and true emissions
spearman <- cor(oos_residuals$proxy_oos, oos_residuals$y_asinh,
                method = "spearman", use = "complete.obs")
pearson  <- cor(oos_residuals$proxy_oos, oos_residuals$y_asinh,
                method = "pearson", use = "complete.obs")
cat("\n  Spearman(proxy, y_asinh):", round(spearman, 4), "\n")
cat("  Pearson(proxy, y_asinh): ", round(pearson, 4), "\n\n")


# =============================================================================
# SECTION 3: Heteroscedasticity — Var(epsilon) vs. proxy level
# =============================================================================
cat("-- Heteroscedasticity: SD(epsilon) by proxy decile -----------------\n")

oos_residuals$proxy_decile <- cut(
  oos_residuals$proxy_oos,
  breaks = quantile(oos_residuals$proxy_oos, probs = seq(0, 1, 0.1),
                     na.rm = TRUE),
  include.lowest = TRUE,
  labels = paste0("D", 1:10)
)

het_table <- oos_residuals %>%
  filter(!is.na(proxy_decile)) %>%
  group_by(proxy_decile) %>%
  summarise(
    n         = n(),
    mean_eps  = mean(epsilon),
    sd_eps    = sd(epsilon),
    mean_proxy = mean(proxy_oos),
    .groups   = "drop"
  )

for (i in seq_len(nrow(het_table))) {
  cat(sprintf("  %3s (mean proxy=%6.2f): SD=%5.3f, mean=%6.3f (n=%d)\n",
              het_table$proxy_decile[i], het_table$mean_proxy[i],
              het_table$sd_eps[i], het_table$mean_eps[i], het_table$n[i]))
}


# =============================================================================
# SECTION 4: Sector dependence
# =============================================================================
cat("\n-- Sector dependence: SD(epsilon) by NACE 2-digit ------------------\n")

sec_table <- oos_residuals %>%
  group_by(nace2d) %>%
  summarise(
    n        = n(),
    sd_eps   = sd(epsilon),
    mean_eps = mean(epsilon),
    .groups  = "drop"
  ) %>%
  arrange(desc(n))

for (i in seq_len(min(nrow(sec_table), 20))) {
  cat(sprintf("  NACE %2s: SD=%5.3f, mean=%6.3f (n=%d)\n",
              sec_table$nace2d[i], sec_table$sd_eps[i],
              sec_table$mean_eps[i], sec_table$n[i]))
}
if (nrow(sec_table) > 20)
  cat(sprintf("  ... and %d more sectors\n", nrow(sec_table) - 20))


# =============================================================================
# SECTION 5: Draw dependence
# =============================================================================
cat("\n-- Draw dependence: SD(epsilon) across draws ------------------------\n")

draw_table <- oos_residuals %>%
  group_by(draw) %>%
  summarise(
    n        = n(),
    sd_eps   = sd(epsilon),
    mean_eps = mean(epsilon),
    .groups  = "drop"
  )

cat("  Across", nrow(draw_table), "draws:\n")
cat("    SD(epsilon) — mean:", round(mean(draw_table$sd_eps), 4),
    " sd:", round(sd(draw_table$sd_eps), 4),
    " range: [", round(min(draw_table$sd_eps), 4), ",",
    round(max(draw_table$sd_eps), 4), "]\n")
cat("    Mean(epsilon) — mean:", round(mean(draw_table$mean_eps), 4),
    " sd:", round(sd(draw_table$mean_eps), 4), "\n")
cat("    N per draw — mean:", round(mean(draw_table$n), 0),
    " range: [", min(draw_table$n), ",", max(draw_table$n), "]\n")


# =============================================================================
# SECTION 6: Summary assessment
# =============================================================================
cat("\n===================================================================\n")
cat("  SUMMARY\n")
cat("===================================================================\n")
cat("Total residuals:", nrow(oos_residuals), "\n")
cat("Overall SD:", round(sd(oos_residuals$epsilon), 4), "\n")
cat("Spearman:", round(spearman, 4), "\n")

# Test for heteroscedasticity: correlation between |epsilon| and proxy
het_cor <- cor(abs(oos_residuals$epsilon), oos_residuals$proxy_oos,
               method = "spearman", use = "complete.obs")
cat("Spearman(|epsilon|, proxy):", round(het_cor, 4),
    if (abs(het_cor) > 0.1) " ** heteroscedastic **" else " (roughly homoscedastic)",
    "\n")

cat("\nFor pooled sampling in prediction intervals, use the full epsilon vector.\n")
cat("===================================================================\n")
