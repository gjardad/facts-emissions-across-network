###############################################################################
# tests/inspect_pareto_checkpoints.R
# Quick inspection of b_loop_pareto_year checkpoint files.
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
library(dplyr)

years <- 2008:2011

for (t in years) {
  f <- file.path(PROC_DATA, sprintf("b_loop_pareto_year_%d.RData", t))
  if (!file.exists(f)) { cat("Missing:", f, "\n"); next }
  load(f)

  cat(sprintf("\n══════════ Year %d ══════════\n", t))

  # --- Allocation flags ---
  cat("\n  Allocation flags:\n")
  flag_tab <- year_flags %>% count(alloc_flag) %>% mutate(pct = round(100 * n / sum(n), 1))
  for (i in seq_len(nrow(flag_tab))) {
    cat(sprintf("    %-25s %6d (%5.1f%%)\n", flag_tab$alloc_flag[i], flag_tab$n[i], flag_tab$pct[i]))
  }

  # --- Stats 2d: basic check ---
  cat(sprintf("\n  stats2d: %d sector-draw obs, %d unique sectors, %d draws\n",
              nrow(year_stats2d), n_distinct(year_stats2d$nace2d), n_distinct(year_stats2d$draw)))

  # Summary of key dispersion stats across sectors (averaged over draws first)
  summary_2d <- year_stats2d %>%
    group_by(nace2d) %>%
    summarise(
      n_firms_mean = mean(n_firms, na.rm = TRUE),
      s1_gini_mean = mean(s1_gini, na.rm = TRUE),
      s1_p90p10_mean = mean(s1_p90p10, na.rm = TRUE),
      cp_var_log_mean = mean(cp_var_log, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_firms_mean))

  cat("\n  Top 10 sectors by avg n_firms:\n")
  print(head(summary_2d, 10), n = 10)

  # --- Dispersion distribution across sectors ---
  cat(sprintf("\n  Gini distribution across sectors: median=%.3f, p10=%.3f, p90=%.3f\n",
              median(summary_2d$s1_gini_mean, na.rm = TRUE),
              quantile(summary_2d$s1_gini_mean, 0.1, na.rm = TRUE),
              quantile(summary_2d$s1_gini_mean, 0.9, na.rm = TRUE)))

  cat(sprintf("  p90/p10 distribution: median=%.1f, p10=%.1f, p90=%.1f\n",
              median(summary_2d$s1_p90p10_mean, na.rm = TRUE),
              quantile(summary_2d$s1_p90p10_mean, 0.1, na.rm = TRUE),
              quantile(summary_2d$s1_p90p10_mean, 0.9, na.rm = TRUE)))

  # --- Cross-draw variability for a few sectors ---
  big_sectors <- summary_2d$nace2d[1:min(3, nrow(summary_2d))]
  cat("\n  Cross-draw variability (SD of Gini) for top sectors:\n")
  for (sec in big_sectors) {
    ginis <- year_stats2d$s1_gini[year_stats2d$nace2d == sec]
    cat(sprintf("    %s: mean=%.3f, sd=%.4f, n_draws=%d\n",
                sec, mean(ginis, na.rm = TRUE), sd(ginis, na.rm = TRUE), length(ginis)))
  }

  # --- Stats 5d ---
  cat(sprintf("\n  stats5d: %d sector-draw obs, %d unique sectors\n",
              nrow(year_stats5d), n_distinct(year_stats5d$nace5d)))
}

cat("\n\n=== DONE ===\n")
