library(dplyr, warn.conflicts = FALSE)

DATA_DIR <- "C:/Users/jota_/Documents/NBB_data/processed"

all_stats <- list()
for (yr in 2005:2021) {
  f <- file.path(DATA_DIR, sprintf("b_loop_pareto_year_%d.RData", yr))
  if (file.exists(f)) {
    load(f)
    all_stats[[as.character(yr)]] <- year_stats5d
  }
}
stats5d <- bind_rows(all_stats)

cat("Total draw-sector-years (5d):", nrow(stats5d), "\n\n")

# Average across draws within each sector-year
sector_year <- stats5d %>%
  group_by(nace5d, year) %>%
  summarise(
    cp_p90p10    = mean(cp_p90p10, na.rm = TRUE),
    cp_p9010_log = mean(cp_p9010_log, na.rm = TRUE),
    cp_var_log   = mean(cp_var_log, na.rm = TRUE),
    n_draws      = n(),
    .groups      = "drop"
  ) %>%
  mutate(nace2d = substr(nace5d, 1, 2))

# Manufacturing only
mfg <- sector_year %>% filter(nace2d >= "10", nace2d <= "33")

cat("══════════════════════════════════════════════════════════════\n")
cat("  MANUFACTURING — NACE 5-digit\n")
cat("══════════════════════════════════════════════════════════════\n\n")

cat("── Overall manufacturing mean (unweighted across 5d sector-years) ──\n")
cat("  p90/p10 ratio:     ", round(mean(mfg$cp_p90p10, na.rm = TRUE), 1), "\n")
cat("  90-10 log diff:    ", round(mean(mfg$cp_p9010_log, na.rm = TRUE), 2), "\n")
cat("  Var(log cp):       ", round(mean(mfg$cp_var_log, na.rm = TRUE), 2), "\n")

cat("\n── By NACE 2-digit parent (avg of 5d sectors within, across years) ──\n")
by_parent <- mfg %>%
  filter(!is.na(cp_p90p10)) %>%
  group_by(nace2d) %>%
  summarise(
    n_5d_sectors = n_distinct(nace5d),
    n_obs        = n(),
    cp_p90p10    = round(mean(cp_p90p10, na.rm = TRUE), 1),
    cp_log_9010  = round(mean(cp_p9010_log, na.rm = TRUE), 2),
    cp_var_log   = round(mean(cp_var_log, na.rm = TRUE), 2),
    .groups      = "drop"
  ) %>%
  arrange(nace2d)
print(as.data.frame(by_parent), row.names = FALSE)

cat("\n── All sectors mean (unweighted across 5d sector-years) ──\n")
cat("  p90/p10 ratio:     ", round(mean(sector_year$cp_p90p10, na.rm = TRUE), 1), "\n")
cat("  90-10 log diff:    ", round(mean(sector_year$cp_p9010_log, na.rm = TRUE), 2), "\n")

cat("\n── Literature benchmarks ──\n")
cat("  Lyubich et al (2018) — US mfg, 6-digit, plant-level:\n")
cat("    90-10 log diff:  2.27  (p90/p10 ~ 9.7x)\n")
cat("  De Lyon & Dechezlepretre (2025) — 3-4 digit:\n")
cat("    p90/p10 carbon:  24x (all), France 14x at 4-digit\n")
