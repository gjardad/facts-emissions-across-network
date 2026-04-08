# Temporary script: compare carbon productivity dispersion with literature
# Lyubich et al (2018): mean within-industry 90-10 log diff = 2.27 (US mfg, plant-level, 6-digit)
# De Lyon & Dechezlepretre (2025): mean within-industry p90/p10 = 24x (9 countries, 3-4 digit)
#   European countries: France ~14x, Netherlands ~18x, UK ~16x, Sweden ~12x

library(dplyr, warn.conflicts = FALSE)

DATA_DIR <- "C:/Users/jota_/Documents/NBB_data/processed"

# Load all year files and combine stats2d
all_stats <- list()
for (yr in 2005:2021) {
  f <- file.path(DATA_DIR, sprintf("b_loop_pareto_year_%d.RData", yr))
  if (file.exists(f)) {
    load(f)
    all_stats[[as.character(yr)]] <- year_stats2d
  }
}
stats2d <- bind_rows(all_stats)

cat("Total draw-sector-years:", nrow(stats2d), "\n\n")

# Average across draws within each sector-year
sector_year <- stats2d %>%
  group_by(nace2d, year) %>%
  summarise(
    cp_p90p10    = mean(cp_p90p10, na.rm = TRUE),
    cp_p9010_log = mean(cp_p9010_log, na.rm = TRUE),
    cp_var_log   = mean(cp_var_log, na.rm = TRUE),
    n_draws      = n(),
    .groups      = "drop"
  )

# ── Manufacturing sectors (NACE 10-33) ──────────────────────────────────────
mfg <- sector_year %>% filter(nace2d >= "10", nace2d <= "33")

cat("══════════════════════════════════════════════════════════════\n")
cat("  MANUFACTURING SECTORS (NACE 10-33)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# By sector (avg across years)
by_sector <- mfg %>%
  filter(!is.na(cp_p90p10)) %>%
  group_by(nace2d) %>%
  summarise(
    n_years      = n(),
    cp_p90p10    = round(mean(cp_p90p10, na.rm = TRUE), 1),
    cp_log_9010  = round(mean(cp_p9010_log, na.rm = TRUE), 2),
    cp_var_log   = round(mean(cp_var_log, na.rm = TRUE), 2),
    .groups      = "drop"
  ) %>%
  arrange(nace2d)

cat("By NACE 2-digit (averaged across years):\n")
print(as.data.frame(by_sector), row.names = FALSE)

cat("\n── Overall manufacturing mean (unweighted across sector-years) ──\n")
cat("  p90/p10 ratio:     ", round(mean(mfg$cp_p90p10, na.rm = TRUE), 1), "\n")
cat("  90-10 log diff:    ", round(mean(mfg$cp_p9010_log, na.rm = TRUE), 2), "\n")
cat("  Var(log cp):       ", round(mean(mfg$cp_var_log, na.rm = TRUE), 2), "\n")

cat("\n── Literature benchmarks ──\n")
cat("  Lyubich et al (2018) — US mfg, 6-digit, plant-level:\n")
cat("    90-10 log diff:  2.27  (≈ p90/p10 of 9.7x)\n")
cat("    SD(log cp):      0.89  (≈ Var(log cp) of 0.79)\n")
cat("  De Lyon & Dechezlepretre (2025) — 9 countries, 3-4 digit:\n")
cat("    p90/p10 carbon:  24x (all), France 14x, NL 18x, UK 16x, Sweden 12x\n")

# ── All sectors ─────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat("  ALL SECTORS\n")
cat("══════════════════════════════════════════════════════════════\n\n")

all_sectors <- sector_year %>%
  filter(!is.na(cp_p90p10)) %>%
  group_by(nace2d) %>%
  summarise(
    n_years      = n(),
    cp_p90p10    = round(mean(cp_p90p10, na.rm = TRUE), 1),
    cp_log_9010  = round(mean(cp_p9010_log, na.rm = TRUE), 2),
    .groups      = "drop"
  ) %>%
  arrange(nace2d)

print(as.data.frame(all_sectors), row.names = FALSE)

cat("\n── Overall mean (all sectors, unweighted) ──\n")
cat("  p90/p10 ratio:     ", round(mean(sector_year$cp_p90p10, na.rm = TRUE), 1), "\n")
cat("  90-10 log diff:    ", round(mean(sector_year$cp_p9010_log, na.rm = TRUE), 2), "\n")
