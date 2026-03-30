###############################################################################
# tests/analyze_flags.R
# Analyze allocation flags from Pareto checkpoint files
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
library(dplyr)

years <- 2008:2011

all_flags <- list()
for (t in years) {
  f <- file.path(PROC_DATA, sprintf("b_loop_pareto_year_%d.RData", t))
  if (!file.exists(f)) { cat("Missing:", f, "\n"); next }
  load(f)
  all_flags[[as.character(t)]] <- year_flags %>% mutate(year = t)
}

flags <- bind_rows(all_flags)

cat("=== OVERALL FLAG DISTRIBUTION (across all draws and years) ===\n")
tab <- flags %>% count(alloc_flag) %>% mutate(pct = round(100 * n / sum(n), 1))
for (i in seq_len(nrow(tab))) {
  cat(sprintf("  %-25s %6d (%5.1f%%)\n", tab$alloc_flag[i], tab$n[i], tab$pct[i]))
}

cat("\n=== PER-YEAR FLAG DISTRIBUTION ===\n")
for (t in years) {
  cat(sprintf("\nYear %d:\n", t))
  tab2 <- flags %>% filter(year == t) %>% count(alloc_flag) %>%
    mutate(pct = round(100 * n / sum(n), 1))
  for (i in seq_len(nrow(tab2))) {
    cat(sprintf("  %-25s %6d (%5.1f%%)\n", tab2$alloc_flag[i], tab2$n[i], tab2$pct[i]))
  }
}

# === MIXED SECTOR ANALYSIS ===
cat("\n=== MIXED SECTOR ANALYSIS ===\n")
mixed <- flags %>% filter(grepl("mixed", alloc_flag))
cat("Total mixed CRF-group x year x draw obs:", nrow(mixed), "\n")

mixed_summary <- mixed %>%
  group_by(crf_group, year) %>%
  summarise(
    n_draws    = n(),
    n_ok       = sum(alloc_flag == "mixed_ok"),
    n_lowered  = sum(alloc_flag == "mixed_lowered"),
    n_capped   = sum(alloc_flag == "mixed_capped"),
    n_no_deploy = sum(alloc_flag == "mixed_no_deploy"),
    pct_capped = round(100 * sum(alloc_flag == "mixed_capped") / n(), 1),
    .groups = "drop"
  )

cat("\nCRF-group x year cells with mixed sectors:", nrow(mixed_summary), "\n")

always_ok <- mixed_summary %>% filter(n_ok == n_draws)
cat("\n1. Threshold NEVER lowered (always mixed_ok):", nrow(always_ok), "\n")

lowered_not_capped <- mixed_summary %>% filter(n_capped == 0, n_lowered > 0)
cat("2. Lowered but NEVER capped:", nrow(lowered_not_capped), "\n")

sometimes_capped <- mixed_summary %>% filter(n_capped > 0, n_capped < n_draws)
cat("3. Capped in SOME draws:", nrow(sometimes_capped), "\n")

always_capped <- mixed_summary %>% filter(n_capped == n_draws)
cat("4. ALWAYS capped (all draws):", nrow(always_capped), "\n")

cat("\n=== DETAIL: Always-capped sector-years ===\n")
if (nrow(always_capped) > 0) {
  print(always_capped %>% select(crf_group, year, n_draws, n_capped), n = 100)
}

cat("\n=== DETAIL: Sometimes-capped sector-years ===\n")
if (nrow(sometimes_capped) > 0) {
  print(sometimes_capped %>%
    select(crf_group, year, n_draws, n_ok, n_lowered, n_capped, pct_capped) %>%
    arrange(desc(pct_capped)), n = 100)
}

cat("\n=== DETAIL: Lowered-not-capped ===\n")
if (nrow(lowered_not_capped) > 0) {
  print(lowered_not_capped %>% select(crf_group, year, n_draws, n_ok, n_lowered), n = 50)
}

cat("\n=== DETAIL: Always-ok ===\n")
if (nrow(always_ok) > 0) {
  print(always_ok %>% select(crf_group, year, n_draws), n = 50)
}

cat("\n=== DONE ===\n")
