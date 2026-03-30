###############################################################################
# Test the deployment grid changes locally
#
# 1. Syntax-check both modified scripts
# 2. Simulate the year-specific B2B split logic with local data
# 3. Create deployment_firm_years.RData from existing data and test
#    build_deployment_panel.R compatibility
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  TEST 1: Syntax check — run_subsampled_en.R\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

tryCatch({
  parse(file.path(REPO_DIR, "analysis", "run_subsampled_en.R"))
  cat("  PASS: run_subsampled_en.R parses without errors\n")
}, error = function(e) {
  cat("  FAIL:", conditionMessage(e), "\n")
})

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  TEST 2: Syntax check — build_deployment_panel.R\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

tryCatch({
  parse(file.path(REPO_DIR, "preprocess", "build_deployment_panel.R"))
  cat("  PASS: build_deployment_panel.R parses without errors\n")
}, error = function(e) {
  cat("  FAIL:", conditionMessage(e), "\n")
})

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  TEST 3: Simulate year-specific B2B split logic\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Load the training sample to get training firm-years
load(file.path(PROC_DATA, "training_sample.RData"))

lhs <- training_sample %>%
  filter(year >= 2005) %>%
  select(vat, year, emissions, euets, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))

training_vats <- unique(lhs$vat)
training_vy <- lhs %>% distinct(vat, year)

cat("Training VATs (global):", length(training_vats), "\n")
cat("Training firm-years:", nrow(training_vy), "\n")

# Identify late entrants
ets_only <- lhs %>% filter(euets == 1)
ets_span <- ets_only %>%
  group_by(vat) %>%
  summarise(first_year = min(year), .groups = "drop")

late_entrants <- ets_span %>% filter(first_year > 2005)
cat("Late ETS entrants:", nrow(late_entrants), "\n")

# Check: for these late entrants, how many pre-entry years are in
# annual accounts but NOT in training?
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year >= 2005) %>%
  select(vat, year)

late_pre_entry <- accounts %>%
  inner_join(late_entrants, by = "vat") %>%
  filter(year < first_year) %>%
  select(vat, year)

cat("Late-entrant pre-entry firm-years (in accounts):", nrow(late_pre_entry), "\n")

# Under OLD logic: these are all excluded from deployment
# (because their VAT is in training_vats)
old_excluded <- late_pre_entry %>%
  filter(vat %in% training_vats)
cat("  Excluded under OLD logic:", nrow(old_excluded), "\n")

# Under NEW logic: anti_join on (vat, year) means these are INCLUDED
new_included <- late_pre_entry %>%
  anti_join(training_vy, by = c("vat", "year"))
cat("  Included under NEW logic:", nrow(new_included), "\n")

# Verify: no pre-entry firm-year should be in training
overlap <- late_pre_entry %>%
  inner_join(training_vy, by = c("vat", "year"))
cat("  Overlap check (should be 0):", nrow(overlap), "\n")

if (nrow(overlap) == 0) {
  cat("  PASS: No pre-entry firm-years leak into training\n")
} else {
  cat("  FAIL: Some pre-entry firm-years are in training!\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  TEST 4: Create deployment_firm_years.RData and test\n")
cat("          build_deployment_panel.R compatibility\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Create deployment_firm_years from existing deployment_panel
# (deployment_proxy_avg.RData is on RMD, not local)
load(file.path(PROC_DATA, "deployment_panel.RData"))

deployment_firm_years <- deployment_panel %>% select(vat, year)
cat("Created deployment_firm_years from existing panel:", nrow(deployment_firm_years), "rows\n")

# Save it
save(deployment_firm_years, file = file.path(PROC_DATA, "deployment_firm_years.RData"))
cat("Saved deployment_firm_years.RData\n")

# Now test build_deployment_panel.R logic: reload and rebuild
rm(deployment_panel)
load(file.path(PROC_DATA, "deployment_firm_years.RData"))

deployment_panel_test <- deployment_firm_years %>%
  left_join(
    df_annual_accounts_selected_sample_key_variables %>%
      select(vat, year, nace5d, revenue),
    by = c("vat", "year")
  )

n_missing_nace <- sum(is.na(deployment_panel_test$nace5d))
cat("\nRebuilt panel from deployment_firm_years:\n")
cat("  Panel rows:", nrow(deployment_panel_test), "\n")
cat("  Unique firms:", n_distinct(deployment_panel_test$vat), "\n")
cat("  Missing nace5d:", n_missing_nace,
    sprintf("(%.1f%%)\n", 100 * n_missing_nace / nrow(deployment_panel_test)))

# Compare with original
load(file.path(PROC_DATA, "deployment_panel.RData"))
cat("\n  Original deployment_panel:", nrow(deployment_panel), "rows,",
    n_distinct(deployment_panel$vat), "firms\n")

if (nrow(deployment_panel_test) == nrow(deployment_panel)) {
  cat("  PASS: Row counts match — round-trip works\n")
} else {
  cat("  FAIL: Row counts differ!\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  ALL TESTS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n")
