###############################################################################
# Evidence on fuel types for small ETS installations in Belgium
# If the smallest ETS firms (just above 20MW) are gas-fired, it's reasonable
# to assume sub-20MW firms are too.
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))

accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% 2005:2021) %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))

eutl <- firm_year_belgian_euets %>%
  filter(!is.na(emissions), emissions > 0, year %in% 2005:2021)

# What does the emission level of smallest ETS firms tell us about fuel type?
# At 20MW, 8760h, 100% capacity: gas=35.3kt, coal=60kt, oil=48.6kt
# So if a firm emits ~20-35kt, it's likely running at or near capacity on gas
# If it emits ~40-60kt, it could be coal or oil at moderate capacity, or gas at very high capacity

cat("═══════════════════════════════════════════════════════════════\n")
cat("  EMISSION LEVELS OF SMALLEST ETS FIRMS\n")
cat("  (evidence on fuel type at the 20MW boundary)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# For each ETS firm: average emissions across years
ets_avg <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d)) %>%
  group_by(vat, nace2d) %>%
  summarise(avg_emissions = mean(emissions), n_years = n(), .groups = "drop")

# Small ETS firms (avg < 50kt) — these are likely near the 20MW threshold
small_ets <- ets_avg %>%
  filter(avg_emissions < 50000) %>%
  arrange(avg_emissions)

cat("ETS firms with avg emissions < 50kt:", nrow(small_ets), "\n\n")

cat("── Distribution of avg emissions for small ETS firms ──\n")
cat(sprintf("  < 5kt:    %d firms\n", sum(small_ets$avg_emissions < 5000)))
cat(sprintf("  5-10kt:   %d firms\n", sum(small_ets$avg_emissions >= 5000 & small_ets$avg_emissions < 10000)))
cat(sprintf("  10-20kt:  %d firms\n", sum(small_ets$avg_emissions >= 10000 & small_ets$avg_emissions < 20000)))
cat(sprintf("  20-35kt:  %d firms\n", sum(small_ets$avg_emissions >= 20000 & small_ets$avg_emissions < 35000)))
cat(sprintf("  35-50kt:  %d firms\n", sum(small_ets$avg_emissions >= 35000 & small_ets$avg_emissions < 50000)))

# The key insight: at 20MW rated thermal input and natural gas:
# - 100% capacity factor: 35.3 kt/year
# - 70% capacity factor: 24.7 kt/year
# - 50% capacity factor: 17.7 kt/year
# So a firm emitting < 25kt on natural gas would need capacity factor < 70%
# A firm emitting < 35kt is consistent with gas at any reasonable capacity factor

cat("\n── Reference: emissions at 20MW by fuel and capacity factor ──\n")
fuels <- data.frame(
  fuel = c("Natural gas", "Gas oil/diesel", "Heavy fuel oil", "Coal"),
  ef_tCO2_per_TJ = c(56.1, 74.1, 77.4, 95.0)
)
cap_factors <- c(0.40, 0.50, 0.60, 0.70, 0.80)

cat("\n  Emissions (kt CO2/year) at 20MW rated thermal input:\n\n")
cat(sprintf("  %-18s", "Fuel / CF"))
for (cf in cap_factors) cat(sprintf("%8.0f%%", cf * 100))
cat("\n")
for (i in seq_len(nrow(fuels))) {
  cat(sprintf("  %-18s", fuels$fuel[i]))
  for (cf in cap_factors) {
    e <- 20 * cf * 8760 * 3.6 / 1000 * fuels$ef_tCO2_per_TJ[i] / 1000
    cat(sprintf("%8.1f", e))
  }
  cat("\n")
}

# NACE breakdown of small ETS firms
cat("\n── NACE 2d of small ETS firms (< 50kt avg) ──\n")
nace_small <- small_ets %>%
  group_by(nace2d) %>%
  summarise(
    n = n(),
    avg = round(mean(avg_emissions)),
    min = round(min(avg_emissions)),
    max = round(max(avg_emissions)),
    .groups = "drop"
  ) %>%
  arrange(desc(n))
print(nace_small, n = 30)

# What fraction of manufacturing energy is gaseous fuels?
# From NID: 45% in 1990, 70% in 2021, so ~55-65% in 2005-2011
# For commercial: 77% currently
# For sub-20MW: even higher because coal requires bulk handling infrastructure

cat("\n── NID evidence summary ──\n")
cat("  Manufacturing gaseous fuel share: 45% (1990) → 55-65% (2005-2011) → 70% (2021)\n")
cat("  Commercial gaseous fuel share: ~77% (current)\n")
cat("  Coal in manufacturing: concentrated in iron/steel (blast furnaces, coke ovens)\n")
cat("    and cement (kiln fuels) — all large ETS installations\n")
cat("  Coal in power: phased out of Flanders by 2017, concentrated in\n")
cat("    large Electrabel plants during 2005-2011\n")
cat("  Flemish region: since 2017 no coal consumed in power plants\n")

cat("\n=== DONE ===\n")
