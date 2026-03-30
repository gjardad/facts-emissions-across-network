###############################################################################
# Identify firms that enter ETS after 2005 and are currently excluded from
# both training and deployment in their pre-entry years.
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

YEARS <- 2005:2021

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))

accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))

eutl <- firm_year_belgian_euets %>%
  filter(!is.na(emissions), emissions > 0, year %in% YEARS)

# ── For each ETS firm: first and last year with observed emissions ──
ets_span <- eutl %>%
  group_by(vat) %>%
  summarise(
    first_ets_year = min(year),
    last_ets_year  = max(year),
    n_ets_years    = n(),
    avg_emissions  = mean(emissions),
    .groups = "drop"
  )

cat("═══════════════════════════════════════════════════════════════\n")
cat("  LATE ETS ENTRANTS: firms entering after 2005\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Total ETS firms (with positive emissions in 2005-2021):", nrow(ets_span), "\n\n")

# Breakdown by first entry year
entry_tab <- ets_span %>%
  count(first_ets_year) %>%
  mutate(cumulative = cumsum(n))
cat("── First ETS year distribution ──\n")
print(entry_tab, n = 30)

# Late entrants: first year > 2005
late <- ets_span %>% filter(first_ets_year > 2005)
cat(sprintf("\n\nLate entrants (first ETS year > 2005): %d firms\n", nrow(late)))
cat(sprintf("  Of which enter after Phase 1 (>2007): %d\n",
            sum(late$first_ets_year > 2007)))
cat(sprintf("  Of which enter after Phase 2 (>2012): %d\n",
            sum(late$first_ets_year > 2012)))

# For late entrants: how many pre-entry firm-years exist in annual accounts?
late_vats <- late$vat
late_in_accounts <- accounts %>%
  filter(vat %in% late_vats, year %in% YEARS) %>%
  left_join(eutl %>% select(vat, year, emissions), by = c("vat", "year")) %>%
  left_join(late %>% select(vat, first_ets_year), by = "vat") %>%
  mutate(pre_entry = year < first_ets_year)

pre_entry <- late_in_accounts %>% filter(pre_entry)
cat(sprintf("\nPre-entry firm-years in annual accounts: %d\n", nrow(pre_entry)))
cat(sprintf("Unique firms with pre-entry years: %d\n", n_distinct(pre_entry$vat)))

# ── NACE breakdown ──
cat("\n── NACE 2d breakdown of pre-entry firm-years ──\n")
nace_breakdown <- pre_entry %>%
  group_by(nace2d) %>%
  summarise(
    n_firms = n_distinct(vat),
    n_firm_years = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_firms))

# Add avg emissions from ETS years for context
nace_emissions <- late %>%
  left_join(accounts %>% distinct(vat, nace2d), by = "vat") %>%
  group_by(nace2d) %>%
  summarise(avg_ets_emissions = round(mean(avg_emissions)), .groups = "drop")

nace_breakdown <- nace_breakdown %>%
  left_join(nace_emissions, by = "nace2d")

print(nace_breakdown, n = 40)

# ── Detailed: per-firm info for late entrants ──
cat("\n── Late entrants detail (by entry year) ──\n\n")
late_detail <- late %>%
  left_join(accounts %>% distinct(vat, nace2d), by = "vat") %>%
  left_join(
    pre_entry %>%
      group_by(vat) %>%
      summarise(pre_entry_years = paste(sort(unique(year)), collapse = ","),
                n_pre_years = n(), .groups = "drop"),
    by = "vat"
  ) %>%
  filter(!is.na(n_pre_years)) %>%
  select(vat, nace2d, first_ets_year, n_ets_years, avg_emissions,
         n_pre_years, pre_entry_years) %>%
  arrange(first_ets_year, nace2d)

print(late_detail, n = 80)

# ── How many are in B2B? (check if they'd actually get a proxy) ──
cat("\n── Are these firms in the B2B selected sample? ──\n")
load(file.path(PROC_DATA, "deployment_panel.RData"))

# Check if any late-entrant pre-entry years are already in deployment
already_deploy <- deployment_panel %>%
  filter(vat %in% late_vats) %>%
  left_join(late %>% select(vat, first_ets_year), by = "vat") %>%
  filter(year < first_ets_year)

cat(sprintf("Late-entrant pre-entry years already in deployment panel: %d\n",
            nrow(already_deploy)))
cat(sprintf("  (These would exist if deploy_vats used year-specific exclusion)\n"))

# Check: are any late-entrant VATs in the current deployment panel at all?
in_deploy <- deployment_panel %>% filter(vat %in% late_vats)
cat(sprintf("\nLate-entrant firm-years in current deployment panel: %d\n",
            nrow(in_deploy)))

# Summary: total missing firm-years
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat(sprintf("Late ETS entrants: %d firms\n", nrow(late)))
cat(sprintf("Pre-entry firm-years (in accounts, not in ETS): %d\n", nrow(pre_entry)))
cat(sprintf("Currently in deployment: %d of those firm-years\n", nrow(already_deploy)))
cat(sprintf("GAP (should be in deployment but aren't): %d firm-years\n",
            nrow(pre_entry) - nrow(already_deploy)))

cat("\n=== DONE ===\n")
