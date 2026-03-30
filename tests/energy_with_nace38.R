###############################################################################
# What happens to energy E_deploy if we add NACE 38 deployment firms?
#
# Currently: CRF "energy" → NACE 35 only
# Proposed:  CRF "energy" → NACE 35 + NACE 38 (waste incinerators)
#
# Two effects:
#   1. More deployment firms to absorb E_deploy → lower per-firm burden
#   2. E_ETS might change if some NACE 38 ETS firms move from "commercial" to "energy"
#      (but we saw only 1-2 NACE 38 ETS firms with small emissions)
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

YEARS <- 2005:2021

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
load(file.path(PROC_DATA, "deployment_panel.RData"))
load(file.path(PROC_DATA, "nir_calibration_targets.RData"))
load(file.path(PROC_DATA, "deployment_proxy_list.RData"))

nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE, colClasses = c(nace2d = "character")
) %>% select(nace2d, crf_group)

accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))

eutl <- firm_year_belgian_euets %>%
  filter(!is.na(emissions), emissions > 0, year %in% YEARS)

# ── Current setup: NACE 35 → energy ──────────────────────────────────────────
ets_with_nace <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d))

deploy_nace <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  filter(!is.na(nace2d))

# ── Compute E_deploy under CURRENT and PROPOSED mappings ─────────────────────

compute_energy_stats <- function(energy_nace_codes, label) {
  # ETS emissions attributed to energy
  e_ets <- ets_with_nace %>%
    filter(nace2d %in% energy_nace_codes) %>%
    group_by(year) %>%
    summarise(E_ETS = sum(emissions), n_ets = n(), .groups = "drop")

  # Deployment firms in energy
  proxy_b1 <- proxy_list[[1]]
  n_deploy <- proxy_b1 %>%
    filter(year %in% YEARS) %>%
    left_join(deploy_nace %>% distinct(vat, nace2d), by = "vat") %>%
    filter(nace2d %in% energy_nace_codes) %>%
    mutate(proxy_pos = pmax(sinh(proxy), 0) > 0) %>%
    group_by(year) %>%
    summarise(
      n_deploy_total = n(),
      n_deploy_pos = sum(proxy_pos),
      .groups = "drop"
    )

  # NIR target for energy (unchanged)
  nir_energy <- nir_targets %>%
    filter(crf_group == "energy", year %in% YEARS) %>%
    mutate(E_NIR = E_NIR_kt * 1000) %>%
    select(year, E_NIR)

  panel <- nir_energy %>%
    left_join(e_ets, by = "year") %>%
    left_join(n_deploy, by = "year") %>%
    mutate(
      E_ETS = coalesce(E_ETS, 0),
      n_ets = coalesce(n_ets, 0L),
      n_deploy_total = coalesce(n_deploy_total, 0L),
      n_deploy_pos = coalesce(n_deploy_pos, 0L),
      E_deploy = pmax(E_NIR - E_ETS, 0),
      per_firm = ifelse(n_deploy_pos > 0, E_deploy / n_deploy_pos, NA),
      feasible_25kt = E_deploy <= n_deploy_pos * 25000
    )

  cat(sprintf("\n── %s ──\n", label))
  cat(sprintf("  NACE codes: %s\n", paste(energy_nace_codes, collapse = ", ")))

  for (t in YEARS) {
    row <- panel %>% filter(year == t)
    if (nrow(row) == 0) next
    cat(sprintf("  %d: E_NIR=%6.0fkt E_ETS=%6.0fkt E_deploy=%6.0fkt | n_ets=%3d n_deploy=%4d | per_firm=%6.0ft | feasible=%s\n",
                t,
                row$E_NIR / 1000, row$E_ETS / 1000, row$E_deploy / 1000,
                row$n_ets, row$n_deploy_pos,
                row$per_firm,
                ifelse(row$feasible_25kt, "YES", "NO")))
  }

  panel
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("  ENERGY E_DEPLOY: CURRENT vs WITH NACE 38\n")
cat("═══════════════════════════════════════════════════════════════\n")

current <- compute_energy_stats("35", "CURRENT: NACE 35 only")
proposed <- compute_energy_stats(c("35", "38"), "PROPOSED: NACE 35 + NACE 38")

# ── Also check: what if we add more NACE codes? ─────────────────────────────
# From NID: CHP joint ventures could be at NACE 19, 20, 24, 10, 11
# But we can't identify which firms have CHP
# Let's just see the deployment firm counts for context

cat("\n\n── Deployment firm counts by NACE for context ──\n\n")
deploy_by_nace <- deploy_nace %>%
  filter(year == 2010) %>%  # representative year
  group_by(nace2d) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

cat("Year 2010 deployment firms by NACE 2d (top 20):\n")
print(head(deploy_by_nace, 20), n = 20)

# ── Compare current vs proposed side by side ─────────────────────────────────
cat("\n\n── Side-by-side comparison (2008-2011) ──\n\n")
comp <- current %>%
  filter(year %in% 2008:2011) %>%
  select(year, n_deploy_pos_old = n_deploy_pos, E_deploy_old = E_deploy,
         per_firm_old = per_firm) %>%
  left_join(
    proposed %>%
      filter(year %in% 2008:2011) %>%
      select(year, n_deploy_pos_new = n_deploy_pos, E_deploy_new = E_deploy,
             per_firm_new = per_firm),
    by = "year"
  ) %>%
  mutate(
    added_firms = n_deploy_pos_new - n_deploy_pos_old,
    pct_more_firms = round(100 * added_firms / n_deploy_pos_old, 1),
    per_firm_reduction_pct = round(100 * (1 - per_firm_new / per_firm_old), 1)
  )

for (i in seq_len(nrow(comp))) {
  r <- comp[i, ]
  cat(sprintf("  %d: deploy firms %d → %d (+%d, +%.1f%%) | per_firm %.0fkt → %.0fkt (%.1f%% reduction)\n",
              r$year,
              r$n_deploy_pos_old, r$n_deploy_pos_new, r$added_firms, r$pct_more_firms,
              r$per_firm_old / 1000, r$per_firm_new / 1000, r$per_firm_reduction_pct))
}

# Also check impact on commercial CRF group (which loses NACE 38 firms)
cat("\n\n── Impact on commercial CRF group (loses NACE 38) ──\n\n")
commercial_current <- compute_energy_stats(
  nace_crf$nace2d[nace_crf$crf_group == "commercial"],
  "COMMERCIAL current"
)
commercial_new_naces <- setdiff(
  nace_crf$nace2d[nace_crf$crf_group == "commercial"], "38"
)
commercial_proposed <- compute_energy_stats(
  commercial_new_naces,
  "COMMERCIAL without NACE 38"
)

cat("\n=== DONE ===\n")
