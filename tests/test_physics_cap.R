###############################################################################
# What if the upper-bound constraint is 25kt instead of min(ETS)?
# Compare with current min_ets_emit and check feasibility.
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

YEARS <- 2008:2011
PHYSICS_CAP <- 25000  # 25 kt CO2/year

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
  filter(year %in% YEARS, !is.na(emissions), emissions > 0)

# ETS summary by CRF
ets_summary <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group)) %>%
  group_by(crf_group, year) %>%
  summarise(
    E_ETS = sum(emissions),
    n_ets = n(),
    min_ets = min(emissions),
    .groups = "drop"
  )

# Deploy counts (draw 1)
deploy_nace <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  select(vat, nace2d) %>%
  distinct(vat, .keep_all = TRUE) %>%
  left_join(nace_crf, by = "nace2d")

proxy_b1 <- proxy_list[[1]]
deploy_counts <- proxy_b1 %>%
  filter(year %in% YEARS) %>%
  left_join(deploy_nace %>% select(vat, crf_group), by = "vat") %>%
  filter(!is.na(crf_group)) %>%
  mutate(proxy_pos = pmax(sinh(proxy), 0) > 0) %>%
  group_by(crf_group, year) %>%
  summarise(n_deploy_pos = sum(proxy_pos), .groups = "drop")

# E_deploy
E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(ets_summary, by = c("crf_group", "year")) %>%
  left_join(deploy_counts, by = c("crf_group", "year")) %>%
  mutate(
    E_ETS = coalesce(E_ETS, 0),
    E_deploy = pmax(E_NIR - E_ETS, 0)
  ) %>%
  filter(!is.na(n_ets), n_ets > 0, !is.na(n_deploy_pos), n_deploy_pos > 0)

# Compare constraints
comparison <- E_deploy_panel %>%
  mutate(
    cap_min_ets = min_ets,
    cap_physics = PHYSICS_CAP,
    # Effective cap: min of physics and min_ets (we don't want deploy > ETS either)
    cap_effective = pmin(cap_physics, min_ets),
    # Feasibility: can E_deploy fit under each cap?
    capacity_min_ets = n_deploy_pos * min_ets,
    capacity_physics = n_deploy_pos * PHYSICS_CAP,
    capacity_effective = n_deploy_pos * cap_effective,
    feasible_min_ets = E_deploy <= capacity_min_ets,
    feasible_physics = E_deploy <= capacity_physics,
    feasible_effective = E_deploy <= capacity_effective,
    # Uniform allocation per firm
    uniform_per_firm = E_deploy / n_deploy_pos
  )

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PHYSICS-BASED CAP (25kt) vs MIN-ETS CAP\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat(sprintf("Feasible with min_ets cap:   %d / %d\n",
            sum(comparison$feasible_min_ets), nrow(comparison)))
cat(sprintf("Feasible with 25kt cap:      %d / %d\n",
            sum(comparison$feasible_physics), nrow(comparison)))
cat(sprintf("Feasible with min(25kt, min_ets): %d / %d\n",
            sum(comparison$feasible_effective), nrow(comparison)))

cat("\n── Detail by CRF group × year ──\n\n")
detail <- comparison %>%
  select(crf_group, year, E_deploy, n_deploy_pos, min_ets,
         uniform_per_firm, feasible_min_ets, feasible_physics) %>%
  mutate(
    E_deploy = round(E_deploy / 1000),  # kt
    min_ets = round(min_ets),
    uniform_per_firm = round(uniform_per_firm),
    headroom_physics = round(PHYSICS_CAP - uniform_per_firm)
  ) %>%
  arrange(crf_group, year)

print(detail, n = 40)

cat("\n── Summary by CRF group ──\n\n")
summary_crf <- comparison %>%
  group_by(crf_group) %>%
  summarise(
    n_cells = n(),
    avg_min_ets = round(mean(min_ets)),
    avg_uniform = round(mean(uniform_per_firm)),
    feasible_min_ets = sum(feasible_min_ets),
    feasible_25kt = sum(feasible_physics),
    .groups = "drop"
  ) %>%
  mutate(
    min_ets_vs_25kt = ifelse(avg_min_ets > PHYSICS_CAP, "min_ets > 25kt",
                              ifelse(avg_min_ets > 5000, "min_ets > 5kt",
                                     "min_ets < 5kt")),
    uniform_vs_25kt = ifelse(avg_uniform > PHYSICS_CAP, "uniform > 25kt",
                              "uniform < 25kt")
  )
print(summary_crf, n = 20)

cat("\n── Key insight: which sectors flip from infeasible to feasible? ──\n\n")
flipped <- comparison %>%
  filter(!feasible_min_ets & feasible_physics) %>%
  select(crf_group, year, E_deploy, n_deploy_pos, min_ets, uniform_per_firm) %>%
  mutate(across(c(E_deploy, uniform_per_firm), round))
cat(sprintf("%d sector-years flip from infeasible to feasible\n\n", nrow(flipped)))
if (nrow(flipped) > 0) print(flipped, n = 40)

cat("\n── Sectors that remain infeasible even at 25kt ──\n\n")
still_infeasible <- comparison %>%
  filter(!feasible_physics) %>%
  select(crf_group, year, E_deploy, n_deploy_pos, uniform_per_firm) %>%
  mutate(E_deploy_kt = round(E_deploy / 1000),
         uniform_kt = round(uniform_per_firm / 1000, 1)) %>%
  select(crf_group, year, E_deploy_kt, n_deploy_pos, uniform_kt)
cat(sprintf("%d sector-years remain infeasible\n\n", nrow(still_infeasible)))
if (nrow(still_infeasible) > 0) print(still_infeasible, n = 40)

cat("\n=== DONE ===\n")
