###############################################################################
# tests/analyze_cap_magnitude.R
#
# Quantify the capping problem:
# 1. For each CRF-group x year: E_deploy, n_deploy, min_ets_emit,
#    and whether uniform allocation would even fit
# 2. Investigate min-ETS emitters: are they stable firms or one-off low years?
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

YEARS <- 2008:2011

# ── Load data ────────────────────────────────────────────────────────────────
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
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)

# ── ETS stats per CRF group x year ──────────────────────────────────────────
ets_crf <- eutl %>%
  filter(!is.na(emissions), emissions > 0) %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group))

ets_summary <- ets_crf %>%
  group_by(crf_group, year) %>%
  summarise(
    n_ets          = n(),
    E_ETS          = sum(emissions),
    min_ets_emit   = min(emissions),
    p05_ets_emit   = quantile(emissions, 0.05),
    p10_ets_emit   = quantile(emissions, 0.10),
    p25_ets_emit   = quantile(emissions, 0.25),
    median_ets_emit = median(emissions),
    .groups = "drop"
  )

# ── E_deploy per CRF group x year ───────────────────────────────────────────
E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(ets_summary %>% select(crf_group, year, E_ETS), by = c("crf_group", "year")) %>%
  mutate(E_ETS = coalesce(E_ETS, 0),
         E_deploy = pmax(E_NIR - E_ETS, 0))

# ── Deployment firm counts per CRF group x year (using draw 1) ──────────────
deploy_nace <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  select(vat, nace2d) %>%
  distinct(vat, .keep_all = TRUE) %>%
  left_join(nace_crf, by = "nace2d")

# Count deployment firms with positive proxy per CRF group (draw 1)
proxy_b1 <- proxy_list[[1]]
deploy_counts <- proxy_b1 %>%
  filter(year %in% YEARS) %>%
  left_join(deploy_nace %>% select(vat, crf_group), by = "vat") %>%
  filter(!is.na(crf_group)) %>%
  mutate(proxy_pos = pmax(sinh(proxy), 0) > 0) %>%
  group_by(crf_group, year) %>%
  summarise(
    n_deploy_total = n(),
    n_deploy_pos   = sum(proxy_pos),
    .groups = "drop"
  )

# ── Combine ──────────────────────────────────────────────────────────────────
combined <- E_deploy_panel %>%
  select(crf_group, year, E_NIR, E_ETS, E_deploy) %>%
  left_join(ets_summary, by = c("crf_group", "year")) %>%
  left_join(deploy_counts, by = c("crf_group", "year")) %>%
  filter(!is.na(n_ets), n_ets > 0, !is.na(n_deploy_pos), n_deploy_pos > 0) %>%
  mutate(
    # Maximum capacity: all deploy firms at min_ets_emit
    max_capacity     = n_deploy_pos * min_ets_emit,
    # Is E_deploy > max_capacity? (mathematically impossible even with uniform)
    impossible       = E_deploy > max_capacity,
    # Ratio: how many times over capacity?
    overcommit_ratio = E_deploy / max_capacity,
    # Uniform allocation per firm
    uniform_per_firm = E_deploy / n_deploy_pos,
    # Ratio of uniform to min_ets_emit
    uniform_vs_min   = uniform_per_firm / min_ets_emit
  )

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CAPPING MAGNITUDE ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("── Mathematically impossible cases (E_deploy > n_deploy * min_ets) ──\n")
cat("   (even uniform allocation would violate the constraint)\n\n")
imp <- combined %>% filter(impossible)
cat(sprintf("  %d of %d mixed sector-years are IMPOSSIBLE\n\n", nrow(imp), nrow(combined)))

if (nrow(imp) > 0) {
  imp_print <- imp %>%
    select(crf_group, year, E_deploy, n_deploy_pos, min_ets_emit,
           max_capacity, overcommit_ratio, uniform_per_firm) %>%
    mutate(across(c(E_deploy, min_ets_emit, max_capacity, uniform_per_firm),
                  ~round(.x, 1)),
           overcommit_ratio = round(overcommit_ratio, 2)) %>%
    arrange(desc(overcommit_ratio))
  print(imp_print, n = 50)
}

cat("\n── Cases where uniform would fit but Pareto concentrates too much ──\n")
feas <- combined %>% filter(!impossible)
cat(sprintf("  %d of %d mixed sector-years are feasible under uniform\n\n", nrow(feas), nrow(combined)))
if (nrow(feas) > 0) {
  feas_print <- feas %>%
    select(crf_group, year, E_deploy, n_deploy_pos, min_ets_emit,
           max_capacity, overcommit_ratio, uniform_per_firm) %>%
    mutate(across(c(E_deploy, min_ets_emit, max_capacity, uniform_per_firm),
                  ~round(.x, 1)),
           overcommit_ratio = round(overcommit_ratio, 2)) %>%
    arrange(desc(overcommit_ratio))
  print(feas_print, n = 50)
}

cat("\n── Summary by CRF group (averaged across years) ──\n\n")
group_avg <- combined %>%
  group_by(crf_group) %>%
  summarise(
    avg_E_deploy       = mean(E_deploy),
    avg_n_ets          = mean(n_ets),
    avg_n_deploy       = mean(n_deploy_pos),
    avg_min_ets        = mean(min_ets_emit),
    avg_overcommit     = mean(overcommit_ratio),
    pct_impossible     = round(100 * mean(impossible), 0),
    .groups = "drop"
  ) %>%
  mutate(across(c(avg_E_deploy, avg_min_ets), ~round(.x, 1)),
         avg_overcommit = round(avg_overcommit, 2)) %>%
  arrange(desc(avg_overcommit))
print(group_avg, n = 20)

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: Investigate min-ETS emitters
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  MIN-ETS EMITTER INVESTIGATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# For each CRF group x year, who is the min emitter?
min_emitters <- ets_crf %>%
  group_by(crf_group, year) %>%
  filter(emissions == min(emissions)) %>%
  ungroup()

cat("── Min emitters: how stable are they over time? ──\n")

# For each min-emitter firm, look at their full time series
min_vats <- unique(min_emitters$vat)
min_histories <- ets_crf %>%
  filter(vat %in% min_vats) %>%
  arrange(vat, year)

# For each min-emitter occurrence, compute ratio to their own median
min_detail <- min_emitters %>%
  left_join(
    ets_crf %>%
      group_by(vat, crf_group) %>%
      summarise(
        n_years     = n(),
        own_median  = median(emissions),
        own_mean    = mean(emissions),
        own_min     = min(emissions),
        own_max     = max(emissions),
        .groups = "drop"
      ),
    by = c("vat", "crf_group")
  ) %>%
  mutate(
    ratio_to_own_median = round(emissions / own_median, 3),
    own_cv              = round(sd(emissions) / own_mean, 3)
  ) %>%
  select(crf_group, year, vat, emissions, own_median, own_mean,
         own_min, own_max, n_years, ratio_to_own_median) %>%
  arrange(crf_group, year)

cat("\nAll min-emitter occurrences with their own emission history:\n")
print(min_detail, n = 60)

# Flag potential outlier years: emission < 20% of own median
cat("\n── Potential outlier min-emitters (emission < 20% of own median) ──\n")
outliers <- min_detail %>% filter(ratio_to_own_median < 0.20)
if (nrow(outliers) > 0) {
  print(outliers, n = 30)
} else {
  cat("  None found.\n")
}

cat("\n── Potential outlier min-emitters (emission < 50% of own median) ──\n")
outliers50 <- min_detail %>% filter(ratio_to_own_median < 0.50)
if (nrow(outliers50) > 0) {
  print(outliers50, n = 30)
} else {
  cat("  None found.\n")
}

# What if we used p05 or p10 of ETS instead of min?
cat("\n\n── Effect of using percentile thresholds instead of min ──\n")
cat("   (how much headroom do we gain?)\n\n")

compare_thresholds <- combined %>%
  mutate(
    cap_at_min = min_ets_emit,
    cap_at_p05 = p05_ets_emit,
    cap_at_p10 = p10_ets_emit,
    cap_at_p25 = p25_ets_emit,
    capacity_min = n_deploy_pos * min_ets_emit,
    capacity_p05 = n_deploy_pos * p05_ets_emit,
    capacity_p10 = n_deploy_pos * p10_ets_emit,
    capacity_p25 = n_deploy_pos * p25_ets_emit,
    feasible_min = E_deploy <= capacity_min,
    feasible_p05 = E_deploy <= capacity_p05,
    feasible_p10 = E_deploy <= capacity_p10,
    feasible_p25 = E_deploy <= capacity_p25
  )

cat(sprintf("  Feasible with cap = min:  %d / %d\n",
            sum(compare_thresholds$feasible_min), nrow(compare_thresholds)))
cat(sprintf("  Feasible with cap = p05:  %d / %d\n",
            sum(compare_thresholds$feasible_p05), nrow(compare_thresholds)))
cat(sprintf("  Feasible with cap = p10:  %d / %d\n",
            sum(compare_thresholds$feasible_p10), nrow(compare_thresholds)))
cat(sprintf("  Feasible with cap = p25:  %d / %d\n",
            sum(compare_thresholds$feasible_p25), nrow(compare_thresholds)))

# ── NACE 2d level analysis ───────────────────────────────────────────────────
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  NACE 2-DIGIT LEVEL ANALYSIS (for approach 1)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# If we set the constraint at NACE 2d level instead of CRF level
ets_nace2d <- ets_crf %>%
  group_by(nace2d, year) %>%
  summarise(
    n_ets        = n(),
    min_ets_emit = min(emissions),
    .groups = "drop"
  )

deploy_nace2d <- proxy_b1 %>%
  filter(year %in% YEARS) %>%
  left_join(deploy_nace %>% select(vat, nace2d, crf_group), by = "vat") %>%
  filter(!is.na(nace2d), !is.na(crf_group)) %>%
  mutate(proxy_pos = pmax(sinh(proxy), 0) > 0) %>%
  group_by(nace2d, crf_group, year) %>%
  summarise(
    n_deploy_pos = sum(proxy_pos),
    .groups = "drop"
  )

# E_deploy is still at CRF level, but we'd need to apportion it to NACE 2d
# For now: just check how many NACE 2d sectors within each CRF group have ETS firms
ets_nace2d_by_crf <- ets_crf %>%
  distinct(crf_group, nace2d, year) %>%
  group_by(crf_group, year) %>%
  summarise(n_nace2d_with_ets = n(),
            nace2d_list = paste(sort(unique(nace2d)), collapse = ", "),
            .groups = "drop")

cat("── NACE 2d sectors with ETS firms, by CRF group ──\n")
for (cg in sort(unique(ets_nace2d_by_crf$crf_group))) {
  sub <- ets_nace2d_by_crf %>% filter(crf_group == cg) %>%
    distinct(nace2d_list) %>% pull(nace2d_list)
  cat(sprintf("  %-12s: %s\n", cg, sub[1]))
}

# For each NACE2d x year with ETS firms: min_ets_emit at NACE2d level
# vs min_ets_emit at CRF level
cat("\n── Min ETS at NACE 2d vs CRF level ──\n")
nace2d_vs_crf <- ets_nace2d %>%
  left_join(nace_crf, by = "nace2d") %>%
  left_join(ets_summary %>% select(crf_group, year, min_ets_emit_crf = min_ets_emit),
            by = c("crf_group", "year")) %>%
  mutate(ratio = round(min_ets_emit / min_ets_emit_crf, 2)) %>%
  filter(year == 2010) %>%   # just show one year
  select(crf_group, nace2d, min_ets_emit, min_ets_emit_crf, n_ets, ratio) %>%
  arrange(crf_group, nace2d)

cat("\n  (showing year 2010 as example)\n")
print(nace2d_vs_crf, n = 40)

cat("\n=== DONE ===\n")
