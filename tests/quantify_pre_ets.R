###############################################################################
# Quantify E_pre_ets by CRF group × year
#
# For each firm entering ETS after 2005:
#   1. Average first two observed emissions
#   2. Backcast using sector-year fixed effects: log(e_hat) = log(e_bar) + (mu_st - mu_sT)
#   3. Assign to CRF groups
#   4. Compare E_deploy before and after subtracting E_pre_ets
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

YEARS <- 2005:2021

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
load(file.path(PROC_DATA, "nir_calibration_targets.RData"))

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

# ── Step 1: Identify late entrants and their anchor emissions ────────────────
ets_span <- eutl %>%
  group_by(vat) %>%
  summarise(first_ets_year = min(year), .groups = "drop")

late <- ets_span %>% filter(first_ets_year > 2005)

# Average of first two observed emissions
anchor <- eutl %>%
  inner_join(late, by = "vat") %>%
  group_by(vat) %>%
  arrange(year) %>%
  slice_head(n = 2) %>%
  summarise(
    e_anchor = mean(emissions),
    anchor_year = min(year),  # first year (for sector-year FE reference)
    n_anchor = n(),
    .groups = "drop"
  )

# Attach NACE
anchor <- anchor %>%
  left_join(
    accounts %>%
      group_by(vat) %>%
      summarise(nace2d = names(which.max(table(nace2d))), .groups = "drop"),
    by = "vat"
  ) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group))

cat("Late entrants with anchor + CRF assignment:", nrow(anchor), "\n\n")

# ── Step 2: Compute sector-year means of log emissions ───────────────────────
sector_year_mu <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d)) %>%
  group_by(nace2d, year) %>%
  summarise(mu = mean(log(emissions)), n_firms = n(), .groups = "drop")

cat("Sector-year means computed for", n_distinct(sector_year_mu$nace2d), "sectors\n\n")

# ── Step 3: Backcast emissions for pre-ETS years ─────────────────────────────
# For each anchor firm, for each year < first_ets_year:
#   log(e_hat) = log(e_anchor) + (mu_{s,t} - mu_{s,T})
#   where T = anchor_year, s = nace2d

pre_ets_list <- list()

for (i in seq_len(nrow(anchor))) {
  v <- anchor$vat[i]
  s <- anchor$nace2d[i]
  T_entry <- anchor$anchor_year[i]
  e_anch <- anchor$e_anchor[i]
  crf <- anchor$crf_group[i]

  # Sector-year mean at anchor year
  mu_T <- sector_year_mu %>%
    filter(nace2d == s, year == T_entry) %>%
    pull(mu)

  if (length(mu_T) == 0) next  # no sector-year mean at entry

  # Pre-entry years where firm exists in accounts
  pre_years <- accounts %>%
    filter(vat == v, year < T_entry, year %in% YEARS) %>%
    pull(year) %>%
    unique()

  if (length(pre_years) == 0) next

  for (t in pre_years) {
    mu_t <- sector_year_mu %>%
      filter(nace2d == s, year == t) %>%
      pull(mu)

    if (length(mu_t) == 0) {
      # No sector-year mean for this year — use flat backcast
      e_hat <- e_anch
    } else {
      e_hat <- e_anch * exp(mu_t - mu_T)
    }

    pre_ets_list[[length(pre_ets_list) + 1]] <- data.frame(
      vat = v, year = t, nace2d = s, crf_group = crf,
      emissions_backcast = e_hat,
      e_anchor = e_anch, anchor_year = T_entry,
      stringsAsFactors = FALSE
    )
  }
}

pre_ets <- bind_rows(pre_ets_list)
cat("Pre-ETS backcasted firm-years:", nrow(pre_ets), "\n\n")

# ── Step 4: Aggregate E_pre_ets by CRF group × year ─────────────────────────
E_pre_ets <- pre_ets %>%
  group_by(crf_group, year) %>%
  summarise(
    E_pre_ets = sum(emissions_backcast),
    n_pre_ets = n(),
    .groups = "drop"
  )

# ── Step 5: Compute original E_deploy and adjusted E_deploy ──────────────────
E_ETS_group <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group)) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions), n_ets = n(), .groups = "drop")

panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(E_ETS_group, by = c("crf_group", "year")) %>%
  left_join(E_pre_ets, by = c("crf_group", "year")) %>%
  mutate(
    E_ETS = coalesce(E_ETS, 0),
    E_pre_ets = coalesce(E_pre_ets, 0),
    n_pre_ets = coalesce(n_pre_ets, 0L),
    E_deploy_old = pmax(E_NIR - E_ETS, 0),
    E_deploy_new = pmax(E_NIR - E_ETS - E_pre_ets, 0),
    reduction = E_deploy_old - E_deploy_new,
    pct_reduction = round(100 * reduction / E_deploy_old, 1)
  )

# ── Report ───────────────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════\n")
cat("  E_DEPLOY REDUCTION FROM PRE-ETS BACKCAST\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# By CRF group, averaged across years
cat("── By CRF group (averaged across all years) ──\n\n")
by_crf <- panel %>%
  filter(E_deploy_old > 0) %>%
  group_by(crf_group) %>%
  summarise(
    avg_E_deploy_old = round(mean(E_deploy_old)),
    avg_E_pre_ets = round(mean(E_pre_ets)),
    avg_E_deploy_new = round(mean(E_deploy_new)),
    avg_n_pre_ets = round(mean(n_pre_ets), 1),
    avg_pct_reduction = round(mean(pct_reduction, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_pct_reduction))
print(by_crf, n = 20)

# By CRF group × year (just 2005-2012 where it matters most)
cat("\n── By CRF group × year (2005-2012 only, where pre-ETS firms exist) ──\n\n")
detail <- panel %>%
  filter(year <= 2012, E_pre_ets > 0) %>%
  select(crf_group, year, E_NIR, E_ETS, E_pre_ets, n_pre_ets,
         E_deploy_old, E_deploy_new, pct_reduction) %>%
  mutate(across(c(E_NIR, E_ETS, E_pre_ets, E_deploy_old, E_deploy_new),
                ~round(.x, 0))) %>%
  arrange(crf_group, year)
print(detail, n = 80)

# Summary
cat("\n── Overall summary ──\n\n")
total_old <- sum(panel$E_deploy_old[panel$year <= 2012])
total_pre <- sum(panel$E_pre_ets[panel$year <= 2012])
cat(sprintf("Total E_deploy (2005-2012): %s tonnes\n",
            format(round(total_old), big.mark = ",")))
cat(sprintf("Total E_pre_ets (2005-2012): %s tonnes\n",
            format(round(total_pre), big.mark = ",")))
cat(sprintf("Reduction: %.1f%%\n", 100 * total_pre / total_old))

# Per-firm detail for energy sector
cat("\n── Pre-ETS firms in energy sector (NACE 35) ──\n\n")
energy_pre <- pre_ets %>%
  filter(crf_group == "energy") %>%
  arrange(vat, year)
if (nrow(energy_pre) > 0) {
  for (v in unique(energy_pre$vat)) {
    sub <- energy_pre[energy_pre$vat == v, ]
    cat(sprintf("  vat=%s | nace=%s | anchor=%d (%.0ft) | backcast years: %s\n",
                substr(v, 1, 12), sub$nace2d[1], sub$anchor_year[1],
                sub$e_anchor[1],
                paste(sub$year, sprintf("(%.0ft)", sub$emissions_backcast),
                      collapse = ", ")))
  }
}

cat("\n=== DONE ===\n")
