###############################################################################
# analysis/table_dispersion_comparison.R
#
# PURPOSE
#   Produce comparison tables: ETS-only vs. literature vs. full-sample
#   carbon productivity dispersion.
#
# OUTPUT
#   Console output + CSV files for downstream formatting
###############################################################################

library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)

DATA_DIR  <- "C:/Users/jota_/Documents/NBB_data/processed"
REPO_DIR  <- "C:/Users/jota_/Documents/facts-emissions-across-network"
OUT_DIR   <- file.path(REPO_DIR, "output", "tables")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

MIN_N <- 3L
YEARS <- 2005:2021

# ── Helper ──────────────────────────────────────────────────────────────────
pct_ratio <- function(x, p_hi, p_lo) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  q <- quantile(x, c(p_lo, p_hi), names = FALSE)
  if (q[1] == 0) return(NA_real_)
  q[2] / q[1]
}

do_stats <- function(df) {
  log_cp <- log(df$cp)
  data.frame(
    n_firms     = nrow(df),
    cp_p90p10   = pct_ratio(df$cp, 0.9, 0.1),
    cp_log_9010 = diff(quantile(log_cp, c(0.1, 0.9))),
    cp_var_log  = var(log_cp)
  )
}

# ═══════════════════════════════════════════════════════════════════════════
# PART A: ETS-only dispersion
# ═══════════════════════════════════════════════════════════════════════════
cat("── Loading data ──\n")

load(file.path(DATA_DIR, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS, !is.na(emissions), emissions > 0) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

load(file.path(DATA_DIR, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

ets <- eutl %>%
  inner_join(accounts, by = c("vat", "year")) %>%
  filter(revenue > 0, nace2d >= "10", nace2d <= "33") %>%
  mutate(cp = revenue / emissions)

# ETS 2-digit
ets_2d <- ets %>%
  group_by(nace2d, year) %>%
  filter(n() >= MIN_N) %>%
  group_modify(~do_stats(.x)) %>%
  ungroup() %>%
  group_by(nace2d) %>%
  summarise(
    ets_n_years     = n(),
    ets_n_firms_avg = round(mean(n_firms), 0),
    ets_p90p10      = mean(cp_p90p10, na.rm = TRUE),
    ets_log9010     = mean(cp_log_9010, na.rm = TRUE),
    ets_var_log     = mean(cp_var_log, na.rm = TRUE),
    .groups = "drop"
  )

# ETS 5-digit
ets_5d <- ets %>%
  group_by(nace5d, year) %>%
  filter(n() >= MIN_N) %>%
  group_modify(~do_stats(.x)) %>%
  ungroup() %>%
  mutate(nace2d = substr(nace5d, 1, 2)) %>%
  group_by(nace2d) %>%
  summarise(
    ets5d_n_sectors  = n_distinct(nace5d),
    ets5d_p90p10     = mean(cp_p90p10, na.rm = TRUE),
    ets5d_log9010    = mean(cp_log_9010, na.rm = TRUE),
    ets5d_var_log    = mean(cp_var_log, na.rm = TRUE),
    .groups = "drop"
  )


# ═══════════════════════════════════════════════════════════════════════════
# PART B: Full-sample dispersion (from Pareto loop output)
# ═══════════════════════════════════════════════════════════════════════════
cat("── Loading Pareto loop results ──\n")

all_2d <- list()
all_5d <- list()
for (yr in YEARS) {
  f <- file.path(DATA_DIR, sprintf("b_loop_pareto_year_%d.RData", yr))
  if (file.exists(f)) {
    load(f)
    all_2d[[as.character(yr)]] <- year_stats2d
    all_5d[[as.character(yr)]] <- year_stats5d
  }
}
full_2d_raw <- bind_rows(all_2d)
full_5d_raw <- bind_rows(all_5d)

# Full 2-digit: average across draws, then across years
full_2d <- full_2d_raw %>%
  filter(nace2d >= "10", nace2d <= "33") %>%
  group_by(nace2d, year) %>%
  summarise(
    n_firms     = mean(n_firms),
    cp_p90p10   = mean(cp_p90p10, na.rm = TRUE),
    cp_log_9010 = mean(cp_p9010_log, na.rm = TRUE),
    cp_var_log  = mean(cp_var_log, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  group_by(nace2d) %>%
  summarise(
    full_n_years     = n(),
    full_n_firms_avg = round(mean(n_firms), 0),
    full_p90p10      = mean(cp_p90p10, na.rm = TRUE),
    full_log9010     = mean(cp_log_9010, na.rm = TRUE),
    full_var_log     = mean(cp_var_log, na.rm = TRUE),
    .groups = "drop"
  )

# Full 5-digit: average across draws, then across sector-years, grouped by 2d parent
full_5d <- full_5d_raw %>%
  mutate(nace2d = substr(nace5d, 1, 2)) %>%
  filter(nace2d >= "10", nace2d <= "33") %>%
  group_by(nace5d, year) %>%
  summarise(
    cp_p90p10   = mean(cp_p90p10, na.rm = TRUE),
    cp_log_9010 = mean(cp_p9010_log, na.rm = TRUE),
    cp_var_log  = mean(cp_var_log, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(nace2d = substr(nace5d, 1, 2)) %>%
  group_by(nace2d) %>%
  summarise(
    full5d_n_sectors = n_distinct(nace5d),
    full5d_p90p10    = mean(cp_p90p10, na.rm = TRUE),
    full5d_log9010   = mean(cp_log_9010, na.rm = TRUE),
    full5d_var_log   = mean(cp_var_log, na.rm = TRUE),
    .groups = "drop"
  )


# ═══════════════════════════════════════════════════════════════════════════
# TABLE 1: Overall summary — Belgium vs. Literature
# ═══════════════════════════════════════════════════════════════════════════

# Compute unweighted means across sectors
ets_2d_overall  <- c(p90p10 = mean(ets_2d$ets_p90p10, na.rm = TRUE),
                     log9010 = mean(ets_2d$ets_log9010, na.rm = TRUE),
                     var_log = mean(ets_2d$ets_var_log, na.rm = TRUE))
ets_5d_overall  <- c(p90p10 = mean(ets_5d$ets5d_p90p10, na.rm = TRUE),
                     log9010 = mean(ets_5d$ets5d_log9010, na.rm = TRUE),
                     var_log = mean(ets_5d$ets5d_var_log, na.rm = TRUE))
full_2d_overall <- c(p90p10 = mean(full_2d$full_p90p10, na.rm = TRUE),
                     log9010 = mean(full_2d$full_log9010, na.rm = TRUE),
                     var_log = mean(full_2d$full_var_log, na.rm = TRUE))
full_5d_overall <- c(p90p10 = mean(full_5d$full5d_p90p10, na.rm = TRUE),
                     log9010 = mean(full_5d$full5d_log9010, na.rm = TRUE),
                     var_log = mean(full_5d$full5d_var_log, na.rm = TRUE))

tab1 <- data.frame(
  Source = c(
    "Belgium ETS-only (NACE 2-digit)",
    "Belgium ETS-only (NACE 5-digit)",
    "",
    "Lyubich et al. (2018)",
    "De Lyon & Dechezlepretre (2025)",
    "",
    "Belgium all firms (NACE 2-digit)",
    "Belgium all firms (NACE 5-digit)"
  ),
  Granularity = c("2-digit", "5-digit", "",
                  "6-digit NAICS", "3-4 digit NACE", "",
                  "2-digit", "5-digit"),
  Unit = c("firm", "firm", "", "plant", "firm (20+ emp.)", "", "firm", "firm"),
  p90_p10 = c(
    round(ets_2d_overall["p90p10"], 1),
    round(ets_5d_overall["p90p10"], 1),
    NA,
    9.7,   # exp(2.27)
    24,    # all-country average
    NA,
    round(full_2d_overall["p90p10"], 1),
    round(full_5d_overall["p90p10"], 1)
  ),
  log_9010 = c(
    round(ets_2d_overall["log9010"], 2),
    round(ets_5d_overall["log9010"], 2),
    NA,
    2.27,
    3.18,  # log(24)
    NA,
    round(full_2d_overall["log9010"], 2),
    round(full_5d_overall["log9010"], 2)
  ),
  stringsAsFactors = FALSE
)

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  TABLE 1: Carbon Productivity Dispersion — Belgium vs. Literature\n")
cat("  Manufacturing sectors only. Mean across sector-years.\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
print(tab1, row.names = FALSE, right = FALSE, na.print = "")
write.csv(tab1, file.path(OUT_DIR, "dispersion_comparison_overview.csv"),
          row.names = FALSE, na = "")


# ═══════════════════════════════════════════════════════════════════════════
# TABLE 2: Sector-by-sector — ETS-only vs. Full sample (NACE 2-digit)
# ═══════════════════════════════════════════════════════════════════════════

# NACE 2-digit labels
nace2d_labels <- c(
  "10" = "Food", "11" = "Beverages", "12" = "Tobacco",
  "13" = "Textiles", "14" = "Wearing apparel", "15" = "Leather",
  "16" = "Wood", "17/18" = "Paper & printing", "19" = "Coke & petroleum",
  "20" = "Chemicals", "21" = "Pharmaceuticals", "22" = "Rubber & plastics",
  "23" = "Non-metallic minerals", "24" = "Basic metals",
  "25" = "Fabricated metals", "26" = "Electronics",
  "27" = "Electrical equipment", "28" = "Machinery",
  "29" = "Motor vehicles", "30" = "Other transport",
  "31" = "Furniture", "32" = "Other manufacturing",
  "33" = "Repair & installation"
)

tab2 <- full_2d %>%
  left_join(ets_2d, by = "nace2d") %>%
  mutate(
    label = nace2d_labels[nace2d],
    inflation_p90p10 = full_p90p10 / ets_p90p10,
    inflation_log    = full_log9010 / ets_log9010
  ) %>%
  transmute(
    NACE     = nace2d,
    Sector   = label,
    ETS_n    = ets_n_firms_avg,
    ETS_yrs  = ets_n_years,
    Full_n   = full_n_firms_avg,
    ETS_p90p10  = round(ets_p90p10, 1),
    Full_p90p10 = round(full_p90p10, 1),
    ETS_log9010  = round(ets_log9010, 2),
    Full_log9010 = round(full_log9010, 2),
    Ratio_p90p10 = round(inflation_p90p10, 1)
  ) %>%
  arrange(NACE)

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  TABLE 2: ETS-only vs. Full Sample — NACE 2-digit\n")
cat("  Carbon productivity p90/p10 and 90-10 log diff.\n")
cat("  Averaged across years. ETS_n = avg firms/year in ETS.\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

print(tab2, row.names = FALSE, right = FALSE, na.print = "—")
write.csv(tab2, file.path(OUT_DIR, "dispersion_ets_vs_full_2d.csv"),
          row.names = FALSE, na = "")


# ═══════════════════════════════════════════════════════════════════════════
# TABLE 3: Same comparison at NACE 5-digit (grouped by 2-digit parent)
# ═══════════════════════════════════════════════════════════════════════════

tab3 <- full_5d %>%
  left_join(ets_5d, by = "nace2d") %>%
  mutate(
    label = nace2d_labels[nace2d],
    inflation_log = full5d_log9010 / ets5d_log9010
  ) %>%
  transmute(
    NACE     = nace2d,
    Sector   = label,
    ETS_5d   = ets5d_n_sectors,
    Full_5d  = full5d_n_sectors,
    ETS_log9010  = round(ets5d_log9010, 2),
    Full_log9010 = round(full5d_log9010, 2),
    Ratio_log    = round(inflation_log, 1)
  ) %>%
  arrange(NACE)

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  TABLE 3: ETS-only vs. Full Sample — NACE 5-digit\n")
cat("  Mean of within-5d-sector 90-10 log diff, by 2-digit parent.\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

print(tab3, row.names = FALSE, right = FALSE, na.print = "—")
write.csv(tab3, file.path(OUT_DIR, "dispersion_ets_vs_full_5d.csv"),
          row.names = FALSE, na = "")

cat("\n── Done. CSVs saved to", OUT_DIR, "──\n")
