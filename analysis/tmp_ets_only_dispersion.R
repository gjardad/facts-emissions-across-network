library(dplyr, warn.conflicts = FALSE)

DATA_DIR  <- "C:/Users/jota_/Documents/NBB_data/processed"
REPO_DIR  <- "C:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

# Load ETS emissions
load(file.path(DATA_DIR, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year >= 2005, year <= 2021, !is.na(emissions), emissions > 0) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

# Load accounts for NACE codes + revenue
load(file.path(DATA_DIR, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year >= 2005, year <= 2021) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

# Merge
ets <- eutl %>%
  inner_join(accounts, by = c("vat", "year")) %>%
  filter(revenue > 0)

cat("ETS firm-years with emissions > 0 and revenue > 0:", nrow(ets), "\n")
cat("Unique firms:", n_distinct(ets$vat), "\n\n")

# Carbon productivity
ets <- ets %>% mutate(cp = revenue / emissions)

# Helper functions
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

MIN_N <- 3L

# ── NACE 2-digit ────────────────────────────────────────────────────────────
stats2d <- ets %>%
  group_by(nace2d, year) %>%
  filter(n() >= MIN_N) %>%
  group_modify(~do_stats(.x)) %>%
  ungroup()

mfg2d <- stats2d %>% filter(nace2d >= "10", nace2d <= "33")

cat("══════════════════════════════════════════════════════════════\n")
cat("  ETS ONLY — NACE 2-digit\n")
cat("══════════════════════════════════════════════════════════════\n\n")

cat("── Manufacturing mean (unweighted across sector-years) ──\n")
cat("  p90/p10 ratio:     ", round(mean(mfg2d$cp_p90p10, na.rm = TRUE), 1), "\n")
cat("  90-10 log diff:    ", round(mean(mfg2d$cp_log_9010, na.rm = TRUE), 2), "\n")
cat("  Var(log cp):       ", round(mean(mfg2d$cp_var_log, na.rm = TRUE), 2), "\n")

cat("\n── By NACE 2-digit (avg across years) ──\n")
by2d <- mfg2d %>%
  group_by(nace2d) %>%
  summarise(
    n_years     = n(),
    n_firms_avg = round(mean(n_firms), 0),
    cp_p90p10   = round(mean(cp_p90p10, na.rm = TRUE), 1),
    cp_log_9010 = round(mean(cp_log_9010, na.rm = TRUE), 2),
    cp_var_log  = round(mean(cp_var_log, na.rm = TRUE), 2),
    .groups     = "drop"
  ) %>%
  arrange(nace2d)
print(as.data.frame(by2d), row.names = FALSE)

# ── NACE 5-digit ────────────────────────────────────────────────────────────
stats5d <- ets %>%
  group_by(nace5d, year) %>%
  filter(n() >= MIN_N) %>%
  group_modify(~do_stats(.x)) %>%
  ungroup() %>%
  mutate(nace2d = substr(nace5d, 1, 2))

mfg5d <- stats5d %>% filter(nace2d >= "10", nace2d <= "33")

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  ETS ONLY — NACE 5-digit\n")
cat("══════════════════════════════════════════════════════════════\n\n")

cat("── Manufacturing mean (unweighted across sector-years) ──\n")
cat("  p90/p10 ratio:     ", round(mean(mfg5d$cp_p90p10, na.rm = TRUE), 1), "\n")
cat("  90-10 log diff:    ", round(mean(mfg5d$cp_log_9010, na.rm = TRUE), 2), "\n")
cat("  Var(log cp):       ", round(mean(mfg5d$cp_var_log, na.rm = TRUE), 2), "\n")

cat("\n── By NACE 2-digit parent (avg of 5d sectors within) ──\n")
by5d_parent <- mfg5d %>%
  group_by(nace2d) %>%
  summarise(
    n_5d_sectors = n_distinct(nace5d),
    n_obs        = n(),
    cp_p90p10    = round(mean(cp_p90p10, na.rm = TRUE), 1),
    cp_log_9010  = round(mean(cp_log_9010, na.rm = TRUE), 2),
    cp_var_log   = round(mean(cp_var_log, na.rm = TRUE), 2),
    .groups      = "drop"
  ) %>%
  arrange(nace2d)
print(as.data.frame(by5d_parent), row.names = FALSE)

cat("\n── Literature benchmarks ──\n")
cat("  Lyubich et al (2018): 90-10 log diff = 2.27, p90/p10 ~ 9.7x\n")
cat("  De Lyon & Dechezlepretre (2025): p90/p10 = 14-24x (3-4 digit)\n")
