library(dplyr, warn.conflicts = FALSE)

DATA_DIR  <- "C:/Users/jota_/Documents/NBB_data/processed"
REPO_DIR  <- "C:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

load(file.path(DATA_DIR, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year >= 2005, year <= 2021, !is.na(emissions), emissions > 0) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

load(file.path(DATA_DIR, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year >= 2005, year <= 2021) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

ets25 <- eutl %>%
  inner_join(accounts, by = c("vat", "year")) %>%
  filter(nace2d == "25", revenue > 0) %>%
  mutate(cp = revenue / emissions)

cat("══════════════════════════════════════════════════════════════\n")
cat("  SECTOR 25 — ETS ONLY — YEAR BY YEAR\n")
cat("══════════════════════════════════════════════════════════════\n\n")

pct_ratio <- function(x, p_hi, p_lo) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  q <- quantile(x, c(p_lo, p_hi), names = FALSE)
  if (q[1] == 0) return(NA_real_)
  q[2] / q[1]
}

by_year <- ets25 %>%
  group_by(year) %>%
  summarise(
    n_firms     = n(),
    min_e       = round(min(emissions)),
    max_e       = round(max(emissions)),
    cp_p90p10   = round(pct_ratio(cp, 0.9, 0.1), 1),
    cp_log_9010 = round(diff(quantile(log(cp), c(0.1, 0.9))), 2),
    .groups     = "drop"
  )

print(as.data.frame(by_year), row.names = FALSE)

cat("\n── Which NACE 5-digit codes are these firms in? ──\n")
ets25 %>%
  count(nace5d) %>%
  arrange(desc(n)) %>%
  print(row.names = FALSE)

cat("\n── Unique firms per year ──\n")
cat("Total unique ETS firms in NACE 25:", n_distinct(ets25$vat), "\n")
