REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
load(file.path(PROC_DATA, "nir_calibration_targets.RData"))

nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE, colClasses = c(nace2d = "character")
) %>% select(nace2d, crf_group)

accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% 2008:2011) %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))

eutl <- firm_year_belgian_euets %>%
  filter(year %in% 2008:2011, !is.na(emissions), emissions > 0)

E_ETS_group <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group)) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions), n_ets = n(), .groups = "drop")

panel <- nir_targets %>%
  filter(year %in% 2008:2011) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(E_ETS_group, by = c("crf_group", "year")) %>%
  mutate(
    E_ETS      = coalesce(E_ETS, 0),
    E_deploy   = pmax(E_NIR - E_ETS, 0),
    pct_deploy = round(100 * E_deploy / E_NIR, 1)
  )

cat("All CRF groups (tonnes):\n\n")
print(
  panel %>%
    select(crf_group, year, E_NIR, E_ETS, E_deploy, n_ets, pct_deploy) %>%
    mutate(across(c(E_NIR, E_ETS, E_deploy), ~round(.x, 0))) %>%
    arrange(crf_group, year),
  n = 60
)
