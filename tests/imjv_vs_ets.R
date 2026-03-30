###############################################################################
# Check if IMJV NACE 38 firms are already in ETS
# Also check overall IMJV vs ETS overlap
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
library(dplyr)

imjv <- read.delim(file.path(RAW_DATA, "IMJV", "imjv_ghg_lucht.tsv"),
                     stringsAsFactors = FALSE)
imjv_co2 <- imjv %>%
  filter(polluent == "CO2") %>%
  mutate(emissions_t = emission_kg / 1000,
         nace2d = substr(nace_code, 1, 2))

# EUTL installation data — has names we can fuzzy-match
inst <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "installation.csv"),
                  stringsAsFactors = FALSE)
be_inst <- inst %>%
  filter(country_id == "BE",
         isAircraftOperator == "False" | isAircraftOperator == FALSE,
         isMaritimeOperator == "False" | isMaritimeOperator == FALSE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  IMJV NACE 38 firms — are they in EU ETS?\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# The NACE 38 firms in IMJV
imjv_38_firms <- imjv_co2 %>%
  filter(nace2d == "38") %>%
  group_by(cbb_number, firm_name) %>%
  summarise(avg_t = round(mean(emissions_t)), .groups = "drop") %>%
  distinct(cbb_number, firm_name, avg_t)

cat("IMJV NACE 38 firms and their names:\n")
for (i in seq_len(nrow(imjv_38_firms))) {
  cat(sprintf("  CBB=%s | avg=%dt | %s\n",
              imjv_38_firms$cbb_number[i],
              imjv_38_firms$avg_t[i],
              imjv_38_firms$firm_name[i]))
}

# Now search for these names in EUTL
cat("\n── Searching EUTL installations for matching names ──\n\n")
keywords_38 <- c("IVAGO", "ISVAG", "INTERGEMEENTELIJK", "RECYPER", "SITA",
                  "VANHEEDE", "DALKIA", "VEOLIA", "Gentse Warmte", "SARPI")

for (kw in keywords_38) {
  matches <- be_inst %>%
    filter(grepl(kw, name, ignore.case = TRUE))
  if (nrow(matches) > 0) {
    for (j in seq_len(nrow(matches))) {
      cat(sprintf("  MATCH '%s': [%s] act=%s nace=%s | %s\n",
                  kw, matches$id[j], matches$activity_id[j],
                  matches$nace_id[j], matches$name[j]))
    }
  } else {
    cat(sprintf("  NO MATCH for '%s'\n", kw))
  }
}

# Also check: the ETS firms we identified as NACE 38 in annual accounts
# (from the previous analysis: 1-2 firms, ~10-25k tonnes)
cat("\n\n── Cross-check: ETS firms classified as NACE 38 in annual accounts ──\n")
load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

accounts <- df_annual_accounts_selected_sample_key_variables %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))

ets_38 <- firm_year_belgian_euets %>%
  filter(!is.na(emissions), emissions > 0) %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(nace2d == "38")

cat("ETS firm-years with NACE 38:", nrow(ets_38), "\n")
cat("Unique ETS firms with NACE 38:", n_distinct(ets_38$vat), "\n")
if (nrow(ets_38) > 0) {
  ets_38_summary <- ets_38 %>%
    group_by(vat) %>%
    summarise(
      n_years = n(),
      avg_emissions = round(mean(emissions)),
      years = paste(range(year), collapse = "-"),
      .groups = "drop"
    )
  cat("\nETS NACE 38 firms:\n")
  for (i in seq_len(nrow(ets_38_summary))) {
    cat(sprintf("  vat=%s | n_years=%d | avg=%dt | %s\n",
                ets_38_summary$vat[i], ets_38_summary$n_years[i],
                ets_38_summary$avg_emissions[i], ets_38_summary$years[i]))
  }
}

cat("\n=== DONE ===\n")
