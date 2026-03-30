###############################################################################
# Investigate IMJV and Climate TRACE for identifying deployment emitters
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
library(dplyr)

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: IMJV (Integraal Milieujaarverslag — Flemish environmental reports)
# ═══════════════════════════════════════════════════════════════════════════════
cat("═══════════════════════════════════════════════════════════════\n")
cat("  IMJV — Flemish Integrated Environmental Annual Reports\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

imjv <- read.delim(file.path(RAW_DATA, "IMJV", "imjv_ghg_lucht.tsv"),
                     stringsAsFactors = FALSE)

cat("Columns:", paste(names(imjv), collapse = ", "), "\n")
cat("Rows:", nrow(imjv), "\n")
cat("Years:", paste(sort(unique(imjv$year)), collapse = ", "), "\n")
cat("Pollutants:", paste(unique(imjv$polluent), collapse = ", "), "\n\n")

# Focus on CO2
imjv_co2 <- imjv %>%
  filter(polluent == "CO2") %>%
  mutate(emissions_t = emission_kg / 1000,  # convert to tonnes
         nace2d = substr(nace_code, 1, 2))

cat("── CO2 records ──\n")
cat("Rows:", nrow(imjv_co2), "\n")
cat("Unique firms:", n_distinct(imjv_co2$cbb_number), "\n\n")

# By NACE 2d sector
cat("── Firms by NACE 2d (CO2 only, unique firms) ──\n")
nace_summary <- imjv_co2 %>%
  group_by(nace2d) %>%
  summarise(
    n_firms = n_distinct(cbb_number),
    n_obs = n(),
    total_emissions_kt = round(sum(emissions_t, na.rm = TRUE) / 1000, 1),
    avg_per_firm_t = round(mean(emissions_t, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  arrange(desc(total_emissions_kt))
print(nace_summary, n = 40)

# NACE 38 specifically
cat("\n── NACE 38 in IMJV ──\n")
imjv_38 <- imjv_co2 %>% filter(nace2d == "38")
cat("CO2 records:", nrow(imjv_38), "\n")
cat("Unique firms:", n_distinct(imjv_38$cbb_number), "\n")

if (nrow(imjv_38) > 0) {
  cat("\nPer-firm summary:\n")
  firm_38 <- imjv_38 %>%
    group_by(cbb_number, firm_name) %>%
    summarise(
      n_years = n(),
      years = paste(range(year), collapse = "-"),
      avg_emissions_t = round(mean(emissions_t), 0),
      total_emissions_t = round(sum(emissions_t), 0),
      .groups = "drop"
    ) %>%
    arrange(desc(total_emissions_t))
  print(firm_38, n = 30)
}

# NACE 35 in IMJV
cat("\n── NACE 35 in IMJV ──\n")
imjv_35 <- imjv_co2 %>% filter(nace2d == "35")
cat("CO2 records:", nrow(imjv_35), "\n")
cat("Unique firms:", n_distinct(imjv_35$cbb_number), "\n")

if (nrow(imjv_35) > 0) {
  firm_35 <- imjv_35 %>%
    group_by(cbb_number, firm_name) %>%
    summarise(
      n_years = n(),
      avg_emissions_t = round(mean(emissions_t), 0),
      total_emissions_t = round(sum(emissions_t), 0),
      .groups = "drop"
    ) %>%
    arrange(desc(total_emissions_t))
  cat("\nTop 20 firms:\n")
  print(head(firm_35, 20), n = 20)
}

# How many IMJV firms are NOT in ETS?
# We can't check this directly without matching CBB to VAT,
# but we can look at the scale of emissions
cat("\n── Emission size distribution in IMJV (CO2, all sectors) ──\n")
firm_avg <- imjv_co2 %>%
  group_by(cbb_number, nace2d) %>%
  summarise(avg_t = mean(emissions_t, na.rm = TRUE), .groups = "drop")

cat(sprintf("  Total firms: %d\n", nrow(firm_avg)))
cat(sprintf("  Firms with avg > 25kt (likely ETS): %d\n", sum(firm_avg$avg_t > 25000)))
cat(sprintf("  Firms with avg 1-25kt (could be non-ETS): %d\n",
            sum(firm_avg$avg_t > 1000 & firm_avg$avg_t <= 25000)))
cat(sprintf("  Firms with avg < 1kt: %d\n", sum(firm_avg$avg_t <= 1000)))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: Climate TRACE — waste incineration
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  CLIMATE TRACE — Waste incineration, Belgium\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

ct_waste_file <- file.path(RAW_DATA, "Climate TRACE", "BEL", "DATA", "waste",
                            "incineration-and-open-burning-of-waste_country_emissions_v5_4_1.csv")
if (file.exists(ct_waste_file)) {
  ct_waste <- read.csv(ct_waste_file, stringsAsFactors = FALSE)
  cat("Waste incineration (country-level):\n")
  cat("Columns:", paste(names(ct_waste), collapse = ", "), "\n")
  # Show CO2 totals by year
  ct_waste_co2 <- ct_waste %>% filter(gas == "co2")
  if (nrow(ct_waste_co2) > 0) {
    cat("\nCO2 by year (tonnes):\n")
    yearly <- ct_waste_co2 %>%
      group_by(start_time) %>%
      summarise(total = sum(emissions_quantity, na.rm = TRUE), .groups = "drop")
    print(yearly, n = 30)
  }
}

# Check if there are facility-level waste files
ct_waste_sources <- file.path(RAW_DATA, "Climate TRACE", "BEL", "DATA", "waste",
                               "incineration-and-open-burning-of-waste_emissions_sources_v5_4_1.csv")
if (file.exists(ct_waste_sources)) {
  cat("\nFacility-level waste incineration data exists!\n")
  ct_ws <- read.csv(ct_waste_sources, stringsAsFactors = FALSE)
  cat("Rows:", nrow(ct_ws), "\n")
  cat("Columns:", paste(names(ct_ws), collapse = ", "), "\n")
} else {
  cat("\nNo facility-level waste incineration data.\n")
}

# Climate TRACE electricity — facility level
ct_elec <- file.path(RAW_DATA, "Climate TRACE", "BEL", "DATA", "power",
                      "electricity-generation_emissions_sources_v5_4_1.csv")
ct_e <- read.csv(ct_elec, stringsAsFactors = FALSE)
cat("\n── Climate TRACE electricity facilities ──\n")
cat("Rows:", nrow(ct_e), "\n")
# Unique facilities
n_fac <- n_distinct(ct_e$source_id)
cat("Unique facilities:", n_fac, "\n")

# Annual totals by source type
cat("\nFacilities by source_type:\n")
fac_types <- ct_e %>%
  filter(gas == "co2") %>%
  group_by(source_type) %>%
  summarise(
    n_fac = n_distinct(source_id),
    total_co2_kt = round(sum(emissions_quantity, na.rm = TRUE) / 1000, 0),
    .groups = "drop"
  ) %>%
  arrange(desc(total_co2_kt))
print(fac_types, n = 20)

# List facility names to check for waste incinerators
cat("\n── Facility names (unique, looking for waste/incineration) ──\n")
fac_names <- ct_e %>% distinct(source_id, source_name, source_type)
waste_fac <- fac_names %>%
  filter(grepl("waste|inciner|afval|verbrand", source_name, ignore.case = TRUE))
if (nrow(waste_fac) > 0) {
  cat("Waste-related facilities in electricity generation:\n")
  print(waste_fac, n = 20)
} else {
  cat("No waste-related facilities found in electricity generation data.\n")
}

cat("\n=== DONE ===\n")
