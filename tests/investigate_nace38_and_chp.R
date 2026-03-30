###############################################################################
# tests/investigate_nace38_and_chp.R
#
# 1. How many NACE 38 firms are in the deployment panel vs ETS?
# 2. What share of 1A1a emissions might come from waste incineration?
# 3. Can we identify CHP firms from the EUTL installation data?
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))
library(dplyr)

YEARS <- 2005:2021

# ── Load data ────────────────────────────────────────────────────────────────
load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
load(file.path(PROC_DATA, "deployment_panel.RData"))

accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = make_nace2d(nace5d))

eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: NACE 38 in deployment and ETS
# ═══════════════════════════════════════════════════════════════════════════════
cat("═══════════════════════════════════════════════════════════════\n")
cat("  NACE 38 (Waste collection, treatment, disposal, recovery)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ETS firms in NACE 38
ets_38 <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(nace2d == "38")

cat("── ETS firms in NACE 38 ──\n")
cat("  Firm-years:", nrow(ets_38), "\n")
cat("  Unique firms:", n_distinct(ets_38$vat), "\n")
if (nrow(ets_38) > 0) {
  cat("  Years present:", paste(sort(unique(ets_38$year)), collapse = ", "), "\n")
  ets_38_summary <- ets_38 %>%
    filter(!is.na(emissions), emissions > 0) %>%
    group_by(year) %>%
    summarise(n = n(), total = sum(emissions), mean = mean(emissions),
              min = min(emissions), max = max(emissions), .groups = "drop")
  cat("\n  Per-year summary:\n")
  print(ets_38_summary, n = 30)
}

# Deployment firms in NACE 38
deploy_38 <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  filter(nace2d == "38")

cat("\n── Deployment firms in NACE 38 ──\n")
cat("  Firm-years:", nrow(deploy_38), "\n")
cat("  Unique firms:", n_distinct(deploy_38$vat), "\n")
if (nrow(deploy_38) > 0) {
  deploy_38_by_year <- deploy_38 %>%
    group_by(year) %>%
    summarise(n = n(), .groups = "drop")
  cat("\n  Per-year counts:\n")
  print(deploy_38_by_year, n = 30)
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: Broader look at NACE sectors near energy/waste boundary
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  ETS FIRMS BY NACE 2D — WHICH SECTORS HAVE ETS FIRMS?\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

ets_by_nace <- eutl %>%
  filter(!is.na(emissions), emissions > 0) %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d)) %>%
  group_by(nace2d) %>%
  summarise(
    n_firm_years = n(),
    n_firms = n_distinct(vat),
    total_emissions = sum(emissions),
    avg_emissions = mean(emissions),
    .groups = "drop"
  ) %>%
  arrange(desc(total_emissions))

cat("All NACE 2d sectors with ETS firms:\n")
print(ets_by_nace, n = 50)

# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: ETS firms in NACE 35 specifically
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  NACE 35 (Electricity, gas, steam, air conditioning)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

ets_35 <- eutl %>%
  filter(!is.na(emissions), emissions > 0) %>%
  left_join(accounts %>% distinct(vat, year, nace2d, nace5d), by = c("vat", "year")) %>%
  filter(nace2d == "35")

cat("ETS firms in NACE 35:\n")
cat("  Firm-years:", nrow(ets_35), "\n")
cat("  Unique firms:", n_distinct(ets_35$vat), "\n")

ets_35_by_year <- ets_35 %>%
  group_by(year) %>%
  summarise(n = n(), total = sum(emissions), .groups = "drop")
cat("\n  Per-year:\n")
print(ets_35_by_year, n = 30)

# Deployment firms in NACE 35
deploy_35 <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  filter(nace2d == "35")

cat("\n  Deployment firms in NACE 35:\n")
cat("  Firm-years:", nrow(deploy_35), "\n")
cat("  Unique firms:", n_distinct(deploy_35$vat), "\n")

deploy_35_by_year <- deploy_35 %>%
  group_by(year) %>%
  summarise(n = n(), .groups = "drop")
cat("\n  Per-year deployment counts:\n")
print(deploy_35_by_year, n = 30)

# NACE 5-digit breakdown for NACE 35
cat("\n  NACE 5d breakdown for ETS firms in NACE 35:\n")
ets_35_5d <- ets_35 %>%
  group_by(nace5d) %>%
  summarise(n_firms = n_distinct(vat), total = sum(emissions), .groups = "drop") %>%
  arrange(desc(total))
print(ets_35_5d, n = 20)

cat("\n  NACE 5d breakdown for deployment firms in NACE 35:\n")
deploy_35_5d <- deploy_35 %>%
  group_by(nace5d) %>%
  summarise(n_firms = n_distinct(vat), .groups = "drop") %>%
  arrange(desc(n_firms))
print(deploy_35_5d, n = 20)

# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: Can we identify CHP from EUTL installation data?
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  CHP IDENTIFICATION FROM EUTL DATA\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Check if installation data has activity type or description
install_path <- file.path(PROC_DATA, "installation.RData")
if (file.exists(install_path)) {
  load(install_path)
  cat("Installation data loaded.\n")
  cat("Columns:", paste(names(installation), collapse = ", "), "\n")
} else {
  # Try raw EUTL data
  install_files <- list.files(file.path(RAW_DATA, "EUTL", "Oct_2024"),
                               pattern = "install", full.names = TRUE)
  cat("Installation .RData not found in processed.\n")
  cat("Raw EUTL files:\n")
  cat(paste(" ", install_files), sep = "\n")

  if (length(install_files) > 0) {
    # Try loading the first one
    f <- install_files[1]
    cat("\nLoading:", f, "\n")
    if (grepl("\\.csv$", f)) {
      inst <- read.csv(f, stringsAsFactors = FALSE)
    } else if (grepl("\\.RData$", f)) {
      load(f)
      inst <- get(ls()[length(ls())])
    }
    if (exists("inst")) {
      cat("Columns:", paste(names(inst), collapse = ", "), "\n")
      cat("Rows:", nrow(inst), "\n")
      # Look for activity type columns
      act_cols <- grep("activ|type|descr|nace|sector", names(inst),
                       ignore.case = TRUE, value = TRUE)
      cat("Activity-related columns:", paste(act_cols, collapse = ", "), "\n")
      if (length(act_cols) > 0) {
        for (col in act_cols) {
          cat(sprintf("\n  Unique values in '%s' (first 20):\n", col))
          vals <- sort(unique(inst[[col]]))
          cat(paste("   ", head(vals, 20)), sep = "\n")
        }
      }
    }
  }
}

cat("\n=== DONE ===\n")
