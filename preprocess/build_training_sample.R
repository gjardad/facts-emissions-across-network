###############################################################################
# preprocess/build_training_sample.R
#
# PURPOSE
#   Construct the canonical training sample for the elastic net proxy and all
#   downstream models.
#
#   The sample contains:
#     - EU ETS firms: observed verified emissions.
#     - Non-ETS firms in sectors with near-complete EU ETS coverage, assigned
#       zero combustion emissions:
#         C19  Oil refining          (100%, CRF 1.A.1.b)
#         C17  Paper                 (~97%, CRF 1.A.2.d)
#         C18  Printing              (~97%, CRF 1.A.2.d)
#         C24  Iron & steel only     (~95%, CRF 1.A.2.a)
#              NACE 24.1, 24.2, 24.3, 24.51, 24.52
#              Excludes non-ferrous metals (24.4x, 24.53, 24.54) where
#              EU ETS covers only ~70% (CRF 1.A.2.b).
#
# INPUT
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   {PROC_DATA}/training_sample.RData
#     Contains: training_sample (data.frame)
#     Columns: vat, year, capital, revenue, wage_bill, fte, value_added,
#              nace5d, emissions, euets, nace2d
#
# RUNS ON: RMD
# COPIED FROM: inferring_emissions/preprocess/build_training_sample.R
###############################################################################

# ── Paths ────────────────────────────────────────────────────────────────────
if (tolower(Sys.info()[["user"]]) == "jardang") {
  REPO_DIR <- "C:/Users/jardang/Documents/facts-emissions-across-network"
} else if (tolower(Sys.info()[["user"]]) == "jota_"){
  REPO_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
                        error = function(e) normalizePath(getwd(), winslash = "/"))
  while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
} else {
  stop("Define REPO_DIR for this user.")
}
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)

# ── Load data ────────────────────────────────────────────────────────────────
load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))

euets <- firm_year_belgian_euets %>%
  select(vat, year, emissions) %>%
  mutate(euets = 1)

# ── Define non-ETS zero-emission sectors ─────────────────────────────────────
ZERO_NACE2D <- c("17/18", "19")
ZERO_C24_NACE4D <- c("2410", "2420", "2431", "2432", "2433", "2434",
                      "2451", "2452")

# ── Build sample ─────────────────────────────────────────────────────────────
training_sample <- df_annual_accounts_selected_sample_key_variables %>%
  left_join(euets, by = c("vat", "year")) %>%
  mutate(
    nace2d = make_nace2d(nace5d),
    nace4d = substr(nace5d, 1, 4),
    euets  = coalesce(euets, 0),
    is_zero_sector = nace2d %in% ZERO_NACE2D |
                     (nace2d == "24" & nace4d %in% ZERO_C24_NACE4D)
  ) %>%
  filter(euets == 1 | is_zero_sector,
         year >= 2005) %>%
  # Non-ETS firms in zero-emission sectors: emissions = 0
  mutate(emissions = if_else(is.na(emissions) & is_zero_sector, 0, emissions)) %>%
  select(-nace4d, -is_zero_sector)

# Drop EU ETS firm-years where emissions is NA (pre-regulation years)
n_na_euets <- sum(is.na(training_sample$emissions))
cat("EU ETS firm-years with NA emissions (not yet regulated):", n_na_euets, "\n")

training_sample <- training_sample %>%
  filter(!is.na(emissions))

cat("Training sample:", nrow(training_sample), "firm-years,",
    n_distinct(training_sample$vat), "firms",
    "(dropped", n_na_euets, "pre-regulation obs)\n")

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n── Sample composition ──\n")
cat("EU ETS firms:", n_distinct(training_sample$vat[training_sample$euets == 1]), "\n")
cat("Non-ETS zero firms:", n_distinct(training_sample$vat[training_sample$euets == 0]), "\n")
cat("Sectors:", paste(sort(unique(training_sample$nace2d)), collapse = ", "), "\n")

# ── Save ─────────────────────────────────────────────────────────────────────
save(training_sample, file = file.path(PROC_DATA, "training_sample.RData"))
cat("\nSaved:", file.path(PROC_DATA, "training_sample.RData"), "\n")
