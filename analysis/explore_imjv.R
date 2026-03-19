###############################################################################
# analysis/explore_imjv.R
#
# PURPOSE
#   Explore the IMJV (Integraal Milieujaarverslag) data — Flanders' mandatory
#   environmental reporting system for firms consuming >0.1 PJ primary energy.
#
# INPUTS
#   {RAW_DATA}/IMJV/imjv_ghg_lucht.tsv
#
# RUNS ON: local 1
###############################################################################

# -- Paths -------------------------------------------------------------------
REPO_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) normalizePath(getwd(), winslash = "/")
)
while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
source(file.path(REPO_DIR, "paths.R"))

library(dplyr)

# -- Load data ---------------------------------------------------------------
imjv <- read.delim(file.path(RAW_DATA, "IMJV", "imjv_ghg_lucht.tsv"),
                   stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

cat("== IMJV Data Overview ==\n\n")
cat("Rows:", nrow(imjv), "\n")
cat("Columns:", paste(names(imjv), collapse = ", "), "\n")
cat("Unique firms (cbb_number):", n_distinct(imjv$cbb_number), "\n")
cat("Year range:", paste(range(imjv$year), collapse = " – "), "\n\n")


# =============================================================================
# 1. Pollutant breakdown
# =============================================================================
cat("== Pollutant breakdown ==\n\n")
poll_summary <- imjv %>%
  group_by(polluent) %>%
  summarise(
    n_obs = n(),
    n_firms = n_distinct(cbb_number),
    total_kg = sum(emission_kg, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(poll_summary))


# =============================================================================
# 2. CO2 emissions distribution
# =============================================================================
cat("\n== CO2 emissions distribution (kg) ==\n\n")
co2 <- imjv %>% filter(polluent == "CO2")
cat("N firm-years:", nrow(co2), "\n")
cat("N unique firms:", n_distinct(co2$cbb_number), "\n\n")

# Convert to tonnes for readability
co2 <- co2 %>% mutate(emission_t = emission_kg / 1000)

cat("In tonnes:\n")
print(summary(co2$emission_t))

cat("\nPercentiles:\n")
cat(sprintf("  p10: %.0f t\n", quantile(co2$emission_t, 0.10)))
cat(sprintf("  p25: %.0f t\n", quantile(co2$emission_t, 0.25)))
cat(sprintf("  p50: %.0f t\n", quantile(co2$emission_t, 0.50)))
cat(sprintf("  p75: %.0f t\n", quantile(co2$emission_t, 0.75)))
cat(sprintf("  p90: %.0f t\n", quantile(co2$emission_t, 0.90)))
cat(sprintf("  p99: %.0f t\n", quantile(co2$emission_t, 0.99)))


# =============================================================================
# 3. Top emitters by mean annual CO2
# =============================================================================
cat("\n== Top 20 emitters by mean annual CO2 (tonnes) ==\n\n")
top_emitters <- co2 %>%
  group_by(cbb_number, firm_name, nace_code) %>%
  summarise(
    mean_annual_t = mean(emission_t),
    years = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_annual_t)) %>%
  head(20)

print(as.data.frame(top_emitters %>%
  select(firm_name, nace_code, mean_annual_t, years)), row.names = FALSE)


# =============================================================================
# 4. Temporal coverage per firm
# =============================================================================
cat("\n== Temporal coverage (years of CO2 data per firm) ==\n\n")
firm_years <- co2 %>%
  group_by(cbb_number) %>%
  summarise(n_years = n_distinct(year), .groups = "drop")

print(table(firm_years$n_years))
cat(sprintf("\nMean years per firm: %.1f\n", mean(firm_years$n_years)))


# =============================================================================
# 5. NACE 2-digit sector breakdown
# =============================================================================
cat("\n== NACE 2-digit sector breakdown (CO2) ==\n\n")
co2 <- co2 %>%
  mutate(nace2d = as.integer(sub("\\..*", "", nace_code)))

nace_summary <- co2 %>%
  group_by(nace2d) %>%
  summarise(
    n_firms = n_distinct(cbb_number),
    n_firmyears = n(),
    total_t = sum(emission_t, na.rm = TRUE),
    mean_firm_t = mean(emission_t, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_t))

print(as.data.frame(nace_summary))

cat(sprintf("\nTotal CO2 across all IMJV firms: %.0f tonnes\n", sum(co2$emission_t)))
cat(sprintf("Mean annual total: %.0f tonnes\n",
            sum(co2$emission_t) / n_distinct(co2$year)))


# =============================================================================
# 6. Year-level totals
# =============================================================================
cat("\n== Annual CO2 totals (tonnes) ==\n\n")
annual <- co2 %>%
  group_by(year) %>%
  summarise(
    n_firms = n_distinct(cbb_number),
    total_t = sum(emission_t, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(annual))
