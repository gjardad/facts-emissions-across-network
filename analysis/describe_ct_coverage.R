###############################################################################
# analysis/describe_ct_coverage.R
#
# PURPOSE
#   Describe Climate TRACE coverage in Belgium:
#     - Which sectors and subsectors are covered?
#     - What is the temporal coverage (year range by sector)?
#     - How many CT sources are NOT matched to any EUTL installation?
#
# INPUTS
#   {PROC_DATA}/ct_eutl_match.RData         (from match_ct_to_eutl.R)
#
# RUNS ON: local 1
###############################################################################

# ── Paths ────────────────────────────────────────────────────────────────────
REPO_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) normalizePath(getwd(), winslash = "/")
)
while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
source(file.path(REPO_DIR, "paths.R"))

library(dplyr)

# ── Load data ────────────────────────────────────────────────────────────────
load(file.path(PROC_DATA, "ct_eutl_match.RData"))  # ct, ct_eutl_match, be_inst


# =============================================================================
# 1. Overall summary
# =============================================================================
cat("== Climate TRACE Belgium — Overall Summary ==\n\n")
cat("Total source-years (CO2 > 0):", nrow(ct), "\n")
cat("Unique sources:", n_distinct(ct$source_id), "\n")
cat("Year range:", paste(range(ct$year), collapse = " – "), "\n")
cat("Sectors:", n_distinct(ct$sector), "\n")
cat("Subsectors:", n_distinct(ct$subsector), "\n\n")


# =============================================================================
# 2. Coverage by sector and subsector
# =============================================================================
cat("== Coverage by sector and subsector ==\n\n")

sector_summary <- ct %>%
  group_by(sector, subsector) %>%
  summarise(
    n_sources = n_distinct(source_id),
    n_source_years = n(),
    year_min = min(year),
    year_max = max(year),
    total_co2 = sum(ct_emissions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sector, subsector)

print(as.data.frame(sector_summary), right = FALSE)


# =============================================================================
# 3. Temporal coverage by sector
# =============================================================================
cat("\n== Temporal coverage by sector ==\n\n")

temporal <- ct %>%
  group_by(sector, year) %>%
  summarise(n_sources = n_distinct(source_id), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = year, values_from = n_sources, values_fill = 0)

print(as.data.frame(temporal), right = FALSE)


# =============================================================================
# 4. EUTL overlap
# =============================================================================
cat("\n== EUTL overlap ==\n\n")

all_sources <- ct %>%
  distinct(source_id, source_name, sector, subsector)

all_sources$matched <- all_sources$source_id %in% ct_eutl_match$source_id

overlap <- all_sources %>%
  group_by(sector, subsector) %>%
  summarise(
    n_total = n(),
    n_matched = sum(matched),
    n_unmatched = sum(!matched),
    pct_matched = round(100 * mean(matched), 1),
    .groups = "drop"
  ) %>%
  arrange(sector, subsector)

cat("By sector/subsector:\n")
print(as.data.frame(overlap), right = FALSE)

cat("\nOverall:\n")
cat("  Total CT sources:", nrow(all_sources), "\n")
cat("  Matched to EUTL:", sum(all_sources$matched), "\n")
cat("  Not matched:", sum(!all_sources$matched), "\n")
cat("  Match rate:", round(100 * mean(all_sources$matched), 1), "%\n")


# =============================================================================
# 5. Unmatched sources — who are they?
# =============================================================================
cat("\n== Unmatched CT sources (sample) ==\n\n")

unmatched <- all_sources %>%
  filter(!matched) %>%
  arrange(sector, subsector, source_name)

# Show up to 30
print(head(as.data.frame(unmatched %>% select(-matched)), 30), right = FALSE)

if (nrow(unmatched) > 30) {
  cat("... and", nrow(unmatched) - 30, "more unmatched sources\n")
}
