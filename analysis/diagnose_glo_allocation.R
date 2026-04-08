###############################################################################
# analysis/diagnose_glo_allocation.R
#
# PURPOSE
#   Diagnostic on the GLO allocation output:
#
#     1. Imputed firm-years with scope1 > 35 kt (35,000 tons; physics-based
#        upper bound).
#     2. Cells where max(scope1[source==imputed]) > min(scope1[source==ets]).
#        (Pre-ETS firms are EXCLUDED from the "min ETS" computation, by
#        instruction -- only true ETS firms count for that comparison.)
#     3. Imputed firm-years exceeding the cell's min true-ETS scope1.
#
#   Reported overall, by CRF group, and by year.
#
# INPUT
#   {PROC_DATA}/allocation_glo/alloc_YYYY.RData
#
# OUTPUT
#   {REPO_DIR}/output/tables/glo_allocation_diagnostic.csv
###############################################################################

# -- Paths --------------------------------------------------------------------
if (tolower(Sys.info()[["user"]]) == "jardang") {
  REPO_DIR <- "C:/Users/jardang/Documents/facts-emissions-across-network"
} else if (tolower(Sys.info()[["user"]]) == "jota_") {
  REPO_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
                        error = function(e) normalizePath(getwd(), winslash = "/"))
  while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
} else {
  stop("Define REPO_DIR for this user.")
}
source(file.path(REPO_DIR, "paths.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

CAP_TONS <- 35000  # 35 kt physics-based bound

cat("===================================================================\n")
cat("  GLO ALLOCATION DIAGNOSTIC\n")
cat("===================================================================\n\n")

# -- Load all years -----------------------------------------------------------
in_dir <- file.path(PROC_DATA, "allocation_glo")
files  <- list.files(in_dir, pattern = "^alloc_\\d{4}\\.RData$", full.names = TRUE)
if (length(files) == 0) stop("No allocation files found in ", in_dir)

all_firms <- list()
for (f in files) {
  load(f)  # year_firms
  all_firms[[length(all_firms) + 1]] <- year_firms
}
firms <- bind_rows(all_firms)
rm(all_firms)

cat("Loaded", length(files), "years,", nrow(firms), "firm-year rows\n\n")

# -- Min ETS per cell (true ETS only; pre_ets excluded) -----------------------
min_ets_cell <- firms %>%
  filter(source == "ets") %>%
  group_by(crf_group, year) %>%
  summarise(min_ets = min(scope1, na.rm = TRUE),
            n_ets   = n(),
            .groups = "drop")

# -- Imputed firms with cell context ------------------------------------------
imputed <- firms %>%
  filter(source == "imputed") %>%
  left_join(min_ets_cell, by = c("crf_group", "year")) %>%
  mutate(
    over_cap        = scope1 > CAP_TONS,
    over_min_ets    = !is.na(min_ets) & scope1 > min_ets
  )

n_imputed_total <- nrow(imputed)
cat("Imputed firm-years:", n_imputed_total, "\n")

# -- Cells with both ETS and imputed ------------------------------------------
cell_summary <- imputed %>%
  filter(!is.na(min_ets)) %>%
  group_by(crf_group, year) %>%
  summarise(
    n_imputed_in_cell = n(),
    max_imputed       = max(scope1, na.rm = TRUE),
    min_ets           = first(min_ets),
    .groups = "drop"
  ) %>%
  mutate(cell_violates = max_imputed > min_ets)

# =============================================================================
# Diagnostic 1: count of imputed > 35 kt
# =============================================================================
cat("\n--- Imputed firm-years with scope1 > 35 kt ---\n")
n_over_cap <- sum(imputed$over_cap)
cat(sprintf("  Total: %d / %d (%.2f%%)\n",
            n_over_cap, n_imputed_total, 100 * n_over_cap / max(n_imputed_total, 1)))

over_cap_by_crf <- imputed %>%
  group_by(crf_group) %>%
  summarise(n_imputed = n(),
            n_over_cap = sum(over_cap),
            pct_over   = 100 * mean(over_cap),
            .groups = "drop") %>%
  arrange(desc(n_over_cap))

over_cap_by_year <- imputed %>%
  group_by(year) %>%
  summarise(n_imputed = n(),
            n_over_cap = sum(over_cap),
            pct_over   = 100 * mean(over_cap),
            .groups = "drop") %>%
  arrange(year)

cat("\n  By CRF group (top 10):\n")
print(head(over_cap_by_crf, 10))
cat("\n  By year:\n")
print(over_cap_by_year)

# =============================================================================
# Diagnostic 2: cells where max imputed > min true-ETS
# =============================================================================
cat("\n--- Cells where max(imputed) > min(true ETS) ---\n")
n_cells_with_both <- nrow(cell_summary)
n_cells_violate   <- sum(cell_summary$cell_violates)
cat(sprintf("  Cells with both ETS and imputed: %d\n", n_cells_with_both))
cat(sprintf("  Of those, cells where max imputed > min ETS: %d (%.2f%%)\n",
            n_cells_violate, 100 * n_cells_violate / max(n_cells_with_both, 1)))

cell_violate_by_crf <- cell_summary %>%
  group_by(crf_group) %>%
  summarise(n_cells          = n(),
            n_cells_violate  = sum(cell_violates),
            pct_cells_violate = 100 * mean(cell_violates),
            .groups = "drop") %>%
  arrange(desc(n_cells_violate))

cell_violate_by_year <- cell_summary %>%
  group_by(year) %>%
  summarise(n_cells          = n(),
            n_cells_violate  = sum(cell_violates),
            pct_cells_violate = 100 * mean(cell_violates),
            .groups = "drop") %>%
  arrange(year)

cat("\n  By CRF group:\n"); print(cell_violate_by_crf)
cat("\n  By year:\n");      print(cell_violate_by_year)

# =============================================================================
# Diagnostic 3: imputed firm-years exceeding min true-ETS in their cell
# =============================================================================
cat("\n--- Imputed firm-years with scope1 > min(true ETS) in cell ---\n")
imputed_with_ets <- imputed %>% filter(!is.na(min_ets))
n_imp_in_mixed_cells <- nrow(imputed_with_ets)
n_imp_violate        <- sum(imputed_with_ets$over_min_ets)
cat(sprintf("  Imputed firms in cells that also have ETS: %d\n",
            n_imp_in_mixed_cells))
cat(sprintf("  Of those, imputed > min ETS: %d (%.2f%%)\n",
            n_imp_violate, 100 * n_imp_violate / max(n_imp_in_mixed_cells, 1)))

firm_violate_by_crf <- imputed_with_ets %>%
  group_by(crf_group) %>%
  summarise(n_imputed       = n(),
            n_over_min_ets  = sum(over_min_ets),
            pct_over_min_ets = 100 * mean(over_min_ets),
            .groups = "drop") %>%
  arrange(desc(n_over_min_ets))

firm_violate_by_year <- imputed_with_ets %>%
  group_by(year) %>%
  summarise(n_imputed       = n(),
            n_over_min_ets  = sum(over_min_ets),
            pct_over_min_ets = 100 * mean(over_min_ets),
            .groups = "drop") %>%
  arrange(year)

cat("\n  By CRF group:\n"); print(firm_violate_by_crf)
cat("\n  By year:\n");      print(firm_violate_by_year)

# =============================================================================
# Save consolidated table
# =============================================================================
out_dir  <- file.path(REPO_DIR, "output", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_path <- file.path(out_dir, "glo_allocation_diagnostic.csv")

overall_row <- data.frame(
  scope                 = "overall",
  group                 = "ALL",
  n_imputed             = n_imputed_total,
  n_over_35kt           = n_over_cap,
  pct_over_35kt         = round(100 * n_over_cap / max(n_imputed_total, 1), 3),
  n_cells_with_both     = n_cells_with_both,
  n_cells_violate       = n_cells_violate,
  pct_cells_violate     = round(100 * n_cells_violate / max(n_cells_with_both, 1), 3),
  n_imp_in_mixed_cells  = n_imp_in_mixed_cells,
  n_imp_over_min_ets    = n_imp_violate,
  pct_imp_over_min_ets  = round(100 * n_imp_violate / max(n_imp_in_mixed_cells, 1), 3),
  stringsAsFactors      = FALSE
)

by_crf_row <- over_cap_by_crf %>%
  full_join(cell_violate_by_crf, by = "crf_group") %>%
  full_join(firm_violate_by_crf %>% rename(n_imputed_mix = n_imputed),
            by = "crf_group") %>%
  transmute(
    scope                 = "by_crf",
    group                 = crf_group,
    n_imputed             = n_imputed,
    n_over_35kt           = n_over_cap,
    pct_over_35kt         = round(pct_over, 3),
    n_cells_with_both     = n_cells,
    n_cells_violate       = n_cells_violate,
    pct_cells_violate     = round(pct_cells_violate, 3),
    n_imp_in_mixed_cells  = n_imputed_mix,
    n_imp_over_min_ets    = n_over_min_ets,
    pct_imp_over_min_ets  = round(pct_over_min_ets, 3)
  )

by_year_row <- over_cap_by_year %>%
  full_join(cell_violate_by_year, by = "year") %>%
  full_join(firm_violate_by_year %>% rename(n_imputed_mix = n_imputed),
            by = "year") %>%
  transmute(
    scope                 = "by_year",
    group                 = as.character(year),
    n_imputed             = n_imputed,
    n_over_35kt           = n_over_cap,
    pct_over_35kt         = round(pct_over, 3),
    n_cells_with_both     = n_cells,
    n_cells_violate       = n_cells_violate,
    pct_cells_violate     = round(pct_cells_violate, 3),
    n_imp_in_mixed_cells  = n_imputed_mix,
    n_imp_over_min_ets    = n_over_min_ets,
    pct_imp_over_min_ets  = round(pct_over_min_ets, 3)
  )

out_table <- bind_rows(overall_row, by_crf_row, by_year_row)
write.csv(out_table, out_path, row.names = FALSE)
cat("\nSaved:", out_path, "\n")
