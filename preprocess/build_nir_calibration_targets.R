###############################################################################
# preprocess/build_nir_calibration_targets.R
#
# PURPOSE
#   Read BEL-CRT Excel files for YEARS 2005–2022, extract CO2 emissions for
#   the CRF categories in scope, map to CRF groups via the NACE–CRF crosswalk,
#   and save nir_calibration_targets.RData.
#
#   Only CO2 is extracted (kt). See preprocess/crosswalks/NACE_CRF_CROSSWALK.md for
#   the rationale and a note on future extension to total GHG (CO2e).
#
# EXTRACTION LOGIC
#   Table1  (fuel combustion):      rows whose CRF code has exactly 4 dots,
#                                   e.g. "1.A.2.g." → compact "1A2g".
#                                   These are the finest level in Table1,
#                                   so no parent rows are included → no
#                                   double-counting.
#   Table2(I) (industrial processes): rows with exactly 2 dots,
#                                   e.g. "2.A." → compact "2A".
#                                   These are the sector aggregates in Table2(I)
#                                   (sub-rows like "2.A.1." are skipped to
#                                   avoid double-counting).
#   "NO" / "IE" / "NA,NO" in the CO2 column are treated as 0 and logged.
#
# INPUTS
#   {RAW_DATA}/NIR/BEL-CRT-2025-V1.0-{YYYY}-*_awaiting submission.xlsx
#   {REPO_DIR}/preprocess/crosswalks/crf_group_definitions.csv
#
# OUTPUT
#   {PROC_DATA}/nir_calibration_targets.RData
#     nir_targets : data frame  (crf_group, year, E_NIR_kt)
#
# RUNS ON: local 1
###############################################################################

# ── Paths ─────────────────────────────────────────────────────────────────────
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

library(dplyr)
library(readxl)

YEARS <- 2005:2022
NIR_DIR   <- file.path(RAW_DATA, "NIR")
CROSS_DIR <- file.path(REPO_DIR, "preprocess", "crosswalks")

cat("══════════════════════════════════════════════════════════\n")
cat("  BUILD NIR CALIBRATION TARGETS\n")
cat("  Years:", min(YEARS), "–", max(YEARS), "\n")
cat("══════════════════════════════════════════════════════════\n\n")


# =============================================================================
# SECTION 1: Load CRF group definitions
# =============================================================================
crf_defs <- read.csv(
  file.path(CROSS_DIR, "crf_group_definitions.csv"),
  stringsAsFactors = FALSE
)
cat("CRF groups loaded:", nrow(crf_defs), "\n")
cat(paste0("  ", crf_defs$crf_group, "  [", crf_defs$crf_pattern, "]",
           collapse = "\n"), "\n\n")


# =============================================================================
# SECTION 2: Helper — read one sheet and return tidy (code, co2_kt) rows
# =============================================================================

#' Parse a BEL-CRT sheet into a tidy data frame.
#'
#' @param path   Path to the xlsx file.
#' @param sheet  Sheet name ("Table1" or "Table2(I)").
#' @param n_dots Target number of dots in the CRF code (4 for Table1, 2 for
#'               Table2(I)).
#' @return data frame with columns: code (compact, e.g. "1A2g"), co2_kt.
read_crt_sheet <- function(path, sheet, n_dots) {

  raw <- read_excel(path, sheet = sheet, col_names = FALSE,
                    .name_repair = "minimal")

  # Identify CO2 column: find the column containing "CO2" in the header rows
  # (rows 1–10 are header/title; "CO2" appears in the column-header row)
  co2_col <- NA_integer_
  for (r in seq_len(min(nrow(raw), 12L))) {
    cols_with_co2 <- which(vapply(raw[r, ], function(x) {
      isTRUE(trimws(as.character(x)) == "CO2")
    }, logical(1L)))
    if (length(cols_with_co2) > 0L) {
      co2_col <- cols_with_co2[1L]
      break
    }
  }
  if (is.na(co2_col))
    stop(sprintf("Could not find 'CO2' column header in sheet '%s' of %s",
                 sheet, basename(path)))

  label_col <- 1L  # CRF label always in column A

  # Extract rows: label col must contain a CRF-style code prefix (starts with a digit)
  labels <- as.character(raw[[label_col]])
  co2_raw <- as.character(raw[[co2_col]])

  # Parse: extract code prefix (everything before first space), count dots
  code_prefix <- sub(" .*", "", trimws(labels))          # e.g. "1.A.2.g."
  n_dots_vec  <- nchar(code_prefix) -
                 nchar(gsub(".", "", code_prefix, fixed = TRUE))
  code_compact <- gsub(".", "", code_prefix, fixed = TRUE)    # e.g. "1A2g"

  # Filter to rows at the target dot level and non-empty labels
  keep <- !is.na(labels) &
          nchar(trimws(labels)) > 0 &
          n_dots_vec == n_dots &
          grepl("^[0-9]", code_prefix)      # must start with a digit

  if (!any(keep)) {
    warning(sprintf("No rows with %d dots in sheet '%s' of %s",
                    n_dots, sheet, basename(path)))
    return(data.frame(code = character(), co2_kt = numeric(),
                      stringsAsFactors = FALSE))
  }

  # Parse CO2: non-numeric entries ("NO", "IE", "NA,NO", etc.) → NA → 0
  co2_num <- suppressWarnings(as.numeric(co2_raw[keep]))

  n_nonnum <- sum(is.na(co2_num))
  if (n_nonnum > 0L)
    message(sprintf("  [%s / %s] %d CO2 value(s) non-numeric → treated as 0",
                    basename(path), sheet, n_nonnum))

  co2_num[is.na(co2_num)] <- 0

  data.frame(
    code   = code_compact[keep],
    co2_kt = co2_num,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# SECTION 3: Map compact CRF codes to groups
# =============================================================================

assign_groups <- function(df_codes) {
  # df_codes: (code, co2_kt)
  if (nrow(df_codes) == 0L)
    return(df_codes %>% mutate(crf_group = character(0)))
  df_codes$crf_group <- rep(NA_character_, nrow(df_codes))

  for (i in seq_len(nrow(crf_defs))) {
    pat   <- crf_defs$crf_pattern[i]
    group <- crf_defs$crf_group[i]
    match <- grepl(pat, df_codes$code)
    # Warn if a code matches multiple groups (should not happen with well-defined patterns)
    already_assigned <- !is.na(df_codes$crf_group) & match
    if (any(already_assigned))
      warning(sprintf(
        "Code(s) %s match group '%s' but were already assigned to '%s'",
        paste(df_codes$code[already_assigned], collapse = ", "),
        group,
        paste(unique(df_codes$crf_group[already_assigned]), collapse = ", ")
      ))
    df_codes$crf_group[match & is.na(df_codes$crf_group)] <- group
  }

  # Codes not matched by any group are dropped (e.g. 1A3 transport)
  n_unmatched <- sum(is.na(df_codes$crf_group))
  if (n_unmatched > 0L)
    message(sprintf("  %d code(s) not matched to any CRF group (excluded): %s",
                    n_unmatched,
                    paste(df_codes$code[is.na(df_codes$crf_group)], collapse = ", ")))

  df_codes[!is.na(df_codes$crf_group), ]
}


# =============================================================================
# SECTION 4: Main loop — one file per year
# =============================================================================
cat("── Reading BEL-CRT files ────────────────────────────────\n")

all_years <- vector("list", length(YEARS))

for (i in seq_along(YEARS)) {
  yr <- YEARS[i]

  # Find file (glob: year appears between two hyphens in the filename)
  pattern <- paste0("BEL-CRT-.*-", yr, "-.*\\.xlsx$")
  files   <- list.files(NIR_DIR, pattern = pattern, full.names = TRUE)

  if (length(files) == 0L) {
    warning(sprintf("No BEL-CRT file found for year %d — skipping", yr))
    next
  }
  if (length(files) > 1L) {
    warning(sprintf("Multiple BEL-CRT files for year %d — using first: %s",
                    yr, basename(files[1L])))
  }
  path <- files[1L]
  cat(sprintf("  %d  %s\n", yr, basename(path)))

  # Read Table1 (4-dot level) and Table2(I) (2-dot level)
  t1 <- tryCatch(
    read_crt_sheet(path, "Table1",    n_dots = 4L),
    error = function(e) { warning(sprintf("Table1 error year %d: %s", yr, e$message)); NULL }
  )
  t2 <- tryCatch(
    read_crt_sheet(path, "Table2(I)", n_dots = 2L),
    error = function(e) { warning(sprintf("Table2(I) error year %d: %s", yr, e$message)); NULL }
  )

  if (is.null(t1) && is.null(t2)) next

  combined <- bind_rows(t1, t2)

  # Assign CRF groups
  combined <- assign_groups(combined)

  # Sum to group level (a group may have multiple matching rows, e.g. metals
  # = 1A1c + 1A2a + 1A2b + 2C)
  year_totals <- combined %>%
    group_by(crf_group) %>%
    summarise(E_NIR_kt = sum(co2_kt, na.rm = TRUE), .groups = "drop") %>%
    mutate(year = yr)

  all_years[[i]] <- year_totals
}

cat("\n")


# =============================================================================
# SECTION 5: Combine, verify, save
# =============================================================================
cat("── Combining years ──────────────────────────────────────\n")

nir_targets <- bind_rows(all_years)

# Every group should appear in every year — flag missing combinations
expected <- expand.grid(
  crf_group = crf_defs$crf_group,
  year      = YEARS,
  stringsAsFactors = FALSE
)
missing_combos <- anti_join(expected, nir_targets, by = c("crf_group", "year"))
if (nrow(missing_combos) > 0L) {
  cat("  WARNING: missing group-year combinations (E_NIR_kt will be treated as 0):\n")
  print(missing_combos)
  # Fill with zeros
  nir_targets <- bind_rows(
    nir_targets,
    missing_combos %>% mutate(E_NIR_kt = 0)
  )
}

cat("  Total rows:", nrow(nir_targets), "\n")
cat("  Groups    :", n_distinct(nir_targets$crf_group), "\n")
cat("  Years     :", n_distinct(nir_targets$year), "\n\n")

# Print summary table
wide <- nir_targets %>%
  filter(year %in% c(2005, 2010, 2015, 2022)) %>%
  tidyr::pivot_wider(names_from = year, values_from = E_NIR_kt) %>%
  arrange(crf_group)
cat("Sample E_NIR_kt (selected years):\n")
print(wide)
cat("\n")


# =============================================================================
# SECTION 6: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "nir_calibration_targets.RData")
save(nir_targets, file = OUT_PATH)

cat("══════════════════════════════════════════════════════════\n")
cat("Saved:", OUT_PATH, "\n")
cat("  Rows:", nrow(nir_targets), "\n")
cat("══════════════════════════════════════════════════════════\n")
