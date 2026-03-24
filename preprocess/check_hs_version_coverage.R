###############################################################################
# preprocess/check_hs_version_coverage.R
#
# PURPOSE
#   Assess the impact of using HS 2002-based ISIC mappings for CN8 codes that
#   may come from later HS revisions (2007, 2012). Two concerns:
#
#   1. NEW CODES: HS6 codes in H3/H4 that don't exist in H2. These would go
#      unmapped in our HS2002-only crosswalk.
#
#   2. REDEFINED CODES: HS6 codes that exist in both H2 and H3/H4 but map to
#      different ISIC Rev 3 codes (meaning the product scope changed, so
#      our H2-based EXIOBASE assignment could be wrong).
#
# INPUT
#   {CROSS_DIR}/hs2_isic3_wits.csv   (HS 2002 -> ISIC Rev 3, 5224 rows)
#   {CROSS_DIR}/hs3_isic3_wits.csv   (HS 2007 -> ISIC Rev 3, 5053 rows)
#   {CROSS_DIR}/hs4_isic3_wits.csv   (HS 2012 -> ISIC Rev 3, 5177 rows)
#
# OUTPUT
#   Console report only (no files written)
#
# RUNS ON: local 1 or RMD
###############################################################################

source("paths.R")
CROSS_DIR <- file.path(REPO_DIR, "preprocess", "crosswalks")

cat("==============================================================\n")
cat("  CHECK HS VERSION COVERAGE\n")
cat("==============================================================\n\n")

# -- Read all three concordances -----------------------------------------------
read_wits <- function(file, label) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  names(df) <- c("hs6", "hs_desc", "isicrev3_4d", "isic_desc")
  df$hs6         <- sprintf("%06d", as.integer(df$hs6))
  df$isicrev3_4d <- sprintf("%04d", as.integer(df$isicrev3_4d))
  cat(sprintf("  %s: %d codes\n", label, nrow(df)))
  df
}

h2 <- read_wits(file.path(CROSS_DIR, "hs2_isic3_wits.csv"), "H2 (HS 2002)")
h3 <- read_wits(file.path(CROSS_DIR, "hs3_isic3_wits.csv"), "H3 (HS 2007)")
h4 <- read_wits(file.path(CROSS_DIR, "hs4_isic3_wits.csv"), "H4 (HS 2012)")
cat("\n")


# =============================================================================
# Problem 1: NEW CODES — in H3/H4 but not in H2
# =============================================================================
cat("--------------------------------------------------------------\n")
cat("  PROBLEM 1: NEW HS6 CODES (not in H2)\n")
cat("--------------------------------------------------------------\n\n")

new_in_h3 <- setdiff(h3$hs6, h2$hs6)
new_in_h4 <- setdiff(h4$hs6, h2$hs6)
# Codes in H4 that are also not in H3 (truly new in HS 2012)
new_in_h4_only <- setdiff(h4$hs6, union(h2$hs6, h3$hs6))

cat(sprintf("  New in H3 (not in H2): %d codes\n", length(new_in_h3)))
cat(sprintf("  New in H4 (not in H2): %d codes\n", length(new_in_h4)))
cat(sprintf("  New in H4 (not in H2 or H3): %d codes\n", length(new_in_h4_only)))
cat(sprintf("  Total unique new codes (H3 or H4, not in H2): %d codes\n\n",
            length(union(new_in_h3, new_in_h4))))

# What ISIC codes do these new codes map to?
new_h3_isic <- unique(h3[h3$hs6 %in% new_in_h3, "isicrev3_4d"])
new_h4_isic <- unique(h4[h4$hs6 %in% new_in_h4, "isicrev3_4d"])
cat("  ISIC codes for new H3 codes:", paste(sort(new_h3_isic), collapse = ", "), "\n")
cat("  ISIC codes for new H4 codes:", paste(sort(new_h4_isic), collapse = ", "), "\n")

# Are these ISIC codes already covered by our H2-based crosswalk?
h2_isic <- unique(h2$isicrev3_4d)
new_isic <- setdiff(union(new_h3_isic, new_h4_isic), h2_isic)
if (length(new_isic) > 0) {
  cat("\n  WARNING: New ISIC codes not covered by H2:\n")
  cat("  ", paste(new_isic, collapse = ", "), "\n")
} else {
  cat("\n  All ISIC codes for new HS6 codes are already in H2 coverage.\n")
  cat("  -> New codes map to same ISIC industries as existing H2 codes.\n")
  cat("  -> Adding H3/H4 concordances would recover these codes into the\n")
  cat("     same EXIOBASE industries our crosswalk already covers.\n")
}

# Show some examples of new codes
cat("\n  Examples of new H3 codes (first 10):\n")
for (code in head(sort(new_in_h3), 10)) {
  row <- h3[h3$hs6 == code, ][1, ]
  cat(sprintf("    %s: %s -> ISIC %s\n", code, row$hs_desc, row$isicrev3_4d))
}

cat("\n  Examples of new H4 codes (first 10):\n")
for (code in head(sort(new_in_h4), 10)) {
  row <- h4[h4$hs6 == code, ][1, ]
  cat(sprintf("    %s: %s -> ISIC %s\n", code, row$hs_desc, row$isicrev3_4d))
}


# =============================================================================
# Problem 2: REDEFINED CODES — same HS6, different ISIC
# =============================================================================
cat("\n--------------------------------------------------------------\n")
cat("  PROBLEM 2: REDEFINED HS6 CODES (same code, different ISIC)\n")
cat("--------------------------------------------------------------\n\n")

# Compare H2 vs H3 for codes that exist in both
shared_h2_h3 <- intersect(h2$hs6, h3$hs6)
h2_map <- h2[, c("hs6", "isicrev3_4d")]
h3_map <- h3[, c("hs6", "isicrev3_4d")]
names(h2_map)[2] <- "isic_h2"
names(h3_map)[2] <- "isic_h3"
comp_h2h3 <- merge(h2_map, h3_map, by = "hs6")
diff_h2h3 <- comp_h2h3[comp_h2h3$isic_h2 != comp_h2h3$isic_h3, ]

cat(sprintf("  H2 vs H3: %d shared codes, %d with different ISIC (%.1f%%)\n",
            length(shared_h2_h3), nrow(diff_h2h3),
            100 * nrow(diff_h2h3) / length(shared_h2_h3)))

# Compare H2 vs H4
shared_h2_h4 <- intersect(h2$hs6, h4$hs6)
h4_map <- h4[, c("hs6", "isicrev3_4d")]
names(h4_map)[2] <- "isic_h4"
comp_h2h4 <- merge(h2_map, h4_map, by = "hs6")
diff_h2h4 <- comp_h2h4[comp_h2h4$isic_h2 != comp_h2h4$isic_h4, ]

cat(sprintf("  H2 vs H4: %d shared codes, %d with different ISIC (%.1f%%)\n",
            length(shared_h2_h4), nrow(diff_h2h4),
            100 * nrow(diff_h2h4) / length(shared_h2_h4)))

# For the redefined codes, check if the NACE division (first 2 digits) also changed
# If the division didn't change, the EXIOBASE mapping is likely the same
if (nrow(diff_h2h3) > 0) {
  diff_h2h3$div_h2 <- substr(diff_h2h3$isic_h2, 1, 2)
  diff_h2h3$div_h3 <- substr(diff_h2h3$isic_h3, 1, 2)
  div_changed_h3 <- diff_h2h3[diff_h2h3$div_h2 != diff_h2h3$div_h3, ]
  cat(sprintf("\n  H2 vs H3 redefined: %d changed ISIC, %d also changed NACE division\n",
              nrow(diff_h2h3), nrow(div_changed_h3)))
  if (nrow(div_changed_h3) > 0) {
    cat("  Codes with changed NACE division (H2 vs H3):\n")
    for (i in seq_len(min(nrow(div_changed_h3), 20))) {
      r <- div_changed_h3[i, ]
      cat(sprintf("    %s: ISIC %s (div %s) -> %s (div %s)\n",
                  r$hs6, r$isic_h2, r$div_h2, r$isic_h3, r$div_h3))
    }
  }
}

if (nrow(diff_h2h4) > 0) {
  diff_h2h4$div_h2 <- substr(diff_h2h4$isic_h2, 1, 2)
  diff_h2h4$div_h4 <- substr(diff_h2h4$isic_h4, 1, 2)
  div_changed_h4 <- diff_h2h4[diff_h2h4$div_h2 != diff_h2h4$div_h4, ]
  cat(sprintf("  H2 vs H4 redefined: %d changed ISIC, %d also changed NACE division\n",
              nrow(diff_h2h4), nrow(div_changed_h4)))
  if (nrow(div_changed_h4) > 0) {
    cat("  Codes with changed NACE division (H2 vs H4):\n")
    for (i in seq_len(min(nrow(div_changed_h4), 20))) {
      r <- div_changed_h4[i, ]
      cat(sprintf("    %s: ISIC %s (div %s) -> %s (div %s)\n",
                  r$hs6, r$isic_h2, r$div_h2, r$isic_h4, r$div_h4))
    }
  }
}

# Read the existing crosswalk to check if EXIOBASE assignment actually changes
cat("\n--------------------------------------------------------------\n")
cat("  EXIOBASE IMPACT OF REDEFINED CODES\n")
cat("--------------------------------------------------------------\n\n")

crosswalk <- read.csv(file.path(CROSS_DIR, "hs6_exiobase_crosswalk.csv"),
                      stringsAsFactors = FALSE)
isic_exio <- unique(crosswalk[, c("isicrev3_4d", "exio_code")])
isic_exio$isicrev3_4d <- sprintf("%04d", isic_exio$isicrev3_4d)

# For redefined codes: would the EXIOBASE industry actually change?
check_exio_impact <- function(diff_df, h2_col, hx_col, label) {
  if (nrow(diff_df) == 0) {
    cat(sprintf("  %s: no redefined codes\n", label))
    return(invisible(NULL))
  }
  diff_df <- merge(diff_df, isic_exio, by.x = h2_col, by.y = "isicrev3_4d", all.x = TRUE)
  names(diff_df)[names(diff_df) == "exio_code"] <- "exio_h2"

  # Build ISIC->EXIO for the hx version using same hardcoded mapping
  diff_df <- merge(diff_df, isic_exio, by.x = hx_col, by.y = "isicrev3_4d", all.x = TRUE)
  names(diff_df)[names(diff_df) == "exio_code"] <- "exio_hx"

  exio_changed <- diff_df[!is.na(diff_df$exio_h2) & !is.na(diff_df$exio_hx) &
                           diff_df$exio_h2 != diff_df$exio_hx, ]
  cat(sprintf("  %s: %d codes changed ISIC, %d would change EXIOBASE industry\n",
              label, nrow(diff_df), nrow(exio_changed)))
  if (nrow(exio_changed) > 0) {
    cat("  Examples:\n")
    for (i in seq_len(min(nrow(exio_changed), 15))) {
      r <- exio_changed[i, ]
      cat(sprintf("    %s: %s -> %s\n", r$hs6, r$exio_h2, r$exio_hx))
    }
  }
}

check_exio_impact(diff_h2h3, "isic_h2", "isic_h3", "H2 vs H3")
check_exio_impact(diff_h2h4, "isic_h2", "isic_h4", "H2 vs H4")


cat("\n==============================================================\n")
cat("  DONE\n")
cat("==============================================================\n")
