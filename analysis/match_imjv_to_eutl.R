###############################################################################
# analysis/match_imjv_to_eutl.R
#
# PURPOSE
#   Match IMJV firms to EUTL installations to assess overlap.
#   Uses firm name similarity, municipality/city, and NACE codes.
#
# INPUTS
#   {RAW_DATA}/IMJV/imjv_co2_lucht.tsv
#   {RAW_DATA}/EUTL/Oct_2024_version/installation.csv
#
# OUTPUTS
#   {PROC_DATA}/imjv_eutl_match.RData
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


# =============================================================================
# STEP 1: Load data
# =============================================================================
cat("== STEP 1: Loading data ==\n")

imjv <- read.delim(file.path(RAW_DATA, "IMJV", "imjv_co2_lucht.tsv"),
                   stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

# Unique IMJV firms
imjv_firms <- imjv %>%
  mutate(nace2d = as.integer(sub("\\..*", "", nace_code))) %>%
  group_by(cbb_number, firm_name, gemeente, postcode, nace_code, nace2d) %>%
  summarise(
    mean_annual_co2_t = mean(emission_kg / 1000, na.rm = TRUE),
    n_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  # Use most recent name/location if a firm appears with multiple
  group_by(cbb_number) %>%
  slice_max(n_years, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("IMJV unique firms:", nrow(imjv_firms), "\n")

# Load EUTL installations
inst <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "installation.csv"),
                 stringsAsFactors = FALSE)
be_inst <- inst %>%
  filter(country_id == "BE",
         isAircraftOperator != "True",
         isMaritimeOperator != "True")

cat("EUTL Belgian installations:", nrow(be_inst), "\n\n")


# =============================================================================
# STEP 2: Hard-coded matches (verified corporate name changes)
# =============================================================================
# These were identified by the algorithmic matcher as city+NACE candidates
# with zero name similarity, then verified via web search as corporate name
# changes, mergers, or acquisitions. See analysis/IMJV_README.md for sources.

hardcoded <- tribble(
  ~imjv_pattern,                              ~eutl_name,
  "BP CHEMBEL",                               "INEOS Aromatics Belgium",
  "OUDEGEM PAPIER",                           "VPK Paper",
  "CYTEC SURFACE SPECIALTIES",                "ALLNEX Belgium",
  "SOLVIN",                                   "INOVYN BELGIUM",
  "SOLVIC",                                   "BASF Antwerpen - 127c",
  "LATEXCO",                                  "Novaya Belgium",
  "MISA ECO",                                 "Nesar",
  "NORTH EUROPEAN SULFERIC ACID REGENERATION", "Nesar",
  "ORRION CHEMICALS REGEN",                   "Nesar",
  "V.B.G.",                                   "VBG nv",
  "ASFALT PRODUCTIE LIMBURG",                 "APL nv"
)

# EDF LUMINUS / S.P.E. entries: match by cbb establishment suffix to specific
# EUTL installations. All share the same enterprise ID (1764822000) but have
# different establishment suffixes mapping to different physical plants.
edf_hardcoded <- tribble(
  ~cbb_suffix, ~eutl_name,
  "0137",      "EDFL - Centrale Harelbeke",
  "0238",      "Electrabel - Langerbrugge",
  "0339",      "Electrabel - Langerbrugge",
  "0743",      "EDFL - Izegem"
)

cat("Hard-coded matches: ", nrow(hardcoded), " by name + ",
    nrow(edf_hardcoded), " EDF/SPE entries\n\n")


# =============================================================================
# STEP 3: Similarity functions (reused from match_ct_to_eutl.R)
# =============================================================================

bigram_jaccard <- function(a, b) {
  make_bigrams <- function(s) {
    s <- tolower(gsub("[^a-z0-9]", "", s))
    if (nchar(s) < 2) return(character(0))
    substring(s, 1:(nchar(s)-1), 2:nchar(s))
  }
  bg_a <- make_bigrams(a)
  bg_b <- make_bigrams(b)
  if (length(bg_a) == 0 || length(bg_b) == 0) return(0)
  length(intersect(bg_a, bg_b)) / length(union(bg_a, bg_b))
}

STOPWORDS <- c("nv", "sa", "bv", "bvba", "sprl", "srl", "cvba",
               "site", "vestiging", "usine", "plant",
               "de", "le", "la", "les", "het", "van", "du", "des",
               "den", "der", "een", "the", "and", "et", "en")

token_overlap <- function(a, b) {
  tokenize <- function(s) {
    toks <- unlist(strsplit(tolower(s), "[^a-z0-9]+"))
    toks <- toks[nchar(toks) >= 3]
    toks <- toks[!toks %in% STOPWORDS]
    unique(toks)
  }
  toks_a <- tokenize(a)
  toks_b <- tokenize(b)
  if (length(toks_a) == 0 || length(toks_b) == 0) return(0)
  length(intersect(toks_a, toks_b)) / min(length(toks_a), length(toks_b))
}

name_similarity <- function(a, b) {
  max(bigram_jaccard(a, b), token_overlap(a, b))
}

# Substring check: does one name contain a key token from the other?
# Catches cases like TEEPAK → ViskoTeepak where "teepak" is embedded
substring_match <- function(a, b) {
  a_clean <- tolower(gsub("[^a-z0-9]", "", a))
  b_clean <- tolower(gsub("[^a-z0-9]", "", b))
  if (nchar(a_clean) < 4 || nchar(b_clean) < 4) return(FALSE)
  grepl(a_clean, b_clean, fixed = TRUE) || grepl(b_clean, a_clean, fixed = TRUE)
}

# Municipality name similarity
city_match <- function(gemeente, eutl_city) {
  if (is.na(gemeente) || is.na(eutl_city) || gemeente == "" || eutl_city == "") return(FALSE)
  g <- tolower(trimws(gemeente))
  c <- tolower(trimws(eutl_city))
  # Exact match or one contains the other
  g == c || grepl(g, c, fixed = TRUE) || grepl(c, g, fixed = TRUE)
}

# NACE compatibility: match at 4-digit or 2-digit level
nace_match_4d <- function(imjv_nace, eutl_nace) {
  if (is.na(imjv_nace) || is.na(eutl_nace) || eutl_nace == "") return(FALSE)
  as.character(imjv_nace) == as.character(eutl_nace)
}

nace_match_2d <- function(imjv_nace, eutl_nace) {
  if (is.na(imjv_nace) || is.na(eutl_nace) || eutl_nace == "") return(FALSE)
  imjv_2d <- sub("\\..*", "", as.character(imjv_nace))
  eutl_2d <- sub("\\..*", "", as.character(eutl_nace))
  imjv_2d == eutl_2d
}


# =============================================================================
# STEP 4: Match each IMJV firm to best EUTL candidate
# =============================================================================
cat("== STEP 3: Matching IMJV firms to EUTL installations ==\n")

# For each IMJV firm, first check hard-coded matches, then fall back to
# algorithmic name-based matching.
# City+NACE matches are flagged for manual review, not auto-accepted.
#
# Auto-accept tiers:
#   1. name_nace:       name sim >= 0.4 AND NACE 2d match
#   2. name:            name sim >= 0.4
#   3. name_city_nace:  name sim >= 0.2 AND city match AND NACE 2d match
#   4. name_city:       name sim >= 0.2 AND city match
#   5. substr_city_nace: substring match AND city match AND NACE 2d match
#
# Flagged for review (not auto-accepted):
#   city+NACE matches with zero name similarity (possible corporate name changes)

matches <- list()
flagged <- list()

for (i in seq_len(nrow(imjv_firms))) {
  f <- imjv_firms[i, ]

  # --- Check hard-coded matches first ---
  hc_hit <- hardcoded %>% filter(imjv_pattern == f$firm_name)
  edf_hit <- NULL
  if (nrow(hc_hit) == 0) {
    # Check EDF/SPE by cbb suffix
    cbb_str <- as.character(f$cbb_number)
    suffix <- substr(cbb_str, nchar(cbb_str) - 3, nchar(cbb_str))
    edf_hit <- edf_hardcoded %>% filter(cbb_suffix == suffix)
    if (nrow(edf_hit) == 0) edf_hit <- NULL
  }

  if (nrow(hc_hit) > 0 || !is.null(edf_hit)) {
    target_name <- if (nrow(hc_hit) > 0) hc_hit$eutl_name[1] else edf_hit$eutl_name[1]
    hc_idx <- which(be_inst$name == target_name)
    if (length(hc_idx) > 0) {
      hc_idx <- hc_idx[1]
      matches[[i]] <- data.frame(
        cbb_number = f$cbb_number,
        imjv_name = f$firm_name,
        imjv_gemeente = f$gemeente,
        imjv_nace = f$nace_code,
        mean_annual_co2_t = f$mean_annual_co2_t,
        eutl_name = be_inst$name[hc_idx],
        eutl_city = be_inst$city[hc_idx],
        eutl_nace = be_inst$nace_id[hc_idx],
        name_sim = 0,
        substr_match = FALSE,
        city_match = TRUE,
        nace_match_2d = TRUE,
        nace_match_4d = TRUE,
        matched = TRUE,
        match_type = "hardcoded",
        stringsAsFactors = FALSE
      )
      if (i %% 20 == 0) cat("  Processed", i, "of", nrow(imjv_firms), "firms\n")
      next
    }
  }

  # --- Algorithmic matching ---
  best_score <- -1
  best_idx <- NA
  best_city <- FALSE
  best_nace2d <- FALSE
  best_nace4d <- FALSE
  best_name_sim <- 0
  best_substr <- FALSE
  best_tier <- "none"

  # Also track best city+NACE candidate for flagging
  best_cn_idx <- NA
  best_cn_nace4d <- FALSE

  for (j in seq_len(nrow(be_inst))) {
    sim <- name_similarity(f$firm_name, be_inst$name[j])
    cm <- city_match(f$gemeente, be_inst$city[j])
    n2d <- nace_match_2d(f$nace_code, be_inst$nace_id[j])
    n4d <- nace_match_4d(f$nace_code, be_inst$nace_id[j])
    sm <- substring_match(f$firm_name, be_inst$name[j])

    # Auto-accept tiers
    tier <- "none"
    score <- 0

    if (sim >= 0.4 && n2d) {
      tier <- "name_nace"
      score <- 6 + sim
    } else if (sim >= 0.4) {
      tier <- "name"
      score <- 5 + sim
    } else if (sim >= 0.2 && cm && n2d) {
      tier <- "name_city_nace"
      score <- 4 + sim
    } else if (sim >= 0.2 && cm) {
      tier <- "name_city"
      score <- 3 + sim
    } else if (sm && cm && n2d) {
      tier <- "substr_city_nace"
      score <- 2.5
    }

    if (score > best_score && tier != "none") {
      best_score <- score
      best_idx <- j
      best_city <- cm
      best_nace2d <- n2d
      best_nace4d <- n4d
      best_name_sim <- sim
      best_substr <- sm
      best_tier <- tier
    }

    # Track city+NACE candidates for flagging (prefer 4d over 2d)
    if (cm && n2d && sim < 0.2) {
      if (is.na(best_cn_idx) || (n4d && !best_cn_nace4d)) {
        best_cn_idx <- j
        best_cn_nace4d <- n4d
      }
    }
  }

  matched <- best_tier != "none"
  match_type <- best_tier

  matches[[i]] <- data.frame(
    cbb_number = f$cbb_number,
    imjv_name = f$firm_name,
    imjv_gemeente = f$gemeente,
    imjv_nace = f$nace_code,
    mean_annual_co2_t = f$mean_annual_co2_t,
    eutl_name = if (!is.na(best_idx)) be_inst$name[best_idx] else NA,
    eutl_city = if (!is.na(best_idx)) be_inst$city[best_idx] else NA,
    eutl_nace = if (!is.na(best_idx)) be_inst$nace_id[best_idx] else NA,
    name_sim = if (!is.na(best_idx)) best_name_sim else NA,
    substr_match = best_substr,
    city_match = best_city,
    nace_match_2d = best_nace2d,
    nace_match_4d = best_nace4d,
    matched = matched,
    match_type = match_type,
    stringsAsFactors = FALSE
  )

  # Flag unmatched firms that have a city+NACE candidate
  if (!matched && !is.na(best_cn_idx)) {
    flagged[[length(flagged) + 1]] <- data.frame(
      imjv_name = f$firm_name,
      imjv_gemeente = f$gemeente,
      imjv_nace = f$nace_code,
      mean_annual_co2_t = f$mean_annual_co2_t,
      eutl_candidate = be_inst$name[best_cn_idx],
      eutl_city = be_inst$city[best_cn_idx],
      eutl_nace = be_inst$nace_id[best_cn_idx],
      nace_match_4d = best_cn_nace4d,
      stringsAsFactors = FALSE
    )
  }

  if (i %% 20 == 0) cat("  Processed", i, "of", nrow(imjv_firms), "firms\n")
}

imjv_eutl <- bind_rows(matches)
flagged_df <- if (length(flagged) > 0) bind_rows(flagged) else data.frame()
cat("Done.\n\n")


# =============================================================================
# STEP 4: Results
# =============================================================================
cat("== Matching results ==\n\n")
cat("Total IMJV firms:", nrow(imjv_eutl), "\n")
cat("Matched to EUTL:", sum(imjv_eutl$matched), "\n")
cat("Not matched:", sum(!imjv_eutl$matched), "\n\n")

cat("By match type:\n")
print(table(imjv_eutl$match_type))

# Matched firms
cat("\n== Matched IMJV firms (sample) ==\n\n")
matched_df <- imjv_eutl %>%
  filter(matched) %>%
  arrange(desc(mean_annual_co2_t))
print(head(as.data.frame(matched_df %>%
  select(imjv_name, eutl_name, name_sim, city_match, mean_annual_co2_t)), 20),
  row.names = FALSE)

# Unmatched firms
cat("\n== Unmatched IMJV firms ==\n\n")
unmatched_df <- imjv_eutl %>%
  filter(!matched) %>%
  arrange(desc(mean_annual_co2_t))

print(as.data.frame(unmatched_df %>%
  select(imjv_name, imjv_gemeente, imjv_nace, mean_annual_co2_t,
         eutl_name, name_sim)), row.names = FALSE)

cat(sprintf("\n== Summary ==\n"))
cat(sprintf("Matched firms: %d (%.0f t/yr mean CO2)\n",
            sum(imjv_eutl$matched),
            mean(imjv_eutl$mean_annual_co2_t[imjv_eutl$matched])))
cat(sprintf("Unmatched firms: %d (%.0f t/yr mean CO2)\n",
            sum(!imjv_eutl$matched),
            mean(imjv_eutl$mean_annual_co2_t[!imjv_eutl$matched])))
cat(sprintf("Total matched CO2: %.0f t/yr\n",
            sum(imjv_eutl$mean_annual_co2_t[imjv_eutl$matched])))
cat(sprintf("Total unmatched CO2: %.0f t/yr\n",
            sum(imjv_eutl$mean_annual_co2_t[!imjv_eutl$matched])))

# NACE breakdown of unmatched
cat("\n== Unmatched by NACE 2-digit ==\n\n")
unmatched_nace <- unmatched_df %>%
  mutate(nace2d = as.integer(sub("\\..*", "", imjv_nace))) %>%
  group_by(nace2d) %>%
  summarise(
    n_firms = n(),
    total_co2_t = sum(mean_annual_co2_t),
    .groups = "drop"
  ) %>%
  arrange(desc(total_co2_t))
print(as.data.frame(unmatched_nace))

# Flagged for manual review: unmatched firms with city+NACE candidates
if (nrow(flagged_df) > 0) {
  cat("\n== FLAGGED FOR MANUAL REVIEW ==\n")
  cat("These unmatched IMJV firms have a city+NACE candidate in EUTL\n")
  cat("but zero name similarity. Could be corporate name changes or\n")
  cat("different firms in the same city/sector. Verify manually.\n\n")
  flagged_df <- flagged_df %>% arrange(desc(mean_annual_co2_t))
  print(as.data.frame(flagged_df), row.names = FALSE)
} else {
  cat("\nNo flagged cases.\n")
}


# =============================================================================
# STEP 5: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "imjv_eutl_match.RData")
save(imjv_eutl, imjv_firms, file = OUT_PATH)
cat("\nSaved to:", OUT_PATH, "\n")
