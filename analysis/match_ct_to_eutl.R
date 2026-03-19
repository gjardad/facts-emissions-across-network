###############################################################################
# analysis/match_ct_to_eutl.R
#
# PURPOSE
#   Match Climate TRACE (CT) facility-level sources to EUTL installations
#   in Belgium using geographic proximity, name similarity, and sector
#   compatibility. See analysis/README.md for full procedure documentation.
#
#   Matching strategy (four-tier):
#     Tier 1: d <= 1 km, name sim >= 0.2, sector compatible
#     Tier 2: d <= 1 km, sector compatible (name waived)
#     Tier 3: 1 < d <= 5 km, name sim >= 0.2, sector compatible
#     Tier 4: 5 < d <= 10 km, name sim >= 0.4
#
#   One-to-one assignment: greedy, sorted by tier then distance.
#
# INPUTS
#   {RAW_DATA}/Climate TRACE/BEL/DATA/   (source-level CSVs)
#   {RAW_DATA}/EUTL/Oct_2024_version/    (installation.csv)
#
# OUTPUTS
#   {PROC_DATA}/ct_eutl_match.RData
#     Contains:
#       ct              - CT source-year panel (CO2, aggregated)
#       ct_eutl_match   - matching table (source_id -> installation_id)
#       be_inst         - EUTL Belgian installations
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

# -- Parameters --------------------------------------------------------------
GEO_MAX_KM      <- 10.0
GEO_TIGHT_KM    <- 1.0
GEO_MED_KM      <- 5.0
NAME_THRESH_LOW  <- 0.2
NAME_THRESH_HIGH <- 0.4


# =============================================================================
# STEP 1: Load Climate TRACE Belgium data
# =============================================================================
cat("== STEP 1: Loading Climate TRACE data ==\n")

ct_dir <- file.path(RAW_DATA, "Climate TRACE", "BEL", "DATA")
ct_files <- list.files(ct_dir, pattern = "emissions_sources_v.*\\.csv$",
                       recursive = TRUE, full.names = TRUE)
ct_files <- ct_files[!grepl("confidence", ct_files)]
cat("Found", length(ct_files), "CT source-level CSV files\n")

needed_cols <- c("source_id", "source_name", "source_type", "lat", "lon",
                 "start_time", "gas", "emissions_quantity",
                 "sector", "subsector")

ct_raw <- lapply(ct_files, function(f) {
  tryCatch({
    df <- read.csv(f, stringsAsFactors = FALSE)
    present <- intersect(needed_cols, names(df))
    if (!all(c("source_id", "lat", "lon", "start_time", "gas",
               "emissions_quantity") %in% present)) return(NULL)
    df <- df[, present, drop = FALSE]
    df$emissions_quantity <- as.numeric(df$emissions_quantity)
    df$lat  <- as.numeric(df$lat)
    df$lon  <- as.numeric(df$lon)
    df$year <- as.integer(substr(df$start_time, 1, 4))
    df$start_time <- NULL
    df
  }, error = function(e) NULL)
})
ct_raw <- bind_rows(ct_raw[!sapply(ct_raw, is.null)])
cat("CT raw rows:", nrow(ct_raw), "\n")

ct <- ct_raw %>%
  filter(gas == "co2") %>%
  group_by(source_id, source_name, source_type, lat, lon, sector, subsector, year) %>%
  summarise(ct_emissions = sum(emissions_quantity, na.rm = TRUE), .groups = "drop") %>%
  filter(ct_emissions > 0)

cat("CT source-years (CO2 > 0):", nrow(ct), "\n")
cat("CT unique sources:", n_distinct(ct$source_id), "\n")
cat("CT year range:", range(ct$year), "\n\n")


# =============================================================================
# STEP 2: Load EUTL installations
# =============================================================================
cat("== STEP 2: Loading EUTL installations ==\n")

inst <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "installation.csv"),
                 stringsAsFactors = FALSE)
be_inst <- inst %>%
  filter(country_id == "BE",
         isAircraftOperator != "True",
         isMaritimeOperator != "True") %>%
  mutate(
    lat_eutl = ifelse(!is.na(latitudeGoogle), latitudeGoogle, latitudeEutl),
    lon_eutl = ifelse(!is.na(longitudeGoogle), longitudeGoogle, longitudeEutl)
  ) %>%
  filter(!is.na(lat_eutl), !is.na(lon_eutl))

cat("EUTL Belgian installations (with coordinates):", nrow(be_inst), "\n\n")


# =============================================================================
# STEP 3: Similarity functions
# =============================================================================

# Haversine distance in km
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

# Metric A: Bigram Jaccard similarity
# J(a, b) = |B(a) ∩ B(b)| / |B(a) ∪ B(b)|
# where B(s) = set of consecutive 2-char substrings of lowercase alphanumeric s
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

# Metric B: Token overlap score
# Tokenize on whitespace/punctuation, lowercase, remove stopwords,
# keep tokens with nchar >= 3. Compute:
# S_token = |shared tokens| / min(|tokens_a|, |tokens_b|)
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

# Combined name similarity: S_name = max(S_bigram, S_token)
name_similarity <- function(a, b) {
  max(bigram_jaccard(a, b), token_overlap(a, b))
}


# =============================================================================
# STEP 4: Sector compatibility map
# =============================================================================
# Maps CT subsectors to compatible EUTL activity_id values.
# activity_id reference (from activity_type.csv):
#   1  = Combustion installations > 20 MW
#   20 = Combustion of fuels
#   21 = Refining of mineral oil
#   24 = Production of pig iron or steel
#   25 = Production or processing of ferrous metals
#   26 = Production of primary aluminium
#   27 = Production of secondary aluminium
#   28 = Production or processing of non-ferrous metals
#   29 = Production of cement clinker
#   30 = Production of lime
#   31 = Manufacture of glass
#   32 = Manufacture of ceramics
#   35 = Production of pulp
#   36 = Production of paper or cardboard
#   38 = Production of nitric acid
#   39 = Production of adipic acid
#   41 = Production of ammonia
#   42 = Production of bulk chemicals
#   43 = Production of hydrogen and synthesis gas

SECTOR_COMPAT <- list(
  "cement"                      = c(29),
  "iron-and-steel"              = c(24, 25),
  "glass"                       = c(31),
  "oil-and-gas-refining"        = c(21),
  "petrochemical-steam-cracking" = c(42, 43),
  "pulp-and-paper"              = c(35, 36),
  "lime"                        = c(30),
  "chemicals"                   = c(38, 39, 41, 42, 43),
  "other-metals"                = c(26, 27, 28),
  # Generic combustion sectors: compatible with any combustion installation
  # electricity-generation also compatible with refining (21) and chemicals
  # (42, 43) because CHP units in industrial complexes (e.g., Antwerp port)
  # are often registered under the co-located facility's EUTL entry
  "electricity-generation"      = c(1, 20, 21, 42, 43),
  "food-beverage-tobacco"       = c(1, 20),
  "textiles-leather-apparel"    = c(1, 20)
)

# For subsectors not in the map, also allow generic combustion
DEFAULT_COMPAT <- c(1, 20)

is_sector_compatible <- function(ct_subsector, eutl_activity_id) {
  compat_ids <- SECTOR_COMPAT[[ct_subsector]]
  if (is.null(compat_ids)) compat_ids <- DEFAULT_COMPAT
  eutl_activity_id %in% compat_ids
}


# =============================================================================
# STEP 5: Build all candidate pairs within 10 km
# =============================================================================
cat("== STEP 3: Building candidate pairs and matching ==\n")

ct_sources <- ct %>%
  distinct(source_id, source_name, lat, lon, sector, subsector) %>%
  filter(!is.na(lat), !is.na(lon))

candidates <- list()
for (i in seq_len(nrow(ct_sources))) {
  s <- ct_sources[i, ]
  dists <- haversine_km(s$lat, s$lon, be_inst$lat_eutl, be_inst$lon_eutl)

  within_range <- which(dists <= GEO_MAX_KM)
  if (length(within_range) == 0) next

  for (j in within_range) {
    d <- dists[j]
    s_name <- name_similarity(s$source_name, be_inst$name[j])
    s_bigram <- bigram_jaccard(s$source_name, be_inst$name[j])
    s_token <- token_overlap(s$source_name, be_inst$name[j])
    compat <- is_sector_compatible(s$subsector, be_inst$activity_id[j])

    # Determine tier
    if (d <= GEO_TIGHT_KM && s_name >= NAME_THRESH_LOW && compat) {
      tier <- 1L
    } else if (d <= GEO_TIGHT_KM && compat) {
      tier <- 2L
    } else if (d <= GEO_MED_KM && s_name >= NAME_THRESH_LOW && compat) {
      tier <- 3L
    } else if (d <= GEO_MAX_KM && s_name >= NAME_THRESH_HIGH) {
      tier <- 4L
    } else {
      next
    }

    candidates[[length(candidates) + 1]] <- data.frame(
      source_id       = s$source_id,
      ct_name         = s$source_name,
      ct_subsector    = s$subsector,
      installation_id = be_inst$id[j],
      eutl_name       = be_inst$name[j],
      eutl_activity   = be_inst$activity_id[j],
      dist_km         = round(d, 4),
      sim_bigram      = round(s_bigram, 4),
      sim_token       = round(s_token, 4),
      sim_name        = round(s_name, 4),
      sector_compat   = compat,
      tier            = tier,
      stringsAsFactors = FALSE
    )
  }
}

all_candidates <- bind_rows(candidates)
cat("Total candidate pairs:", nrow(all_candidates), "\n")
cat("By tier:\n")
print(table(all_candidates$tier))


# =============================================================================
# STEP 6: One-to-one greedy assignment
# =============================================================================
cat("\n== STEP 4: One-to-one greedy assignment ==\n")

# Sort by tier (ascending = better), then distance (ascending = closer)
all_candidates <- all_candidates %>%
  arrange(tier, dist_km)

claimed_ct   <- character(0)
claimed_eutl <- integer(0)
matches      <- list()

for (i in seq_len(nrow(all_candidates))) {
  row <- all_candidates[i, ]
  if (row$source_id %in% claimed_ct) next
  if (row$installation_id %in% claimed_eutl) next

  claimed_ct   <- c(claimed_ct, row$source_id)
  claimed_eutl <- c(claimed_eutl, row$installation_id)
  matches[[length(matches) + 1]] <- row
}

ct_eutl_match <- bind_rows(matches)

# Label match type for readability
ct_eutl_match <- ct_eutl_match %>%
  mutate(match_type = case_when(
    tier == 1L ~ "geo_name_sector",
    tier == 2L ~ "geo_sector",
    tier == 3L ~ "med_name_sector",
    tier == 4L ~ "extended_name"
  ))


# =============================================================================
# STEP 7: Diagnostics
# =============================================================================
cat("\n-- Matching results --\n")
cat("CT sources total:", nrow(ct_sources), "\n")
cat("CT sources matched:", nrow(ct_eutl_match), "\n")
cat("CT sources unmatched:", nrow(ct_sources) - nrow(ct_eutl_match), "\n\n")

cat("By tier / match type:\n")
print(table(ct_eutl_match$match_type))

cat("\nDistance distribution (km):\n")
print(summary(ct_eutl_match$dist_km))

cat("\nName similarity distribution (combined):\n")
print(summary(ct_eutl_match$sim_name))

cat("\nBigram Jaccard distribution:\n")
print(summary(ct_eutl_match$sim_bigram))

cat("\nToken overlap distribution:\n")
print(summary(ct_eutl_match$sim_token))

# Show matches by CT sector
cat("\nMatched CT sectors:\n")
matched_sectors <- ct %>%
  filter(source_id %in% ct_eutl_match$source_id) %>%
  distinct(source_id, sector, subsector) %>%
  count(sector, subsector, name = "n_sources")
print(as.data.frame(matched_sectors))

cat("\nUnmatched CT sectors:\n")
unmatched_sectors <- ct %>%
  filter(!source_id %in% ct_eutl_match$source_id) %>%
  distinct(source_id, sector, subsector) %>%
  count(sector, subsector, name = "n_sources")
print(as.data.frame(unmatched_sectors))

# Flag ambiguous cases: CT sources with multiple EUTL candidates within 2 km
ambiguous <- all_candidates %>%
  group_by(source_id) %>%
  filter(n() > 1, min(dist_km) <= 2) %>%
  arrange(source_id, tier, dist_km)

if (nrow(ambiguous) > 0) {
  cat("\n-- Ambiguous cases (multiple candidates within 2 km) --\n")
  print(as.data.frame(ambiguous %>%
    select(ct_name, eutl_name, dist_km, sim_name, tier)), row.names = FALSE)
}


# =============================================================================
# STEP 8: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "ct_eutl_match.RData")
save(ct, ct_eutl_match, be_inst, file = OUT_PATH)
cat("\nSaved to:", OUT_PATH, "\n")
cat("  ct:", nrow(ct), "source-years\n")
cat("  ct_eutl_match:", nrow(ct_eutl_match), "matched sources\n")
cat("  be_inst:", nrow(be_inst), "EUTL installations\n")
