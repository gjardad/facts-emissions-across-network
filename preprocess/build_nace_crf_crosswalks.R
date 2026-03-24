###############################################################################
# preprocess/build_nace_crf_crosswalks.R
#
# PURPOSE
#   Hard-code the two lookup tables used by the NIR calibration pipeline:
#     1. nace_crf_crosswalk.csv      — NACE 2-digit → crf_group (82 rows)
#     2. crf_group_definitions.csv   — crf_group → regex pattern for CRF codes
#
#   The content is domain knowledge derived from:
#     - Annex I to the AEA Manual (2015): CRF–NACE Rev. 2 correspondence
#     - NIR Belgium sub2025, pp. 94–95, 141–142 (CHP treatment)
#   See preprocess/crosswalks/NACE_CRF_CROSSWALK.md for full rationale.
#
# OUTPUT
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#   {REPO_DIR}/preprocess/crosswalks/crf_group_definitions.csv
#
# RUNS ON: local 1 or RMD (no data dependency — pure hard-coded lookup)
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

CROSS_DIR <- file.path(REPO_DIR, "preprocess", "crosswalks")
if (!dir.exists(CROSS_DIR)) dir.create(CROSS_DIR, recursive = TRUE)

cat("══════════════════════════════════════════════════════════\n")
cat("  BUILD NACE–CRF CROSSWALKS\n")
cat("══════════════════════════════════════════════════════════\n\n")


# =============================================================================
# 1. NACE 2-digit → CRF group
# =============================================================================
# Excluded: NACE 34 (does not exist), 40 (does not exist), 44, 48, 54, 57, 67,
#           76, 83 (do not exist), 49–51 (transport — excluded by design),
#           97–99 (households / extraterritorial).
# See NACE_CRF_CROSSWALK.md §Exclusions for rationale.

nace_crf <- rbind(
  # Agriculture, forestry, fishing
  data.frame(nace2d = "01", nace_label = "Crop and animal production",              crf_group = "agriculture"),
  data.frame(nace2d = "02", nace_label = "Forestry and logging",                    crf_group = "agriculture"),
  data.frame(nace2d = "03", nace_label = "Fishing and aquaculture",                 crf_group = "agriculture"),
  # Mining
  data.frame(nace2d = "05", nace_label = "Mining of coal and lignite",              crf_group = "mfg_other"),
  data.frame(nace2d = "06", nace_label = "Extraction of crude petroleum and natural gas", crf_group = "mfg_other"),
  data.frame(nace2d = "07", nace_label = "Mining of metal ores",                   crf_group = "mfg_other"),
  data.frame(nace2d = "08", nace_label = "Other mining and quarrying",              crf_group = "mfg_other"),
  data.frame(nace2d = "09", nace_label = "Mining support service activities",       crf_group = "mfg_other"),
  # Food, beverages, tobacco
  data.frame(nace2d = "10", nace_label = "Manufacture of food products",            crf_group = "food"),
  data.frame(nace2d = "11", nace_label = "Manufacture of beverages",                crf_group = "food"),
  data.frame(nace2d = "12", nace_label = "Manufacture of tobacco products",         crf_group = "food"),
  # Textiles, leather, wood
  data.frame(nace2d = "13", nace_label = "Manufacture of textiles",                 crf_group = "mfg_other"),
  data.frame(nace2d = "14", nace_label = "Manufacture of wearing apparel",          crf_group = "mfg_other"),
  data.frame(nace2d = "15", nace_label = "Manufacture of leather and related products", crf_group = "mfg_other"),
  data.frame(nace2d = "16", nace_label = "Manufacture of wood and wood products",   crf_group = "mfg_other"),
  # Paper and printing
  data.frame(nace2d = "17", nace_label = "Manufacture of paper and paper products", crf_group = "paper"),
  data.frame(nace2d = "18", nace_label = "Printing and reproduction of recorded media", crf_group = "paper"),
  # Refining
  data.frame(nace2d = "19", nace_label = "Manufacture of coke and refined petroleum products", crf_group = "refining"),
  # Chemicals and pharma
  data.frame(nace2d = "20", nace_label = "Manufacture of chemicals and chemical products", crf_group = "chemicals"),
  data.frame(nace2d = "21", nace_label = "Manufacture of basic pharmaceutical products",  crf_group = "chemicals"),
  # Rubber, plastics, non-metallic minerals
  data.frame(nace2d = "22", nace_label = "Manufacture of rubber and plastic products",    crf_group = "mfg_other"),
  data.frame(nace2d = "23", nace_label = "Manufacture of non-metallic mineral products",  crf_group = "minerals"),
  # Metals (C24 + C25 both assigned to metals — see NACE_CRF_CROSSWALK.md Assumption B)
  data.frame(nace2d = "24", nace_label = "Manufacture of basic metals",                   crf_group = "metals"),
  data.frame(nace2d = "25", nace_label = "Manufacture of fabricated metal products",      crf_group = "metals"),
  # Electronics, electrical, machinery, transport equipment, other mfg
  data.frame(nace2d = "26", nace_label = "Manufacture of computer electronic and optical products", crf_group = "mfg_other"),
  data.frame(nace2d = "27", nace_label = "Manufacture of electrical equipment",           crf_group = "mfg_other"),
  data.frame(nace2d = "28", nace_label = "Manufacture of machinery and equipment n.e.c.", crf_group = "mfg_other"),
  data.frame(nace2d = "29", nace_label = "Manufacture of motor vehicles trailers and semi-trailers", crf_group = "mfg_other"),
  data.frame(nace2d = "30", nace_label = "Manufacture of other transport equipment",      crf_group = "mfg_other"),
  data.frame(nace2d = "31", nace_label = "Manufacture of furniture",                      crf_group = "mfg_other"),
  data.frame(nace2d = "32", nace_label = "Other manufacturing",                           crf_group = "mfg_other"),
  data.frame(nace2d = "33", nace_label = "Repair and installation of machinery and equipment", crf_group = "mfg_other"),
  # Energy
  data.frame(nace2d = "35", nace_label = "Electricity gas steam and air conditioning supply", crf_group = "energy"),
  # Water, waste, remediation → commercial
  data.frame(nace2d = "36", nace_label = "Water collection treatment and supply",         crf_group = "commercial"),
  data.frame(nace2d = "37", nace_label = "Sewage",                                        crf_group = "commercial"),
  data.frame(nace2d = "38", nace_label = "Waste collection treatment disposal and recovery", crf_group = "commercial"),
  data.frame(nace2d = "39", nace_label = "Remediation activities and other waste management", crf_group = "commercial"),
  # Construction
  data.frame(nace2d = "41", nace_label = "Construction of buildings",                     crf_group = "mfg_other"),
  data.frame(nace2d = "42", nace_label = "Civil engineering",                             crf_group = "mfg_other"),
  data.frame(nace2d = "43", nace_label = "Specialised construction activities",           crf_group = "mfg_other"),
  # Trade
  data.frame(nace2d = "45", nace_label = "Wholesale and retail trade and repair of motor vehicles", crf_group = "commercial"),
  data.frame(nace2d = "46", nace_label = "Wholesale trade except of motor vehicles",      crf_group = "commercial"),
  data.frame(nace2d = "47", nace_label = "Retail trade except of motor vehicles",         crf_group = "commercial"),
  # Transport support (49–51 excluded; 52–53 included)
  data.frame(nace2d = "52", nace_label = "Warehousing and support activities for transportation", crf_group = "commercial"),
  data.frame(nace2d = "53", nace_label = "Postal and courier activities",                 crf_group = "commercial"),
  # Accommodation and food service
  data.frame(nace2d = "55", nace_label = "Accommodation",                                 crf_group = "commercial"),
  data.frame(nace2d = "56", nace_label = "Food and beverage service activities",          crf_group = "commercial"),
  # Information and communication
  data.frame(nace2d = "58", nace_label = "Publishing activities",                         crf_group = "commercial"),
  data.frame(nace2d = "59", nace_label = "Motion picture video and television programme production", crf_group = "commercial"),
  data.frame(nace2d = "60", nace_label = "Programming and broadcasting activities",       crf_group = "commercial"),
  data.frame(nace2d = "61", nace_label = "Telecommunications",                            crf_group = "commercial"),
  data.frame(nace2d = "62", nace_label = "Computer programming consultancy and related activities", crf_group = "commercial"),
  data.frame(nace2d = "63", nace_label = "Information service activities",                crf_group = "commercial"),
  # Finance and insurance
  data.frame(nace2d = "64", nace_label = "Financial service activities except insurance and pension funding", crf_group = "commercial"),
  data.frame(nace2d = "65", nace_label = "Insurance reinsurance and pension funding",     crf_group = "commercial"),
  data.frame(nace2d = "66", nace_label = "Activities auxiliary to financial services and insurance", crf_group = "commercial"),
  # Real estate
  data.frame(nace2d = "68", nace_label = "Real estate activities",                        crf_group = "commercial"),
  # Professional, scientific and technical
  data.frame(nace2d = "69", nace_label = "Legal and accounting activities",               crf_group = "commercial"),
  data.frame(nace2d = "70", nace_label = "Activities of head offices and management consultancy", crf_group = "commercial"),
  data.frame(nace2d = "71", nace_label = "Architectural and engineering activities and technical testing", crf_group = "commercial"),
  data.frame(nace2d = "72", nace_label = "Scientific research and development",           crf_group = "commercial"),
  data.frame(nace2d = "73", nace_label = "Advertising and market research",               crf_group = "commercial"),
  data.frame(nace2d = "74", nace_label = "Other professional scientific and technical activities", crf_group = "commercial"),
  data.frame(nace2d = "75", nace_label = "Veterinary activities",                         crf_group = "commercial"),
  # Administrative and support
  data.frame(nace2d = "77", nace_label = "Rental and leasing activities",                 crf_group = "commercial"),
  data.frame(nace2d = "78", nace_label = "Employment activities",                         crf_group = "commercial"),
  data.frame(nace2d = "79", nace_label = "Travel agency tour operator and other reservation activities", crf_group = "commercial"),
  data.frame(nace2d = "80", nace_label = "Security and investigation activities",         crf_group = "commercial"),
  data.frame(nace2d = "81", nace_label = "Services to buildings and landscape activities", crf_group = "commercial"),
  data.frame(nace2d = "82", nace_label = "Office administrative support and other business activities", crf_group = "commercial"),
  # Public administration, education, health, social
  data.frame(nace2d = "84", nace_label = "Public administration and defence compulsory social security", crf_group = "commercial"),
  data.frame(nace2d = "85", nace_label = "Education",                                     crf_group = "commercial"),
  data.frame(nace2d = "86", nace_label = "Human health activities",                       crf_group = "commercial"),
  data.frame(nace2d = "87", nace_label = "Residential care activities",                   crf_group = "commercial"),
  data.frame(nace2d = "88", nace_label = "Social work activities without accommodation",  crf_group = "commercial"),
  # Arts, entertainment, recreation, other services
  data.frame(nace2d = "90", nace_label = "Creative arts and entertainment activities",    crf_group = "commercial"),
  data.frame(nace2d = "91", nace_label = "Libraries archives museums and other cultural activities", crf_group = "commercial"),
  data.frame(nace2d = "92", nace_label = "Gambling and betting activities",               crf_group = "commercial"),
  data.frame(nace2d = "93", nace_label = "Sports activities and amusement and recreation activities", crf_group = "commercial"),
  data.frame(nace2d = "94", nace_label = "Activities of membership organisations",        crf_group = "commercial"),
  data.frame(nace2d = "95", nace_label = "Repair of computers and personal and household goods", crf_group = "commercial"),
  data.frame(nace2d = "96", nace_label = "Other personal service activities",             crf_group = "commercial"),
  stringsAsFactors = FALSE
)

cat("nace_crf_crosswalk: ", nrow(nace_crf), "rows\n")
cat("Groups:", paste(sort(unique(nace_crf$crf_group)), collapse = ", "), "\n\n")

# Sanity: no duplicate NACE codes
stopifnot(!any(duplicated(nace_crf$nace2d)))


# =============================================================================
# 2. CRF group definitions (regex patterns for matching compact CRF codes)
# =============================================================================
# Compact CRF codes have dots removed: "1.A.2.g." → "1A2g".
# Each pattern is matched via grepl() against these compact codes.
# See NACE_CRF_CROSSWALK.md for full documentation of scope and assumptions.

crf_defs <- data.frame(
  crf_group   = c("energy", "refining", "metals", "chemicals",
                  "paper", "food", "minerals", "mfg_other",
                  "commercial", "agriculture"),
  crf_pattern = c(
    "^1A1a",           # Public electricity and heat production incl. main-activity CHP
    "^1A1b",           # Petroleum refining (verification only — no deployment firms in NACE 19)
    "^1A1c|^1A2a|^1A2b|^2C",  # Metals: combustion (incl. 1A1c solid fuels absorbed) + all metal-industry processes
    "^1A2c|^2B",       # Chemicals and pharma: combustion + all chemical-industry processes
    "^1A2d",           # Pulp paper and print: combustion only
    "^1A2e",           # Food beverages and tobacco: combustion only
    "^1A2f|^2A",       # Non-metallic minerals: combustion + all mineral-industry processes
    "^1A2g|^2E",       # Other manufacturing and construction: combustion + electronics-industry processes
    "^1A4a",           # Commercial and institutional: stationary combustion and off-road mobile
    "^1A4c"            # Agriculture forestry and fishing: stationary combustion off-road mobile and fishing
  ),
  description = c(
    "Public electricity and heat production incl. main-activity CHP",
    "Petroleum refining (verification only — no deployment firms in NACE 19)",
    "Metals: combustion (incl. 1A1c solid fuels absorbed) + all metal-industry processes",
    "Chemicals and pharma: combustion + all chemical-industry processes",
    "Pulp paper and print: combustion only",
    "Food beverages and tobacco: combustion only",
    "Non-metallic minerals: combustion + all mineral-industry processes",
    "Other manufacturing and construction: combustion + electronics-industry processes",
    "Commercial and institutional: stationary combustion and off-road mobile",
    "Agriculture forestry and fishing: stationary combustion off-road mobile and fishing"
  ),
  stringsAsFactors = FALSE
)

cat("crf_group_definitions:", nrow(crf_defs), "rows\n")
cat("Groups:", paste(crf_defs$crf_group, collapse = ", "), "\n\n")


# =============================================================================
# 3. Write
# =============================================================================
write.csv(nace_crf,  file.path(CROSS_DIR, "nace_crf_crosswalk.csv"),    row.names = FALSE)
write.csv(crf_defs,  file.path(CROSS_DIR, "crf_group_definitions.csv"), row.names = FALSE)

cat("══════════════════════════════════════════════════════════\n")
cat("Saved to:", CROSS_DIR, "\n")
cat("  nace_crf_crosswalk.csv   —", nrow(nace_crf), "rows\n")
cat("  crf_group_definitions.csv —", nrow(crf_defs), "rows\n")
cat("══════════════════════════════════════════════════════════\n")
