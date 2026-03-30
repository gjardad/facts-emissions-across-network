###############################################################################
# Search EUTL for waste incinerators under any NACE code
###############################################################################
REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))

inst <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "installation.csv"),
                  stringsAsFactors = FALSE)
be <- inst[inst$country_id == "BE" &
           (inst$isAircraftOperator == "False" | inst$isAircraftOperator == FALSE) &
           (inst$isMaritimeOperator == "False" | inst$isMaritimeOperator == FALSE), ]

# Search for waste/incineration/afval/verbrand/stoom keywords
kw <- "waste|inciner|afval|verbrand|stoom|huisvuil|IVAGO|ISVAG|IMOG|MIROM|IDM|IGEAN|IVM|IVOO|IVBO"
matches <- be[grepl(kw, be$name, ignore.case = TRUE), ]

cat("EUTL installations matching waste/incinerator keywords:\n\n")
for (i in seq_len(nrow(matches))) {
  cat(sprintf("  [%s] act=%s nace=%s | %s\n",
              matches$id[i], matches$activity_id[i],
              matches$nace_id[i], matches$name[i]))
}

# Also: check compliance data for these installations
comp <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "compliance.csv"),
                  stringsAsFactors = FALSE)

# Get Belgian installations with activity 1 or 20 but classified 35.11
# that might be waste incinerators based on their emissions pattern
cat("\n\n── All NACE 35.11 installations — checking for waste-to-energy ──\n")
be_35 <- be[grepl("^35", as.character(be$nace_id)), ]
cat("Total NACE 35.x installations:", nrow(be_35), "\n\n")

# Show names that could be waste-to-energy
wte_kw <- "biostoom|biomass|stoom|waste|afval|verbrand|huisvuil|energie.*duffel|2valorise|sleco"
be_wte <- be_35[grepl(wte_kw, be_35$name, ignore.case = TRUE), ]
cat("Possible waste-to-energy among NACE 35 installations:\n")
for (i in seq_len(nrow(be_wte))) {
  cat(sprintf("  [%s] nace=%s | %s\n", be_wte$id[i], be_wte$nace_id[i], be_wte$name[i]))
}

# Check the account data to find company names behind installations
acct <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "account.csv"),
                  stringsAsFactors = FALSE)
be_acct <- acct[grepl("^BE_", acct$installation_id), ]
cat("\n\n── Account holders for Belgian installations ──\n")
cat("Rows:", nrow(be_acct), "\n")
cat("Columns:", paste(names(be_acct), collapse = ", "), "\n\n")

# Search account holders for waste companies
if ("companyRegistrationNumber" %in% names(be_acct) || "name" %in% names(be_acct)) {
  name_col <- if ("name" %in% names(be_acct)) "name" else "accountHolderName"
  waste_acct <- be_acct[grepl("IVAGO|ISVAG|INTERGEMEENTELIJK|afval|waste|huisvuil|IMOG|MIROM",
                               be_acct[[name_col]], ignore.case = TRUE), ]
  cat(sprintf("Account holders matching waste keywords (in '%s'):\n", name_col))
  if (nrow(waste_acct) > 0) {
    for (i in seq_len(nrow(waste_acct))) {
      cat(sprintf("  inst=%s | %s\n", waste_acct$installation_id[i], waste_acct[[name_col]][i]))
    }
  } else {
    cat("  None found.\n")
    # Try other columns
    cat("\nFirst few account holder names:\n")
    sample_names <- head(unique(be_acct[[name_col]]), 20)
    cat(paste("  ", sample_names), sep = "\n")
  }
}
