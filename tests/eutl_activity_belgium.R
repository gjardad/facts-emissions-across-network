REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
library(dplyr)

inst <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "installation.csv"),
                  stringsAsFactors = FALSE)
act <- read.csv(file.path(RAW_DATA, "EUTL", "Oct_2024_version", "activity_type.csv"),
                 stringsAsFactors = FALSE)

be_inst <- inst %>%
  filter(country_id == "BE",
         isAircraftOperator == "False" | isAircraftOperator == FALSE,
         isMaritimeOperator == "False" | isMaritimeOperator == FALSE) %>%
  left_join(act, by = c("activity_id" = "id")) %>%
  mutate(nace_id = as.character(nace_id),
         nace20_id = as.character(nace20_id))

cat("Belgian installations:", nrow(be_inst), "\n\n")

# Activity 1 = Combustion > 20MW — these are the power/CHP plants
comb <- be_inst %>% filter(activity_id == 1)
cat("── Activity 1 (Combustion > 20MW): ", nrow(comb), "installations ──\n")
cat("NACE codes:\n")
print(table(comb$nace_id, useNA = "ifany"))

# Activity 20 = Combustion of fuels (newer numbering)
comb20 <- be_inst %>% filter(activity_id == 20)
cat("\n── Activity 20 (Combustion of fuels): ", nrow(comb20), "installations ──\n")
cat("NACE codes:\n")
print(table(comb20$nace_id, useNA = "ifany"))

# Search for CHP/cogen in names
chp_kw <- "CHP|cogen|warmte|stoom|WKK|kraftwarm|chaleur|combined.*heat"
chp <- be_inst %>%
  filter(grepl(chp_kw, name, ignore.case = TRUE))
cat("\n── Installations with CHP/cogen keywords in name:", nrow(chp), "──\n")
if (nrow(chp) > 0) {
  for (i in seq_len(nrow(chp))) {
    cat(sprintf("  [%s] act=%d nace=%s name=%s\n",
                chp$id[i], chp$activity_id[i], chp$nace_id[i], chp$name[i]))
  }
}

# Belgian NACE 35 installations
cat("\n── NACE 35.xx installations ──\n")
be_35 <- be_inst %>% filter(grepl("^35", nace_id))
cat("Count:", nrow(be_35), "\n")
for (i in seq_len(nrow(be_35))) {
  cat(sprintf("  [%s] act=%d nace=%s name=%s\n",
              be_35$id[i], be_35$activity_id[i], be_35$nace_id[i], be_35$name[i]))
}

# Belgian NACE 38 installations
cat("\n── NACE 38.xx installations ──\n")
be_38 <- be_inst %>% filter(grepl("^38", nace_id))
cat("Count:", nrow(be_38), "\n")
for (i in seq_len(nrow(be_38))) {
  cat(sprintf("  [%s] act=%d nace=%s name=%s\n",
              be_38$id[i], be_38$activity_id[i], be_38$nace_id[i], be_38$name[i]))
}

cat("\n=== DONE ===\n")
