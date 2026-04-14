###############################################################################
# analysis/validate_imjv.R
#
# PURPOSE
#   Compare imputed emissions (from b_allocation_glo.R) against IMJV
#   self-reported emissions for non-ETS Flemish firms.
#
#   Reports:
#     - Extensive margin: share of IMJV emitters correctly classified
#     - Intensive margin: APD among correctly classified emitters
#
# INPUT
#   {RAW_DATA}/IMJV/imjv_co2_lucht.tsv          (IMJV CO2 emissions)
#   {RAW_DATA}/IMJV/crosswalk/imjv_cbb_vat_ano.csv  (CBB -> vat_ano)
#   {PROC_DATA}/allocation_glo_balanced/alloc_YYYY.RData
#
# RUNS ON: local 1 or RMD
###############################################################################

# -- Paths -------------------------------------------------------------------
REPO_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) normalizePath(getwd(), winslash = "/")
)
while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
source(file.path(REPO_DIR, "paths.R"))

library(dplyr)

ALLOC_DIR <- file.path(PROC_DATA, "allocation_glo_balanced")

cat("===================================================================\n")
cat("  IMJV vs imputed emissions validation\n")
cat("===================================================================\n\n")


# =============================================================================
# 1. Load IMJV emissions
# =============================================================================
cat("-- Loading IMJV data ------------------------------------------------\n")

imjv <- read.delim(
  file.path(RAW_DATA, "IMJV", "imjv_co2_lucht.tsv"),
  stringsAsFactors = FALSE,
  colClasses = c(cbb_number = "character", cbb_current = "character")
)
cat("  Raw IMJV rows:", nrow(imjv), "\n")

# Load crosswalk: CBB -> vat_ano
xwalk <- read.csv(
  file.path(RAW_DATA, "IMJV", "crosswalk", "imjv_cbb_vat_ano.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(cbb = "character", vat_ano = "character")
)
cat("  Crosswalk rows:", nrow(xwalk), "\n")

# CBBs are already 14-digit character strings; just trim whitespace
imjv$cbb_number <- trimws(imjv$cbb_number)
xwalk$cbb       <- trimws(xwalk$cbb)

# Each CBB appears once in the crosswalk → one-to-many is fine (one cbb → one vat_ano)
imjv <- imjv %>%
  left_join(xwalk %>% select(cbb, vat_ano), by = c("cbb_number" = "cbb"))

n_no_vat <- sum(is.na(imjv$vat_ano))
if (n_no_vat > 0)
  cat("  WARNING:", n_no_vat, "IMJV rows with no vat_ano match\n")

# Aggregate IMJV to (vat_ano, year): sum across sites of the same firm
imjv_firm <- imjv %>%
  filter(!is.na(vat_ano)) %>%
  group_by(vat_ano, year) %>%
  summarise(
    imjv_co2_t  = sum(emission_kg, na.rm = TRUE) / 1000,  # kg -> tonnes
    n_sites     = n_distinct(cbb_number),
    firm_name   = first(firm_name),
    .groups     = "drop"
  )
cat("  IMJV firm-years (vat_ano x year):", nrow(imjv_firm), "\n")
cat("  Unique IMJV firms:", n_distinct(imjv_firm$vat_ano), "\n\n")


# =============================================================================
# 2. Load allocation files and stack
# =============================================================================
cat("-- Loading allocation files -----------------------------------------\n")

alloc_files <- list.files(ALLOC_DIR, pattern = "^alloc_\\d{4}\\.RData$",
                          full.names = TRUE)
cat("  Found", length(alloc_files), "allocation files\n")

alloc_all <- list()
for (f in alloc_files) {
  load(f)  # loads year_firms
  alloc_all[[length(alloc_all) + 1]] <- year_firms
}
alloc <- bind_rows(alloc_all)
rm(alloc_all, year_firms)
cat("  Total allocation firm-years:", nrow(alloc), "\n")
cat("  Sources:", paste(names(table(alloc$source)), collapse = ", "), "\n\n")


# =============================================================================
# 3. Identify non-ETS IMJV firms
# =============================================================================
cat("-- Identifying non-ETS IMJV firms ----------------------------------\n")

# Non-ETS defined as: not matched to any EUTL installation by EITHER
# (a) BvD/KBO identifier matching, or (b) name+address matching.
# The curated list of 38 non-ETS CBBs (= 33 unique firms) was pre-computed
# from the union of both matching methods. See IMJV_README.md for details.
non_ets_cbbs <- read.csv(
  file.path(RAW_DATA, "IMJV", "crosswalk", "non_ets_cbbs.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(cbb = "character", vat_ano = "character")
)
non_ets_vats <- unique(non_ets_cbbs$vat_ano)

imjv_firm <- imjv_firm %>%
  mutate(is_ets = !(vat_ano %in% non_ets_vats))

n_ets     <- n_distinct(imjv_firm$vat_ano[imjv_firm$is_ets])
n_non_ets <- n_distinct(imjv_firm$vat_ano[!imjv_firm$is_ets])
cat("  IMJV firms ETS (by BvD or name matching):", n_ets, "\n")
cat("  IMJV firms non-ETS:", n_non_ets, "\n\n")

# Keep only non-ETS firm-years within the allocation year range
alloc_years <- sort(unique(alloc$year))
imjv_non_ets <- imjv_firm %>%
  filter(!is_ets, year %in% alloc_years)

cat("  Non-ETS IMJV firm-years in allocation year range (",
    min(alloc_years), "-", max(alloc_years), "):", nrow(imjv_non_ets), "\n")
cat("  Unique non-ETS firms in range:", n_distinct(imjv_non_ets$vat_ano), "\n\n")


# =============================================================================
# 4. Match to allocation
# =============================================================================
cat("-- Matching to allocation -------------------------------------------\n")

# For each non-ETS IMJV (vat_ano, year), check if it appears in allocation
# as source == "imputed" with scope1 > 0
alloc_imputed <- alloc %>%
  filter(source == "imputed") %>%
  select(vat, year, scope1, p_i, proxy_mean_i, crf_group, rank_in_cell,
         n_cell, n_imputed_in_cell) %>%
  rename(vat_ano = vat, imputed_co2_t = scope1)

matched <- imjv_non_ets %>%
  left_join(alloc_imputed, by = c("vat_ano", "year"))

# Classification status
matched <- matched %>%
  mutate(
    in_alloc        = !is.na(imputed_co2_t),
    classified_emit = in_alloc & imputed_co2_t > 0,
    imjv_positive   = imjv_co2_t > 0
  )


# =============================================================================
# 5. Extensive margin
# =============================================================================
cat("\n")
cat("===================================================================\n")
cat("  EXTENSIVE MARGIN\n")
cat("===================================================================\n\n")

# Among IMJV firm-years with positive emissions, how many did we classify?
ext <- matched %>% filter(imjv_positive)

tp <- sum(ext$classified_emit)
fn <- sum(!ext$classified_emit)
n_total <- nrow(ext)

cat("  IMJV firm-years with positive CO2 (non-ETS, in alloc years):", n_total, "\n")
cat("  Correctly classified as emitters (TP):", tp, "\n")
cat("  Missed (FN):", fn, "\n")
cat("  Detection rate (TP / total):", round(tp / n_total, 3), "\n\n")

# Also check: among IMJV firm-years with zero emissions, how many did we
# incorrectly classify as emitters?
ext_zero <- matched %>% filter(!imjv_positive)
if (nrow(ext_zero) > 0) {
  fp_in_zero <- sum(ext_zero$classified_emit)
  cat("  IMJV firm-years with zero CO2:", nrow(ext_zero), "\n")
  cat("  Incorrectly classified as emitters:", fp_in_zero, "\n\n")
}

# Breakdown by firm: ever detected vs never detected
firm_ext <- matched %>%
  filter(imjv_positive) %>%
  group_by(vat_ano, firm_name) %>%
  summarise(
    n_years          = n(),
    n_detected       = sum(classified_emit),
    mean_imjv_co2_t  = mean(imjv_co2_t),
    .groups = "drop"
  ) %>%
  mutate(ever_detected = n_detected > 0)

cat("  --- Firm-level summary ---\n")
cat("  Non-ETS IMJV firms with positive CO2:", nrow(firm_ext), "\n")
cat("  Ever detected in at least one year:", sum(firm_ext$ever_detected), "\n")
cat("  Never detected:", sum(!firm_ext$ever_detected), "\n\n")

cat("  Never-detected firms:\n")
never <- firm_ext %>%
  filter(!ever_detected) %>%
  arrange(-mean_imjv_co2_t)
for (i in seq_len(nrow(never))) {
  cat(sprintf("    %-45s  mean CO2: %8.1f t/yr  (%d years)\n",
              never$firm_name[i], never$mean_imjv_co2_t[i], never$n_years[i]))
}

# Diagnostic: proxy values for missed emitters (FN)
cat("\n  --- Proxy diagnostics for missed emitters (FN) ---\n")

# Get all imputed rows from allocation (including those below threshold)
# and the full proxy summary to check p_i and proxy_mean
alloc_all_sources <- alloc %>%
  filter(source == "imputed") %>%
  select(vat, year, scope1, p_i, proxy_mean_i) %>%
  rename(vat_ano = vat)

fn_detail <- matched %>%
  filter(imjv_positive, !classified_emit) %>%
  select(vat_ano, year, firm_name, imjv_co2_t)

# Check if FN firms appear anywhere in the allocation (as imputed)
fn_in_alloc <- fn_detail %>%
  left_join(alloc_all_sources, by = c("vat_ano", "year"))

cat(sprintf("  FN firm-years: %d\n", nrow(fn_in_alloc)))
cat(sprintf("  Of which appear in allocation as 'imputed': %d\n",
            sum(!is.na(fn_in_alloc$scope1))))
cat(sprintf("  Of which do NOT appear in allocation at all: %d\n\n",
            sum(is.na(fn_in_alloc$scope1))))

# For those not in allocation, they might not be in deployment_panel at all
# (no proxy was ever computed). Let's check if they have any proxy data
# by looking at the allocation for ANY year.
fn_vats <- unique(fn_detail$vat_ano)
fn_ever_in_alloc <- alloc %>%
  filter(vat %in% fn_vats) %>%
  distinct(vat, source)
cat("  FN firms ever appearing in allocation (any year, any source):\n")
if (nrow(fn_ever_in_alloc) > 0) {
  print(as.data.frame(table(fn_ever_in_alloc$source)))
} else {
  cat("    None\n")
}

# Summarise by firm: are these firms simply absent from the deployment panel?
fn_firm_summary <- fn_detail %>%
  left_join(alloc_all_sources, by = c("vat_ano", "year")) %>%
  group_by(vat_ano, firm_name) %>%
  summarise(
    n_years       = n(),
    mean_imjv     = mean(imjv_co2_t),
    n_in_alloc    = sum(!is.na(scope1)),
    mean_p_i      = mean(p_i, na.rm = TRUE),
    mean_proxy    = mean(proxy_mean_i, na.rm = TRUE),
    mean_imputed  = mean(scope1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(-mean_imjv)

cat("\n  FN firm-level detail:\n")
cat(sprintf("  %-40s %10s %5s %5s %8s %10s %10s\n",
            "Firm", "IMJV t/yr", "Yrs", "InAl", "mean_p_i", "mean_proxy", "imp_if_any"))
cat("  ", strrep("-", 95), "\n")
for (i in seq_len(nrow(fn_firm_summary))) {
  r <- fn_firm_summary[i, ]
  cat(sprintf("  %-40s %10.1f %5d %5d %8.3f %10.4f %10.1f\n",
              substr(r$firm_name, 1, 40),
              r$mean_imjv, r$n_years, r$n_in_alloc,
              ifelse(is.nan(r$mean_p_i), NA, r$mean_p_i),
              ifelse(is.nan(r$mean_proxy), NA, r$mean_proxy),
              ifelse(is.nan(r$mean_imputed), NA, r$mean_imputed)))
}


# =============================================================================
# 6. Intensive margin
# =============================================================================
cat("\n")
cat("===================================================================\n")
cat("  INTENSIVE MARGIN (among correctly classified emitters)\n")
cat("===================================================================\n\n")

int <- matched %>%
  filter(imjv_positive, classified_emit) %>%
  mutate(
    apd = abs(imputed_co2_t - imjv_co2_t) / ((imputed_co2_t + imjv_co2_t) / 2),
    log_ratio = log(imputed_co2_t / imjv_co2_t)
  )

if (nrow(int) > 0) {
  cat("  Matched firm-years:", nrow(int), "\n")
  cat("  Unique firms:", n_distinct(int$vat_ano), "\n\n")

  cat("  APD (absolute percentage deviation):\n")
  cat("    Mean:   ", round(mean(int$apd), 3), "\n")
  cat("    Median: ", round(median(int$apd), 3), "\n")
  cat("    P25:    ", round(quantile(int$apd, 0.25), 3), "\n")
  cat("    P75:    ", round(quantile(int$apd, 0.75), 3), "\n\n")

  cat("  Log ratio ln(imputed / IMJV):\n")
  cat("    Mean:   ", round(mean(int$log_ratio), 3),
      " (imputed is on average", round(exp(mean(int$log_ratio)), 1), "x IMJV)\n")
  cat("    Median: ", round(median(int$log_ratio), 3),
      " (", round(exp(median(int$log_ratio)), 1), "x)\n\n")

  cat("  Correlation (imputed vs IMJV):\n")
  cat("    Spearman rho:", round(cor(int$imputed_co2_t, int$imjv_co2_t,
                                      method = "spearman"), 3), "\n")
  cat("    Pearson r:   ", round(cor(int$imputed_co2_t, int$imjv_co2_t,
                                      method = "pearson"), 3), "\n\n")

  # Breakdown by firm
  int_firm <- int %>%
    group_by(vat_ano, firm_name) %>%
    summarise(
      n_years      = n(),
      mean_imjv    = mean(imjv_co2_t),
      mean_imputed = mean(imputed_co2_t),
      mean_apd     = mean(apd),
      .groups = "drop"
    ) %>%
    arrange(-mean_imjv)

  cat("  --- Firm-level detail (sorted by mean IMJV emissions) ---\n")
  cat(sprintf("  %-40s %10s %10s %8s %5s\n",
              "Firm", "IMJV t/yr", "Imp t/yr", "APD", "Yrs"))
  cat("  ", strrep("-", 80), "\n")
  for (i in seq_len(nrow(int_firm))) {
    cat(sprintf("  %-40s %10.1f %10.1f %8.3f %5d\n",
                substr(int_firm$firm_name[i], 1, 40),
                int_firm$mean_imjv[i],
                int_firm$mean_imputed[i],
                int_firm$mean_apd[i],
                int_firm$n_years[i]))
  }
} else {
  cat("  No matched firm-years with positive emissions on both sides.\n")
}

cat("\n===================================================================\n")
cat("  Done.\n")
cat("===================================================================\n")
