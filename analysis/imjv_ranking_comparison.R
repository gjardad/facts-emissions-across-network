###############################################################################
# analysis/imjv_ranking_comparison.R
#
# PURPOSE
#   For the 56 IMJV TP firm-years (non-ETS firms correctly classified as
#   emitters), compare their q_i ranking to their revenue ranking within
#   CRF-year, NACE 2d-year, and NACE 5d-year cells.
#
# INPUT
#   {RAW_DATA}/IMJV/imjv_co2_lucht.tsv
#   {RAW_DATA}/IMJV/crosswalk/imjv_cbb_vat_ano.csv
#   {RAW_DATA}/IMJV/crosswalk/non_ets_cbbs.csv
#   {PROC_DATA}/allocation_glo_balanced/alloc_YYYY.RData
#   {PROC_DATA}/extensive_margin_calibration.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData  (FULL, not downsampled)
#   {PROC_DATA}/deployment_panel.RData
#   {REPO_DIR}/preprocess/crosswalks/nace_crf_crosswalk.csv
#
# OUTPUT
#   Printed summary + {PROC_DATA}/imjv_ranking_comparison.csv
#
# RUNS ON: RMD (needs full annual accounts and deployment panel)
###############################################################################

# -- Paths -------------------------------------------------------------------
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

cat("===================================================================\n")
cat("  IMJV ranking comparison: q_i vs revenue\n")
cat("===================================================================\n\n")


# =============================================================================
# 1. Load logistic calibration
# =============================================================================
load(file.path(PROC_DATA, "extensive_margin_calibration.RData"))
WEIGHT_SCHEME <- "balanced"
q_fn   <- ext_q_function[[WEIGHT_SCHEME]]
q_star <- q_star_pooled[[WEIGHT_SCHEME]]
cat("q_star =", round(q_star, 4), "\n\n")


# =============================================================================
# 2. Load IMJV + crosswalk → non-ETS firm-years
# =============================================================================
cat("-- Loading IMJV data ------------------------------------------------\n")

imjv <- read.delim(file.path(RAW_DATA, "IMJV", "imjv_co2_lucht.tsv"),
                    stringsAsFactors = FALSE,
                    colClasses = c(cbb_number = "character", cbb_current = "character"))
xwalk <- read.csv(file.path(RAW_DATA, "IMJV", "crosswalk", "imjv_cbb_vat_ano.csv"),
                   stringsAsFactors = FALSE,
                   colClasses = c(cbb = "character", vat_ano = "character"))
imjv <- imjv %>%
  left_join(xwalk %>% select(cbb, vat_ano), by = c("cbb_number" = "cbb")) %>%
  filter(!is.na(vat_ano))

imjv_vy <- imjv %>%
  group_by(vat_ano, year) %>%
  summarise(imjv_co2_t = sum(emission_kg, na.rm = TRUE) / 1000, .groups = "drop")

vat_names <- imjv %>%
  group_by(vat_ano) %>%
  summarise(firm_name = first(firm_name), .groups = "drop")

non_ets <- read.csv(file.path(RAW_DATA, "IMJV", "crosswalk", "non_ets_cbbs.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(cbb = "character", vat_ano = "character"))
non_ets_vats <- unique(non_ets$vat_ano)


# =============================================================================
# 3. Load allocation → get q_i for all imputed firms
# =============================================================================
cat("-- Loading allocation files -----------------------------------------\n")

alloc_list <- list()
for (yr in 2005:2021) {
  load(file.path(PROC_DATA, "allocation_glo_balanced", sprintf("alloc_%d.RData", yr)))
  alloc_list[[as.character(yr)]] <- year_firms
  rm(year_firms)
}
alloc <- bind_rows(alloc_list)
rm(alloc_list)

alloc_imp <- alloc %>%
  filter(source == "imputed") %>%
  mutate(q_i = q_fn(p_i, proxy_mean_i))

cat("  Imputed firm-years:", nrow(alloc_imp), "\n")


# =============================================================================
# 4. Load NACE codes from deployment panel
# =============================================================================
cat("-- Loading deployment panel (for NACE codes) -----------------------\n")

load(file.path(PROC_DATA, "deployment_panel.RData"))
deploy_nace <- deployment_panel %>%
  select(vat, year, nace5d) %>%
  mutate(nace2d = substr(nace5d, 1, 2))
rm(deployment_panel)
cat("  Deployment panel rows:", nrow(deploy_nace), "\n")

# NACE-CRF crosswalk
nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(nace2d = "character")
) %>% select(nace2d, crf_group)


# =============================================================================
# 5. Load revenue from annual accounts (FULL version on RMD)
# =============================================================================
cat("-- Loading annual accounts (for revenue) ---------------------------\n")

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  select(vat, year, revenue, nace5d) %>%
  mutate(nace2d = substr(nace5d, 1, 2))
rm(df_annual_accounts_selected_sample_key_variables)
cat("  Annual accounts rows:", nrow(accounts), "\n\n")


# =============================================================================
# 6. Join NACE and revenue to imputed firms, compute rankings
# =============================================================================
cat("-- Computing rankings -----------------------------------------------\n")

# Add NACE to imputed firms
alloc_imp <- alloc_imp %>%
  left_join(deploy_nace, by = c("vat", "year"))

# Add revenue
alloc_imp <- alloc_imp %>%
  left_join(accounts %>% select(vat, year, revenue), by = c("vat", "year"))

# Compute within-cell rankings for q_i and revenue
# CRF-year
alloc_imp <- alloc_imp %>%
  group_by(crf_group, year) %>%
  mutate(
    n_crf_yr       = n(),
    q_rank_crf     = rank(-q_i, ties.method = "average"),
    q_pctile_crf   = 1 - (q_rank_crf - 1) / n_crf_yr,
    rev_rank_crf   = rank(-revenue, ties.method = "average", na.last = "keep"),
    n_rev_crf      = sum(!is.na(revenue)),
    rev_pctile_crf = ifelse(!is.na(revenue),
                            1 - (rank(-revenue[!is.na(revenue)], ties.method = "average")[
                              cumsum(!is.na(revenue))[seq_len(n())] * as.integer(!is.na(revenue))
                            ] - 1) / n_rev_crf,
                            NA_real_)
  ) %>%
  ungroup()

# Simpler revenue percentile: rank among non-NA revenue firms
alloc_imp <- alloc_imp %>%
  group_by(crf_group, year) %>%
  mutate(
    rev_pctile_crf = {
      r <- rep(NA_real_, n())
      ok <- !is.na(revenue)
      if (sum(ok) > 1) {
        ranks <- rank(-revenue[ok], ties.method = "average")
        r[ok] <- 1 - (ranks - 1) / sum(ok)
      }
      r
    }
  ) %>%
  ungroup()

# NACE 2d-year
alloc_imp <- alloc_imp %>%
  group_by(nace2d, year) %>%
  mutate(
    n_nace2d_yr       = n(),
    q_pctile_nace2d   = {
      r <- rank(-q_i, ties.method = "average")
      1 - (r - 1) / n()
    },
    rev_pctile_nace2d = {
      r <- rep(NA_real_, n())
      ok <- !is.na(revenue)
      if (sum(ok) > 1) {
        ranks <- rank(-revenue[ok], ties.method = "average")
        r[ok] <- 1 - (ranks - 1) / sum(ok)
      }
      r
    }
  ) %>%
  ungroup()

# NACE 5d-year
alloc_imp <- alloc_imp %>%
  group_by(nace5d, year) %>%
  mutate(
    n_nace5d_yr       = n(),
    q_pctile_nace5d   = {
      r <- rank(-q_i, ties.method = "average")
      1 - (r - 1) / n()
    },
    rev_pctile_nace5d = {
      r <- rep(NA_real_, n())
      ok <- !is.na(revenue)
      if (sum(ok) > 1) {
        ranks <- rank(-revenue[ok], ties.method = "average")
        r[ok] <- 1 - (ranks - 1) / sum(ok)
      }
      r
    }
  ) %>%
  ungroup()


# =============================================================================
# 7. Extract IMJV TP firm-years and report
# =============================================================================

tp <- imjv_vy %>%
  filter(vat_ano %in% non_ets_vats, year >= 2005, year <= 2021, imjv_co2_t > 0) %>%
  inner_join(
    alloc_imp %>% select(vat, year, crf_group, scope1, q_i, revenue,
                          n_crf_yr, q_pctile_crf, rev_pctile_crf,
                          nace2d, n_nace2d_yr, q_pctile_nace2d, rev_pctile_nace2d,
                          nace5d, n_nace5d_yr, q_pctile_nace5d, rev_pctile_nace5d),
    by = c("vat_ano" = "vat", "year")
  ) %>%
  filter(scope1 > 0) %>%
  left_join(vat_names, by = "vat_ano") %>%
  mutate(
    gap_crf   = q_pctile_crf - rev_pctile_crf,
    gap_nace2d = q_pctile_nace2d - rev_pctile_nace2d,
    gap_nace5d = q_pctile_nace5d - rev_pctile_nace5d
  )

cat("\n===================================================================\n")
cat("  RESULTS\n")
cat("===================================================================\n\n")

cat("TP firm-years:", nrow(tp), "\n")
cat("With non-NA revenue:", sum(!is.na(tp$revenue)), "\n\n")

tp_rev <- tp %>% filter(!is.na(revenue))

cat("--- CRF-year rankings ---\n")
cat("q_i percentile:     median =", round(median(tp_rev$q_pctile_crf), 3), "\n")
cat("Revenue percentile: median =", round(median(tp_rev$rev_pctile_crf), 3), "\n")
cat("Gap (q - rev):      median =", round(median(tp_rev$gap_crf), 3), "\n")
cat("Spearman(q_pct, rev_pct):  ", round(cor(tp_rev$q_pctile_crf, tp_rev$rev_pctile_crf,
                                               method = "spearman"), 3), "\n\n")

cat("--- NACE 2d-year rankings ---\n")
tp_n2 <- tp_rev %>% filter(!is.na(q_pctile_nace2d), !is.na(rev_pctile_nace2d))
cat("n =", nrow(tp_n2), "\n")
cat("q_i percentile:     median =", round(median(tp_n2$q_pctile_nace2d), 3), "\n")
cat("Revenue percentile: median =", round(median(tp_n2$rev_pctile_nace2d), 3), "\n")
cat("Gap (q - rev):      median =", round(median(tp_n2$gap_nace2d), 3), "\n")
cat("Spearman(q_pct, rev_pct):  ", round(cor(tp_n2$q_pctile_nace2d, tp_n2$rev_pctile_nace2d,
                                               method = "spearman"), 3), "\n\n")

cat("--- NACE 5d-year rankings ---\n")
tp_n5 <- tp_rev %>% filter(!is.na(q_pctile_nace5d), !is.na(rev_pctile_nace5d), n_nace5d_yr >= 10)
cat("n =", nrow(tp_n5), "(cells with >= 10 firms)\n")
if (nrow(tp_n5) > 0) {
  cat("q_i percentile:     median =", round(median(tp_n5$q_pctile_nace5d), 3), "\n")
  cat("Revenue percentile: median =", round(median(tp_n5$rev_pctile_nace5d), 3), "\n")
  cat("Gap (q - rev):      median =", round(median(tp_n5$gap_nace5d), 3), "\n")
  cat("Spearman(q_pct, rev_pct):  ", round(cor(tp_n5$q_pctile_nace5d, tp_n5$rev_pctile_nace5d,
                                                 method = "spearman"), 3), "\n\n")
}

cat("--- Is q_i ranking above or below revenue ranking? ---\n")
cat("CRF-year:   q > rev in", sum(tp_rev$gap_crf > 0), "of", nrow(tp_rev), "firm-years (",
    round(100 * mean(tp_rev$gap_crf > 0), 1), "%)\n")
if (nrow(tp_n2) > 0)
  cat("NACE2d-year: q > rev in", sum(tp_n2$gap_nace2d > 0), "of", nrow(tp_n2), "firm-years (",
      round(100 * mean(tp_n2$gap_nace2d > 0), 1), "%)\n")

cat("\n--- Firm-level detail (CRF-year, sorted by gap) ---\n")
cat(sprintf("%-28s %5s %8s %7s %7s %7s %6s\n",
            "Firm", "Year", "CRF", "q_pct", "rev_pct", "gap", "ratio"))
cat(strrep("-", 78), "\n")
tp_detail <- tp_rev %>% arrange(gap_crf)
for (i in seq_len(nrow(tp_detail))) {
  r <- tp_detail[i, ]
  cat(sprintf("%-28s %5d %8s %7.3f %7.3f %7.3f %6.2f\n",
              substr(r$firm_name, 1, 28), r$year, r$crf_group,
              r$q_pctile_crf, r$rev_pctile_crf, r$gap_crf,
              r$scope1 / r$imjv_co2_t))
}

cat("\n--- Cutoff comparison: q_i vs revenue (CRF-year) ---\n")
cat("If we classify as emitter all firms in the top X% of their CRF-year cell,\n")
cat("how many of the 56 IMJV TP firm-years would we detect?\n\n")

cutoffs <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
cat(sprintf("%8s %12s %12s\n", "Top X%", "q_i", "Revenue"))
cat(strrep("-", 36), "\n")
for (x in cutoffs) {
  n_q   <- sum(tp$q_pctile_crf >= (1 - x))
  n_rev <- sum(tp_rev$rev_pctile_crf >= (1 - x), na.rm = TRUE)
  denom_q   <- nrow(tp)
  denom_rev <- nrow(tp_rev)
  cat(sprintf("%7.0f%% %5d / %-4d %5d / %-4d\n",
              x * 100, n_q, denom_q, n_rev, denom_rev))
}
cat("\n")

# Save
out <- tp %>%
  select(firm_name, vat_ano, year, crf_group, nace2d, nace5d, imjv_co2_t, scope1,
         q_i, revenue, q_pctile_crf, rev_pctile_crf, gap_crf,
         q_pctile_nace2d, rev_pctile_nace2d, gap_nace2d,
         q_pctile_nace5d, rev_pctile_nace5d, gap_nace5d,
         n_crf_yr, n_nace2d_yr, n_nace5d_yr)

write.csv(out, file.path(PROC_DATA, "imjv_ranking_comparison.csv"), row.names = FALSE)
cat("\nSaved:", file.path(PROC_DATA, "imjv_ranking_comparison.csv"), "\n")
