###############################################################################
# analysis/table_dispersion_comparison_v2.R
#
# PURPOSE
#   Compare within-sector carbon productivity dispersion:
#     - Belgium ETS-only (directly observed)
#     - Belgium all firms (observed + Pareto-imputed)
#     - Lyubich et al. (2018) — US plants, 6-digit NAICS
#     - De Lyon & Dechezlepretre (2025) — EU firms (20+ emp), 3-4 digit NACE
#
# INPUT
#   {PROC_DATA}/firm_year_belgian_euets.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#   {PROC_DATA}/b_loop_scope1_dispersion_pareto.RData
#     (from b_loop_scope1_dispersion_v2.R)
#
# OUTPUT
#   {OUT_DIR}/dispersion_overview.tex          (Panel: ETS / literature / all)
#   {OUT_DIR}/dispersion_by_sector.tex         (Sector-by-sector ETS vs full)
#   {OUT_DIR}/dispersion_comparison_overview.csv
#   {OUT_DIR}/dispersion_ets_vs_full_2d.csv
#   {OUT_DIR}/dispersion_ets_vs_full_5d.csv
#
# RUNS ON: local
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
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr, warn.conflicts = FALSE)

OUT_DIR <- file.path(REPO_DIR, "output", "tables")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

MIN_N <- 3L
YEARS <- 2005:2021

# -- Helpers ------------------------------------------------------------------
pct_ratio <- function(x, p_hi, p_lo) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  q <- quantile(x, c(p_lo, p_hi), names = FALSE)
  if (q[1] == 0) return(NA_real_)
  q[2] / q[1]
}

do_stats <- function(df) {
  log_cp <- log(df$cp)
  data.frame(
    n_firms     = nrow(df),
    cp_p90p10   = pct_ratio(df$cp, 0.9, 0.1),
    cp_log_9010 = diff(quantile(log_cp, c(0.1, 0.9))),
    cp_var_log  = var(log_cp)
  )
}

fmt  <- function(x, digits = 2) {
  ifelse(is.na(x), "---", formatC(round(x, digits), format = "f", digits = digits))
}
fmt1 <- function(x) fmt(x, 1)
fmt0 <- function(x) ifelse(is.na(x), "---", as.character(round(x, 0)))


# =============================================================================
# PART A: ETS-only dispersion (directly from EUTL microdata)
# =============================================================================
cat("-- Loading data ------------------------------------------------------\n")

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS, !is.na(emissions), emissions > 0) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

ets <- eutl %>%
  inner_join(accounts, by = c("vat", "year")) %>%
  filter(revenue > 0, nace2d >= "10", nace2d <= "33") %>%
  mutate(cp = revenue / emissions)

# ETS 2-digit: average across years per sector
ets_2d <- ets %>%
  group_by(nace2d, year) %>%
  filter(n() >= MIN_N) %>%
  group_modify(~do_stats(.x)) %>%
  ungroup() %>%
  group_by(nace2d) %>%
  summarise(
    ets_n_years     = n(),
    ets_n_firms_avg = round(mean(n_firms), 0),
    ets_p90p10      = mean(cp_p90p10, na.rm = TRUE),
    ets_log9010     = mean(cp_log_9010, na.rm = TRUE),
    .groups = "drop"
  )

# ETS 5-digit: aggregated to 2-digit parent
ets_5d <- ets %>%
  group_by(nace5d, year) %>%
  filter(n() >= MIN_N) %>%
  group_modify(~do_stats(.x)) %>%
  ungroup() %>%
  mutate(nace2d = substr(nace5d, 1, 2)) %>%
  group_by(nace2d) %>%
  summarise(
    ets5d_n_sectors = n_distinct(nace5d),
    ets5d_log9010   = mean(cp_log_9010, na.rm = TRUE),
    .groups = "drop"
  )


# =============================================================================
# PART B: Full-sample dispersion (from scope 1 dispersion output)
# =============================================================================
cat("-- Loading scope 1 dispersion results --------------------------------\n")

load(file.path(PROC_DATA, "b_loop_scope1_dispersion_pareto.RData"))
# all_stats2d, all_stats5d (draw-level)

# Full 2-digit: average across draws, then across years
full_2d <- all_stats2d %>%
  filter(nace2d >= "10", nace2d <= "33") %>%
  group_by(nace2d, year) %>%
  summarise(
    n_firms     = mean(n_firms),
    cp_p90p10   = mean(cp_p90p10, na.rm = TRUE),
    cp_log_9010 = mean(cp_p9010_log, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  group_by(nace2d) %>%
  summarise(
    full_n_firms_avg = round(mean(n_firms), 0),
    full_p90p10      = mean(cp_p90p10, na.rm = TRUE),
    full_log9010     = mean(cp_log_9010, na.rm = TRUE),
    .groups = "drop"
  )

# Full 5-digit: aggregated to 2-digit parent
full_5d <- all_stats5d %>%
  mutate(nace2d = substr(nace5d, 1, 2)) %>%
  filter(nace2d >= "10", nace2d <= "33") %>%
  group_by(nace5d, year) %>%
  summarise(
    cp_log_9010 = mean(cp_p9010_log, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(nace2d = substr(nace5d, 1, 2)) %>%
  group_by(nace2d) %>%
  summarise(
    full5d_n_sectors = n_distinct(nace5d),
    full5d_log9010   = mean(cp_log_9010, na.rm = TRUE),
    .groups = "drop"
  )

rm(all_stats2d, all_stats5d, stats2d_summary, stats5d_summary)


# =============================================================================
# TABLE 1 (CSV): Overall summary -- Belgium vs. Literature
# =============================================================================

ets_2d_mean  <- c(p90p10  = mean(ets_2d$ets_p90p10, na.rm = TRUE),
                  log9010 = mean(ets_2d$ets_log9010, na.rm = TRUE))
ets_5d_mean  <- c(log9010 = mean(ets_5d$ets5d_log9010, na.rm = TRUE))
full_2d_mean <- c(p90p10  = mean(full_2d$full_p90p10, na.rm = TRUE),
                  log9010 = mean(full_2d$full_log9010, na.rm = TRUE))
full_5d_mean <- c(log9010 = mean(full_5d$full5d_log9010, na.rm = TRUE))

tab1 <- data.frame(
  Source = c(
    "Belgium ETS-only (NACE 2-digit)",
    "Belgium ETS-only (NACE 5-digit)",
    "",
    "Lyubich et al. (2018)",
    "De Lyon & Dechezlepretre (2025)",
    "",
    "Belgium all firms (NACE 2-digit)",
    "Belgium all firms (NACE 5-digit)"
  ),
  Granularity = c("2-digit", "5-digit", "",
                  "6-digit NAICS", "3-4 digit NACE", "",
                  "2-digit", "5-digit"),
  Unit = c("firm", "firm", "", "plant", "firm (20+ emp.)", "", "firm", "firm"),
  p90_p10 = c(
    round(ets_2d_mean["p90p10"], 1),
    NA,
    NA,
    9.7,
    24,
    NA,
    round(full_2d_mean["p90p10"], 1),
    NA
  ),
  log_9010 = c(
    round(ets_2d_mean["log9010"], 2),
    round(ets_5d_mean["log9010"], 2),
    NA,
    2.27,
    3.18,
    NA,
    round(full_2d_mean["log9010"], 2),
    round(full_5d_mean["log9010"], 2)
  ),
  stringsAsFactors = FALSE
)

cat("\n")
cat("======================================================================\n")
cat("  TABLE 1: Carbon Productivity Dispersion -- Belgium vs. Literature\n")
cat("  Manufacturing sectors only. Mean across sector-years.\n")
cat("======================================================================\n\n")
print(tab1, row.names = FALSE, right = FALSE, na.print = "")
write.csv(tab1, file.path(OUT_DIR, "dispersion_comparison_overview.csv"),
          row.names = FALSE, na = "")


# =============================================================================
# TABLE 1 (LaTeX): dispersion_overview.tex
# =============================================================================

tex1 <- c(
  "\\begin{tabular}{l l l cc}",
  "\\toprule",
  " & Granularity & Unit & $p_{90}/p_{10}$ & $\\ln p_{90} - \\ln p_{10}$ \\\\",
  "\\midrule",
  "\\addlinespace[2pt]",
  "\\multicolumn{5}{l}{\\textit{Panel A: Observed emissions (ETS firms only)}} \\\\",
  "\\addlinespace[2pt]",
  sprintf("Belgium & NACE 2-digit & Firm & %s & %s \\\\",
          fmt1(ets_2d_mean["p90p10"]), fmt(ets_2d_mean["log9010"])),
  sprintf("Belgium & NACE 5-digit & Firm & --- & %s \\\\",
          fmt(ets_5d_mean["log9010"])),
  "\\addlinespace[6pt]",
  "\\multicolumn{5}{l}{\\textit{Panel B: Literature benchmarks}} \\\\",
  "\\addlinespace[2pt]",
  "Lyubich et al.\\ (2018) & 6-digit NAICS & Plant & 9.7 & 2.27 \\\\",
  "De Lyon \\& Dechezlepr\\^{e}tre (2025) & 3--4 digit NACE & Firm ($\\geq$20 emp.) & 24.0 & 3.18 \\\\",
  "\\addlinespace[6pt]",
  "\\multicolumn{5}{l}{\\textit{Panel C: Observed + imputed emissions (all firms)}} \\\\",
  "\\addlinespace[2pt]",
  sprintf("Belgium & NACE 2-digit & Firm & %s & %s \\\\",
          fmt1(full_2d_mean["p90p10"]), fmt(full_2d_mean["log9010"])),
  sprintf("Belgium & NACE 5-digit & Firm & --- & %s \\\\",
          fmt(full_5d_mean["log9010"])),
  "\\bottomrule",
  "\\end{tabular}"
)

writeLines(tex1, file.path(OUT_DIR, "dispersion_overview.tex"))
cat("-- Wrote dispersion_overview.tex\n")


# =============================================================================
# TABLE 2 (CSV + LaTeX): Sector-by-sector, ETS vs full (NACE 2-digit)
# =============================================================================

nace2d_labels <- c(
  "10" = "Food products", "11" = "Beverages", "12" = "Tobacco",
  "13" = "Textiles", "14" = "Wearing apparel", "15" = "Leather",
  "16" = "Wood products", "17/18" = "Paper \\& printing",
  "19" = "Coke \\& refined petroleum",
  "20" = "Chemicals", "21" = "Pharmaceuticals",
  "22" = "Rubber \\& plastics", "23" = "Non-metallic minerals",
  "24" = "Basic metals", "25" = "Fabricated metals",
  "26" = "Computer \\& electronics", "27" = "Electrical equipment",
  "28" = "Machinery n.e.c.", "29" = "Motor vehicles",
  "30" = "Other transport eq.", "31" = "Furniture",
  "32" = "Other manufacturing", "33" = "Repair \\& installation"
)

# Plain-text labels for CSV
nace2d_labels_csv <- c(
  "10" = "Food", "11" = "Beverages", "12" = "Tobacco",
  "13" = "Textiles", "14" = "Wearing apparel", "15" = "Leather",
  "16" = "Wood", "17/18" = "Paper & printing", "19" = "Coke & petroleum",
  "20" = "Chemicals", "21" = "Pharmaceuticals", "22" = "Rubber & plastics",
  "23" = "Non-metallic minerals", "24" = "Basic metals",
  "25" = "Fabricated metals", "26" = "Electronics",
  "27" = "Electrical equipment", "28" = "Machinery",
  "29" = "Motor vehicles", "30" = "Other transport",
  "31" = "Furniture", "32" = "Other manufacturing",
  "33" = "Repair & installation"
)

# -- CSV table ----------------------------------------------------------------
tab2_csv <- full_2d %>%
  left_join(ets_2d, by = "nace2d") %>%
  mutate(
    label = nace2d_labels_csv[nace2d],
    inflation_p90p10 = full_p90p10 / ets_p90p10,
    inflation_log    = full_log9010 / ets_log9010
  ) %>%
  transmute(
    NACE     = nace2d,
    Sector   = label,
    ETS_n    = ets_n_firms_avg,
    ETS_yrs  = ets_n_years,
    Full_n   = full_n_firms_avg,
    ETS_p90p10  = round(ets_p90p10, 1),
    Full_p90p10 = round(full_p90p10, 1),
    ETS_log9010  = round(ets_log9010, 2),
    Full_log9010 = round(full_log9010, 2),
    Ratio_p90p10 = round(inflation_p90p10, 1)
  ) %>%
  arrange(NACE)

cat("\n")
cat("======================================================================\n")
cat("  TABLE 2: ETS-only vs. Full Sample -- NACE 2-digit\n")
cat("======================================================================\n\n")
print(tab2_csv, row.names = FALSE, right = FALSE, na.print = "---")
write.csv(tab2_csv, file.path(OUT_DIR, "dispersion_ets_vs_full_2d.csv"),
          row.names = FALSE, na = "")

# -- LaTeX table --------------------------------------------------------------
tab_data <- full_2d %>%
  left_join(ets_2d, by = "nace2d") %>%
  left_join(ets_5d, by = "nace2d") %>%
  left_join(full_5d, by = "nace2d") %>%
  mutate(label = nace2d_labels[nace2d]) %>%
  arrange(nace2d)

tex2 <- c(
  "\\begin{tabular}{l r r cc cc}",
  "\\toprule",
  " & & & \\multicolumn{2}{c}{ETS only} & \\multicolumn{2}{c}{All firms} \\\\",
  "\\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
  " & $\\bar{N}_{\\text{ETS}}$ & $\\bar{N}_{\\text{all}}$ & 2-digit & 5-digit & 2-digit & 5-digit \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(tab_data))) {
  r <- tab_data[i, ]
  tex2 <- c(tex2, sprintf(
    "%s & %s & %s & %s & %s & %s & %s \\\\",
    r$label,
    fmt0(r$ets_n_firms_avg),
    fmt0(r$full_n_firms_avg),
    fmt(r$ets_log9010),
    fmt(r$ets5d_log9010),
    fmt(r$full_log9010),
    fmt(r$full5d_log9010)
  ))
}

tex2 <- c(tex2,
  "\\midrule",
  sprintf(
    "Mean & & & %s & %s & %s & %s \\\\",
    fmt(mean(tab_data$ets_log9010, na.rm = TRUE)),
    fmt(mean(tab_data$ets5d_log9010, na.rm = TRUE)),
    fmt(mean(tab_data$full_log9010, na.rm = TRUE)),
    fmt(mean(tab_data$full5d_log9010, na.rm = TRUE))
  ),
  "\\addlinespace[3pt]",
  "\\multicolumn{7}{l}{\\textit{Literature benchmarks}} \\\\",
  "\\addlinespace[2pt]",
  "Lyubich et al.\\ (2018) & & & & 2.27 & & \\\\",
  "De Lyon \\& Dechezlepr\\^{e}tre (2025) & & & 3.18 & & & \\\\",
  "\\bottomrule",
  "\\end{tabular}"
)

writeLines(tex2, file.path(OUT_DIR, "dispersion_by_sector.tex"))
cat("-- Wrote dispersion_by_sector.tex\n")


# =============================================================================
# TABLE 3 (CSV): NACE 5-digit comparison
# =============================================================================

tab3 <- full_5d %>%
  left_join(ets_5d, by = "nace2d") %>%
  mutate(
    label = nace2d_labels_csv[nace2d],
    inflation_log = full5d_log9010 / ets5d_log9010
  ) %>%
  transmute(
    NACE     = nace2d,
    Sector   = label,
    ETS_5d   = ets5d_n_sectors,
    Full_5d  = full5d_n_sectors,
    ETS_log9010  = round(ets5d_log9010, 2),
    Full_log9010 = round(full5d_log9010, 2),
    Ratio_log    = round(inflation_log, 1)
  ) %>%
  arrange(NACE)

cat("\n")
cat("======================================================================\n")
cat("  TABLE 3: ETS-only vs. Full Sample -- NACE 5-digit\n")
cat("======================================================================\n\n")
print(tab3, row.names = FALSE, right = FALSE, na.print = "---")
write.csv(tab3, file.path(OUT_DIR, "dispersion_ets_vs_full_5d.csv"),
          row.names = FALSE, na = "")

cat("\n-- Done. Tables saved to", OUT_DIR, "--\n")
