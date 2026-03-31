###############################################################################
# analysis/table_dispersion_comparison_latex.R
#
# PURPOSE
#   Produce LaTeX tables comparing carbon productivity dispersion:
#     Table 1: Belgium (ETS-only + all firms) vs. literature benchmarks
#     Table 2: Sector-by-sector ETS-only vs. full sample (NACE 2-digit)
#
# OUTPUT
#   {REPO_DIR}/output/tables/dispersion_overview.tex
#   {REPO_DIR}/output/tables/dispersion_by_sector.tex
###############################################################################

library(dplyr, warn.conflicts = FALSE)

DATA_DIR  <- "C:/Users/jota_/Documents/NBB_data/processed"
REPO_DIR  <- "C:/Users/jota_/Documents/facts-emissions-across-network"
OUT_DIR   <- file.path(REPO_DIR, "output", "tables")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

MIN_N <- 3L
YEARS <- 2005:2021

# ── Helpers ─────────────────────────────────────────────────────────────────
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

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "---", formatC(round(x, digits), format = "f", digits = digits))
}

fmt1 <- function(x) fmt(x, 1)
fmt0 <- function(x) ifelse(is.na(x), "---", as.character(round(x, 0)))

# ═══════════════════════════════════════════════════════════════════════════
# DATA
# ═══════════════════════════════════════════════════════════════════════════
cat("── Loading data ──\n")

load(file.path(DATA_DIR, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS, !is.na(emissions), emissions > 0) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

load(file.path(DATA_DIR, "annual_accounts_selected_sample_key_variables.RData"))
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

# ── ETS 2-digit ─────────────────────────────────────────────────────────────
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

# ── ETS 5-digit (aggregated to 2-digit parent) ─────────────────────────────
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

# ── Full sample from Pareto loop ───────────────────────────────────────────
cat("── Loading Pareto loop results ──\n")

all_2d <- list()
all_5d <- list()
for (yr in YEARS) {
  f <- file.path(DATA_DIR, sprintf("b_loop_pareto_year_%d.RData", yr))
  if (file.exists(f)) {
    load(f)
    all_2d[[as.character(yr)]] <- year_stats2d
    all_5d[[as.character(yr)]] <- year_stats5d
  }
}

full_2d <- bind_rows(all_2d) %>%
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

full_5d <- bind_rows(all_5d) %>%
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


# ═══════════════════════════════════════════════════════════════════════════
# TABLE 1: Overview — Belgium vs. Literature
# ═══════════════════════════════════════════════════════════════════════════

ets_2d_mean  <- c(mean(ets_2d$ets_p90p10, na.rm = TRUE),
                  mean(ets_2d$ets_log9010, na.rm = TRUE))
ets_5d_mean  <- c(mean(ets_5d$ets5d_log9010, na.rm = TRUE))  # p90p10 less meaningful
full_2d_mean <- c(mean(full_2d$full_p90p10, na.rm = TRUE),
                  mean(full_2d$full_log9010, na.rm = TRUE))
full_5d_mean <- c(mean(full_5d$full5d_log9010, na.rm = TRUE))

tex1 <- c(
  "\\begin{tabular}{l l l cc}",
  "\\toprule",
  " & Granularity & Unit & $p_{90}/p_{10}$ & $\\ln p_{90} - \\ln p_{10}$ \\\\",
  "\\midrule",
  "\\addlinespace[2pt]",
  "\\multicolumn{5}{l}{\\textit{Panel A: Observed emissions (ETS firms only)}} \\\\",
  "\\addlinespace[2pt]",
  sprintf("Belgium & NACE 2-digit & Firm & %s & %s \\\\",
          fmt1(ets_2d_mean[1]), fmt(ets_2d_mean[2])),
  sprintf("Belgium & NACE 5-digit & Firm & --- & %s \\\\",
          fmt(ets_5d_mean[1])),
  "\\addlinespace[6pt]",
  "\\multicolumn{5}{l}{\\textit{Panel B: Literature benchmarks}} \\\\",
  "\\addlinespace[2pt]",
  "Lyubich et al.\\ (2018) & 6-digit NAICS & Plant & 9.7 & 2.27 \\\\",
  "De Lyon \\& Dechezlepr\\^{e}tre (2025) & 3--4 digit NACE & Firm ($\\geq$20 emp.) & 24.0 & 3.18 \\\\",
  "\\addlinespace[6pt]",
  "\\multicolumn{5}{l}{\\textit{Panel C: Observed + imputed emissions (all firms)}} \\\\",
  "\\addlinespace[2pt]",
  sprintf("Belgium & NACE 2-digit & Firm & %s & %s \\\\",
          fmt1(full_2d_mean[1]), fmt(full_2d_mean[2])),
  sprintf("Belgium & NACE 5-digit & Firm & --- & %s \\\\",
          fmt(full_5d_mean[1])),
  "\\bottomrule",
  "\\end{tabular}"
)

writeLines(tex1, file.path(OUT_DIR, "dispersion_overview.tex"))
cat("── Wrote dispersion_overview.tex ──\n")


# ═══════════════════════════════════════════════════════════════════════════
# TABLE 2: Sector-by-sector (NACE 2-digit)
# ═══════════════════════════════════════════════════════════════════════════

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

# Overall mean row
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
cat("── Wrote dispersion_by_sector.tex ──\n")
cat("── Done ──\n")
