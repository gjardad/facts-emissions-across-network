###############################################################################
# analysis/diagnose_rowsum_ge1.R
#
# PURPOSE
#   Identify WHY some rows of the A matrix have row sums >= 1.
#   For each buyer j, row_sum_j = domestic_inputs_j / cost_vec_j.
#   cost_vec_j = domestic_inputs_j + wage_bill_j + imports_j + carbon_j + 1e-6
#   So row_sum >= 1 only when the non-domestic-input components are ~zero.
#
#   This script profiles those firms: how many are there, what do they look
#   like, and where does the gap come from?
#
# RUNS ON: RMD or local 1
###############################################################################

# ── Paths ────────────────────────────────────────────────────────────────────
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

YEARS <- 2005:2021

CARBON_PRICE <- c(
  "2005" = 25.29, "2006" = 21.53, "2007" =  0.86, "2008" = 25.74,
  "2009" = 18.41, "2010" = 18.98, "2011" = 18.08, "2012" =  9.49,
  "2013" =  5.94, "2014" =  7.89, "2015" =  8.52, "2016" =  5.92,
  "2017" =  6.63, "2018" = 18.55, "2019" = 27.84, "2020" = 27.94,
  "2021" = 62.25
)

# ── Load data ────────────────────────────────────────────────────────────────
cat("── Loading data ─────────────────────────────────────────────\n")

load(file.path(PROC_DATA, "b2b_selected_sample.RData"))
b2b <- df_b2b_selected_sample %>% filter(year %in% YEARS)
rm(df_b2b_selected_sample)

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, wage_bill) %>%
  mutate(wage_bill = pmax(coalesce(wage_bill, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

load(file.path(PROC_DATA, "firm_year_total_imports.RData"))
if ("vat_ano" %in% names(firm_year_total_imports))
  firm_year_total_imports <- rename(firm_year_total_imports, vat = vat_ano)
imports <- firm_year_total_imports %>%
  filter(year %in% YEARS) %>%
  select(vat, year, total_imports)
rm(firm_year_total_imports)

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

cat("  Data loaded.\n\n")

# ── Year-by-year diagnostic ─────────────────────────────────────────────────
summary_rows <- list()

for (t in YEARS) {
  b2b_t      <- b2b %>% filter(year == t)
  accounts_t <- accounts %>% filter(year == t)
  imports_t  <- imports %>% filter(year == t)
  eutl_t     <- eutl %>% filter(year == t)

  # Aggregate B2B: total purchases by each buyer j
  buyer_inputs <- b2b_t %>%
    group_by(vat_j_ano) %>%
    summarise(domestic_inputs = sum(corr_sales_ij, na.rm = TRUE), .groups = "drop") %>%
    filter(domestic_inputs > 0) %>%
    rename(vat = vat_j_ano)

  n_buyers <- nrow(buyer_inputs)

  # Merge with other cost components
  firm <- buyer_inputs %>%
    left_join(accounts_t %>% select(vat, wage_bill), by = "vat") %>%
    left_join(imports_t  %>% select(vat, total_imports), by = "vat") %>%
    left_join(eutl_t     %>% select(vat, emissions), by = "vat") %>%
    mutate(
      wage_bill     = pmax(coalesce(wage_bill, 0), 0),
      total_imports = pmax(coalesce(total_imports, 0), 0),
      carbon_cost   = pmax(coalesce(emissions, 0), 0) * CARBON_PRICE[as.character(t)],
      non_dom       = wage_bill + total_imports + carbon_cost,
      cost_vec      = domestic_inputs + non_dom + 1e-6,
      rowsum        = domestic_inputs / cost_vec
    )

  # --- Profile firms at two thresholds ---
  classify_high <- function(firm_df, threshold) {
    high <- firm_df %>% filter(rowsum >= threshold)
    if (nrow(high) == 0) {
      return(data.frame(
        n_high = 0L, pct_buyers = 0,
        not_in_accounts = 0L, in_acc_zero_wage = 0L,
        no_imports = 0L, no_wage_imp_carb = 0L,
        dom_inputs_high_M = 0, pct_dom_vol = 0
      ))
    }
    high <- high %>%
      mutate(
        has_wage    = wage_bill > 0,
        has_imports = total_imports > 0,
        has_carbon  = carbon_cost > 0,
        in_accounts = vat %in% accounts_t$vat,
        in_imports  = vat %in% imports_t$vat
      )
    data.frame(
      n_high            = nrow(high),
      pct_buyers        = round(100 * nrow(high) / n_buyers, 2),
      not_in_accounts   = sum(!high$in_accounts),
      in_acc_zero_wage  = sum(high$in_accounts & !high$has_wage),
      no_imports        = sum(!high$has_imports),
      no_wage_imp_carb  = sum(!high$has_wage & !high$has_imports & !high$has_carbon),
      dom_inputs_high_M = round(sum(high$domestic_inputs) / 1e6, 1),
      pct_dom_vol       = round(100 * sum(high$domestic_inputs) /
                                  sum(firm$domestic_inputs), 2)
    )
  }

  stats_999 <- classify_high(firm, 0.999)
  stats_95  <- classify_high(firm, 0.95)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    year            = t,
    n_buyers        = n_buyers,
    dom_inputs_all_M = round(sum(firm$domestic_inputs) / 1e6, 1),
    # --- threshold 0.95 ---
    n_rs_ge95           = stats_95$n_high,
    pct_buyers_95       = stats_95$pct_buyers,
    not_in_acc_95       = stats_95$not_in_accounts,
    in_acc_zero_wage_95 = stats_95$in_acc_zero_wage,
    no_imports_95       = stats_95$no_imports,
    no_wage_imp_carb_95 = stats_95$no_wage_imp_carb,
    pct_dom_vol_95      = stats_95$pct_dom_vol,
    # --- threshold 0.999 ---
    n_rs_ge999           = stats_999$n_high,
    pct_buyers_999       = stats_999$pct_buyers,
    not_in_acc_999       = stats_999$not_in_accounts,
    in_acc_zero_wage_999 = stats_999$in_acc_zero_wage,
    no_imports_999       = stats_999$no_imports,
    no_wage_imp_carb_999 = stats_999$no_wage_imp_carb,
    pct_dom_vol_999      = stats_999$pct_dom_vol
  )

  if (t == YEARS[1]) {
    # Print a few example firms in year 1 for inspection
    example <- firm %>%
      filter(rowsum >= 0.95) %>%
      mutate(
        in_accounts = vat %in% accounts_t$vat,
        in_imports  = vat %in% imports_t$vat
      ) %>%
      arrange(desc(domestic_inputs)) %>%
      head(15) %>%
      select(vat, domestic_inputs, wage_bill, total_imports, carbon_cost,
             rowsum, in_accounts, in_imports) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
    cat(sprintf("── Year %d: example firms with row sum >= 0.95 ──\n", t))
    print(example)
    cat("\n")
  }
}

# ── Summary table ────────────────────────────────────────────────────────────
summary_df <- bind_rows(summary_rows)

cat("── Summary by year ────────────────────────────────────────────\n")
cat("  Suffix _95 = row sum >= 0.95 | _999 = row sum >= 0.999\n")
cat("  not_in_acc     = # not found in annual accounts at all\n")
cat("  in_acc_zero_wage = found in accounts but wage_bill = 0\n")
cat("  no_imports     = # with zero imports\n")
cat("  no_wage_imp_carb = # with ALL three non-domestic components zero\n")
cat("  pct_dom_vol    = % of total domestic input volume from these firms\n\n")

cat("--- Threshold: row sum >= 0.95 ---\n")
print(
  as.data.frame(summary_df %>%
    select(year, n_buyers, dom_inputs_all_M,
           n_rs_ge95, pct_buyers_95, not_in_acc_95,
           in_acc_zero_wage_95, no_imports_95,
           no_wage_imp_carb_95, pct_dom_vol_95)),
  row.names = FALSE
)

cat("\n--- Threshold: row sum >= 0.999 ---\n")
print(
  as.data.frame(summary_df %>%
    select(year, n_buyers,
           n_rs_ge999, pct_buyers_999, not_in_acc_999,
           in_acc_zero_wage_999, no_imports_999,
           no_wage_imp_carb_999, pct_dom_vol_999)),
  row.names = FALSE
)
