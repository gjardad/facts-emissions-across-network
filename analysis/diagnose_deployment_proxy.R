###############################################################################
# analysis/diagnose_deployment_proxy.R
#
# PURPOSE
#   Diagnose the distribution of proxy_avg among deployment (non-ETS) firms.
#   Key question: does the 34.5% predicted-emitter rate reflect many
#   near-zero proxy firms (benign) or genuinely widespread positive proxies
#   (concerning)? And how concentrated is the proxy signal within sectors —
#   i.e., would the NIR residual effectively concentrate on a few high-proxy
#   firms regardless of the classification threshold?
#
# INPUT
#   {PROC_DATA}/deployment_proxy_avg.RData
#   {PROC_DATA}/annual_accounts_selected_sample_key_variables.RData
#
# OUTPUT
#   Console diagnostics + tables/diagnose_proxy_*.csv
#
# RUNS ON: local 1 (or RMD)
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
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)

dir.create(file.path(REPO_DIR, "tables"), showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
load(file.path(PROC_DATA, "deployment_proxy_avg.RData"))
load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))

# Attach NACE 2-digit (modal sector across years per firm)
nace <- df_annual_accounts_selected_sample_key_variables %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  group_by(vat) %>%
  summarise(nace2d = names(which.max(table(nace2d))), .groups = "drop")

proxy <- proxy_avg %>%
  left_join(nace, by = "vat") %>%
  mutate(
    emitter    = proxy_avg > 0,
    sinh_proxy = sinh(proxy_avg)
  )

cat("══════════════════════════════════════════════════════════\n")
cat("  DEPLOYMENT PROXY DIAGNOSTICS\n")
cat("══════════════════════════════════════════════════════════\n\n")


# =============================================================================
# 1. Overall summary
# =============================================================================
cat("── 1. Overall ──────────────────────────────────────────────\n")
cat("Total firm-years:         ", nrow(proxy), "\n")
cat("Predicted emitters:       ", sum(proxy$emitter),
    sprintf("(%.1f%%)\n", 100 * mean(proxy$emitter)))
cat("Unique firms:             ", n_distinct(proxy$vat), "\n")
cat("Unique predicted emitters:", n_distinct(proxy$vat[proxy$emitter]), "\n\n")

cat("Distribution of proxy_avg among predicted emitters:\n")
q <- quantile(proxy$proxy_avg[proxy$emitter],
              probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99))
print(round(q, 4))
cat("\n")


# =============================================================================
# 2. By NACE 2-digit: predicted emitter rate + proxy quantiles
# =============================================================================
cat("── 2. Predicted emitter rate by NACE 2-digit ───────────────\n")

by_sector <- proxy %>%
  group_by(nace2d) %>%
  summarise(
    n_firmyears   = n(),
    n_emitters    = sum(emitter),
    pct_emitters  = 100 * mean(emitter),
    p25_proxy_pos = quantile(proxy_avg[emitter], 0.25, na.rm = TRUE),
    p50_proxy_pos = quantile(proxy_avg[emitter], 0.50, na.rm = TRUE),
    p75_proxy_pos = quantile(proxy_avg[emitter], 0.75, na.rm = TRUE),
    p90_proxy_pos = quantile(proxy_avg[emitter], 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_emitters))

print(by_sector, n = 50)
write.csv(by_sector,
          file.path(REPO_DIR, "tables", "diagnose_proxy_by_sector.csv"),
          row.names = FALSE)
cat("\nSaved: tables/diagnose_proxy_by_sector.csv\n\n")


# =============================================================================
# 3. NIR concentration: what share of sinh(proxy) goes to top X% of firms
#    within each sector-year? This tells us how concentrated imputed emissions
#    would be regardless of how many near-zero firms are classified as emitters.
# =============================================================================
cat("── 3. NIR concentration: share of sinh(proxy) by top-X firms ──\n")

concentration <- proxy %>%
  filter(emitter) %>%
  group_by(nace2d, year) %>%
  arrange(desc(proxy_avg), .by_group = TRUE) %>%
  mutate(
    rank         = row_number(),
    n            = n(),
    sinh_total   = sum(sinh_proxy),
    sinh_cumshare = cumsum(sinh_proxy) / sinh_total,
    pctile_rank  = rank / n
  ) %>%
  summarise(
    n_emitters          = first(n),
    top10pct_share      = sinh_cumshare[which.min(abs(pctile_rank - 0.10))],
    top25pct_share      = sinh_cumshare[which.min(abs(pctile_rank - 0.25))],
    top50pct_share      = sinh_cumshare[which.min(abs(pctile_rank - 0.50))],
    .groups = "drop"
  )

# Summarise across sector-years (median concentration by sector)
conc_by_sector <- concentration %>%
  group_by(nace2d) %>%
  summarise(
    n_sy              = n(),
    med_n_emitters    = median(n_emitters),
    med_top10pct_share = median(top10pct_share, na.rm = TRUE),
    med_top25pct_share = median(top25pct_share, na.rm = TRUE),
    med_top50pct_share = median(top50pct_share, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(med_top10pct_share))

cat("Median share of sinh(proxy) held by top 10%/25%/50% of predicted emitters,\n")
cat("by NACE 2-digit (higher = more concentrated = near-zero firms matter less):\n\n")
print(conc_by_sector, n = 50)
write.csv(conc_by_sector,
          file.path(REPO_DIR, "tables", "diagnose_proxy_concentration.csv"),
          row.names = FALSE)
cat("\nSaved: tables/diagnose_proxy_concentration.csv\n\n")

# Overall concentration (pooled across all sector-years)
cat("── Overall concentration (pooled) ──────────────────────────\n")
overall_conc <- proxy %>%
  filter(emitter) %>%
  group_by(nace2d, year) %>%
  arrange(desc(proxy_avg), .by_group = TRUE) %>%
  mutate(
    n             = n(),
    sinh_total    = sum(sinh_proxy),
    sinh_cumshare = cumsum(sinh_proxy) / sinh_total,
    pctile_rank   = row_number() / n
  ) %>%
  ungroup() %>%
  summarise(
    med_top10pct_share = median(
      sapply(split(., interaction(.$nace2d, .$year)), function(d)
        d$sinh_cumshare[which.min(abs(d$pctile_rank - 0.10))]),
      na.rm = TRUE),
    med_top25pct_share = median(
      sapply(split(., interaction(.$nace2d, .$year)), function(d)
        d$sinh_cumshare[which.min(abs(d$pctile_rank - 0.25))]),
      na.rm = TRUE),
    med_top50pct_share = median(
      sapply(split(., interaction(.$nace2d, .$year)), function(d)
        d$sinh_cumshare[which.min(abs(d$pctile_rank - 0.50))]),
      na.rm = TRUE)
  )
cat("Median (across sector-years) share of sinh(proxy) held by:\n")
cat("  Top 10% of emitters:", round(100 * overall_conc$med_top10pct_share, 1), "%\n")
cat("  Top 25% of emitters:", round(100 * overall_conc$med_top25pct_share, 1), "%\n")
cat("  Top 50% of emitters:", round(100 * overall_conc$med_top50pct_share, 1), "%\n\n")


# =============================================================================
# 4. Near-zero proxy: how many predicted emitters have negligible proxy values?
#    Use training sample proxy distribution as a reference for "negligible".
# =============================================================================
cat("── 4. Near-zero proxy among predicted emitters ─────────────\n")

thresholds <- c(0.001, 0.01, 0.1, 0.5, 1, 5)
emitters   <- proxy %>% filter(emitter)
cat("Among", nrow(emitters), "predicted emitter firm-years:\n")
for (thr in thresholds) {
  n_below <- sum(emitters$proxy_avg < thr)
  cat(sprintf("  proxy_avg < %-6g : %6d  (%.1f%%)\n",
              thr, n_below, 100 * n_below / nrow(emitters)))
}
cat("\n")

cat("══════════════════════════════════════════════════════════\n")
cat("Outputs saved to tables/\n")
cat("══════════════════════════════════════════════════════════\n")
