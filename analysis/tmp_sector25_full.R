library(dplyr, warn.conflicts = FALSE)

DATA_DIR  <- "C:/Users/jota_/Documents/NBB_data/processed"
REPO_DIR  <- "C:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

# Load accounts for NACE codes + revenue
load(file.path(DATA_DIR, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year >= 2005, year <= 2021) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

# Load deployment panel to know who is deployment
load(file.path(DATA_DIR, "deployment_panel.RData"))
deploy_vats <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  filter(nace2d == "25") %>%
  distinct(vat, year)
rm(deployment_panel)

cat("Deployment firm-years in NACE 25:", nrow(deploy_vats), "\n")
cat("Unique deployment firms:", n_distinct(deploy_vats$vat), "\n\n")

# CRF group for sector 25
nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE, colClasses = c(nace2d = "character")
)
crf25 <- nace_crf$crf_group[nace_crf$nace2d == "25"]
cat("NACE 25 maps to CRF group:", crf25, "\n")
cat("Other NACE sectors in same CRF group:",
    paste(nace_crf$nace2d[nace_crf$crf_group == crf25 & nace_crf$nace2d != "25"],
          collapse = ", "), "\n\n")

# Load one year file to see imputed emissions for sector 25
# We need to reconstruct from the draw-level data
# Let's use year 2019 as example and look at draw-level stats
load(file.path(DATA_DIR, "b_loop_pareto_year_2019.RData"))

cat("══════════════════════════════════════════════════════════════\n")
cat("  SECTOR 25 — FULL SAMPLE (2019, across draws)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

s25_2019 <- year_stats2d %>% filter(nace2d == "25")
cat("Draws with stats for NACE 25:", nrow(s25_2019), "\n")
cat("n_firms per draw: mean =", round(mean(s25_2019$n_firms)),
    ", min =", min(s25_2019$n_firms), ", max =", max(s25_2019$n_firms), "\n")
cat("cp_p90p10: mean =", round(mean(s25_2019$cp_p90p10, na.rm=TRUE), 1),
    ", sd =", round(sd(s25_2019$cp_p90p10, na.rm=TRUE), 1),
    ", min =", round(min(s25_2019$cp_p90p10, na.rm=TRUE), 1),
    ", max =", round(max(s25_2019$cp_p90p10, na.rm=TRUE), 1), "\n")
cat("cp_p9010_log: mean =", round(mean(s25_2019$cp_p9010_log, na.rm=TRUE), 2),
    ", sd =", round(sd(s25_2019$cp_p9010_log, na.rm=TRUE), 2), "\n\n")

# Now reconstruct the actual imputed emissions for one draw
# We need the proxy list + the allocation logic
# Instead, let's look at what we can infer from the stats

# Load proxy list — just draw 1 for sector 25
load(file.path(DATA_DIR, "deployment_proxy_list.RData"))

# Draw 1, year 2019
px1 <- proxy_list[[1]]
px1_2019 <- px1[px1$year == 2019, ]

# Match to deployment firms in NACE 25
deploy25_2019 <- deploy_vats %>% filter(year == 2019)
px1_25 <- px1_2019 %>% filter(vat %in% deploy25_2019$vat)

cat("── Proxy distribution for NACE 25 deployment firms (draw 1, 2019) ──\n")
cat("Total deployment firms:", nrow(deploy25_2019), "\n")
cat("With proxy value:", nrow(px1_25), "\n")
cat("With proxy > 0:", sum(px1_25$proxy > 0), "\n")
cat("With proxy = 0:", sum(px1_25$proxy == 0), "\n")
cat("With proxy < 0:", sum(px1_25$proxy < 0), "\n\n")

pos_proxy <- px1_25$proxy[px1_25$proxy > 0]
cat("── Among positive-proxy firms ──\n")
cat("N:", length(pos_proxy), "\n")

# Threshold at 62nd percentile
thr <- quantile(sinh(pos_proxy), probs = 0.62, names = FALSE)
sinh_vals <- sinh(pos_proxy)
n_above <- sum(sinh_vals > thr)
n_below <- sum(sinh_vals <= thr)
cat("sinh(proxy) quantiles:\n")
print(round(quantile(sinh_vals, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)), 4))
cat("\nAfter 62% threshold:\n")
cat("  Firms surviving (emit > 0):", n_above, "\n")
cat("  Firms zeroed out:", n_below, "\n")
cat("  Ratio surviving:", round(n_above / length(pos_proxy) * 100, 1), "%\n")

# Ratio of max to min surviving proxy
surviving <- sinh_vals[sinh_vals > thr]
cat("\n── Surviving firms' sinh(proxy) range ──\n")
cat("  min:", round(min(surviving), 4), "\n")
cat("  max:", round(max(surviving), 4), "\n")
cat("  max/min ratio:", round(max(surviving) / min(surviving), 1), "\n")
print(round(quantile(surviving, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), 4))

rm(proxy_list)
gc(verbose = FALSE)

# NIR calibration: how much E_deploy goes to sector 25's CRF group?
load(file.path(DATA_DIR, "nir_calibration_targets.RData"))
nir25 <- nir_targets %>% filter(crf_group == crf25)
cat("\n── NIR targets for CRF group '", crf25, "' ──\n", sep = "")
print(as.data.frame(nir25 %>% filter(year >= 2017, year <= 2021) %>%
        mutate(E_NIR_kt = round(E_NIR_kt, 1))), row.names = FALSE)
