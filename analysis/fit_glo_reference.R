###############################################################################
# analysis/fit_glo_reference.R
#
# PURPOSE
#   Fit a single GLO (Generalized Logistic) distribution to the within-cell
#   shape of log emissions, pooled across ALL sectors and years.
#
#   Reference distribution construction (per inferring_emissions/analysis/
#   active/threshold_p_levels_eval_all.R:84, with the mixed-sector exclusion
#   removed -- for deployment we want all sectors):
#
#     1. Take training firm-years with y > 0 (true emitters), all CRF groups.
#     2. log_y = log(y)
#     3. tilde = log_y - mean(log_y) within year
#     4. d     = tilde  - mean(tilde)  within primary_crf_group
#     5. Pool {d_i}, fit GLO via L-moments (samlmu + pelglo).
#
#   Under a location-shift assumption (within-cell shapes differ only in
#   location), the L-moments lambda2, tau3, tau4 of pooled {d_i} consistently
#   estimate the L-moments of the common within-cell shape. The fitted GLO
#   parameters are reused unchanged across all (CRF group x year) cells in
#   the deployment allocation.
#
# INPUT
#   {PROC_DATA}/repeated_cv_proxy_crf_asinh.RData
#     repeated_cv_proxy_panel : data.frame with vat, year, y,
#                               primary_crf_group, ...
#
# OUTPUT
#   {PROC_DATA}/glo_reference_params.RData
#     glo_par : named numeric vector (xi, alpha, k)
#     lmoms   : sample L-moments (length 4)
#     n_ref   : number of firm-years used in the fit
#
# RUNS ON: local 1 or RMD (training-sample-derived file is on both).
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
source(file.path(REPO_DIR, "utils", "glo_fit.R"))

suppressPackageStartupMessages(library(dplyr))

cat("===================================================================\n")
cat("  FIT GLO REFERENCE DISTRIBUTION (all sectors)\n")
cat("===================================================================\n\n")

# -- Load training panel ------------------------------------------------------
load(file.path(PROC_DATA, "repeated_cv_proxy_crf_asinh.RData"))

panel <- repeated_cv_proxy_panel
cat("Loaded repeated_cv_proxy_panel:", nrow(panel), "firm-years\n")
cat("  Columns:", paste(names(panel), collapse = ", "), "\n\n")

# -- Build reference distribution --------------------------------------------
# All sectors, true emitters only.
ref_df <- panel %>%
  filter(!is.na(y), y > 0, !is.na(primary_crf_group)) %>%
  mutate(log_y = log(y)) %>%
  group_by(year)              %>% mutate(tilde = log_y - mean(log_y)) %>% ungroup() %>%
  group_by(primary_crf_group) %>% mutate(d     = tilde - mean(tilde)) %>% ungroup()

cat("Reference distribution:\n")
cat("  n firm-years (y > 0):", nrow(ref_df), "\n")
cat("  Distinct CRF groups :", n_distinct(ref_df$primary_crf_group), "\n")
cat("  Year range          :", min(ref_df$year), "-", max(ref_df$year), "\n\n")

# -- Fit GLO ------------------------------------------------------------------
glo_par <- fit_glo(ref_df$d)
if (is.null(glo_par)) stop("GLO fit failed.")

lmoms <- attr(glo_par, "lmoms")
n_ref <- attr(glo_par, "n")

cat("L-moments (lambda1, lambda2, tau3, tau4):\n"); print(round(lmoms, 4))
cat("\nGLO parameters (xi, alpha, k):\n"); print(round(glo_par, 4))
cat("\nn used for fit:", n_ref, "\n")

# Strip attributes before saving for clean object
glo_par <- as.numeric(glo_par)
names(glo_par) <- c("xi", "alpha", "k")

# -- Save ---------------------------------------------------------------------
out_path <- file.path(PROC_DATA, "glo_reference_params.RData")
save(glo_par, lmoms, n_ref, file = out_path)
cat("\nSaved:", out_path, "\n")
