###############################################################################
# analysis/build_summary_objects.R
#
# PURPOSE
#   Collapse the 200 CV repeats in repeated_cv_proxy_crf_asinh.RData into
#   per firm-year summary objects:
#     p_i          = mean(proxy_ir > 0) across 200 repeats
#     proxy_mean_i = mean(proxy_ir)     across 200 repeats
#
#   These are the inputs to the extensive- and intensive-margin calibration
#   models (see deployment_prediction_interval_step_by_step.md, Step 1).
#
# INPUT
#   {PROC_DATA}/repeated_cv_proxy_crf_asinh.RData
#     proxy_matrix            : N x 200 matrix of held-out proxy values
#     repeated_cv_proxy_panel : firm-year panel with vat, year, nace2d, y,
#                                emit, primary_crf_group, ...
#
# OUTPUT
#   {PROC_DATA}/training_summary_objects.RData
#     training_summary : data.frame with columns
#       vat, year, nace2d, primary_crf_group, D_i, y_i,
#       p_i, proxy_mean_i
#
# RUNS ON: local 1
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

cat("===================================================================\n")
cat("  BUILD TRAINING SUMMARY OBJECTS\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load
# =============================================================================

load(file.path(PROC_DATA, "repeated_cv_proxy_crf_asinh.RData"))

cat("Loaded repeated_cv_proxy_crf_asinh.RData\n")
cat("  proxy_matrix:", nrow(proxy_matrix), "firm-years x", ncol(proxy_matrix), "repeats\n")
cat("  panel cols:  ", paste(names(repeated_cv_proxy_panel), collapse = ", "), "\n\n")

stopifnot(nrow(proxy_matrix) == nrow(repeated_cv_proxy_panel))


# =============================================================================
# SECTION 2: Compute summary objects
# =============================================================================

cat("Computing p_i and proxy_mean_i across 200 repeats...\n")

p_i          <- rowMeans(proxy_matrix > 0)
proxy_mean_i <- rowMeans(proxy_matrix)

cat("  p_i summary:\n")
print(summary(p_i))
cat("\n  proxy_mean_i summary:\n")
print(summary(proxy_mean_i))


# =============================================================================
# SECTION 3: Assemble training panel
# =============================================================================

training_summary <- data.frame(
  vat               = repeated_cv_proxy_panel$vat,
  year              = repeated_cv_proxy_panel$year,
  nace2d            = repeated_cv_proxy_panel$nace2d,
  primary_crf_group = repeated_cv_proxy_panel$primary_crf_group,
  D_i               = as.integer(repeated_cv_proxy_panel$emit),
  y_i               = repeated_cv_proxy_panel$y,
  p_i               = p_i,
  proxy_mean_i      = proxy_mean_i,
  stringsAsFactors  = FALSE
)

cat("\nTraining summary panel:\n")
cat("  Rows:", nrow(training_summary), "\n")
cat("  True emitters (D_i = 1):", sum(training_summary$D_i == 1), "\n")
cat("  CRF groups:\n")
print(table(training_summary$primary_crf_group))


# =============================================================================
# SECTION 4: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA, "training_summary_objects.RData")
save(training_summary, file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat("  ", nrow(training_summary), "firm-years\n")
cat("===================================================================\n")
