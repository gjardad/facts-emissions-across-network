###############################################################################
# analysis/calibrate_extensive_margin.R
#
# PURPOSE
#   Calibrate p(x) = P(true emitter | proxy percentile = x) via logistic
#   regression on held-out firms in mixed CRF categories.
#
#   For each of M = 200 CV repeats, fit logistic regression of true emitter
#   status on proxy-derived explanatory variables among held-out firms in
#   mixed CRFs (paper, refinery, metals). Average coefficients across
#   repeats to obtain a single calibration function p_hat(x).
#
#   Two specifications:
#     Model A: raw proxy + CRF fixed effects
#     Model B: within-CRF percentile of year-demeaned proxy (pooled across years)
#
# INPUT
#   {PROC_DATA}/repeated_cv_proxy_crf_asinh.RData
#     proxy_matrix           : N x 200 matrix of proxy values
#     repeated_cv_proxy_panel: data.frame with vat, year, emit, primary_crf_group
#     crf_group_map          : NACE 2d -> CRF group mapping
#
# OUTPUT
#   {PROC_DATA}/extensive_margin_calibration.RData
#     beta_avg_A, beta_avg_B  : averaged coefficients for Model A and B
#     beta_all_A, beta_all_B  : M x p matrices of per-repeat coefficients
#     p_hat_A(), p_hat_B()    : calibrated probability functions
#     validation              : comparison with observed FPR/FNR
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
cat("  EXTENSIVE MARGIN CALIBRATION\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load data
# =============================================================================

load(file.path(PROC_DATA, "repeated_cv_proxy_crf_asinh.RData"))

cat("Loaded repeated_cv_proxy_crf_asinh.RData\n")
cat("  proxy_matrix:", nrow(proxy_matrix), "firm-years x", ncol(proxy_matrix), "repeats\n")
cat("  Panel columns:", paste(names(repeated_cv_proxy_panel), collapse = ", "), "\n\n")

# Mixed CRF categories: only these have both emitters and non-emitters in training
MIXED_CRFS <- c("paper", "refining", "metals")

panel <- repeated_cv_proxy_panel
M     <- ncol(proxy_matrix)

# Restrict to mixed CRFs
mixed_idx <- which(panel$primary_crf_group %in% MIXED_CRFS)
cat("Firms in mixed CRFs:", length(mixed_idx), "out of", nrow(panel), "\n")
cat("  CRF breakdown:\n")
for (crf in MIXED_CRFS) {
  n_crf   <- sum(panel$primary_crf_group[mixed_idx] == crf)
  n_emit  <- sum(panel$emit[mixed_idx][panel$primary_crf_group[mixed_idx] == crf])
  cat(sprintf("    %-10s: %5d firm-years, %4d emitters (%.1f%%)\n",
              crf, n_crf, n_emit, 100 * n_emit / n_crf))
}
cat("\n")


# =============================================================================
# SECTION 2: Fit logistic regressions per repeat
# =============================================================================

cat("Fitting logistic regressions for M =", M, "repeats...\n\n")

# Storage for coefficients
# Model A: intercept + CRF FEs + beta_proxy
# Model B: intercept + beta_percentile
coef_names_A <- NULL
beta_all_A   <- vector("list", M)
beta_all_B   <- matrix(NA_real_, nrow = M, ncol = 2,
                         dimnames = list(NULL, c("intercept", "percentile")))

t0 <- Sys.time()

for (m in seq_len(M)) {

  # --- Extract proxy for this repeat, restricted to mixed CRFs ---------------
  proxy_m <- proxy_matrix[mixed_idx, m]
  emit_m  <- panel$emit[mixed_idx]
  crf_m   <- panel$primary_crf_group[mixed_idx]
  year_m  <- panel$year[mixed_idx]

  df_m <- data.frame(
    emit  = emit_m,
    proxy = proxy_m,
    crf   = factor(crf_m),
    year  = year_m
  )

  # --- Model A: raw proxy + CRF FE ------------------------------------------
  fit_A <- glm(emit ~ proxy + crf, data = df_m, family = binomial(link = "logit"))

  if (m == 1L) {
    coef_names_A <- names(coef(fit_A))
    beta_all_A   <- matrix(NA_real_, nrow = M, ncol = length(coef_names_A),
                            dimnames = list(NULL, coef_names_A))
  }
  beta_all_A[m, ] <- coef(fit_A)

  # --- Model B: within-CRF percentile of year-demeaned proxy -----------------
  # Step (a): demean by CRF-year
  df_m$proxy_demeaned <- ave(df_m$proxy, df_m$crf, df_m$year, FUN = function(x) x - mean(x))

  # Step (b): within-CRF percentile, pooling across years
  df_m$percentile <- ave(df_m$proxy_demeaned, df_m$crf,
                          FUN = function(x) rank(x, ties.method = "average") / length(x))

  # Step (c): fit logistic
  fit_B <- glm(emit ~ percentile, data = df_m, family = binomial(link = "logit"))
  beta_all_B[m, ] <- coef(fit_B)

  if (m %% 50 == 0) {
    elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
    cat(sprintf("  Repeat %3d / %d (%.1fs elapsed)\n", m, M, elapsed))
  }
}

total_time <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat(sprintf("\nAll %d repeats done in %.1fs\n\n", M, total_time))


# =============================================================================
# SECTION 3: Average coefficients across repeats
# =============================================================================

beta_avg_A <- colMeans(beta_all_A)
beta_avg_B <- colMeans(beta_all_B)
beta_sd_A  <- apply(beta_all_A, 2, sd)
beta_sd_B  <- apply(beta_all_B, 2, sd)

cat("-- Model A (raw proxy + CRF FE) ------------------------------------\n")
for (i in seq_along(beta_avg_A)) {
  cat(sprintf("  %-20s: avg = %8.4f  sd = %6.4f\n",
              names(beta_avg_A)[i], beta_avg_A[i], beta_sd_A[i]))
}

cat("\n-- Model B (within-CRF demeaned percentile) -------------------------\n")
for (i in seq_along(beta_avg_B)) {
  cat(sprintf("  %-20s: avg = %8.4f  sd = %6.4f\n",
              names(beta_avg_B)[i], beta_avg_B[i], beta_sd_B[i]))
}
cat("\n")


# =============================================================================
# SECTION 4: Build p_hat functions from averaged coefficients
# =============================================================================

# Model A: needs proxy value + CRF group
# Returns P(emitter) for given proxy values and CRF groups
p_hat_A <- function(proxy, crf, beta = beta_avg_A) {
  mm <- model.matrix(~ proxy + factor(crf))
  # Align columns with coefficient names
  eta <- mm %*% beta[match(colnames(mm), names(beta))]
  as.numeric(plogis(eta))
}

# Model B: needs percentile (within-CRF of year-demeaned proxy)
p_hat_B <- function(percentile, beta = beta_avg_B) {
  eta <- beta["intercept"] + beta["percentile"] * percentile
  as.numeric(plogis(eta))
}


# =============================================================================
# SECTION 5: Validation — reproduce aggregate FPR / FNR
# =============================================================================

cat("-- Validation: aggregate FPR / FNR ----------------------------------\n\n")

emit_val <- panel$emit[mixed_idx]

# Per-repeat FPR/TPR (matches inferring_emissions table methodology)
fpr_vec <- numeric(M)
tpr_vec <- numeric(M)
for (m in seq_len(M)) {
  proxy_m <- proxy_matrix[mixed_idx, m]
  pred_m  <- as.integer(proxy_m > 0)
  tp <- sum(pred_m == 1 & emit_val == 1)
  fp <- sum(pred_m == 1 & emit_val == 0)
  tn <- sum(pred_m == 0 & emit_val == 0)
  fn <- sum(pred_m == 0 & emit_val == 1)
  tpr_vec[m] <- tp / (tp + fn)
  fpr_vec[m] <- fp / (fp + tn)
}

obs_TPR <- mean(tpr_vec)
obs_FPR <- mean(fpr_vec)

cat(sprintf("  Per-repeat FPR/TPR averaged across M = %d repeats (mixed CRFs):\n", M))
cat(sprintf("    avg TPR = %.4f (sd = %.4f)\n", obs_TPR, sd(tpr_vec)))
cat(sprintf("    avg FPR = %.4f (sd = %.4f)\n\n", obs_FPR, sd(fpr_vec)))

# Model A implied probabilities (using proxy_mean for a single validation pass)
proxy_val <- panel$proxy_mean[mixed_idx]
crf_val   <- panel$primary_crf_group[mixed_idx]
year_val  <- panel$year[mixed_idx]

df_val <- data.frame(proxy = proxy_val, crf = factor(crf_val), year = year_val)
mm_val <- model.matrix(~ proxy + crf, data = df_val)
eta_val_A <- mm_val %*% beta_avg_A[match(colnames(mm_val), names(beta_avg_A))]
phat_val_A <- plogis(eta_val_A)

# Model B implied probabilities
df_val$proxy_demeaned <- ave(df_val$proxy, df_val$crf, df_val$year,
                              FUN = function(x) x - mean(x))
df_val$percentile <- ave(df_val$proxy_demeaned, df_val$crf,
                          FUN = function(x) rank(x, ties.method = "average") / length(x))
phat_val_B <- plogis(beta_avg_B["intercept"] + beta_avg_B["percentile"] * df_val$percentile)

implied_TPR_A <- mean(phat_val_A[emit_val == 1])
implied_FPR_A <- mean(phat_val_A[emit_val == 0])
implied_TPR_B <- mean(phat_val_B[emit_val == 1])
implied_FPR_B <- mean(phat_val_B[emit_val == 0])

cat("  Model A implied (avg p_hat among true emitters / non-emitters):\n")
cat(sprintf("    E[p_hat | emitter] = %.4f  E[p_hat | non-emitter] = %.4f\n",
            implied_TPR_A, implied_FPR_A))
cat("\n  Model B implied:\n")
cat(sprintf("    E[p_hat | emitter] = %.4f  E[p_hat | non-emitter] = %.4f\n\n",
            implied_TPR_B, implied_FPR_B))


# =============================================================================
# SECTION 6: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA, "extensive_margin_calibration.RData")

validation <- list(
  obs_TPR = obs_TPR, obs_FPR = obs_FPR,
  tpr_vec = tpr_vec, fpr_vec = fpr_vec,
  implied_TPR_A = implied_TPR_A, implied_FPR_A = implied_FPR_A,
  implied_TPR_B = implied_TPR_B, implied_FPR_B = implied_FPR_B
)

save(
  beta_avg_A, beta_avg_B,
  beta_sd_A, beta_sd_B,
  beta_all_A, beta_all_B,
  MIXED_CRFS,
  validation,
  file = OUT_PATH
)

cat("===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat("  Model A:", length(beta_avg_A), "coefficients (raw proxy + CRF FE)\n")
cat("  Model B:", length(beta_avg_B), "coefficients (demeaned percentile)\n")
cat("  Validation metrics included\n")
cat("===================================================================\n")
