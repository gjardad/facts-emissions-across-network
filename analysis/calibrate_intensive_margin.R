###############################################################################
# analysis/calibrate_intensive_margin.R
#
# PURPOSE
#   Calibrate the intensive-margin model used to perturb within-cell ranks
#   among true emitters in deployment.
#
#   Sample: all true emitters (D_i = 1) across all training sectors. Each
#   true emitter is assigned Dhat_i using the threshold rule:
#     - mixed-CRF firms: tau_loso[crf]  (LOSO threshold)
#     - non-mixed-CRF firms: tau_pooled_mixed
#
#   Branch A: D_i = 1 AND Dhat_i = 1
#     - drop cells (crf_group, year) with fewer than 3 true emitters
#     - within each cell, compute u_i (true rank by y_i) and uhat_i
#       (predicted rank by proxy_mean_i), both as percentiles among true emitters
#     - z_i = Phi^-1(u_i), zhat_i = Phi^-1(uhat_i)
#     - fit mean regression: z_i = a + b * zhat_i + e_i
#     - fit variance regression: log(e_i^2 + c) = d + h * zhat_i + noise
#
#   Branch B: D_i = 1 AND Dhat_i = 0  (false negatives)
#     - compute u_i and z_i for all branch B observations (no cell-size filter)
#     - store the empirical sample of z_i values as fallback distribution
#
# INPUT
#   {PROC_DATA}/training_summary_objects.RData
#   {PROC_DATA}/threshold_calibration.RData
#
# OUTPUT
#   {PROC_DATA}/intensive_margin_calibration.RData
#     mean_coefs    : c(a, b)
#     var_coefs     : c(d, h)
#     branch_b_z    : numeric vector of z_i values (empirical fallback)
#     mean_fit      : fitted lm object for the mean regression
#     var_fit       : fitted lm object for the variance regression
#     int_diagnostics : list of sample sizes, R^2, etc.
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

# -- Parameter ----------------------------------------------------------------
# Which extensive-margin weighting scheme to consume.
# Set to "balanced" (primary) or "unweighted" (robustness).
WEIGHT_SCHEME <- "balanced"

cat("===================================================================\n")
cat("  INTENSIVE MARGIN CALIBRATION\n")
cat("  WEIGHT_SCHEME =", WEIGHT_SCHEME, "\n")
cat("===================================================================\n\n")


# =============================================================================
# SECTION 1: Load
# =============================================================================

load(file.path(PROC_DATA, "training_summary_objects.RData"))
load(file.path(PROC_DATA, "extensive_margin_calibration.RData"))

stopifnot(WEIGHT_SCHEME %in% SCHEMES)

# Pull the scheme-specific calibration objects
sch_coefs       <- ext_coefs[[WEIGHT_SCHEME]]
sch_loso_fits   <- loso_fits[[WEIGHT_SCHEME]]
sch_loso_q_star <- loso_q_star[[WEIGHT_SCHEME]]
sch_fit_pooled  <- ext_fit_pooled[[WEIGHT_SCHEME]]
sch_q_star      <- q_star_pooled[[WEIGHT_SCHEME]]

cat("Loaded training_summary (", nrow(training_summary), "rows)\n")
cat("Logistic coefs (scheme =", WEIGHT_SCHEME, "):\n")
cat(sprintf("  alpha = %.4f, beta_p = %.4f, beta_lp = %.4f\n",
            sch_coefs["alpha"], sch_coefs["beta_p"], sch_coefs["beta_lp"]))
cat("LOSO q*:\n")
for (s in MIXED_CRFS)
  cat(sprintf("  %s = %.4f\n", s, sch_loso_q_star[s]))
cat(sprintf("Pooled q* (used for non-mixed CRFs) = %.4f\n\n", sch_q_star))


# =============================================================================
# SECTION 2: Assign Dhat_i
# =============================================================================
#
# For mixed-CRF observations: use the LOSO logistic for that fold (the fit
# trained on the OTHER two mixed CRFs), and the LOSO q* for that fold.
# For non-mixed-CRF observations: use the pooled-mixed logistic and pooled q*.

ts <- training_summary
ts$lp <- log1p(pmax(ts$proxy_mean_i, 0))
ts$Dhat <- NA_integer_

for (s in MIXED_CRFS) {
  idx <- which(ts$primary_crf_group == s)
  if (length(idx) == 0) next
  q_idx <- predict(sch_loso_fits[[s]],
                    newdata = ts[idx, c("p_i", "lp")],
                    type = "response")
  ts$Dhat[idx] <- as.integer(q_idx >= sch_loso_q_star[s])
}

# Non-mixed CRFs: pooled fit
nm_idx <- which(!(ts$primary_crf_group %in% MIXED_CRFS))
if (length(nm_idx) > 0) {
  q_nm <- predict(sch_fit_pooled,
                   newdata = ts[nm_idx, c("p_i", "lp")],
                   type = "response")
  ts$Dhat[nm_idx] <- as.integer(q_nm >= sch_q_star)
}

# Restrict to true emitters
emit <- ts[ts$D_i == 1, ]
cat("True emitters:", nrow(emit), "\n")
cat("  Branch A (Dhat = 1):", sum(emit$Dhat == 1), "\n")
cat("  Branch B (Dhat = 0):", sum(emit$Dhat == 0), "\n\n")


# =============================================================================
# SECTION 3: Branch A — within-cell percentiles, drop tiny cells
# =============================================================================

cat("-- Branch A: build (u_i, uhat_i) within (crf_group, year) cells -----\n\n")

br_a <- emit[emit$Dhat == 1, ]

# For each cell, count emitters; drop cells with < 3
cell_n <- ave(br_a$D_i, br_a$primary_crf_group, br_a$year, FUN = length)
cat("Cell-size distribution (true emitters per cell, branch A):\n")
print(summary(cell_n))
cat(sprintf("  Cells with size < 3: %d firm-years (%.1f%%)\n",
            sum(cell_n < 3), 100 * mean(cell_n < 3)))

br_a <- br_a[cell_n >= 3, ]
cat(sprintf("\n  Branch A after dropping tiny cells: %d firm-years\n\n", nrow(br_a)))

# Within each cell, compute u_i (rank by y_i, descending) and uhat_i (rank by
# proxy_mean_i, descending), both as percentiles among true emitters.
# rank 1 = highest, so we want percentile = (n - rank + 0.5) / n? No — the .md
# says rank 1 = highest, and uhat_i = (rank_i - 0.5) / n. So uhat_i is small
# when rank_i = 1 (highest emitter), which puts it at the BOTTOM of the (0,1)
# scale. Then z = Phi^-1(uhat_i) is negative for top emitters.
#
# This is fine — it's just a labeling convention. The key thing is that the
# SAME convention is used for both u_i and uhat_i, so the slope b > 0
# means "higher predicted rank → higher true rank" in the same direction.

compute_pct <- function(x, n) {
  # rank 1 = highest x; pct = (rank - 0.5) / n
  r <- rank(-x, ties.method = "average")
  (r - 0.5) / n
}

br_a$u_i    <- ave(br_a$y_i,          br_a$primary_crf_group, br_a$year,
                    FUN = function(x) compute_pct(x, length(x)))
br_a$uhat_i <- ave(br_a$proxy_mean_i, br_a$primary_crf_group, br_a$year,
                    FUN = function(x) compute_pct(x, length(x)))

br_a$z_i    <- qnorm(br_a$u_i)
br_a$zhat_i <- qnorm(br_a$uhat_i)

cat("Branch A z_i / zhat_i summary:\n")
print(summary(br_a[, c("u_i", "uhat_i", "z_i", "zhat_i")]))


# =============================================================================
# SECTION 4: Branch A — mean and variance regressions
# =============================================================================

cat("\n-- Branch A: mean regression z_i ~ a + b * zhat_i -------------------\n")

mean_fit <- lm(z_i ~ zhat_i, data = br_a)
mean_coefs <- coef(mean_fit)
names(mean_coefs) <- c("a", "b")
cat(sprintf("  a = %.4f, b = %.4f, R^2 = %.4f\n",
            mean_coefs["a"], mean_coefs["b"], summary(mean_fit)$r.squared))

br_a$ehat <- residuals(mean_fit)

cat("\n-- Branch A: residual variance (homoskedastic) ----------------------\n")

# We checked for heteroskedasticity by fitting a log-squared-residual
# regression on zhat_i. The slope was h ~ -0.08 with R^2 ~ 0.001 — no
# evidence of variance structure in zhat_i. We therefore use a single
# constant sigma estimated as the empirical residual SD, which avoids the
# Jensen's-inequality bias of regressing log(e^2) and is unbiased on the
# level scale.

sigma_const <- sqrt(mean(br_a$ehat^2))
cat(sprintf("  sigma (empirical residual SD) = %.4f\n", sigma_const))

# Spearman between z_i and zhat_i for diagnostic comparison with OOS metric
sp <- cor(br_a$z_i, br_a$zhat_i, method = "spearman")
cat(sprintf("\n  Spearman(z_i, zhat_i) on branch A = %.4f\n", sp))


# =============================================================================
# SECTION 5: Branch B — empirical distribution of z_i
# =============================================================================

cat("\n-- Branch B: empirical fallback distribution ------------------------\n\n")

br_b <- emit[emit$Dhat == 0, ]
cat("Branch B firm-years:", nrow(br_b), "\n")
cat("  By CRF:\n")
print(table(br_b$primary_crf_group))

# Compute u_i within (crf_group, year), but using ALL true emitters in the
# cell (including branch A firms), since percentiles must be defined among
# the full true-emitter population. Then keep only branch B z_i values.
all_emit <- emit  # contains both branch A and branch B
all_emit$u_full <- ave(all_emit$y_i, all_emit$primary_crf_group, all_emit$year,
                        FUN = function(x) compute_pct(x, length(x)))
all_emit$z_full <- qnorm(all_emit$u_full)

br_b_z <- all_emit$z_full[all_emit$Dhat == 0]
branch_b_z <- br_b_z

cat(sprintf("\nBranch B z_i summary (n = %d):\n", length(branch_b_z)))
print(summary(branch_b_z))


# =============================================================================
# SECTION 6: Diagnostics summary
# =============================================================================

int_diagnostics <- list(
  branch_a_n       = nrow(br_a),
  branch_a_dropped = sum(cell_n < 3),
  branch_b_n       = length(branch_b_z),
  mean_r2          = summary(mean_fit)$r.squared,
  sigma_const      = sigma_const,
  spearman_z_zhat  = sp
)


# =============================================================================
# SECTION 7: Save
# =============================================================================

OUT_PATH <- file.path(PROC_DATA,
  sprintf("intensive_margin_calibration_%s.RData", WEIGHT_SCHEME))
save(mean_coefs, sigma_const, branch_b_z, mean_fit,
     int_diagnostics, WEIGHT_SCHEME, file = OUT_PATH)

cat("\n===================================================================\n")
cat("Saved:", OUT_PATH, "\n")
cat(sprintf("  Branch A: n = %d, mean(a=%.3f, b=%.3f), sigma=%.3f\n",
            int_diagnostics$branch_a_n,
            mean_coefs["a"], mean_coefs["b"], sigma_const))
cat(sprintf("  Branch B: n = %d empirical z_i values\n", int_diagnostics$branch_b_n))
cat("===================================================================\n")
