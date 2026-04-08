###############################################################################
# analysis/calibrate_pooled_youden_tau.R
#
# PURPOSE
#   Calibrate the deployment threshold tau* on p_i = share of CV repeats with
#   positive proxy. Pool the three mixed CRF groups (paper, refining, metals)
#   and pick tau by maximizing Youden's J:
#
#       J(tau) = TPR(tau) - FPR(tau)   on (p_i, true emit) pairs
#
#   The same tau* is then applied at deployment time to all CRF groups.
#
# INPUT
#   {PROC_DATA}/repeated_cv_proxy_crf_asinh.RData
#     proxy_matrix            : N x M matrix of training-CV proxies
#     repeated_cv_proxy_panel : data.frame with vat, year, emit,
#                               primary_crf_group, ...
#
# OUTPUT
#   {PROC_DATA}/pooled_youden_tau.RData
#     tau_star : numeric, optimal threshold on p_i
#     J_star   : numeric, Youden's J at tau_star
#     n_pooled : integer, n firm-years used in calibration
#     tpr      : TPR at tau_star
#     fpr      : FPR at tau_star
#
# RUNS ON: local 1 or RMD.
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
cat("  POOLED-MIXED YOUDEN tau CALIBRATION\n")
cat("===================================================================\n\n")

MIXED_CRFS <- c("paper", "refining", "metals")

# -- Load training panel ------------------------------------------------------
load(file.path(PROC_DATA, "repeated_cv_proxy_crf_asinh.RData"))

stopifnot(nrow(proxy_matrix) == nrow(repeated_cv_proxy_panel))
cat("Loaded:", nrow(proxy_matrix), "firm-years x", ncol(proxy_matrix), "repeats\n\n")

# -- Compute p_i --------------------------------------------------------------
n_obs <- rowSums(!is.na(proxy_matrix))
n_pos <- rowSums(proxy_matrix > 0, na.rm = TRUE)
p_i   <- ifelse(n_obs > 0, n_pos / n_obs, NA_real_)

panel <- repeated_cv_proxy_panel
panel$p_i <- p_i

# -- Restrict to mixed CRFs ---------------------------------------------------
mixed <- panel[panel$primary_crf_group %in% MIXED_CRFS &
                 !is.na(panel$p_i) & !is.na(panel$emit), ]

cat("Mixed-CRF pooled sample:\n")
for (g in MIXED_CRFS) {
  sub <- mixed[mixed$primary_crf_group == g, ]
  cat(sprintf("  %-10s : %5d firm-years, %4d emitters (%.1f%%)\n",
              g, nrow(sub), sum(sub$emit == 1),
              100 * mean(sub$emit == 1)))
}
cat(sprintf("  pooled     : %5d firm-years, %4d emitters\n\n",
            nrow(mixed), sum(mixed$emit == 1)))

# -- Maximize Youden's J ------------------------------------------------------
emit <- as.integer(mixed$emit == 1)
p    <- mixed$p_i
P    <- sum(emit == 1)
N    <- sum(emit == 0)
if (P == 0 || N == 0) stop("Need both emitters and non-emitters in pooled sample.")

# Sort descending by p; sweep cumulative TPR/FPR
ord  <- order(p, decreasing = TRUE)
p_o  <- p[ord]
e_o  <- emit[ord]
tp_cum <- cumsum(e_o == 1)
fp_cum <- cumsum(e_o == 0)
tpr    <- tp_cum / P
fpr    <- fp_cum / N
j_vec  <- tpr - fpr

# Evaluate only at the last index per distinct p value (i.e. take >=)
last_at_value <- !duplicated(p_o, fromLast = TRUE)
p_eval <- p_o[last_at_value]
j_eval <- j_vec[last_at_value]
tpr_eval <- tpr[last_at_value]
fpr_eval <- fpr[last_at_value]

best_idx <- which.max(j_eval)
tau_star <- p_eval[best_idx]
J_star   <- j_eval[best_idx]
tpr_star <- tpr_eval[best_idx]
fpr_star <- fpr_eval[best_idx]
n_pooled <- nrow(mixed)

cat(sprintf("tau_star = %.4f\n", tau_star))
cat(sprintf("J_star   = %.4f  (TPR = %.4f, FPR = %.4f)\n",
            J_star, tpr_star, fpr_star))
cat(sprintf("n_pooled = %d\n", n_pooled))

# -- Save ---------------------------------------------------------------------
tpr <- tpr_star
fpr <- fpr_star
out_path <- file.path(PROC_DATA, "pooled_youden_tau.RData")
save(tau_star, J_star, n_pooled, tpr, fpr, file = out_path)
cat("\nSaved:", out_path, "\n")
