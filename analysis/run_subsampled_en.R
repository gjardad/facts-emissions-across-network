###############################################################################
# analysis/run_subsampled_en.R
#
# PURPOSE
#   Compute an averaged deployment proxy for non-ETS firms via B repeated
#   elastic net fits on stratified subsamples of the training sample.
#
#   For each draw b = 1, ..., B:
#     1. Leave-20%-out of training firms (stratified by NACE 2-digit × emitter)
#     2. Fit elastic net on the 80% subsample
#     3. Apply supplier coefficients to ALL deployment firms (non-ETS, not in
#        training sample)
#   Average and record SD of proxy values across B draws.
#
#   Design closely follows inferring_emissions/analysis/active/
#   build_repeated_cv_proxy_asinh.R, which runs M=200 K-fold CV to evaluate
#   out-of-sample performance for training firms. This script instead computes
#   the deployment proxy for the full population of non-ETS firms.
#
# SPECIFICATION
#   LHS: asinh(emissions) [EU ETS verified, or 0 for confirmed-zero sectors]
#   RHS: asinh(B2B purchases from each eligible supplier) [penalized]
#        log(revenue) + year FEs + NACE 2-digit sector FEs [unpenalized]
#   Eligible suppliers: those with >= MIN_LHS_BUYERS distinct LHS buyers in
#                       the training subsample.
#
# INPUT
#   {PROC_DATA}/b2b_selected_sample.RData   — buyer-seller-year B2B panel
#   {PROC_DATA}/training_sample.RData        — EU ETS + zero-emission non-ETS
#
# OUTPUT
#   {PROC_DATA}/deployment_proxy_avg.RData
#     proxy_avg  : data.frame (vat, year, proxy_avg, proxy_sd) for all
#                  deployment firm-years observed in B2B data
#     proxy_list : list of length B, each a data.frame (vat, year, proxy)
#                  for one subsample draw — used in uncertainty propagation
#
# RUNS ON: RMD (requires b2b_selected_sample with deployment firm purchases)
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
library(Matrix)
library(glmnet)
library(doParallel)
library(foreach)

# ── Parameters ───────────────────────────────────────────────────────────────
B_DRAWS        <- 50L    # number of subsample draws
SUBSAMPLE_FRAC <- 0.8    # fraction of training firms kept in each subsample
K_INNER        <- 10L    # inner folds for cv.glmnet lambda tuning
MIN_LHS_BUYERS <- 5L     # sellers must supply >= MIN_LHS_BUYERS LHS buyers
ALPHA          <- 0.5    # elastic net mixing parameter (0.5 = true EN)
BASE_SEED      <- 2026L  # draw b uses seed = BASE_SEED + b
N_CORES        <- 40L    # parallel workers (set to 1 for sequential)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  SUBSAMPLED EN DEPLOYMENT PROXY\n")
cat("  B =", B_DRAWS, "draws | subsample =", SUBSAMPLE_FRAC,
    "| alpha =", ALPHA, "| cores =", N_CORES, "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")


# =============================================================================
# STEP 1: Load data
# =============================================================================
cat("Loading B2B selected sample...\n")
load(file.path(PROC_DATA, "b2b_selected_sample.RData"))
b2b <- df_b2b_selected_sample
rm(df_b2b_selected_sample)

cat("Loading training sample...\n")
load(file.path(PROC_DATA, "training_sample.RData"))

lhs <- training_sample %>%
  filter(year >= 2005) %>%
  mutate(
    y           = emissions,
    log_revenue = log(pmax(revenue, 1e-12))
  ) %>%
  select(vat, year, y, log_revenue, nace2d, euets) %>%
  arrange(vat, year)
rm(training_sample)

lhs$emit <- as.integer(lhs$y > 0)
cat("Training LHS:", nrow(lhs), "firm-years,",
    length(unique(lhs$vat)), "unique firms\n")
cat("  Emitter firm-years:", sum(lhs$emit), "\n\n")


# =============================================================================
# STEP 2: Assign primary sector to each training firm
# =============================================================================
firm_sector <- lhs %>%
  group_by(vat) %>%
  summarise(primary_nace2d = names(which.max(table(nace2d))), .groups = "drop")

lhs <- lhs %>%
  left_join(firm_sector, by = "vat")

training_vats <- unique(lhs$vat)
cat("Sectors in training sample:", length(unique(lhs$primary_nace2d)), "\n\n")


# =============================================================================
# STEP 3: Split B2B into training-buyer and deployment-buyer panels
# =============================================================================
cat("Splitting B2B into training and deployment panels (years >= 2005)...\n")
b2b_filtered <- b2b %>%
  filter(year >= 2005)

deploy_vats <- setdiff(unique(b2b_filtered$vat_j_ano), training_vats)
cat("  Training buyers:", length(training_vats), "\n")
cat("  Deployment buyers:", length(deploy_vats), "\n")

b2b_lhs    <- b2b_filtered[b2b_filtered$vat_j_ano %in% training_vats, ]
b2b_deploy <- b2b_filtered[b2b_filtered$vat_j_ano %in% deploy_vats, ]
rm(b2b, b2b_filtered)

cat("  B2B training rows:", nrow(b2b_lhs),
    "| B2B deployment rows:", nrow(b2b_deploy), "\n\n")


# =============================================================================
# STEP 4: Helper functions
# =============================================================================

# Extract non-zero supplier coefficients from a cv.glmnet fit
extract_suppliers <- function(fit, n_ctrl, eligible_sellers, s = "lambda.min") {
  co           <- coef(fit, s = s)
  supplier_idx <- (n_ctrl + 2):length(co)  # skip intercept + controls
  vals         <- co[supplier_idx]
  data.frame(
    vat_i_ano = eligible_sellers,
    coef      = as.numeric(vals),
    stringsAsFactors = FALSE
  ) %>%
    filter(coef != 0) %>%
    arrange(desc(abs(coef)))
}

# Stratified leave-(1-frac)-out subsample of training firms
draw_subsample_vats <- function(lhs, frac, seed) {
  set.seed(seed)
  lhs %>%
    distinct(vat, primary_nace2d, emit) %>%
    group_by(primary_nace2d, emit) %>%
    slice_sample(prop = frac) %>%
    ungroup() %>%
    pull(vat)
}


# =============================================================================
# STEP 5: Single-draw function
# =============================================================================
run_one_subsample <- function(lhs, b2b_lhs, b2b_deploy,
                               alpha, K_inner, min_buyers,
                               subsample_frac, seed) {

  # --- 5.1 Stratified subsample of training firms ---
  subsample_vats <- draw_subsample_vats(lhs, subsample_frac, seed)
  train_lhs      <- lhs[lhs$vat %in% subsample_vats, ]
  train_lhs$row_idx <- seq_len(nrow(train_lhs))

  # --- 5.2 Eligible sellers: must supply >= min_buyers distinct LHS buyers ---
  b2b_train    <- b2b_lhs[b2b_lhs$vat_j_ano %in% train_lhs$vat, ]
  seller_counts <- b2b_train %>%
    distinct(vat_i_ano, vat_j_ano) %>%
    count(vat_i_ano, name = "n_lhs_buyers")
  eligible_sellers <- seller_counts %>%
    filter(n_lhs_buyers >= min_buyers) %>%
    pull(vat_i_ano) %>%
    sort()

  if (length(eligible_sellers) == 0) {
    return(data.frame(vat = character(0), year = integer(0),
                      proxy = numeric(0), stringsAsFactors = FALSE))
  }

  seller_map <- data.frame(
    vat_i_ano = eligible_sellers,
    col_idx   = seq_along(eligible_sellers),
    stringsAsFactors = FALSE
  )

  # --- 5.3 Sparse B2B supplier matrix (training firms × eligible sellers) ---
  b2b_agg <- b2b_train %>%
    filter(vat_i_ano %in% eligible_sellers) %>%
    group_by(vat_i_ano, vat_j_ano, year) %>%
    summarise(sales = sum(corr_sales_ij, na.rm = TRUE), .groups = "drop") %>%
    filter(sales != 0) %>%
    inner_join(train_lhs %>% select(vat, year, row_idx),
               by = c("vat_j_ano" = "vat", "year" = "year")) %>%
    inner_join(seller_map, by = "vat_i_ano")

  n_rows <- nrow(train_lhs)
  n_cols <- length(eligible_sellers)

  X_asinh <- sparseMatrix(
    i = b2b_agg$row_idx, j = b2b_agg$col_idx,
    x = asinh(b2b_agg$sales), dims = c(n_rows, n_cols)
  )
  colnames(X_asinh) <- eligible_sellers
  rm(b2b_agg, b2b_train)

  # --- 5.4 Control matrix: log_revenue + year FEs + sector FEs ---
  year_dummies <- model.matrix(~ factor(year), data = train_lhs)[, -1, drop = FALSE]
  colnames(year_dummies) <- paste0("yr_", sort(unique(train_lhs$year))[-1])

  train_sectors  <- sort(unique(train_lhs$nace2d))
  sector_dummies <- model.matrix(~ factor(nace2d, levels = train_sectors),
                                  data = train_lhs)[, -1, drop = FALSE]
  colnames(sector_dummies) <- paste0("sec_", train_sectors[-1])

  X_controls <- cbind(log_revenue = train_lhs$log_revenue,
                       year_dummies, sector_dummies)
  n_controls <- ncol(X_controls)

  X_full <- cbind(Matrix(X_controls, sparse = TRUE), X_asinh)
  pf     <- c(rep(0, n_controls), rep(1, n_cols))  # 0 = unpenalized
  rm(X_asinh, X_controls, year_dummies, sector_dummies)

  # Drop zero-variance columns
  col_vars <- colMeans(X_full^2) - colMeans(X_full)^2
  zero_var <- which(col_vars == 0)
  if (length(zero_var) > 0) {
    X_full <- X_full[, -zero_var]
    pf     <- pf[-zero_var]
    n_controls <- sum(pf == 0)
    eligible_sellers_enet <- colnames(X_full)[which(pf == 1)]
  } else {
    eligible_sellers_enet <- eligible_sellers
  }

  # --- 5.5 Inner folds for lambda tuning (grouped by firm) ---
  set.seed(seed + 1000L)
  firms_b        <- unique(train_lhs$vat)
  inner_fold_map <- sample(rep(1:K_inner, length.out = length(firms_b)))
  names(inner_fold_map) <- firms_b
  inner_foldid <- unname(inner_fold_map[train_lhs$vat])

  # --- 5.6 Fit cv.glmnet ---
  fit <- cv.glmnet(
    x = X_full, y = asinh(train_lhs$y),
    family = "gaussian", alpha = alpha,
    penalty.factor = pf, foldid = inner_foldid,
    standardize = TRUE
  )
  rm(X_full)

  # --- 5.7 Extract positive supplier coefficients at lambda.min ---
  coef_lookup <- extract_suppliers(fit, n_controls, eligible_sellers_enet,
                                    "lambda.min")
  rm(fit)

  if (nrow(coef_lookup) == 0) {
    return(data.frame(vat = character(0), year = integer(0),
                      proxy = numeric(0), stringsAsFactors = FALSE))
  }

  # --- 5.8 Apply supplier coefficients to deployment firms ---
  # Proxy = sum of coef_s * asinh(B2B purchase from seller s), floored at 0
  proxy <- b2b_deploy %>%
    inner_join(coef_lookup, by = "vat_i_ano") %>%
    group_by(vat_j_ano, year) %>%
    summarise(proxy = sum(coef * asinh(corr_sales_ij), na.rm = TRUE),
              .groups = "drop") %>%
    rename(vat = vat_j_ano) %>%
    mutate(proxy = pmax(proxy, 0))

  proxy
}


# =============================================================================
# STEP 6: Main loop over B draws
# =============================================================================
use_parallel <- N_CORES > 1

if (use_parallel) {
  cat("\n═══ Starting", B_DRAWS, "draws in PARALLEL (", N_CORES, "cores) ═══\n\n")

  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(dplyr)
    library(Matrix)
    library(glmnet)
  })

  t0_all <- Sys.time()

  all_results <- foreach(
    b = seq_len(B_DRAWS),
    .packages = c("dplyr", "Matrix", "glmnet"),
    .export = c("lhs", "b2b_lhs", "b2b_deploy",
                "ALPHA", "K_INNER", "MIN_LHS_BUYERS", "SUBSAMPLE_FRAC",
                "BASE_SEED", "extract_suppliers", "draw_subsample_vats",
                "run_one_subsample")
  ) %dopar% {
    t0_b <- Sys.time()
    seed_b <- BASE_SEED + b
    proxy_b <- run_one_subsample(lhs, b2b_lhs, b2b_deploy,
                                  ALPHA, K_INNER, MIN_LHS_BUYERS,
                                  SUBSAMPLE_FRAC, seed_b)
    list(proxy = proxy_b,
         time  = as.numeric(difftime(Sys.time(), t0_b, units = "mins")))
  }

  stopCluster(cl)

  total_time  <- round(difftime(Sys.time(), t0_all, units = "mins"), 1)
  proxy_list  <- lapply(all_results, `[[`, "proxy")
  draw_timing <- sapply(all_results, `[[`, "time")
  rm(all_results)

  cat(sprintf("═══ All %d draws complete in %.1f min (avg %.1f min/draw) ═══\n\n",
              B_DRAWS, total_time, mean(draw_timing)))

} else {
  cat("\n═══ Starting", B_DRAWS, "draws SEQUENTIALLY ═══\n\n")

  t0_all      <- Sys.time()
  proxy_list  <- vector("list", B_DRAWS)
  draw_timing <- numeric(B_DRAWS)

  for (b in seq_len(B_DRAWS)) {
    t0_b   <- Sys.time()
    seed_b <- BASE_SEED + b
    cat(sprintf("── Draw %d / %d ", b, B_DRAWS))

    proxy_list[[b]] <- run_one_subsample(lhs, b2b_lhs, b2b_deploy,
                                          ALPHA, K_INNER, MIN_LHS_BUYERS,
                                          SUBSAMPLE_FRAC, seed_b)
    gc()

    elapsed_b      <- as.numeric(difftime(Sys.time(), t0_b, units = "mins"))
    draw_timing[b] <- elapsed_b
    avg_time       <- mean(draw_timing[1:b])
    remaining      <- round(avg_time * (B_DRAWS - b), 1)
    cat(sprintf("(%.1f min, est. %.0f min remaining) ──\n", elapsed_b, remaining))
  }

  total_time <- round(difftime(Sys.time(), t0_all, units = "mins"), 1)
  cat(sprintf("\n═══ All %d draws complete in %.1f min (avg %.1f min/draw) ═══\n\n",
              B_DRAWS, total_time, mean(draw_timing)))
}


# =============================================================================
# STEP 7: Average proxy across B draws
# =============================================================================
cat("Averaging proxy across B draws...\n")

# Canonical deployment firm-year grid: union of all (vat, year) pairs
all_deploy_fy <- b2b_deploy %>%
  distinct(vat_j_ano, year) %>%
  rename(vat = vat_j_ano)

# For each draw, fill missing firm-years with 0 (no supplier matches → proxy = 0)
proxy_matrix <- do.call(cbind, lapply(proxy_list, function(px) {
  all_deploy_fy %>%
    left_join(px, by = c("vat", "year")) %>%
    mutate(proxy = coalesce(proxy, 0)) %>%
    pull(proxy)
}))

proxy_avg_vec <- rowMeans(proxy_matrix)
proxy_sd_vec  <- apply(proxy_matrix, 1, sd)
rm(proxy_matrix)

proxy_avg <- all_deploy_fy %>%
  mutate(
    proxy_avg = proxy_avg_vec,
    proxy_sd  = proxy_sd_vec
  )

cat("Deployment proxy panel:", nrow(proxy_avg), "firm-years,",
    n_distinct(proxy_avg$vat), "firms\n")
cat("  Predicted emitters (proxy_avg > 0):",
    sum(proxy_avg$proxy_avg > 0),
    sprintf("(%.1f%%)\n", 100 * mean(proxy_avg$proxy_avg > 0)))
cat("  Mean proxy SD:", round(mean(proxy_avg$proxy_sd, na.rm = TRUE), 4), "\n\n")


# =============================================================================
# STEP 8: Save
# =============================================================================
OUT_PATH <- file.path(PROC_DATA, "deployment_proxy_avg.RData")

save(proxy_avg, proxy_list, file = OUT_PATH)

cat("══════════════════════════════════════════════\n")
cat("Saved to:", OUT_PATH, "\n")
cat("  proxy_avg :", nrow(proxy_avg), "rows x", ncol(proxy_avg), "cols\n")
cat("  proxy_list:", length(proxy_list), "draws\n")
cat("  B =", B_DRAWS, "| alpha =", ALPHA, "| seed =", BASE_SEED, "\n")
cat("══════════════════════════════════════════════\n")
