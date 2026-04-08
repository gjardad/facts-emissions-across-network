###############################################################################
# analysis/run_subsampled_en.R
#
# PURPOSE
#   Compute B deployment proxy vectors for non-ETS firms via repeated
#   elastic net fits on stratified subsamples of the training sample.
#
#   For each draw b = 1, ..., B:
#     1. Leave-20%-out of training firms (stratified by NACE 2-digit × emitter)
#     2. Fit elastic net on the 80% subsample
#     3. Apply supplier coefficients to ALL deployment firms (non-ETS, not in
#        training sample)
#
#   The B proxy vectors are the primary output. Prediction intervals are
#   constructed downstream by nesting K perturbations within each draw b
#   (see analysis/uncertainty_propagation.R).
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
#   {PROC_DATA}/deployment_firm_years.RData
#     deployment_firm_years : data.frame (vat, year) defining the deployment grid
#   {PROC_DATA}/deployment_proxy_list.RData
#     proxy_list : list of length B, each a data.frame (vat, year, proxy)
#                  for one subsample draw
#   {PROC_DATA}/deployment_draw_params.RData
#     draw_params : list of length B, each containing gpa_params
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
source(file.path(REPO_DIR, "utils", "calibration_pareto.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)
library(Matrix)
library(glmnet)
library(doParallel)
library(foreach)

# ── Parameters ───────────────────────────────────────────────────────────────
B_DRAWS        <- 200L   # number of subsample draws
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
# Year-specific split: a firm that enters ETS in year T is in training for
# years >= T but in deployment for years < T. This avoids dropping pre-entry
# years of late ETS entrants from both training and deployment.
cat("Splitting B2B into training and deployment panels (years >= 2005)...\n")
b2b_filtered <- b2b %>%
  filter(year >= 2005)

# Build year-specific training firm-year set
training_vy <- lhs %>%
  distinct(vat, year)

b2b_lhs <- b2b_filtered %>%
  inner_join(training_vy, by = c("vat_j_ano" = "vat", "year" = "year"))

b2b_deploy <- b2b_filtered %>%
  anti_join(training_vy, by = c("vat_j_ano" = "vat", "year" = "year"))

# For reporting: permanent training VATs vs late-entrant deployment firm-years
deploy_vats <- unique(b2b_deploy$vat_j_ano)
late_entrant_vy <- b2b_deploy %>%
  filter(vat_j_ano %in% training_vats) %>%
  distinct(vat_j_ano, year)

cat("  Training buyers (ever in ETS/zero-sector):", length(training_vats), "\n")
cat("  Deployment buyers (total unique):", length(deploy_vats), "\n")
cat("  Late-entrant pre-ETS firm-years in deployment:",
    nrow(late_entrant_vy), "\n")

rm(b2b, b2b_filtered, training_vy, late_entrant_vy)

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
    return(list(
      proxy      = data.frame(vat = character(0), year = integer(0),
                               proxy = numeric(0), stringsAsFactors = FALSE),
      gpa_params = NULL,
      threshold_pct = NA_real_
    ))
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

  # --- 5.6 Stage B: Fit cv.glmnet on full subsample (deployment coefficients) -
  y_vec <- asinh(train_lhs$y)

  fit <- cv.glmnet(
    x = X_full, y = y_vec,
    family = "gaussian", alpha = alpha,
    penalty.factor = pf, foldid = inner_foldid,
    standardize = TRUE
  )

  # --- 5.7 Extract supplier coefficients at lambda.min (for deployment) ---
  coef_lookup <- extract_suppliers(fit, n_controls, eligible_sellers_enet,
                                    "lambda.min")

  if (nrow(coef_lookup) == 0) {
    rm(fit, X_full)
    return(list(
      proxy      = data.frame(vat = character(0), year = integer(0),
                               proxy = numeric(0), stringsAsFactors = FALSE),
      gpa_params = NULL,
      threshold_pct = NA_real_
    ))
  }

  # --- 5.8 Apply supplier coefficients to deployment firms ---
  proxy <- b2b_deploy %>%
    inner_join(coef_lookup, by = "vat_i_ano") %>%
    group_by(vat_j_ano, year) %>%
    summarise(proxy = sum(coef * asinh(corr_sales_ij), na.rm = TRUE),
              .groups = "drop") %>%
    rename(vat = vat_j_ano) %>%
    mutate(proxy = pmax(proxy, 0))

  # --- 5.8b GPA reference distribution from subsample's ETS emitters ---
  ets_emit <- train_lhs[train_lhs$euets == 1 & train_lhs$y > 0, ]
  ref_dist   <- build_reference_dist(ets_emit$y, ets_emit$year, ets_emit$nace2d)
  gpa_params <- fit_gpa(ref_dist)

  rm(fit, X_full)

  list(proxy = proxy, gpa_params = gpa_params)
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
    library(lmom)
  })

  t0_all <- Sys.time()

  all_results <- foreach(
    b = seq_len(B_DRAWS),
    .packages = c("dplyr", "Matrix", "glmnet", "lmom"),
    .export = c("lhs", "b2b_lhs", "b2b_deploy",
                "ALPHA", "K_INNER", "MIN_LHS_BUYERS", "SUBSAMPLE_FRAC",
                "BASE_SEED", "extract_suppliers", "draw_subsample_vats",
                "run_one_subsample",
                "build_reference_dist", "fit_gpa")
  ) %dopar% {
    t0_b <- Sys.time()
    seed_b <- BASE_SEED + b
    res_b <- run_one_subsample(lhs, b2b_lhs, b2b_deploy,
                                ALPHA, K_INNER, MIN_LHS_BUYERS,
                                SUBSAMPLE_FRAC, seed_b)
    list(proxy      = res_b$proxy,
         gpa_params = res_b$gpa_params,
         time       = as.numeric(difftime(Sys.time(), t0_b, units = "mins")))
  }

  stopCluster(cl)

  total_time  <- round(difftime(Sys.time(), t0_all, units = "mins"), 1)
  proxy_list  <- lapply(all_results, `[[`, "proxy")
  draw_params <- lapply(all_results, function(r) list(gpa_params = r$gpa_params))
  draw_timing <- sapply(all_results, `[[`, "time")
  rm(all_results)

  cat(sprintf("═══ All %d draws complete in %.1f min (avg %.1f min/draw) ═══\n\n",
              B_DRAWS, total_time, mean(draw_timing)))

} else {
  cat("\n═══ Starting", B_DRAWS, "draws SEQUENTIALLY ═══\n\n")

  t0_all      <- Sys.time()
  proxy_list  <- vector("list", B_DRAWS)
  draw_params <- vector("list", B_DRAWS)
  draw_timing <- numeric(B_DRAWS)

  for (b in seq_len(B_DRAWS)) {
    t0_b   <- Sys.time()
    seed_b <- BASE_SEED + b
    cat(sprintf("── Draw %d / %d ", b, B_DRAWS))

    res_b <- run_one_subsample(lhs, b2b_lhs, b2b_deploy,
                                ALPHA, K_INNER, MIN_LHS_BUYERS,
                                SUBSAMPLE_FRAC, seed_b)
    proxy_list[[b]]  <- res_b$proxy
    draw_params[[b]] <- list(gpa_params = res_b$gpa_params)
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
# STEP 7: Build deployment firm-year grid and save
# =============================================================================

# Canonical deployment firm-year grid: union of all (vat, year) pairs
deployment_firm_years <- b2b_deploy %>%
  distinct(vat_j_ano, year) %>%
  rename(vat = vat_j_ano)

cat("Deployment firm-years:", nrow(deployment_firm_years), "obs,",
    n_distinct(deployment_firm_years$vat), "firms\n")

# Quick summary of emitter rates across draws
emit_counts <- sapply(proxy_list, function(px) sum(px$proxy > 0))
cat("  Emitters per draw: mean =", round(mean(emit_counts)),
    " sd =", round(sd(emit_counts)), "\n\n")


# =============================================================================
# STEP 8: Save
# =============================================================================

save(deployment_firm_years, file = file.path(PROC_DATA, "deployment_firm_years.RData"))
save(proxy_list, file = file.path(PROC_DATA, "deployment_proxy_list.RData"))
save(draw_params, file = file.path(PROC_DATA, "deployment_draw_params.RData"))

gpa_ok <- sapply(draw_params, function(x) !is.null(x$gpa_params))

cat("══════════════════════════════════════════════\n")
cat("Saved:\n")
cat("  deployment_firm_years.RData  —", nrow(deployment_firm_years), "firm-years\n")
cat("  deployment_proxy_list.RData  —", length(proxy_list), "draws\n")
cat("  deployment_draw_params.RData —", length(draw_params), "draws\n")
cat("    GPA fitted:", sum(gpa_ok), "/", length(draw_params), "draws\n")
cat("  B =", B_DRAWS, "| alpha =", ALPHA, "| seed =", BASE_SEED, "\n")
cat("══════════════════════════════════════════════\n")
