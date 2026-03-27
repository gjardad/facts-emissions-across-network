###############################################################################
# tests/test_subsampled_en.R
# Quick local test of run_one_subsample with B=2 draws on downsampled data.
# NOT a unit test — just checks the code runs without errors.
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "calibration_pareto.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)
library(Matrix)
library(glmnet)

B_DRAWS        <- 2L
SUBSAMPLE_FRAC <- 0.8
K_INNER        <- 5L
MIN_LHS_BUYERS <- 5L
ALPHA          <- 0.5
BASE_SEED      <- 2026L

cat("Loading data...\n")
load(file.path(PROC_DATA, "b2b_selected_sample.RData"))
b2b <- df_b2b_selected_sample
rm(df_b2b_selected_sample)

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
cat("Training LHS:", nrow(lhs), "firm-years\n")
cat("Sectors:", paste(sort(unique(lhs$nace2d)), collapse = ", "), "\n")

# Assign primary sector
firm_sector <- lhs %>%
  group_by(vat) %>%
  summarise(primary_nace2d = names(which.max(table(nace2d))), .groups = "drop")
lhs <- lhs %>% left_join(firm_sector, by = "vat")
training_vats <- unique(lhs$vat)

# Split B2B
b2b_filtered <- b2b %>% filter(year >= 2005)
deploy_vats <- setdiff(unique(b2b_filtered$vat_j_ano), training_vats)
cat("Training buyers:", length(training_vats),
    "| Deployment buyers:", length(deploy_vats), "\n")

b2b_lhs    <- b2b_filtered[b2b_filtered$vat_j_ano %in% training_vats, ]
b2b_deploy <- b2b_filtered[b2b_filtered$vat_j_ano %in% deploy_vats, ]
rm(b2b, b2b_filtered)

MIXED_SECTORS <- c("17/18", "19", "24")

# ── Helper functions (copied from run_subsampled_en.R) ──────────────────────

extract_suppliers <- function(fit, n_ctrl, eligible_sellers, s = "lambda.min") {
  co           <- coef(fit, s = s)
  supplier_idx <- (n_ctrl + 2):length(co)
  vals         <- co[supplier_idx]
  data.frame(
    vat_i_ano = eligible_sellers,
    coef      = as.numeric(vals),
    stringsAsFactors = FALSE
  ) %>%
    filter(coef != 0) %>%
    arrange(desc(abs(coef)))
}

youden_threshold_pct <- function(proxy_vals, is_emit, n_grid = 200) {
  pos_vals <- proxy_vals[proxy_vals > 0]
  if (sum(is_emit) < 3 || sum(!is_emit) < 3 || length(pos_vals) < 5) return(0)
  pct_grid <- seq(0, 0.95, length.out = n_grid)
  best_j   <- -Inf
  best_pct <- 0
  for (pct in pct_grid) {
    tau <- quantile(pos_vals, probs = pct, names = FALSE)
    pred_pos <- (proxy_vals > tau)
    tp  <- sum(pred_pos & is_emit)
    fp  <- sum(pred_pos & !is_emit)
    fn  <- sum(!pred_pos & is_emit)
    tn  <- sum(!pred_pos & !is_emit)
    tpr <- tp / max(tp + fn, 1)
    fpr <- fp / max(fp + tn, 1)
    j   <- tpr - fpr
    if (j > best_j) { best_j <- j; best_pct <- pct }
  }
  best_pct
}

draw_subsample_vats <- function(lhs, frac, seed) {
  set.seed(seed)
  lhs %>%
    distinct(vat, primary_nace2d, emit) %>%
    group_by(primary_nace2d, emit) %>%
    slice_sample(prop = frac) %>%
    ungroup() %>%
    pull(vat)
}

# ── run_one_subsample (from run_subsampled_en.R, with diagnostics) ──────────

run_one_subsample <- function(lhs, b2b_lhs, b2b_deploy,
                               alpha, K_inner, min_buyers,
                               subsample_frac, seed) {
  subsample_vats <- draw_subsample_vats(lhs, subsample_frac, seed)
  train_lhs      <- lhs[lhs$vat %in% subsample_vats, ]
  train_lhs$row_idx <- seq_len(nrow(train_lhs))

  b2b_train    <- b2b_lhs[b2b_lhs$vat_j_ano %in% train_lhs$vat, ]
  seller_counts <- b2b_train %>%
    distinct(vat_i_ano, vat_j_ano) %>%
    count(vat_i_ano, name = "n_lhs_buyers")
  eligible_sellers <- seller_counts %>%
    filter(n_lhs_buyers >= min_buyers) %>%
    pull(vat_i_ano) %>%
    sort()

  cat("  Eligible sellers:", length(eligible_sellers), "\n")

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

  year_dummies <- model.matrix(~ factor(year), data = train_lhs)[, -1, drop = FALSE]
  colnames(year_dummies) <- paste0("yr_", sort(unique(train_lhs$year))[-1])
  train_sectors  <- sort(unique(train_lhs$nace2d))
  sector_dummies <- model.matrix(~ factor(nace2d, levels = train_sectors),
                                  data = train_lhs)[, -1, drop = FALSE]
  colnames(sector_dummies) <- paste0("sec_", train_sectors[-1])
  X_controls <- cbind(log_revenue = train_lhs$log_revenue, year_dummies, sector_dummies)
  n_controls <- ncol(X_controls)
  X_full <- cbind(Matrix(X_controls, sparse = TRUE), X_asinh)
  pf     <- c(rep(0, n_controls), rep(1, n_cols))
  rm(X_asinh, X_controls, year_dummies, sector_dummies)

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

  set.seed(seed + 1000L)
  firms_b        <- unique(train_lhs$vat)
  inner_fold_map <- sample(rep(1:K_inner, length.out = length(firms_b)))
  names(inner_fold_map) <- firms_b
  inner_foldid <- unname(inner_fold_map[train_lhs$vat])

  y_vec <- asinh(train_lhs$y)
  cat("  Fitting cv.glmnet (", n_rows, "obs x", ncol(X_full), "cols)...\n")
  fit <- cv.glmnet(
    x = X_full, y = y_vec,
    family = "gaussian", alpha = alpha,
    penalty.factor = pf, foldid = inner_foldid,
    standardize = TRUE
  )
  cat("  lambda.min =", fit$lambda.min, "\n")

  coef_lookup <- extract_suppliers(fit, n_controls, eligible_sellers_enet, "lambda.min")
  cat("  Non-zero suppliers:", nrow(coef_lookup), "\n")

  if (nrow(coef_lookup) == 0) {
    rm(fit, X_full)
    return(list(
      proxy = data.frame(vat = character(0), year = integer(0),
                          proxy = numeric(0), stringsAsFactors = FALSE),
      gpa_params = NULL,
      threshold_pct = NA_real_
    ))
  }

  # Stage B: deployment proxy
  proxy <- b2b_deploy %>%
    inner_join(coef_lookup, by = "vat_i_ano") %>%
    group_by(vat_j_ano, year) %>%
    summarise(proxy = sum(coef * asinh(corr_sales_ij), na.rm = TRUE), .groups = "drop") %>%
    rename(vat = vat_j_ano) %>%
    mutate(proxy = pmax(proxy, 0))
  cat("  Deployment proxy rows:", nrow(proxy), "\n")

  # Stage A: Inner-CV threshold
  mixed_idx <- which(train_lhs$nace2d %in% MIXED_SECTORS)
  threshold_pct <- NA_real_
  cat("  Mixed-sector firm-years:", length(mixed_idx), "\n")

  if (length(mixed_idx) > 0) {
    mixed_df <- data.frame(
      vat   = train_lhs$vat[mixed_idx],
      year  = train_lhs$year[mixed_idx],
      emit  = train_lhs$y[mixed_idx] > 0,
      sec   = train_lhs$nace2d[mixed_idx],
      fold  = inner_foldid[mixed_idx],
      proxy = 0,
      stringsAsFactors = FALSE
    )

    for (k in seq_len(K_inner)) {
      rows_k <- which(mixed_df$fold == k)
      if (length(rows_k) == 0) next
      train_rows <- which(inner_foldid != k)

      fit_k <- glmnet(
        x = X_full[train_rows, , drop = FALSE],
        y = y_vec[train_rows],
        family = "gaussian", alpha = alpha,
        penalty.factor = pf, lambda = fit$lambda.min,
        standardize = TRUE
      )

      coef_k <- extract_suppliers(fit_k, n_controls, eligible_sellers_enet,
                                   s = fit$lambda.min)
      cat("    Fold", k, ":", nrow(coef_k), "suppliers,",
          length(rows_k), "mixed held-out\n")
      if (nrow(coef_k) == 0) next

      fold_vats <- mixed_df$vat[rows_k]
      proxy_k <- b2b_lhs[b2b_lhs$vat_j_ano %in% fold_vats, ] %>%
        inner_join(coef_k, by = "vat_i_ano") %>%
        group_by(vat_j_ano, year) %>%
        summarise(proxy = sum(coef * asinh(corr_sales_ij), na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(proxy = pmax(proxy, 0))

      key_df    <- paste(mixed_df$vat[rows_k], mixed_df$year[rows_k])
      key_proxy <- paste(proxy_k$vat_j_ano, proxy_k$year)
      m <- match(key_df, key_proxy)
      mixed_df$proxy[rows_k] <- ifelse(is.na(m), 0, proxy_k$proxy[m])
    }

    cat("  Mixed-sector proxy > 0:", sum(mixed_df$proxy > 0),
        "/", nrow(mixed_df), "\n")

    pcts <- numeric(0)
    for (sec in unique(mixed_df$sec)) {
      sec_mask <- mixed_df$sec == sec
      pct_s <- youden_threshold_pct(mixed_df$proxy[sec_mask], mixed_df$emit[sec_mask])
      cat("    Sector", sec, ": Youden pct =", pct_s, "\n")
      if (!is.na(pct_s) && pct_s > 0) pcts <- c(pcts, pct_s)
    }
    if (length(pcts) > 0) threshold_pct <- mean(pcts)
  }
  rm(fit, X_full)

  ets_emit <- train_lhs[train_lhs$euets == 1 & train_lhs$y > 0, ]
  ref_dist   <- build_reference_dist(ets_emit$y, ets_emit$year, ets_emit$nace2d)
  gpa_params <- fit_gpa(ref_dist)

  cat("  threshold_pct =", threshold_pct,
      "| GPA fitted:", !is.null(gpa_params), "\n")
  list(proxy = proxy, gpa_params = gpa_params, threshold_pct = threshold_pct)
}

# ── Run 2 draws ─────────────────────────────────────────────────────────────
cat("\n=== Draw 1 ===\n")
res1 <- run_one_subsample(lhs, b2b_lhs, b2b_deploy,
                           ALPHA, K_INNER, MIN_LHS_BUYERS,
                           SUBSAMPLE_FRAC, BASE_SEED + 1)

cat("\n=== Draw 2 ===\n")
res2 <- run_one_subsample(lhs, b2b_lhs, b2b_deploy,
                           ALPHA, K_INNER, MIN_LHS_BUYERS,
                           SUBSAMPLE_FRAC, BASE_SEED + 2)

cat("\n=== SUCCESS ===\n")
cat("Draw 1: proxy rows =", nrow(res1$proxy),
    "| threshold =", res1$threshold_pct, "\n")
cat("Draw 2: proxy rows =", nrow(res2$proxy),
    "| threshold =", res2$threshold_pct, "\n")
