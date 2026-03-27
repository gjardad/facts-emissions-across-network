###############################################################################
# tests/test_pareto_redistribution.R
# Quick local test of the rewritten Pareto redistribution with combined ranking.
# Runs B=2 draws, 2 years, sequential. Checks no errors and prints diagnostics.
###############################################################################

REPO_DIR <- "c:/Users/jota_/Documents/facts-emissions-across-network"
source(file.path(REPO_DIR, "paths.R"))
source(file.path(REPO_DIR, "utils", "calibration_pareto.R"))
source(file.path(REPO_DIR, "utils", "sector_conventions.R"))

library(dplyr)
library(lmom)

# ── Load data (same as b_loop_scope1_dispersion_pareto.R Sections 1-3) ──────
YEARS <- 2010:2011  # just 2 years for testing

load(file.path(PROC_DATA, "deployment_proxy_list.RData"))
B <- min(2L, length(proxy_list))
proxy_list <- proxy_list[1:B]
cat("B =", B, "draws,", length(YEARS), "years\n")

load(file.path(PROC_DATA, "firm_year_belgian_euets.RData"))
eutl <- firm_year_belgian_euets %>%
  filter(year %in% YEARS) %>%
  select(vat, year, emissions)
rm(firm_year_belgian_euets)

load(file.path(PROC_DATA, "annual_accounts_selected_sample_key_variables.RData"))
accounts <- df_annual_accounts_selected_sample_key_variables %>%
  filter(year %in% YEARS) %>%
  select(vat, year, nace5d, revenue) %>%
  mutate(nace2d  = make_nace2d(nace5d),
         revenue = pmax(coalesce(revenue, 0), 0))
rm(df_annual_accounts_selected_sample_key_variables)

load(file.path(PROC_DATA, "deployment_panel.RData"))

nace_crf <- read.csv(
  file.path(REPO_DIR, "preprocess", "crosswalks", "nace_crf_crosswalk.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(nace2d = "character")
) %>%
  select(nace2d, crf_group)

deploy_nace <- deployment_panel %>%
  mutate(nace2d = make_nace2d(nace5d)) %>%
  select(vat, nace2d, nace5d) %>%
  distinct(vat, .keep_all = TRUE) %>%
  left_join(nace_crf, by = "nace2d")
rm(deployment_panel)

# NIR targets
load(file.path(PROC_DATA, "nir_calibration_targets.RData"))

E_ETS_group <- eutl %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  left_join(nace_crf, by = "nace2d") %>%
  filter(!is.na(crf_group), year %in% YEARS) %>%
  group_by(crf_group, year) %>%
  summarise(E_ETS = sum(emissions, na.rm = TRUE), .groups = "drop")

E_deploy_panel <- nir_targets %>%
  filter(year %in% YEARS) %>%
  mutate(E_NIR = E_NIR_kt * 1000) %>%
  left_join(E_ETS_group, by = c("crf_group", "year")) %>%
  mutate(E_ETS    = coalesce(E_ETS, 0),
         E_deploy = pmax(E_NIR - E_ETS, 0))

# GPA — use global fallback for test
eutl_with_nace <- eutl %>%
  filter(emissions > 0) %>%
  left_join(accounts %>% distinct(vat, year, nace2d), by = c("vat", "year")) %>%
  filter(!is.na(nace2d))
ref_dist <- build_reference_dist(eutl_with_nace$emissions,
                                  eutl_with_nace$year,
                                  eutl_with_nace$nace2d)
gpa_global <- fit_gpa(ref_dist)
cat("GPA fitted:", !is.null(gpa_global), "\n")

FALLBACK_THRESHOLD <- 0.35
MIN_N_STATS <- 3L

# ── Simplified sequential run (no parallel) ──────────────────────────────────
all_flags <- list()

for (t in YEARS) {
  cat(sprintf("\n=== Year %d ===\n", t))

  eutl_t     <- eutl[eutl$year == t, ]
  accounts_t <- accounts[accounts$year == t, ]
  E_dep_t    <- E_deploy_panel[E_deploy_panel$year == t, c("crf_group", "E_deploy")]

  proxy_t <- lapply(proxy_list, function(px) {
    sub <- px[px$year == t, c("vat", "proxy")]
    setNames(sub$proxy, sub$vat)
  })

  ok_ets <- !is.na(eutl_t$emissions) & eutl_t$emissions > 0
  ets_vats_t      <- eutl_t$vat[ok_ets]
  ets_emissions_t <- setNames(eutl_t$emissions[ok_ets], eutl_t$vat[ok_ets])

  # Firm-level ETS with CRF group
  ets_crf_t <- eutl_t[ok_ets, ] %>%
    left_join(accounts_t %>% distinct(vat, nace2d), by = "vat") %>%
    left_join(nace_crf, by = "nace2d") %>%
    filter(!is.na(crf_group)) %>%
    select(vat, emissions, crf_group)
  cat("  ETS firms with CRF:", nrow(ets_crf_t), "\n")

  for (b in 1:B) {
    cat(sprintf("  Draw %d:\n", b))
    proxy_b     <- proxy_t[[b]]
    gpa_b       <- gpa_global
    threshold_b <- FALLBACK_THRESHOLD

    deploy_b <- data.frame(
      vat       = names(proxy_b),
      proxy_val = as.numeric(proxy_b),
      stringsAsFactors = FALSE
    ) %>%
      filter(proxy_val > 0) %>%
      left_join(deploy_nace, by = "vat") %>%
      filter(!is.na(crf_group)) %>%
      left_join(E_dep_t, by = "crf_group") %>%
      mutate(E_deploy = coalesce(E_deploy, 0))

    all_crf_groups <- union(unique(deploy_b$crf_group),
                             unique(ets_crf_t$crf_group))

    n_mixed <- 0; n_pure_d <- 0; n_pure_e <- 0
    n_lowered <- 0; n_capped <- 0

    for (cg in all_crf_groups) {
      cg_deploy <- deploy_b[deploy_b$crf_group == cg, ]
      cg_ets    <- ets_crf_t[ets_crf_t$crf_group == cg, ]
      n_ets     <- nrow(cg_ets)
      E_cg      <- if (nrow(cg_deploy) > 0) cg_deploy$E_deploy[1] else {
        ed <- E_dep_t$E_deploy[E_dep_t$crf_group == cg]
        if (length(ed) > 0) ed[1] else 0
      }
      min_ets_emit <- if (n_ets > 0) min(cg_ets$emissions) else Inf
      is_mixed <- n_ets > 0

      if (nrow(cg_deploy) == 0) { n_pure_e <- n_pure_e + 1; next }
      if (E_cg <= 0) next

      proxy_levels_full <- pmax(sinh(cg_deploy$proxy_val), 0)
      thr_current       <- threshold_b
      flag              <- if (is_mixed) "mixed_ok" else "pure_deploy"

      for (iter in 1:20) {
        proxy_levels <- proxy_levels_full
        pos_mask <- proxy_levels > 0
        if (sum(pos_mask) >= 2 && thr_current > 0) {
          tau <- quantile(proxy_levels[pos_mask], probs = thr_current, names = FALSE)
          proxy_levels[proxy_levels <= tau] <- 0
        }

        emit_mask     <- proxy_levels > 0
        n_deploy_emit <- sum(emit_mask)

        if (n_deploy_emit == 0) {
          if (thr_current > 1e-4) { thr_current <- thr_current / 2; next }
          break
        }

        n_total    <- n_ets + n_deploy_emit
        emit_idx   <- which(emit_mask)
        emit_proxy <- proxy_levels[emit_idx]
        ranks      <- rank(emit_proxy, ties.method = "average")

        pw <- pareto_weights(ranks, n_total, gpa_b)

        if (is.null(pw)) {
          flag <- "gpa_fallback_sinh"
          total_sinh <- sum(sinh(cg_deploy$proxy_val[emit_idx]))
          if (total_sinh > 0) {
            emissions_b <- E_cg * sinh(cg_deploy$proxy_val[emit_idx]) / total_sinh
          } else {
            emissions_b <- rep(E_cg / n_deploy_emit, n_deploy_emit)
          }
        } else {
          emissions_b <- E_cg * pw / sum(pw)
        }

        if (is_mixed && max(emissions_b) >= min_ets_emit) {
          if (thr_current > 1e-4) {
            thr_current <- thr_current / 2
            flag <- "mixed_lowered"
            next
          }
          flag <- "mixed_capped"
          cap  <- min_ets_emit * (1 - 1e-6)
          over <- which(emissions_b > cap)
          freed <- sum(emissions_b[over] - cap)
          emissions_b[over] <- cap
          under <- which(emissions_b <= cap)
          if (length(under) > 0) {
            emissions_b[under] <- emissions_b[under] +
              freed * emissions_b[under] / sum(emissions_b[under])
          }
        }
        break
      }

      if (flag == "mixed_ok") n_mixed <- n_mixed + 1
      if (flag == "pure_deploy") n_pure_d <- n_pure_d + 1
      if (flag == "mixed_lowered") n_lowered <- n_lowered + 1
      if (flag == "mixed_capped") n_capped <- n_capped + 1

      all_flags[[length(all_flags) + 1]] <- data.frame(
        crf_group = cg, alloc_flag = flag, year = t, draw = b,
        n_ets = n_ets, n_deploy_emit = n_deploy_emit,
        max_deploy_emit = if (n_deploy_emit > 0) max(emissions_b) else NA,
        min_ets_emit = if (is.finite(min_ets_emit)) min_ets_emit else NA,
        stringsAsFactors = FALSE
      )
    }

    cat(sprintf("    pure_ets=%d  pure_deploy=%d  mixed_ok=%d  lowered=%d  capped=%d\n",
                n_pure_e, n_pure_d, n_mixed, n_lowered, n_capped))
  }
}

flags_df <- bind_rows(all_flags)
cat("\n=== FLAG SUMMARY ===\n")
print(table(flags_df$alloc_flag))

# Check constraint: for mixed sectors, max deploy emission < min ETS emission
mixed <- flags_df[grepl("mixed", flags_df$alloc_flag) &
                    !is.na(flags_df$max_deploy_emit) &
                    !is.na(flags_df$min_ets_emit), ]
if (nrow(mixed) > 0) {
  violations <- mixed[mixed$max_deploy_emit >= mixed$min_ets_emit, ]
  cat("\nConstraint violations (max_deploy >= min_ets):", nrow(violations),
      "out of", nrow(mixed), "mixed-sector cells\n")
  if (nrow(violations) > 0) print(violations)
} else {
  cat("\nNo mixed-sector cells to check constraint\n")
}

cat("\n=== SUCCESS ===\n")
