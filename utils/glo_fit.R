###############################################################################
# utils/glo_fit.R
#
# Self-contained GLO (Generalized Logistic) helpers used by the deployment
# allocation pipeline.
#
#   fit_glo(d_vec)        : fit GLO by L-moments to a numeric vector
#   glo_weights(ranks,    : evaluate quaglo at plotting positions and return
#               n_total,    un-normalized exp(w - max(w)) weights for the
#               glo_par)    given ranks within a cell of size n_total
#
# References:
#   inferring_emissions/analysis/active/threshold_p_levels_eval_all.R:84
###############################################################################

suppressPackageStartupMessages(library(lmom))

#' Fit a GLO by method of L-moments.
#'
#' @param d_vec Numeric vector (typically doubly-demeaned log emissions).
#' @return Named numeric vector of GLO parameters (xi, alpha, k), or NULL if
#'   the fit fails.
fit_glo <- function(d_vec) {
  d_vec <- d_vec[is.finite(d_vec)]
  if (length(d_vec) < 4) return(NULL)
  lmoms <- lmom::samlmu(sort(d_vec), nmom = 4)
  par <- tryCatch(lmom::pelglo(lmoms[1:3]), error = function(e) NULL)
  if (is.null(par)) return(NULL)
  attr(par, "lmoms") <- lmoms
  attr(par, "n")     <- length(d_vec)
  par
}

#' Compute un-normalized GLO rank weights for a cell.
#'
#' For ranks 1..n_total in a cell of size n_total, compute plotting positions
#' p_k = (rank_k - 0.5) / n_total, evaluate w_k = quaglo(p_k, glo_par), and
#' return exp(w_k - max(w)). The caller is responsible for normalizing across
#' whatever subset of ranks the shares should sum to 1 over.
#'
#' Rank 1 receives the highest weight (highest GLO quantile).
#'
#' @param ranks    Integer vector of ranks (in 1..n_total).
#' @param n_total  Total number of firms in the cell.
#' @param glo_par  GLO parameter vector returned by fit_glo().
#' @return Numeric vector of un-normalized weights, same length as `ranks`.
glo_weights <- function(ranks, n_total, glo_par) {
  if (n_total < 1L) return(numeric(0))
  # Convert "rank 1 = highest" to a plotting position near 1.
  p_k <- (n_total - ranks + 0.5) / n_total
  w_k <- lmom::quaglo(p_k, glo_par)
  if (any(!is.finite(w_k))) return(rep(NA_real_, length(ranks)))
  exp(w_k - max(w_k))
}
