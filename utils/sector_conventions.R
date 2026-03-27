###############################################################################
# utils/sector_conventions.R
#
# Sector recoding conventions used project-wide.
#
# RULE: NACE 17 (paper) and 18 (printing) are ALWAYS treated as a single
#       sector "17/18". See THIS_PROJECT.md § Sector Conventions.
###############################################################################

#' Derive nace2d from nace5d with project conventions applied.
#' NACE 17 and 18 are recoded to "17/18".
#'
#' @param nace5d Character vector of 5-digit NACE codes.
#' @return Character vector of 2-digit NACE codes (with "17/18" convention).
make_nace2d <- function(nace5d) {
  n2 <- substr(nace5d, 1, 2)
  n2[n2 %in% c("17", "18")] <- "17/18"
  n2
}
