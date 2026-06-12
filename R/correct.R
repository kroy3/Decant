## correct.R
## Subtract expected ambient counts with a HARD guarantee that the method can
## only ever remove or leave counts -- never add them. This is the property the
## April 2026 benchmark found scAR and CellClear violating (they fabricated
## counts). It is cheap to guarantee and worth advertising.

#' Remove ambient counts from each cell.
#'
#' corrected_gc = clamp( observed_gc - rho_c * T_c * soup_g , 0 , observed_gc )
#'
#' @param observed genes x cells count matrix.
#' @param soup length-G ambient probability vector.
#' @param rho length-cells contamination fractions (see estimate_rho()).
#' @return genes x cells corrected matrix. Guaranteed 0 <= corrected <= observed.
correct_counts <- function(observed, soup, rho) {
  T_c <- colSums(observed)
  ## expected ambient counts: outer(soup, rho * T_c)
  expected <- outer(soup, rho * T_c)
  corrected <- observed - expected
  corrected[corrected < 0] <- 0
  ## never remove more than was there (clamp upper bound) -- redundant given the
  ## floor, but makes the guarantee explicit and survives future edits.
  over <- corrected > observed
  corrected[over] <- observed[over]
  corrected
}

#' Estimate per-cell contamination fraction rho.
#'
#' v0 deliberately keeps rho estimation simple and SHARED across methods in the
#' benchmark, so that comparisons isolate the variable under test (the ambient
#' PROFILE), not rho-estimation noise. Here: for each cell, find genes the soup
#' expresses but the cell's own program does not, and attribute their counts to
#' soup. Replace with a SoupX-style marker estimator before any real use.
#'
#' @param observed genes x cells.
#' @param soup length-G ambient profile.
#' @param quantile_floor genes below this soup-quantile are ignored as estimators.
estimate_rho <- function(observed, soup, low_expr_q = 0.5) {
  T_c <- colSums(observed); T_c[T_c == 0] <- 1
  ## "soup-diagnostic" genes: clearly present in soup
  diag_genes <- which(soup >= quantile(soup[soup > 0], low_expr_q))
  obs_frac <- sweep(observed[diag_genes, , drop = FALSE], 2, T_c, "/")
  ## if a cell's observed fraction on these genes tracks the soup fraction,
  ## the ratio approximates rho. Use a robust median ratio, capped to [0,1).
  sp <- soup[diag_genes]
  ratio <- apply(obs_frac, 2, function(of) {
    use <- sp > 0
    median(of[use] / sp[use])
  })
  rho <- pmin(pmax(ratio, 0), 0.95)
  rho[is.na(rho)] <- 0
  rho
}
