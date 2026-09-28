## correct.R
## Remove expected ambient counts with a HARD guarantee that the method can only
## ever remove or leave counts -- never add them. This is the property the April
## 2026 benchmark found scAR and CellClear violating (they fabricated counts).
##
## All rules operate on non-zero entries only (a zero has nothing to remove), so
## they are sparse-safe and never build a dense genes x cells intermediate.
##
##  "subtract"     : v0 rule. corrected = clamp(x - rho*T*soup, 0, x). Every gene
##                   clamped at 0 silently drops the rest of its expected removal,
##                   so it systematically UNDER-removes.
##  "redistribute" : SoupX-style. Remove exactly rho*T counts per cell (or all of
##                   it, if less), re-spreading the shortfall from clamped genes
##                   over genes that still have counts, in soup proportion.
##  "bayes"        : posterior-expected native counts under the mixture
##                   x ~ (1-rho)*phi_k + rho*soup, with phi_k each cluster's
##                   soup-subtracted profile. Removes contamination where the
##                   soup explains the counts, not blindly by soup rank.
## Which rule is the default is decided by gate_correction(), not by taste.

#' Remove ambient counts from each cell (mass-conserving).
#'
#' @param observed genes x cells counts (base matrix or dgCMatrix).
#' @param soup length-G ambient probability vector (one soup for all cells;
#'   for per-sample soups use [decant()]).
#' @param rho length-cells contamination fractions.
#' @param method one of "subtract", "redistribute", "bayes" (see Details).
#' @param clusters cluster labels, required for method = "bayes".
#' @param max_iter iterations for "redistribute".
#' @param drop_zeros for sparse input, drop entries corrected to exactly 0.
#'   FALSE keeps the input's sparsity pattern (explicit zeros).
#' @return corrected counts, same class as `observed` (dense or dgCMatrix).
#'   Guaranteed 0 <= corrected <= observed element-wise.
#' @details See the header of `R/correct.R` for the three rules and
#'   [gate_correction()] for how they compare against ground truth.
#' @export
correct_counts <- function(observed, soup, rho,
                           method = c("subtract", "redistribute", "bayes"),
                           clusters = NULL, max_iter = 50, drop_zeros = TRUE) {
  method <- match.arg(method)
  observed <- .as_counts(observed, "observed")
  if (length(soup) != nrow(observed)) stop("length(soup) must equal nrow(observed)", call. = FALSE)
  if (length(rho) == 1L) rho <- rep(rho, ncol(observed))
  if (length(rho) != ncol(observed)) stop("length(rho) must equal ncol(observed)", call. = FALSE)
  soup <- soup / sum(soup)
  T_c <- .col_sums(observed)
  nz <- .nz(observed)

  v <- switch(method,
    subtract     = .corr_subtract(nz, soup, rho, T_c),
    redistribute = .corr_redistribute(nz, soup, rho, T_c, max_iter),
    bayes = {
      if (is.null(clusters)) stop("method = 'bayes' needs `clusters`", call. = FALSE)
      .corr_bayes(observed, nz, soup, rho, T_c, clusters)
    })
  ## the guarantee, enforced (not assumed) on every path
  v <- pmin(pmax(v, 0), nz$v)
  .nz_set(observed, nz, v, drop = drop_zeros)
}

.corr_subtract <- function(nz, soup, rho, T_c) {
  nz$v - (rho * T_c)[nz$j] * soup[nz$i]
}

.corr_redistribute <- function(nz, soup, rho, T_c, max_iter) {
  N <- length(T_c)
  target <- pmin(rho * T_c, T_c)
  rem <- numeric(length(nz$v))
  active <- nz$v > 0
  for (it in seq_len(max_iter)) {
    done <- .tabulate_sum(rem, nz$j, N)
    left <- target - done
    w <- ifelse(active, soup[nz$i], 0)
    wsum <- .tabulate_sum(w, nz$j, N)
    need <- left > 1e-8 & wsum > 0
    if (!any(need)) break
    add <- ifelse(need[nz$j], left[nz$j] * w / pmax(wsum[nz$j], 1e-300), 0)
    rem <- pmin(rem + add, nz$v)
    active <- rem < nz$v
  }
  nz$v - rem
}

.corr_bayes <- function(observed, nz, soup, rho, T_c, clusters) {
  ks <- sort(unique(clusters))
  Xk <- .cluster_sums(observed, clusters, ks)
  amb_k <- as.numeric(tapply(rho * T_c, factor(clusters, levels = ks), sum))
  phi <- pmax(Xk - outer(soup, amb_k), 0)
  phi <- sweep(phi, 2, pmax(colSums(phi), 1e-300), "/")
  kj <- match(clusters, ks)[nz$j]
  a <- (1 - rho[nz$j]) * phi[cbind(nz$i, kj)]
  b <- rho[nz$j] * soup[nz$i]
  nz$v * a / pmax(a + b, 1e-300)
}

## sum of x within groups 1..n (fast, no factor overhead)
.tabulate_sum <- function(x, g, n) {
  out <- numeric(n)
  s <- rowsum(x, g, reorder = FALSE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

#' Legacy v0 per-cell rho estimator (superseded).
#'
#' Median of observed/soup ratios over soup-rich genes. Because those genes
#' include ones the cell expresses natively, it overestimates rho 1.5-4x on
#' ground truth (see [gate_rho()]). Kept only so earlier benchmark numbers
#' remain reproducible; use [estimate_rho_cluster()] instead.
#'
#' @param observed genes x cells.
#' @param soup length-G ambient profile.
#' @param low_expr_q genes below this soup-quantile are ignored as estimators.
#' @return per-cell rho.
#' @export
estimate_rho <- function(observed, soup, low_expr_q = 0.5) {
  observed <- .as_counts(observed, "observed")
  T_c <- .col_sums(observed); T_c[T_c == 0] <- 1
  diag_genes <- which(soup >= stats::quantile(soup[soup > 0], low_expr_q))
  obs_frac <- as.matrix(observed[diag_genes, , drop = FALSE])
  obs_frac <- sweep(obs_frac, 2, T_c, "/")
  sp <- soup[diag_genes]
  ratio <- apply(obs_frac, 2, function(of) {
    use <- sp > 0
    stats::median(of[use] / sp[use])
  })
  rho <- pmin(pmax(ratio, 0), 0.95)
  rho[is.na(rho)] <- 0
  rho
}
