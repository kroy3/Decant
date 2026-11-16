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
##  "bayes_cell"   : "bayes" made cell-aware. "bayes" scales every cell of a
##                   cluster by the same ratio, so a cell whose count of a gene
##                   is far too high to be ambient (a doublet, a rare expresser)
##                   loses most of it, and it inflates the cluster's native
##                   profile so that the other cells keep ambient counts. Here
##                   an entry above the upper Poisson bound of its own ambient
##                   expectation loses only that expectation, and the cluster
##                   profile is built from the remaining entries.
## Which rule is the default is decided by gate_correction(), not by taste.

#' Remove ambient counts from each cell (mass-conserving).
#'
#' @param observed genes x cells counts (base matrix or dgCMatrix).
#' @param soup length-G ambient probability vector (one soup for all cells;
#'   for per-sample soups use [decant()]).
#' @param rho length-cells contamination fractions.
#' @param method one of "subtract", "redistribute", "bayes", "bayes_cell"
#'   (see Details).
#' @param clusters cluster labels, required for the two bayes methods.
#' @param max_iter iterations for "redistribute".
#' @param cap_q for "bayes_cell": an entry is treated as not explainable by
#'   ambient RNA when its count exceeds the `1 - cap_q` Poisson quantile of
#'   its ambient expectation.
#' @param cap_test for "bayes_cell": how the "not explainable by ambient"
#'   test is applied. "entry" tests each entry on its own; "cluster" corrects
#'   it for the number of cells of its cluster it is applied to (Bonferroni):
#'   an entry is set aside only if a count that high would be unexpected
#'   among all those cells. With "entry", any single count of a gene whose
#'   expected ambient count is below `cap_q` is set aside and kept, although
#'   thousands of cells make some such counts expected; on the 20k species
#'   mixture this kept most ambient counts of rarely contaminating genes.
#'   "cluster" removes more but also keeps less of the cells' own rare
#'   counts: on held-out mixtures removal 0.519 vs 0.444, retain 0.974 vs
#'   0.986 on the 10k mixture (`validation/BENCHMARK.md`, section 1e), so
#'   it is not the default.
#' @param drop_zeros for sparse input, drop entries corrected to exactly 0.
#'   FALSE keeps the input's sparsity pattern (explicit zeros).
#' @return corrected counts, same class as `observed` (dense or dgCMatrix).
#'   Guaranteed 0 <= corrected <= observed element-wise.
#' @details See the header of `R/correct.R` for the rules and
#'   [gate_correction()] for how they compare against ground truth.
#' @export
correct_counts <- function(observed, soup, rho,
                           method = c("subtract", "redistribute", "bayes", "bayes_cell"),
                           clusters = NULL, max_iter = 50, drop_zeros = TRUE, cap_q = 1e-3,
                           cap_test = c("entry", "cluster")) {
  method <- match.arg(method); cap_test <- match.arg(cap_test)
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
    },
    bayes_cell = {
      if (is.null(clusters)) stop("method = 'bayes_cell' needs `clusters`", call. = FALSE)
      .corr_bayes_cell(observed, nz, soup, rho, T_c, clusters, cap_q, cap_test)
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

## Cell-aware bayes. lambda = rho_c * T_c * soup_g is an entry's expected ambient
## count. An entry above the (1 - cap_q) Poisson quantile of lambda cannot be
## ambient alone: it keeps everything but lambda. The cluster's native profile
## is then estimated from the other entries only, against the ambient expected
## in those same cells, and applied to them as in "bayes".
.corr_bayes_cell <- function(observed, nz, soup, rho, T_c, clusters, cap_q, cap_test = "entry") {
  ks <- sort(unique(clusters)); K <- length(ks); G <- nrow(observed)
  kj <- match(clusters, ks)[nz$j]
  amb_c <- rho * T_c
  lam <- amb_c[nz$j] * soup[nz$i]
  high <- if (cap_test == "entry") {
    nz$v > stats::qpois(1 - cap_q, lam)
  } else {
    ## P(count >= observed | ambient only), times the cells it is tested in
    n_k <- tabulate(match(clusters, ks), K)
    stats::ppois(nz$v - 1, lam, lower.tail = FALSE) * n_k[kj] < cap_q
  }
  gk <- function(w) {                       # sum w over entries, by (gene, cluster)
    m <- matrix(0, G, K)
    s <- rowsum(w, (kj - 1L) * G + nz$i, reorder = FALSE)
    m[as.integer(rownames(s))] <- s[, 1]
    m
  }
  ## counts and expected ambient among the entries that ambient could explain.
  ## Zero entries are always in that set, so start from the cluster totals.
  amb_k <- as.numeric(tapply(amb_c, factor(clusters, levels = ks), sum))
  X_low <- gk(ifelse(high, 0, nz$v))
  E_low <- outer(soup, amb_k) - soup * gk(ifelse(high, amb_c[nz$j], 0))
  native <- pmax(X_low - E_low, 0)
  ## native rate per native UMI, among the same cells: their native depth is
  ## the cluster's minus that of the cells set aside for this gene
  nat_c <- (1 - rho) * T_c
  nat_k <- as.numeric(tapply(nat_c, factor(clusters, levels = ks), sum))
  D_low <- pmax(matrix(nat_k, G, K, byrow = TRUE) - gk(ifelse(high, nat_c[nz$j], 0)), 1e-300)
  phi <- native / D_low
  a <- (1 - rho[nz$j]) * phi[cbind(nz$i, kj)]
  b <- rho[nz$j] * soup[nz$i]
  ifelse(high, nz$v - pmin(lam, nz$v), nz$v * a / pmax(a + b, 1e-300))
}
