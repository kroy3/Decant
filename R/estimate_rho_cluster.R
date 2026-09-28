## estimate_rho_cluster.R
## Replacement for the v0 marker-ratio rho estimator, which took a median of
## observed/soup ratios over ALL soup-rich genes -- including genes the cell
## expresses natively -- and so overestimated rho 1.5-4x (worst at the low rho
## typical of scRNA-seq). See gate_rho().
##
## Identifiability argument (the same one SoupX rests on): for a gene a cluster
## does NOT express, every count it shows is ambient, so its count is
## Poisson(rho * T_k * soup_g). Those "absent" genes pin rho; every expressed
## gene can only sit ABOVE that line. So:
##   1. pool counts per cluster (per-cell counts on absent genes are too sparse);
##   2. start from the lowest decile of observed/expected ratios among
##      well-powered genes;
##   3. iterate: absent set A = genes NOT significantly above rho*expected
##      (one-sided Poisson test); rho = sum_A counts / sum_A expected;
##   4. per cell, shrink toward the cluster rho with a Gamma-Poisson empirical
##      Bayes posterior mean over the same absent set.
## When almost no gene is cleanly absent from a cluster, rho is genuinely not
## identifiable from counts; the estimator then biases UP, and says so via the
## per-cluster diagnostics (low_power) rather than hiding it.
##
## KNOWN FAILURE MODE (not fixed, detected): if the empty droplets contain
## material that is absent from the cells (debris, cell types filtered out of
## the cell set), those genes sit far BELOW rho*expected, the fixed point
## collapses toward 0, and rho is badly UNDERestimated (10% junk in empties took
## a true 0.2 to 0.005 in testing). A two-sided variant is robust to this but
## biased up on clean data, so it is used only as a tripwire: when the two
## disagree sharply across the DATASET (cell-weighted; per-cluster comparisons
## were too noisy and false-alarmed on clean data) the fit is flagged
## `soup_mismatch`.

#' Cluster-pooled, absent-gene rho estimator with per-cell shrinkage.
#'
#' @param observed genes x cells counts (base matrix or dgCMatrix).
#' @param soup length-G ambient profile (e.g. from [ambient_global()]).
#' @param clusters length-cells cluster labels. Prefer a moderate number of
#'   reasonably large clusters: very small clusters lose the power to reject
#'   weakly-expressed genes and bias rho upward.
#' @param alpha one-sided Poisson level for calling a gene "above soup"
#'   (i.e. natively expressed) in a cluster.
#' @param min_absent_soup a cluster whose absent set covers less than this
#'   fraction of soup mass is flagged `low_power`.
#' @param max_iter iteration cap for the absent-set fixed point.
#' @return numeric vector of per-cell rho, with attribute `"clusters"`: a
#'   data.frame of per-cluster diagnostics (n_cells, rho, n_absent_genes,
#'   absent_soup_frac, absent_counts, low_power, borrowed, rho_two_sided), and
#'   attribute `"soup_mismatch"` (logical). Low-power clusters take the
#'   evidence-weighted rho of the other clusters (`borrowed = TRUE`).
#'   `soup_mismatch = TRUE` means the soup profile looks inconsistent with the
#'   cells and rho is probably UNDERestimated -- inspect the empties or supply
#'   rho yourself.
#' @export
estimate_rho_cluster <- function(observed, soup, clusters, alpha = 0.01,
                                 min_absent_soup = 0.02, max_iter = 100) {
  observed <- .as_counts(observed, "observed")
  N <- ncol(observed)
  if (length(clusters) != N) stop("length(clusters) must equal ncol(observed)", call. = FALSE)
  if (length(soup) != nrow(observed)) stop("length(soup) must equal nrow(observed)", call. = FALSE)
  soup <- soup / sum(soup)
  T_c <- .col_sums(observed)
  ks <- sort(unique(clusters))
  Xk <- .cluster_sums(observed, clusters, ks)
  rho_c <- numeric(N)
  diag <- vector("list", length(ks))
  absent <- vector("list", length(ks))

  ## pass 1: cluster-level rho from the absent-gene fixed point
  for (i in seq_along(ks)) {
    idx <- which(clusters == ks[i])
    X <- Xk[, i]
    e <- soup * sum(T_c[idx])
    ok <- e > 0
    rho <- rho2 <- 0; A <- rep(FALSE, length(X))
    if (sum(X) > 0 && any(ok)) {
      r <- X / pmax(e, 1e-300)
      pw <- ok & e >= stats::quantile(e[ok], 0.5)
      lo <- pw & r <= stats::quantile(r[pw], 0.1)
      rho0 <- sum(X[lo]) / sum(e[lo])
      fit <- .absent_fixed_point(X, e, ok, rho0, alpha, max_iter, two_sided = FALSE)
      rho <- fit$rho; A <- fit$A
      rho2 <- .absent_fixed_point(X, e, ok, rho0, alpha, max_iter, two_sided = TRUE)$rho
    }
    absent[[i]] <- which(A)
    SA <- sum(soup[absent[[i]]])
    diag[[i]] <- data.frame(cluster = ks[i], n_cells = length(idx), rho = rho,
                            n_absent_genes = length(absent[[i]]), absent_soup_frac = SA,
                            absent_counts = sum(X[absent[[i]]]),
                            low_power = SA < min_absent_soup,
                            borrowed = FALSE,
                            rho_two_sided = rho2,
                            stringsAsFactors = FALSE)
  }
  diag <- do.call(rbind, diag)

  ## A cluster with (almost) no absent genes -- typically the fragile type that
  ## dominates the soup itself -- cannot identify its own rho. Rather than keep
  ## an arbitrary value, it borrows the evidence-weighted rho of the powered
  ## clusters as its prior (assumes rho does not depend strongly on cell type;
  ## flagged, not hidden). Its own per-cell evidence, if any, still counts.
  weak <- diag$low_power
  if (any(weak) && any(!weak)) {
    w <- diag$absent_counts[!weak]
    diag$rho[weak] <- if (sum(w) > 0) sum(diag$rho[!weak] * w) / sum(w) else mean(diag$rho[!weak])
    diag$borrowed[weak] <- TRUE
  }

  ## pass 2: per-cell empirical-Bayes shrinkage toward the cluster rho
  ## (Gamma prior, Poisson counts on that cluster's absent genes)
  for (i in seq_along(ks)) {
    idx <- which(clusters == ks[i])
    gA <- absent[[i]]; rho <- diag$rho[i]
    if (length(gA) > 0 && rho > 0) {
      y <- .col_sums(observed[gA, idx, drop = FALSE])
      n <- T_c[idx] * diag$absent_soup_frac[i]
      raw <- y / pmax(n, 1e-12)
      v_between <- if (length(idx) > 1) stats::var(raw) - mean(rho / pmax(n, 1e-12)) else 0
      v_between <- max(v_between, (0.05 * rho)^2, 1e-10)
      a <- rho^2 / v_between
      rho_c[idx] <- (a + y) / (a / rho + n)
    } else {
      rho_c[idx] <- rho
    }
  }
  out <- pmin(pmax(rho_c, 0), 0.99)
  names(out) <- colnames(observed)
  attr(out, "clusters") <- diag
  r1 <- stats::weighted.mean(diag$rho, diag$n_cells)
  r2 <- stats::weighted.mean(diag$rho_two_sided, diag$n_cells)
  attr(out, "soup_mismatch") <- r1 < 0.5 * r2 && r2 - r1 > 0.05
  out
}

## Absent-set fixed point: rho = sum_A X / sum_A e, where A = genes consistent
## with pure soup at the current rho (one-sided: not significantly ABOVE;
## two-sided: also not significantly BELOW).
.absent_fixed_point <- function(X, e, ok, rho, alpha, max_iter, two_sided) {
  A <- rep(FALSE, length(X))
  for (it in seq_len(max_iter)) {
    A <- ok & stats::ppois(X - 1, rho * e, lower.tail = FALSE) > alpha
    if (two_sided) A <- A & stats::ppois(X, rho * e) > alpha
    if (!any(A)) break
    new <- sum(X[A]) / sum(e[A])
    if (abs(new - rho) < 1e-8) { rho <- new; break }
    rho <- new
  }
  list(rho = rho, A = A)
}

#' Default clustering used when the caller supplies none.
#'
#' k-means on log-normalised top-variance genes. Sparse-safe (only an
#' n_hvg x cells block is densified). For real analyses pass your own
#' Seurat / SingleCellExperiment clusters instead.
#'
#' @param observed genes x cells counts.
#' @param k number of clusters.
#' @param n_hvg number of top-variance genes used.
#' @param seed RNG seed for k-means.
#' @return integer cluster labels.
#' @export
quick_labels <- function(observed, k = 6, n_hvg = 200, seed = 1) {
  observed <- .as_counts(observed, "observed")
  set.seed(seed)
  L <- .lognorm_hvg(observed, n_hvg)
  k <- min(k, ncol(observed) - 1L)
  stats::kmeans(t(L), centers = k, nstart = 5, iter.max = 50)$cluster
}
