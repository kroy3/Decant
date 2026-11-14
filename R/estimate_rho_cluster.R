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
##   3. iterate: absent set A = well-powered genes (expected ambient count
##      >= min_expected) NOT significantly above rho*expected (one-sided
##      Poisson test); rho = sum_A counts / sum_A expected;
##   4. per cell, shrink toward the cluster rho with a Gamma-Poisson empirical
##      Bayes posterior mean over the same absent set.
## When almost no gene is cleanly absent from a cluster, rho is genuinely not
## identifiable from counts; the estimator then biases UP, and says so via the
## per-cluster diagnostics (low_power) rather than hiding it.
##
## KNOWN FAILURE MODE: empties that do not match what the cells were bathed in.
##  * Genes in the empties that no cell contains at all are detected and
##    dropped from the soup (see `floor_rho`); a large share raises
##    `soup_mismatch`.
##  * Genes the cells do contain but the empties over-represent (debris of a
##    cell type that is under-sampled among the cells) are NOT detected. They
##    sit below rho*expected, the fixed point settles on them, and rho is
##    UNDERestimated (10% such junk took a true 0.2 to 0.005 in simulation).
##    A two-sided variant of the fit used to serve as a tripwire for this. On
##    real PBMC data it raised a false alarm (it runs away in the clusters
##    that dominate the soup), and with the coverage requirement it no longer
##    separates the junk case from a true rho of zero. It is still reported
##    per cluster (rho_two_sided) but no longer raises a warning.

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
#' @param min_expected a gene can join the absent set only if its expected
#'   ambient count in the cluster, at the current rho, is at least this. Below
#'   it the Poisson test cannot tell a weakly expressed gene from an absent
#'   one, and on real data (tens of thousands of low-count genes) letting
#'   those genes in drove rho to several times its true value.
#' @param floor_rho the smallest contamination considered possible. Genes
#'   whose total count in the cells is significantly below this fraction of
#'   their soup expectation are treated as soup-only and dropped from the soup.
#' @return numeric vector of per-cell rho, with attribute `"clusters"`: a
#'   data.frame of per-cluster diagnostics (n_cells, rho, n_absent_genes,
#'   absent_soup_frac, absent_counts, low_power, borrowed, rho_two_sided), and
#'   attributes `"soup_mismatch"` (logical), `"soup_only"` (logical per gene)
#'   and `"soup_only_frac"` (their share of the soup). Low-power clusters take the
#'   evidence-weighted rho of the other clusters (`borrowed = TRUE`).
#'   `soup_mismatch = TRUE` means more than 5% of the soup was in genes no
#'   cell contains; those genes are dropped, but the empties deserve a look.
#' @export
estimate_rho_cluster <- function(observed, soup, clusters, alpha = 0.01,
                                 min_absent_soup = 0.02, max_iter = 100, min_expected = 5,
                                 floor_rho = 0.005) {
  observed <- .as_counts(observed, "observed")
  N <- ncol(observed)
  if (length(clusters) != N) stop("length(clusters) must equal ncol(observed)", call. = FALSE)
  if (length(soup) != nrow(observed)) stop("length(soup) must equal nrow(observed)", call. = FALSE)
  soup <- soup / sum(soup)
  T_c <- .col_sums(observed)
  ks <- sort(unique(clusters))
  Xk <- .cluster_sums(observed, clusters, ks)
  ## Soup-only genes. Ambient RNA lands in cells too, so a gene at soup
  ## fraction s must show at least rho * s in the cells as a whole. A gene
  ## whose total count is significantly below even 0.5% contamination is in
  ## the empties but not in what the cells were bathed in (debris, cell types
  ## lost before capture). Left in, such genes sit far below the line and drag
  ## the absent-gene fixed point toward zero, so they are removed from the
  ## soup and the rest is renormalised.
  e_all <- soup * sum(T_c)
  soup_only <- floor_rho * e_all >= min_expected &
    stats::ppois(rowSums(Xk), floor_rho * e_all) < alpha
  soup_only_frac <- sum(soup[soup_only])
  if (any(soup_only)) { soup[soup_only] <- 0; soup <- soup / sum(soup) }
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
      ## start among well-covered genes when there are enough of them (real
      ## data: most genes have almost no soup, and the upper half by soup is
      ## still mostly zeros); otherwise the upper half by soup
      pw <- ok & e >= 4 * min_expected
      if (sum(pw) < 50) pw <- ok & e >= stats::quantile(e[ok], 0.5)
      lo <- pw & r <= stats::quantile(r[pw], 0.1)
      rho0 <- sum(X[lo]) / sum(e[lo])
      fit <- .absent_fixed_point(X, e, ok, rho0, alpha, max_iter, two_sided = FALSE, min_expected)
      rho <- fit$rho; A <- fit$A
      ## two-sided variant, reported as a diagnostic only
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
  ## an arbitrary value, it borrows the evidence-weighted median rho of the
  ## powered clusters as its prior (assumes rho does not depend strongly on cell type;
  ## flagged, not hidden). Its own per-cell evidence, if any, still counts.
  weak <- diag$low_power
  if (any(weak) && any(!weak)) {
    ## weighted MEDIAN, by how much soup each powered cluster's absent set
    ## covers: a count-weighted mean let one small cluster whose fit had run
    ## away (high rho, hence many counts) set the value for a large one
    w <- (diag$n_cells * diag$absent_soup_frac)[!weak]
    diag$rho[weak] <- if (sum(w) > 0) .weighted_median(diag$rho[!weak], w) else stats::median(diag$rho[!weak])
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
  attr(out, "soup_mismatch") <- soup_only_frac > 0.05
  attr(out, "soup_only") <- soup_only
  attr(out, "soup_only_frac") <- soup_only_frac
  out
}

.weighted_median <- function(x, w) {
  o <- order(x); cw <- cumsum(w[o]) / sum(w)
  x[o][which(cw >= 0.5)[1]]
}

## Absent-set fixed point: rho = sum_A X / sum_A e, where A = genes consistent
## with pure soup at the current rho (one-sided: not significantly ABOVE;
## two-sided: also not significantly BELOW).
.absent_fixed_point <- function(X, e, ok, rho, alpha, max_iter, two_sided, min_expected = 0,
                                n_top = 50) {
  A <- rep(FALSE, length(X))
  ## the best-covered soup genes stay eligible whatever rho is, so that small
  ## data sets and a true rho near zero still have a candidate set
  top <- ok & e >= sort(e[ok], decreasing = TRUE)[min(n_top, sum(ok))]
  for (it in seq_len(max_iter)) {
    cons <- ok & stats::ppois(X - 1, rho * e, lower.tail = FALSE) > alpha
    if (two_sided) cons <- cons & stats::ppois(X, rho * e) > alpha
    ## when rho is so low that no gene reaches min_expected (true rho near
    ## zero, or very small data), coverage cannot be asked for: every
    ## consistent gene counts, as before
    covered <- rho * e >= min_expected
    new_A <- if (any(covered[ok])) cons & (covered | top) else cons
    if (!any(new_A)) new_A <- cons
    if (!any(new_A)) break
    A <- new_A
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
