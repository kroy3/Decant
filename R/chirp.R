## chirp.R -- ambient-aware ligand-receptor calling.
##
## The failure this addresses: a ligand made in abundance by one cell type
## spills into the ambient pool and then shows up, at low level, in every
## droplet. A caller that asks "is the ligand detected in cluster k?" then
## reports every cluster as a sender. chirp() asks instead whether a cluster
## has more ligand (or receptor) than its measured ambient RNA can explain,
## and scores pairs on the excess only.
##
## Status: experimental. Validated on simulated ground truth only
## (tests/testthat/test-chirp.R). Calls are only as good as rho: if rho is
## underestimated, ambient counts leak into the "native" excess. `margin`
## inflates the ambient expectation to absorb that. On simulations (2 seeds,
## rho 0.05-0.3) margin = 2 gave no false calls with rho at 2/3 of its value;
## margin = 1.5 gave 17-100. The cost: pairs whose weaker side is expressed at
## 5-20 per 10k are mostly missed (recall 0.1-0.5); >= 20 per 10k, all found.

#' Ambient-aware ligand-receptor inference between clusters ("chirp").
#'
#' For each gene and cluster, compares the observed cluster total with the
#' count ambient RNA is expected to contribute, `sum_c rho_c * T_c * soup_g`
#' (per sample when `soup` has one column per sample). A gene counts as
#' natively expressed in a cluster when the observed total is significantly
#' above `margin` times that expectation (Poisson upper tail, Benjamini-Hochberg
#' across genes x clusters). Pair scores use native expression, the observed
#' total minus the expected ambient count, per 10,000 native UMIs.
#'
#' @param observed genes x cells raw counts (matrix or dgCMatrix), uncorrected.
#' @param pairs data.frame with columns `ligand` and `receptor` (gene names
#'   matching `rownames(observed)`). A multi-subunit receptor or ligand can be
#'   given as genes joined by `"+"`; the complex is scored by its least
#'   expressed subunit and called only if every subunit is called.
#' @param fit optional `decant` object; supplies `soup`, `rho` and `clusters`
#'   when those are not given.
#' @param soup ambient profile: a length-genes vector, or genes x samples matrix
#'   with `sample_of`.
#' @param rho per-cell (or scalar) contamination fraction.
#' @param clusters length-cells cluster labels.
#' @param sample_of length-cells sample index into the columns of `soup`
#'   (required when `soup` has more than one column).
#' @param alpha FDR level for the native-expression call.
#' @param margin multiplier on the expected ambient count before testing.
#'   The default, 2, kept false calls at zero in simulations with rho
#'   underestimated by a third; it does not protect if rho is half its true
#'   value. Lower it for more sensitivity only if rho is trusted.
#' @param min_native minimum native expression, in counts per 10,000 native
#'   UMIs, for a call.
#' @return data.frame, one row per pair x sender x receiver, sorted by score:
#'   `ligand`, `receptor`, `sender`, `receiver`, `score` (geometric mean of
#'   native ligand in sender and native receptor in receiver), `call`
#'   (both sides natively expressed), `naive_score` (same on raw counts),
#'   `naive_call` (both sides detected in at least 10% of cells, the usual
#'   rule), `ambient_share_ligand` (expected ambient share of the ligand
#'   counts in the sender) and `p_ligand`, `p_receptor` (BH-adjusted).
#'   Rows with `naive_call & !call` are the interactions ambient RNA would
#'   have produced. Attribute `"native"`: genes x clusters native expression.
#' @examples
#' sim <- simulate_experiment(n_genes = 200, n_types = 4, n_cells = 300, seed = 1)
#' g <- rownames(sim$observed)
#' pairs <- data.frame(ligand = g[1:5], receptor = g[6:10])
#' res <- chirp(sim$observed, pairs, soup = ambient_global(sim$empty),
#'              rho = sim$rho_true, clusters = sim$labels)
#' head(res[res$call, ])
#' @export
chirp <- function(observed, pairs, fit = NULL, soup = NULL, rho = NULL,
                  clusters = NULL, sample_of = NULL, alpha = 0.05,
                  margin = 2, min_native = 1) {
  observed <- .as_counts(observed, "observed")
  if (!is.null(fit)) {
    if (!inherits(fit, "decant")) stop("`fit` must be a decant object", call. = FALSE)
    if (is.null(soup))     soup     <- fit$ambient
    if (is.null(rho))      rho      <- fit$rho
    if (is.null(clusters)) clusters <- fit$clusters
  }
  if (is.null(soup) || is.null(rho) || is.null(clusters))
    stop("supply `soup`, `rho` and `clusters`, or a decant `fit`", call. = FALSE)
  G <- nrow(observed); N <- ncol(observed)
  genes <- rownames(observed)
  if (is.null(genes)) stop("observed needs gene rownames", call. = FALSE)
  if (!is.data.frame(pairs) || !all(c("ligand", "receptor") %in% names(pairs)))
    stop("`pairs` must be a data.frame with columns ligand and receptor", call. = FALSE)
  if (length(clusters) != N) stop("length(clusters) must equal ncol(observed)", call. = FALSE)
  rho <- as.numeric(rho)
  if (length(rho) == 1L) rho <- rep(rho, N)
  if (length(rho) != N || anyNA(rho) || any(rho < 0 | rho >= 1))
    stop("rho must be scalar or length ncol(observed), in [0, 1)", call. = FALSE)
  if (margin < 1) stop("margin must be >= 1", call. = FALSE)

  soup <- as.matrix(soup)
  if (nrow(soup) != G) stop("soup needs nrow(observed) genes", call. = FALSE)
  soup <- sweep(soup, 2, colSums(soup), "/")
  if (ncol(soup) == 1L) {
    sample_idx <- rep(1L, N)
  } else {
    if (is.null(sample_of) || length(sample_of) != N)
      stop("`sample_of` (length ncol(observed)) is required for a per-sample soup", call. = FALSE)
    sample_idx <- if (is.numeric(sample_of)) as.integer(sample_of)
                  else match(sample_of, colnames(soup))
    if (anyNA(sample_idx) || any(sample_idx < 1L | sample_idx > ncol(soup)))
      stop("`sample_of` does not map onto the columns of soup", call. = FALSE)
  }

  ks <- sort(unique(clusters))
  K <- length(ks)
  T_c <- .col_sums(observed)
  X <- .cluster_sums(observed, clusters, ks)                       # genes x K
  D <- .cluster_sums(methods::as(observed > 0, "dMatrix"), clusters, ks)
  n_k <- tabulate(match(clusters, ks), K)
  ## expected ambient counts per gene x cluster
  w <- rho * T_c
  W <- matrix(0, ncol(soup), K)                                    # samples x K
  for (s in unique(sample_idx)) {
    i <- sample_idx == s
    W[s, ] <- vapply(ks, function(k) sum(w[i & clusters == k]), 0)
  }
  A <- soup %*% W                                                  # genes x K
  native_umis <- vapply(ks, function(k) sum((1 - rho[clusters == k]) * T_c[clusters == k]), 0)
  native <- sweep(pmax(X - A, 0), 2, pmax(native_umis, 1), "/") * 1e4
  raw    <- sweep(X, 2, pmax(vapply(ks, function(k) sum(T_c[clusters == k]), 0), 1), "/") * 1e4

  ## native-expression test: X ~ Poisson(margin * A) under "ambient only"
  lam <- margin * A
  p <- stats::ppois(X - 1, lam, lower.tail = FALSE)
  p[X == 0] <- 1
  p[lam == 0 & X > 0] <- 0
  p_adj <- matrix(stats::p.adjust(as.numeric(p), "BH"), G, K)
  called <- p_adj < alpha & native >= min_native
  detected <- sweep(D, 2, pmax(n_k, 1), "/") >= 0.1
  dimnames(native) <- dimnames(raw) <- dimnames(p_adj) <- dimnames(called) <-
    dimnames(detected) <- dimnames(A) <- dimnames(X) <- list(genes, as.character(ks))

  ## complexes: "A+B" -> least expressed subunit, all subunits called
  side <- function(spec) {
    sub <- strsplit(spec, "+", fixed = TRUE)[[1]]
    miss <- setdiff(sub, genes)
    if (length(miss)) return(NULL)
    agg <- function(M, f) apply(M[sub, , drop = FALSE], 2, f)
    list(native = agg(native, min), raw = agg(raw, min),
         called = agg(called, all), detected = agg(detected, all),
         p = agg(p_adj, max), amb = colSums(A[sub, , drop = FALSE]) /
           pmax(colSums(X[sub, , drop = FALSE]), 1))
  }

  out <- list()
  for (r in seq_len(nrow(pairs))) {
    L <- side(as.character(pairs$ligand[r])); R <- side(as.character(pairs$receptor[r]))
    if (is.null(L) || is.null(R)) next
    grid <- expand.grid(s = seq_len(K), t = seq_len(K))
    out[[length(out) + 1L]] <- data.frame(
      ligand = as.character(pairs$ligand[r]), receptor = as.character(pairs$receptor[r]),
      sender = as.character(ks)[grid$s], receiver = as.character(ks)[grid$t],
      score = sqrt(L$native[grid$s] * R$native[grid$t]),
      call = L$called[grid$s] & R$called[grid$t],
      naive_score = sqrt(L$raw[grid$s] * R$raw[grid$t]),
      naive_call = L$detected[grid$s] & R$detected[grid$t],
      ambient_share_ligand = pmin(L$amb[grid$s], 1),
      p_ligand = L$p[grid$s], p_receptor = R$p[grid$t],
      stringsAsFactors = FALSE)
  }
  skipped <- nrow(pairs) - length(out)
  if (skipped > 0)
    warning(skipped, " pair(s) skipped: gene(s) not in rownames(observed)", call. = FALSE)
  res <- if (length(out)) do.call(rbind, out) else
    data.frame(ligand = character(), receptor = character(), sender = character(),
               receiver = character(), score = numeric(), call = logical(),
               naive_score = numeric(), naive_call = logical(),
               ambient_share_ligand = numeric(), p_ligand = numeric(), p_receptor = numeric())
  res <- res[order(!res$call, -res$score), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "native") <- native
  res
}
