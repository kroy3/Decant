## decant.R -- the assembled pipeline.
## Composes ONLY the components that passed their gate. Each capability switches
## on when the data that justifies it is present, and the mass-conservation
## guarantee holds end to end. The DE module is intentionally excluded by default
## (it failed its gate); it is reachable as decant_de_experimental() with a loud
## warning, not wired into the main path.
##
## Sparse inputs (dgCMatrix) stay sparse throughout; no dense genes x cells
## intermediate is built on the default (non-splice) path.

#' Run Decant.
#'
#' @param obs_total genes x cells total counts (base matrix or dgCMatrix).
#' @param empties EITHER a genes x droplets matrix (single sample) OR a list of
#'   such matrices, one per sample, which triggers hierarchical pooling (Gap 3).
#'   See [split_droplets()] to obtain these from a raw 10x matrix.
#' @param sample_of length-cells sample index (required if empties is a list).
#' @param clusters optional length-cells cluster labels (e.g. from Seurat/SCE).
#'   Used by the rho estimator. If NULL, k-means with `k` clusters is run.
#' @param k number of clusters when `clusters` is NULL.
#' @param rho optional user-supplied per-cell (or scalar) contamination
#'   fraction; skips estimation entirely.
#' @param correction correction rule, see [correct_counts()]. The default,
#'   "bayes_cell", has the lowest error on simulated ground truth
#'   ([gate_correction()]). On real PBMC data (`validation/benchmark_pbmc_six_tools.R`) it
#'   removes 88 to 95% of marker counts that ambient RNA could explain in cell
#'   types that do not express them, and leaves most of the counts that are
#'   too high in a cell to be ambient (likely doublets), which "bayes" removed.
#' @param cap_test for "bayes_cell": "entry" or "cluster", how the test that
#'   sets aside counts too high to be ambient is applied; see
#'   [correct_counts()].
#' @param obs_unspliced,obs_spliced genes x cells layers; if both given with
#'   empty layers, correction is splice-aware (Gap 1).
#' @param empties_unspliced,empties_spliced layer-resolved empties (matrix or
#'   list per sample), required for splice-aware mode.
#' @param allelic optional list(own, other, donor, abund) for genotype-aware rho
#'   (Gap 5); overrides count-based rho when present.
#' @param basis optional genes x K program basis; if given, returns a lysis
#'   diagnostic (Gap 4).
#' @param verbose print a short run summary.
#' @return an object of class `decant`: a list with $corrected (same class as
#'   the input; + $corrected_unspliced/$corrected_spliced in splice mode),
#'   $rho, $rho_diagnostics (per-cluster, from [estimate_rho_cluster()]),
#'   $ambient (genes x samples), $clusters, $lysis (if basis), $modules.
#' @export
decant <- function(obs_total, empties, sample_of = NULL, clusters = NULL, k = 8,
                   rho = NULL, correction = c("bayes_cell", "bayes", "subtract", "redistribute", "bayes_absent"),
                   cap_test = c("entry", "cluster"),
                   obs_unspliced = NULL, obs_spliced = NULL,
                   empties_unspliced = NULL, empties_spliced = NULL,
                   allelic = NULL, basis = NULL, verbose = FALSE) {
  correction <- match.arg(correction); cap_test <- match.arg(cap_test)
  obs_total <- .as_counts(obs_total, "obs_total")
  modules <- character(0)
  multi <- is.list(empties) && !is.data.frame(empties)
  G <- nrow(obs_total); N <- ncol(obs_total)
  if (N < 2) stop("need at least 2 cells", call. = FALSE)

  ## ---- ambient profile (hierarchical pooling if multi-sample) [GAP 3] ----
  if (multi) {
    if (is.null(sample_of)) stop("`sample_of` is required when `empties` is a list", call. = FALSE)
    if (length(sample_of) != N) stop("length(sample_of) must equal ncol(obs_total)", call. = FALSE)
    empties <- lapply(empties, .as_counts, what = "empties")
    if (any(vapply(empties, nrow, 1L) != G)) stop("every empties matrix needs nrow(obs_total) genes", call. = FALSE)
    sample_idx <- if (is.numeric(sample_of)) as.integer(sample_of)
                  else match(sample_of, if (is.null(names(empties))) unique(sample_of) else names(empties))
    if (anyNA(sample_idx) || any(sample_idx < 1L | sample_idx > length(empties)))
      stop("`sample_of` does not map onto `empties`", call. = FALSE)
    pp <- pool_soup(empties)
    soup_by_sample <- pp$pooled
    soup_umis <- vapply(empties, function(e) sum(.col_sums(e)), 1)
    modules <- c(modules, "hierarchical_ambient")
  } else {
    empties <- .as_counts(empties, "empties")
    if (nrow(empties) != G) stop("empties needs nrow(obs_total) genes", call. = FALSE)
    soup_by_sample <- matrix(ambient_global(empties), ncol = 1)
    soup_umis <- sum(.col_sums(empties))
    sample_idx <- rep(1L, N)
  }
  rownames(soup_by_sample) <- rownames(obs_total)

  ## ---- clusters (needed by the rho estimator and bayes correction) ----
  need_clusters <- (is.null(rho) && is.null(allelic)) || correction %in% c("bayes", "bayes_cell", "bayes_absent")
  if (is.null(clusters) && need_clusters) {
    clusters <- quick_labels(obs_total, k = k)
    modules <- c(modules, "kmeans_clusters")
  }
  if (!is.null(clusters) && length(clusters) != N)
    stop("length(clusters) must equal ncol(obs_total)", call. = FALSE)

  ## ---- per-cell rho: user > allelic [GAP 5] > cluster/absent-gene ----
  rho_diag <- NULL; mismatch <- FALSE
  if (!is.null(rho)) {
    if (length(rho) == 1L) rho <- rep(rho, N)
    if (length(rho) != N || any(rho < 0 | rho > 1)) stop("`rho` must be in [0,1], length 1 or ncol(obs_total)", call. = FALSE)
    modules <- c(modules, "user_rho")
  } else if (!is.null(allelic)) {
    rho <- estimate_rho_allelic(allelic)
    modules <- c(modules, "allelic_rho")
  } else {
    ## cell-weighted mean soup (never expand to genes x cells)
    w <- tabulate(sample_idx, nbins = ncol(soup_by_sample)) / N
    soup_mean <- as.numeric(soup_by_sample %*% w)
    rho <- estimate_rho_cluster(obs_total, soup_mean, clusters)
    rho_diag <- attr(rho, "clusters")
    mismatch <- isTRUE(attr(rho, "soup_mismatch"))
    so <- attr(rho, "soup_only"); so_frac <- attr(rho, "soup_only_frac")
    unidentifiable <- isTRUE(attr(rho, "unidentifiable"))
    rho <- as.numeric(rho)
    modules <- c(modules, "cluster_rho")
    if (any(so)) {
      ## correct with the same soup rho was estimated against
      soup_by_sample[so, ] <- 0
      soup_by_sample <- sweep(soup_by_sample, 2, colSums(soup_by_sample), "/")
      modules <- c(modules, "soup_only_genes_dropped")
    }
    if (mismatch)
      warning(sprintf("%.0f%% of the soup is in %d gene(s) the cells do not contain ", 100 * so_frac, sum(so)),
              "(empties may hold debris or cell types absent from the cell set). ",
              "Those genes were dropped from the soup, but other genes may be ",
              "affected too, in which case rho is UNDERestimated. ",
              "See $rho_diagnostics; consider supplying `rho`.", call. = FALSE)
    if (unidentifiable) {
      warning("rho could not be identified: no cluster has genes that are cleanly ",
              "absent, or most cells are in clusters whose estimate is implausibly high ",
              "(> 0.5). This happens ",
              "when the soup resembles the cells, e.g. a single cell type. The estimate ",
              "(median ", sprintf("%.2f", stats::median(rho)), ") is NOT reliable and correction may ",
              "remove real signal. Supply `rho` from an independent source (genotypes, a ",
              "species-mixing control, or a value from comparable data).", call. = FALSE)
    } else if (any(rho_diag$low_power)) {
      n_imp <- sum(rho_diag$implausible)
      warning(sum(rho_diag$low_power), " cluster(s) have almost no genes that are ",
              "cleanly absent", if (n_imp) sprintf(" or an implausible rho above 0.5 (%d)", n_imp) else "",
              ", so rho is weakly identified there; they borrow rho from the other ",
              "clusters. See $rho_diagnostics.", call. = FALSE)
    }
  }

  ## ---- correction ----
  splice_mode <- !is.null(obs_unspliced) && !is.null(obs_spliced) &&
                 !is.null(empties_unspliced) && !is.null(empties_spliced)
  out <- list(rho = rho, rho_diagnostics = rho_diag, soup_mismatch = mismatch,
              rho_unidentifiable = if (is.null(rho_diag)) FALSE else unidentifiable,
              ambient = soup_by_sample,
              clusters = clusters)

  if (splice_mode) {
    eu <- if (is.list(empties_unspliced)) do.call(cbind, empties_unspliced) else empties_unspliced
    es <- if (is.list(empties_spliced))   do.call(cbind, empties_spliced)   else empties_spliced
    obs_unspliced <- .as_counts(obs_unspliced, "obs_unspliced")
    obs_spliced <- .as_counts(obs_spliced, "obs_spliced")
    sa <- correct_splice_aware(obs_unspliced, obs_spliced, eu, es, rho, drop_zeros = FALSE)
    .check_mass(obs_unspliced, sa$unspliced)
    .check_mass(obs_spliced, sa$spliced)
    out$corrected_unspliced <- .drop0(sa$unspliced)
    out$corrected_spliced   <- .drop0(sa$spliced)
    out$corrected <- out$corrected_unspliced + out$corrected_spliced
    modules <- c(modules, "splice_aware_correction")
  } else if (ncol(soup_by_sample) == 1L) {
    out$corrected <- correct_counts(obs_total, soup_by_sample[, 1], rho,
                                    method = correction, clusters = clusters,
                                    drop_zeros = FALSE, cap_test = cap_test, soup_umis = soup_umis[1])
    modules <- c(modules, paste0("correction_", correction))
  } else {
    ## per-sample soup: correct each sample's cells with that sample's soup
    out$corrected <- obs_total
    if (!.is_sparse(out$corrected)) storage.mode(out$corrected) <- "double"
    for (s in unique(sample_idx)) {
      idx <- which(sample_idx == s)
      out$corrected[, idx] <- correct_counts(obs_total[, idx, drop = FALSE],
                                             soup_by_sample[, s], rho[idx],
                                             method = correction,
                                             clusters = clusters[idx],
                                             drop_zeros = FALSE, cap_test = cap_test, soup_umis = soup_umis[s])
    }
    modules <- c(modules, paste0("correction_", correction))
  }

  ## ---- lysis diagnostic [GAP 4] ----
  if (!is.null(basis)) {
    st <- ambient_structured(if (multi) do.call(cbind, empties) else empties, basis)
    out$lysis <- st$lysis
    modules <- c(modules, "lysis_diagnostic")
  }

  ## hard guarantee check, element-wise (cheap slot comparison when sparse)
  if (!splice_mode) {
    .check_mass(obs_total, out$corrected)
    out$corrected <- .drop0(out$corrected)
  }
  out$modules <- modules
  class(out) <- "decant"
  if (verbose) print(out)
  out
}

#' Print a decant result.
#' @param x a `decant` object.
#' @param ... unused.
#' @method print decant
#' @export
print.decant <- function(x, ...) {
  cat("<decant> ", nrow(x$corrected), " genes x ", ncol(x$corrected), " cells",
      if (.is_sparse(x$corrected)) " (sparse)" else "", "\n", sep = "")
  cat("  modules : ", paste(x$modules, collapse = ", "), "\n", sep = "")
  q <- stats::quantile(x$rho, c(0.1, 0.5, 0.9))
  cat(sprintf("  rho     : median %.3f (10%%-90%%: %.3f-%.3f)\n", q[2], q[1], q[3]))
  if (isTRUE(x$soup_mismatch))
    cat("  WARNING : soup/cell mismatch (soup-only genes dropped; rho may be UNDERestimated)\n")
  if (isTRUE(x$rho_unidentifiable))
    cat("  WARNING : rho NOT identifiable from these data; supply `rho`\n")
  else if (!is.null(x$rho_diagnostics) && any(x$rho_diagnostics$low_power))
    cat("  WARNING : ", sum(x$rho_diagnostics$low_power),
        " cluster(s) without usable evidence for rho; they borrow it from the others\n", sep = "")
  invisible(x)
}

#' EXPERIMENTAL and OFF by default: failed its gate (made ambient-driven DE
#' false positives worse via covariate collinearity). Kept only so the negative
#' result is reproducible. Superseded by [ambient_de()] / [DecantDE()].
#' @param ... passed to [de_ambient_aware()].
#' @export
decant_de_experimental <- function(...) {
  warning("decant_de_experimental FAILED its benchmark gate (higher false-positive ",
          "rate than naive). Provided for reproducibility only; do not use.")
  de_ambient_aware(...)
}
