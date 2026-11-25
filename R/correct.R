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
##  "bayes_absent" : "bayes_cell" with an absence posterior. The rules above
##                   remove a gene's EXPECTED ambient count and keep the rest
##                   as native. But ambient levels scatter 2 to 4 times more
##                   than Poisson from gene to gene, so a gene a cluster does
##                   not express usually sits above its expectation and keeps
##                   the excess (20k species mixture: at most 0.77 of foreign
##                   signal removable even with the true rho). Here each gene
##                   in each cluster gets a posterior probability of being
##                   absent, from a two-component Poisson-lognormal mixture
##                   over genes of cluster count against ambient expectation
##                   (ambient-only genes scatter around the expectation;
##                   expressed genes lie far above it). An absent gene's
##                   native rate is zero, so all its counts are removed.
##                   Counts too high to be ambient are still set aside
##                   (cap_test = "hybrid": too high anywhere in the cluster,
##                   or too high in the cell itself and at least 5 UMIs).
## Which rule is the default is decided by gate_correction(), not by taste.

#' Remove ambient counts from each cell (mass-conserving).
#'
#' @param observed genes x cells counts (base matrix or dgCMatrix).
#' @param soup length-G ambient probability vector (one soup for all cells;
#'   for per-sample soups use [decant()]).
#' @param rho length-cells contamination fractions.
#' @param method one of "subtract", "redistribute", "bayes", "bayes_cell",
#'   "bayes_absent" (see Details).
#' @param soup_umis for "bayes_absent": total UMIs the soup profile was
#'   estimated from. If given, the absence test compares each gene's count
#'   in a cluster with its count in the empty droplets (a binomial given
#'   their total), which carries the soup's own sampling noise; otherwise
#'   the soup is taken as known.
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
#'   it is not the default. "hybrid" (used by "bayes_absent") sets an entry
#'   aside if it passes the "cluster" test, or the "entry" test with at
#'   least 5 UMIs: a single UMI against a tiny expectation is weak evidence
#'   on its own, but five are not, so doublet-level counts stay protected.
#' @param drop_zeros for sparse input, drop entries corrected to exactly 0.
#'   FALSE keeps the input's sparsity pattern (explicit zeros).
#' @return corrected counts, same class as `observed` (dense or dgCMatrix).
#'   Guaranteed 0 <= corrected <= observed element-wise.
#' @details See the header of `R/correct.R` for the rules and
#'   [gate_correction()] for how they compare against ground truth.
#' @export
correct_counts <- function(observed, soup, rho,
                           method = c("subtract", "redistribute", "bayes", "bayes_cell", "bayes_absent"),
                           clusters = NULL, max_iter = 50, drop_zeros = TRUE, cap_q = 1e-3,
                           cap_test = c("entry", "cluster", "hybrid"), soup_umis = NULL) {
  method <- match.arg(method); cap_test <- match.arg(cap_test)
  observed <- .as_counts(observed, "observed")
  if (length(soup) != nrow(observed)) stop("length(soup) must equal nrow(observed)", call. = FALSE)
  if (length(rho) == 1L) rho <- rep(rho, ncol(observed))
  if (length(rho) != ncol(observed)) stop("length(rho) must equal ncol(observed)", call. = FALSE)
  soup <- soup / sum(soup)
  T_c <- .col_sums(observed)
  where <- paste0("correct_counts(method = \"", method, "\")")
  ## The cell-aware rules run blockwise over columns: for a sparse input no
  ## vector the length of all non-zero entries is built except the result
  ## (hgmm 10k, 55M non-zeros: peak 4.7 GB above the input before).
  if (method %in% c("bayes_cell", "bayes_absent") && .is_sparse(observed)) {
    if (is.null(clusters)) stop("method = '", method, "' needs `clusters`", call. = FALSE)
    v <- .corr_bayes_cell(observed, soup, rho, T_c, clusters, cap_q,
                          if (method == "bayes_absent") "hybrid" else cap_test,
                          absent = method == "bayes_absent", soup_umis = soup_umis)
    observed@x <- .assert_bounds(v, observed@x, where = where)
    return(if (drop_zeros) Matrix::drop0(observed) else observed)
  }
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
      .corr_bayes_cell(observed, soup, rho, T_c, clusters, cap_q, cap_test)
    },
    bayes_absent = {
      if (is.null(clusters)) stop("method = 'bayes_absent' needs `clusters`", call. = FALSE)
      .corr_bayes_cell(observed, soup, rho, T_c, clusters, cap_q, "hybrid",
                       absent = TRUE, soup_umis = soup_umis)
    })
  ## The guarantee holds by construction on every path (each rule removes a
  ## non-negative share of an entry, at most the entry itself). It is checked
  ## here, O(nnz), and a violation is an error, not something clipped away:
  ## only floating-point rounding (a relative 1e-10) is snapped back.
  v <- .assert_bounds(v, nz$v, where = where)
  .nz_set(observed, nz, v, drop = drop_zeros)
}

## O(nnz) runtime check of invariant I1: 0 <= corrected <= observed, finite.
## Violations beyond rounding raise a classed error (`decant_invariant_violation`).
.assert_bounds <- function(v, x, where = "correction", rel_tol = 1e-10) {
  if (!anyNA(v) && all(v >= 0) && all(v <= x)) return(v)     # fast path: exact
  tol <- rel_tol * pmax(abs(x), 1)
  bad <- !is.finite(v) | v < -tol | v > x + tol
  if (any(bad)) {
    i <- which(bad)[1]
    stop(structure(class = c("decant_invariant_violation", "error", "condition"), list(
      message = sprintf(paste0("internal error in %s: %d corrected value(s) outside [0, observed] ",
                               "(first: corrected %g vs observed %g)"), where, sum(bad), v[i], x[i]),
      call = NULL, n_violations = sum(bad))))
  }
  pmin(pmax(v, 0), x)            # rounding only: every value is within tol
}

## "subtract" is defined as clamp(x - expected ambient, 0, x): the clamp at 0 is
## the rule itself (and why it under-removes), not a safety net.
.corr_subtract <- function(nz, soup, rho, T_c) {
  pmax(nz$v - (rho * T_c)[nz$j] * soup[nz$i], 0)
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
.corr_bayes_cell <- function(observed, soup, rho, T_c, clusters, cap_q, cap_test = "entry",
                             absent = FALSE, soup_umis = NULL, absent_min_genes = 50,
                             absent_cap = 0.10, block_nnz = 4e6) {
  ks <- sort(unique(clusters)); K <- length(ks); G <- nrow(observed); N <- ncol(observed)
  cl_idx <- match(clusters, ks)
  n_k <- tabulate(cl_idx, K)
  amb_c <- rho * T_c
  nat_c <- (1 - rho) * T_c
  t_entry <- cap_q * (1 + 1e-12) + 64 * .Machine$double.eps
  ## Column blocks of at most ~block_nnz non-zero entries. Values come back in
  ## the order of .nz(observed) (column-major), which concatenating blocks of
  ## consecutive columns preserves.
  nnz_col <- if (.is_sparse(observed)) diff(observed@p) else colSums(observed != 0)
  blk <- .column_blocks(nnz_col, block_nnz)
  off <- c(0, cumsum(nnz_col))
  ## per block: non-zeros, their expected ambient count, and whether each is
  ## too high to be ambient. P(count >= observed | ambient only) serves both
  ## tests: "entry" = above the (1 - cap_q) Poisson quantile (pu <= cap_q,
  ## the 64-epsilon fuzz matching qpois()); "cluster" = pu times the cells of
  ## the cluster < cap_q.
  block <- function(cols, test = TRUE) {
    sub <- observed[, cols, drop = FALSE]
    nzb <- .nz(sub)
    j <- cols[nzb$j]
    lam <- amb_c[j] * soup[nzb$i]
    if (!test) return(list(i = nzb$i, j = j, v = nzb$v, lam = lam))
    pu <- .ppois_upper(nzb$v, lam, decide = c(0.5 * cap_q / max(n_k), 2 * t_entry))
    high <- switch(cap_test, entry = pu <= t_entry, cluster = pu * n_k[cl_idx[j]] < cap_q,
                   hybrid = pu * n_k[cl_idx[j]] < cap_q | (nzb$v >= 5 & pu <= t_entry))
    list(i = nzb$i, j = j, v = nzb$v, lam = lam, high = high, sub = sub, idx = nzb$idx, cols = cols)
  }
  ## pass 1: sums over entries by (gene, cluster), as a sparse product with
  ## the cells x clusters indicator
  ind <- Matrix::sparseMatrix(i = seq_len(N), j = cl_idx, x = 1, dims = c(N, K))
  by_gk <- function(b, w) {
    M <- b$sub
    if (.is_sparse(M)) M@x <- as.numeric(w) else { M[] <- 0; M[b$idx] <- w }
    m <- as.matrix(M %*% ind[b$cols, , drop = FALSE]); dimnames(m) <- NULL; m
  }
  X_low <- E_high <- N_high <- X_all <- matrix(0, G, K)
  high_all <- logical(off[N + 1])
  for (cols in blk) {
    b <- block(cols)
    if (length(b$v)) high_all[(off[cols[1]] + 1):off[cols[length(cols)] + 1]] <- b$high
    X_low  <- X_low  + by_gk(b, b$v * !b$high)
    E_high <- E_high + by_gk(b, amb_c[b$j] * b$high)
    N_high <- N_high + by_gk(b, nat_c[b$j] * b$high)
    if (absent) X_all <- X_all + by_gk(b, b$v)
  }
  ## counts and expected ambient among the entries that ambient could explain.
  ## Zero entries are always in that set, so start from the cluster totals.
  amb_k <- as.numeric(tapply(amb_c, factor(clusters, levels = ks), sum))
  native <- pmax(X_low - (outer(soup, amb_k) - soup * E_high), 0)
  if (absent) {
    ## absence posterior per gene and cluster; an absent gene has no native rate
    ## judged on the cluster's FULL counts: for an expressed gene nearly
    ## every entry is "high", and the few entries left would look ambient
    W <- vapply(seq_len(K), function(k)
      as.numeric(.absence_posterior(X_all[, k], soup, amb_k[k], soup_umis, min_info_genes = absent_min_genes)),
      numeric(G))
    W <- .absence_degeneracy_cap(matrix(W, G, K), X_all, absent_cap)
    native <- native * (1 - W)
  }
  ## native rate per native UMI, among the same cells: their native depth is
  ## the cluster's minus that of the cells set aside for this gene
  nat_k <- as.numeric(tapply(nat_c, factor(clusters, levels = ks), sum))
  phi <- native / pmax(matrix(nat_k, G, K, byrow = TRUE) - N_high, 1e-300)
  ## pass 2: allocate each entry
  out <- numeric(off[N + 1])
  for (cols in blk) {
    b <- block(cols, test = FALSE)
    if (!length(b$v)) next
    rng <- (off[cols[1]] + 1):off[cols[length(cols)] + 1]
    high <- high_all[rng]
    a <- (1 - rho[b$j]) * phi[cbind(b$i, cl_idx[b$j])]
    bb <- rho[b$j] * soup[b$i]
    ## a + b = 0: no native rate and no soup (a gene absent from the cluster
    ## that the empties never showed); with the absence posterior it is removed
    o <- b$v * a / pmax(a + bb, 1e-300)
    if (absent) o[a + bb <= 0] <- 0
    o[high] <- b$v[high] - pmin(b$lam[high], b$v[high])
    out[rng] <- o
  }
  out
}

## Split columns 1..N into consecutive blocks holding about `target` non-zero
## entries each (at least one column per block).
.column_blocks <- function(nnz_col, target) {
  N <- length(nnz_col)
  if (!N) return(list())
  id <- floor(cumsum(as.numeric(nnz_col)) / max(target, 1))
  split(seq_len(N), cumsum(c(TRUE, diff(id) != 0)))
}

## P(X >= x) for X ~ Poisson(lam), x >= 1 integer counts. Most non-zero
## entries of droplet data are 1, where this is 1 - exp(-lam) (computed as
## -expm1(-lam) for accuracy); the rest go through ppois().
## If `decide = c(lo, hi)` is given, the result is only used to compare with
## thresholds inside [lo, hi]; then, for x >= 2, the bounds
## exp(-lam) lam^x / x! <= P(X >= x) <= lam^x / x! settle most entries
## exactly (below lo or above hi) and ppois() runs only on the rest.
.ppois_upper <- function(x, lam, decide = NULL) {
  out <- -expm1(-lam)
  o <- which(x != 1)
  if (!length(o)) return(out)
  if (!is.null(decide)) {
    xo <- x[o]; lo <- lam[o]
    lub <- xo * log(lo) - lgamma(xo + 1)              # log upper bound
    below <- xo >= 2 & lub < log(decide[1])
    above <- xo >= 2 & (lub - lo) > log(decide[2])     # log lower bound
    out[o[below]] <- exp(lub[below])
    out[o[above]] <- exp(lub[above] - lo[above])
    o <- o[!(below | above)]
  }
  if (length(o)) out[o] <- stats::ppois(x[o] - 1, lam[o], lower.tail = FALSE)
  out
}

## Posterior probability that each gene is absent from a cluster (all of its
## counts ambient), from a two-component mixture over genes. A gene's cluster
## count X is compared with its count C in the empty droplets (soup * M):
## given their total, X ~ Binomial(X + C, p), with
##   p = A r / (A r + M),   r = exp(m + s * z),  z ~ N(0, 1),
## A the cluster's expected ambient UMIs and M the soup's. r is the gene's
## rate in the cluster relative to pure ambient: (m, s) = (m_A, s_A) for
## ambient-only genes (r near 1, with the gene-to-gene scatter of ambient
## levels) and (m_B, s_B) for expressed genes (r far above 1). Conditioning
## on the total removes the gene's unknown soup rate and carries the sampling
## noise of the soup itself, so a gene the empties never showed is not
## treated as ambient on the strength of a guess. Without M (soup_umis) the
## soup is taken as known: X ~ Poisson(soup * A * r).
## The likelihood is integrated over z by Gauss-Hermite quadrature; (m, s)
## come from posterior-weighted moments of log ratios among genes with
## enough counts, with the sampling part of the variance removed. m_A is
## free, so a cluster rho that is too low moves the ambient component
## instead of turning absent genes into expressed ones.
.absence_posterior <- function(X, soup, A, M = NULL, n_iter = 30, min_count = 5, n_nodes = 15,
                               min_info_genes = 50) {
  w <- numeric(length(X))
  C <- if (is.null(M)) NULL else round(soup * M)
  ok <- if (is.null(M)) soup > 0 else (X + C) > 0
  if (A <= 0 || sum(ok & X > 0) < 20) return(w)
  x <- round(X[ok]); e <- soup[ok] * A; cc <- if (is.null(M)) NULL else C[ok]
  ## Information guard. The ambient component's spread s_A is estimated from
  ## genes with at least `min_count` expected ambient counts in this cluster
  ## (and in the empties). With n such genes the log-SD has a standard error
  ## of about 1 / sqrt(2 n): ~10% at n = 50, the precision the 2-SD separation
  ## rule below takes for granted. A cluster with fewer (in practice a
  ## cluster of a handful of cells) cannot tell an absent gene from a weakly
  ## expressed one, and its fit was seen to call nearly every gene absent
  ## (17,543 genes in a 9-cell brain cluster; simulated 1-5 cell clusters:
  ## 93% of the genes called absent were expressed). It keeps the plain rule.
  n_info <- sum(e >= min_count & (if (is.null(M)) TRUE else cc >= min_count))
  if (n_info < min_info_genes) {
    attr(w, "fit") <- c(m_A = NA, s_A = NA, m_B = NA, s_B = NA, used = FALSE, n_info = n_info)
    return(w)
  }
  gh <- .gauss_hermite(n_nodes)
  ## With the empties (binomial form) the success probability does not depend
  ## on the gene, so the likelihood depends on (x, c) only: it is evaluated
  ## once per distinct pair (tens of thousands of genes share a few thousand).
  if (!is.null(M)) {
    key <- x * (max(cc) + 1) + cc
    uk <- unique(key); grp <- match(key, uk); first <- match(uk, key)
    ux <- x[first]; ucc <- cc[first]
  }
  loglik <- function(m, s) {
    xx <- if (is.null(M)) x else ux
    mx <- NULL; L <- vector("list", n_nodes)
    for (q in seq_len(n_nodes)) {
      r <- exp(m + s * gh$x[q])
      L[[q]] <- (if (is.null(M)) stats::dpois(x, e * r, log = TRUE)
                 else stats::dbinom(ux, ux + ucc, A * r / (A * r + M), log = TRUE)) + log(gh$w[q])
      mx <- if (is.null(mx)) L[[q]] else pmax(mx, L[[q]])
    }
    acc <- numeric(length(xx))
    for (q in seq_len(n_nodes)) acc <- acc + exp(L[[q]] - mx)
    out <- mx + log(acc)
    if (is.null(M)) out else out[grp]
  }
  cnt <- if (is.null(M)) rep(Inf, length(x)) else cc
  ## start from the two modes of the log ratios themselves (two-means on
  ## well-covered genes), not from a test against the ambient expectation:
  ## when ambient genes sit well above it (rho too low, or plain scatter with
  ## large counts) such a test rejects most of them, the expressed component
  ## starts as a blend of both groups, and EM does not recover (hgmm 10k).
  lr_all <- log(pmax(x, 0.5) / pmax(e, 1e-300))
  cov <- e > 0 & x >= min_count & cnt >= min_count
  if (sum(cov) < 20) cov <- e > 0 & x >= 1
  ctr <- stats::quantile(lr_all[cov], c(0.1, 0.9), names = FALSE)
  for (i in 1:50) {
    lowg <- abs(lr_all[cov] - ctr[1]) <= abs(lr_all[cov] - ctr[2])
    if (!any(lowg) || all(lowg)) break
    nc <- c(mean(lr_all[cov][lowg]), mean(lr_all[cov][!lowg]))
    if (max(abs(nc - ctr)) < 1e-6) break
    ctr <- nc
  }
  wA <- as.numeric(e > 0 & lr_all < mean(ctr))
  moments <- function(wt, use) {
    u <- use & wt > 1e-6 & e > 0
    if (sum(wt[u]) < 5) return(NULL)
    lr <- log(pmax(x[u], 0.5) / e[u]); ww <- wt[u] / sum(wt[u])
    m <- sum(ww * lr); v <- sum(ww * (lr - m)^2) - sum(ww * (1 / pmax(x[u], 1) + 1 / pmax(cnt[u], 1)))
    c(m = m, s = sqrt(max(v, 1e-4)))
  }
  pA <- c(m = 0, s = 0.5); pB <- NULL
  for (it in seq_len(n_iter)) {
    a <- moments(wA, e >= min_count & cnt >= min_count); b <- moments(1 - wA, x >= min_count & cnt >= min_count)
    if (is.null(b)) b <- moments(1 - wA, x >= min_count)
    if (is.null(b)) { wA[] <- 1; break }          # nothing looks expressed
    if (!is.null(a)) pA <- a
    pB <- b
    if (pB[["m"]] <= pA[["m"]]) break              # components not separable
    piA <- min(max(mean(wA), 1e-4), 1 - 1e-4)
    lA <- log(piA) + loglik(pA[["m"]], pA[["s"]]); lB <- log(1 - piA) + loglik(pB[["m"]], pB[["s"]])
    new <- 1 / (1 + exp(pmin(lB - lA, 700)))
    new[!is.finite(new)] <- 0
    if (max(abs(new - wA)) < 1e-4) { wA <- new; break }
    wA <- new
  }
  ## The posterior is used only when the two groups are well apart: at least
  ## two standard deviations each between the component means. In a small or
  ## heavily contaminated cluster (expressed genes only a few times above
  ## ambient) they overlap and no classification is reliable; the cluster then
  ## keeps the plain rule (w = 0).
  sep_ok <- !is.null(pB) && (pB[["m"]] - pA[["m"]]) >= 2 * (pA[["s"]] + pB[["s"]])
  if (sep_ok) w[ok] <- wA
  attr(w, "fit") <- c(m_A = pA[["m"]], s_A = pA[["s"]], m_B = if (is.null(pB)) NA else pB[["m"]],
                      s_B = if (is.null(pB)) NA else pB[["s"]], used = sep_ok, n_info = n_info)
  w
}

## Degeneracy cap across clusters. Among genes with counts in a cluster, the
## share called absent (posterior > 0.5) may exceed that of the OTHER clusters
## whose call applied by at most `cap` (default 0.10); a cluster above it keeps
## the plain rule (its column set to 0). With no other cluster to compare to,
## nothing is changed (the information guard in .absence_posterior() is then
## the only protection).
.absence_degeneracy_cap <- function(W, X, cap = 0.10) {
  if (!is.finite(cap) || ncol(W) < 2) return(W)
  used <- colSums(W) > 0
  frac <- vapply(seq_len(ncol(W)), function(k) {
    h <- X[, k] > 0; if (!any(h)) 0 else mean(W[h, k] > 0.5) }, 0)
  for (k in which(used)) {
    ref <- frac[used & seq_along(frac) != k]
    if (length(ref) && frac[k] > max(ref) + cap) W[, k] <- 0
  }
  W
}

## Gauss-Hermite nodes and weights for a standard normal (Golub-Welsch)
.gauss_hermite <- function(n) {
  J <- matrix(0, n, n); i <- seq_len(n - 1)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i)
  ev <- eigen(J, symmetric = TRUE)
  list(x = ev$values, w = ev$vectors[1, ]^2)
}
