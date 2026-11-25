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
#' @param max_rho the largest contamination treated as identifiable. A
#'   cluster whose estimate exceeds it is flagged (`implausible`) and borrows
#'   rho from the other clusters, like a low-power cluster. Estimates above
#'   about 0.5 mean the cells would be mostly soup, which points to a soup
#'   that resembles the cells (a single cell type, or a cluster that dominates
#'   the soup) rather than to real contamination: on a mouse-only cell line
#'   with 0.4% true contamination the estimator returned 0.5 to 0.94.
#' @param overdispersed test genes against a negative binomial whose extra
#'   gene-to-gene variance is estimated from the data, instead of a Poisson.
#'   On real data, genes a cluster does not express scatter around rho by
#'   about 25-30% (soup and capture differ slightly gene by gene), far more
#'   than Poisson noise. A Poisson test then rejects the absent genes that
#'   fall above the current rho and keeps those below, and the fixed point
#'   settles near the bottom of their spread: on a human/mouse mixture this
#'   halved rho. The variance is estimated from the lower tail of the absent
#'   genes, which cannot contain expressed genes. With no extra variance the
#'   test reduces to the Poisson one.
#' @param weak_cosine,weak_kappa thresholds above which a cluster's profile
#'   counts as close to the soup (`collinear`): the weighted cosine between
#'   the cluster's raw profile and the soup (0.9: a two-profile fit's
#'   variance inflation 1 / (1 - c^2) then exceeds 5) and the Gram condition
#'   number (30, i.e. a cosine of 0.935, so it only matters if `weak_cosine`
#'   is raised).
#' @param weak_ci_width a cluster whose rho has a 95% interval wider than
#'   this is weakly identified (0.05: wider than typical contamination
#'   itself, 0.01-0.1).
#' @param weak_rule "evidence" (default): a cluster is `weak` when the
#'   estimator's own evidence is thin (interval wider than `weak_ci_width`,
#'   almost no absent genes, or an implausible own estimate). "any": also
#'   when it is `collinear`. A collinear cluster whose absent genes still
#'   cover enough soup is usually estimated well, so "any" raises false
#'   alarms on soup-dominant cell types.
#' @return numeric vector of per-cell rho, with attribute `"clusters"`: a
#'   data.frame of per-cluster diagnostics (n_cells, rho, n_absent_genes,
#'   absent_soup_frac, absent_counts, low_power, borrowed, rho_two_sided,
#'   dispersion: the extra squared coefficient of variation used,
#'   implausible: rho above `max_rho`; identifiability: cosine, kappa, js,
#'   fisher_umi, n_eff, rho_se, ci_width, collinear, and weak / weak_reason /
#'   action_taken), attribute `"unidentifiable"` (no
#'   cluster had usable evidence, or more than half the cells are in
#'   implausible clusters; rho should then be supplied), and
#'   attributes `"soup_mismatch"` (logical), `"soup_only"` (logical per gene)
#'   and `"soup_only_frac"` (their share of the soup). Low-power clusters take the
#'   evidence-weighted rho of the other clusters (`borrowed = TRUE`).
#'   `soup_mismatch = TRUE` means more than 5% of the soup was in genes no
#'   cell contains; those genes are dropped, but the empties deserve a look.
#' @export
estimate_rho_cluster <- function(observed, soup, clusters, alpha = 0.01,
                                 min_absent_soup = 0.02, max_iter = 100, min_expected = 5,
                                 floor_rho = 0.005, overdispersed = TRUE, max_rho = 0.5,
                                 weak_cosine = 0.9, weak_kappa = 30, weak_ci_width = 0.05,
                                 weak_rule = c("evidence", "any")) {
  weak_rule <- match.arg(weak_rule)
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
    rho <- rho2 <- phi <- 0; A <- rep(FALSE, length(X))
    if (sum(X) > 0 && any(ok)) {
      r <- X / pmax(e, 1e-300)
      ## start among well-covered genes when there are enough of them (real
      ## data: most genes have almost no soup, and the upper half by soup is
      ## still mostly zeros); otherwise the upper half by soup
      pw <- ok & e >= 4 * min_expected
      if (sum(pw) < 50) pw <- ok & e >= stats::quantile(e[ok], 0.5)
      lo <- pw & r <= stats::quantile(r[pw], 0.1)
      rho0 <- sum(X[lo]) / sum(e[lo])
      ## with sparse counts the lowest decile is all zeros; rho = 0 is an
      ## absorbing state of the fixed point (any gene with a count is then
      ## "significantly above" soup), so start no lower than floor_rho
      rho0 <- max(rho0, floor_rho)
      fit <- .absent_fixed_point(X, e, ok, rho0, alpha, max_iter, two_sided = FALSE, min_expected,
                                 overdispersed = overdispersed)
      rho <- fit$rho; A <- fit$A; phi <- fit$phi
      ## two-sided variant, reported as a diagnostic only
      rho2 <- .absent_fixed_point(X, e, ok, rho0, alpha, max_iter, two_sided = TRUE)$rho
    }
    absent[[i]] <- which(A)
    SA <- sum(soup[absent[[i]]])
    idf <- .identifiability(X, soup, rho, e, A, phi)
    diag[[i]] <- data.frame(cluster = ks[i], n_cells = length(idx), rho = rho,
                            n_absent_genes = length(absent[[i]]), absent_soup_frac = SA,
                            absent_counts = sum(X[absent[[i]]]),
                            low_power = SA < min_absent_soup,
                            borrowed = FALSE,
                            rho_two_sided = rho2,
                            dispersion = phi,
                            implausible = rho > max_rho,
                            idf,
                            stringsAsFactors = FALSE)
  }
  diag <- do.call(rbind, diag)

  ## A cluster with (almost) no absent genes -- typically the fragile type that
  ## dominates the soup itself -- cannot identify its own rho. Rather than keep
  ## an arbitrary value, it borrows the evidence-weighted median rho of the
  ## powered clusters as its prior (assumes rho does not depend strongly on cell type;
  ## flagged, not hidden). Its own per-cell evidence, if any, still counts.
  ## an implausibly high rho is treated like a low-power cluster
  diag$low_power <- diag$low_power | diag$implausible
  weak <- diag$low_power
  if (any(weak) && any(!weak)) {
    ## weighted MEDIAN, by how much soup each powered cluster's absent set
    ## covers: a count-weighted mean let one small cluster whose fit had run
    ## away (high rho, hence many counts) set the value for a large one
    w <- (diag$n_cells * diag$absent_soup_frac)[!weak]
    diag$rho[weak] <- if (sum(w) > 0) .weighted_median(diag$rho[!weak], w) else stats::median(diag$rho[!weak])
    diag$borrowed[weak] <- TRUE
  }

  ## weak identifiability: the reasons, and what was done about it
  wi <- .weak_identifiability(diag, weak_cosine, weak_kappa, weak_ci_width, weak_rule)
  diag$collinear <- wi$collinear
  diag$weak_reason <- wi$reason
  diag$weak <- nzchar(diag$weak_reason)
  diag$action_taken <- ifelse(diag$borrowed, "rho borrowed from the other clusters",
                              ifelse(diag$weak, "none: rho kept, treat it with caution", ""))

  ## pass 2: per-cell empirical-Bayes shrinkage toward the cluster rho
  ## (Gamma prior, Poisson counts on that cluster's absent genes)
  for (i in seq_along(ks)) {
    idx <- which(clusters == ks[i])
    gA <- absent[[i]]; rho <- diag$rho[i]
    ## an implausible cluster's "absent" genes are not absent: they cannot
    ## serve as per-cell evidence, so its cells keep the borrowed rho
    if (length(gA) > 0 && rho > 0 && !diag$implausible[i]) {
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
  ## rho is not estimable from these counts when no cluster has usable
  ## evidence, or when most cells sit in clusters whose estimate was
  ## implausible (the soup resembles the cells, so the remaining clusters'
  ## estimates are not trustworthy either)
  imp_share <- sum(diag$n_cells[diag$implausible]) / sum(diag$n_cells)
  attr(out, "unidentifiable") <- all(diag$low_power) || imp_share > 0.5
  out
}

## Identifiability diagnostics for one cluster (Objective C).
## X: cluster gene sums; soup: profile; rho: cluster estimate; e: soup *
## cluster UMIs; A: absent set; phi: extra squared CV used by the test.
##  * cosine: weighted cosine between the cluster's RAW profile p and the soup,
##    weights 1/q with q = (p + soup) / 2 (the multinomial Fisher metric at the
##    midpoint; weighting by 1/p alone would let any soup gene with a sampling
##    zero in the cluster set it to 0). Near 1, the cluster looks like the soup:
##    with p = (1 - rho) theta + rho soup, theta is then ~ soup and rho cannot
##    be separated from expression. The soup-subtracted profile is not used:
##    subtracting rho * soup biases the similarity down by construction.
##  * kappa = (1 + |c|) / (1 - |c|): condition number of the 2 x 2 Gram matrix
##    of the two unit profiles; it grows like the variance inflation of rho in
##    a two-profile fit.
##  * js: Jensen-Shannon divergence (nats) between p and the soup.
##  * fisher_umi: Fisher information about rho per UMI, sum (theta - soup)^2 / m,
##    m = (1 - rho) theta + rho soup, theta = soup-subtracted profile. It
##    treats theta as known, so it is an upper bound on what the data carry.
##    n_eff: (sum w)^2 / sum w^2 of its per-gene terms: how many genes carry it.
##  * rho_se, ci_width: sampling SD of the absent-gene estimator
##    rho = sum_A X / sum_A e under the negative binomial used by the test,
##    Var = (rho sum_A e + phi sum_A (rho e)^2) / (sum_A e)^2, and the width of
##    the 95% interval. NA when no gene is absent.
.identifiability <- function(X, soup, rho, e, A, phi) {
  tot <- sum(X)
  if (tot <= 0) return(data.frame(cosine = NA_real_, kappa = NA_real_, js = NA_real_, fisher_umi = NA_real_,
                                  n_eff = NA_real_, rho_se = NA_real_, ci_width = NA_real_))
  p <- X / tot
  mid <- (p + soup) / 2; u <- mid > 0
  cosine <- sum(p[u] * soup[u] / mid[u]) / sqrt(sum(p[u]^2 / mid[u]) * sum(soup[u]^2 / mid[u]))
  kappa <- (1 + cosine) / (1 - cosine)
  kl <- function(a) { v <- a > 0; sum(a[v] * log(a[v] / mid[v])) }
  js <- 0.5 * kl(p) + 0.5 * kl(soup)
  th <- pmax(X - rho * e, 0); th <- if (sum(th) > 0) th / sum(th) else p
  m <- (1 - rho) * th + rho * soup; v <- m > 0
  w <- (th[v] - soup[v])^2 / m[v]
  EA <- sum(e[A])
  se <- if (EA > 0) sqrt((rho * EA + phi * sum((rho * e[A])^2)) / EA^2) else NA_real_
  data.frame(cosine = cosine, kappa = kappa, js = js, fisher_umi = sum(w),
             n_eff = if (sum(w) > 0) sum(w)^2 / sum(w^2) else 0,
             rho_se = se, ci_width = 2 * stats::qnorm(0.975) * se)
}

## Which clusters are weakly identified, and why (Objective C). Two kinds of
## reason: the estimator's own evidence is thin (wide interval, almost no
## absent genes, implausible own estimate), and the cluster's profile is close
## to the soup's (cosine / kappa: "collinear"). With rule "evidence" (default)
## only the first kind makes a cluster weak: in simulation a cluster at cosine
## 0.93 whose absent genes still covered 4.7% of the soup was estimated at
## 0.100 against a true 0.105 (95% CI width 0.009). Rule "any" also flags
## collinear clusters.
.weak_identifiability <- function(d, max_cosine = 0.9, max_kappa = 30, max_ci_width = 0.05,
                                  rule = c("evidence", "any")) {
  rule <- match.arg(rule)
  coll <- (d$cosine > max_cosine) %in% TRUE | (d$kappa > max_kappa) %in% TRUE
  why <- vapply(seq_len(nrow(d)), function(i) paste(c(
    if (isTRUE(d$ci_width[i] > max_ci_width)) sprintf("95%% CI width %.3f > %g", d$ci_width[i], max_ci_width),
    if (isTRUE(d$implausible[i])) "own estimate above max_rho",
    if (isTRUE(d$low_power[i]) && !isTRUE(d$implausible[i])) "almost no absent genes",
    if (rule == "any" && coll[i]) sprintf("profile close to the soup (cosine %.2f, kappa %.0f)", d$cosine[i], d$kappa[i])),
    collapse = "; "), "")
  list(reason = why, collinear = coll)
}

.weighted_median <- function(x, w) {
  o <- order(x); cw <- cumsum(w[o]) / sum(w)
  x[o][which(cw >= 0.5)[1]]
}

## Absent-set fixed point: rho = sum_A X / sum_A e, where A = genes consistent
## with pure soup at the current rho (one-sided: not significantly ABOVE;
## two-sided: also not significantly BELOW).
.absent_fixed_point <- function(X, e, ok, rho, alpha, max_iter, two_sided, min_expected = 0,
                                n_top = 50, overdispersed = FALSE, min_disp_expected = 50) {
  A <- rep(FALSE, length(X)); phi <- 0
  ## the best-covered soup genes stay eligible whatever rho is, so that small
  ## data sets and a true rho near zero still have a candidate set
  top <- ok & e >= sort(e[ok], decreasing = TRUE)[min(n_top, sum(ok))]
  for (it in seq_len(max_iter)) {
    if (overdispersed) phi <- .lower_half_dispersion(X, e, A, rho, min_disp_expected)
    cons <- ok & if (phi > 0)
      stats::pnbinom(X - 1, size = 1 / phi, mu = rho * e, lower.tail = FALSE) > alpha
    else stats::ppois(X - 1, rho * e, lower.tail = FALSE) > alpha
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
  list(rho = rho, A = A, phi = phi)
}

## Extra (beyond Poisson) squared coefficient of variation of absent genes'
## counts around rho * e, estimated from the LOWER tail of the current absent
## set's log ratios (counts / soup expectation): the spread between its 10th
## and 25th percentiles, scaled to a normal SD. Expressed genes can only lie
## above the absent genes, so the lower tail is free of them. Two worse
## choices were tried and rejected:
##  * spread measured against the current rho: when rho is too low only the
##    far tail lies below it, the spread looks tiny and rho stays low;
##  * median minus lower quartile: weakly expressed genes that enter the
##    absent set lift the median, inflate the spread, let in more of them,
##    and rho runs high (gate_rho failed in its stress setting).
## Only genes with at least `min_n` expected counts are used, and their mean
## Poisson variance (1 / mu on the log scale) is subtracted.
.lower_half_dispersion <- function(X, e, A, rho, min_n = 50) {
  mu <- rho * e
  use <- A & mu >= min_n & X > 0
  if (sum(use) < 20) return(0)
  lr <- log(X[use] / e[use])
  q <- stats::quantile(lr, c(0.10, 0.25), names = FALSE)
  s <- (q[2] - q[1]) / (stats::qnorm(0.25) - stats::qnorm(0.10))
  max(s^2 - mean(1 / mu[use]), 0)
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
