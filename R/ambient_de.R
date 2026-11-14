## ambient_de.R
## Contamination-aware differential expression.
##
## Correcting counts and then testing does not stop ambient leakage: when the
## soup differs between samples (fragile disease tissue lyses more), the
## residual after subtraction still tracks the condition, and an abundant cell
## type's disease genes are reported as DE in rare cell types (see
## experiments/contamination_aware_de and gate_ambient_de()).
##
## Here the soup is not subtracted; it is MODELLED. For target-cell-type
## pseudobulk y_gs in sample s:
##     y_gs ~ NB( mu_gs = N_s * exp(x_s' beta_g) + E_gs )
## where E_gs = rho_s * T_s * soup_gs is a KNOWN additive ambient expectation:
##   * soup_gs comes from sample s's own empty droplets, shrunk per gene by
##     empirical Bayes only as far as its sampling noise warrants;
##   * rho_s is re-estimated from the target type's own pseudobulk (absent-gene
##     estimator), so errors in upstream cell-level rho do not propagate;
##   * the soup's sampling variance is added to the NB variance;
##   * N_s carries TMM composition factors computed on the soup-removed counts.
## beta is fitted by vectorised Fisher scoring and tested with a
## quasi-likelihood F-test whose residual variances are moderated by empirical
## Bayes (as edgeR's glmQLFTest does), which keeps small-sample tests calibrated.
## Outlier samples are down-weighted per gene (Huber weights), which keeps
## large-cohort tests calibrated when a few samples carry extreme counts.

#' Empirical-Bayes soup profiles per sample.
#'
#' Each sample's soup is shrunk, gene by gene, toward the cross-sample mean by
#' `w = tau2 / (tau2 + v)`, where `v` is the multinomial sampling variance of
#' that sample's estimate and `tau2` the estimated between-sample variance of
#' the true soups. Well-measured samples keep their own soup; poorly measured
#' ones borrow. (A fixed-rule pooling was tried and rejected: it erased real
#' sample-specific soup differences and inflated false positives.)
#'
#' @param empties list of genes x droplets empty-droplet matrices, one per
#'   sample (same genes, same order).
#' @param design optional sample-level design matrix. Each sample's soup is
#'   then shrunk toward the design-fitted soup instead of the overall mean, so
#'   real soup differences between conditions are not shrunk away.
#' @return list(soup = genes x samples, var = genes x samples posterior
#'   variance, umis = total soup UMIs per sample).
#' @export
soup_eb <- function(empties, design = NULL) {
  empties <- lapply(empties, .as_counts, what = "empties")
  if (length(empties) < 2) stop("soup_eb() needs at least 2 samples", call. = FALSE)
  if (length(unique(vapply(empties, nrow, 1L))) != 1) stop("all empties need the same genes", call. = FALSE)
  raw <- vapply(empties, function(e) { p <- .row_sums(e); p / sum(p) }, numeric(nrow(empties[[1]])))
  M <- vapply(empties, function(e) sum(.col_sums(e)), 1)
  v <- sweep(raw * (1 - raw), 2, M, "/")
  S <- length(empties)
  if (is.null(design)) design <- matrix(1, S, 1)
  design <- as.matrix(design)
  if (nrow(design) != S) stop("nrow(design) must equal the number of samples", call. = FALSE)
  ## prior mean: the design-fitted soup, so that a soup that differs by
  ## condition is not shrunk toward a condition-blind average
  H <- design %*% solve(crossprod(design), t(design))
  m <- raw %*% H
  df <- max(S - ncol(design), 1)
  tau2 <- pmax(rowSums((raw - m)^2) / df - rowMeans(v), 0)
  w <- tau2 / (tau2 + v + 1e-300)
  post <- w * raw + (1 - w) * m
  post <- sweep(post, 2, colSums(post), "/")
  dimnames(post) <- list(rownames(empties[[1]]), names(empties))
  list(soup = post, var = w * v, umis = M)
}

#' Contamination-aware differential expression on pseudobulk counts.
#'
#' @param counts genes x samples pseudobulk counts of ONE cell type.
#' @param empties list (one per sample, same order as `counts` columns) of
#'   genes x droplets empty-droplet matrices for the same genes; OR supply
#'   `soup` directly.
#' @param design sample-level design matrix (rows = samples).
#' @param coef column of `design` to test (index or name).
#' @param soup optional precomputed [soup_eb()] result.
#' @param min_count genes are tested only if their soup-removed pseudobulk
#'   passes edgeR's `filterByExpr` at this `min.count`.
#' @param naive also run standard edgeR QL on the raw pseudobulk and report
#'   which of its hits the contamination-aware test does not support.
#' @param robust down-weight outlier samples per gene (Huber weights on
#'   Pearson residuals, with the dispersion re-estimated under those weights).
#'   Guards against the false positives NB tests produce in large cohorts
#'   when a few samples carry extreme counts.
#' @param joint_rho refine each sample's ambient fraction on genes that look
#'   absent from the cell type across ALL samples pooled. A single sample
#'   cannot tell a weakly expressed gene from an absent one, which biases its
#'   ambient fraction upward; the pooled test can.
#' @param soup_design shrink each sample's soup toward the design-fitted soup
#'   (see [soup_eb()]) instead of the overall mean.
#' @param max_ambient_frac genes whose counts in this cell type are, on
#'   average, more than this fraction ambient are not tested (they are listed
#'   in attribute `"ambient_dominated"`). Their native signal is too small
#'   for the ambient estimate's error.
#' @param lfc test against `|log2 fold change| <= lfc` instead of zero (as
#'   limma's TREAT). In cohorts of a hundred or more samples per group,
#'   residual biases well under 1% become statistically detectable; a small
#'   threshold such as `log2(1.1)` asks for effects that matter.
#' @param huber_k Huber tuning constant, in Pearson-residual units. The
#'   default 2 kept null genes calibrated with and without outliers in
#'   [gate_ambient_de_cohort()]; smaller values over-trim clean data and
#'   inflate false positives, larger ones let outliers through.
#' @param robust_iter number of re-weighting rounds.
#' @param maxit,tol Fisher-scoring controls.
#' @return data.frame per tested gene: logFC (native-expression log2 fold
#'   change), F, PValue, FDR, ambient_frac (mean fraction of the gene's counts
#'   attributed to ambient), and when `naive = TRUE` naive_PValue, naive_FDR
#'   and `ambient_driven` (naive FDR < 0.05 but contamination-aware FDR >=
#'   0.05). Attribute `"samples"`: per-sample ambient fraction and library
#'   factors. Attribute `"ambient_dominated"`: genes skipped as mostly soup.
#' @export
ambient_de <- function(counts, empties = NULL, design, coef = ncol(design),
                       soup = NULL, min_count = 10, naive = TRUE, robust = TRUE, joint_rho = TRUE,
                       soup_design = TRUE, lfc = 0, max_ambient_frac = 0.75,
                       huber_k = 2, robust_iter = 2, maxit = 50, tol = 1e-8) {
  for (p in c("edgeR", "limma"))
    if (!requireNamespace(p, quietly = TRUE)) stop("ambient_de() needs the ", p, " package", call. = FALSE)
  Y <- as.matrix(counts); storage.mode(Y) <- "double"
  S <- ncol(Y)
  design <- as.matrix(design)
  if (nrow(design) != S) stop("nrow(design) must equal ncol(counts)", call. = FALSE)
  if (S - ncol(design) < 1) stop("no residual degrees of freedom", call. = FALSE)
  if (is.character(coef)) coef <- match(coef, colnames(design))
  if (is.na(coef) || coef < 1 || coef > ncol(design)) stop("`coef` not found in design", call. = FALSE)
  if (is.null(soup)) {
    if (is.null(empties) || length(empties) != S)
      stop("supply one empties matrix per sample (or `soup`)", call. = FALSE)
    soup <- soup_eb(empties, if (soup_design) design else NULL)
  }
  if (!identical(dim(soup$soup), dim(Y))) stop("soup and counts dimensions differ", call. = FALSE)

  ## ---- ambient load per sample, self-calibrated on the target pseudobulk ----
  T_s <- colSums(Y)
  rho <- vapply(seq_len(S), function(s)
    as.numeric(estimate_rho_cluster(Y[, s, drop = FALSE], soup$soup[, s], clusters = 1L)), 1)
  rho_var <- rep(0, S)
  if (joint_rho) {
    jr <- .joint_rho(Y, soup$soup, rho)
    rho <- jr$rho; rho_var <- jr$var
  }
  E <- sweep(soup$soup, 2, rho * T_s, "*")
  ## Var(E) = (rho T)^2 Var(soup) + (T soup)^2 Var(rho)
  VE <- sweep(soup$var, 2, (rho * T_s)^2, "*") + sweep(soup$soup^2, 2, T_s^2 * rho_var, "*")
  native <- pmax(Y - E, 0)

  keep <- edgeR::filterByExpr(native, design = design, min.count = min_count)
  ## A gene whose counts are mostly soup has too little native signal: a few
  ## percent error in E is a large error in what is left. It is not tested.
  dominated <- keep & rowMeans(E / pmax(Y, 1)) > max_ambient_frac
  keep <- keep & !dominated
  if (sum(keep) < 2) stop("fewer than 2 genes pass filtering", call. = FALSE)
  nf <- .norm_lib_sizes(native[keep, , drop = FALSE])
  N <- (1 - rho) * T_s * nf

  ## ---- NB dispersion from the soup-removed counts ----
  dge <- edgeR::DGEList(round(native[keep, , drop = FALSE]), norm.factors = nf)
  phi <- edgeR::estimateDisp(dge, design)$tagwise.dispersion

  y <- Y[keep, , drop = FALSE]; e <- E[keep, , drop = FALSE]; ve <- VE[keep, , drop = FALSE]
  full <- .fit_ambient_nb(y, e, ve, N, design, phi, maxit, tol)
  w <- NULL
  if (robust) {
    ## Outlier samples: a few extreme counts shift a group mean while the NB
    ## dispersion (a likelihood, not a moment, estimate) understates their
    ## variance, so NB tests call the gene (Li et al. 2022). Down-weight
    ## observations by Huber weights on the full model's Pearson residuals,
    ## re-estimate the dispersion with those weights, and iterate (the scheme
    ## of edgeR::estimateGLMRobustDisp). Weights come from the FULL model, so
    ## a real condition effect is fitted, not down-weighted.
    for (i in seq_len(robust_iter)) {
      w <- pmin(huber_k / pmax(abs(full$pearson), 1e-8), 1)
      dge$weights <- w
      phi <- edgeR::estimateDisp(dge, design)$tagwise.dispersion
      full <- .fit_ambient_nb(y, e, ve, N, design, phi, maxit, tol, w)
    }
  }
  red <- .fit_ambient_nb(y, e, ve, N, design[, -coef, drop = FALSE], phi, maxit, tol, w)

  ## ---- quasi-likelihood F-test with EB-moderated residual variance ----
  ## down-weighting removes information: the residual df is the total weight
  ## minus the parameters, per gene
  df_res <- if (is.null(w)) rep(S - ncol(design), nrow(y)) else pmax(rowSums(w) - ncol(design), 1)
  s2 <- pmax(full$deviance / df_res, 1e-8)
  sq <- limma::squeezeVar(s2, df_res, covariate = log(rowMeans(native[keep, , drop = FALSE]) + 1))
  lr <- pmax(red$deviance - full$deviance, 0)
  Fstat <- lr / sq$var.post
  df_tot <- df_res + sq$df.prior
  pval <- stats::pf(Fstat, 1, df_tot, lower.tail = FALSE)
  if (lfc > 0) {
    ## TREAT-style interval null, with the standard error implied by F
    b <- abs(full$beta[, coef]) / log(2)
    se <- b / sqrt(pmax(Fstat, 1e-12))
    pval <- stats::pt((b - lfc) / se, df_tot, lower.tail = FALSE) +
      stats::pt((b + lfc) / se, df_tot, lower.tail = FALSE)
    pval <- pmin(pval, 1)
  }

  out <- data.frame(gene = rownames(Y)[keep], logFC = full$beta[, coef] / log(2),
                    F = Fstat, PValue = pval, FDR = stats::p.adjust(pval, "BH"),
                    ambient_frac = rowMeans(e / pmax(y, 1)),
                    stringsAsFactors = FALSE, row.names = NULL)
  if (naive) {
    d0 <- .norm_lib_sizes(edgeR::DGEList(Y[keep, , drop = FALSE]))
    d0 <- edgeR::estimateDisp(d0, design)
    t0 <- edgeR::glmQLFTest(edgeR::glmQLFit(d0, design), coef = coef)$table
    out$naive_PValue <- t0$PValue
    out$naive_FDR <- stats::p.adjust(t0$PValue, "BH")
    out$ambient_driven <- out$naive_FDR < 0.05 & out$FDR >= 0.05
  }
  out <- out[order(out$PValue), ]
  attr(out, "ambient_dominated") <- rownames(Y)[dominated]
  attr(out, "samples") <- data.frame(sample = colnames(Y) %||% seq_len(S), ambient_frac = rho,
                                     soup_umis = soup$umis, norm_factor = nf,
                                     stringsAsFactors = FALSE)
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## Joint absent-gene refinement of per-sample rho. The absent set is the genes
## whose counts, pooled over all samples, are not significantly above the
## pooled ambient expectation; each sample's rho is then its own count over
## its own expectation on that set. Pooling gives the power to exclude weakly
## expressed genes that a single sample's test lets through. A sample with
## too few counts on the joint set keeps its per-sample estimate. Returns the
## refined rho and its sampling variance.
.joint_rho <- function(Y, soup, rho, alpha = 0.01, max_iter = 20, min_counts = 50,
                       min_expected = 5) {
  e1 <- sweep(soup, 2, colSums(Y), "*")
  start <- unname(rho); rho <- start; A <- rep(FALSE, nrow(Y))
  none <- list(rho = start, var = rep(0, length(start)))
  for (it in seq_len(max_iter)) {
    m <- as.numeric(e1 %*% rho)
    A <- m >= min_expected & stats::ppois(rowSums(Y) - 1, m, lower.tail = FALSE) > alpha
    if (!any(A)) return(none)
    new <- unname(colSums(Y[A, , drop = FALSE]) / pmax(colSums(e1[A, , drop = FALSE]), 1e-12))
    done <- max(abs(new - rho)) < 1e-6
    rho <- new
    if (done) break
  }
  cnt <- colSums(Y[A, , drop = FALSE]); den <- pmax(colSums(e1[A, , drop = FALSE]), 1e-12)
  v <- unname(cnt / den^2)                        # Poisson variance of count / expectation
  weak <- unname(cnt < min_counts)
  rho[weak] <- start[weak]; v[weak] <- 0
  list(rho = pmin(pmax(rho, 0), 0.99), var = v)
}

## TMM factors. edgeR >= 4.0 renamed calcNormFactors() to normLibSizes() and
## prints a message on every call to the old name; use whichever exists.
.norm_lib_sizes <- function(x, ...) {
  f <- if ("normLibSizes" %in% getNamespaceExports("edgeR")) edgeR::normLibSizes
       else edgeR::calcNormFactors
  f(x, ...)
}

## Vectorised Fisher scoring for y ~ NB(mu = N * exp(X beta) + e), all genes at
## once. Variance V = mu + phi * m^2 + ve (moment-matched NB, with the soup's
## sampling variance ve). `w` are optional genes x samples observation weights.
## Returns beta (genes x p), the (weighted) NB deviance and Pearson residuals.
.fit_ambient_nb <- function(y, e, ve, N, X, phi, maxit, tol, w = NULL) {
  G <- nrow(y); S <- ncol(y); p <- ncol(X)
  if (is.null(w)) w <- matrix(1, G, S)
  logN <- log(N)
  eta_off <- matrix(logN, G, S, byrow = TRUE)
  ## start: log-linear fit to the soup-removed counts
  z0 <- log(pmax(y - e, 0.5)) - eta_off
  beta <- matrix(t(qr.solve(X, t(z0))), G, p)
  native <- function(beta, idx = seq_len(G))
    exp(beta %*% t(X) + eta_off[idx, , drop = FALSE])
  ## NB log-likelihood with the size held FIXED (alternating scheme): every
  ## step and line search then optimise one objective, and the fixed point
  ## solves the quasi-likelihood equations sum_s (y - mu)/V * m * x_s = 0
  ## with V = mu + phi m^2 + ve.
  ll_fixed <- function(beta, size, idx = seq_len(G)) {
    mu <- native(beta, idx) + e[idx, , drop = FALSE]
    rowSums(w[idx, , drop = FALSE] *
              stats::dnbinom(y[idx, , drop = FALSE], size = size[idx, , drop = FALSE], mu = mu, log = TRUE))
  }
  size_at <- function(beta) {
    m <- native(beta); mu <- m + e
    mu^2 / pmax(phi * m^2 + ve, 1e-10)
  }
  for (it in seq_len(maxit)) {
    size <- size_at(beta)
    ll <- ll_fixed(beta, size)
    m <- native(beta); mu <- m + e
    V <- mu + mu^2 / size                       # = mu + phi m^2 + ve
    score <- (w * (y - mu) / V * m) %*% X       # G x p
    iw <- w * m^2 / V
    info <- array(0, c(G, p, p))
    for (j in seq_len(p)) for (k in j:p) {
      info[, j, k] <- info[, k, j] <- iw %*% (X[, j] * X[, k])
    }
    step <- .batched_solve(info, score)
    lam <- rep(1, G); new <- beta + step; ll_new <- ll_fixed(new, size)
    for (h in 1:30) {
      bad <- which(!is.finite(ll_new) | ll_new < ll - 1e-12)
      if (!length(bad)) break
      lam[bad] <- lam[bad] / 2
      new[bad, ] <- beta[bad, , drop = FALSE] + lam[bad] * step[bad, , drop = FALSE]
      ll_new[bad] <- ll_fixed(new[bad, , drop = FALSE], size, bad)
    }
    bad <- !is.finite(ll_new) | ll_new < ll - 1e-12
    new[bad, ] <- beta[bad, ]
    delta <- max(abs(new - beta))
    beta <- new
    if (delta < tol) break
  }
  size <- size_at(beta)
  mu <- native(beta) + e
  sat <- stats::dnbinom(y, size = size, mu = pmax(y, 1e-10), log = TRUE)
  fit <- stats::dnbinom(y, size = size, mu = mu, log = TRUE)
  list(beta = beta, deviance = pmax(2 * rowSums(w * (sat - fit)), 0), loglik = rowSums(w * fit),
       pearson = (y - mu) / sqrt(mu + mu^2 / size), iterations = it)
}

## Solve A[g,,] x = b[g,] for every g at once (A symmetric positive definite,
## small p), by vectorised Gaussian elimination with a small ridge.
.batched_solve <- function(A, b) {
  G <- dim(A)[1]; p <- dim(A)[2]
  for (j in seq_len(p)) A[, j, j] <- A[, j, j] + 1e-8
  for (j in seq_len(p)) {
    piv <- A[, j, j]
    if (j < p) for (i in (j + 1):p) {
      f <- A[, i, j] / piv
      A[, i, ] <- A[, i, , drop = FALSE][, 1, ] - f * A[, j, , drop = FALSE][, 1, ]
      b[, i] <- b[, i] - f * b[, j]
    }
  }
  x <- matrix(0, G, p)
  for (j in p:1) {
    s <- b[, j]
    if (j < p) for (k in (j + 1):p) s <- s - A[, j, k] * x[, k]
    x[, j] <- s / A[, j, j]
  }
  x
}
