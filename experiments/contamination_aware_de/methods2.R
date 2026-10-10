## v2: self-calibrated ambient load + pooled soup + soup sampling variance.
caware2_test <- function(Y, soup, M, cond, genes, disp) {
  ## Y: genes x samples target-type pseudobulk; soup: genes x samples (pooled);
  ## M: total UMIs in each sample's empties (soup precision).
  S <- ncol(Y); T_s <- colSums(Y)
  ## (1) re-estimate ambient fraction per sample from the target pseudobulk itself
  rho_pb <- vapply(seq_len(S), function(s)
    as.numeric(estimate_rho_cluster(Y[, s, drop = FALSE], soup[, s], clusters = 1L)), 1)
  E <- sweep(soup, 2, rho_pb * T_s, "*")
  ## (2) soup sampling variance of E (multinomial over the empties' UMIs)
  VE <- sweep(soup * (1 - soup), 2, (rho_pb * T_s)^2 / M, "*")
  N <- (1 - rho_pb) * T_s; N <- N / mean(N)
  rownames(E) <- rownames(VE) <- rownames(Y)
  x <- as.numeric(cond) - 1
  nll <- function(par, y, e, ve, phi, b1fixed = NULL) {
    b1 <- if (is.null(b1fixed)) par[2] else b1fixed
    mn <- N * exp(par[1] + b1 * x); mu <- mn + e
    size <- mu^2 / pmax(phi * mn^2 + ve, 1e-8)          # moment-matched NB
    -sum(dnbinom(y, size = size, mu = mu, log = TRUE))
  }
  p <- vapply(seq_along(genes), function(i) {
    g <- genes[i]; y <- Y[g, ]; e <- E[g, ]; ve <- VE[g, ]; phi <- disp[i]
    start <- log(max(sum(pmax(y - e, 0.5)), 1) / sum(N))
    f1 <- optim(c(start, 0), nll, y = y, e = e, ve = ve, phi = phi, method = "BFGS")
    f0 <- optim(start, function(b) nll(b, y, e, ve, phi, b1fixed = 0), method = "BFGS")
    stats::pchisq(max(2 * (f0$value - f1$value), 0), 1, lower.tail = FALSE)
  }, numeric(1))
  names(p) <- genes; attr(p, "rho_pb") <- rho_pb; p
}

## Empirical-Bayes soup: per gene, shrink each sample's soup toward the
## cross-sample mean by w = tau2 / (tau2 + v), with v the multinomial sampling
## variance of that sample's estimate and tau2 the between-sample variance of
## the TRUE soups (method of moments). Rich samples keep their own soup; poor
## ones borrow. Returns the posterior mean and posterior variance.
eb_soup <- function(empties) {
  raw <- vapply(empties, function(e) { p <- Matrix::rowSums(e); p / sum(p) }, numeric(nrow(empties[[1]])))
  M <- vapply(empties, sum, 1)
  v <- sweep(raw * (1 - raw), 2, M, "/")                      # sampling variance
  m <- rowMeans(raw)
  tau2 <- pmax(apply(raw, 1, stats::var) - rowMeans(v), 0)     # between-sample true variance
  w <- tau2 / (tau2 + v + 1e-30)
  post <- w * raw + (1 - w) * m
  post <- sweep(post, 2, colSums(post), "/")
  list(soup = post, var = w * v)                               # posterior variance approx
}

caware3_test <- function(Y, empties, cond, genes, disp) {
  eb <- eb_soup(empties)
  S <- ncol(Y); T_s <- colSums(Y)
  rho_pb <- vapply(seq_len(S), function(s)
    as.numeric(estimate_rho_cluster(Y[, s, drop = FALSE], eb$soup[, s], clusters = 1L)), 1)
  E <- sweep(eb$soup, 2, rho_pb * T_s, "*")
  VE <- sweep(eb$var, 2, (rho_pb * T_s)^2, "*")
  N <- (1 - rho_pb) * T_s; N <- N / mean(N)
  rownames(E) <- rownames(VE) <- rownames(Y)
  x <- as.numeric(cond) - 1
  nll <- function(par, y, e, ve, phi, b1fixed = NULL) {
    b1 <- if (is.null(b1fixed)) par[2] else b1fixed
    mn <- N * exp(par[1] + b1 * x); mu <- mn + e
    size <- mu^2 / pmax(phi * mn^2 + ve, 1e-8)
    -sum(dnbinom(y, size = size, mu = mu, log = TRUE))
  }
  p <- vapply(seq_along(genes), function(i) {
    g <- genes[i]; y <- Y[g, ]; e <- E[g, ]; ve <- VE[g, ]; phi <- disp[i]
    start <- log(max(sum(pmax(y - e, 0.5)), 1) / sum(N))
    f1 <- optim(c(start, 0), nll, y = y, e = e, ve = ve, phi = phi, method = "BFGS")
    f0 <- optim(start, function(b) nll(b, y, e, ve, phi, b1fixed = 0), method = "BFGS")
    stats::pchisq(max(2 * (f0$value - f1$value), 0), 1, lower.tail = FALSE)
  }, numeric(1))
  names(p) <- genes; attr(p, "rho_pb") <- rho_pb; p
}


## v4 = v3 + composition normalisation: TMM factors computed on the
## soup-removed native estimate, so strong DE in a few genes does not shift all.
caware4_test <- function(Y, empties, cond, genes, disp) {
  eb <- eb_soup(empties)
  S <- ncol(Y); T_s <- colSums(Y)
  rho_pb <- vapply(seq_len(S), function(s)
    as.numeric(estimate_rho_cluster(Y[, s, drop = FALSE], eb$soup[, s], clusters = 1L)), 1)
  E <- sweep(eb$soup, 2, rho_pb * T_s, "*"); rownames(E) <- rownames(Y)
  VE <- sweep(eb$var, 2, (rho_pb * T_s)^2, "*"); rownames(VE) <- rownames(Y)
  native <- pmax(Y - E, 0)
  nf <- edgeR::calcNormFactors(native[genes, , drop = FALSE])
  N <- (1 - rho_pb) * T_s * nf; N <- N / mean(N)
  x <- as.numeric(cond) - 1
  nll <- function(par, y, e, ve, phi, b1fixed = NULL) {
    b1 <- if (is.null(b1fixed)) par[2] else b1fixed
    mn <- N * exp(par[1] + b1 * x); mu <- mn + e
    size <- mu^2 / pmax(phi * mn^2 + ve, 1e-8)
    -sum(dnbinom(y, size = size, mu = mu, log = TRUE))
  }
  p <- vapply(seq_along(genes), function(i) {
    g <- genes[i]; y <- Y[g, ]; e <- E[g, ]; ve <- VE[g, ]; phi <- disp[i]
    start <- log(max(sum(pmax(y - e, 0.5)), 1) / sum(N))
    f1 <- optim(c(start, 0), nll, y = y, e = e, ve = ve, phi = phi, method = "BFGS")
    f0 <- optim(start, function(b) nll(b, y, e, ve, phi, b1fixed = 0), method = "BFGS")
    stats::pchisq(max(2 * (f0$value - f1$value), 0), 1, lower.tail = FALSE)
  }, numeric(1))
  names(p) <- genes; attr(p, "rho_pb") <- rho_pb; p
}
