suppressMessages({ library(edgeR) })

pb <- function(M, sample, keep) {                    # pseudobulk: genes x samples
  S <- sort(unique(sample))
  sapply(S, function(s) Matrix::rowSums(M[, keep & sample == s, drop = FALSE]))
}

edger_test <- function(Y, cond, genes) {
  d <- DGEList(Y[genes, ]); d <- calcNormFactors(d)
  design <- model.matrix(~cond)
  d <- estimateDisp(d, design)
  fit <- glmQLFit(d, design)
  res <- glmQLFTest(fit, coef = 2)
  p <- res$table$PValue; names(p) <- genes; p
}

## Contamination-aware NB: mu_s = N_s * exp(b0 + b1 x_s) + E_s, E_s known.
## Dispersion borrowed from edgeR on the corrected counts (trended+tagwise).
caware_test <- function(Y, E, N, cond, genes, disp) {
  x <- as.numeric(cond) - 1
  nll <- function(par, y, e, n, phi, b1fixed = NULL) {
    b1 <- if (is.null(b1fixed)) par[2] else b1fixed
    mu <- n * exp(par[1] + b1 * x) + e
    -sum(dnbinom(y, size = 1 / phi, mu = mu, log = TRUE))
  }
  p <- vapply(seq_along(genes), function(i) {
    g <- genes[i]; y <- Y[g, ]; e <- E[g, ]; n <- N; phi <- disp[i]
    start <- log(max(sum(pmax(y - e, 0.5)), 1) / sum(n))
    f1 <- optim(c(start, 0), nll, y = y, e = e, n = n, phi = phi, method = "BFGS")
    f0 <- optim(start, function(b) nll(b, y, e, n, phi, b1fixed = 0), method = "BFGS")
    stat <- max(2 * (f0$value - f1$value), 0)
    stats::pchisq(stat, 1, lower.tail = FALSE)
  }, numeric(1))
  names(p) <- genes; p
}

run_all <- function(st, target = "R", rho_scale = 1, clusters = c("true", "kmeans")) {
  clusters <- match.arg(clusters)
  pkg <- if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(Sys.getenv("DECANT_ROOT", "../.."), quiet = TRUE)
  keepR <- st$label == target
  S <- sort(unique(st$sample))
  Y0 <- pb(st$X, st$sample, keepR)
  ## per-sample decant: rho per cell, corrected counts, soup per sample
  corr <- st$X; rho <- numeric(ncol(st$X)); soup <- matrix(0, nrow(st$X), length(S))
  for (s in S) {
    idx <- which(st$sample == s)
    cl <- if (clusters == "true") st$label[idx] else quick_labels(st$X[, idx], k = 4, seed = s)
    r <- suppressWarnings(decant(st$X[, idx], st$empties[[s]], clusters = cl))
    corr[, idx] <- as.matrix(r$corrected); rho[idx] <- pmin(r$rho * rho_scale, 0.99); soup[, s] <- r$ambient[, 1]
  }
  Y1 <- pb(corr, st$sample, keepR)
  T_c <- colSums(st$X)
  amb_load <- sapply(S, function(s) sum((rho * T_c)[keepR & st$sample == s]))   # expected ambient UMIs
  E <- sweep(soup, 2, amb_load, "*"); rownames(E) <- rownames(st$X)
  N <- sapply(S, function(s) sum(((1 - rho) * T_c)[keepR & st$sample == s]))
  N <- N / mean(N)
  genes <- rownames(Y0)[filterByExpr(DGEList(Y0), model.matrix(~st$cond))]
  ## dispersion for caware: edgeR on corrected counts
  d1 <- estimateDisp(calcNormFactors(DGEList(Y1[genes, ])), model.matrix(~st$cond))
  list(genes = genes,
       naive = edger_test(Y0, st$cond, genes),
       corrected = edger_test(Y1, st$cond, genes),
       caware = caware_test(Y0, E, N, st$cond, genes, d1$tagwise.dispersion),
       caware2 = caware2_test(Y0, pool_soup(st$empties)$pooled,
                              vapply(st$empties, sum, 1), st$cond, genes, d1$tagwise.dispersion),
       caware3 = caware3_test(Y0, st$empties, st$cond, genes, d1$tagwise.dispersion),
       caware4 = caware4_test(Y0, st$empties, st$cond, genes, d1$tagwise.dispersion),
       caware2_nopool = caware2_test(Y0, pool_soup(st$empties)$independent,
                              vapply(st$empties, sum, 1), st$cond, genes, d1$tagwise.dispersion),
       rho_est = tapply(rho[keepR], st$sample[keepR], mean))
}

score <- function(st, r, alpha = 0.05, fdr = 0.1) {
  gi <- as.integer(sub("g", "", r$genes))
  is_de <- gi %in% st$de_rare; is_leak <- gi %in% st$de_abund
  is_sh <- gi %in% st$shared
  do.call(rbind, lapply(c("naive", "corrected", "caware", "caware3", "caware4"), function(m) {
    p <- r[[m]]; q <- p.adjust(p, "BH")
    data.frame(method = m,
               fp_leak = mean(p[is_leak] < alpha),         # A's disease genes called DE in R
               fp_other_null = mean(p[!is_de & !is_leak] < alpha),
               power = mean(p[is_de & !is_sh] < alpha),
               power_shared = if (any(is_sh)) mean(p[is_sh] < alpha) else NA_real_,
               n_disc = sum(q < fdr),
               fdp = if (sum(q < fdr)) mean(!is_de[q < fdr]) else 0,
               leak_tested = sum(is_leak))
  }))
}
