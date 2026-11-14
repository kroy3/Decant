skip_if_not_installed("edgeR")
skip_if_not_installed("limma")

pb_R <- function(st) {
  keep <- st$label == "R"; S <- sort(unique(st$sample))
  Y <- vapply(S, function(i) rowSums(st$X[, keep & st$sample == i, drop = FALSE]), numeric(nrow(st$X)))
  rownames(Y) <- rownames(st$X); Y
}

test_that("ambient_de stops leakage that naive edgeR reports, and stays calibrated", {
  st <- simulate_ambient_study(seed = 3, G = 600, n_cells = 500)
  res <- ambient_de(pb_R(st), st$empties, model.matrix(~st$cond))
  gi <- as.integer(sub("g", "", res$gene))
  leak <- gi %in% st$de_abund; de <- gi %in% st$de_rare
  expect_equal(sum(res$FDR[leak] < 0.1), 0)                     # nothing leaked
  expect_lt(mean(res$PValue[!leak & !de] < 0.05), 0.1)           # calibrated nulls
  expect_gt(mean(res$FDR[de] < 0.1), 0.6)                        # real signal found
  expect_true(any(res$ambient_driven))                           # naive hits flagged
  expect_true(all(c("logFC", "F", "PValue", "FDR", "ambient_frac") %in% names(res)))
  s <- attr(res, "samples")
  expect_true(mean(s$ambient_frac[st$cond == "dis"]) > mean(s$ambient_frac[st$cond == "ctrl"]))
})

test_that("soup_eb keeps well-measured samples and shrinks noisy ones", {
  set.seed(1)
  base <- rgamma(200, 0.5); base <- base / sum(base)
  rich <- lapply(1:3, function(i) rmultinom(5000, 20, base))
  poor <- list(rmultinom(5, 20, base))
  sb <- soup_eb(c(rich, poor))
  expect_equal(dim(sb$soup), c(200, 4))
  expect_equal(unname(colSums(sb$soup)), rep(1, 4))
  err <- colSums(abs(sb$soup - base))
  raw_poor <- rowSums(poor[[1]]) / sum(poor[[1]])
  expect_lt(err[4], sum(abs(raw_poor - base)))                   # shrinkage helped the poor one
})

test_that("Fisher scoring solves the quasi-likelihood equations, close to full ML", {
  set.seed(2)
  G <- 30; S <- 8; X <- cbind(1, rep(0:1, each = 4)); N <- runif(S, 0.8, 1.2)
  e <- matrix(rpois(G * S, 20), G); y <- e + matrix(rpois(G * S, 50), G)
  phi <- rep(0.05, G); ve <- matrix(0, G, S)
  f <- .fit_ambient_nb(y, e, ve, N, X, phi, 100, 1e-10)
  ## quasi-score sum_s (y - mu) / V * m * x_s is zero at the solution
  m <- sweep(exp(f$beta %*% t(X)), 2, N, "*"); mu <- m + e; V <- mu + phi * m^2
  score <- ((y - mu) / V * m) %*% X
  expect_lt(max(abs(score)), 1e-4)
  ## and the solution is within a whisker of the full-likelihood optimum
  nll <- function(b, i) { m <- N * exp(drop(X %*% b)); mu <- m + e[i, ]
    -sum(dnbinom(y[i, ], size = mu^2 / (phi[i] * m^2), mu = mu, log = TRUE)) }
  for (i in 1:5) {
    o <- optim(f$beta[i, ] + 0.3, nll, i = i, method = "BFGS", control = list(reltol = 1e-12))
    expect_lt(abs(-o$value - f$loglik[i]), 0.02)
    expect_lt(max(abs(o$par - f$beta[i, ])), 0.05)
  }
})

test_that("input validation", {
  st <- simulate_ambient_study(seed = 1, G = 300, n_cells = 200, n_per = 3)
  Y <- pb_R(st); d <- model.matrix(~st$cond)
  expect_error(ambient_de(Y, st$empties[-1], d), "one empties matrix per sample")
  expect_error(ambient_de(Y, st$empties, d[-1, ]), "nrow\\(design\\)")
  expect_error(ambient_de(Y, st$empties, d, coef = "nope"), "coef")
})

test_that("observation weights enter the fit: a down-weighted outlier stops moving beta", {
  set.seed(4)
  G <- 20; S <- 40; X <- cbind(1, rep(0:1, each = 20)); N <- rep(1, S)
  e <- matrix(5, G, S); y <- e + matrix(rnbinom(G * S, mu = 100, size = 20), G)
  phi <- rep(0.05, G); ve <- matrix(0, G, S)
  clean <- .fit_ambient_nb(y, e, ve, N, X, phi, 100, 1e-10)
  y2 <- y; y2[, S] <- y2[, S] * 20                     # one extreme sample in group 2
  hit <- .fit_ambient_nb(y2, e, ve, N, X, phi, 100, 1e-10)
  w <- matrix(1, G, S); w[, S] <- 1e-6
  wtd <- .fit_ambient_nb(y2, e, ve, N, X, phi, 100, 1e-10, w)
  expect_gt(median(abs(hit$beta[, 2] - clean$beta[, 2])), 0.3)
  expect_lt(median(abs(wtd$beta[, 2] - clean$beta[, 2])), 0.05)
  expect_equal(dim(wtd$pearson), c(G, S))
})

test_that("robust weighting keeps a cohort with outlier samples calibrated", {
  st <- simulate_ambient_cohort(n_per = 25, G = 800, out_p = 0.01, seed = 2)
  d <- model.matrix(~st$cond)
  null <- function(res) {
    gi <- as.integer(sub("g", "", res$gene))
    mean(res$PValue[!gi %in% c(st$de_rare, st$de_abund)] < 0.05)
  }
  rb <- suppressWarnings(ambient_de(st$Y, st$empties, d, naive = FALSE))
  nr <- suppressWarnings(ambient_de(st$Y, st$empties, d, naive = FALSE, robust = FALSE))
  expect_lt(null(rb), 0.085)
  expect_gt(null(nr), 0.15)                                      # the failure being guarded
  gi <- as.integer(sub("g", "", rb$gene))
  expect_equal(sum(rb$FDR[gi %in% st$de_abund] < 0.1), 0)
})

test_that("simulate_ambient_cohort returns a consistent study", {
  st <- simulate_ambient_cohort(n_per = 4, G = 300, seed = 1)
  expect_equal(dim(st$Y), c(300, 8))
  expect_length(st$empties, 8)
  expect_equal(nrow(st$empties[[1]]), 300)
  expect_true(mean(st$rho_true[st$cond == "dis"]) > mean(st$rho_true[st$cond == "ctrl"]))
  expect_error(simulate_ambient_cohort(frac = c(A = 1)), "must name types")
})

test_that("joint refinement removes the upward bias of per-sample rho", {
  st <- simulate_ambient_cohort(n_per = 10, G = 800, seed = 3)
  sp <- soup_eb(st$empties)
  r0 <- vapply(seq_len(ncol(st$Y)), function(s)
    as.numeric(estimate_rho_cluster(st$Y[, s, drop = FALSE], sp$soup[, s], clusters = 1L)), 1)
  r1 <- .joint_rho(st$Y, sp$soup, r0)$rho
  expect_lt(abs(mean(r1 - st$rho_true)), abs(mean(r0 - st$rho_true)) / 2)
  expect_lt(abs(mean(r1 - st$rho_true)), 0.005)
  ## with no usable joint set, the per-sample estimate is returned unchanged
  expect_equal(.joint_rho(st$Y, sp$soup, r0, min_counts = Inf)$rho, r0)
})

test_that("soup_eb with a design does not shrink away a condition difference", {
  set.seed(5)
  a <- rgamma(300, 0.5); a <- a / sum(a); b <- a; b[1:30] <- b[1:30] * 3; b <- b / sum(b)
  emp <- c(lapply(1:6, function(i) rmultinom(1, 3000, a)), lapply(1:6, function(i) rmultinom(1, 3000, b)))
  X <- cbind(1, rep(0:1, each = 6))
  gap <- function(sb) sum(sb$soup[1:30, 7:12]) / sum(sb$soup[1:30, 1:6])
  truth <- sum(b[1:30]) / sum(a[1:30])
  expect_lt(abs(gap(soup_eb(emp, X)) - truth), abs(gap(soup_eb(emp)) - truth))
  expect_error(soup_eb(emp, X[-1, ]), "number of samples")
})

test_that("lfc threshold tests an interval null", {
  st <- simulate_ambient_cohort(n_per = 12, G = 600, fc_rare = 3, seed = 4)
  d <- model.matrix(~st$cond)
  r0 <- suppressWarnings(ambient_de(st$Y, st$empties, d, naive = FALSE))
  r1 <- suppressWarnings(ambient_de(st$Y, st$empties, d, naive = FALSE, lfc = log2(1.2)))
  r1 <- r1[match(r0$gene, r1$gene), ]
  sig <- r0$PValue < 0.05
  expect_true(all(r1$PValue[sig] >= r0$PValue[sig] - 1e-12))
  de <- as.integer(sub("g", "", r1$gene)) %in% st$de_rare
  expect_gt(mean(r1$FDR[de] < 0.1), 0.6)                         # 3-fold genes survive
  small <- abs(r1$logFC) < log2(1.2)
  expect_equal(sum(r1$FDR[small] < 0.1), 0)                      # sub-threshold effects do not
})

test_that("a confounding batch is adjusted for through the design", {
  st <- simulate_ambient_cohort(n_per = 20, G = 800, n_batch = 200, seed = 6)
  fp <- function(d) {
    r <- suppressWarnings(ambient_de(st$Y, st$empties, d, naive = FALSE))
    mean(r$PValue[as.integer(sub("g", "", r$gene)) %in% st$batch_genes] < 0.05)
  }
  adjusted <- fp(model.matrix(~st$batch + st$cond))
  unadjusted <- fp(model.matrix(~st$cond))
  expect_lt(adjusted, 0.15)
  expect_gt(unadjusted, 2 * adjusted)                            # the confounding is real
  expect_equal(length(st$batch), ncol(st$Y))
})

test_that("ambient-dominated genes are skipped and reported", {
  st <- simulate_ambient_cohort(n_per = 8, G = 600, seed = 7)
  d <- model.matrix(~st$cond)
  r <- suppressWarnings(ambient_de(st$Y, st$empties, d, naive = FALSE))
  expect_true(all(r$ambient_frac <= 0.75))
  all_in <- suppressWarnings(ambient_de(st$Y, st$empties, d, naive = FALSE, max_ambient_frac = 1))
  expect_length(attr(all_in, "ambient_dominated"), 0)
  expect_setequal(setdiff(all_in$gene, r$gene), attr(r, "ambient_dominated"))
})
