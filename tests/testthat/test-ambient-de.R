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
