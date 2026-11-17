test_that("cap_test = 'cluster' removes rare ambient singletons but keeps doublet-level counts", {
  set.seed(1)
  G <- 20; N <- 10000
  ## gene 1 is rare in the soup: about 0.0005 expected ambient counts per cell
  ## (below cap_q, where the per-entry test sets every count aside), about 5
  ## over the 10000 cells, so a few ambient singletons are expected
  soup <- c(2.5e-5, rep((1 - 2.5e-5) / (G - 1), G - 1))
  X <- matrix(stats::rpois(G * N, 20), G, N)
  X[1, ] <- 0                                         # the cells do not express gene 1
  X[1, 1:4] <- 1                                      # ambient singletons
  X[1, 6] <- 40                                       # far above ambient: a doublet
  rho <- rep(0.05, N); cl <- rep("a", N)
  e <- correct_counts(X, soup, rho, method = "bayes_cell", clusters = cl, cap_test = "entry")
  k <- correct_counts(X, soup, rho, method = "bayes_cell", clusters = cl, cap_test = "cluster")
  ## per-entry test: each singleton is "too high to be ambient" and kept
  expect_true(all(e[1, 1:4] > 0.99))
  ## corrected for the cells it is applied to: the singletons are ambient
  expect_true(all(k[1, 1:4] < 0.01))
  ## the doublet-level count keeps all but its ambient expectation under both
  lam <- rho[6] * sum(X[, 6]) * soup[1]
  expect_equal(k[1, 6], 40 - lam, tolerance = 1e-8)
  expect_equal(e[1, 6], 40 - lam, tolerance = 1e-8)
  ## never adds counts
  expect_true(all(k <= X + 1e-12) && all(k >= 0))
  ## the default is unchanged
  expect_equal(correct_counts(X, soup, rho, method = "bayes_cell", clusters = cl), e)
})

test_that("bayes_absent removes a gene a cluster does not express and keeps one it does", {
  set.seed(2)
  G <- 400; N <- 600; cl <- rep(c("a", "b"), each = N / 2)
  ## genes 1-200 expressed by cluster a only, 201-400 by cluster b only
  prof <- cbind(a = c(rep(1, 200), rep(0, 200)), b = c(rep(0, 200), rep(1, 200))) * stats::rexp(G, 1)
  prof <- sweep(prof, 2, colSums(prof), "/")
  soup <- rowMeans(prof); rho <- 0.08; depth <- 4000
  scatter <- exp(stats::rnorm(G, 0, 0.5))                       # ambient levels scatter gene to gene
  amb <- soup * scatter / sum(soup * scatter)
  X <- sapply(seq_len(N), function(j) stats::rpois(G, depth * ((1 - rho) * prof[, cl[j]] + rho * amb)))
  M <- 2e6
  empt_soup <- as.numeric(stats::rmultinom(1, M, soup)) / M     # soup as estimated from empties
  fa <- 201:400; ja <- cl == "a"                                # genes foreign to cluster a
  cell <- correct_counts(X, empt_soup, rho, method = "bayes_cell", clusters = cl)
  absn <- correct_counts(X, empt_soup, rho, method = "bayes_absent", clusters = cl, soup_umis = M)
  rem <- function(Y) 1 - sum(Y[fa, ja]) / sum(X[fa, ja])
  expect_gt(rem(absn), 0.95)                                    # foreign counts removed
  expect_gt(rem(absn), rem(cell) + 0.1)                         # more than the expectation-based rule
  expect_gt(sum(absn[1:200, ja]) / sum(X[1:200, ja]), 0.9)      # the cluster's own genes kept
  expect_true(all(absn <= X + 1e-9) && all(absn >= 0))          # never adds counts
  ## a doublet-level count of a foreign gene keeps all but its ambient share
  X2 <- X; X2[250, 1] <- 60
  ab2 <- correct_counts(X2, empt_soup, rho, method = "bayes_absent", clusters = cl, soup_umis = M)
  expect_gt(ab2[250, 1], 55)
})
