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
