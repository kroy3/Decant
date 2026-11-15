## The one property that survived every experiment: correction can never
## fabricate counts. Checked for every rule, dense and sparse.
test_that("no correction rule ever adds or negates counts", {
  sim  <- small_sim()
  soup <- ambient_global(sim$empty)
  cl   <- quick_labels(sim$observed, k = 4)
  rho  <- as.numeric(estimate_rho_cluster(sim$observed, soup, cl))
  for (m in c("subtract", "redistribute", "bayes", "bayes_cell")) for (sp in c(FALSE, TRUE)) {
    obs  <- if (sp) as_sparse(sim$observed) else sim$observed
    corr <- as.matrix(correct_counts(obs, soup, rho, method = m, clusters = cl))
    expect_true(all(corr >= 0), info = paste(m, sp))
    expect_true(all(corr <= sim$observed + 1e-9), info = paste(m, sp))
  }
})

test_that("guarantee holds even for absurd rho", {
  sim  <- small_sim()
  soup <- ambient_global(sim$empty)
  for (m in c("subtract", "redistribute")) {
    corr <- correct_counts(sim$observed, soup, rep(0.99, ncol(sim$observed)), method = m)
    expect_true(all(corr >= 0 & corr <= sim$observed))
  }
})

test_that("redistribute removes exactly rho * T when counts allow", {
  sim  <- small_sim(rho_mean = 0.1)
  soup <- ambient_global(sim$empty)
  rho  <- rep(0.1, ncol(sim$observed))
  corr <- correct_counts(sim$observed, soup, rho, method = "redistribute")
  removed <- colSums(sim$observed) - colSums(corr)
  expect_equal(removed, 0.1 * colSums(sim$observed), tolerance = 1e-6)
})

test_that("bayes_cell keeps a count that ambient cannot explain and clears the rest", {
  set.seed(3)
  G <- 200; N <- 300
  soup <- c(rep(20, 5), rep(1, G - 5)); soup <- soup / sum(soup)          # genes 1-5 dominate the soup
  cl <- rep(1:2, each = N / 2)
  prof <- matrix(1, G, 2); prof[1:5, 1] <- 60; prof[1:5, 2] <- 0           # only cluster 1 expresses them
  prof <- sweep(prof, 2, colSums(prof), "/")
  rho <- rep(0.1, N); T_c <- rep(3000, N)
  obs <- vapply(seq_len(N), function(c) stats::rpois(G, T_c[c] * ((1 - rho[c]) * prof[, cl[c]] + rho[c] * soup)),
                numeric(G))
  dbl <- N                                                                  # one doublet in cluster 2
  obs[1:5, dbl] <- obs[1:5, dbl] + 120
  T_c <- colSums(obs)
  old <- correct_counts(obs, soup, rho, "bayes", clusters = cl)
  new <- correct_counts(obs, soup, rho, "bayes_cell", clusters = cl)
  others <- setdiff(which(cl == 2), dbl)
  ## the doublet keeps its foreign counts under bayes_cell, not under bayes
  expect_gt(sum(new[1:5, dbl]) / sum(obs[1:5, dbl]), 0.8)
  expect_lt(sum(old[1:5, dbl]) / sum(obs[1:5, dbl]), 0.6)
  ## and the other cells of the cluster are cleared at least as well
  expect_lt(sum(new[1:5, others]), 0.15 * sum(obs[1:5, others]))
  expect_lte(sum(new[1:5, others]), sum(old[1:5, others]))
  ## native expression in cluster 1 is touched no more than by bayes
  expect_gt(sum(new[1:5, cl == 1]) / sum(obs[1:5, cl == 1]), 0.9)
  expect_error(correct_counts(obs, soup, rho, "bayes_cell"), "needs `clusters`")
})
