## The one property that survived every experiment: correction can never
## fabricate counts. Checked for every rule, dense and sparse.
test_that("no correction rule ever adds or negates counts", {
  sim  <- small_sim()
  soup <- ambient_global(sim$empty)
  cl   <- quick_labels(sim$observed, k = 4)
  rho  <- as.numeric(estimate_rho_cluster(sim$observed, soup, cl))
  for (m in c("subtract", "redistribute", "bayes")) for (sp in c(FALSE, TRUE)) {
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
