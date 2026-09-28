test_that("cluster rho estimator is calibrated on ground truth", {
  for (rm in c(0.05, 0.2)) {
    sim  <- small_sim(rho_mean = rm)
    soup <- ambient_global(sim$empty)
    r <- estimate_rho_cluster(sim$observed, soup, quick_labels(sim$observed, k = 4))
    expect_lt(abs(mean(r) - mean(sim$rho_true)), 0.25 * mean(sim$rho_true))
    ## per-cell resolution is limited when a cluster has few absent genes (here
    ## the soup-dominating type at rho = 0.2 borrows its prior), so the bar for
    ## per-cell correlation is lower than for the cluster-level mean.
    expect_gt(cor(as.numeric(r), sim$rho_true), 0.5)
    d <- attr(r, "clusters")
    expect_setequal(names(d), c("cluster", "n_cells", "rho", "n_absent_genes",
                                "absent_soup_frac", "absent_counts", "low_power",
                                "borrowed", "rho_two_sided"))
    expect_false(attr(r, "soup_mismatch"))
  }
})

test_that("negative control: no contamination -> rho ~ 0", {
  sim <- small_sim(rho_mean = 0)
  r <- estimate_rho_cluster(sim$observed, ambient_global(sim$empty),
                            quick_labels(sim$observed, k = 4))
  expect_lt(mean(r), 0.01)
})

test_that("cluster estimator beats legacy estimator", {
  sim  <- small_sim(rho_mean = 0.05)
  soup <- ambient_global(sim$empty)
  new <- as.numeric(estimate_rho_cluster(sim$observed, soup, quick_labels(sim$observed, k = 4)))
  old <- estimate_rho(sim$observed, soup)
  rmse <- function(x) sqrt(mean((x - sim$rho_true)^2))
  expect_lt(rmse(new), rmse(old))
})

test_that("input validation", {
  sim <- small_sim()
  soup <- ambient_global(sim$empty)
  expect_error(estimate_rho_cluster(sim$observed, soup, 1:3), "clusters")
  expect_error(estimate_rho_cluster(sim$observed, soup[-1], rep(1, 400)), "soup")
})

test_that("junk in the empties does not fail silently", {
  ## Known failure mode: material in the empties that no cell contains makes
  ## rho collapse. The estimator must at least flag it.
  sim <- simulate_experiment(n_genes = 600, n_cells = 800, n_empty = 2000,
                             rho_mean = 0.2, seed = 3)
  set.seed(9)
  G <- nrow(sim$empty)
  junk <- stats::rgamma(G, 0.05); junk[sample(G, 20)] <- junk[sample(G, 20)] + 50
  emp <- sim$empty + vapply(seq_len(ncol(sim$empty)), function(i)
    stats::rmultinom(1, round(0.2 * sum(sim$empty[, i])), junk)[, 1], numeric(G))
  expect_warning(decant(sim$observed, emp, k = 6), "UNDERestimated")
})
