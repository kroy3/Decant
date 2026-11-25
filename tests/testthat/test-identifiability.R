## Objective C: identifiability diagnostics and the classed warning.

test_that("well-separated clusters are not flagged; diagnostics are reported", {
  sim <- small_sim(rho_mean = 0.1)
  soup <- ambient_global(sim$empty)
  r <- estimate_rho_cluster(sim$observed, soup, quick_labels(sim$observed, k = 4))
  d <- attr(r, "clusters")
  for (col in c("cosine", "kappa", "js", "fisher_umi", "n_eff", "rho_se", "ci_width",
                "collinear", "weak", "weak_reason", "action_taken")) expect_true(col %in% names(d), info = col)
  expect_true(all(d$cosine >= 0 & d$cosine <= 1 + 1e-12))
  expect_true(all(d$kappa >= 1))
  expect_true(all(d$ci_width > 0))
  expect_false(any(d$weak))
  fit <- expect_no_warning(decant(as_sparse(sim$observed), as_sparse(sim$empty), k = 4))
  expect_s3_class(fit$identifiability, "data.frame")
  expect_false(any(fit$identifiability$weak))
})

test_that("a soup that resembles one cluster is flagged with a classed warning", {
  ## one type makes up most cells and nearly all of the soup
  sim <- simulate_experiment(n_genes = 600, n_types = 4, rho_mean = 0.1, seed = 3,
                             type_prob = c(0.1, 0.1, 0.1, 1), soup_bias = 10)
  cl <- quick_labels(sim$observed, k = 4, seed = 1)
  w <- NULL
  fit <- withCallingHandlers(
    decant(as_sparse(sim$observed), as_sparse(sim$empty), clusters = cl),
    decant_weak_identifiability = function(c) { w <<- c; invokeRestart("muffleWarning") })
  expect_s3_class(w, "decant_weak_identifiability")
  expect_s3_class(w, "decant_warning")
  idt <- w$identifiability
  weak <- idt[idt$weak, ]
  expect_equal(nrow(weak), 1)
  expect_gt(weak$cosine, 0.9)
  expect_gt(weak$kappa, 30)
  expect_true(weak$collinear)
  expect_true(nzchar(weak$action_taken))
  expect_match(conditionMessage(w), "independent source")
  ## the dominant cluster's own estimate is not used: the cluster-level value
  ## borrowed from the others is close to the truth
  dom <- cl == weak$cluster
  expect_lt(abs(weak$rho - mean(sim$rho_true[dom])), 0.02)
  expect_equal(fit$identifiability, idt)
})

test_that("a collinear cluster that the estimator still pins down is reported, not flagged", {
  sim <- small_sim(rho_mean = 0.1)
  d <- attr(estimate_rho_cluster(sim$observed, ambient_global(sim$empty), sim$labels), "clusters")
  dom <- which.max(d$cosine)
  expect_true(d$collinear[dom]); expect_false(d$weak[dom])
  d2 <- attr(estimate_rho_cluster(sim$observed, ambient_global(sim$empty), sim$labels, weak_rule = "any"), "clusters")
  expect_true(d2$weak[dom])
})

test_that("a single cell type is reported as unidentifiable", {
  sim <- simulate_experiment(n_genes = 600, n_types = 4, rho_mean = 0.1, seed = 3,
                             type_prob = c(0, 0, 0, 1), soup_bias = 10)
  expect_warning(fit <- decant(as_sparse(sim$observed), as_sparse(sim$empty), k = 4),
                 class = "decant_unidentifiable")
  expect_true(fit$rho_unidentifiable)
  expect_true(all(fit$identifiability$weak))
})

test_that("rho's interval is the sampling spread of the absent-gene estimator", {
  ## at a known rho, the reported SE matches the spread across replicate data sets
  est <- se <- numeric(0)
  for (s in 1:12) {
    sim <- small_sim(rho_mean = 0.1, seed = 100 + s)
    soup <- ambient_global(sim$empty)
    d <- attr(estimate_rho_cluster(sim$observed, soup, sim$labels), "clusters")
    est <- c(est, d$rho[2]); se <- c(se, d$rho_se[2])
  }
  ## same order of magnitude (rho also varies with the simulated cells' own rho)
  expect_gt(stats::sd(est) / mean(se), 0.3)
  expect_lt(stats::sd(est) / mean(se), 5)
})
