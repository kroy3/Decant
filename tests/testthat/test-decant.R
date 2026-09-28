test_that("single-sample pipeline runs and reports modules", {
  sim <- small_sim()
  res <- decant(sim$observed, sim$empty, k = 4)
  expect_s3_class(res, "decant")
  expect_true(all(c("kmeans_clusters", "cluster_rho", "correction_subtract") %in% res$modules))
  expect_equal(dim(res$corrected), dim(sim$observed))
  expect_length(res$rho, ncol(sim$observed))
  expect_output(print(res), "decant")
  ## correction actually improves on doing nothing
  none <- score_correction(sim$observed, sim$observed, sim$truth)$l1_error
  ours <- score_correction(sim$observed, res$corrected, sim$truth)$l1_error
  expect_lt(ours, none)
})

test_that("user-supplied clusters and rho are respected", {
  sim <- small_sim()
  res <- decant(sim$observed, sim$empty, clusters = sim$labels)
  expect_false("kmeans_clusters" %in% res$modules)
  res2 <- decant(sim$observed, sim$empty, rho = 0.1)
  expect_true("user_rho" %in% res2$modules)
  expect_equal(res2$rho, rep(0.1, ncol(sim$observed)))
})

test_that("multi-sample and splice-aware modes run", {
  sim <- simulate_multimodal(n_genes = 300, n_types = 4, n_cells = 300, n_samples = 3,
                             empty_per_sample = c(800, 800, 100), seed = 2)
  emp_tot <- lapply(sim$empties, function(e) e$unspliced + e$spliced)
  res <- decant(sim$obs_total, emp_tot, sample_of = sim$sample_of, k = 4,
                obs_unspliced = sim$obs_unspliced, obs_spliced = sim$obs_spliced,
                empties_unspliced = lapply(sim$empties, `[[`, "unspliced"),
                empties_spliced = lapply(sim$empties, `[[`, "spliced"))
  expect_true(all(c("hierarchical_ambient", "splice_aware_correction") %in% res$modules))
  expect_true(all(res$corrected_unspliced <= sim$obs_unspliced))
  expect_equal(ncol(res$ambient), 3)

  res_ns <- decant(sim$obs_total, emp_tot, sample_of = sim$sample_of, k = 4)
  expect_true(all(res_ns$corrected <= sim$obs_total))
})

test_that("bad inputs fail loudly", {
  sim <- small_sim()
  expect_error(decant(sim$observed, list(sim$empty)), "sample_of")
  expect_error(decant(sim$observed, sim$empty[-1, ]), "genes")
  expect_error(decant(sim$observed, sim$empty, rho = 2), "rho")
  expect_error(decant(sim$observed, sim$empty, clusters = 1:3), "clusters")
})
