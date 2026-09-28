test_that("new simulator knobs leave default output unchanged (reproducibility)", {
  ## reference values computed with the v0.1 simulator code
  s <- simulate_experiment(n_cells = 200, n_empty = 300, seed = 7)
  expect_equal(c(sum(s$observed), sum(s$truth), sum(s$empty)), c(620425, 497229, 7396))
  m <- simulate_multimodal(n_cells = 200, n_samples = 2, empty_per_sample = c(100, 50), seed = 3)
  expect_equal(c(sum(m$obs_unspliced), sum(m$truth_spliced)), c(387769, 176162))
})

test_that("DE benchmark is calibrated on null genes", {
  sim <- simulate_de(seed = 1)
  p <- de_naive(sim)
  expect_lt(mean(p[sim$null_idx] < 0.05, na.rm = TRUE), 0.10)
})
