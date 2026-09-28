test_that("sparse and dense inputs give identical results", {
  sim  <- small_sim()
  sp   <- as_sparse(sim$observed)
  soup <- ambient_global(as_sparse(sim$empty))
  expect_equal(soup, ambient_global(sim$empty))
  cl <- quick_labels(sim$observed, k = 4)
  expect_equal(quick_labels(sp, k = 4), cl)
  r_d <- estimate_rho_cluster(sim$observed, soup, cl)
  r_s <- estimate_rho_cluster(sp, soup, cl)
  expect_equal(as.numeric(r_d), as.numeric(r_s))
  for (m in c("subtract", "redistribute", "bayes")) {
    cd <- correct_counts(sim$observed, soup, as.numeric(r_d), method = m, clusters = cl)
    cs <- correct_counts(sp, soup, as.numeric(r_d), method = m, clusters = cl)
    expect_s4_class(cs, "dgCMatrix")
    expect_equal(as.matrix(cs), cd, ignore_attr = TRUE, info = m)
  }
})

test_that("decant keeps sparse input sparse", {
  sim <- small_sim()
  res <- decant(as_sparse(sim$observed), as_sparse(sim$empty), k = 4)
  expect_s4_class(res$corrected, "dgCMatrix")
  res_d <- decant(sim$observed, sim$empty, k = 4)
  expect_equal(as.matrix(res$corrected), res_d$corrected, ignore_attr = TRUE)
})

test_that("splice-aware correction: sparse == dense, and never fabricates", {
  sim <- simulate_multimodal(n_genes = 200, n_types = 4, n_cells = 150, n_samples = 2,
                             empty_per_sample = c(300, 300), seed = 5)
  eu <- do.call(cbind, lapply(sim$empties, `[[`, "unspliced"))
  es <- do.call(cbind, lapply(sim$empties, `[[`, "spliced"))
  rho <- rep(0.2, ncol(sim$obs_total))
  d <- correct_splice_aware(sim$obs_unspliced, sim$obs_spliced, eu, es, rho)
  s <- correct_splice_aware(as_sparse(sim$obs_unspliced), as_sparse(sim$obs_spliced),
                            as_sparse(eu), as_sparse(es), rho)
  expect_equal(as.matrix(s$unspliced), d$unspliced, ignore_attr = TRUE)
  expect_equal(as.matrix(s$spliced), d$spliced, ignore_attr = TRUE)
  expect_true(all(d$unspliced >= 0 & d$unspliced <= sim$obs_unspliced))
  expect_true(all(d$spliced >= 0 & d$spliced <= sim$obs_spliced))
})
