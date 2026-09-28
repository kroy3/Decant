test_that("gene_species parses Cell Ranger multi-genome prefixes", {
  expect_equal(gene_species(c("GRCh38_ENSG1", "mm10___ENSMUSG1")), c("GRCh38", "mm10"))
  expect_error(gene_species(c("A_1", "B_1", "C_1")), "two genome")
})

test_that("species_mix_truth recovers the known rho from cross-species reads", {
  s <- sim_species_mix()
  tr <- species_mix_truth(s$cells, s$empties, gene_species(rownames(s$cells)))
  expect_equal(tr$species, s$species)
  expect_false(any(tr$doublet))
  ## unbiased and tightly correlated with the simulated truth
  expect_lt(abs(mean(tr$rho_true) - mean(s$rho)), 0.01)
  expect_gt(cor(tr$rho_true, s$rho), 0.8)
})

test_that("hidden-species truth matches contamination among the kept genes", {
  s <- sim_species_mix(rho_mean = 0.15)
  sp <- gene_species(rownames(s$cells))
  h <- species_mix_hidden(s$cells, s$empties, sp)
  expect_true(all(gene_species(c(rownames(h$cells), "mm10___x"))[seq_len(nrow(h$cells))] == h$species))
  ## decant on the hidden view should land near the hidden truth
  r <- suppressWarnings(decant(h$cells, h$empties, k = 3))$rho
  expect_lt(abs(mean(r) - mean(h$rho_true)), 0.25 * mean(h$rho_true))
})

test_that("score_species_mix: perfect correction removes all cross-species counts", {
  s <- sim_species_mix()
  sp <- gene_species(rownames(s$cells))
  tr <- species_mix_truth(s$cells, s$empties, sp)
  ## oracle-ish correction: zero out every other-species count
  corr <- s$cells
  for (c in seq_len(ncol(corr))) corr[sp != tr$species[c], c] <- 0
  sc <- score_species_mix(s$cells, s$empties, sp, corr, rho_hat = tr$rho_true, truth = tr)
  expect_equal(sc$cross_removed, 1)
  expect_equal(sc$own_removed_ratio, 0)      # removed nothing of its own species
  expect_equal(sc$rho_rmse, 0)
  none <- score_species_mix(s$cells, s$empties, sp, s$cells, truth = tr)
  expect_equal(none$cross_removed, 0)
})
