## chirp(): ambient-aware ligand-receptor calls against simulated ground truth.
## The soup is dominated by one cell type; its marker genes are used as
## ligands, so a detection-based caller sees every cluster as a sender.

chirp_setup <- function(rho_mean = 0.15, seed = 7) {
  sim <- simulate_experiment(n_genes = 400, n_types = 5, n_cells = 800, n_empty = 3000,
                             rho_mean = rho_mean, soup_bias = 3, seed = seed)
  prof <- sim$profiles
  dom <- which.max(sim$lysis_true)
  top <- function(k, n) {
    s <- prof[, k] / rowSums(prof)
    names(sort(-(s * prof[, k])))[seq_len(n)]
  }
  others <- setdiff(seq_len(ncol(prof)), dom)
  pairs <- expand.grid(ligand = top(dom, 8),
                       receptor = unlist(lapply(others, top, n = 2)),
                       stringsAsFactors = FALSE)
  list(sim = sim, pairs = pairs, soup = ambient_global(sim$empty))
}

## true expression (per 10k) of the weaker side of each row
chirp_truth <- function(res, prof) {
  e <- prof * 1e4
  g <- rownames(prof)
  pmin(e[cbind(match(res$ligand, g), as.integer(res$sender))],
       e[cbind(match(res$receptor, g), as.integer(res$receiver))])
}

test_that("chirp removes ambient-driven senders that detection calls keep", {
  x <- chirp_setup()
  rho <- as.numeric(estimate_rho_cluster(x$sim$observed, x$soup, x$sim$labels))
  res <- chirp(x$sim$observed, x$pairs, soup = x$soup, rho = rho, clusters = x$sim$labels)
  m <- chirp_truth(res, x$sim$profiles)
  absent <- m < 0.5; strong <- m >= 20
  ## detection rule: calls most absent pairs (ligand is ambient everywhere)
  expect_gt(mean(res$naive_call[absent]), 0.5)
  ## chirp: no absent pair called, every strong pair found
  expect_equal(sum(res$call & absent), 0)
  expect_true(all(res$call[strong]))
  ## native expression tracks the truth for markers in their own cluster
  nat <- attr(res, "native")
  g <- unique(c(x$pairs$ligand, x$pairs$receptor))
  hi <- x$sim$profiles[g, ] * 1e4 > 50
  expect_lt(max(abs(nat[g, ][hi] / (x$sim$profiles[g, ][hi] * 1e4) - 1)), 0.15)
})

test_that("chirp with rho = 0 scores raw expression", {
  x <- chirp_setup()
  res <- chirp(x$sim$observed, x$pairs[1:4, ], soup = x$soup, rho = 0, clusters = x$sim$labels)
  expect_equal(res$score, res$naive_score)
})

test_that("chirp is not fooled when rho is underestimated by a third", {
  x <- chirp_setup()
  rho <- as.numeric(estimate_rho_cluster(x$sim$observed, x$soup, x$sim$labels))
  res <- chirp(x$sim$observed, x$pairs, soup = x$soup, rho = rho * 2 / 3,
               clusters = x$sim$labels)
  expect_equal(sum(res$call & chirp_truth(res, x$sim$profiles) < 0.5), 0)
})

test_that("chirp accepts sparse input, complexes, per-sample soup and a decant fit", {
  x <- chirp_setup()
  sp <- as_sparse(x$sim$observed)
  rho <- as.numeric(estimate_rho_cluster(x$sim$observed, x$soup, x$sim$labels))
  a <- chirp(x$sim$observed, x$pairs[1:3, ], soup = x$soup, rho = rho, clusters = x$sim$labels)
  b <- chirp(sp, x$pairs[1:3, ], soup = x$soup, rho = rho, clusters = x$sim$labels)
  expect_equal(a, b)
  ## same soup in two "samples" gives the same answer
  two <- cbind(s1 = x$soup, s2 = x$soup)
  smp <- rep(c("s1", "s2"), length.out = ncol(sp))
  c2 <- chirp(sp, x$pairs[1:3, ], soup = two, rho = rho, clusters = x$sim$labels,
              sample_of = smp)
  expect_equal(a$score, c2$score)
  ## complex: called only if every subunit is
  rr <- unique(x$pairs$receptor)[c(1, 3)]   # markers of two different types
  cx <- data.frame(ligand = x$pairs$ligand[1], receptor = paste(rr, collapse = "+"))
  r <- chirp(sp, cx, soup = x$soup, rho = rho, clusters = x$sim$labels)
  one <- function(rc) {
    z <- chirp(sp, data.frame(ligand = cx$ligand, receptor = rc), soup = x$soup,
               rho = rho, clusters = x$sim$labels)
    z[order(z$sender, z$receiver), ]
  }
  r <- r[order(r$sender, r$receiver), ]
  u <- one(rr[1]); v <- one(rr[2])
  expect_equal(r$call, u$call & v$call)
  expect_equal(r$score, pmin(u$score, v$score))
  expect_lt(sum(r$call), sum(u$call))
  expect_warning(chirp(sp, data.frame(ligand = "nope", receptor = "gene1"),
                       soup = x$soup, rho = rho, clusters = x$sim$labels), "skipped")
  fit <- decant(sp, x$sim$empty, clusters = x$sim$labels)
  f <- chirp(sp, x$pairs[1:3, ], fit = fit)
  expect_equal(nrow(f), 3 * 5 * 5)
})
