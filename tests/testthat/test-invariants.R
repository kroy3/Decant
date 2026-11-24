## Invariant harness (I1, I2): properties that must hold on EVERY input, not
## just on benchmark data. Random negative-binomial matrices over a grid of
## contamination levels, plus the degenerate shapes real data produces
## (all-zero cells, one-cell clusters, clusters that express a single gene).

nb_matrix <- function(G, N, n_types, rho, seed, size = 2, empty_n = 400) {
  set.seed(seed)
  prof <- matrix(stats::rgamma(G * n_types, 0.3), G, n_types)
  prof <- sweep(prof, 2, colSums(prof), "/")
  soup <- as.numeric(prof %*% rep(1 / n_types, n_types))
  lab <- sample(seq_len(n_types), N, replace = TRUE)
  lib <- round(stats::rlnorm(N, log(1500), 0.6))
  mu <- vapply(seq_len(N), function(c) lib[c] * ((1 - rho) * prof[, lab[c]] + rho * soup), numeric(G))
  obs <- matrix(stats::rnbinom(G * N, size = size, mu = mu), G, N)
  empty <- vapply(seq_len(empty_n), function(i) stats::rmultinom(1, 30, soup)[, 1], numeric(G))
  dimnames(obs) <- list(paste0("g", seq_len(G)), paste0("c", seq_len(N)))
  rownames(empty) <- rownames(obs)
  list(obs = obs, empty = empty, soup = soup, lab = lab)
}

all_rules <- c("subtract", "redistribute", "bayes", "bayes_cell", "bayes_absent")

expect_bounded <- function(corr, obs, info) {
  corr <- as.matrix(corr)
  expect_true(all(is.finite(corr)), info = info)
  expect_true(all(corr >= 0), info = info)
  expect_true(all(corr <= obs), info = info)      # exact: no tolerance
}

test_that("every rule stays within [0, observed] on NB data across rho", {
  for (rho in c(0, 0.05, 0.3, 0.6, 0.9)) {
    d <- nb_matrix(150, 120, 3, rho, seed = round(rho * 100) + 1)
    for (m in all_rules) for (sp in c(FALSE, TRUE)) {
      obs <- if (sp) as_sparse(d$obs) else d$obs
      corr <- correct_counts(obs, d$soup, rep(rho, ncol(obs)), method = m, clusters = d$lab,
                             soup_umis = sum(d$empty))
      expect_bounded(corr, d$obs, paste(m, rho, sp))
    }
  }
})

test_that("decant() end to end stays bounded and returns rho in [0, 1)", {
  for (rho in c(0, 0.05, 0.3, 0.6, 0.9)) {
    d <- nb_matrix(150, 120, 3, rho, seed = 10 + round(rho * 100))
    fit <- suppressWarnings(decant(as_sparse(d$obs), as_sparse(d$empty), clusters = d$lab))
    expect_bounded(fit$corrected, d$obs, paste("decant", rho))
    expect_true(all(fit$rho >= 0 & fit$rho < 1), info = paste("rho", rho))
  }
})

test_that("degenerate shapes: all-zero cells, one-cell and single-gene clusters", {
  d <- nb_matrix(120, 60, 3, 0.1, seed = 7)
  obs <- d$obs
  obs[, 1:3] <- 0                              # all-zero cells
  lab <- d$lab
  lab[4] <- 99                                 # a one-cell cluster
  sg <- 5:8; obs[, sg] <- 0; obs[1, sg] <- c(40, 3, 17, 9); lab[sg] <- 77   # a cluster with one gene
  for (m in all_rules) for (sp in c(FALSE, TRUE)) {
    o <- if (sp) as_sparse(obs) else obs
    for (rho in c(0, 0.3, 0.9)) {
      corr <- correct_counts(o, d$soup, rep(rho, ncol(o)), method = m, clusters = lab,
                             soup_umis = sum(d$empty))
      expect_bounded(corr, obs, paste(m, sp, rho))
      expect_true(all(as.matrix(corr)[, 1:3] == 0), info = paste(m, "zero cells"))
    }
  }
  fit <- suppressWarnings(decant(as_sparse(obs), as_sparse(d$empty), clusters = lab))
  expect_bounded(fit$corrected, obs, "decant degenerate")
  expect_true(all(is.finite(fit$rho)) && all(fit$rho >= 0 & fit$rho < 1))
})

test_that("a violation of the bound is an error, never silently clipped", {
  expect_error(Decant:::.assert_bounds(c(1, 2.5), c(1, 2)), class = "decant_invariant_violation")
  expect_error(Decant:::.assert_bounds(c(-0.1, 1), c(1, 2)), class = "decant_invariant_violation")
  expect_error(Decant:::.assert_bounds(c(NaN, 1), c(1, 2)), class = "decant_invariant_violation")
  ## floating-point rounding is not a violation
  expect_equal(Decant:::.assert_bounds(c(2 + 1e-13, -1e-14), c(2, 1)), c(2, 0))
})

## I2: high-count non-ambient signal (doublets, rare expressers) survives the
## DEFAULT correction. Foreign-lineage markers are spiked into a few cells of
## another cluster at doublet level; >= 90% of the spike must be kept.
test_that("default correction keeps doublet-level foreign counts (I2)", {
  sim <- simulate_experiment(n_genes = 600, n_types = 6, n_cells = 1200, n_empty = 3000,
                             rho_mean = 0.1, seed = 11)
  obs <- sim$observed
  mk <- 1:15                                    # the simulator's type-1 markers
  host <- which(sim$labels == 2)[1:3]
  set.seed(1)
  spike <- matrix(stats::rpois(15 * length(host), 10), 15)
  obs[mk, host] <- obs[mk, host] + spike
  ## retention of the SPIKE itself: corrected with it minus corrected without
  ## it, at the same rho (so only the allocation is compared)
  kept_share <- function(rule) {
    f1 <- decant(as_sparse(obs), as_sparse(sim$empty), clusters = sim$labels, correction = rule)
    f0 <- decant(as_sparse(sim$observed), as_sparse(sim$empty), clusters = sim$labels,
                 correction = rule, rho = f1$rho)
    sum(as.matrix(f1$corrected)[mk, host] - as.matrix(f0$corrected)[mk, host]) / sum(spike)
  }
  for (rule in c("bayes_absent", "bayes_cell"))
    expect_gte(kept_share(rule), 0.9, label = paste(rule, "doublet retention"))
  ## the test can fail: the cluster-scaled "bayes" rule loses part of the spike
  expect_lt(kept_share("bayes"), 0.9)
})

## I2: a rare expresser (a gene a cluster does not express, except in two
## cells, at a level far above ambient) keeps its counts in those cells.
test_that("default correction keeps rare high expressers (I2)", {
  sim <- simulate_experiment(n_genes = 600, n_types = 6, n_cells = 1200, n_empty = 3000,
                             rho_mean = 0.15, seed = 12)
  obs <- sim$observed
  g <- 3 * 15 + 1                                # a type-4 marker
  host <- which(sim$labels == 1)[1:2]
  obs[g, host] <- obs[g, host] + 60
  fit <- decant(as_sparse(obs), as_sparse(sim$empty), clusters = sim$labels)
  expect_gte(sum(as.matrix(fit$corrected)[g, host]) / 120, 0.9)
})

## Tiny-cluster degeneracy (seen on real brain data: a 9-cell cluster had
## 17,543 genes called absent). A cluster too small to estimate the ambient
## spread must not have its genes called absent; it keeps the bayes_cell rule.
test_that("tiny clusters do not get a degenerate absence call", {
  sim <- simulate_experiment(n_genes = 4000, n_types = 6, n_cells = 3000, n_empty = 3000,
                             rho_mean = 0.1, seed = 1, overdispersion = 0.3)
  for (ntiny in c(1, 5)) {
    keep <- sim$labels != 6 | cumsum(sim$labels == 6) <= ntiny
    obs <- sim$observed[, keep]; tru <- sim$truth[, keep]; cl <- sim$labels[keep]
    soup <- ambient_global(sim$empty); su <- sum(sim$empty)
    rho <- as.numeric(suppressWarnings(estimate_rho_cluster(obs, soup, cl)))
    amb <- sum(rho[cl == 6] * colSums(obs)[cl == 6])
    X <- rowSums(obs[, cl == 6, drop = FALSE])
    w <- Decant:::.absence_posterior(X, soup, amb, su)
    ## at most 1% of the genes with counts are called absent
    expect_lte(mean(w[X > 0] > 0.5), 0.01, label = paste("absent share,", ntiny, "cells"))
    tk <- cl == 6
    ret <- function(m) {
      cc <- correct_counts(obs, soup, rho, m, clusters = cl, soup_umis = su)
      sum(pmin(cc[, tk], tru[, tk])) / sum(tru[, tk])
    }
    expect_gte(ret("bayes_absent"), ret("bayes_cell") - 1e-6)
  }
})

test_that("degeneracy cap drops a cluster far above the others", {
  W <- cbind(c(1, 1, 0, 0), c(1, 1, 1, 1), c(1, 0, 0, 0))
  X <- matrix(1, 4, 3)
  out <- Decant:::.absence_degeneracy_cap(W, X, cap = 0.10)
  expect_equal(out[, 2], rep(0, 4))          # 100% vs others' 50% / 25%
  expect_equal(out[, c(1, 3)], W[, c(1, 3)])
  expect_equal(Decant:::.absence_degeneracy_cap(W, X, cap = Inf), W)
})
