## module_hier.R  (GAP 3)
## Partial-pooling of per-sample soup estimates. Samples with few empties get
## shrunk toward the global mean (where SoupX's own docs say single-sample
## estimation breaks); samples with rich empties keep their own estimate.

#' Empirical-Bayes-style shrinkage of per-sample soup profiles.
#' @param empties_total list of genes x droplets matrices, one per sample.
#' @return list with $pooled (genes x samples) and $independent (genes x samples).
pool_soup <- function(empties_total) {
  S <- length(empties_total)
  ## independent per-sample estimate + the count mass behind it
  indep <- vapply(empties_total, function(e) { p <- rowSums(e); p / sum(p) },
                  numeric(nrow(empties_total[[1]])))
  mass <- vapply(empties_total, function(e) sum(e), numeric(1))
  global <- rowMeans(indep)

  ## between-sample variance per gene (signal) vs within-sample sampling noise.
  ## shrink weight lambda_s grows with that sample's mass relative to a tau set
  ## from the average mass (more empties -> trust the sample's own estimate).
  tau <- median(mass)
  lambda <- mass / (mass + tau)                       # in (0,1), per sample
  pooled <- indep
  for (s in seq_len(S))
    pooled[, s] <- lambda[s] * indep[, s] + (1 - lambda[s]) * global
  pooled <- sweep(pooled, 2, colSums(pooled), "/")
  list(pooled = pooled, independent = indep, lambda = lambda, mass = mass)
}

#' GATE: for empty-poor samples, is the pooled soup closer to truth than the
#' independent estimate? For empty-rich samples, are they about equal (pooling
#' should not hurt)?
gate_hier <- function(seeds = 1:3) {
  poor_gain <- c(); rich_gain <- c()
  for (s in seeds) {
    sim <- simulate_multimodal(n_cells = 1500, n_samples = 6,
                               empty_per_sample = c(4000, 4000, 4000, 200, 200, 200),
                               between_sample_var = 0.15, seed = s)
    emp_tot <- lapply(sim$empties, function(e) e$unspliced + e$spliced)
    pp <- pool_soup(emp_tot)
    l1 <- function(p, t) sum(abs(p / sum(p) - t / sum(t)))
    for (j in seq_len(sim$n_samples)) {
      e_indep <- l1(pp$independent[, j], sim$soup_true[[j]])
      e_pool  <- l1(pp$pooled[, j], sim$soup_true[[j]])
      if (sim$empty_per_sample[j] <= 300) poor_gain <- c(poor_gain, e_indep - e_pool)
      else rich_gain <- c(rich_gain, e_indep - e_pool)
    }
  }
  cat(sprintf("  empty-POOR samples | mean soup-err reduction from pooling: %+.3f  %s\n",
              mean(poor_gain), if (mean(poor_gain) > 0) "POOLING HELPS" else "pooling hurts"))
  cat(sprintf("  empty-RICH samples | mean soup-err reduction from pooling: %+.3f  %s\n",
              mean(rich_gain), if (abs(mean(rich_gain)) < 0.02) "~neutral (good)" else
                if (mean(rich_gain) < 0) "pooling hurts rich (bad)" else "pooling helps"))
}
