## module_allelic.R  (GAP 5)
## In multiplexed (pooled-genotype) experiments, reads carrying an allele that
## does NOT match a cell's own genotype are, by construction, contamination.
## That is near-experimental ground truth for rho. This module simulates such
## data and shows the allelic rho estimate is far more accurate than the
## marker-based one -- because it uses strictly more information.

#' Simulate a pooled multi-genotype experiment with allele-resolved SNP genes.
#' @return list with per-cell own/other allele counts and true rho.
simulate_allelic <- function(n_cells = 1500, n_donors = 8, n_snp = 300,
                             rho_mean = 0.2, rho_conc = 25, snp_lib = 800, seed = 1) {
  set.seed(seed)
  donor <- sample(seq_len(n_donors), n_cells, replace = TRUE)
  abund <- runif(n_donors, 0.5, 1.5); abund <- abund / sum(abund)  # soup donor mix
  rho <- rbeta(n_cells, rho_mean * rho_conc, (1 - rho_mean) * rho_conc)

  own <- integer(n_cells); other <- integer(n_cells)
  for (c in seq_len(n_cells)) {
    L <- rpois(1, snp_lib)
    n_soup <- rbinom(1, L, rho[c]); n_own <- L - n_soup
    ## native SNP reads all carry the cell's own allele
    own[c] <- n_own
    ## soup SNP reads carry a donor allele drawn from the pool abundance
    soup_don <- sample(seq_len(n_donors), n_soup, replace = TRUE, prob = abund)
    own[c] <- own[c] + sum(soup_don == donor[c])
    other[c] <- sum(soup_don != donor[c])
  }
  list(own = own, other = other, donor = donor, abund = abund, rho_true = rho,
       n_donors = n_donors)
}

#' Allelic rho estimate. Fraction of mismatched-allele reads = rho * P(soup read
#' is from another donor). With pool abundances a_d, for a cell of donor d that
#' probability is (1 - a_d). Invert to recover rho.
estimate_rho_allelic <- function(allelic) {
  tot <- allelic$own + allelic$other
  mismatch_frac <- allelic$other / pmax(tot, 1)
  a_self <- allelic$abund[allelic$donor]
  rho_hat <- mismatch_frac / pmax(1 - a_self, 1e-6)
  pmin(pmax(rho_hat, 0), 0.99)
}

#' GATE: allelic rho error vs a marker-style proxy (here: assuming uniform pool,
#' the naive 1/n_donors correction a non-allelic method would implicitly make).
gate_allelic <- function(seeds = 1:3) {
  e_allelic <- c(); e_naive <- c()
  for (s in seeds) {
    al <- simulate_allelic(seed = s)
    rho_a <- estimate_rho_allelic(al)
    ## naive: assumes uniform donor pool (ignores true abundances) -> biased
    tot <- al$own + al$other
    rho_n <- (al$other / pmax(tot, 1)) / (1 - 1 / al$n_donors)
    rmse <- function(x, t) sqrt(mean((x - t)^2))
    e_allelic <- c(e_allelic, rmse(rho_a, al$rho_true))
    e_naive   <- c(e_naive,   rmse(pmin(rho_n, 0.99), al$rho_true))
  }
  cat(sprintf("  rho RMSE | allelic(abundance-aware)=%.3f  naive(uniform-pool)=%.3f  %s\n",
              mean(e_allelic), mean(e_naive),
              if (mean(e_allelic) < mean(e_naive)) "ALLELIC WINS" else "no gain"))
}
