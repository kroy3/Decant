## simulate_groundtruth.R
## Synthetic data with KNOWN clean counts, known soup composition, and a pool
## of empty droplets, so any correction method can be scored against truth.
##
## IMPORTANT HONESTY NOTE (read before trusting any number this produces):
## Synthetic data can be rigged to favour whatever method you are selling.
## The ONLY purpose of this simulator is sufficiency testing and negative
## controls -- e.g. "does the structured-soup idea help when the soup is
## biased, and does it correctly do NOTHING when the soup is uniform?".
## It is NOT evidence that anything beats CellBender. Real claims require
## species-mixing / genotype-mixing data with experimental ground truth.

## Draw K cell-type expression profiles (each a probability vector over genes).
## A block of "marker" genes per type makes cross-type contamination visible.
.make_profiles <- function(n_genes, n_types, n_markers = 15, marker_strength = 40,
                           dirichlet_alpha = 0.3) {
  prof <- matrix(0, nrow = n_genes, ncol = n_types)
  for (k in seq_len(n_types)) {
    base <- rgamma(n_genes, shape = dirichlet_alpha, rate = 1)   # shared-ish background
    idx <- ((k - 1) * n_markers + 1):(k * n_markers)
    idx <- idx[idx <= n_genes]
    base[idx] <- base[idx] + marker_strength                     # type-specific markers
    prof[, k] <- base / sum(base)
  }
  rownames(prof) <- paste0("gene", seq_len(n_genes))
  colnames(prof) <- paste0("type", seq_len(n_types))
  prof
}

#' Simulate a contaminated droplet experiment with known ground truth.
#'
#' @param n_genes,n_types,n_cells dimensions of the experiment.
#' @param soup_bias controls how non-uniform the soup is. 0 = every type lyses
#'   equally (soup == global average of cells, SoupX's assumption holds exactly);
#'   higher values = a few fragile types dominate the soup (the regime where a
#'   structured model should help, and where global models should fail).
#' @param rho_mean mean contamination fraction per cell (snRNA-seq tends to run
#'   higher than scRNA-seq; try 0.1 vs 0.3).
#' @param n_empty number of empty droplets profiled (the soup observation).
#' @return list with $observed, $truth (clean own counts), $empty, $soup_true,
#'   $lysis_true, $rho_true, $labels. Matrices are genes x cells.
simulate_experiment <- function(n_genes = 600, n_types = 6, n_cells = 1200,
                                soup_bias = 3, rho_mean = 0.2, rho_conc = 25,
                                n_empty = 4000, lib_mean = 3000, seed = 1) {
  set.seed(seed)
  prof <- .make_profiles(n_genes, n_types)

  ## lysis weights: which types bleed into the soup. soup_bias=0 -> uniform.
  fragility <- exp(soup_bias * scale(seq_len(n_types))[, 1])
  lysis_w <- fragility / sum(fragility)
  soup_true <- as.numeric(prof %*% lysis_w)
  soup_true <- soup_true / sum(soup_true)

  ## cells
  labels <- sample(seq_len(n_types), n_cells, replace = TRUE)
  libs   <- round(rlnorm(n_cells, log(lib_mean), 0.4))
  rho    <- rbeta(n_cells, rho_mean * rho_conc, (1 - rho_mean) * rho_conc)

  observed <- matrix(0L, n_genes, n_cells)
  truth    <- matrix(0L, n_genes, n_cells)
  for (c in seq_len(n_cells)) {
    T_c <- libs[c]
    n_own  <- rbinom(1, T_c, 1 - rho[c])
    n_soup <- T_c - n_own
    own  <- rmultinom(1, n_own,  prof[, labels[c]])[, 1]
    soup <- rmultinom(1, n_soup, soup_true)[, 1]
    truth[, c]    <- own              # what a perfect method must recover
    observed[, c] <- own + soup       # what the sequencer actually reports
  }

  ## empty droplets = direct (noisy) observations of the soup
  empty_libs <- rpois(n_empty, 25)
  empty <- vapply(empty_libs, function(L) rmultinom(1, L, soup_true)[, 1],
                  numeric(n_genes))

  dn <- list(rownames(prof), paste0("cell", seq_len(n_cells)))
  dimnames(observed) <- dn; dimnames(truth) <- dn
  rownames(empty) <- rownames(prof)

  list(observed = observed, truth = truth, empty = empty,
       soup_true = soup_true, lysis_true = lysis_w, rho_true = rho,
       labels = labels, profiles = prof)
}
