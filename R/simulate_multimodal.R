## simulate_multimodal.R
## Extends the ground-truth simulator to emit the ORTHOGONAL signals the new
## modules need: separate spliced/unspliced layers and multi-sample structure.
## Allelic data has a different shape and gets its own simulator (module_allelic.R).
##
## Same honesty caveat as before: this is for sufficiency tests and negative
## controls only. A module that cannot beat baseline on data DESIGNED to contain
## its signal is dead; a module that wins here still must be re-tested on real
## species/genotype-mixing data before any claim.

.profiles <- function(n_genes, n_types, n_markers = 15, marker_strength = 40, alpha = 0.3) {
  prof <- matrix(0, n_genes, n_types)
  for (k in seq_len(n_types)) {
    base <- rgamma(n_genes, alpha, 1)
    idx <- ((k - 1) * n_markers + 1):(k * n_markers); idx <- idx[idx <= n_genes]
    base[idx] <- base[idx] + marker_strength
    prof[, k] <- base / sum(base)
  }
  prof
}

#' Multi-sample, splice-resolved contaminated experiment with ground truth.
#'
#' @param n_samples number of samples; each gets its own soup that shares a global
#'   structure but deviates (between_sample_var), and its own empty-droplet pool.
#' @param empty_per_sample vector (length n_samples) of empties per sample; set
#'   some low to stress hierarchical pooling.
#' @param splice_distinct in [0,1]: how differently ambient vs native split into
#'   unspliced/spliced. 0 = identical (negative control for the splice module).
#' @return list of layers (genes x cells) for spliced/unspliced observed and
#'   truth, per-sample empties (with splice layers), sample/label vectors, and
#'   the true soup profiles.
simulate_multimodal <- function(n_genes = 600, n_types = 6, n_cells = 1500,
                                n_samples = 6, empty_per_sample = NULL,
                                soup_bias = 3, rho_mean = 0.2, rho_conc = 25,
                                between_sample_var = 0.15, splice_distinct = 0.6,
                                lib_mean = 3000, seed = 1) {
  set.seed(seed)
  if (is.null(empty_per_sample))
    empty_per_sample <- rep(c(3000, 300), length.out = n_samples)  # mix of rich/poor
  prof <- .profiles(n_genes, n_types)

  ## native vs ambient splice fractions (prob a molecule is UNSPLICED/intronic).
  ## nuclei native skews unspliced (high); ambient (cytoplasmic) skews spliced (low).
  pu_native <- rbeta(n_genes, 4, 2)                              # ~0.67 mean, unspliced-rich
  pu_ambient <- pu_native * (1 - splice_distinct)               # ambient shifted to spliced
  ## when splice_distinct=0, pu_ambient==pu_native (no usable signal)

  ## global lysis weights -> global soup; per-sample deviation
  fragility <- exp(soup_bias * scale(seq_len(n_types))[, 1])
  lysis_g <- fragility / sum(fragility)

  sample_of <- sort(sample(seq_len(n_samples), n_cells, replace = TRUE))
  labels <- sample(seq_len(n_types), n_cells, replace = TRUE)
  libs <- round(rlnorm(n_cells, log(lib_mean), 0.4))
  rho <- rbeta(n_cells, rho_mean * rho_conc, (1 - rho_mean) * rho_conc)

  Su <- matrix(0L, n_genes, n_cells); Ss <- matrix(0L, n_genes, n_cells)   # observed
  Tu <- matrix(0L, n_genes, n_cells); Ts <- matrix(0L, n_genes, n_cells)   # truth (native)

  soup_true <- vector("list", n_samples)
  for (s in seq_len(n_samples)) {
    dev <- rgamma(n_types, 1 / between_sample_var, 1 / between_sample_var)  # mean 1, variable
    ly_s <- lysis_g * dev; ly_s <- ly_s / sum(ly_s)
    sp <- as.numeric(prof %*% ly_s); sp <- sp / sum(sp)
    soup_true[[s]] <- sp
  }

  split_counts <- function(total_per_gene, pu) {
    u <- rbinom(length(total_per_gene), total_per_gene, pu)
    list(u = u, s = total_per_gene - u)
  }

  for (c in seq_len(n_cells)) {
    s <- sample_of[c]; T_c <- libs[c]
    n_own <- rbinom(1, T_c, 1 - rho[c]); n_soup <- T_c - n_own
    own  <- rmultinom(1, n_own,  prof[, labels[c]])[, 1]
    soup <- rmultinom(1, n_soup, soup_true[[s]])[, 1]
    own_sp  <- split_counts(own,  pu_native)
    soup_sp <- split_counts(soup, pu_ambient)
    Tu[, c] <- own_sp$u;            Ts[, c] <- own_sp$s
    Su[, c] <- own_sp$u + soup_sp$u; Ss[, c] <- own_sp$s + soup_sp$s
  }

  ## empties per sample, splice-resolved (observe the ambient splice signature)
  empties <- vector("list", n_samples)
  for (s in seq_len(n_samples)) {
    M <- empty_per_sample[s]; L <- rpois(M, 25)
    eu <- matrix(0L, n_genes, M); es <- matrix(0L, n_genes, M)
    for (j in seq_len(M)) {
      tot <- rmultinom(1, L[j], soup_true[[s]])[, 1]
      sp <- split_counts(tot, pu_ambient)
      eu[, j] <- sp$u; es[, j] <- sp$s
    }
    empties[[s]] <- list(unspliced = eu, spliced = es)
  }

  list(obs_unspliced = Su, obs_spliced = Ss,
       truth_unspliced = Tu, truth_spliced = Ts,
       obs_total = Su + Ss, truth_total = Tu + Ts,
       empties = empties, sample_of = sample_of, labels = labels,
       rho_true = rho, soup_true = soup_true,
       pu_native = pu_native, pu_ambient = pu_ambient,
       n_samples = n_samples, empty_per_sample = empty_per_sample)
}
