## A synthetic two-species mixture with KNOWN per-cell rho, used to check that
## species_mix_truth() recovers rho from cross-species reads alone.
sim_species_mix <- function(n_genes = 400, n_cells = 600, n_empty = 3000,
                            rho_mean = 0.1, human_frac = 0.6, seed = 1) {
  set.seed(seed)
  G2 <- n_genes / 2
  genes <- c(sprintf("GRCh38_G%04d", seq_len(G2)), sprintf("mm10___G%04d", seq_len(G2)))
  mk <- function(k, sp) {                       # k types per species
    sapply(seq_len(k), function(j) {
      p <- numeric(n_genes); idx <- if (sp == "h") seq_len(G2) else G2 + seq_len(G2)
      b <- stats::rgamma(G2, 0.3, 1); b[((j - 1) * 15 + 1):(j * 15)] <- b[((j - 1) * 15 + 1):(j * 15)] + 30
      p[idx] <- b / sum(b); p })
  }
  prof <- cbind(mk(3, "h"), mk(3, "m"))
  sp_of_type <- rep(c("GRCh38", "mm10"), each = 3)
  type <- ifelse(stats::runif(n_cells) < human_frac, sample(1:3, n_cells, TRUE), sample(4:6, n_cells, TRUE))
  soup <- rowMeans(prof[, type]); soup <- soup / sum(soup)
  rho <- stats::rbeta(n_cells, rho_mean * 25, (1 - rho_mean) * 25)
  lib <- round(stats::rlnorm(n_cells, log(3000), 0.4))
  obs <- sapply(seq_len(n_cells), function(c) {
    no <- stats::rbinom(1, lib[c], 1 - rho[c])
    stats::rmultinom(1, no, prof[, type[c]])[, 1] + stats::rmultinom(1, lib[c] - no, soup)[, 1] })
  emp <- sapply(stats::rpois(n_empty, 25), function(L) stats::rmultinom(1, L, soup)[, 1])
  dimnames(obs) <- list(genes, paste0("c", seq_len(n_cells))); rownames(emp) <- genes
  list(cells = obs, empties = emp, rho = rho, species = sp_of_type[type], genes = genes)
}
