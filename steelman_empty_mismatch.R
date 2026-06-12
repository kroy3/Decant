#!/usr/bin/env Rscript
## steelman_empty_mismatch.R
## Fair test of the structured-soup idea in the ONE regime where it could win:
## the empty-droplet pool is contaminated by an "off-manifold" technical profile
## (a junk vector that is NOT a mixture of real cell programs). The in-cell
## contaminant is still the true biological soup. Question: does constraining the
## soup to the span of cell programs project the junk away and beat the global
## estimate that swallows it whole?

for (f in list.files("Decant/R", pattern="\\.R$", full.names=TRUE)) source(f)

run_one <- function(junk_frac, seed) {
  sim <- simulate_experiment(soup_bias = 3, rho_mean = 0.2, n_empty = 4000, seed = seed)
  G <- nrow(sim$observed)

  ## off-manifold junk: a spiky technical profile unrelated to cell programs
  set.seed(seed + 7)
  junk <- rgamma(G, 0.05); junk[sample(G, 20)] <- junk[sample(G, 20)] + 50
  junk <- junk / sum(junk)

  ## corrupt the EMPTY observations only (in-cell contaminant stays = soup_true)
  empty_corrupt <- sim$empty
  ne <- ncol(empty_corrupt)
  add <- vapply(seq_len(ne), function(i) {
    L <- round(junk_frac * sum(empty_corrupt[, i]))
    if (L < 1) return(numeric(G))
    rmultinom(1, L, junk)[, 1]
  }, numeric(G))
  empty_corrupt <- empty_corrupt + add

  soup0   <- ambient_global(empty_corrupt)
  rho_hat <- estimate_rho(sim$observed, soup0)

  ## global: uses corrupted empties directly
  soup_g <- ambient_global(empty_corrupt)
  m_g <- score_correction(sim$observed, correct_counts(sim$observed, soup_g, rho_hat), sim$truth)
  eg  <- soup_profile_error(soup_g, sim$soup_true)

  ## structured: projects corrupted empties onto cell-program basis
  basis <- cluster_basis(sim$observed, k = 6)
  st <- ambient_structured(empty_corrupt, basis)
  m_s <- score_correction(sim$observed, correct_counts(sim$observed, st$soup, rho_hat), sim$truth)
  es  <- soup_profile_error(st$soup, sim$soup_true)

  data.frame(junk_frac, 
             soupErr_global = eg, soupErr_struct = es,
             sigDestroy_global = m_g$signal_destroyed, sigDestroy_struct = m_s$signal_destroyed,
             preserv_global = m_g$preservation, preserv_struct = m_s$preservation)
}

res <- do.call(rbind, lapply(c(0, 0.3, 0.7, 1.5), function(jf)
          do.call(rbind, lapply(1:3, function(s) run_one(jf, 100 + s)))))
agg <- aggregate(. ~ junk_frac, data = res, FUN = mean)

cat("Empty-droplet junk contamination -> does structuring rescue the soup?\n\n")
cat(sprintf("%-9s | %-13s %-13s | %-15s %-15s\n",
            "junkFrac","soupErr_glob","soupErr_str","sigDestroy_glob","sigDestroy_str"))
cat(strrep("-", 78), "\n")
for (i in seq_len(nrow(agg))) {
  cat(sprintf("%-9.1f | %-13.3f %-13.3f | %-15.3f %-15.3f\n",
      agg$junk_frac[i], agg$soupErr_global[i], agg$soupErr_struct[i],
      agg$sigDestroy_global[i], agg$sigDestroy_struct[i]))
}
cat("\nIf soupErr_str stays low while soupErr_glob climbs with junk, the structured\n")
cat("idea is alive but RE-SCOPED: its value is robustness to bad empties, not\n")
cat("biased-soup recovery. If both climb together, the idea is dead.\n")
