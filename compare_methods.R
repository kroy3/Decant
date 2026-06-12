#!/usr/bin/env Rscript
## compare_methods.R
## Honest 4-way comparison on ground-truth simulations. No method is tuned to
## win. Whatever the table says is the result.
##
##  - none      : no correction (floor / sanity check)
##  - SoupX-like: global soup from empties + shared rho + subtract
##  - DecontX-like: Bayesian-mixture EM, no empties, infers from clusters
##  - Decant    : structured soup (mixture of cell programs) + shared rho
##
## Scored identically. rho swept to cover scRNA-like (0.1) and snRNA-like (0.3).

for (f in list.files("Decant/R", pattern="\\.R$", full.names=TRUE)) source(f)

eval_all <- function(rho_mean, soup_bias = 3, seed = 1) {
  sim <- simulate_experiment(soup_bias = soup_bias, rho_mean = rho_mean,
                             n_empty = 4000, seed = seed)
  obs <- sim$observed; truth <- sim$truth

  ## shared pieces
  soup0   <- ambient_global(sim$empty)
  rho_hat <- estimate_rho(obs, soup0)
  labels  <- quick_labels(obs, k = 6, seed = seed)

  res <- list()

  ## none
  res$none <- score_correction(obs, obs, truth)

  ## SoupX-like (global soup + shared rho)
  cg <- correct_counts(obs, ambient_global(sim$empty), rho_hat)
  res$soupx <- score_correction(obs, cg, truth)

  ## DecontX-like (EM, its own theta, no empties)
  dx <- decontx_em(obs, labels, iters = 25)
  res$decontx <- score_correction(obs, dx$corrected, truth)

  ## Decant (structured soup + shared rho)
  basis <- cluster_basis(obs, k = 6, seed = seed)
  st <- ambient_structured(sim$empty, basis)
  cs <- correct_counts(obs, st$soup, rho_hat)
  res$decant <- score_correction(obs, cs, truth)

  do.call(rbind, lapply(names(res), function(nm) {
    m <- res[[nm]]
    data.frame(rho_mean = rho_mean, method = nm,
               sensitivity = m$sensitivity,
               residual_contam = m$residual_contam_frac,
               signal_destroyed = m$signal_destroyed_frac,
               preservation = m$preservation_cosine,
               fabricated = m$counts_fabricated)
  }))
}

grid <- expand.grid(rho = c(0.1, 0.3), seed = 1:4)
all <- do.call(rbind, Map(function(r, s) eval_all(r, seed = s), grid$rho, grid$seed))

agg <- aggregate(cbind(sensitivity, residual_contam, signal_destroyed, preservation, fabricated) ~
                   method + rho_mean, data = all, FUN = mean)
agg <- agg[order(agg$rho_mean, -agg$sensitivity), ]

lab <- c(none="none", soupx="SoupX-like", decontx="DecontX-like", decant="Decant")
cat("\n4-way benchmark on ground-truth simulation (mean over 4 seeds)\n")
cat("higher sensitivity & preservation = better; lower residual & signalDestroyed = better\n\n")
cat(sprintf("%-5s %-13s | %-11s %-12s %-13s %-11s %-5s\n",
            "rho","method","sensitivity","residContam","signalDestroy","preserv","fab"))
cat(strrep("-", 78), "\n")
last <- NA
for (i in seq_len(nrow(agg))) {
  if (!is.na(last) && last != agg$rho_mean[i]) cat("\n")
  cat(sprintf("%-5.1f %-13s | %-11.3f %-12.3f %-13.3f %-11.3f %-5.0f\n",
      agg$rho_mean[i], lab[agg$method[i]], agg$sensitivity[i], agg$residual_contam[i],
      agg$signal_destroyed[i], agg$preservation[i], agg$fabricated[i]))
  last <- agg$rho_mean[i]
}
