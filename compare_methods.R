#!/usr/bin/env Rscript
## compare_methods.R
## Honest 4-way comparison on ground-truth simulations. No method is tuned to
## win. Whatever the table says is the result.
##
##  - none        : no correction (floor / sanity check)
##  - v0          : global soup + legacy marker-ratio rho + subtract (Decant 0.1)
##  - DecontX-like: Bayesian-mixture EM, no empties, infers from clusters
##  - Decant      : decant() as shipped (global soup + cluster rho + subtract)
##
## NOTE: single-sample Decant is SoupX's model (global soup, subtract) with a
## different rho estimator. The real SoupX autoEstCont is NOT reimplemented
## here; "beats v0" is not "beats SoupX". Run the official package for that.
##
## Scored identically. rho swept to cover scRNA-like (0.1) and snRNA-like (0.3).

## Run from the repo root. Uses the installed package, or the source tree via pkgload.
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)

eval_all <- function(rho_mean, soup_bias = 3, seed = 1) {
  sim <- simulate_experiment(soup_bias = soup_bias, rho_mean = rho_mean,
                             n_empty = 4000, seed = seed)
  obs <- sim$observed; truth <- sim$truth

  ## shared pieces
  soup0   <- ambient_global(sim$empty)
  labels  <- quick_labels(obs, k = 6, seed = seed)

  res <- list()

  ## none
  res$none <- score_correction(obs, obs, truth)

  ## v0 (global soup + legacy rho)
  cg <- correct_counts(obs, soup0, estimate_rho(obs, soup0))
  res$v0 <- score_correction(obs, cg, truth)

  ## DecontX-like (EM, its own theta, no empties)
  dx <- decontx_em(obs, labels, iters = 25)
  res$decontx <- score_correction(obs, dx$corrected, truth)

  ## Decant as shipped
  res$decant <- score_correction(obs, decant(obs, sim$empty, clusters = labels)$corrected, truth)

  do.call(rbind, lapply(names(res), function(nm) {
    m <- res[[nm]]
    data.frame(rho_mean = rho_mean, method = nm,
               sensitivity = m$sensitivity,
               residual_contam = m$residual_contam_frac,
               signal_destroyed = m$signal_destroyed_frac,
               preservation = m$preservation_cosine,
               l1_error = m$l1_error,
               fabricated = m$counts_fabricated)
  }))
}

grid <- expand.grid(rho = c(0.1, 0.3), seed = 1:4)
all <- do.call(rbind, Map(function(r, s) eval_all(r, seed = s), grid$rho, grid$seed))

agg <- aggregate(cbind(sensitivity, residual_contam, signal_destroyed, preservation, l1_error, fabricated) ~
                   method + rho_mean, data = all, FUN = mean)
agg <- agg[order(agg$rho_mean, agg$l1_error), ]

lab <- c(none="none", v0="v0 (legacy rho)", decontx="DecontX-like", decant="Decant")
cat("\n4-way benchmark on ground-truth simulation (mean over 4 seeds)\n")
cat("higher sensitivity & preservation = better; lower residual, signalDestroyed, L1err = better\n\n")
cat(sprintf("%-5s %-15s | %-11s %-12s %-13s %-9s %-7s %-5s\n",
            "rho","method","sensitivity","residContam","signalDestroy","preserv","L1err","fab"))
cat(strrep("-", 88), "\n")
last <- NA
for (i in seq_len(nrow(agg))) {
  if (!is.na(last) && last != agg$rho_mean[i]) cat("\n")
  cat(sprintf("%-5.1f %-15s | %-11.3f %-12.3f %-13.3f %-9.3f %-7.3f %-5.0f\n",
      agg$rho_mean[i], lab[agg$method[i]], agg$sensitivity[i], agg$residual_contam[i],
      agg$signal_destroyed[i], agg$preservation[i], agg$l1_error[i], agg$fabricated[i]))
  last <- agg$rho_mean[i]
}
