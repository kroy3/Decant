## module_diagnose.R  (GAP 4)
## The structured-soup model FAILED as a corrector (proven earlier). Its honest
## residual value is diagnostic: recover per-program lysis weights so a user can
## see which populations are seeding the soup. Tested with an oracle basis to
## isolate the recovery mechanism from clustering error (which degrades it in
## practice -- stated, not hidden).

#' GATE: with a clean (oracle) program basis, do recovered lysis weights track
#' the true lysis weights? Reported as rank correlation; clustering noise lowers
#' this in real use.
gate_diagnose <- function(seeds = 1:4) {
  rhos <- c()
  for (s in seeds) {
    sim <- simulate_experiment(soup_bias = 3, n_empty = 4000, seed = s)
    st <- ambient_structured(sim$empty, sim$profiles)     # oracle basis = true type means
    rhos <- c(rhos, suppressWarnings(cor(st$lysis, sim$lysis_true, method = "spearman")))
  }
  cat(sprintf("  lysis-weight recovery (oracle basis) | mean Spearman rho = %.3f  %s\n",
              mean(rhos), if (mean(rhos) > 0.6) "DIAGNOSTIC USABLE (as QC, not correction)"
                          else "weak"))
  cat("  note: with k-means basis instead of oracle, expect this to drop; it is a\n")
  cat("        QC signal, not a correction method.\n")
}
