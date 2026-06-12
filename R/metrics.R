## metrics.R
## Score a corrected matrix against known-clean truth. Mirrors the axes used by
## the published benchmarks: contamination removed, endogenous signal preserved,
## and an explicit fabrication audit.

#' @param observed,corrected,truth genes x cells matrices (truth = clean own counts).
#' @return named list of scalar metrics.
score_correction <- function(observed, corrected, truth) {
  injected  <- observed - truth                # the soup counts that were added (>=0)
  removed   <- observed - corrected            # what the method took out (>=0 by guarantee)
  injected[injected < 0] <- 0

  tot_injected <- sum(injected)
  tot_truth    <- sum(truth)

  ## of the contamination that was actually present, how much got removed
  ## (capped per element so over-removal does not count as a win)
  true_removed <- sum(pmin(removed, injected))
  sensitivity  <- if (tot_injected > 0) true_removed / tot_injected else NA_real_

  ## endogenous damage: counts removed beyond the contamination = destroyed signal
  over_removed <- sum(pmax(removed - injected, 0))
  signal_destroyed_frac <- over_removed / tot_truth

  ## residual contamination left behind
  residual_contam_frac <- (tot_injected - true_removed) / tot_injected

  ## per-gene preservation of the true profile (mean across cells of cosine sim
  ## between corrected and true cell vectors)
  cos <- function(a, b) {
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2)); if (d == 0) return(NA_real_); sum(a * b) / d
  }
  preservation <- mean(vapply(seq_len(ncol(truth)),
                              function(c) cos(corrected[, c], truth[, c]),
                              numeric(1)), na.rm = TRUE)

  ## FABRICATION AUDIT: must be exactly 0 for a mass-conserving method.
  fabricated <- sum(pmax(corrected - observed, 0))

  list(
    sensitivity            = sensitivity,            # higher better (0-1)
    residual_contam_frac   = residual_contam_frac,   # lower better
    signal_destroyed_frac  = signal_destroyed_frac,  # lower better
    preservation_cosine    = preservation,           # higher better (~1)
    counts_fabricated      = fabricated               # MUST be 0
  )
}

#' Soup-profile recovery error: how close the estimated ambient profile is to
#' the true soup. This is the mechanism the structured method is meant to fix,
#' so report it directly. L1 distance between probability vectors (0 = perfect).
soup_profile_error <- function(soup_hat, soup_true) {
  sum(abs(soup_hat / sum(soup_hat) - soup_true / sum(soup_true)))
}
