## benchmark.R
## Run the global-soup baseline against the structured-soup method on the same
## simulated data, holding rho estimation FIXED so the comparison isolates the
## ambient-profile model. Sweep soup_bias to expose the regime where each wins.

#' @param soup_bias_grid vector of soup-bias values to test (0 = uniform soup).
#' @param reps replicate simulations per setting (different seeds).
#' @return data.frame, one row per (method, soup_bias, rep).
run_benchmark <- function(soup_bias_grid = c(0, 1, 2, 4),
                          reps = 3, rho_mean = 0.2,
                          n_genes = 600, n_types = 6, n_cells = 1200, n_empty = 4000) {
  out <- list(); row <- 1L
  for (sb in soup_bias_grid) {
    for (r in seq_len(reps)) {
      sim <- simulate_experiment(n_genes = n_genes, n_types = n_types,
                                 n_cells = n_cells, soup_bias = sb,
                                 rho_mean = rho_mean, n_empty = n_empty,
                                 seed = 100 * r + sb)

      ## SHARED rho, estimated once from a neutral (global) soup so neither
      ## method is advantaged on the rho axis.
      soup0 <- ambient_global(sim$empty)
      rho_hat <- estimate_rho(sim$observed, soup0)

      ## --- baseline: global soup ---
      soup_g <- ambient_global(sim$empty)
      corr_g <- correct_counts(sim$observed, soup_g, rho_hat)
      m_g <- score_correction(sim$observed, corr_g, sim$truth)
      m_g$soup_err <- soup_profile_error(soup_g, sim$soup_true)
      m_g$method <- "global (SoupX/DecontX assumption)"

      ## --- structured soup ---
      basis <- cluster_basis(sim$observed, k = n_types)
      st <- ambient_structured(sim$empty, basis)
      corr_s <- correct_counts(sim$observed, st$soup, rho_hat)
      m_s <- score_correction(sim$observed, corr_s, sim$truth)
      m_s$soup_err <- soup_profile_error(st$soup, sim$soup_true)
      m_s$method <- "structured (Decant)"

      for (m in list(m_g, m_s)) {
        out[[row]] <- data.frame(method = m$method, soup_bias = sb, rep = r,
                                 sensitivity = m$sensitivity,
                                 residual_contam = m$residual_contam_frac,
                                 signal_destroyed = m$signal_destroyed_frac,
                                 preservation = m$preservation_cosine,
                                 soup_err = m$soup_err,
                                 fabricated = m$counts_fabricated,
                                 stringsAsFactors = FALSE)
        row <- row + 1L
      }
    }
  }
  do.call(rbind, out)
}

#' Aggregate to mean +/- sd per (method, soup_bias).
summarise_benchmark <- function(df) {
  agg <- aggregate(cbind(sensitivity, residual_contam, signal_destroyed,
                         preservation, soup_err, fabricated) ~ method + soup_bias,
                   data = df, FUN = function(x) c(mean = mean(x), sd = sd(x)))
  agg
}
