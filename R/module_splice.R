## module_splice.R  (GAP 1)
## Decontaminate the spliced AND unspliced layers using the ambient's splice
## signature, instead of correcting the total and splitting by each cell's own
## observed ratio. Uniquely outputs corrected layers (the RNA-velocity payoff),
## which total-only methods never provide.

#' @param obs_u,obs_s genes x cells observed unspliced / spliced counts.
#' @param emp_u,emp_s genes x droplets empty-droplet unspliced / spliced counts.
#' @param rho per-cell contamination fraction (estimate upstream).
#' @return list with $unspliced, $spliced (corrected, mass-conserving) and the
#'   per-gene ambient unspliced fraction used.
correct_splice_aware <- function(obs_u, obs_s, emp_u, emp_s, rho) {
  amb_u <- rowSums(emp_u); amb_s <- rowSums(emp_s)
  amb_tot <- amb_u + amb_s
  p_amb <- amb_tot / sum(amb_tot)                 # ambient gene distribution
  phi_amb <- amb_u / pmax(amb_tot, 1e-9)          # ambient unspliced fraction per gene

  T_c <- colSums(obs_u) + colSums(obs_s)
  exp_tot <- outer(p_amb, rho * T_c)              # expected ambient total per gene/cell
  exp_u <- exp_tot * phi_amb                       # split by AMBIENT signature, not cell's
  exp_s <- exp_tot * (1 - phi_amb)

  cu <- obs_u - exp_u; cu[cu < 0] <- 0; cu <- pmin(cu, obs_u)
  cs <- obs_s - exp_s; cs[cs < 0] <- 0; cs <- pmin(cs, obs_s)
  list(unspliced = cu, spliced = cs, phi_amb = phi_amb)
}

#' Baseline to beat: correct the TOTAL, then split the removal by each cell's own
#' observed unspliced/spliced ratio (what you get if you decontaminate total and
#' naively apportion to layers).
correct_total_then_split <- function(obs_u, obs_s, emp_u, emp_s, rho) {
  p_amb <- (rowSums(emp_u) + rowSums(emp_s)); p_amb <- p_amb / sum(p_amb)
  obs_t <- obs_u + obs_s
  T_c <- colSums(obs_t)
  exp_tot <- outer(p_amb, rho * T_c)
  removed <- pmin(exp_tot, obs_t)
  frac_u <- obs_u / pmax(obs_t, 1e-9)             # split by the CELL's ratio
  cu <- obs_u - removed * frac_u; cu[cu < 0] <- 0
  cs <- obs_s - removed * (1 - frac_u); cs[cs < 0] <- 0
  list(unspliced = cu, spliced = cs)
}

#' GATE: does splice-aware beat total-then-split on per-layer recovery, and does
#' it correctly NOT help when ambient is not splice-distinct (negative control)?
gate_splice <- function(seeds = 1:3) {
  for (sd in c(0, 0.6)) {                          # splice_distinct: control vs signal
    errs_aware <- c(); errs_base <- c()
    for (s in seeds) {
      sim <- simulate_multimodal(n_cells = 1000, n_samples = 3,
                                 empty_per_sample = c(2500, 2500, 2500),
                                 splice_distinct = sd, seed = s)
      eu <- do.call(cbind, lapply(sim$empties, `[[`, "unspliced"))
      es <- do.call(cbind, lapply(sim$empties, `[[`, "spliced"))
      p0 <- (rowSums(eu) + rowSums(es)); p0 <- p0 / sum(p0)
      rho <- estimate_rho(sim$obs_total, p0)
      a <- correct_splice_aware(sim$obs_unspliced, sim$obs_spliced, eu, es, rho)
      b <- correct_total_then_split(sim$obs_unspliced, sim$obs_spliced, eu, es, rho)
      ## error on the UNSPLICED layer (the one velocity cares about)
      l1 <- function(x, t) sum(abs(x - t)) / sum(t)
      errs_aware <- c(errs_aware, l1(a$unspliced, sim$truth_unspliced))
      errs_base  <- c(errs_base,  l1(b$unspliced, sim$truth_unspliced))
    }
    cat(sprintf("  splice_distinct=%.1f | unspliced L1 err  aware=%.3f  baseline=%.3f  %s\n",
                sd, mean(errs_aware), mean(errs_base),
                if (mean(errs_aware) < mean(errs_base) - 1e-3) "AWARE WINS"
                else if (abs(mean(errs_aware) - mean(errs_base)) <= 1e-3) "tie (expected at 0)"
                else "aware loses"))
  }
}
