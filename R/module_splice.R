## module_splice.R  (GAP 1)
## Decontaminate the spliced AND unspliced layers using the ambient's splice
## signature, instead of correcting the total and splitting by each cell's own
## observed ratio. Uniquely outputs corrected layers (the RNA-velocity payoff),
## which total-only methods never provide.

#' Splice-aware layer decontamination.
#'
#' @param obs_u,obs_s genes x cells observed unspliced / spliced counts.
#' @param emp_u,emp_s genes x droplets empty-droplet unspliced / spliced counts.
#' @param rho per-cell contamination fraction (estimate upstream).
#' @param drop_zeros for sparse input, drop entries corrected to exactly 0.
#' @return list with $unspliced, $spliced (corrected, mass-conserving) and the
#'   per-gene ambient unspliced fraction used.
#' @export
correct_splice_aware <- function(obs_u, obs_s, emp_u, emp_s, rho, drop_zeros = TRUE) {
  obs_u <- .as_counts(obs_u, "obs_unspliced"); obs_s <- .as_counts(obs_s, "obs_spliced")
  if (!identical(dim(obs_u), dim(obs_s))) stop("unspliced/spliced layers differ in shape", call. = FALSE)
  amb_u <- .row_sums(.as_counts(emp_u, "empties_unspliced"))
  amb_s <- .row_sums(.as_counts(emp_s, "empties_spliced"))
  amb_tot <- amb_u + amb_s
  p_amb <- amb_tot / sum(amb_tot)                 # ambient gene distribution
  phi_amb <- amb_u / pmax(amb_tot, 1e-9)          # ambient unspliced fraction per gene
  T_c <- .col_sums(obs_u) + .col_sums(obs_s)

  ## Work on entries where EITHER layer is non-zero (a zero total has nothing
  ## to remove). Keys are doubles: G * N can exceed .Machine$integer.max.
  G <- nrow(obs_u)
  nu <- .nz(obs_u); ns <- .nz(obs_s)
  ku <- (nu$j - 1) * G + nu$i; ks <- (ns$j - 1) * G + ns$i
  key <- sort(unique(c(ku, ks)))
  mu <- match(ku, key); ms <- match(ks, key)
  u <- numeric(length(key)); u[mu] <- nu$v
  s <- numeric(length(key)); s[ms] <- ns$v
  i <- (key - 1) %% G + 1; j <- (key - 1) %/% G + 1

  ## Clamp the TOTAL removal once (clamping each layer separately loses mass
  ## twice and made the negative control fail), split it by the AMBIENT splice
  ## signature, and route any overflow to the layer that still has counts.
  tot_rm <- pmin(p_amb[i] * rho[j] * T_c[j], u + s)
  ru <- pmin(tot_rm * phi_amb[i], u)
  rs <- pmin(tot_rm * (1 - phi_amb[i]), s)
  left <- tot_rm - ru - rs
  add_u <- pmin(left, u - ru); ru <- ru + add_u
  rs <- rs + pmin(left - add_u, s - rs)

  list(unspliced = .nz_set(obs_u, nu, pmin(pmax(u - ru, 0), u)[mu], drop = drop_zeros),
       spliced   = .nz_set(obs_s, ns, pmin(pmax(s - rs, 0), s)[ms], drop = drop_zeros),
       phi_amb = phi_amb)
}

#' Baseline to beat: correct the TOTAL, then split the removal by each cell's own
#' observed unspliced/spliced ratio (what you get if you decontaminate total and
#' naively apportion to layers).
#' @inheritParams correct_splice_aware
#' @return list with $unspliced, $spliced.
#' @export
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
#' @param seeds replicate seeds.
#' @export
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
      rho <- as.numeric(estimate_rho_cluster(sim$obs_total, p0,
                                             quick_labels(sim$obs_total, k = 6, seed = s)))
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
