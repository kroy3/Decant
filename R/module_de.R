## module_de.R  (GAP 2)
## Differential ambient across conditions creates FALSE-POSITIVE DE genes: a gene
## can look condition-different purely because the soup composition differs. The
## standard pipeline (subtract a point estimate, then test) does not account for
## this. This module instead carries the ambient contribution INTO the DE model
## as a covariate, so ambient-driven differences are absorbed rather than called.

#' Self-contained pseudobulk DE simulator with three gene classes:
#'   null      : no difference anywhere
#'   true_de   : differs in NATIVE expression between conditions (real signal)
#'   ambient_de: native identical, but SOUP differs between conditions (a trap)
simulate_de <- function(n_genes = 400, n_true = 40, n_amb = 40, n_rep = 5,
                        rho_bar = 0.2, fc = 2.5, amb_fc = 3, seed = 1) {
  set.seed(seed)
  base <- rgamma(n_genes, 0.5, 1); base <- base / sum(base)
  true_idx <- 1:n_true
  amb_idx  <- (n_true + 1):(n_true + n_amb)
  cond <- rep(c("A", "B"), each = n_rep)
  S <- length(cond)

  ## native profile per condition: true_de genes up in B
  nat_A <- base; nat_B <- base
  nat_B[true_idx] <- nat_B[true_idx] * fc; nat_B <- nat_B / sum(nat_B)
  nat_A <- nat_A / sum(nat_A)

  ## soup per condition: ambient_de genes enriched in B's soup (the trap)
  soup_A <- base; soup_B <- base
  soup_B[amb_idx] <- soup_B[amb_idx] * amb_fc
  soup_A <- soup_A / sum(soup_A); soup_B <- soup_B / sum(soup_B)

  depth <- 2e5
  obs <- matrix(0, n_genes, S); amb_load <- matrix(0, n_genes, S)
  for (s in seq_len(S)) {
    nat <- if (cond[s] == "A") nat_A else nat_B
    soup <- if (cond[s] == "A") soup_A else soup_B
    lam <- depth * ((1 - rho_bar) * nat + rho_bar * soup)
    obs[, s] <- rpois(n_genes, lam)
    amb_load[, s] <- rho_bar * soup            # estimated ambient fraction (known here)
  }
  list(obs = obs, amb_load = amb_load, cond = factor(cond),
       true_idx = true_idx, amb_idx = amb_idx,
       null_idx = setdiff(seq_len(n_genes), c(true_idx, amb_idx)))
}

.logcpm <- function(m) log1p(sweep(m, 2, colSums(m), "/") * 1e6)

#' Naive: subtract ambient point estimate, then test condition per gene.
de_naive <- function(sim) {
  corrected <- sim$obs - sweep(sim$amb_load, 2, colSums(sim$obs), "*")
  corrected[corrected < 0] <- 0
  y <- .logcpm(corrected)
  apply(y, 1, function(g) tryCatch(
    summary(lm(g ~ sim$cond))$coefficients[2, 4], error = function(e) NA))
}

#' Decant: test condition on OBSERVED data with the ambient load as a covariate,
#' so condition-correlated ambient is absorbed instead of mistaken for signal.
de_ambient_aware <- function(sim) {
  y <- .logcpm(sim$obs)
  al <- .logcpm(sim$amb_load + 1e-9)
  vapply(seq_len(nrow(y)), function(i) tryCatch(
    summary(lm(y[i, ] ~ sim$cond + al[i, ]))$coefficients["sim$condB", 4],
    error = function(e) NA_real_), numeric(1))
}

#' GATE: false-positive rate on ambient-trap genes (should be low for Decant,
#' inflated for naive) while preserving power on true_de genes.
gate_de <- function(seeds = 1:4, alpha = 0.05) {
  fp_n <- c(); fp_d <- c(); pw_n <- c(); pw_d <- c()
  for (s in seeds) {
    sim <- simulate_de(seed = s)
    pn <- de_naive(sim); pd <- de_ambient_aware(sim)
    fp_n <- c(fp_n, mean(pn[sim$amb_idx] < alpha, na.rm = TRUE))
    fp_d <- c(fp_d, mean(pd[sim$amb_idx] < alpha, na.rm = TRUE))
    pw_n <- c(pw_n, mean(pn[sim$true_idx] < alpha, na.rm = TRUE))
    pw_d <- c(pw_d, mean(pd[sim$true_idx] < alpha, na.rm = TRUE))
  }
  cat(sprintf("  ambient-trap FALSE POSITIVE rate | naive=%.2f  Decant=%.2f  %s\n",
              mean(fp_n), mean(fp_d),
              if (mean(fp_d) < mean(fp_n) - 0.05) "DECANT CONTROLS FP" else "no gain"))
  cat(sprintf("  true-DE POWER                    | naive=%.2f  Decant=%.2f  %s\n",
              mean(pw_n), mean(pw_d),
              if (mean(pw_d) > mean(pw_n) - 0.1) "power preserved" else "POWER LOST (bad)"))
}
