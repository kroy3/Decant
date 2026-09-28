## gates_core.R
## Falsification gates for the two components EVERY run depends on: the rho
## estimator and the correction rule. Same rule as the modules: a component
## becomes the default only by winning here, and the gates run on both the
## original (easy) simulator and a stress configuration designed to break the
## absent-gene assumption. Clusters come from k-means, never from oracle labels.

.gate_settings <- function() {
  list(
    easy = list(),
    stress = list(n_genes = 1000, n_types = 8, shared_frac = 0.5, overdispersion = 0.3,
                  type_prob = c(rep(1, 6), 0.05, 0.05), lib_sdlog = 0.6)
  )
}

.gate_sim <- function(setting, rho_mean, seed) {
  args <- c(.gate_settings()[[setting]], list(rho_mean = rho_mean, seed = seed))
  do.call(simulate_experiment, args)
}

#' GATE: is the cluster/absent-gene rho estimator calibrated, and better than
#' the legacy marker-ratio estimator?
#'
#' Pass requires, in every (setting, rho) cell: lower RMSE than legacy, and
#' |mean bias| < max(25% of true rho, 0.015). Negative control: at rho = 0 the
#' estimator must not hallucinate contamination (mean estimate < 0.01). The
#' soup-mismatch tripwire must not fire on any of these clean simulations.
#'
#' @param rho_grid true mean contamination levels to test.
#' @param seeds replicate seeds.
#' @param verbose print the table.
#' @return invisibly, a data.frame of results with attribute `"pass"`.
#' @export
gate_rho <- function(rho_grid = c(0, 0.02, 0.05, 0.1, 0.2, 0.3), seeds = 1:3,
                     verbose = TRUE) {
  rows <- list()
  for (setting in names(.gate_settings())) for (rm in rho_grid) for (s in seeds) {
    sim <- .gate_sim(setting, rm, s)
    soup <- ambient_global(sim$empty)
    cl <- quick_labels(sim$observed, k = length(sim$lysis_true), seed = s)
    r_fit <- estimate_rho_cluster(sim$observed, soup, cl)
    r_new <- as.numeric(r_fit)
    r_old <- estimate_rho(sim$observed, soup)
    t <- sim$rho_true
    rows[[length(rows) + 1]] <- data.frame(
      setting = setting, rho = rm, seed = s, true_mean = mean(t),
      new_mean = mean(r_new), new_rmse = sqrt(mean((r_new - t)^2)),
      new_cor = if (stats::sd(t) > 0) stats::cor(r_new, t) else NA_real_,
      old_mean = mean(r_old), old_rmse = sqrt(mean((r_old - t)^2)),
      false_alarm = isTRUE(attr(r_fit, "soup_mismatch")))
  }
  df <- do.call(rbind, rows)
  agg <- stats::aggregate(cbind(true_mean, new_mean, new_rmse, new_cor, old_mean, old_rmse,
                                false_alarm) ~
                            setting + rho, data = df, FUN = mean, na.action = stats::na.pass)
  agg$bias_ok <- abs(agg$new_mean - agg$true_mean) < pmax(0.25 * agg$true_mean, 0.015)
  agg$beats_legacy <- agg$new_rmse < agg$old_rmse
  neg <- agg$rho == 0
  agg$neg_ok <- ifelse(neg, agg$new_mean < 0.01, NA)
  ## the soup-mismatch tripwire must stay silent on clean data
  agg$no_false_alarm <- agg$false_alarm == 0
  pass <- all(agg$bias_ok) && all(agg$beats_legacy) && all(agg$neg_ok, na.rm = TRUE) &&
          all(agg$no_false_alarm)
  if (verbose) {
    cat(sprintf("  %-7s %5s | %6s | %-22s | %-15s | %s\n", "setting", "rho", "true",
                "cluster: mean rmse cor", "legacy: mean rmse", "ok"))
    for (i in seq_len(nrow(agg))) with(agg[i, ], cat(sprintf(
      "  %-7s %5.2f | %6.3f | %6.3f %6.3f %6s | %6.3f %6.3f | %s\n",
      setting, rho, true_mean, new_mean, new_rmse,
      if (is.na(new_cor)) "  -" else sprintf("%.2f", new_cor),
      old_mean, old_rmse,
      if (bias_ok && beats_legacy && !isFALSE(neg_ok) && no_false_alarm) "ok" else "FAIL")))
    cat(sprintf("  => %s\n", if (pass) "PASS: cluster estimator calibrated and beats legacy"
                             else "FAIL: see rows marked FAIL"))
  }
  attr(agg, "pass") <- pass
  invisible(agg)
}

#' GATE: which mass-conserving correction rule is closest to truth?
#'
#' Scores "subtract" (v0), "redistribute" (SoupX-style) and "bayes" with the
#' ESTIMATED rho (not oracle), on easy and stress simulations. The default in
#' [decant()] is whichever has the lowest mean `l1_error` here.
#'
#' @param rho_grid true mean contamination levels.
#' @param seeds replicate seeds.
#' @param verbose print the table.
#' @return invisibly, a data.frame with attribute `"winner"`.
#' @export
gate_correction <- function(rho_grid = c(0.05, 0.2), seeds = 1:3, verbose = TRUE) {
  methods <- c("subtract", "redistribute", "bayes")
  rows <- list()
  for (setting in names(.gate_settings())) for (rm in rho_grid) for (s in seeds) {
    sim <- .gate_sim(setting, rm, s)
    soup <- ambient_global(sim$empty)
    cl <- quick_labels(sim$observed, k = length(sim$lysis_true), seed = s)
    rho <- as.numeric(estimate_rho_cluster(sim$observed, soup, cl))
    for (m in methods) {
      sc <- score_correction(sim$observed,
                             correct_counts(sim$observed, soup, rho, method = m, clusters = cl),
                             sim$truth)
      rows[[length(rows) + 1]] <- data.frame(
        setting = setting, rho = rm, seed = s, method = m,
        l1_error = sc$l1_error, sensitivity = sc$sensitivity,
        signal_destroyed = sc$signal_destroyed_frac, fabricated = sc$counts_fabricated)
    }
  }
  df <- do.call(rbind, rows)
  agg <- stats::aggregate(cbind(l1_error, sensitivity, signal_destroyed, fabricated) ~
                            setting + rho + method, data = df, FUN = mean)
  overall <- tapply(agg$l1_error, agg$method, mean)
  winner <- names(which.min(overall))
  if (verbose) {
    cat(sprintf("  %-7s %5s %-13s | %8s %8s %9s %4s\n", "setting", "rho", "method",
                "L1err", "sens", "destroyed", "fab"))
    agg <- agg[order(agg$setting, agg$rho, agg$l1_error), ]
    for (i in seq_len(nrow(agg))) with(agg[i, ], cat(sprintf(
      "  %-7s %5.2f %-13s | %8.3f %8.3f %9.3f %4.0f\n",
      setting, rho, method, l1_error, sensitivity, signal_destroyed, fabricated)))
    cat("  mean L1 error across settings: ",
        paste(sprintf("%s=%.3f", names(overall), overall), collapse = "  "), "\n", sep = "")
    cat(sprintf("  => lowest error: %s\n", winner))
  }
  attr(agg, "winner") <- winner
  attr(agg, "overall") <- overall
  invisible(agg)
}
