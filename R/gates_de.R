## gates_de.R
## Falsification gate for contamination-aware DE. Each method is run with its
## own standard gene filter; leakage is counted over ALL of the abundant type's
## disease genes (a gene a method declines to test counts as not called).

.de_case_scores <- function(st, alpha = 0.05, fdr = 0.1) {
  keep <- st$label == "R"; S <- sort(unique(st$sample))
  Y <- vapply(S, function(i) Matrix::rowSums(st$X[, keep & st$sample == i, drop = FALSE]),
              numeric(nrow(st$X)))
  rownames(Y) <- rownames(st$X)
  design <- stats::model.matrix(~st$cond)
  res <- ambient_de(Y, st$empties, design, naive = FALSE)
  kn <- edgeR::filterByExpr(edgeR::DGEList(Y), design)
  d0 <- edgeR::estimateDisp(.norm_lib_sizes(edgeR::DGEList(Y[kn, ])), design)
  pn <- edgeR::glmQLFTest(edgeR::glmQLFit(d0, design), coef = 2)$table$PValue
  names(pn) <- rownames(Y)[kn]
  gi <- seq_len(nrow(Y))
  is_de <- gi %in% st$de_rare; is_leak <- gi %in% st$de_abund; is_sh <- gi %in% st$shared
  one <- function(pv, m) {
    p <- stats::setNames(rep(1, nrow(Y)), rownames(Y)); p[names(pv)] <- pv
    q <- p; q[] <- 1; q[names(pv)] <- stats::p.adjust(pv, "BH")
    tested <- rownames(Y) %in% names(pv)
    data.frame(method = m, leak_called = sum(q[is_leak] < fdr),
               null_fp = mean(p[tested & !is_de & !is_leak] < alpha),
               power = mean(q[is_de & !is_sh] < fdr),
               power_shared = if (any(is_sh)) mean(q[is_sh] < fdr) else NA_real_,
               fdp = if (any(q < fdr)) mean(!is_de[q < fdr]) else 0)
  }
  rbind(one(stats::setNames(res$PValue, res$gene), "ambient_de"), one(pn, "naive_edgeR"))
}

#' GATE: does contamination-aware DE stop ambient leakage while staying
#' calibrated and powerful?
#'
#' Simulated 6 vs 6 case/control studies in which disease tissue carries more
#' ambient RNA. Tests DE in a rare cell type; the abundant type's disease
#' genes leak in through the soup. Pass requires, averaged over seeds, in
#' every case: at most 1 leaked gene called, null false-positive rate <= 0.075,
#' false discovery proportion <= 0.15, and power within 0.1 of naive edgeR
#' (including on genes truly DE in both types).
#'
#' @param seeds replicate seeds.
#' @param verbose print the table.
#' @return invisibly, a data.frame with attribute `"pass"`.
#' @export
gate_ambient_de <- function(seeds = 1:2, verbose = TRUE) {
  for (p in c("edgeR", "limma"))
    if (!requireNamespace(p, quietly = TRUE)) stop("gate_ambient_de() needs ", p, call. = FALSE)
  cases <- list(base = list(), neg_control = list(rho_dis = 0.05),
                shared_true_de = list(n_shared = 20))
  rows <- list()
  for (nm in names(cases)) for (s in seeds) {
    st <- do.call(simulate_ambient_study, c(cases[[nm]], list(seed = s)))
    rows[[length(rows) + 1]] <- cbind(case = nm, .de_case_scores(st))
  }
  df <- do.call(rbind, rows)
  agg <- stats::aggregate(cbind(leak_called, null_fp, power, fdp) ~ case + method, df, mean)
  sh <- stats::aggregate(power_shared ~ case + method, df[!is.na(df$power_shared), ], mean)
  agg <- merge(agg, sh, all.x = TRUE)
  a <- agg[agg$method == "ambient_de", ]; n <- agg[agg$method == "naive_edgeR", ]
  n <- n[match(a$case, n$case), ]
  a$ok <- a$leak_called <= 1 & a$null_fp <= 0.075 & a$fdp <= 0.15 &
    a$power >= n$power - 0.1 & (is.na(a$power_shared) | a$power_shared >= n$power_shared - 0.1)
  pass <- all(a$ok)
  if (verbose) {
    cat(sprintf("  %-15s %-12s | %8s %8s %6s %7s %6s\n", "case", "method", "leaked", "null_fp",
                "power", "shared", "FDP"))
    for (i in seq_len(nrow(agg))) with(agg[i, ], cat(sprintf(
      "  %-15s %-12s | %8.1f %8.3f %6.2f %7s %6.2f\n", case, method, leak_called, null_fp, power,
      if (is.na(power_shared)) "-" else sprintf("%.2f", power_shared), fdp)))
    cat(sprintf("  => %s\n", if (pass) "PASS: no leakage, calibrated, power kept"
                             else "FAIL: see ambient_de rows"))
  }
  attr(agg, "pass") <- pass
  invisible(agg)
}

## Vectorised Wilcoxon rank-sum test per row (normal approximation with tie
## and continuity corrections), for the large-cohort gate's comparators.
.wilcox_rows <- function(X, g) {
  g <- as.integer(factor(g)) == 2L; n1 <- sum(g); n0 <- sum(!g); n <- n0 + n1
  R <- t(apply(X, 1, rank))
  W <- rowSums(R[, g, drop = FALSE]) - n1 * (n1 + 1) / 2
  ties <- apply(R, 1, function(r) { t <- table(r); sum(t^3 - t) })
  v <- n0 * n1 / 12 * ((n + 1) - ties / (n * (n - 1)))
  d <- W - n0 * n1 / 2
  z <- (d - sign(d) * 0.5) / sqrt(v)
  stats::setNames(2 * stats::pnorm(-abs(z)), rownames(X))
}

.cohort_scores <- function(st, alpha = 0.05, fdr = 0.1) {
  Y <- st$Y
  design <- if (length(st$batch_genes)) stats::model.matrix(~st$batch + st$cond)
            else stats::model.matrix(~st$cond)
  cf <- ncol(design)
  pv <- list()
  for (rb in c(TRUE, FALSE)) {
    r <- ambient_de(Y, st$empties, design, coef = cf, naive = FALSE, robust = rb)
    pv[[if (rb) "ambient_de" else "ambient_de_nonrobust"]] <- stats::setNames(r$PValue, r$gene)
  }
  kn <- edgeR::filterByExpr(edgeR::DGEList(Y), design)
  d0 <- edgeR::estimateDisp(.norm_lib_sizes(edgeR::DGEList(Y[kn, ])), design)
  pv$naive_edgeR <- stats::setNames(
    edgeR::glmQLFTest(edgeR::glmQLFit(d0, design, robust = TRUE), coef = cf)$table$PValue,
    rownames(Y)[kn])
  pv$naive_wilcoxon <- .wilcox_rows(edgeR::cpm(.norm_lib_sizes(edgeR::DGEList(Y[kn, ]))), st$cond)
  ## subtract the same ambient estimate ambient_de uses, then rank-test
  sp <- soup_eb(st$empties, design)
  rho <- vapply(seq_len(ncol(Y)), function(s)
    as.numeric(estimate_rho_cluster(Y[, s, drop = FALSE], sp$soup[, s], clusters = 1L)), 1)
  rho <- .joint_rho(Y, sp$soup, rho)$rho
  nat <- pmax(Y - sweep(sp$soup, 2, rho * colSums(Y), "*"), 0)
  ks <- edgeR::filterByExpr(nat, design = design)
  pv$subtract_wilcoxon <- .wilcox_rows(edgeR::cpm(.norm_lib_sizes(edgeR::DGEList(nat[ks, ]))), st$cond)

  gi <- seq_len(nrow(Y)); gn <- rownames(Y)
  is_de <- gi %in% st$de_rare; is_leak <- gi %in% st$de_abund
  do.call(rbind, lapply(names(pv), function(m) {
    p <- stats::setNames(rep(1, nrow(Y)), gn); p[names(pv[[m]])] <- pv[[m]]
    q <- p; q[] <- 1; q[names(pv[[m]])] <- stats::p.adjust(pv[[m]], "BH")
    tested <- gn %in% names(pv[[m]])
    data.frame(method = m, leak_called = sum(q[is_leak] < fdr),
               null_fp = mean(p[tested & !is_de & !is_leak] < alpha),
               power = mean(q[is_de] < fdr),
               fdp = if (any(q < fdr)) mean(!is_de[q < fdr]) else 0)
  }))
}

#' GATE: does contamination-aware DE stay calibrated in large cohorts whose
#' counts are not negative binomial?
#'
#' Simulated case/control cohorts ([simulate_ambient_cohort()]) with
#' `n_per` samples per group, in four cases: clean NB-like counts, per-gene
#' outlier samples, heavy-tailed biological noise, and a batch that confounds
#' condition and shifts 200 genes (the model-based methods adjust for it; a
#' rank test cannot). Compared against
#' ambient_de without outlier down-weighting, naive edgeR (robust QL), a naive
#' Wilcoxon rank-sum test, and a Wilcoxon test on ambient-subtracted counts.
#' Pass requires of ambient_de, averaged over seeds, in every case: at most 1
#' leaked gene called, null false-positive rate <= 0.075, false discovery
#' proportion <= 0.15, and power no more than 0.1 below the best comparator
#' that itself passes the first three criteria.
#'
#' @param seeds replicate seeds.
#' @param n_per samples per condition.
#' @param verbose print the table.
#' @return invisibly, a data.frame with attribute `"pass"`.
#' @export
gate_ambient_de_cohort <- function(seeds = 1:4, n_per = 50, verbose = TRUE) {
  for (p in c("edgeR", "limma"))
    if (!requireNamespace(p, quietly = TRUE)) stop("gate_ambient_de_cohort() needs ", p, call. = FALSE)
  cases <- list(clean = list(), outliers = list(out_p = 0.01), heavy_tail = list(heavy_df = 3),
                batch = list(n_batch = 200))
  rows <- list()
  for (nm in names(cases)) for (s in seeds) {
    st <- do.call(simulate_ambient_cohort, c(cases[[nm]], list(n_per = n_per, seed = s)))
    rows[[length(rows) + 1]] <- cbind(case = nm, suppressWarnings(.cohort_scores(st)))
  }
  df <- do.call(rbind, rows)
  agg <- stats::aggregate(cbind(leak_called, null_fp, power, fdp) ~ case + method, df, mean)
  agg$valid <- agg$leak_called <= 1 & agg$null_fp <= 0.075 & agg$fdp <= 0.15
  a <- agg[agg$method == "ambient_de", ]
  best <- vapply(a$case, function(cs) {
    o <- agg[agg$case == cs & agg$method != "ambient_de" & agg$valid, ]
    if (nrow(o)) max(o$power) else 0
  }, 1)
  a$ok <- a$valid & a$power >= best - 0.1
  pass <- all(a$ok)
  if (verbose) {
    agg <- agg[order(agg$case, agg$method), ]
    cat(sprintf("  %-11s %-21s | %7s %8s %6s %6s %6s\n", "case", "method", "leaked", "null_fp",
                "power", "FDP", "valid"))
    for (i in seq_len(nrow(agg))) with(agg[i, ], cat(sprintf(
      "  %-11s %-21s | %7.1f %8.3f %6.2f %6.2f %6s\n", case, method, leak_called, null_fp, power,
      fdp, if (valid) "yes" else "NO")))
    cat(sprintf("  => %s\n", if (pass) "PASS: calibrated under outliers and heavy tails, power kept"
                             else "FAIL: see ambient_de rows"))
  }
  attr(agg, "pass") <- pass
  invisible(agg)
}
