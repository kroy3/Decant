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
  d0 <- edgeR::estimateDisp(edgeR::calcNormFactors(edgeR::DGEList(Y[kn, ])), design)
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
