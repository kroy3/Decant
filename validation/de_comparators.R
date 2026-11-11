#!/usr/bin/env Rscript
## validation/de_comparators.R
## Does correcting counts with an existing tool, then testing, stop ambient
## leakage into differential expression? Compared against modelling the soup
## inside the test (ambient_de).
##
## Run from the repo root:
##   Rscript validation/de_comparators.R            # seeds 1 2
##   Rscript validation/de_comparators.R 1 2 3 4    # chosen seeds
##
## Arms, all testing the rare cell type R with the same scoring as
## gate_ambient_de():
##   naive_edgeR         edgeR QL on the raw pseudobulk
##   ambient_de          Decant's contamination-aware test
##   <Tool> + edgeR      each sample's cells corrected by <Tool> (official
##                       package, see validation/comparators.R), pseudobulked,
##                       then edgeR QL. A tool that is not installed is
##                       reported as "not installed", never substituted.
##
## Caveats that favour the comparators, stated so nobody has to find them:
##   * every corrector is given the TRUE cell-type labels as clusters;
##   * each sample is corrected with its own empty droplets.
## Caveat that favours ambient_de:
##   * the data come from simulate_ambient_study(), which is built around the
##     additive-soup mechanism the test assumes. This shows what correct-then-
##     test can and cannot do under that mechanism. It is not real-data
##     evidence against any tool.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
for (p in c("edgeR", "limma"))
  if (!requireNamespace(p, quietly = TRUE)) stop("this script needs ", p, call. = FALSE)
source(file.path("validation", "comparators.R"))

OUT <- file.path("validation", "results")
CASES <- list(base = list(), neg_control = list(rho_dis = 0.05),
              shared_true_de = list(n_shared = 20))
ALPHA <- 0.05; FDR <- 0.1
## Decant's non-default correction rule, as its own arm: the default
## `subtract` clamps at zero without redistributing, so it removes less than
## rho of each cell's counts.
COMPARATORS <- append(COMPARATORS, list(
  `Decant (redistribute)` = function(cells, empties, clusters)
    run_decant(cells, empties, clusters, correction = "redistribute")), after = 1)
norm_lib <- utils::getFromNamespace(".norm_lib_sizes", "Decant")

pseudobulk <- function(X, keep, sample) {
  S <- sort(unique(sample))
  Y <- vapply(S, function(i) Matrix::rowSums(X[, keep & sample == i, drop = FALSE]),
              numeric(nrow(X)))
  rownames(Y) <- rownames(X); Y
}

edger_ql <- function(Y, design) {
  k <- edgeR::filterByExpr(edgeR::DGEList(Y), design)
  d <- edgeR::estimateDisp(norm_lib(edgeR::DGEList(Y[k, , drop = FALSE])), design)
  stats::setNames(edgeR::glmQLFTest(edgeR::glmQLFit(d, design), coef = 2)$table$PValue,
                  rownames(Y)[k])
}

## identical definitions to .de_case_scores(): an untested gene counts as
## not called; null_fp is over tested genes that are neither true DE nor leak.
score <- function(pv, st, genes) {
  gi <- seq_along(genes)
  is_de <- gi %in% st$de_rare; is_leak <- gi %in% st$de_abund; is_sh <- gi %in% st$shared
  p <- stats::setNames(rep(1, length(genes)), genes); p[names(pv)] <- pv
  q <- p; q[] <- 1; q[names(pv)] <- stats::p.adjust(pv, "BH")
  tested <- genes %in% names(pv)
  data.frame(status = "ok", leak_called = sum(q[is_leak] < FDR),
             leak_total = sum(is_leak),
             null_fp = mean(p[tested & !is_de & !is_leak] < ALPHA),
             power = mean(q[is_de & !is_sh] < FDR),
             power_shared = if (any(is_sh)) mean(q[is_sh] < FDR) else NA_real_,
             fdp = if (any(q < FDR)) mean(!is_de[q < FDR]) else 0)
}
blank <- function(status) data.frame(status = status, leak_called = NA_real_, leak_total = NA_real_,
                                     null_fp = NA_real_, power = NA_real_,
                                     power_shared = NA_real_, fdp = NA_real_)

run_case <- function(st) {
  genes <- rownames(st$X); keep <- st$label == "R"
  design <- stats::model.matrix(~st$cond)
  Y <- pseudobulk(st$X, keep, st$sample)
  rows <- list(
    naive_edgeR = score(edger_ql(Y, design), st, genes),
    ambient_de = {
      r <- ambient_de(Y, st$empties, design, naive = FALSE)
      score(stats::setNames(r$PValue, r$gene), st, genes)
    })
  for (m in names(COMPARATORS)) {
    arm <- paste0(m, " + edgeR")
    S <- sort(unique(st$sample)); corrected <- vector("list", length(S)); status <- "ok"
    notes <- character(0)
    for (s in S) {
      idx <- st$sample == s
      r <- COMPARATORS[[m]](st$X[, idx, drop = FALSE], st$empties[[s]], st$label[idx])
      if (is.null(r)) { status <- "not installed"; break }
      if (!is.null(r$error)) { status <- paste("error:", r$error); break }
      if (!is.null(r$note) && nzchar(r$note)) notes <- union(notes, r$note)
      corrected[[s]] <- Matrix::rowSums(as.matrix(r$corrected)[, keep[idx], drop = FALSE])
    }
    rows[[arm]] <- if (status != "ok") blank(status) else {
      Yc <- do.call(cbind, corrected); rownames(Yc) <- genes
      cbind(score(edger_ql(Yc, design), st, genes), note = paste(notes, collapse = "; "))
    }
    if (is.null(rows[[arm]]$note)) rows[[arm]]$note <- ""
  }
  rows <- lapply(rows, function(r) { if (is.null(r$note)) r$note <- ""; r })
  cbind(method = names(rows), do.call(rbind, rows))
}

args <- commandArgs(trailingOnly = TRUE)
seeds <- if (length(args)) as.integer(args) else 1:2
all <- list()
for (nm in names(CASES)) for (s in seeds) {
  cat(sprintf("case %-15s seed %d\n", nm, s))
  st <- do.call(simulate_ambient_study, c(CASES[[nm]], list(seed = s)))
  all[[length(all) + 1]] <- cbind(case = nm, seed = s, run_case(st))
}
df <- do.call(rbind, all); rownames(df) <- NULL
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(df, file.path(OUT, "de_comparators_per_seed.csv"), row.names = FALSE)

## mean over seeds; an arm that did not run in every seed is reported with its status
num <- c("leak_called", "leak_total", "null_fp", "power", "power_shared", "fdp")
key <- unique(df[, c("case", "method")])
agg <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
  d <- df[df$case == key$case[i] & df$method == key$method[i], ]
  ok <- all(d$status == "ok")
  v <- if (ok) lapply(d[num], function(x) if (all(is.na(x))) NA_real_ else mean(x)) else
    stats::setNames(as.list(rep(NA_real_, length(num))), num)
  cbind(key[i, ], status = if (ok) "ok" else d$status[d$status != "ok"][1],
        n_seeds = nrow(d), as.data.frame(v),
        note = paste(setdiff(unique(d$note), ""), collapse = "; "))
}))
utils::write.csv(agg, file.path(OUT, "de_comparators.csv"), row.names = FALSE)
print(agg, digits = 3, row.names = FALSE)
cat("\nwrote", file.path(OUT, "de_comparators.csv"), "\n")
