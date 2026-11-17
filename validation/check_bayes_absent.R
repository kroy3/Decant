#!/usr/bin/env Rscript
## validation/check_bayes_absent.R
## The pre-stated checks for correction = "bayes_absent" besides the fresh
## held-out mixtures (benchmark_removal_heldout.R, DECANT_HELDOUT=fresh):
##   1. simulation gate: gate_correction()'s settings, bayes_cell vs bayes_absent;
##   2. PBMC 4k: marker counts that are and are not ambient-compatible, and
##      native markers (the data and metrics of benchmark_pbmc_six_tools.R);
##   3. Smart-seq2 negative control (the data and score of
##      benchmark_negative_control.R): retain must stay >= 0.99.
## The data loading is taken verbatim from those scripts, so the inputs are
## identical. Results: validation/results/bayes_absent_checks_*.csv
## (BENCHMARK.md section 1g). Run from the repo root.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
OUT_DIR <- file.path("validation", "results")
RULES <- c("bayes_cell", "bayes_absent")
`%||%` <- function(a, b) if (is.null(a)) b else a
lines_between <- function(file, from, to) {
  s <- readLines(file); a <- grep(from, s, fixed = TRUE)[1]; b <- grep(to, s, fixed = TRUE)[1]
  stopifnot(!is.na(a), !is.na(b), b > a); s[a:(b - 1)]
}

## ---- 1. simulation gate ----
rows <- list()
for (setting in names(.gate_settings())) for (rm in c(0.05, 0.2)) for (s in 1:3) {
  sim <- .gate_sim(setting, rm, s)
  soup <- ambient_global(sim$empty); su <- sum(.col_sums(sim$empty))
  cl <- quick_labels(sim$observed, k = length(sim$lysis_true), seed = s)
  rho <- as.numeric(estimate_rho_cluster(sim$observed, soup, cl))
  for (m in RULES) {
    sc <- score_correction(sim$observed, correct_counts(sim$observed, soup, rho, method = m, clusters = cl,
                           soup_umis = if (m == "bayes_absent") su), sim$truth)
    rows[[length(rows) + 1]] <- data.frame(setting = setting, rho = rm, seed = s, method = m,
      l1_error = sc$l1_error, sensitivity = sc$sensitivity,
      signal_destroyed = sc$signal_destroyed_frac, fabricated = sc$counts_fabricated)
  }
}
gate <- stats::aggregate(cbind(l1_error, sensitivity, signal_destroyed, fabricated) ~ setting + rho + method,
                         do.call(rbind, rows), mean)
cat("== simulation gate\n"); print(gate, digits = 3, row.names = FALSE)

## ---- 2. PBMC 4k ----
source(file.path("validation", "comparators.R"))
mode <- "score"
eval(parse(text = lines_between(file.path("validation", "benchmark_pbmc_six_tools.R"),
                                "CACHE <- file.path(", "if (identical(mode, \"export\"))")))
soup <- soup_counts / sum(soup_counts); names(soup) <- rownames(cells)
T_c <- Matrix::colSums(cells); AMBIENT_BOUND <- 0.15
sets <- list(immunoglobulin = list(genes = c("IGKC", "IGHM", "IGLC2"), home = "B"),
             monocyte = list(genes = c("LYZ", "S100A8", "S100A9"), home = "MNP"),
             nk_cytotoxic = list(genes = c("GNLY", "NKG7"), home = c("NK", "T_CD8")))
pbmc <- do.call(rbind, lapply(RULES, function(m) {
  corr <- run_decant(cells, empties, md$Cluster, m)$corrected
  do.call(rbind, lapply(names(sets), function(s) {
    g <- sets[[s]]$genes; home <- md$Annotation %in% sets[[s]]$home
    X <- as.matrix(cells[g, !home]); C <- as.matrix(corr[g, !home])
    inc <- X > 0 & stats::ppois(X - 1, outer(soup[g], AMBIENT_BOUND * T_c[!home]), lower.tail = FALSE) < 0.001
    data.frame(method = m, set = s,
               compatible_removed = 1 - sum(C[!inc]) / sum(X[!inc]),
               incompatible_removed = 1 - sum(C[inc]) / sum(X[inc]),
               native_removed = 1 - sum(corr[g, home]) / sum(cells[g, home]),
               total_removed = 1 - sum(corr) / sum(cells),
               counts_added = sum(pmax(corr - cells, 0)))
  }))
}))
cat("\n== PBMC 4k\n"); print(pbmc, digits = 3, row.names = FALSE)
rm(cells, empties, md); invisible(gc())

## ---- 3. Smart-seq2 negative control ----
nc_file <- file.path("validation", "benchmark_negative_control.R")
nc <- readLines(nc_file)
a <- grep('"Decant" = function', nc, fixed = TRUE); b <- grep('"DecontX" = function', nc, fixed = TRUE)
e <- grep("^res <- do.call\\(rbind, rows\\)", nc)
stopifnot(length(a) == 1, length(b) == 1, length(e) == 1, a < b, b < e)
nc <- c(nc[seq_len(a - 1)],
        '    "bayes_cell" = function() run_decant(cells, soup, cl, "bayes_cell"),',
        '    "bayes_absent" = function() run_decant(cells, soup, cl, "bayes_absent"))',
        nc[(b + 1):e])
eval(parse(text = nc))
negctrl <- stats::aggregate(cbind(retain, removed) ~ method, res, mean)
cat("\n== Smart-seq2 negative control (mean over donors)\n"); print(negctrl, digits = 3, row.names = FALSE)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(gate, file.path(OUT_DIR, "bayes_absent_checks_gate.csv"), row.names = FALSE)
utils::write.csv(pbmc, file.path(OUT_DIR, "bayes_absent_checks_pbmc.csv"), row.names = FALSE)
utils::write.csv(res, file.path(OUT_DIR, "bayes_absent_checks_negative_control.csv"), row.names = FALSE)
