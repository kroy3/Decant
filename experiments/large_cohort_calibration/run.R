#!/usr/bin/env Rscript
## experiments/large_cohort_calibration/run.R
## Sample-size sweep behind gate_ambient_de_cohort(): does contamination-aware
## DE stay calibrated as cohorts grow and counts stop being negative binomial?
## Run from the repo root:  Rscript experiments/large_cohort_calibration/run.R
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)

cases <- list(clean = list(), outliers = list(out_p = 0.01), heavy_tail = list(heavy_df = 3),
              batch = list(n_batch = 200))
## ambient_de against a minimum effect of 10% (TREAT-style), scored like the rest
lfc_scores <- function(st, alpha = 0.05, fdr = 0.1) {
  design <- if (length(st$batch_genes)) stats::model.matrix(~st$batch + st$cond)
            else stats::model.matrix(~st$cond)
  r <- ambient_de(st$Y, st$empties, design, naive = FALSE, lfc = log2(1.1))
  gi <- as.integer(sub("g", "", r$gene)); de <- gi %in% st$de_rare; lk <- gi %in% st$de_abund
  data.frame(method = "ambient_de_lfc1.1", leak_called = sum(r$FDR[lk] < fdr),
             null_fp = mean(r$PValue[!de & !lk] < alpha),
             power = sum(r$FDR[de] < fdr) / length(st$de_rare),
             fdp = if (any(r$FDR < fdr)) mean(!de[r$FDR < fdr]) else 0)
}
rows <- list()
for (n in c(6, 15, 50, 100, 200)) for (nm in names(cases)) for (seed in 1:3) {
  st <- do.call(simulate_ambient_cohort, c(cases[[nm]], list(n_per = n, seed = seed)))
  sc <- suppressWarnings(rbind(Decant:::.cohort_scores(st), lfc_scores(st)))
  rows[[length(rows) + 1]] <- cbind(n_per = n, case = nm, seed = seed, sc)
  message(n, " ", nm, " ", seed)
}
df <- do.call(rbind, rows)
agg <- stats::aggregate(cbind(leak_called, null_fp, power, fdp) ~ method + case + n_per, df, mean)
out <- file.path("experiments", "large_cohort_calibration")
utils::write.csv(df, file.path(out, "results_per_seed.csv"), row.names = FALSE)
utils::write.csv(agg, file.path(out, "results.csv"), row.names = FALSE)
print(agg, digits = 3, row.names = FALSE)
