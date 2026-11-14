#!/usr/bin/env Rscript
## experiments/large_cohort_calibration/run.R
## Sample-size sweep behind gate_ambient_de_cohort(): does contamination-aware
## DE stay calibrated as cohorts grow and counts stop being negative binomial?
## Run from the repo root:  Rscript experiments/large_cohort_calibration/run.R
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)

cases <- list(clean = list(), outliers = list(out_p = 0.01), heavy_tail = list(heavy_df = 3))
rows <- list()
for (n in c(6, 15, 50, 100)) for (nm in names(cases)) for (seed in 1:3) {
  st <- do.call(simulate_ambient_cohort, c(cases[[nm]], list(n_per = n, seed = seed)))
  rows[[length(rows) + 1]] <- cbind(n_per = n, case = nm, seed = seed,
                                    suppressWarnings(Decant:::.cohort_scores(st)))
  message(n, " ", nm, " ", seed)
}
df <- do.call(rbind, rows)
agg <- stats::aggregate(cbind(leak_called, null_fp, power, fdp) ~ method + case + n_per, df, mean)
out <- file.path("experiments", "large_cohort_calibration")
utils::write.csv(df, file.path(out, "results_per_seed.csv"), row.names = FALSE)
utils::write.csv(agg, file.path(out, "results.csv"), row.names = FALSE)
print(agg, digits = 3, row.names = FALSE)
