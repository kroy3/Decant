#!/usr/bin/env Rscript
## run_demo.R -- sources Decant directly (no install needed) and runs the
## benchmark. The point is the NEGATIVE CONTROL as much as the positive result:
##   soup_bias = 0  -> soup is uniform, the global assumption is exactly true,
##                     structured must NOT beat it (else the model is overfitting).
##   soup_bias > 0  -> soup is biased toward fragile types, global is wrong,
##                     structured should recover the true soup and correct better.

pkg <- file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), "Decant", "R")
if (!dir.exists(pkg)) pkg <- "Decant/R"
for (f in list.files(pkg, pattern = "\\.R$", full.names = TRUE)) source(f)

cat("Running benchmark (global soup vs structured soup) ...\n\n")
df  <- run_benchmark(soup_bias_grid = c(0, 1, 2, 4), reps = 3, rho_mean = 0.2)

agg <- aggregate(cbind(soup_err, sensitivity, signal_destroyed, preservation, fabricated) ~
                   method + soup_bias, data = df, FUN = mean)
agg <- agg[order(agg$soup_bias, agg$method), ]

fmt <- function(x) formatC(x, digits = 3, format = "f")
cat(sprintf("%-34s %5s | %8s %8s %10s %8s %6s\n",
            "method", "bias", "soupErr", "sens", "sigDestroy", "preserv", "fab"))
cat(strrep("-", 92), "\n")
for (i in seq_len(nrow(agg))) {
  cat(sprintf("%-34s %5.0f | %8s %8s %10s %8s %6.0f\n",
              agg$method[i], agg$soup_bias[i],
              fmt(agg$soup_err[i]), fmt(agg$sensitivity[i]),
              fmt(agg$signal_destroyed[i]), fmt(agg$preservation[i]),
              agg$fabricated[i]))
}

cat("\nReading the table:\n")
cat(" soupErr   = L1 error of the estimated soup profile vs truth (lower better)\n")
cat(" sens      = fraction of real contamination removed (higher better)\n")
cat(" sigDestroy= endogenous counts wrongly removed, as frac of true (lower better)\n")
cat(" preserv   = cosine sim of corrected vs clean cells (higher better)\n")
cat(" fab       = counts fabricated; mass-conservation guarantee => must be 0\n")
