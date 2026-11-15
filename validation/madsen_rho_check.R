#!/usr/bin/env Rscript
## validation/madsen_rho_check.R
## On the low-complexity human/mouse mixture, compare each tool's per-cell
## contamination estimate with the contamination implied by cross-species
## reads. For a human cell, mouse reads can only be ambient, so
##   rho_true = (mouse fraction of the cell's UMIs) / (mouse fraction of the soup)
## and symmetrically for mouse cells. This is an estimate of the truth (it
## assumes the soup is the empty-droplet profile), but it uses no tool's output.
## Needs the cache written by `Rscript validation/madsen_benchmark.R low_complexity`.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))
source(file.path("validation", "comparators.R"))
p <- readRDS(file.path("validation", "data", "madsen", "low_complexity_prep.rds"))
g <- ifelse(grepl("^GRCh38", rownames(p$cells)), "human", "mouse")
soup <- Matrix::rowSums(p$empties); soup_h <- sum(soup[g == "human"]) / sum(soup)
tot <- Matrix::colSums(p$cells)
other <- ifelse(p$sp == "human", Matrix::colSums(p$cells[g == "mouse", ]), Matrix::colSums(p$cells[g == "human", ])) / tot
rho_true <- pmin(1, other / ifelse(p$sp == "human", 1 - soup_h, soup_h))

rows <- list(data.frame(method = "cross-species (truth)", rho_median = stats::median(rho_true),
                        rho_mean = mean(rho_true), spearman_with_truth = 1))
for (nm in c("Decant", "SoupX", "DecontX")) {
  invisible(gc())
  r <- COMPARATORS[[nm]](p$cells, p$empties, p$cl)
  rr <- rep_len(as.numeric(r$rho), ncol(p$cells))   # SoupX reports one value per channel
  rows[[length(rows) + 1]] <- data.frame(method = nm, rho_median = stats::median(rr), rho_mean = mean(rr),
    spearman_with_truth = if (stats::sd(rr) > 0) stats::cor(rr, rho_true, method = "spearman") else NA)
  rm(r)
}
res <- do.call(rbind, rows)
cat(sprintf("soup: %.1f%% human UMIs\n", 100 * soup_h)); print(res, digits = 3, row.names = FALSE)
utils::write.csv(res, file.path("validation", "results", "madsen_low_complexity_rho.csv"), row.names = FALSE)
