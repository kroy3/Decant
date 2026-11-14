#!/usr/bin/env Rscript
## validation/benchmark_soupx.R
## Decant against the official SoupX on (A) real PBMC data without ground
## truth and (B) simulated data with ground truth.
## Run from the repo root:  Rscript validation/benchmark_soupx.R
##
## Both tools get the same cells, the same soup or empties, and the same
## clusters. SoupX is run with its defaults (autoEstCont, adjustCounts); if
## autoEstCont finds no markers, its own advice (soupQuantile = 0.5) is
## followed once and noted. Caveat for (B): the simulator follows the additive
## ambient model both tools assume, and was written alongside Decant.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))
if (!requireNamespace("SoupX", quietly = TRUE)) stop("this script needs SoupX", call. = FALSE)
source(file.path("validation", "comparators.R"))
OUT <- file.path("validation", "results"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------------- (A) real PBMC 4k ----------------
CACHE <- file.path("validation", "data", "SoupX")
if (!dir.exists(CACHE))
  system2("git", c("clone", "-q", "--depth", "1", "https://github.com/constantAmateur/SoupX", CACHE))
e <- new.env()
load(file.path(CACHE, "data", "PBMC_sc.RData"), e); load(file.path(CACHE, "data", "PBMC_metaData.RData"), e)
sc <- e$PBMC_sc; md <- e$PBMC_metaData[colnames(e$PBMC_sc$toc), ]
ok <- md$Annotation != "?"
cells <- methods::as(methods::as(sc$toc, "CsparseMatrix"), "dgCMatrix")[, ok]; md <- md[ok, ]
empties <- Matrix::Matrix(matrix(sc$soupProfile$counts, ncol = 1,
                                 dimnames = list(rownames(cells), "pooled")), sparse = TRUE)
soup <- empties[, 1] / sum(empties); T_c <- Matrix::colSums(cells)

t_dec <- system.time(dec <- suppressWarnings(decant(cells, empties, clusters = md$Cluster)))[["elapsed"]]
## SoupX on its own SoupChannel (its soup profile is the same counts)
sx <- sc; sx$toc <- sx$toc[, ok]; sx$metaData <- sx$metaData[colnames(cells), , drop = FALSE]
sx <- SoupX::setClusters(sx, stats::setNames(as.character(md$Cluster), colnames(cells)))
t_sx <- system.time({
  sx <- suppressWarnings(SoupX::autoEstCont(sx, doPlot = FALSE, verbose = FALSE))
  sx_corr <- SoupX::adjustCounts(sx, verbose = 0)
})[["elapsed"]]
sx_corr <- methods::as(methods::as(sx_corr, "CsparseMatrix"), "dgCMatrix")

sets <- list(immunoglobulin = list(genes = c("IGKC", "IGHM", "IGLC2"), home = "B"),
             monocyte = list(genes = c("LYZ", "S100A8", "S100A9"), home = "MNP"),
             nk_cytotoxic = list(genes = c("GNLY", "NKG7"), home = c("NK", "T_CD8")))
real <- function(name, corr, rho, secs) {
  row <- data.frame(method = name, rho_median = stats::median(rho), rho_mean = mean(rho),
                    total_removed = 1 - sum(corr) / sum(cells),
                    counts_added = sum(pmax(corr - cells, 0)),
                    non_integer = mean(abs(corr@x - round(corr@x)) > 1e-8), seconds = secs)
  for (nm in names(sets)) {
    g <- sets[[nm]]$genes; home <- md$Annotation %in% sets[[nm]]$home
    row[[paste0(nm, "_foreign_removed")]] <- 1 - sum(corr[g, !home]) / sum(cells[g, !home])
    row[[paste0(nm, "_native_removed")]] <- 1 - sum(corr[g, home]) / sum(cells[g, home])
  }
  row
}
tabA <- rbind(real("Decant", dec$corrected, dec$rho, t_dec),
              real("SoupX 1.6.2", sx_corr, sx$metaData$rho, t_sx))
cat("\n(A) real PBMC 4k, no ground truth\n"); print(t(tabA), digits = 3)
utils::write.csv(tabA, file.path(OUT, "benchmark_soupx_pbmc.csv"), row.names = FALSE)

## ---------------- (B) simulated, with ground truth ----------------
rows <- list()
for (setting in c("easy", "stress")) for (rm in c(0, 0.05, 0.1, 0.2)) for (seed in 1:3) {
  sim <- Decant:::.gate_sim(setting, rm, seed)
  cl <- quick_labels(sim$observed, k = length(sim$lysis_true), seed = seed)
  for (m in c("Decant", "SoupX")) {
    r <- if (m == "Decant") run_decant(sim$observed, sim$empty, cl) else run_soupx(sim$observed, sim$empty, cl)
    if (!is.null(r$error)) { rows[[length(rows) + 1]] <- data.frame(setting, rho = rm, seed, method = m,
        rho_est = NA, rho_rmse = NA, l1_error = NA, sensitivity = NA, signal_destroyed = NA, note = r$error); next }
    s <- score_correction(sim$observed, as.matrix(r$corrected), sim$truth)
    rows[[length(rows) + 1]] <- data.frame(setting, rho = rm, seed, method = m, rho_est = mean(r$rho),
      rho_rmse = sqrt(mean((r$rho - sim$rho_true)^2)), l1_error = s$l1_error,
      sensitivity = s$sensitivity, signal_destroyed = s$signal_destroyed_frac, note = r$note)
  }
}
dfB <- do.call(rbind, rows)
utils::write.csv(dfB, file.path(OUT, "benchmark_soupx_sim_per_seed.csv"), row.names = FALSE)
tabB <- stats::aggregate(cbind(rho_est, rho_rmse, l1_error, sensitivity, signal_destroyed) ~ method + rho + setting,
                         dfB, mean, na.action = stats::na.pass)
tabB$failed_runs <- stats::aggregate(is.na(l1_error) ~ method + rho + setting, dfB, sum)[[4]]
tabB$fallback_runs <- stats::aggregate(grepl("soupQuantile", note) ~ method + rho + setting, dfB, sum)[[4]]
cat("\n(B) simulated, mean of 3 seeds\n"); print(tabB, digits = 3, row.names = FALSE)
utils::write.csv(tabB, file.path(OUT, "benchmark_soupx_sim.csv"), row.names = FALSE)
