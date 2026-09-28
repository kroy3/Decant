#!/usr/bin/env Rscript
## validation/species_mix.R
## Experimental-ground-truth validation on public 10x human/mouse mixtures.
## Run from the repo root:  Rscript validation/species_mix.R [dataset ...]
##
## Comparators are the OFFICIAL packages, used only if installed:
##   SoupX (CRAN)             install.packages("SoupX")
##   DecontX (Bioconductor)   BiocManager::install("celda")
## A missing comparator is reported as "not run", never replaced by the
## in-repo reimplementations (those are not fair stand-ins).
##
## Two tests per dataset:
##   full    all cells, both genomes. Comparable to published benchmarks, but
##           EASY for absent-gene estimators (the other genome is a perfect
##           absent set).
##   hidden  majority-species cells, that species' genes only. The estimator
##           must work as on single-species data; truth still comes from the
##           cross-species reads it cannot see. This is the test that matters.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)

DATASETS <- list(
  hgmm_1k_v3  = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_1k_v3/hgmm_1k_v3",
  hgmm_5k_v3  = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_5k_v3/hgmm_5k_v3",
  hgmm_10k_v3 = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_10k_v3/hgmm_10k_v3"
)
CACHE <- file.path("validation", "data")
OUT   <- file.path("validation", "results")

fetch <- function(name, which) {
  dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
  tgz <- file.path(CACHE, sprintf("%s_%s_feature_bc_matrix.tar.gz", name, which))
  dest <- file.path(CACHE, name, which)
  if (!dir.exists(dest)) {
    if (!file.exists(tgz))
      utils::download.file(sprintf("%s_%s_feature_bc_matrix.tar.gz", DATASETS[[name]], which),
                           tgz, mode = "wb", quiet = TRUE)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    utils::untar(tgz, exdir = dest)
  }
  list.dirs(dest)[grepl("feature_bc_matrix$", list.dirs(dest))][1]
}

## ---- methods: each returns list(corrected, rho) or NULL if unavailable ----
run_decant <- function(cells, empties, clusters) {
  res <- suppressWarnings(decant(cells, empties, clusters = clusters))
  list(corrected = res$corrected, rho = res$rho,
       note = if (isTRUE(res$soup_mismatch)) "soup_mismatch warning" else "")
}

run_soupx <- function(cells, empties, clusters) {
  if (!requireNamespace("SoupX", quietly = TRUE)) return(NULL)
  tryCatch({
    tod <- cbind(cells, empties)
    sc <- SoupX::SoupChannel(tod, cells, calcSoupProfile = FALSE)
    sc <- SoupX::estimateSoup(sc, soupRange = c(1, 100))
    sc <- SoupX::setClusters(sc, stats::setNames(as.character(clusters), colnames(cells)))
    sc <- suppressWarnings(SoupX::autoEstCont(sc, doPlot = FALSE, verbose = FALSE))
    list(corrected = SoupX::adjustCounts(sc, verbose = 0),
         rho = sc$metaData$rho, note = "")
  }, error = function(e) list(error = conditionMessage(e)))
}

run_decontx <- function(cells, empties, clusters) {
  if (!requireNamespace("celda", quietly = TRUE)) return(NULL)
  tryCatch({
    d <- celda::decontX(x = cells, z = clusters, background = empties, verbose = FALSE)
    list(corrected = d$decontXcounts, rho = d$contamination, note = "")
  }, error = function(e) list(error = conditionMessage(e)))
}

METHODS <- list(Decant = run_decant, SoupX = run_soupx, DecontX = run_decontx)

evaluate <- function(label, cells, empties, species, truth_rho = NULL, full = TRUE) {
  clusters <- quick_labels(cells, k = 8)
  truth <- if (full) species_mix_truth(cells, empties, species) else NULL
  rows <- list()
  for (m in names(METHODS)) {
    t0 <- Sys.time()
    r <- METHODS[[m]](cells, empties, clusters)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(r)) { rows[[m]] <- data.frame(test = label, method = m, status = "not installed"); next }
    if (!is.null(r$error)) { rows[[m]] <- data.frame(test = label, method = m, status = paste("error:", r$error)); next }
    sc <- if (full) {
      score_species_mix(cells, empties, species, r$corrected, r$rho, truth = truth)
    } else {
      rr <- as.numeric(r$rho)
      data.frame(n_cells = length(truth_rho), rho_true_mean = mean(truth_rho),
                 rho_hat_mean = mean(rr), rho_bias = mean(rr) - mean(truth_rho),
                 rho_rmse = sqrt(mean((rr - truth_rho)^2)),
                 rho_cor = if (stats::sd(rr) > 0) stats::cor(rr, truth_rho) else NA_real_)
    }
    rows[[m]] <- cbind(data.frame(test = label, method = m, status = "ok",
                                  seconds = round(secs, 1), note = r$note), sc)
  }
  do.call(rbind, lapply(rows, function(x) { x[setdiff(ALLCOLS, names(x))] <- NA; x[ALLCOLS] }))
}
ALLCOLS <- c("test", "method", "status", "seconds", "note", "n_cells", "rho_true_mean",
             "rho_hat_mean", "rho_bias", "rho_rmse", "rho_cor", "cross_removed",
             "own_removed_ratio")

args <- commandArgs(trailingOnly = TRUE)
todo <- if (length(args)) args else names(DATASETS)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
all <- list()
for (ds in todo) {
  cat("==", ds, "==\n")
  raw <- read_10x_counts(fetch(ds, "raw"))
  filt <- read_10x_counts(fetch(ds, "filtered"))
  d <- split_droplets(raw, cells = colnames(filt))
  sp <- gene_species(rownames(d$cells))
  tr <- species_mix_truth(d$cells, d$empties, sp)
  cat(sprintf("  %d cells (%d doublets excluded), %d empties; soup share: %s\n",
              nrow(tr), sum(tr$doublet), ncol(d$empties),
              paste(sprintf("%s=%.2f", names(attr(tr, "soup_share")), attr(tr, "soup_share")), collapse = " ")))
  keep <- !tr$doublet
  res_full <- evaluate("full", d$cells[, keep], d$empties, sp)
  h <- species_mix_hidden(d$cells, d$empties, sp, truth = tr)
  res_hid <- evaluate(paste0("hidden (", h$species, " only)"), h$cells, h$empties, NULL,
                      truth_rho = h$rho_true, full = FALSE)
  res <- cbind(dataset = ds, rbind(res_full, res_hid))
  print(res[, c("test", "method", "status", "rho_true_mean", "rho_hat_mean", "rho_rmse",
                "rho_cor", "cross_removed", "own_removed_ratio")], digits = 3, row.names = FALSE)
  all[[ds]] <- res
}
res <- do.call(rbind, all)
utils::write.csv(res, file.path(OUT, "species_mix.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT, "species_mix.csv"), "\n")
