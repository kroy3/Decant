#!/usr/bin/env Rscript
## validation/benchmark_removal_heldout.R
## Retain and removal scores (Cargnelli et al. 2026, the same per-gene scorer
## as benchmark_species_mixture.R) on the 10x hgmm 1k, 5k and 10k v3
## mixtures, which are held out from the design of Decant's correction.
## Run from the repo root:
##   Rscript validation/benchmark_removal_heldout.R [tag] [decant correction] [cap_test]
## `tag` names the Decant variant in the output (default "current").
## SoupX and DecontX results are cached per dataset (they do not depend on
## Decant); delete validation/data/<dataset>/removal_comparators.rds to rerun.
##
## Inputs per dataset, built once and cached (removal_inputs.rds): cells
## called by Cell Ranger, cross-species doublets excluded; empties 1 to 100
## UMIs; clusters by the same Seurat recipe as the 20k mixture (resolution 0.1).

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))
source(file.path("validation", "comparators.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

src <- readLines(file.path("validation", "benchmark_species_mixture.R"))
eval(parse(text = src[grep("^seurat_clusters <- function", src):(grep("^run_tools <- function", src) - 1)]))

args <- commandArgs(trailingOnly = TRUE)
TAG <- if (length(args) >= 1) args[1] else "current"
CORR <- if (length(args) >= 2) args[2] else "bayes_cell"
CAPT <- if (length(args) >= 3) args[3] else NULL
DATASETS <- c("hgmm_1k_v3", "hgmm_5k_v3", "hgmm_10k_v3")
CACHE <- file.path("validation", "data"); OUT <- file.path("validation", "results")

mtx_dir <- function(ds, which) {
  d <- list.dirs(file.path(CACHE, ds, which)); d <- d[grepl("feature_bc_matrix$", d)][1]
  if (is.na(d)) {
    utils::untar(file.path(CACHE, sprintf("%s_%s_feature_bc_matrix.tar.gz", ds, which)),
                 exdir = file.path(CACHE, ds, which))
    d <- list.dirs(file.path(CACHE, ds, which)); d <- d[grepl("feature_bc_matrix$", d)][1]
  }
  d
}

## per-gene retain/removal, per species, averaged over species (as the paper)
scores <- function(cells, corrected, gene_sp, cell_sp) {
  res <- sapply(unique(cell_sp), function(s) {
    cc <- cell_sp == s
    x_end <- Matrix::rowSums(cells[gene_sp == s, cc, drop = FALSE]); y_end <- Matrix::rowSums(corrected[gene_sp == s, cc, drop = FALSE])
    x_amb <- Matrix::rowSums(cells[gene_sp != s, cc, drop = FALSE]); y_amb <- Matrix::rowSums(corrected[gene_sp != s, cc, drop = FALSE])
    k <- x_end > 0; a <- x_amb > 0
    c(retain = mean(pmax(0, 1 - abs(y_end[k] - x_end[k]) / x_end[k])),
      removal = mean(pmin(1, abs(x_amb[a] - y_amb[a]) / x_amb[a])))
  })
  rowMeans(res)
}

inputs <- function(ds) {
  f <- file.path(CACHE, ds, "removal_inputs.rds")
  if (file.exists(f)) return(readRDS(f))
  raw <- read_10x_counts(mtx_dir(ds, "raw")); filt <- read_10x_counts(mtx_dir(ds, "filtered"))
  d <- split_droplets(raw, cells = colnames(filt)); rm(raw, filt); invisible(gc())
  gsp <- gene_species(rownames(d$cells))
  tr <- species_mix_truth(d$cells, d$empties, gsp)
  keep <- !tr$doublet
  cells <- d$cells[, keep]
  x <- list(cells = cells, empties = d$empties, gene_sp = gsp, cell_sp = tr$species[keep],
            clusters = seurat_clusters(cells, resolution = 0.1))
  saveRDS(x, f); x
}

rows <- list()
for (ds in DATASETS) {
  x <- inputs(ds)
  cat(sprintf("== %s: %d cells, %d empties, %d clusters\n", ds, ncol(x$cells), ncol(x$empties), length(unique(x$clusters))))
  cf <- file.path(CACHE, ds, "removal_comparators.rds")
  comp <- if (file.exists(cf)) readRDS(cf) else {
    cmp <- list()
    for (m in c("SoupX", "DecontX")) {
      t0 <- Sys.time()
      r <- (if (m == "SoupX") run_soupx else run_decontx)(x$cells, x$empties, x$clusters)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      sc <- if (!is.null(r$corrected)) { cm <- methods::as(r$corrected, "CsparseMatrix"); dimnames(cm) <- dimnames(x$cells)
                                         scores(x$cells, cm, x$gene_sp, x$cell_sp) } else c(retain = NA, removal = NA)
      cmp[[m]] <- data.frame(dataset = ds, method = m, retain = sc[["retain"]], removal = sc[["removal"]],
                             seconds = secs, note = r$note %||% r$error %||% "not installed")
      invisible(gc())
    }
    cmp <- do.call(rbind, cmp); saveRDS(cmp, cf); cmp
  }
  t0 <- Sys.time()
  res <- suppressWarnings(do.call(decant, c(list(x$cells, x$empties, clusters = x$clusters, correction = CORR),
                                             if (!is.null(CAPT)) list(cap_test = CAPT))))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cm <- res$corrected; dimnames(cm) <- dimnames(x$cells)
  stopifnot(all(cm <= x$cells))   # never adds counts
  sc <- scores(x$cells, cm, x$gene_sp, x$cell_sp)
  rows[[ds]] <- rbind(comp, data.frame(dataset = ds, method = paste0("Decant (", TAG, ")"),
                                       retain = sc[["retain"]], removal = sc[["removal"]], seconds = secs,
                                       note = paste(CORR, CAPT %||% "")))
  print(rows[[ds]], digits = 3, row.names = FALSE)
  rm(x, res, cm); invisible(gc())
}
out <- do.call(rbind, rows)
utils::write.csv(out, file.path(OUT, sprintf("removal_heldout_%s.csv", TAG)), row.names = FALSE)
cat("\nMean over held-out mixtures:\n")
print(stats::aggregate(cbind(retain, removal) ~ method, out, mean), digits = 3, row.names = FALSE)
