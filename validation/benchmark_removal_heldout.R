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
## Two sets. "v3": hgmm 1k, 5k, 10k v3 (seen while bayes_absent was designed:
## 10k was diagnosed, 1k and 5k results were looked at). "fresh": three 10x
## mixtures never used for any design decision, with other chemistries; set
## DECANT_HELDOUT=fresh to run them. They are the judge of bayes_absent.
SETS <- list(
  v3 = c("hgmm_1k_v3", "hgmm_5k_v3", "hgmm_10k_v3"),
  fresh = c("hgmm_6k", "hgmm_12k", "10k_hgmm_3p_nextgem_Chromium_X"),
  ## unused until the hybrid cap was frozen (preregistration/bayes_absent_hybrid.md)
  unused = c("hgmm_1k", "hgmm_1k_v2", "hgmm_100"))
URLS <- c(hgmm_6k = "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_6k/hgmm_6k_%s_gene_bc_matrices.tar.gz",
          hgmm_12k = "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_12k/hgmm_12k_%s_gene_bc_matrices.tar.gz",
          `10k_hgmm_3p_nextgem_Chromium_X` = "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_hgmm_3p_nextgem_Chromium_X/10k_hgmm_3p_nextgem_Chromium_X_%s_feature_bc_matrix.tar.gz",
          hgmm_1k = "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_1k/hgmm_1k_%s_gene_bc_matrices.tar.gz",
          hgmm_100 = "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_100/hgmm_100_%s_gene_bc_matrices.tar.gz",
          hgmm_1k_v2 = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_1k_v2/hgmm_1k_v2_%s_feature_bc_matrix.tar.gz")
SET <- Sys.getenv("DECANT_HELDOUT", "v3")
DATASETS <- SETS[[SET]]
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
## counts of one matrix kind ("raw" or "filtered") for any dataset. Cell
## Ranger 2 writes one matrix per genome (same barcodes); they are stacked.
read_kind <- function(ds, which) {
  if (!ds %in% names(URLS)) return(read_10x_counts(mtx_dir(ds, which)))
  dest <- file.path(CACHE, ds, which)
  if (!dir.exists(dest)) {
    tgz <- file.path(CACHE, paste0(ds, "_", which, ".tar.gz"))
    if (!file.exists(tgz)) utils::download.file(sprintf(URLS[[ds]], which), tgz, mode = "wb", quiet = TRUE)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE); utils::untar(tgz, exdir = dest)
  }
  dd <- list.dirs(dest); has <- vapply(dd, function(x) any(grepl("^matrix.mtx", list.files(x))), TRUE)
  ms <- lapply(sort(dd[has]), read_10x_counts)
  if (length(ms) == 1) return(ms[[1]])
  ## same barcodes in every genome's matrix: reindex columns through the
  ## triplet form (a sparse sub-assignment over ~700k barcodes exhausts memory)
  bc <- Reduce(union, lapply(ms, colnames))
  do.call(rbind, lapply(ms, function(m) {
    if (identical(colnames(m), bc)) return(m)
    t <- methods::as(m, "TsparseMatrix")
    Matrix::sparseMatrix(i = t@i + 1L, j = match(colnames(m), bc)[t@j + 1L], x = t@x,
                         dims = c(nrow(m), length(bc)), dimnames = list(rownames(m), bc))
  }))
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
  dir.create(file.path(CACHE, ds), recursive = TRUE, showWarnings = FALSE)
  raw <- read_kind(ds, "raw"); filt <- read_kind(ds, "filtered")
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
  ## single-sample decant() reads the empties only through their per-gene
  ## totals (ambient_global), so keep just those: identical output, and a
  ## 5k-cell channel with 1.1M empty droplets fits in memory
  x$empties <- Matrix::Matrix(Matrix::rowSums(x$empties), ncol = 1, sparse = TRUE,
                              dimnames = list(rownames(x$cells), NULL))
  invisible(gc())
  t0 <- Sys.time()
  res <- suppressWarnings(do.call(decant, c(list(x$cells, x$empties, clusters = x$clusters, correction = CORR),
                                             if (!is.null(CAPT)) list(cap_test = CAPT))))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cm <- res$corrected; dimnames(cm) <- dimnames(x$cells)
  stopifnot(max(cm - x$cells) <= 1e-9)   # never adds counts (sparse check)
  sc <- scores(x$cells, cm, x$gene_sp, x$cell_sp)
  rows[[ds]] <- rbind(comp, data.frame(dataset = ds, method = paste0("Decant (", TAG, ")"),
                                       retain = sc[["retain"]], removal = sc[["removal"]], seconds = secs,
                                       note = paste(CORR, CAPT %||% "")))
  print(rows[[ds]], digits = 3, row.names = FALSE)
  rm(x, res, cm); invisible(gc())
}
out <- do.call(rbind, rows)
utils::write.csv(out, file.path(OUT, sprintf("removal_heldout%s_%s.csv", if (SET == "v3") "" else paste0("_", SET), TAG)), row.names = FALSE)
cat("\nMean over held-out mixtures:\n")
print(stats::aggregate(cbind(retain, removal) ~ method, out, mean), digits = 3, row.names = FALSE)
