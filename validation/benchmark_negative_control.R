#!/usr/bin/env Rscript
## validation/benchmark_negative_control.R
## Negative control from Cargnelli et al. (2026): Smart-seq2 human
## islets (Segerstolpe et al. 2016, ArrayExpress E-MTAB-5061). Each cell was
## sorted into its own well, so there is no ambient RNA to remove and the best
## possible correction is no correction. The score is the paper's endogenous
## retain score with every gene counted as endogenous: per donor, over genes
## with input count x_g > 0 (summed over the donor's cells),
##   retain = mean_g max(0, 1 - |y_g - x_g| / x_g)
## and then the mean over donors.
##
## What follows the paper (AmbientRemovalBenchmark, scripts 03-05):
##   * donors HP1506401, HP1507101, HP1508501T2D, each run separately;
##   * cells = wells with 100,000-500,000 counts and 2,500-10,500 detected genes;
##   * "empty" boundaries set explicitly at 1,000-100,000 counts: SoupX gets
##     soupRange = c(1000, 100000) over all wells, and autoEstCont(forceAccept = TRUE);
##   * DecontX gets the cells and, as background, every other well of the donor
##     (the paper passes the raw matrix; decontX drops the cell barcodes from it);
##   * clusters: Seurat on all three donors, Harmony by donor, 10 PCs, resolution 0.3.
## Decant gets the same soup as SoupX: the wells with 1,000-100,000 counts.
##
## What differs from the paper:
##   * counts are the authors' processed RefSeq count matrix from EBI
##     (pancreas_refseq_rpkms_counts_3514sc.txt, count columns), not the
##     paper's own STARsolo re-alignment. Thresholds are the paper's but are
##     applied to these counts, so the cell and soup well sets differ somewhat.
##   * Seurat defaults are used where the paper chose interactively.
##
## Also reported: Decant with the soup set to the pooled expression of the
## donor's cells (sum over cells, normalised), the profile SoupX "reduced"
## uses; it needs no non-cell wells at all.
##
## Usage (repo root): Rscript validation/benchmark_negative_control.R

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))
source(file.path("validation", "comparators.R"))
OUT <- file.path("validation", "results")
DATA <- file.path("validation", "data", "segerstolpe")
SAMPLES <- c("HP1506401", "HP1507101", "HP1508501T2D")
PUBLISHED <- data.frame(method = c("SoupX", "DecontX", "CellBender", "scAR"),
                        retain = c(0.988, 0.833, 0.797, 0.401))

f <- file.path(DATA, "pancreas_refseq_rpkms_counts_3514sc.txt")
if (!file.exists(f)) {
  dir.create(DATA, recursive = TRUE, showWarnings = FALSE)
  utils::download.file(paste0("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/061/E-MTAB-5061/Files/",
                              basename(f)), f, mode = "wb")
}
## Layout: header "#samples" + 3514 well names; rows: symbol, RefSeq id,
## 3514 RPKM columns, then 3514 count columns in the same order.
wells <- strsplit(readLines(f, n = 1), "\t")[[1]][-1]
nw <- length(wells)
tab <- data.table::fread(f, skip = 1, header = FALSE, sep = "\t", showProgress = FALSE)
stopifnot(ncol(tab) == 2 + 2 * nw)
keep_cols <- which(sub("_.*", "", wells) %in% SAMPLES)
counts <- as.matrix(tab[, 2 + nw + keep_cols, with = FALSE])
dimnames(counts) <- list(make.unique(tab[[1]]), wells[keep_cols])
rm(tab); invisible(gc())
frac <- mean(abs(counts - round(counts)) > 1e-8)
cat(sprintf("non-integer count entries: %.4f%% (rounded)\n", 100 * frac))
counts <- Matrix::Matrix(round(counts), sparse = TRUE)
donor <- sub("_.*", "", colnames(counts))

tot <- Matrix::colSums(counts); nfeat <- Matrix::colSums(counts > 0)
is_cell <- tot >= 1e5 & tot <= 5e5 & nfeat >= 2500 & nfeat <= 10500
is_soup <- !is_cell & tot >= 1000 & tot <= 1e5

## clusters on all donors' cells, Harmony by donor
cells_all <- counts[Matrix::rowSums(counts[, is_cell]) > 0, is_cell]
so <- suppressWarnings(Seurat::CreateSeuratObject(cells_all))
so$donor <- donor[is_cell]
so <- Seurat::NormalizeData(so, verbose = FALSE)
so <- Seurat::FindVariableFeatures(so, verbose = FALSE)
so <- Seurat::ScaleData(so, verbose = FALSE)
so <- Seurat::RunPCA(so, verbose = FALSE)
red <- "pca"
if (requireNamespace("harmony", quietly = TRUE)) {
  so <- suppressWarnings(harmony::RunHarmony(so, group.by.vars = "donor", verbose = FALSE)); red <- "harmony"
}
so <- Seurat::FindNeighbors(so, dims = 1:10, reduction = red, verbose = FALSE)
cl_all <- stats::setNames(as.character(Seurat::FindClusters(so, resolution = 0.3, verbose = FALSE)$seurat_clusters),
                          colnames(cells_all))
cat(sprintf("clusters: %d (reduction: %s)\n", length(unique(cl_all)), red))

retain_score <- function(x, y) {
  xs <- Matrix::rowSums(x); ys <- Matrix::rowSums(y); k <- xs > 0
  mean(pmax(0, 1 - abs(ys[k] - xs[k]) / xs[k]))
}

rows <- list()
for (s in SAMPLES) {
  in_s <- donor == s
  cells <- counts[, in_s & is_cell]; soup <- counts[, in_s & is_soup]
  bg <- counts[, in_s & !is_cell]; bg <- bg[, Matrix::colSums(bg) > 0, drop = FALSE]
  g <- Matrix::rowSums(cells) > 0
  cells <- cells[g, ]; soup <- soup[g, , drop = FALSE]; bg <- bg[g, , drop = FALSE]
  cl <- cl_all[colnames(cells)]
  cat(sprintf("%s: cells %d, soup wells (1k-100k) %d, other non-cell wells %d, clusters %d\n",
              s, ncol(cells), ncol(soup), ncol(bg), length(unique(cl))))
  pooled <- Matrix::Matrix(Matrix::rowSums(cells), ncol = 1, sparse = TRUE,
                           dimnames = list(rownames(cells), "pooled"))
  runs <- list(
    "Decant" = function() run_decant(cells, soup, cl),
    "Decant (soup = pooled cells)" = function() run_decant(cells, pooled, cl),
    "SoupX" = function() run_soupx(cells, soup, cl, soup_range = c(1000, 1e5), force_accept = TRUE),
    "DecontX" = function() run_decontx(cells, bg, cl))
  for (nm in names(runs)) {
    t0 <- Sys.time(); r <- runs[[nm]]()
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(r) || !is.null(r$error)) {
      rows[[length(rows) + 1]] <- data.frame(sample = s, method = nm, retain = NA, removed = NA,
        rho_mean = NA, rho_median = NA, seconds = secs, note = if (is.null(r)) "not installed" else r$error)
      next
    }
    corr <- methods::as(methods::as(r$corrected, "CsparseMatrix"), "dgCMatrix"); dimnames(corr) <- dimnames(cells)
    rho <- if (is.null(r$rho)) NA else r$rho
    rows[[length(rows) + 1]] <- data.frame(sample = s, method = nm, retain = retain_score(cells, corr),
      removed = 1 - sum(corr) / sum(cells), rho_mean = mean(rho), rho_median = stats::median(rho),
      seconds = secs, note = r$note %||% "")
  }
}
res <- do.call(rbind, rows)
cat("\nPer donor:\n"); print(res, digits = 3, row.names = FALSE)
summ <- stats::aggregate(cbind(retain, removed, rho_mean, rho_median) ~ method, res, mean, na.action = stats::na.pass)
summ <- summ[match(unique(res$method), summ$method), ]
cat("\nMean over donors:\n"); print(summ, digits = 3, row.names = FALSE)
cat("\nPublished (Cargnelli et al. 2026, Smart-seq2 retain score):\n"); print(PUBLISHED, row.names = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, file.path(OUT, "negative_control_scores_per_donor.csv"), row.names = FALSE)
utils::write.csv(summ, file.path(OUT, "negative_control_scores.csv"), row.names = FALSE)
