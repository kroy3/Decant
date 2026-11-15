#!/usr/bin/env Rscript
## validation/madsen_benchmark.R
## Score Decant with the two metrics of Cargnelli, Nielsen & Madsen (2026),
## "Benchmarking computational decontamination of ambient RNA"
## (bioRxiv 10.64898/2026.01.13.699237; code: github.com/madsen-lab/AmbientRemovalBenchmark),
## so that Decant's numbers can be set beside their Supplementary Table S2.
##
## Their definitions (01_functions.R, calc_*_UMIs_per_feature), reproduced in
## madsen_scores() below. In a two-species mixture, for the cells of one
## species, "endogenous" genes are that species' genes and "ambient" genes the
## other species'. Per gene g with input count x_g > 0 (summed over the
## species' cells of one sample) and corrected count y_g:
##   endogenous retain score = mean_g max(0, 1 - |y_g - x_g| / x_g)
##   ambient removal score   = mean_g min(1, |x_g - y_g| / x_g)
## Means are unweighted over genes; human and mouse cells are scored
## separately and averaged.
##
## Usage (from the repo root):
##   Rscript validation/madsen_benchmark.R selftest
##       synthetic two-species mixture with known contamination; checks that
##       the scoring and the pipeline run, and prints scores for Decant, SoupX
##       and DecontX. NOT comparable to the paper.
##   Rscript validation/madsen_benchmark.R low_complexity
##       the paper's low-complexity species-mixing dataset. Needs
##       validation/data/madsen/20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5
##       (cf.10xgenomics.com; downloaded automatically when reachable).
##
## Preprocessing follows the paper for this dataset: emptyDrops (lower = 500,
## FDR <= 0.05) to call cells, cells with >= 75% of UMIs from one species
## kept, Seurat clusters at resolution 0.1. Their other QC steps (in
## 04_preprocessing_and_quality_control.qmd) are not reproduced, so numbers
## will differ somewhat from theirs even for the same tool.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))
source(file.path("validation", "comparators.R"))
OUT <- file.path("validation", "results")

## Their Supplementary Table S2, low-complexity species mixing column, for reference.
PUBLISHED_LOW <- data.frame(
  method = c("SoupX", "SoupX (reduced)", "FastCAR", "DecontX", "DecontX (reduced)", "scCDC",
             "CellBender", "scAR", "CellClear"),
  retain = c(0.988, 0.986, 1.000, 0.980, 0.990, 1.000, 0.996, 0.470, 0.876),
  removal = c(0.422, 0.493, 0.000, 0.624, 0.533, 0.002, 0.729, 0.959, 0.023))

genome_of <- function(genes) {
  g <- ifelse(grepl("^GRCh38", genes), "human", ifelse(grepl("^mm10", genes), "mouse", NA))
  if (anyNA(g)) stop("genes without a GRCh38 / mm10 prefix: ", paste(utils::head(genes[is.na(g)], 3), collapse = ", "))
  g
}

## per sample and species, then averaged over species, then over samples
madsen_scores <- function(input, corrected, cell_species, sample = rep(1, ncol(input))) {
  gs <- genome_of(rownames(input)); res <- NULL
  for (s in unique(sample)) for (sp in c("human", "mouse")) {
    cc <- sample == s & cell_species == sp
    if (sum(cc) < 2) next
    x_end <- Matrix::rowSums(input[gs == sp, cc, drop = FALSE]); y_end <- Matrix::rowSums(corrected[gs == sp, cc, drop = FALSE])
    x_amb <- Matrix::rowSums(input[gs != sp, cc, drop = FALSE]); y_amb <- Matrix::rowSums(corrected[gs != sp, cc, drop = FALSE])
    k <- x_end > 0; a <- x_amb > 0
    res <- rbind(res, data.frame(sample = s, species = sp,
      retain = mean(pmax(0, 1 - abs(y_end[k] - x_end[k]) / x_end[k])),
      removal = mean(pmin(1, abs(x_amb[a] - y_amb[a]) / x_amb[a]))))
  }
  by_sample <- stats::aggregate(cbind(retain, removal) ~ sample, res, mean)
  c(retain = mean(by_sample$retain), removal = mean(by_sample$removal))
}

seurat_clusters <- function(cells, resolution = 0.1) {
  so <- suppressWarnings(Seurat::CreateSeuratObject(cells))
  so <- Seurat::NormalizeData(so, verbose = FALSE)
  so <- Seurat::FindVariableFeatures(so, verbose = FALSE)
  so <- Seurat::ScaleData(so, verbose = FALSE)
  so <- Seurat::RunPCA(so, npcs = 10, verbose = FALSE)
  so <- Seurat::FindNeighbors(so, dims = 1:10, verbose = FALSE)
  as.character(Seurat::FindClusters(so, resolution = resolution, verbose = FALSE)$seurat_clusters)
}

run_tools <- function(cells, empties, clusters, cell_species) {
  tools <- list(Decant = run_decant, SoupX = run_soupx, DecontX = run_decontx)
  do.call(rbind, lapply(names(tools), function(nm) {
    invisible(gc()); message("running ", nm)
    t0 <- Sys.time(); r <- tools[[nm]](cells, empties, clusters)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(r)) return(data.frame(method = nm, retain = NA, removal = NA, seconds = NA, note = "not installed"))
    if (!is.null(r$error)) return(data.frame(method = nm, retain = NA, removal = NA, seconds = secs, note = r$error))
    corr <- methods::as(methods::as(r$corrected, "CsparseMatrix"), "dgCMatrix")
    dimnames(corr) <- dimnames(cells)
    sc <- madsen_scores(cells, corr, cell_species)
    message(sprintf("  %s: retain %.3f, removal %.3f, %.0f s", nm, sc[["retain"]], sc[["removal"]], secs))
    data.frame(method = nm, retain = sc[["retain"]], removal = sc[["removal"]], seconds = secs, note = r$note %||% "")
  }))
}

mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode)) mode <- "selftest"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

if (mode == "selftest") {
  ## same generator as tests/testthat/helper-species.R
  source(file.path("tests", "testthat", "helper-species.R"))
  rows <- list()
  for (rm in c(0.05, 0.1, 0.2)) for (seed in 1:2) {
    sim <- sim_species_mix(n_genes = 600, n_cells = 800, n_empty = 4000, rho_mean = rm, seed = seed)
    sp <- ifelse(sim$species == "GRCh38", "human", "mouse")
    rownames(sim$cells) <- rownames(sim$empties) <- sub("^mm10___", "mm10_", rownames(sim$cells))
    colnames(sim$empties) <- paste0("e", seq_len(ncol(sim$empties)))
    cl <- quick_labels(sim$cells, k = 6, seed = seed)
    r <- run_tools(Matrix::Matrix(sim$cells, sparse = TRUE), Matrix::Matrix(sim$empties, sparse = TRUE), cl, sp)
    rows[[length(rows) + 1]] <- cbind(rho = rm, seed = seed, r)
  }
  df <- do.call(rbind, rows)
  tab <- stats::aggregate(cbind(retain, removal) ~ method + rho, df, mean)
  cat("Synthetic two-species mixture (self-test, not comparable to the paper), mean of 2 seeds:\n")
  print(tab, digits = 3, row.names = FALSE)
  utils::write.csv(tab, file.path(OUT, "madsen_selftest.csv"), row.names = FALSE)
}

if (mode == "low_complexity") {
  DATA <- file.path("validation", "data", "madsen")
  f <- file.path(DATA, "20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5")
  if (!file.exists(f)) {
    dir.create(DATA, recursive = TRUE, showWarnings = FALSE)
    utils::download.file(paste0("https://cf.10xgenomics.com/samples/cell-exp/6.1.0/20k_hgmm_3p_HT_nextgem_Chromium_X/",
                                basename(f)), f, mode = "wb")
  }
  ## emptyDrops on ~3M barcodes is the slow step; cache the preprocessed inputs
  prep_f <- file.path(DATA, "low_complexity_prep.rds")
  if (file.exists(prep_f)) {
    prep <- readRDS(prep_f)
  } else {
    raw <- Seurat::Read10X_h5(f); if (is.list(raw)) raw <- raw[["Gene Expression"]]
    raw <- raw[Matrix::rowSums(raw) > 0, ]
    set.seed(1)
    ed <- DropletUtils::emptyDrops(raw, lower = 500)
    is_cell <- !is.na(ed$FDR) & ed$FDR <= 0.05
    cells <- raw[, is_cell]
    gs <- genome_of(rownames(cells))
    hfrac <- Matrix::colSums(cells[gs == "human", ]) / Matrix::colSums(cells)
    keep <- hfrac >= 0.75 | hfrac <= 0.25
    cells <- cells[, keep]; sp <- ifelse(hfrac[keep] >= 0.75, "human", "mouse")
    ## empties: SoupX's default range, among droplets emptyDrops did not call
    tot <- Matrix::colSums(raw)
    empties <- raw[, !is_cell & tot > 0 & tot <= 100]
    cl <- seurat_clusters(cells, resolution = 0.1)
    prep <- list(cells = cells, empties = empties, sp = sp, cl = cl)
    saveRDS(prep, prep_f)
  }
  cells <- prep$cells; empties <- prep$empties; sp <- prep$sp; cl <- prep$cl
  cat(sprintf("cells %d (human %d, mouse %d), empties %d, clusters %d\n",
              ncol(cells), sum(sp == "human"), sum(sp == "mouse"), ncol(empties), length(unique(cl))))
  r <- run_tools(cells, empties, cl, sp)
  print(r, digits = 3, row.names = FALSE)
  cat("\nPublished (Cargnelli et al. 2026, Table S2, low complexity):\n"); print(PUBLISHED_LOW, row.names = FALSE)
  utils::write.csv(r, file.path(OUT, "madsen_low_complexity.csv"), row.names = FALSE)
}
