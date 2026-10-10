#!/usr/bin/env Rscript
## validation/kidney_mixture.R
## Per-cell rho validation against genotype ground truth.
##
## Data: Janssen et al. (2023) Genome Biology 24:140, Zenodo 10.5281/zenodo.15006097.
## Each 10x channel pools kidney cells from CAST/EiJ (M. m. castaneus) with
## BL6 and 129S1 (M. m. domesticus). In CAST cells, reads carrying domesticus
## alleles at ~32k homozygous SNPs must be background, so the per-cell
## background fraction is measured, not modelled. Decant never sees the SNPs.
##
## Caveat: their estimate counts ALL background (ambient + barcode swapping /
## chimeras), while Decant models ambient only. The paper finds ambient is the
## majority, so expect Decant slightly BELOW truth on average, not above.
##
## Run from the repo root:  Rscript validation/kidney_mixture.R [replicate ...]
## Needs network access to zenodo.org, or the zips unpacked under
## validation/data/kidney/<replicate>/ by hand.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))

RECORD <- "15006097"
CACHE  <- file.path("validation", "data", "kidney")
OUT    <- file.path("validation", "results")

## ---- fetch: list the record's zips and unpack each into its own folder ----
fetch_all <- function() {
  dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
  have <- list.dirs(CACHE, recursive = FALSE, full.names = FALSE)
  if (length(have)) return(have)
  meta <- jsonlite::fromJSON(sprintf("https://zenodo.org/api/records/%s", RECORD))
  files <- meta$files
  for (i in seq_len(nrow(files))) {
    key <- files$key[i]; url <- files$links$self[i]
    if (!grepl("\\.zip$", key)) next
    rep <- sub("\\.zip$", "", key)
    zip <- file.path(CACHE, key)
    if (!file.exists(zip)) utils::download.file(url, zip, mode = "wb", quiet = TRUE)
    utils::unzip(zip, exdir = file.path(CACHE, rep))
    unlink(zip)
  }
  list.dirs(CACHE, recursive = FALSE, full.names = FALSE)
}

find1 <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("no file matching '", pattern, "' under ", dir, call. = FALSE)
  f[1]
}

## barcodes differ between Cell Ranger, Seurat (prefixes/suffixes) and the
## truth table; match on the 16-nt 10x barcode core
core <- function(x) sub(".*([ACGT]{16}).*", "\\1", x)

read_truth <- function(path) {
  x <- readRDS(path)
  if (is.numeric(x) && !is.null(names(x))) return(stats::setNames(x, core(names(x))))
  x <- as.data.frame(x)
  bc <- names(x)[vapply(x, function(v) is.character(v) && mean(grepl("[ACGT]{16}", v)) > 0.9, NA)]
  bc <- if (length(bc)) x[[bc[1]]] else rownames(x)
  num <- names(x)[vapply(x, is.numeric, NA)]
  pick <- grep("cont|rho|binom|frac", num, ignore.case = TRUE, value = TRUE)
  pick <- if (length(pick)) pick[1] else num[1]
  if (is.null(pick) || is.na(pick)) stop("no numeric contamination column in ", path, call. = FALSE)
  message("    truth column: ", pick, " (", nrow(x), " cells)")
  stats::setNames(x[[pick]], core(bc))
}

read_clusters <- function(path) {
  obj <- readRDS(path)
  md <- obj@meta.data
  col <- intersect(c("celltype", "cell_type", "CellType", "seurat_clusters"), names(md))
  cl <- if (length(col)) md[[col[1]]] else as.character(SeuratObject::Idents(obj))
  message("    clusters from: ", if (length(col)) col[1] else "Idents", " (",
          length(unique(cl)), " groups)")
  stats::setNames(as.character(cl), core(colnames(obj)))
}

## ---- one replicate ----
run_rep <- function(rep) {
  dir <- file.path(CACHE, rep)
  message("== ", rep)
  raw <- Seurat::Read10X_h5(find1(dir, "raw_feature_bc_matrix\\.h5$"))
  if (is.list(raw)) raw <- raw[["Gene Expression"]]
  colnames(raw) <- core(colnames(raw))
  truth <- read_truth(find1(dir, "perCell.*CAST.*\\.RDS$"))
  cl_all <- read_clusters(find1(dir, "^seurat\\.RDS$"))

  ## Decant sees every called cell (all strains), clustered as in the paper
  cells <- intersect(names(cl_all), colnames(raw))
  sp <- split_droplets(raw, cells = cells)
  res <- suppressWarnings(decant(sp$cells, sp$empties, clusters = cl_all[colnames(sp$cells)]))
  est <- stats::setNames(res$rho, colnames(sp$cells))

  ## score on CAST cells, where truth exists
  common <- intersect(names(truth), names(est))
  t <- truth[common]; e <- est[common]
  ok <- is.finite(t) & is.finite(e); t <- t[ok]; e <- e[ok]
  cl <- cl_all[names(t)]
  per_cl <- stats::aggregate(cbind(truth = t, decant = e) ~ cluster,
                             data.frame(t, e, cluster = cl), stats::median)
  data.frame(replicate = rep, n_cells = length(t),
             truth_median = stats::median(t), decant_median = stats::median(e),
             ratio = stats::median(e) / stats::median(t),
             mae = mean(abs(e - t)),
             spearman_cell = suppressWarnings(stats::cor(e, t, method = "spearman")),
             spearman_cluster = suppressWarnings(stats::cor(per_cl$decant, per_cl$truth,
                                                            method = "spearman")),
             soup_mismatch = isTRUE(res$soup_mismatch))
}

main <- function(reps) {
  avail <- fetch_all()
  if (!length(reps)) reps <- avail
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  tab <- do.call(rbind, lapply(reps, run_rep))
  print(tab, digits = 3, row.names = FALSE)
  utils::write.csv(tab, file.path(OUT, "kidney_mixture.csv"), row.names = FALSE)
  invisible(tab)
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
