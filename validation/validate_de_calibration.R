#!/usr/bin/env Rscript
## validation/validate_de_calibration.R
## Calibration of ambient_de() on REAL case/control single-nucleus data.
## Run from the repo root:
##   Rscript validation/validate_de_calibration.R [dataset] [config ...]
## Configs: default, tagwise (dispersion before the fix), no_joint_rho,
## no_robust, neither (default: default tagwise).
##
## Disease labels are shuffled between donors, so no gene can be truly DE
## between the shuffled groups, apart from what the remaining true donors in
## each arm carry. The same 24 shuffles (seed 7, at most 2 true disease donors
## in the disease arm) are used for every method. A calibrated test reports a
## hit in few shuffles and about 1% of genes at p < 0.01.
## Naive edgeR on the raw pseudobulk is run alongside as the reference.
##
## Cell types come from a fixed recipe (Seurat LogNormalize, 2000 HVG, 30 PCs,
## Louvain resolution 0.5, clusters labelled by marker means). Seurat is
## needed for this script only.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages({ library(Matrix); library(rhdf5) })

DATASETS <- list(
  ## Sporadic early-onset AD, prefrontal cortex, 10x Multiome (GEX used).
  ## UT04 and UT09 have "sEOAD" in their titles but "Control" in the
  ## disease-state field; the field matches the stated 5 vs 4 design.
  gse272082_pfc = list(gse = "GSE272082",
                       ctrl = c("NIH01", "NIH04", "NIH28", "UT04", "UT09"),
                       case = c("NIH10", "NIH13", "NIH16", "UT2105"),
                       types = c("Microglia", "Oligodendrocyte")),
  ## Pure null: nine control donors of GSE190015 (occipital cortex, one
  ## batch; Control3 to Control11 = TWCC-JC1 to JC9). No disease in either
  ## arm, so every hit is a false positive. Only raw matrices are deposited;
  ## cells are barcodes with >= 500 UMIs.
  gse190015_ctrl = list(gse = "GSE190015",
                        ctrl = paste0("TWCC-JC", 1:9, "-JC", 1:9, "-lib1"), case = character(0),
                        filtered = FALSE, cell_min_umi = 500,
                        types = c("Microglia", "Oligodendrocyte"))
)
CONFIGS <- list(default = list(), tagwise = list(dispersion = "tagwise"), no_joint_rho = list(joint_rho = FALSE),
                no_robust = list(robust = FALSE), neither = list(robust = FALSE, joint_rho = FALSE))
MARKERS <- list(Microglia = c("P2RY12", "CSF1R", "CX3CR1"), Oligodendrocyte = c("MBP", "PLP1", "MOBP"),
                Astrocyte = c("GFAP", "AQP4", "SLC1A2"), OPC = c("PDGFRA", "VCAN"),
                Excitatory = c("SLC17A7", "SATB2"), Inhibitory = c("GAD1", "GAD2"))

args <- commandArgs(trailingOnly = TRUE)
ds <- if (length(args)) args[1] else "gse272082_pfc"
cfgs <- if (length(args) > 1) args[-1] else c("default", "tagwise")
D <- DATASETS[[ds]]; if (is.null(D)) stop("unknown dataset ", ds)
CACHE <- file.path("validation", "data", ds); dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path("validation", "results"); dir.create(OUT, showWarnings = FALSE)
samples <- c(D$ctrl, D$case)

## ---- download (GEO per-sample supplementary files) ----
geo_file <- function(sample, kind) {
  fl <- file.path(CACHE, "filelist.txt")
  if (!file.exists(fl))
    utils::download.file(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/filelist.txt",
                                 substr(D$gse, 1, nchar(D$gse) - 3), D$gse), fl, quiet = TRUE)
  nm <- utils::read.delim(fl, stringsAsFactors = FALSE)$Name
  f <- grep(sprintf("_%s_%s_feature_bc_matrix.h5$", sample, kind), nm, value = TRUE)
  if (length(f) != 1) stop("no unique ", kind, " file for ", sample)
  dest <- file.path(CACHE, f)
  if (!file.exists(dest)) {
    gsm <- sub("_.*", "", f)
    utils::download.file(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%snnn/%s/suppl/%s",
                                 substr(gsm, 1, nchar(gsm) - 3), gsm, f), dest, mode = "wb", quiet = TRUE)
  }
  dest
}
read_h5_gex <- function(f) {
  g <- function(p) as.vector(h5read(f, p))
  m <- sparseMatrix(i = g("matrix/indices") + 1, p = g("matrix/indptr"), x = as.numeric(g("matrix/data")),
                    dims = g("matrix/shape"), dimnames = list(NULL, g("matrix/barcodes")))
  ty <- g("matrix/features/feature_type"); nm <- g("matrix/features/name"); H5close()
  m <- m[ty == "Gene Expression", ]; rownames(m) <- make.unique(nm[ty == "Gene Expression"]); m
}

## ---- cells, empties, cell types, pseudobulks (cached) ----
pbf <- file.path(CACHE, "pseudobulk.rds")
if (!file.exists(pbf)) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("this script needs Seurat for cell typing")
  cells <- list(); empties <- list()
  for (s in samples) {
    raw <- read_h5_gex(geo_file(s, "raw"))
    if (isFALSE(D$filtered)) {
      d <- split_droplets(raw, cell_min_umi = D$cell_min_umi)
    } else {
      bc <- as.vector(h5read(geo_file(s, "filtered"), "matrix/barcodes")); H5close()
      d <- split_droplets(raw, cells = bc)
    }
    colnames(d$cells) <- paste0(s, "_", colnames(d$cells))
    ## the soup needs only per-gene totals over the empty droplets: keep a
    ## one-column matrix of row sums (soup_eb gives the same result)
    cells[[s]] <- d$cells
    empties[[s]] <- Matrix::Matrix(Matrix::rowSums(d$empties), ncol = 1, sparse = TRUE,
                                   dimnames = list(rownames(d$empties), NULL))
    cat(s, ": cells", ncol(d$cells), " empties", ncol(d$empties), "\n")
    rm(raw, d); invisible(gc(verbose = FALSE))
  }
  ## the full matrix waits on disk while Seurat holds its own copy
  X <- do.call(cbind, cells); rm(cells)
  xf <- tempfile(fileext = ".rds"); saveRDS(X, xf, compress = FALSE)
  ## same filters as CreateSeuratObject(min.cells = 10) then
  ## subset(nFeature_RNA >= 300), without Seurat's intermediate copies
  gk <- Matrix::rowSums(X > 0) >= 10
  ck <- Matrix::colSums(X[gk, ] > 0) >= 300
  so <- Seurat::CreateSeuratObject(X[gk, ck]); rm(X); invisible(gc(verbose = FALSE))
  so <- Seurat::NormalizeData(so, verbose = FALSE)
  so <- Seurat::FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
  so <- Seurat::ScaleData(so, verbose = FALSE); so <- Seurat::RunPCA(so, npcs = 30, verbose = FALSE)
  so <- Seurat::FindNeighbors(so, dims = 1:30, verbose = FALSE)
  so <- Seurat::FindClusters(so, resolution = 0.5, random.seed = 0, verbose = FALSE)
  avg <- suppressMessages(Seurat::AverageExpression(so, features = unlist(MARKERS),
                                                    group.by = "seurat_clusters", verbose = FALSE)$RNA)
  sc <- sapply(MARKERS, function(g) colMeans(log1p(as.matrix(avg[g, , drop = FALSE]))))
  z <- scale(sc)
  lab <- colnames(sc)[apply(z, 1, which.max)]; lab[apply(z, 1, max) < 1] <- "Other"
  names(lab) <- sub("^g", "", rownames(sc))
  ct <- stats::setNames(lab[as.character(so$seurat_clusters)], colnames(so))
  print(table(sub("_.*", "", names(ct)), ct))
  ## pseudobulks use the typed cells (>= 300 genes)
  rm(so, avg); invisible(gc(verbose = FALSE))
  X <- readRDS(xf); unlink(xf)
  pb <- lapply(stats::setNames(D$types, D$types), function(t) sapply(samples, function(s) {
    cs <- names(ct)[ct == t & startsWith(names(ct), paste0(s, "_"))]
    Matrix::rowSums(X[, cs, drop = FALSE])
  }))
  saveRDS(list(pb = pb, empties = empties, celltypes = ct), pbf)
}
P <- readRDS(pbf)

## ---- shuffled-label calibration ----
is_case <- samples %in% D$case
combos <- utils::combn(length(samples), if (length(D$case)) length(D$case) else length(samples) %/% 2)
ov <- apply(combos, 2, function(i) sum(is_case[i]))
set.seed(7); pick <- sample(which(ov <= 2), 24)
rows <- list()
for (t in D$types) {
  Y <- P$pb[[t]]
  E <- lapply(P$empties, function(e) e[rownames(Y), , drop = FALSE])
  soup <- soup_eb(E)   # condition-blind soup, shared by all shuffles
  for (k in cfgs) for (j in c(if (any(is_case)) 0, pick)) {
    grp <- if (j == 0) is_case else seq_along(samples) %in% combos[, j]
    d <- cbind(1, as.numeric(grp))
    r <- suppressWarnings(do.call(ambient_de, c(list(Y, design = d, soup = soup), CONFIGS[[k]])))
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, cell_type = t, config = k,
      labels = if (j == 0) "true" else "shuffled", case_arm = paste(samples[grp], collapse = ","),
      true_cases_in_arm = sum(is_case[grp]), genes_tested = nrow(r),
      hits = sum(r$FDR < 0.05), frac_p01 = mean(r$PValue < 0.01),
      naive_hits = sum(r$naive_FDR < 0.05), naive_frac_p01 = mean(r$naive_PValue < 0.01))
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(OUT, sprintf("de_calibration_%s.csv", ds)), row.names = FALSE)
sh <- res[res$labels == "shuffled", ]
summ <- do.call(rbind, lapply(split(sh, list(sh$cell_type, sh$config), drop = TRUE), function(x) data.frame(
  cell_type = x$cell_type[1], config = x$config[1], shuffles = nrow(x),
  shuffles_with_hits = sum(x$hits > 0), median_hits = stats::median(x$hits), max_hits = max(x$hits),
  median_frac_p01 = round(stats::median(x$frac_p01), 3),
  naive_shuffles_with_hits = sum(x$naive_hits > 0), naive_median_frac_p01 = round(stats::median(x$naive_frac_p01), 3))))
rownames(summ) <- NULL; print(summ)
if (any(res$labels == "true")) {
  cat("\nTrue labels:\n"); print(res[res$labels == "true", c("cell_type", "config", "hits", "naive_hits", "frac_p01")])
}
