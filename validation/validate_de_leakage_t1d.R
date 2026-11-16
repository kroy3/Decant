#!/usr/bin/env Rscript
## validation/validate_de_leakage_t1d.R
## Real-data test of ambient-driven differential expression, pre-registered in
## validation/preregistration/t1d_leakage.md (read it first).
## Run from the repo root:  Rscript validation/validate_de_leakage_t1d.R
##
## GSE273597: whole-pancreas snRNA-seq of non-diabetic (ND) and type 1
## diabetic (T1D) donors. T1D destroys beta cells, so beta-cell RNA in the
## ambient soup differs between the groups. A gene that a cell type does not
## express (by an external, ambient-free reference: Segerstolpe Smart-seq2
## islets) but that a test calls DE in that cell type is an ambient-driven
## false positive.
##
## Needs: the raw h5 files of GSE273597 in validation/data/gse273597 (fetched
## if missing), the Segerstolpe counts in validation/data/segerstolpe
## (benchmark_negative_control.R fetches them), edgeR, limma, and optionally
## SoupX and DecontX (validation/comparators.R).

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages({ library(Matrix); library(rhdf5) })
source(file.path("validation", "comparators.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

DIR <- file.path("validation", "data", "gse273597"); dir.create(DIR, recursive = TRUE, showWarnings = FALSE)
SEG <- file.path("validation", "data", "segerstolpe")
OUT <- file.path("validation", "results"); dir.create(OUT, showWarnings = FALSE)
LAB <- utils::read.csv(file.path("validation", "t1d_donor_labels.csv"), colClasses = "character")
TYPES <- c("acinar", "ductal", "alpha")
MK <- list(beta = "INS", alpha = "GCG", delta = "SST",
           acinar = c("PRSS1", "CPA1", "CTRB2", "CELA3A"), ductal = c("KRT19", "CFTR"))
MIN_CELLS <- 20

## Optional argument: a tag for the output files, so that reruns with later
## code do not overwrite the pre-registered results (t1d_leakage_*.csv).
TAG <- commandArgs(trailingOnly = TRUE)[1]
SFX <- if (is.na(TAG)) "" else paste0("_", TAG)

## ---------- truth: genes foreign to each tested type (Segerstolpe) ----------
truth_sets <- function() {
  f <- file.path(SEG, "pancreas_refseq_rpkms_counts_3514sc.txt")
  sd <- utils::read.delim(file.path(SEG, "E-MTAB-5061.sdrf.txt"), check.names = FALSE, stringsAsFactors = FALSE)
  ok <- sd[["Characteristics [disease]"]] == "normal" & sd[["Characteristics [submitted single cell quality]"]] == "OK"
  lab <- stats::setNames(sub(" cell$", "", sd[["Characteristics [inferred cell type]"]][ok]), sd[["Source Name"]][ok])
  wells <- strsplit(readLines(f, n = 1), "\t")[[1]][-1]; nw <- length(wells)
  use <- which(wells %in% names(lab))
  tab <- data.table::fread(f, skip = 1, header = FALSE, sep = "\t", showProgress = FALSE,
                           select = c(1, 2 + nw + use))
  sym <- tab[[1]]; cnt <- as.matrix(tab[, -1, with = FALSE]); rm(tab)
  cnt <- rowsum(round(cnt), sym)                      # one row per symbol
  ct <- lab[wells[use]]
  big <- names(which(table(ct) >= 10)); big <- setdiff(big, c("co-expression", "unclassified endocrine"))
  cpm <- sapply(big, function(t) { y <- rowSums(cnt[, ct == t, drop = FALSE]); y / sum(y) * 1e6 })
  det <- sapply(big, function(t) rowMeans(cnt[, ct == t, drop = FALSE] > 0))
  cat("Segerstolpe reference types:", paste(sprintf("%s (%d)", big, as.vector(table(ct)[big])), collapse = ", "), "\n")
  lapply(stats::setNames(TYPES, TYPES), function(t) {
    oth <- setdiff(big, t)
    mo <- apply(cpm[, oth, drop = FALSE], 1, max)
    rownames(cpm)[mo >= 100 & cpm[, t] <= mo / 100]   # Amendment 2: no detection rule
  })
}

## ---------- per library: cells, cell types, corrected pseudobulks ----------
read_h5_gex <- function(f) {
  g <- function(p) as.vector(h5read(f, p))
  m <- sparseMatrix(i = g("matrix/indices") + 1, p = g("matrix/indptr"), x = as.numeric(g("matrix/data")),
                    dims = g("matrix/shape"), dimnames = list(NULL, g("matrix/barcodes")))
  ty <- g("matrix/features/feature_type"); nm <- g("matrix/features/name"); H5close()
  m <- m[ty == "Gene Expression", ]; rownames(m) <- make.unique(nm[ty == "Gene Expression"]); m
}
marker_types <- function(X) {
  lib <- Matrix::colSums(X)
  cpm <- t(t(as.matrix(X[unlist(MK), , drop = FALSE])) / lib) * 1e4
  sc <- sapply(MK, function(g) colMeans(log1p(cpm[g, , drop = FALSE])))
  srt <- t(apply(sc, 1, sort, decreasing = TRUE))
  ifelse(srt[, 1] - srt[, 2] > 1, colnames(sc)[apply(sc, 1, which.max)], "unclear")
}
lib_file <- function(lib) {
  fl <- file.path(DIR, "filelist.txt")
  if (!file.exists(fl)) utils::download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273597/suppl/filelist.txt", fl, quiet = TRUE)
  nm <- utils::read.delim(fl, stringsAsFactors = FALSE)$Name
  f <- grep(sprintf("_%s_raw_feature_bc_matrix.h5$", lib), nm, value = TRUE)
  dest <- file.path(DIR, f)
  if (!file.exists(dest)) {
    gsm <- sub("_.*", "", f)
    utils::download.file(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%snnn/%s/suppl/%s",
                                 substr(gsm, 1, nchar(gsm) - 3), gsm, f), dest, mode = "wb", quiet = TRUE)
  }
  dest
}
process_lib <- function(lib) {
  cache <- file.path(DIR, paste0("pb_", lib, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  d <- suppressWarnings(split_droplets(read_h5_gex(lib_file(lib)), cell_min_umi = 500))
  X <- d$cells[, Matrix::colSums(d$cells > 0) >= 300, drop = FALSE]; E <- d$empties; rm(d)
  ty <- marker_types(X)
  pbf <- function(M) sapply(TYPES, function(t) Matrix::rowSums(M[, ty == t, drop = FALSE]))
  out <- list(raw = pbf(X), n = table(factor(ty, levels = c(names(MK), "unclear"))),
              soup = Matrix::rowSums(E), soup_umis = sum(E))
  t0 <- Sys.time()
  dc <- suppressWarnings(decant(X, E, clusters = ty))
  out$decant <- pbf(dc$corrected); rm(dc)
  sx <- run_soupx(X, E, ty)
  out$soupx <- if (!is.null(sx$corrected)) pbf(sx$corrected); out$soupx_note <- sx$note %||% sx$error %||% "not installed"
  dx <- run_decontx(X, E, ty)
  out$decontx <- if (!is.null(dx$corrected)) pbf(dx$corrected); out$decontx_note <- dx$note %||% dx$error %||% "not installed"
  cat(sprintf("%s: %d cells (acinar %d, ductal %d, alpha %d, beta %d)  %.0f s  soupx: %s  decontx: %s\n",
              lib, ncol(X), out$n["acinar"], out$n["ductal"], out$n["alpha"], out$n["beta"],
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              if (is.null(out$soupx)) out$soupx_note else "ok", if (is.null(out$decontx)) out$decontx_note else "ok"))
  saveRDS(out, cache); out
}

## ---------- tests ----------
## The standard pipeline: filterByExpr on the counts, TMM, QL F-test.
edger_ql <- function(Y, design) {
  Y <- round(Y); keep <- edgeR::filterByExpr(Y, design = design)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(Y[keep, , drop = FALSE]))
  dge <- edgeR::estimateDisp(dge, design)
  t <- edgeR::glmQLFTest(edgeR::glmQLFit(dge, design), coef = ncol(design))$table
  stats::setNames(stats::p.adjust(t$PValue, "BH"), rownames(t))
}

truth <- truth_sets()
cat("foreign genes:", paste(sprintf("%s %d", names(truth), lengths(truth)), collapse = ", "), "\n")

libs <- LAB$library
P <- lapply(stats::setNames(libs, libs), process_lib)

## pool libraries by donor
donors <- unique(LAB$donor)
pool <- function(field) lapply(stats::setNames(donors, donors), function(dn) {
  ls <- LAB$library[LAB$donor == dn]; x <- P[[ls[1]]][[field]]
  if (is.null(x)) return(NULL)
  for (l in ls[-1]) { y <- P[[l]][[field]]; if (is.null(y)) return(NULL); x <- x + y }
  x
})
D <- list(raw = pool("raw"), decant = pool("decant"), soupx = pool("soupx"), decontx = pool("decontx"),
          soup = pool("soup"), n = pool("n"))
cond <- stats::setNames(LAB$condition[match(donors, LAB$donor)], donors)
sub <- stats::setNames(LAB$subtype[match(donors, LAB$donor)], donors)

contrasts <- list(primary = list(ctrl = names(cond)[cond == "Control"], case = names(cond)[cond == "T1D"]),
                  late = list(ctrl = names(cond)[cond == "Control"], case = names(sub)[sub == "T1D_late"]))
rows <- list(); genes <- list(); shuf <- list()
for (cn in names(contrasts)) for (t in TYPES) {
  ds <- c(contrasts[[cn]]$ctrl, contrasts[[cn]]$case)
  ds <- ds[vapply(ds, function(dn) D$n[[dn]][t] >= MIN_CELLS, TRUE)]
  is_case <- ds %in% contrasts[[cn]]$case
  if (sum(is_case) < 2 || sum(!is_case) < 2) next
  design <- cbind(1, as.numeric(is_case))
  mat <- function(m) { x <- sapply(ds, function(dn) D[[m]][[dn]][, t]); if (is.list(x)) NULL else x }
  Y <- mat("raw")
  E <- lapply(ds, function(dn) Matrix::Matrix(D$soup[[dn]], ncol = 1, sparse = TRUE, dimnames = list(names(D$soup[[dn]]), NULL)))
  fdr <- list(naive_edgeR = edger_ql(Y, design))
  r <- ambient_de(Y, empties = E, design = design, naive = FALSE)
  fdr$ambient_de <- stats::setNames(r$FDR, r$gene)
  dominated <- attr(r, "ambient_dominated")
  for (m in c("decant", "soupx", "decontx")) { Ym <- mat(m); if (!is.null(Ym)) fdr[[paste0(m, "_edgeR")]] <- edger_ql(Ym, design) }
  fg <- intersect(truth[[t]], rownames(Y))
  for (m in names(fdr)) {
    called <- names(fdr[[m]])[fdr[[m]] < 0.05]
    rows[[length(rows) + 1]] <- data.frame(contrast = cn, cell_type = t, method = m,
      n_ctrl = sum(!is_case), n_case = sum(is_case), foreign_genes = length(fg),
      foreign_called = sum(called %in% fg), other_called = sum(!called %in% fg),
      foreign_tested = sum(names(fdr[[m]]) %in% fg),
      foreign_skipped_as_ambient = if (m == "ambient_de") sum(fg %in% dominated) else NA)
    if (length(called)) genes[[length(genes) + 1]] <- data.frame(contrast = cn, cell_type = t, method = m,
      gene = called, foreign = called %in% fg, FDR = fdr[[m]][called])
  }
  ## calibration (primary contrast): shuffled labels
  if (cn == "primary") {
    set.seed(7); n <- length(ds); k <- sum(is_case)
    for (j in 1:24) {
      g <- seq_len(n) %in% sample(n, k); dj <- cbind(1, as.numeric(g))
      rj <- suppressWarnings(ambient_de(Y, empties = E, design = dj, naive = FALSE))
      nj <- edger_ql(Y, dj)
      shuf[[length(shuf) + 1]] <- data.frame(cell_type = t, shuffle = j, true_cases_in_arm = sum(is_case[g]),
        ambient_de_hits = sum(rj$FDR < 0.05), naive_hits = sum(nj < 0.05))
    }
  }
}
res <- do.call(rbind, rows); G <- do.call(rbind, genes); S <- do.call(rbind, shuf)
utils::write.csv(res, file.path(OUT, paste0("t1d_leakage_summary", SFX, ".csv")), row.names = FALSE)
utils::write.csv(G, file.path(OUT, paste0("t1d_leakage_calls", SFX, ".csv")), row.names = FALSE)
utils::write.csv(S, file.path(OUT, paste0("t1d_leakage_shuffles", SFX, ".csv")), row.names = FALSE)
options(width = 200)
print(res)
cat("\nShuffled labels (primary contrast): shuffles with a hit, of 24\n")
print(stats::aggregate(cbind(ambient_de = ambient_de_hits > 0, naive = naive_hits > 0) ~ cell_type, S, sum))
cat("\nForeign genes called in acinar, primary contrast:\n")
print(G[G$contrast == "primary" & G$cell_type == "acinar" & G$foreign, c("method", "gene", "FDR")])
