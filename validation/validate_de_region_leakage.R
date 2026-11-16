#!/usr/bin/env Rscript
## validation/validate_de_region_leakage.R
## Pre-registered in validation/preregistration/region_leakage.md (read it
## first). Paired brain-region comparison in GSE272082: does a standard
## pseudobulk test call genes foreign to microglia (or oligodendrocytes) DE
## between regions because the regions' soups differ, and does ambient_de stop it?
## Run from the repo root:  Rscript validation/validate_de_region_leakage.R
##
## Needs the GSE272082 GEX h5 files in validation/data/gse272082_pfc (fetched
## if missing), the Darmanis 2015 reference (validation/data/darmanis2015,
## built by validate_de_region_leakage.R's truth step from GSE67835), edgeR,
## limma, and optionally SoupX and DecontX.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages({ library(Matrix); library(rhdf5) })
source(file.path("validation", "comparators.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

DIR <- file.path("validation", "data", "gse272082_pfc")
REF <- file.path("validation", "data", "darmanis2015")
OUT <- file.path("validation", "results"); dir.create(OUT, showWarnings = FALSE)
DONORS <- list(D01 = c(PFC = "NIH01", EC = "NIH02", HIP = "NIH03"), D04 = c(PFC = "NIH04", EC = "NIH05", HIP = "NIH06"),
               D10 = c(PFC = "NIH10", EC = "NIH11", HIP = "NIH12"), D13 = c(PFC = "NIH13", EC = "NIH14", HIP = "NIH15"),
               D16 = c(PFC = "NIH16", EC = "NIH17", HIP = "NIH18"), D28 = c(PFC = "NIH28", EC = "NIH29", HIP = "NIH30"))
TYPES <- c("microglia", "oligodendrocyte")
MK <- list(microglia = c("P2RY12", "CSF1R", "CX3CR1"), oligodendrocyte = c("MBP", "PLP1", "MOBP"),
           astrocyte = c("GFAP", "AQP4", "SLC1A2"), OPC = c("PDGFRA", "VCAN"),
           excitatory = c("SLC17A7", "SATB2"), inhibitory = c("GAD1", "GAD2"))
REF_TYPE <- c(microglia = "microglia", oligodendrocyte = "oligodendrocytes")
MIN_CELLS <- 20

## ---------- truth: genes foreign to each tested type (Darmanis 2015) ----------
truth_sets <- function() {
  f <- file.path(REF, "darmanis_adult.rds")
  if (!file.exists(f)) {
    dir.create(file.path(REF, "raw"), recursive = TRUE, showWarnings = FALSE)
    tarf <- file.path(REF, "RAW.tar")
    if (!file.exists(tarf)) utils::download.file(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67835/suppl/GSE67835_RAW.tar", tarf, mode = "wb")
    utils::untar(tarf, exdir = file.path(REF, "raw"))
    gf <- file.path(REF, "gsm.txt")
    if (!file.exists(gf)) utils::download.file(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835&targ=gsm&form=text&view=brief", gf, quiet = TRUE)
    g <- readLines(gf)
    ids <- sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", g, value = TRUE))
    ct <- sub(".*cell type: ", "", grep("cell type:", g, value = TRUE))
    fs <- list.files(file.path(REF, "raw"), full.names = TRUE)
    read1 <- function(x) { d <- utils::read.delim(x, header = FALSE, sep = "\t", strip.white = TRUE); stats::setNames(d[[2]], trimws(d[[1]])) }
    M <- sapply(fs, read1); colnames(M) <- sub("_.*", "", basename(fs))
    M <- M[!grepl("^__|no_feature|ambiguous|alignment_not_unique|too_low_aQual|not_aligned", rownames(M)), ]
    type <- stats::setNames(ct, ids)[colnames(M)]
    adult <- type %in% c("neurons", "astrocytes", "oligodendrocytes", "OPC", "microglia", "endothelial")
    saveRDS(list(counts = M[, adult], type = type[adult]), f)
  }
  d <- readRDS(f); M <- d$counts; type <- d$type
  cpm <- sapply(sort(unique(type)), function(t) { y <- rowSums(M[, type == t, drop = FALSE]); y / sum(y) * 1e6 })
  lapply(REF_TYPE, function(t) {
    mo <- apply(cpm[, setdiff(colnames(cpm), t), drop = FALSE], 1, max)
    rownames(cpm)[mo >= 100 & cpm[, t] <= mo / 100]
  })
}

## ---------- per library ----------
geo_file <- function(sample, kind) {
  fl <- file.path(DIR, "filelist_full.txt")
  if (!file.exists(fl)) utils::download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272082/suppl/filelist.txt", fl, quiet = TRUE)
  nm <- utils::read.delim(fl, stringsAsFactors = FALSE)$Name
  f <- grep(sprintf("_%s_%s_feature_bc_matrix.h5$", sample, kind), nm, value = TRUE)
  dest <- file.path(DIR, f)
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
marker_types <- function(X) {
  lib <- Matrix::colSums(X)
  cpm <- t(t(as.matrix(X[unlist(MK), , drop = FALSE])) / lib) * 1e4
  sc <- sapply(MK, function(g) colMeans(log1p(cpm[g, , drop = FALSE])))
  srt <- t(apply(sc, 1, sort, decreasing = TRUE))
  ifelse(srt[, 1] - srt[, 2] > 1, colnames(sc)[apply(sc, 1, which.max)], "unassigned")
}
process_lib <- function(s) {
  cache <- file.path(DIR, paste0("region_pb_", s, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  raw <- read_h5_gex(geo_file(s, "raw"))
  bc <- as.vector(h5read(geo_file(s, "filtered"), "matrix/barcodes")); H5close()
  d <- split_droplets(raw, cells = bc); rm(raw)
  X <- d$cells[, Matrix::colSums(d$cells > 0) >= 300, drop = FALSE]
  E1 <- Matrix::Matrix(Matrix::rowSums(d$empties), ncol = 1, sparse = TRUE, dimnames = list(rownames(X), NULL))
  rm(d); invisible(gc())
  ty <- marker_types(X)
  pbf <- function(M) sapply(TYPES, function(t) Matrix::rowSums(M[, ty == t, drop = FALSE]))
  out <- list(raw = pbf(X), n = table(factor(ty, levels = c(names(MK), "unassigned"))), soup = as.numeric(E1[, 1]))
  names(out$soup) <- rownames(X)
  dc <- suppressWarnings(decant(X, E1, clusters = ty)); out$decant <- pbf(dc$corrected); rm(dc)
  sx <- run_soupx(X, E1, ty); out$soupx <- if (!is.null(sx$corrected)) pbf(sx$corrected)
  out$soupx_note <- sx$note %||% sx$error %||% "not installed"
  dx <- run_decontx(X, E1, ty); out$decontx <- if (!is.null(dx$corrected)) pbf(dx$corrected)
  out$decontx_note <- dx$note %||% dx$error %||% "not installed"
  cat(sprintf("%s: %d cells (microglia %d, oligodendrocyte %d); soupx %s; decontx %s\n", s, ncol(X),
              out$n["microglia"], out$n["oligodendrocyte"], if (is.null(out$soupx)) out$soupx_note else "ok",
              if (is.null(out$decontx)) out$decontx_note else "ok"))
  saveRDS(out, cache); invisible(gc()); out
}

edger_ql <- function(Y, design) {
  Y <- round(Y); keep <- edgeR::filterByExpr(Y, design = design)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(Y[keep, , drop = FALSE]))
  dge <- edgeR::estimateDisp(dge, design)
  t <- edgeR::glmQLFTest(edgeR::glmQLFit(dge, design), coef = ncol(design))$table
  stats::setNames(stats::p.adjust(t$PValue, "BH"), rownames(t))
}

truth <- truth_sets()
cat("foreign genes:", paste(sprintf("%s %d", names(truth), lengths(truth)), collapse = ", "), "\n")
P <- lapply(stats::setNames(unlist(DONORS), unlist(DONORS)), process_lib)

contrasts <- list(HIP_vs_PFC = "HIP", EC_vs_PFC = "EC")
rows <- list(); calls <- list(); shuf <- list()
for (cn in names(contrasts)) for (t in TYPES) {
  reg <- contrasts[[cn]]
  ok <- vapply(DONORS, function(d) P[[d[["PFC"]]]]$n[t] >= MIN_CELLS && P[[d[[reg]]]]$n[t] >= MIN_CELLS, TRUE)
  dn <- names(DONORS)[ok]
  if (length(dn) < 3) { cat(cn, t, ": fewer than 3 donors qualify\n"); next }
  libs <- unname(c(sapply(DONORS[dn], `[[`, "PFC"), sapply(DONORS[dn], `[[`, reg)))
  donor <- factor(rep(dn, 2)); region <- rep(c(0, 1), each = length(dn))
  design <- cbind(stats::model.matrix(~ donor), region = region)
  mat <- function(m) { x <- lapply(libs, function(l) P[[l]][[m]]); if (any(vapply(x, is.null, TRUE))) NULL else sapply(x, function(z) z[, t]) }
  Y <- mat("raw"); colnames(Y) <- libs
  E <- lapply(libs, function(l) Matrix::Matrix(P[[l]]$soup, ncol = 1, sparse = TRUE, dimnames = list(names(P[[l]]$soup), NULL)))
  fdr <- list(naive_edgeR = edger_ql(Y, design))
  r <- suppressWarnings(ambient_de(Y, empties = E, design = design, naive = FALSE))
  fdr$ambient_de <- stats::setNames(r$FDR, r$gene); dominated <- attr(r, "ambient_dominated")
  for (m in c("decant", "soupx", "decontx")) { Ym <- mat(m); if (!is.null(Ym)) fdr[[paste0(m, "_edgeR")]] <- edger_ql(Ym, design) }
  fg <- intersect(truth[[t]], rownames(Y))
  for (m in names(fdr)) {
    called <- names(fdr[[m]])[fdr[[m]] < 0.05]
    rows[[length(rows) + 1]] <- data.frame(contrast = cn, cell_type = t, method = m, donors = length(dn),
      foreign_genes = length(fg), foreign_called = sum(called %in% fg), other_called = sum(!called %in% fg),
      foreign_skipped_as_ambient = if (m == "ambient_de") sum(fg %in% dominated) else NA)
    if (length(called)) calls[[length(calls) + 1]] <- data.frame(contrast = cn, cell_type = t, method = m,
      gene = called, foreign = called %in% fg, FDR = fdr[[m]][called])
  }
  ## calibration: region labels swapped within donors (donor 1 fixed; the
  ## mirror of a relabeling is the same test)
  nd <- length(dn)
  for (k in seq_len(2^(nd - 1) - 1)) {
    sw <- c(FALSE, as.logical(bitwAnd(k, 2^(0:(nd - 2))) > 0))
    rj <- ifelse(rep(sw, 2), 1 - region, region)
    dj <- cbind(stats::model.matrix(~ donor), region = rj)
    a <- tryCatch(suppressWarnings(ambient_de(Y, empties = E, design = dj, naive = FALSE)),
                  error = function(e) conditionMessage(e))
    shuf[[length(shuf) + 1]] <- data.frame(contrast = cn, cell_type = t, relabeling = k, donors_swapped = sum(sw),
      ambient_de_hits = if (is.character(a)) NA else sum(a$FDR < 0.05), naive_hits = sum(edger_ql(Y, dj) < 0.05),
      error = if (is.character(a)) a else "")
  }
}
res <- do.call(rbind, rows); C <- do.call(rbind, calls); S <- do.call(rbind, shuf)
utils::write.csv(res, file.path(OUT, "region_leakage_summary.csv"), row.names = FALSE)
utils::write.csv(C, file.path(OUT, "region_leakage_calls.csv"), row.names = FALSE)
utils::write.csv(S, file.path(OUT, "region_leakage_relabelings.csv"), row.names = FALSE)
options(width = 200); print(res)
cat("\nRelabelings with a hit:\n")
print(stats::aggregate(cbind(ambient_de = ambient_de_hits > 0, naive = naive_hits > 0, n = 1) ~ contrast + cell_type, S, sum, na.action = na.pass))
if (any(S$error != "")) { cat("\nambient_de errors in relabelings:\n"); print(S[S$error != "", c("contrast", "cell_type", "relabeling", "error")]) }
