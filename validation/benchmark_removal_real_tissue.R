#!/usr/bin/env Rscript
## validation/benchmark_removal_real_tissue.R
## Pre-registered in validation/preregistration/real_tissue_removal.md (read it
## first). Ambient removal in real brain tissue (GSE272082, six PFC libraries),
## scored against genes foreign to each cell type in an independent C1
## reference (Darmanis 2015). Compares decant() bayes_absent (default) and
## bayes_cell, SoupX and DecontX on identical cells, empties and labels.
## Run from the repo root:  Rscript validation/benchmark_removal_real_tissue.R
## Results: validation/results/removal_real_tissue_*.csv

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages({ library(Matrix); library(rhdf5) })
source(file.path("validation", "comparators.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

## data access and the fixed cell-type recipe, taken verbatim from the region test
src <- readLines(file.path("validation", "validate_de_region_leakage.R"))
grab <- function(from, to) src[grep(from, src)[1]:(grep(to, src)[1] - 1)]
eval(parse(text = grab("^DIR <- file.path", "^DONORS <- list")))
eval(parse(text = grab("^MK <- list", "^REF_TYPE <- ")))
eval(parse(text = grab("^geo_file <- function", "^process_lib <- function")))

LIBS <- c("NIH01", "NIH04", "NIH10", "NIH13", "NIH16", "NIH28")
SCORED <- c(microglia = "microglia", oligodendrocyte = "oligodendrocytes", astrocyte = "astrocytes")
METHODS <- c("bayes_absent", "bayes_cell", "SoupX", "DecontX")

## ---------- truth (Darmanis 2015, C1) ----------
ref <- readRDS(file.path(REF, "darmanis_adult.rds"))
cpm <- sapply(sort(unique(ref$type)), function(t) { y <- rowSums(ref$counts[, ref$type == t, drop = FALSE]); y / sum(y) * 1e6 })
sets <- lapply(SCORED, function(t) {
  others <- apply(cpm[, setdiff(colnames(cpm), t), drop = FALSE], 1, max)
  list(foreign = rownames(cpm)[others >= 100 & cpm[, t] <= others / 100],
       native = rownames(cpm)[cpm[, t] >= 100 & cpm[, t] >= 10 * others])
})
cat("truth sets:", paste(sprintf("%s foreign %d native %d", names(sets),
    sapply(sets, function(s) length(s$foreign)), sapply(sets, function(s) length(s$native))), collapse = "; "), "\n")

## ---------- per library ----------
score <- function(X, Y, ty, soup_p, T_c) {
  do.call(rbind, lapply(names(SCORED), function(t) {
    j <- which(ty == t); if (length(j) < 20) return(NULL)
    f <- intersect(sets[[t]]$foreign, rownames(X)); nn <- intersect(sets[[t]]$native, rownames(X))
    Xf <- X[f, j, drop = FALSE]; Yf <- Y[f, j, drop = FALSE]
    ## doublet-level foreign entries: improbable under ambient alone at rho bound 0.15
    lam <- outer(soup_p[f], 0.15 * T_c[j])
    xf <- as.matrix(Xf); hi <- xf > 0 & stats::ppois(xf - 1, lam, lower.tail = FALSE) < 0.001
    data.frame(type = t, cells = length(j),
               foreign_removed = 1 - sum(Yf) / sum(Xf),
               native_removed = 1 - sum(Y[nn, j]) / sum(X[nn, j]),
               doublet_level_kept = if (any(hi)) sum(as.matrix(Yf)[hi]) / sum(xf[hi]) else NA,
               doublet_level_entries = sum(hi))
  }))
}

rows <- list(); fire <- list()
for (s in LIBS) {
  raw <- read_h5_gex(geo_file(s, "raw"))
  bc <- as.vector(h5read(geo_file(s, "filtered"), "matrix/barcodes")); H5close()
  d <- split_droplets(raw, cells = bc); rm(raw)
  X <- d$cells[, Matrix::colSums(d$cells > 0) >= 300, drop = FALSE]
  soup_c <- Matrix::rowSums(d$empties); rm(d); invisible(gc())
  E1 <- Matrix::Matrix(soup_c, ncol = 1, sparse = TRUE, dimnames = list(rownames(X), NULL))
  ty <- marker_types(X); T_c <- Matrix::colSums(X); soup_p <- soup_c / sum(soup_c)
  cat(sprintf("== %s: %d cells; %s\n", s, ncol(X), paste(sprintf("%s %d", names(table(ty)), table(ty)), collapse = ", ")))
  for (m in METHODS) {
    t0 <- Sys.time()
    r <- switch(m,
      bayes_absent = run_decant(X, E1, ty, "bayes_absent"),
      bayes_cell = run_decant(X, E1, ty, "bayes_cell"),
      SoupX = run_soupx(X, E1, ty), DecontX = run_decontx(X, E1, ty))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(r$corrected)) { cat(sprintf("  %s: %s\n", m, r$error %||% r$note %||% "not installed")); next }
    Y <- methods::as(r$corrected, "CsparseMatrix"); dimnames(Y) <- dimnames(X)
    stopifnot(max(Y - X) <= 1e-6 || m %in% c("SoupX", "DecontX"))
    sc <- score(X, Y, ty, soup_p, T_c)
    if (!is.null(sc)) rows[[length(rows) + 1]] <- cbind(library = s, method = m, sc, seconds = secs)
    if (m == "bayes_absent") {
      ## was the gene-level call used in each cluster? (same inputs as decant())
      for (k in sort(unique(ty))) {
        j <- ty == k
        w <- .absence_posterior(Matrix::rowSums(X[, j, drop = FALSE]), soup_p, sum(r$rho[j] * T_c[j]), sum(soup_c))
        f <- attr(w, "fit")
        fire[[length(fire) + 1]] <- data.frame(library = s, cluster = k, cells = sum(j), used = as.logical(f[["used"]]),
          genes_absent = sum(w > 0.5), m_A = f[["m_A"]], s_A = f[["s_A"]], m_B = f[["m_B"]], s_B = f[["s_B"]])
      }
    }
    rm(r, Y); invisible(gc())
  }
  print(do.call(rbind, rows)[do.call(rbind, rows)$library == s, ], digits = 3, row.names = FALSE)
  rm(X, E1); invisible(gc())
}

res <- do.call(rbind, rows); fires <- do.call(rbind, fire)
dir.create(file.path("validation", "results"), showWarnings = FALSE)
utils::write.csv(res, file.path("validation", "results", "removal_real_tissue_by_library.csv"), row.names = FALSE)
utils::write.csv(fires, file.path("validation", "results", "removal_real_tissue_absence_calls.csv"), row.names = FALSE)
agg <- stats::aggregate(cbind(foreign_removed, native_removed, doublet_level_kept) ~ type + method, res, mean)
utils::write.csv(agg, file.path("validation", "results", "removal_real_tissue_summary.csv"), row.names = FALSE)
cat("\n== Mean over libraries\n"); print(agg[order(agg$type, agg$method), ], digits = 3, row.names = FALSE)
cat("\n== bayes_absent gene-level call used, by cluster type\n")
print(stats::aggregate(cbind(used, genes_absent) ~ cluster, fires, mean), digits = 3, row.names = FALSE)

## pre-stated outcome (microglia), as written in the pre-registration
a <- agg[agg$type == "microglia" & agg$method == "bayes_absent", ]; b <- agg[agg$type == "microglia" & agg$method == "bayes_cell", ]
d_f <- a$foreign_removed - b$foreign_removed; d_n <- a$native_removed - b$native_removed
verdict <- if (d_n > 0.01) "HARMFUL" else if (d_f >= 0.05) "GAIN HOLDS" else if (abs(d_f) <= 0.01) "NO GAIN (falls back)" else "IN BETWEEN"
cat(sprintf("\nmicroglia: foreign removed %+.3f, native removed %+.3f (bayes_absent - bayes_cell) => %s\n", d_f, d_n, verdict))
