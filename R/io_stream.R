## io_stream.R
## Streaming input for large raw droplet matrices. A raw 10x matrix holds
## millions of barcodes (hgmm 10k: 6.8M barcodes, 85M non-zeros; reading it
## whole with read_10x_counts() peaks at 4.3 GB). Decant only needs, from the
## empty droplets, O(genes) statistics: per-gene UMI sums (the soup profile),
## and for dispersion the per-gene sums of squares and detection counts. So
## the empties never need to exist as a matrix.

#' Summary statistics of empty droplets (a streamed soup).
#'
#' An O(genes) summary that [decant()], [ambient_global()] and [pool_soup()]
#' accept in place of an empty-droplet matrix. Build it from a matrix with
#' [soup_summary()] or, without ever loading the empties, with
#' [read_10x_streaming()].
#'
#' @param empty genes x droplets matrix (dense or sparse).
#' @return an object of class `decant_soup`: a list with `gene_sums`,
#'   `gene_sumsq` (sum of squared counts), `gene_detected` (droplets with a
#'   count), `total_umis`, `n_droplets` and `genes` (row names).
#' @export
soup_summary <- function(empty) {
  if (inherits(empty, "decant_soup")) return(empty)
  empty <- .as_counts(empty, "empty")
  sq <- empty; nz <- empty
  if (.is_sparse(empty)) { sq@x <- empty@x^2; nz@x <- as.numeric(empty@x != 0) }
  else { sq <- empty^2; nz <- (empty != 0) * 1 }
  .new_soup(.row_sums(empty), .row_sums(sq), .row_sums(nz), ncol(empty), rownames(empty))
}

.new_soup <- function(sums, sumsq, detected, n_droplets, genes = NULL) {
  structure(list(gene_sums = as.numeric(sums), gene_sumsq = as.numeric(sumsq),
                 gene_detected = as.numeric(detected), total_umis = sum(sums),
                 n_droplets = n_droplets, genes = genes), class = "decant_soup")
}

#' Print a soup summary.
#' @param x a `decant_soup`.
#' @param ... unused.
#' @method print decant_soup
#' @export
print.decant_soup <- function(x, ...) {
  cat(sprintf("<decant_soup> %d genes, %d droplets, %s UMIs\n", length(x$gene_sums), x$n_droplets,
              format(x$total_umis, big.mark = ",")))
  invisible(x)
}

## empty-droplet input: a matrix, or a decant_soup summary
.is_soup <- function(x) inherits(x, "decant_soup")
.soup_nrow <- function(x) if (.is_soup(x)) length(x$gene_sums) else nrow(x)
.soup_gene_sums <- function(x) if (.is_soup(x)) x$gene_sums else .row_sums(.as_counts(x, "empties"))
.soup_umis <- function(x) if (.is_soup(x)) x$total_umis else sum(.col_sums(x))

#' Per-gene negative-binomial dispersion of the empty droplets, from moments.
#'
#' For each gene, the method-of-moments estimate
#' `alpha = max(0, (var - mean) / mean^2)` of per-droplet counts, then shrunk
#' toward a trend in log mean (a running median over genes ordered by mean,
#' weighted by detection), so that rarely detected genes borrow strength.
#' This ignores droplet-to-droplet differences in library size, which inflate
#' `alpha`; in empty droplets (a narrow UMI range) the effect is small.
#'
#' @param soup a [soup_summary()] (or an empty-droplet matrix).
#' @param min_detected genes detected in fewer droplets take the trend value.
#' @return numeric vector of per-gene dispersions (NA where the gene has no
#'   counts), with attribute `"raw"` (unshrunk estimates).
#' @export
soup_dispersion <- function(soup, min_detected = 20) {
  soup <- soup_summary(soup)
  n <- soup$n_droplets
  m <- soup$gene_sums / n
  v <- (soup$gene_sumsq - n * m^2) / max(n - 1, 1)
  raw <- ifelse(m > 0, pmax(0, (v - m) / m^2), NA_real_)
  ok <- is.finite(raw) & soup$gene_detected >= min_detected
  out <- raw
  if (sum(ok) >= 10) {
    o <- which(ok)[order(m[ok])]
    k <- max(3L, 2L * (min(length(o), 201L) %/% 2L) - 1L)       # odd window
    tr <- stats::runmed(raw[o], k, endrule = "constant")
    ## trend at every gene with counts, by interpolation in log mean
    f <- stats::approxfun(log(m[o]), tr, rule = 2, ties = mean)
    trend <- ifelse(m > 0, f(log(pmax(m, 1e-300))), NA_real_)
    w <- soup$gene_detected / (soup$gene_detected + min_detected)
    out <- ifelse(ok, w * raw + (1 - w) * trend, trend)
  }
  attr(out, "raw") <- raw
  out
}

#' Read a raw 10x matrix in a streaming pass, keeping only the cells.
#'
#' Reads `matrix.mtx(.gz)` in chunks, twice: first to total each barcode,
#' then to keep the cells' entries and reduce the empty droplets to a
#' [soup_summary()]. Neither the raw matrix nor the empty-droplet matrix is
#' ever held in memory; peak memory is the cells plus one chunk. Cells and
#' empties are defined as in [split_droplets()].
#'
#' @param dir a Cell Ranger `raw_feature_bc_matrix` directory (matrix.mtx,
#'   features.tsv or genes.tsv, barcodes.tsv; optionally gzipped).
#' @param cells optional character vector of cell barcodes (recommended).
#' @param cell_min_umi used only when `cells` is NULL.
#' @param empty_umi_range inclusive UMI range for empty droplets.
#' @param gene_column,feature_type as in [read_10x_counts()].
#' @param chunk_lines matrix entries read per chunk.
#' @return list(cells = genes x cells dgCMatrix, soup = a `decant_soup`,
#'   summary = the same counts as [split_droplets()]).
#' @export
read_10x_streaming <- function(dir, cells = NULL, cell_min_umi = 500, empty_umi_range = c(1, 100),
                               gene_column = 1, feature_type = "Gene Expression",
                               chunk_lines = 2e6) {
  pick <- function(stems) {
    for (s in stems) for (ext in c("", ".gz")) {
      f <- file.path(dir, paste0(s, ext)); if (file.exists(f)) return(f)
    }
    stop("none of ", paste(stems, collapse = "/"), " found in ", dir, call. = FALSE)
  }
  mtx <- pick("matrix.mtx")
  feats <- utils::read.delim(pick(c("features.tsv", "genes.tsv")), header = FALSE,
                             stringsAsFactors = FALSE)
  ## barcodes: raw matrices have millions; keep only a count, the cells'
  ## names and (without `cells`) nothing else, reading them in chunks
  bc_file <- pick("barcodes.tsv")
  read_barcodes <- function(keep) {          # keep(names) -> logical; returns list(n, keep, names)
    con <- if (grepl("\\.gz$", bc_file)) gzfile(bc_file, "r") else file(bc_file, "r")
    on.exit(close(con))
    n <- 0; k <- list(); nm <- list()
    repeat {
      b <- readLines(con, n = 1e6)
      if (!length(b)) break
      b <- sub("\t.*$", "", b)
      kk <- keep(b); k[[length(k) + 1]] <- kk; nm[[length(nm) + 1]] <- b[kk]; n <- n + length(b)
    }
    list(n = n, keep = unlist(k, use.names = FALSE), names = unlist(nm, use.names = FALSE))
  }
  bc <- read_barcodes(if (is.null(cells)) function(b) rep(FALSE, length(b)) else function(b) b %in% cells)
  n_bc <- bc$n
  gene_keep <- if (!is.null(feature_type) && ncol(feats) >= 3) feats[, 3] == feature_type else rep(TRUE, nrow(feats))
  ## Each pass opens the file, skips the header, then scans triplets chunk by
  ## chunk. The loops run in this function's own frame, so the large vectors
  ## filled below are modified in place (assigning into them from a closure
  ## or an environment copies the whole vector on every chunk).
  open_mtx <- function() {
    con <- if (grepl("\\.gz$", mtx)) gzfile(mtx, "r") else file(mtx, "r")
    repeat { l <- readLines(con, n = 1); if (!length(l) || !startsWith(l, "%")) break }
    dims <- as.numeric(strsplit(trimws(l), "[[:space:]]+")[[1]])
    if (dims[1] != nrow(feats) || dims[2] != n_bc) {
      close(con); stop("matrix dimensions do not match features/barcodes", call. = FALSE)
    }
    con
  }
  next_chunk <- function(con) {
    tr <- scan(con, what = list(0L, 0L, 0), nmax = chunk_lines, quiet = TRUE)
    if (!length(tr[[1]])) return(NULL)
    keep <- gene_keep[tr[[1]]]
    if (all(keep)) tr else lapply(tr, `[`, keep)
  }
  ## pass 1: barcode totals and non-zero counts (selected features only)
  tot <- nnzb <- numeric(n_bc)
  con <- open_mtx()
  while (!is.null(tr <- next_chunk(con))) {
    s <- rowsum(cbind(tr[[3]], 1), tr[[2]], reorder = FALSE); jj <- as.integer(rownames(s))
    tot[jj] <- tot[jj] + s[, 1]; nnzb[jj] <- nnzb[jj] + s[, 2]
  }
  close(con)
  is_cell <- if (!is.null(cells)) {
    miss <- setdiff(cells, bc$names)
    if (length(miss)) stop(length(miss), " cell barcodes not found in `raw`", call. = FALSE)
    bc$keep
  } else tot >= cell_min_umi
  cell_names <- if (!is.null(cells)) bc$names else {
    at <- 0
    read_barcodes(function(b) { k <- is_cell[at + seq_along(b)]; at <<- at + length(b); k })$names
  }
  rm(bc); invisible(gc())
  is_empty <- !is_cell & tot >= empty_umi_range[1] & tot <= empty_umi_range[2]
  if (sum(is_empty) < 100)
    warning("only ", sum(is_empty), " empty droplets: the soup estimate will be noisy. ",
            "Use a raw (unfiltered) matrix.", call. = FALSE)
  ## pass 2: the cells' entries go straight into compressed-column slots (the
  ## column pointers are known from pass 1), and the empties are reduced to
  ## O(G) sums. Nothing the size of the raw matrix, or of the empties, exists.
  gkeep_idx <- cumsum(gene_keep); gkeep_idx[!gene_keep] <- 0L
  G <- sum(gene_keep); nc <- sum(is_cell)
  col_of <- cumsum(is_cell); col_of[!is_cell] <- 0L
  pp <- c(0, cumsum(nnzb[is_cell]))
  ci <- integer(pp[nc + 1]); cx <- numeric(pp[nc + 1]); fill <- numeric(nc)
  unsorted <- integer(0)
  s1 <- s2 <- dt <- numeric(G)
  con <- open_mtx()
  while (!is.null(tr <- next_chunk(con))) {
    i <- gkeep_idx[tr[[1]]]; j <- tr[[2]]; v <- tr[[3]]
    c_ <- which(is_cell[j])
    if (length(c_)) {
      jc <- col_of[j[c_]]; o <- order(jc, method = "radix"); jc <- jc[o]; c_ <- c_[o]
      r <- rle(jc)
      pos <- pp[jc] + fill[jc] + sequence(r$lengths)
      ci[pos] <- i[c_] - 1L; cx[pos] <- v[c_]
      fill[r$values] <- fill[r$values] + r$lengths
      ## rows must increase within a column: note columns where they do not
      later <- pos > pp[jc] + 1
      if (any(later)) {
        bad <- ci[pos[later]] <= ci[pos[later] - 1]
        if (any(bad)) unsorted <- unique(c(unsorted, jc[later][bad]))
      }
    }
    e_ <- is_empty[j]
    if (any(e_)) {
      r <- rowsum(cbind(v[e_], v[e_]^2, 1), i[e_], reorder = FALSE); gg <- as.integer(rownames(r))
      s1[gg] <- s1[gg] + r[, 1]; s2[gg] <- s2[gg] + r[, 2]; dt[gg] <- dt[gg] + r[, 3]
    }
  }
  close(con); rm(tr, i, j, v)
  for (k in unsorted) {
    r <- (pp[k] + 1):pp[k + 1]; o <- order(ci[r]); ci[r] <- ci[r][o]; cx[r] <- cx[r][o]
  }
  gn <- make.unique(feats[, gene_column])[gene_keep]
  cellm <- new("dgCMatrix", i = ci, p = as.integer(pp), x = cx, Dim = c(as.integer(G), as.integer(nc)),
               Dimnames = list(gn, cell_names))
  rm(ci, cx)
  soup <- .new_soup(s1, s2, dt, sum(is_empty), gn)
  list(cells = cellm, soup = soup,
       summary = c(n_cells = sum(is_cell), n_empties = sum(is_empty),
                   n_dropped = sum(!is_cell & !is_empty), soup_umis = sum(s1)))
}
