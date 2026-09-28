## io.R
## Getting real data in. Ambient correction needs the RAW (unfiltered) droplet
## matrix, because the empty droplets ARE the soup measurement. Cell Ranger
## writes it as raw_feature_bc_matrix/{matrix.mtx,features.tsv,barcodes.tsv}[.gz].

#' Read a 10x Genomics matrix directory into a sparse matrix.
#'
#' @param dir path to a Cell Ranger `raw_feature_bc_matrix` (or filtered)
#'   directory containing matrix.mtx, features.tsv (or genes.tsv) and
#'   barcodes.tsv, optionally gzipped.
#' @param gene_column which column of features.tsv to use as row names
#'   (1 = Ensembl ID, 2 = symbol). IDs are the safe default: symbols repeat.
#' @param feature_type keep only features of this type (column 3 of
#'   features.tsv), e.g. "Gene Expression". NULL keeps all.
#' @return genes x barcodes dgCMatrix.
#' @export
read_10x_counts <- function(dir, gene_column = 1, feature_type = "Gene Expression") {
  pick <- function(stems) {
    for (s in stems) for (ext in c("", ".gz")) {
      f <- file.path(dir, paste0(s, ext)); if (file.exists(f)) return(f)
    }
    stop("none of ", paste(stems, collapse = "/"), " found in ", dir, call. = FALSE)
  }
  mtx <- Matrix::readMM(pick("matrix.mtx"))
  feats <- utils::read.delim(pick(c("features.tsv", "genes.tsv")), header = FALSE,
                             stringsAsFactors = FALSE)
  bcs <- utils::read.delim(pick("barcodes.tsv"), header = FALSE,
                           stringsAsFactors = FALSE)[, 1]
  if (nrow(mtx) != nrow(feats) || ncol(mtx) != length(bcs))
    stop("matrix dimensions do not match features/barcodes", call. = FALSE)
  mtx <- .as_counts(mtx)
  dimnames(mtx) <- list(make.unique(feats[, gene_column]), bcs)
  if (!is.null(feature_type) && ncol(feats) >= 3) mtx <- mtx[feats[, 3] == feature_type, , drop = FALSE]
  mtx
}

#' Split a raw droplet matrix into cells and empty droplets.
#'
#' Cells are either the barcodes you pass (e.g. from Cell Ranger's filtered
#' matrix, EmptyDrops or CellBender calls -- recommended) or, failing that,
#' droplets with at least `cell_min_umi` counts. Empties are droplets whose
#' total lies in `empty_umi_range` (SoupX's default soup range is 0-100 UMIs)
#' and that are not cells. Droplets in between are ambiguous and are dropped
#' from both sets rather than risk calling a cell "soup".
#'
#' @param raw genes x barcodes raw matrix (see [read_10x_counts()]).
#' @param cells optional character vector of cell barcodes.
#' @param cell_min_umi UMI threshold used only when `cells` is NULL.
#' @param empty_umi_range inclusive UMI range for empty droplets. The lower
#'   bound defaults to 1 because zero-count barcodes carry no information.
#' @return list(cells = genes x cells, empties = genes x droplets, summary).
#' @export
split_droplets <- function(raw, cells = NULL, cell_min_umi = 500,
                           empty_umi_range = c(1, 100)) {
  raw <- .as_counts(raw, "raw")
  tot <- .col_sums(raw)
  is_cell <- if (!is.null(cells)) {
    if (is.null(colnames(raw))) stop("`raw` has no barcodes (colnames)", call. = FALSE)
    miss <- setdiff(cells, colnames(raw))
    if (length(miss)) stop(length(miss), " cell barcodes not found in `raw`", call. = FALSE)
    colnames(raw) %in% cells
  } else tot >= cell_min_umi
  is_empty <- !is_cell & tot >= empty_umi_range[1] & tot <= empty_umi_range[2]
  if (sum(is_empty) < 100)
    warning("only ", sum(is_empty), " empty droplets: the soup estimate will be noisy. ",
            "Use a raw (unfiltered) matrix.", call. = FALSE)
  list(cells = raw[, is_cell, drop = FALSE], empties = raw[, is_empty, drop = FALSE],
       summary = c(n_cells = sum(is_cell), n_empties = sum(is_empty),
                   n_dropped = sum(!is_cell & !is_empty),
                   soup_umis = sum(tot[is_empty])))
}
