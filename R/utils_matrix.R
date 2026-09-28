## utils_matrix.R
## Matrix plumbing shared by every module. Real droplet data is ~30k genes x
## 10^4-10^5 cells and >90% zeros, so nothing on the main path may build a dense
## genes x cells intermediate (the v0 correction did, via outer(), which is
## ~2.4 GB for 30k x 10k). Everything here works on base matrices AND on
## Matrix::dgCMatrix, and corrections touch only the non-zero entries -- which
## is sufficient, because a zero count can never have anything removed from it.

.is_sparse <- function(x) inherits(x, "sparseMatrix")

## Coerce any supported count container to either a base numeric matrix or a
## dgCMatrix (never silently densify a sparse input).
.as_counts <- function(x, what = "counts") {
  if (.is_sparse(x)) {
    x <- methods::as(x, "CsparseMatrix")
    x <- methods::as(x, "generalMatrix")
    return(methods::as(x, "dMatrix"))
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x))
    stop(what, " must be a numeric matrix or a Matrix::sparseMatrix (genes x barcodes)",
         call. = FALSE)
  if (anyNA(x)) stop(what, " contains NA values", call. = FALSE)
  x
}

.col_sums <- function(x) as.numeric(Matrix::colSums(x))
.row_sums <- function(x) as.numeric(Matrix::rowSums(x))

## Non-zero entries as (row, col, value) triplets.
.nz <- function(x) {
  if (.is_sparse(x)) {
    list(i = x@i + 1L, j = rep.int(seq_len(ncol(x)), diff(x@p)), v = x@x)
  } else {
    idx <- which(x != 0)
    G <- nrow(x)
    list(i = (idx - 1L) %% G + 1L, j = (idx - 1L) %/% G + 1L, v = x[idx], idx = idx)
  }
}

## Write new values back into the non-zero slots returned by .nz(). With
## drop = FALSE a sparse result keeps the input's exact pattern (explicit
## zeros), which lets .check_mass() compare values slot-for-slot.
.nz_set <- function(x, nz, v, drop = TRUE) {
  if (.is_sparse(x)) {
    x@x <- as.numeric(v)
    if (drop) Matrix::drop0(x) else x
  } else {
    storage.mode(x) <- "double"
    x[nz$idx] <- v
    x
  }
}

## The guarantee, checked element-wise: 0 <= corrected <= observed. When both
## are sparse with the same pattern this is a cheap slot comparison; otherwise
## fall back to a (memory-hungry) sparse subtraction.
.check_mass <- function(observed, corrected, tol = 1e-6) {
  ok <- if (.is_sparse(observed) && .is_sparse(corrected) &&
            identical(observed@p, corrected@p) && identical(observed@i, corrected@i)) {
    all(corrected@x >= 0) && all(corrected@x <= observed@x + tol)
  } else if (!.is_sparse(observed) && !.is_sparse(corrected)) {
    all(corrected >= 0) && all(corrected <= observed + tol)
  } else {
    d <- methods::as(observed - corrected, "CsparseMatrix")
    all(d@x >= -tol) && all(.nz(corrected)$v >= 0)
  }
  if (!isTRUE(ok)) stop("internal error: mass-conservation guarantee violated", call. = FALSE)
  invisible(TRUE)
}

.drop0 <- function(x) if (.is_sparse(x)) Matrix::drop0(x) else x

## Per-cluster gene sums: genes x K dense (K is small).
.cluster_sums <- function(x, clusters, levels = sort(unique(clusters))) {
  ind <- Matrix::sparseMatrix(i = seq_along(clusters), j = match(clusters, levels),
                              x = 1, dims = c(length(clusters), length(levels)))
  out <- as.matrix(x %*% ind)
  colnames(out) <- as.character(levels)
  out
}

## log1p(CP10k) on the top-variance genes, returned DENSE but only n_hvg x N.
.lognorm_hvg <- function(x, n_hvg = 200) {
  cs <- .col_sums(x); cs[cs == 0] <- 1
  if (.is_sparse(x)) {
    N <- ncol(x)
    L <- x; L@x <- log1p(x@x / cs[rep.int(seq_len(N), diff(x@p))] * 1e4)
    mu <- .row_sums(L) / N
    L@x <- L@x^2
    v <- .row_sums(L) / N - mu^2
    hvg <- order(v, decreasing = TRUE)[seq_len(min(n_hvg, nrow(x)))]
    L <- x[hvg, , drop = FALSE]
    L@x <- log1p(L@x / cs[rep.int(seq_len(N), diff(L@p))] * 1e4)
    as.matrix(L)
  } else {
    logn <- log1p(sweep(x, 2, cs, "/") * 1e4)
    v <- apply(logn, 1, var)
    hvg <- order(v, decreasing = TRUE)[seq_len(min(n_hvg, nrow(logn)))]
    logn[hvg, , drop = FALSE]
  }
}
