## seurat.R
## Run Decant on a Seurat object. Most single-cell analyses live in Seurat, but a
## Seurat object never contains the empty droplets that ambient correction
## needs, and its barcodes and feature names usually no longer match the raw
## Cell Ranger output. This wrapper exists mainly to get that bookkeeping
## right, and to fail loudly instead of silently mis-aligning cells or genes:
##
##  * features: Seurat replaces "_" with "-" and reads symbols (made unique),
##    not Ensembl IDs, so the raw matrix is read and renamed the same way;
##  * barcodes: merge()/integration add a "sample_" prefix or "_1" suffix,
##    so several standard rewrites are tried and every cell must match;
##  * split v5 layers (counts.1, counts.2, ...) are joined before use;
##  * normalised data is refused: Decant needs raw UMI counts.
##
## The corrected counts go into a NEW assay; the original assay is untouched.

#' Run Decant on a Seurat object.
#'
#' @param object a Seurat object holding raw UMI counts.
#' @param raw the raw (unfiltered) droplet matrix: a path to a Cell Ranger
#'   `raw_feature_bc_matrix` directory or a genes x barcodes matrix. For a
#'   multi-sample object, a named list of these, one per sample, with names
#'   matching the values of `meta.data[[sample_col]]`. Either `raw` or
#'   `empties` is required.
#' @param empties alternatively, the empty-droplet matrix (or named list per
#'   sample) directly.
#' @param assay assay holding the raw counts.
#' @param clusters cluster labels: the name of a meta.data column, a vector
#'   with one label per cell, or NULL to use `Idents(object)`. k-means is used
#'   only if the identities are a single level.
#' @param sample_col meta.data column naming each cell's sample; required
#'   when `raw`/`empties` is a list.
#' @param new_assay name of the assay that receives the corrected counts.
#' @param set_default make `new_assay` the default assay.
#' @param empty_umi_range UMI range defining empty droplets (see
#'   [split_droplets()]).
#' @param ... passed to [decant()] (e.g. `correction`, `rho`, `k`).
#' @return the Seurat object with: a new assay (`new_assay`) of corrected
#'   counts, meta.data column `decant_rho`, and `Misc(object, "decant")` holding
#'   the per-cluster rho diagnostics, the soup profile(s), the modules that ran
#'   and the soup-mismatch flag. Run `NormalizeData()` on the new assay before
#'   downstream use.
#' @examples
#' \dontrun{
#' obj <- RunDecant(obj, raw = "sample/outs/raw_feature_bc_matrix")
#' obj <- NormalizeData(obj, assay = "Decant")
#' }
#' @export
RunDecant <- function(object, raw = NULL, empties = NULL, assay = "RNA",
                      clusters = NULL, sample_col = NULL, new_assay = "Decant",
                      set_default = FALSE, empty_umi_range = c(1, 100), ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("RunDecant() needs the SeuratObject package", call. = FALSE)
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object", call. = FALSE)
  if (is.null(raw) == is.null(empties))
    stop("supply exactly one of `raw` (raw droplet matrix) or `empties`", call. = FALSE)
  if (!assay %in% SeuratObject::Assays(object)) stop("assay '", assay, "' not found", call. = FALSE)
  if (new_assay %in% SeuratObject::Assays(object))
    stop("assay '", new_assay, "' already exists; choose another `new_assay`", call. = FALSE)

  counts <- .seurat_counts(object[[assay]])
  cells <- colnames(object)

  ## ---- samples ----
  src <- if (!is.null(raw)) raw else empties
  multi <- is.list(src) && !is.data.frame(src) && !inherits(src, "Matrix")
  if (multi) {
    if (is.null(sample_col) || !sample_col %in% colnames(object@meta.data))
      stop("a list of raw/empties needs `sample_col`, a meta.data column naming each cell's sample",
           call. = FALSE)
    samp <- as.character(object@meta.data[[sample_col]])
    if (is.null(names(src)) || !all(unique(samp) %in% names(src)))
      stop("names of the raw/empties list must cover every value of meta.data$", sample_col,
           ": missing ", paste(setdiff(unique(samp), names(src)), collapse = ", "), call. = FALSE)
    src <- src[unique(samp)]
  } else {
    samp <- rep("all", length(cells)); src <- list(all = src)
  }

  ## ---- empties per sample, aligned to the object's features ----
  emp <- .seurat_empties(src, cells, samp, rownames(counts), is_raw = !is.null(raw),
                         empty_umi_range = empty_umi_range)

  ## ---- clusters ----
  cl <- if (is.null(clusters)) {
    id <- SeuratObject::Idents(object)
    if (nlevels(droplevels(id)) > 1) as.character(id) else NULL
  } else if (is.character(clusters) && length(clusters) == 1 &&
             clusters %in% colnames(object@meta.data)) {
    as.character(object@meta.data[[clusters]])
  } else {
    if (length(clusters) != length(cells)) stop("`clusters` must name a meta.data column or have one label per cell", call. = FALSE)
    as.character(clusters)
  }

  res <- if (multi) {
    decant(counts, emp, sample_of = match(samp, names(emp)), clusters = cl, ...)
  } else {
    decant(counts, emp[[1]], clusters = cl, ...)
  }

  ## ---- write back ----
  corrected <- res$corrected
  dimnames(corrected) <- dimnames(counts)
  object[[new_assay]] <- SeuratObject::CreateAssayObject(counts = methods::as(corrected, "CsparseMatrix"))
  object[[new_assay]] <- .assay_with_key(object[[new_assay]], new_assay)
  object$decant_rho <- res$rho
  SeuratObject::Misc(object, "decant") <- list(
    rho_diagnostics = res$rho_diagnostics, soup_mismatch = res$soup_mismatch,
    ambient = res$ambient, modules = res$modules, source_assay = assay)
  if (set_default) SeuratObject::DefaultAssay(object) <- new_assay
  object
}

## Empty-droplet matrices per sample, aligned to `genes`. `src` is a named
## list (one entry per sample) of raw matrices / Cell Ranger paths
## (is_raw = TRUE: empties are the non-cell droplets in empty_umi_range) or of
## empty-droplet matrices (is_raw = FALSE).
.seurat_empties <- function(src, cells, samp, genes, is_raw, empty_umi_range) {
  emp <- lapply(names(src), function(s) {
    x <- src[[s]]
    if (is.character(x)) x <- read_10x_counts(x, gene_column = 2)
    x <- .as_counts(x, if (is_raw) "raw" else "empties")
    rownames(x) <- .seurat_feature_names(rownames(x))
    x <- .align_features(x, genes, s)
    if (!is_raw) return(x)
    idx <- .match_barcodes(cells[samp == s], colnames(x), s)
    tot <- .col_sums(x)
    is_emp <- !(seq_len(ncol(x)) %in% idx) &
      tot >= empty_umi_range[1] & tot <= empty_umi_range[2]
    if (sum(is_emp) < 100)
      warning("sample '", s, "': only ", sum(is_emp), " empty droplets; the soup ",
              "estimate will be noisy. Is `raw` the UNFILTERED matrix?", call. = FALSE)
    x[, is_emp, drop = FALSE]
  })
  names(emp) <- names(src)
  emp
}

## Raw counts from an Assay or Assay5, joining split v5 layers.
.seurat_counts <- function(a) {
  if (inherits(a, "Assay5")) {
    lay <- grep("^counts", SeuratObject::Layers(a), value = TRUE)
    if (length(lay) == 0) stop("assay has no 'counts' layer", call. = FALSE)
    if (length(lay) > 1) a <- SeuratObject::JoinLayers(a, layers = "counts")
    m <- SeuratObject::LayerData(a, layer = "counts")
  } else {
    m <- SeuratObject::GetAssayData(a, slot = "counts")
  }
  if (prod(dim(m)) == 0) stop("the counts layer is empty", call. = FALSE)
  v <- if (inherits(m, "sparseMatrix")) m@x else m[m != 0]
  if (any(v < 0) || any(abs(v - round(v)) > 1e-8))
    stop("the counts layer is not raw UMI counts (found negative or non-integer ",
         "values). Decant must run on raw counts, before NormalizeData/SCTransform.",
         call. = FALSE)
  .as_counts(m)
}

## Mirror Seurat's feature renaming: "_" is not allowed and becomes "-".
.seurat_feature_names <- function(x) make.unique(gsub("_", "-", x, fixed = TRUE))

.align_features <- function(x, genes, sample) {
  miss <- setdiff(genes, rownames(x))
  if (length(miss) > 0)
    stop("sample '", sample, "': ", length(miss), " of the object's features are not in the ",
         "raw matrix (e.g. ", paste(utils::head(miss, 3), collapse = ", "), "). ",
         "Was the object built from the same Cell Ranger run, with gene symbols?",
         call. = FALSE)
  x[genes, , drop = FALSE]
}

## Map object barcodes onto raw barcodes. merge(add.cell.ids=) adds a
## "sample_" prefix and integration adds a "_1" suffix; try the standard
## rewrites and keep the one that matches every cell. Partial matches error.
.match_barcodes <- function(obj_bcs, raw_bcs, sample) {
  rewrites <- list(
    identity = function(b) b,
    strip_prefix = function(b) sub("^.*_", "", b),
    strip_suffix = function(b) sub("_[0-9]+$", "", b),
    strip_both = function(b) sub("^.*_", "", sub("_[0-9]+$", "", b)))
  hits <- vapply(rewrites, function(f) sum(f(obj_bcs) %in% raw_bcs), numeric(1))
  best <- names(which.max(hits))
  idx <- match(rewrites[[best]](obj_bcs), raw_bcs)
  if (anyNA(idx) || anyDuplicated(idx))
    stop("sample '", sample, "': could not match ", sum(is.na(idx)), " of ", length(obj_bcs),
         " cell barcodes to the raw matrix (e.g. ",
         paste(utils::head(obj_bcs[is.na(idx)], 3), collapse = ", "),
         "). For a merged object, pass one raw matrix per sample with `sample_col`.",
         call. = FALSE)
  idx
}

.assay_with_key <- function(a, name) {
  SeuratObject::Key(a) <- paste0(tolower(gsub("[^A-Za-z0-9]", "", name)), "_")
  a
}
