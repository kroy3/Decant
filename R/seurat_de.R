## seurat_de.R
## Contamination-aware differential expression from a Seurat object.

#' Contamination-aware DE for one cell type of a multi-sample Seurat object.
#'
#' Pseudobulks the chosen cell type per sample, measures each sample's soup
#' from its own empty droplets, and tests a sample-level design with
#' [ambient_de()]. Unlike correcting counts and then testing, this stops an
#' abundant cell type's condition-dependent genes from being reported as DE
#' in the target type when ambient RNA differs between samples.
#'
#' @param object a Seurat object with raw UMI counts and cells from several
#'   samples.
#' @param cell_type the cell type to test (a value of `celltype_col`).
#' @param sample_col meta.data column naming each cell's sample.
#' @param formula sample-level design, e.g. `~ condition` or
#'   `~ batch + condition`. Every variable must be a meta.data column that is
#'   constant within each sample.
#' @param raw named list, one entry per sample (names = values of
#'   `sample_col`), of raw droplet matrices or Cell Ranger
#'   `raw_feature_bc_matrix` paths. Or supply `empties` instead.
#' @param empties named list of empty-droplet matrices per sample.
#' @param celltype_col meta.data column of cell-type labels; NULL uses
#'   `Idents(object)`.
#' @param coef design column to test; default the last one (the last term of
#'   `formula`).
#' @param assay assay holding raw counts.
#' @param min_cells samples with fewer target cells are dropped (with a
#'   warning).
#' @param empty_umi_range UMI range defining empty droplets.
#' @param ... passed to [ambient_de()].
#' @return a data.frame as from [ambient_de()], with `ambient_driven` marking
#'   genes a standard edgeR pseudobulk test calls DE but the contamination-
#'   aware test does not. Attributes `"samples"` (per-sample ambient fraction,
#'   soup depth, cells) and `"design"`.
#' @examples
#' \dontrun{
#' res <- DecantDE(obj, cell_type = "Microglia", sample_col = "donor",
#'                 formula = ~ diagnosis, raw = raw_paths_by_donor)
#' head(res)
#' subset(res, ambient_driven)   # standard-pipeline hits explained by soup
#' }
#' @export
DecantDE <- function(object, cell_type, sample_col, formula, raw = NULL, empties = NULL,
                     celltype_col = NULL, coef = NULL, assay = "RNA", min_cells = 10,
                     empty_umi_range = c(1, 100), ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("DecantDE() needs the SeuratObject package", call. = FALSE)
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object", call. = FALSE)
  if (is.null(raw) == is.null(empties))
    stop("supply exactly one of `raw` or `empties` (named lists, one per sample)", call. = FALSE)
  md <- object@meta.data
  if (!sample_col %in% colnames(md)) stop("meta.data has no column '", sample_col, "'", call. = FALSE)
  if (!inherits(formula, "formula")) stop("`formula` must be a formula such as ~ condition", call. = FALSE)

  ct <- if (is.null(celltype_col)) as.character(SeuratObject::Idents(object)) else {
    if (!celltype_col %in% colnames(md)) stop("meta.data has no column '", celltype_col, "'", call. = FALSE)
    as.character(md[[celltype_col]])
  }
  if (!cell_type %in% ct) stop("cell type '", cell_type, "' not found", call. = FALSE)
  samp <- as.character(md[[sample_col]])

  ## ---- samples with enough target cells ----
  n_target <- table(factor(samp[ct == cell_type], levels = unique(samp)))
  use <- names(n_target)[n_target >= min_cells]
  if (length(use) < length(n_target))
    warning("dropping ", length(n_target) - length(use), " sample(s) with < ", min_cells,
            " '", cell_type, "' cells: ", paste(setdiff(names(n_target), use), collapse = ", "),
            call. = FALSE)

  ## ---- sample-level design ----
  vars <- all.vars(formula)
  miss <- setdiff(vars, colnames(md))
  if (length(miss)) stop("formula variables not in meta.data: ", paste(miss, collapse = ", "), call. = FALSE)
  sdat <- lapply(vars, function(v) {
    vals <- tapply(as.character(md[[v]]), samp, function(x) unique(x))
    bad <- names(vals)[lengths(vals) != 1]
    if (length(bad)) stop("'", v, "' is not constant within sample(s): ",
                          paste(utils::head(bad, 3), collapse = ", "), call. = FALSE)
    unlist(vals)[use]
  })
  sdat <- as.data.frame(stats::setNames(sdat, vars), stringsAsFactors = TRUE)
  rownames(sdat) <- use
  design <- stats::model.matrix(formula, sdat)
  if (qr(design)$rank < ncol(design)) stop("the sample-level design is not of full rank", call. = FALSE)
  if (nrow(design) - ncol(design) < 2)
    stop("too few samples for the design (need >= 2 residual degrees of freedom)", call. = FALSE)
  if (is.null(coef)) coef <- ncol(design)

  ## ---- pseudobulk of the target type ----
  counts <- .seurat_counts(object[[assay]])
  sel <- ct == cell_type & samp %in% use
  Y <- .cluster_sums(counts[, sel, drop = FALSE], factor(samp[sel], levels = use), use)

  ## ---- empties per sample ----
  src <- if (!is.null(raw)) raw else empties
  if (is.null(names(src)) || !all(use %in% names(src)))
    stop("raw/empties must be a list named by sample, covering: ",
         paste(setdiff(use, names(src)), collapse = ", "), call. = FALSE)
  emp <- .seurat_empties(src[use], colnames(object), samp, rownames(counts),
                         is_raw = !is.null(raw), empty_umi_range = empty_umi_range)

  res <- ambient_de(Y, emp, design, coef = coef, ...)
  smp <- attr(res, "samples")
  smp$sample <- use; smp$n_cells <- as.integer(n_target[use])
  attr(res, "samples") <- smp
  attr(res, "design") <- design
  attr(res, "cell_type") <- cell_type
  res
}
