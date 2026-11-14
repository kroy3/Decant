## validation/benchmark_export.R
## Write one dataset in the layout validation/py/run_python_tools.py reads, and
## read the Python tools' outputs back.

## 10x v3 MTX folder with all droplets (cells first, then empties).
export_dataset <- function(dir, cells, empties, clusters) {
  dir.create(file.path(dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  if (is.null(colnames(empties))) colnames(empties) <- paste0("EMPTY", seq_len(ncol(empties)), "-1")
  full <- cbind(methods::as(cells, "CsparseMatrix"), methods::as(empties, "CsparseMatrix"))
  gz <- function(lines, name) { con <- gzfile(file.path(dir, "raw", name), "w"); writeLines(lines, con); close(con) }
  mm <- file.path(dir, "raw", "matrix.mtx"); Matrix::writeMM(full, mm)
  gz(readLines(mm), "matrix.mtx.gz"); unlink(mm)
  gz(paste(rownames(full), rownames(full), "Gene Expression", sep = "\t"), "features.tsv.gz")
  gz(colnames(full), "barcodes.tsv.gz")
  writeLines(colnames(cells), file.path(dir, "cells.txt"))
  writeLines(as.character(ncol(full)), file.path(dir, "n_droplets.txt"))
  utils::write.csv(data.frame(barcode = colnames(cells), cluster = as.character(clusters)),
                   file.path(dir, "clusters.csv"), row.names = FALSE)
  invisible(dir)
}

## list(corrected, rho, note, seconds, version) or list(error = )
read_python_tool <- function(dir, tool, cells) {
  js <- file.path(dir, "out", paste0(tool, ".json"))
  if (!file.exists(js)) return(NULL)
  meta <- jsonlite::fromJSON(js)
  if (!is.null(meta$error)) return(list(error = meta$error, seconds = meta$seconds))
  X <- methods::as(methods::as(Matrix::readMM(file.path(dir, "out", paste0(tool, ".mtx"))), "CsparseMatrix"), "dgCMatrix")
  dimnames(X) <- dimnames(cells)
  list(corrected = X, rho = if (length(meta$rho)) as.numeric(meta$rho) else NULL, note = meta$note,
       seconds = meta$seconds, version = meta$version)
}
