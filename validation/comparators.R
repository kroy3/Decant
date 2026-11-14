## validation/comparators.R
## Wrappers around the OFFICIAL competing tools, for benchmarks only.
##
## Contract (same as validation/species_mix.R): every wrapper takes
##   (cells, empties, clusters)   genes x droplets count matrices + cell labels
## and returns one of
##   NULL                          the tool is not installed  -> "not installed"
##   list(error = "<message>")     the tool was run and failed -> "error: ..."
##   list(corrected, rho, note)    corrected genes x cells counts; rho is the
##                                 per-cell contamination estimate, or NULL
##                                 when the tool does not report one.
## A missing tool is never replaced by an in-repo reimplementation.
##
## Install:
##   SoupX       install.packages("SoupX")
##   DecontX     BiocManager::install("decontX")      # or the older "celda"
##   FastCAR     remotes::install_github("LungCellAtlas/FastCAR")
##   CellBender  pip install cellbender   (a Python CLI; also needs the R
##               package hdf5r to read its output)
##
## Status of each wrapper is recorded in validation/README.md. Do not quote a
## number from a wrapper listed there as "not yet run".

.as_dgc <- function(x) {
  if (methods::is(x, "dgCMatrix")) return(x)
  methods::as(methods::as(methods::as(Matrix::Matrix(as.matrix(x), sparse = TRUE),
                                      "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

## Droplets need unique names for every tool that binds cells and empties.
.name_droplets <- function(cells, empties) {
  if (is.null(colnames(cells))) colnames(cells) <- paste0("cell", seq_len(ncol(cells)))
  if (is.null(colnames(empties)) || anyDuplicated(colnames(empties)) ||
      any(colnames(empties) %in% colnames(cells)))
    colnames(empties) <- paste0("empty", seq_len(ncol(empties)))
  if (is.null(rownames(empties))) rownames(empties) <- rownames(cells)
  list(cells = cells, empties = empties)
}

run_decant <- function(cells, empties, clusters, correction = "bayes") {
  tryCatch({
    res <- suppressWarnings(decant(cells, empties, clusters = clusters, correction = correction))
    list(corrected = res$corrected, rho = res$rho,
         note = if (isTRUE(res$soup_mismatch)) "soup_mismatch warning" else "")
  }, error = function(e) list(error = conditionMessage(e)))
}

run_soupx <- function(cells, empties, clusters) {
  if (!requireNamespace("SoupX", quietly = TRUE)) return(NULL)
  tryCatch({
    d <- .name_droplets(.as_dgc(cells), .as_dgc(empties))
    tod <- cbind(d$cells, d$empties)
    sc <- SoupX::SoupChannel(tod, d$cells, calcSoupProfile = FALSE)
    ## SoupX's default soupRange is c(0, 100); the upper bound is raised only
    ## if the supplied empties are larger than that, so that the droplets the
    ## caller has already declared empty are the ones used.
    hi <- max(100, max(Matrix::colSums(d$empties)) + 1)
    sc <- SoupX::estimateSoup(sc, soupRange = c(0, hi))
    sc <- SoupX::setClusters(sc, stats::setNames(as.character(clusters), colnames(d$cells)))
    ## Defaults first. If autoEstCont finds no marker genes it says to "reduce
    ## tfidfMin or soupQuantile"; follow that once, with one fixed setting
    ## (never tuned to the truth), and record that it happened.
    est <- function(...) suppressWarnings(SoupX::autoEstCont(sc, doPlot = FALSE, verbose = FALSE, ...))
    note <- ""
    sc <- tryCatch(est(), error = function(e) {
      if (!grepl("No plausible marker genes", conditionMessage(e))) stop(e)
      note <<- "defaults found no markers; soupQuantile=0.5"
      est(soupQuantile = 0.5)
    })
    list(corrected = SoupX::adjustCounts(sc, verbose = 0),
         rho = sc$metaData$rho, note = note)
  }, error = function(e) list(error = conditionMessage(e)))
}

run_decontx <- function(cells, empties, clusters) {
  pkg <- if (requireNamespace("decontX", quietly = TRUE)) "decontX"
         else if (requireNamespace("celda", quietly = TRUE)) "celda" else NULL
  if (is.null(pkg)) return(NULL)
  tryCatch({
    d <- .name_droplets(.as_dgc(cells), .as_dgc(empties))
    f <- utils::getFromNamespace("decontX", pkg)
    r <- f(x = d$cells, z = as.character(clusters), background = d$empties, verbose = FALSE)
    list(corrected = r$decontXcounts, rho = r$contamination, note = pkg)
  }, error = function(e) list(error = conditionMessage(e)))
}

## FastCAR removes, for each gene seen in more than `chance_cutoff` of the
## empty droplets, the highest count that gene reaches in any empty droplet,
## from every cell. It reports no per-cell contamination fraction.
## `empty_cutoff` is FastCAR's emptyDropletCutoff: droplets with fewer UMIs
## are treated as empty. Defaults are the values in the FastCAR README.
run_fastcar <- function(cells, empties, clusters, empty_cutoff = 100, chance_cutoff = 0.05) {
  if (!requireNamespace("FastCAR", quietly = TRUE)) return(NULL)
  tryCatch({
    d <- .name_droplets(.as_dgc(cells), .as_dgc(empties))
    full <- cbind(d$cells, d$empties)
    amb <- FastCAR::determine.background.to.remove(full, d$cells, empty_cutoff, chance_cutoff)
    list(corrected = FastCAR::remove.background(d$cells, amb), rho = NULL,
         note = sprintf("%d genes corrected", sum(amb > 0)))
  }, error = function(e) list(error = conditionMessage(e)))
}

## CellBender is a Python CLI. The droplets are written as a Cell Ranger v3
## MTX folder, `cellbender remove-background` is run on it, and the output h5
## is read back for the supplied cell barcodes.
## NOT YET RUN: this wrapper was written against the CellBender v0.3
## documentation and has not been executed. Treat its first results as a test
## of the wrapper, not of CellBender.
run_cellbender <- function(cells, empties, clusters, epochs = 150, cuda = FALSE,
                           extra_args = character(0)) {
  exe <- Sys.which("cellbender")
  if (!nzchar(exe) || !requireNamespace("hdf5r", quietly = TRUE)) return(NULL)
  tryCatch({
    d <- .name_droplets(.as_dgc(cells), .as_dgc(empties))
    full <- cbind(d$cells, d$empties)
    dir <- tempfile("cellbender_in_"); dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    gz <- function(lines, name) {
      con <- gzfile(file.path(dir, name), "w"); on.exit(close(con))
      writeLines(lines, con)
    }
    mm <- file.path(dir, "matrix.mtx")
    Matrix::writeMM(full, mm)
    gz(readLines(mm), "matrix.mtx.gz"); unlink(mm)
    gz(paste(rownames(full), rownames(full), "Gene Expression", sep = "\t"), "features.tsv.gz")
    gz(colnames(full), "barcodes.tsv.gz")
    out <- file.path(dir, "out.h5")
    args <- c("remove-background", "--input", dir, "--output", out,
              "--expected-cells", ncol(d$cells),
              "--total-droplets-included", ncol(full),
              "--epochs", epochs, if (cuda) "--cuda", extra_args)
    status <- system2(exe, args, stdout = FALSE, stderr = FALSE)
    if (status != 0 || !file.exists(out)) stop("cellbender exited with status ", status)
    h <- hdf5r::H5File$new(out, mode = "r"); on.exit(h$close_all(), add = TRUE)
    m <- h[["matrix"]]
    shape <- m[["shape"]]$read()
    X <- Matrix::sparseMatrix(i = m[["indices"]]$read() + 1L, p = m[["indptr"]]$read(),
                              x = as.numeric(m[["data"]]$read()), dims = shape)
    colnames(X) <- m[["barcodes"]]$read()
    rownames(X) <- m[["features"]][["id"]]$read()
    miss <- setdiff(colnames(d$cells), colnames(X))
    if (length(miss)) stop(length(miss), " cell barcodes missing from the CellBender output")
    list(corrected = X[rownames(d$cells), colnames(d$cells), drop = FALSE], rho = NULL,
         note = "wrapper not yet validated")
  }, error = function(e) list(error = conditionMessage(e)))
}

COMPARATORS <- list(Decant = run_decant, SoupX = run_soupx, DecontX = run_decontx,
                    FastCAR = run_fastcar, CellBender = run_cellbender)
