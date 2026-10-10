skip_if_not_installed("SeuratObject")
suppressPackageStartupMessages(library(SeuratObject))

## A raw droplet matrix (cells + empties) whose gene names contain "_", which
## Seurat renames to "-", plus the Seurat object a user would have built.
make_case <- function(seed = 1, prefix = NULL) {
  sim <- simulate_experiment(n_genes = 300, n_types = 4, n_cells = 300, n_empty = 1500,
                             rho_mean = 0.1, seed = seed)
  genes <- paste0("GENE_", seq_len(nrow(sim$observed)))
  rownames(sim$observed) <- rownames(sim$empty) <- genes
  bc <- paste0("BC", seed, "_", seq_len(ncol(sim$observed)), "-1")
  bc <- sub("_", "x", bc)                      # raw barcodes contain no "_"
  colnames(sim$observed) <- bc
  colnames(sim$empty) <- paste0("EMPTY", seed, "x", seq_len(ncol(sim$empty)))
  raw <- methods::as(cbind(sim$observed, sim$empty), "CsparseMatrix")
  cells <- methods::as(sim$observed, "CsparseMatrix")
  if (!is.null(prefix)) colnames(cells) <- paste0(prefix, "_", colnames(cells))
  obj <- suppressWarnings(CreateSeuratObject(cells))
  obj$type <- paste0("t", sim$labels)
  Idents(obj) <- "type"
  list(sim = sim, raw = raw, obj = obj)
}

test_that("RunDecant matches decant() exactly and writes a new assay", {
  cs <- make_case()
  obj <- RunDecant(cs$obj, raw = cs$raw)
  expect_true("Decant" %in% Assays(obj))
  expect_equal(DefaultAssay(obj), "RNA")
  expect_length(obj$decant_rho, ncol(obj))
  expect_false(is.null(Misc(obj, "decant")$rho_diagnostics))

  ## same answer as the matrix API on the same cells, empties and clusters
  ref <- decant(cs$sim$observed, cs$sim$empty, clusters = paste0("t", cs$sim$labels))
  got <- LayerData(obj[["Decant"]], layer = "counts")
  expect_equal(unname(as.matrix(got)), unname(as.matrix(ref$corrected)))
  expect_equal(unname(obj$decant_rho), ref$rho)
  ## original counts untouched, never fabricated
  orig <- LayerData(obj[["RNA"]], layer = "counts")
  expect_true(all(as.matrix(got) <= as.matrix(orig)))
})

test_that("feature renaming (_ -> -) and filtered features are handled", {
  cs <- make_case()
  expect_true(all(grepl("^GENE-", rownames(cs$obj))))
  keep <- rownames(cs$obj)[seq(1, nrow(cs$obj), by = 2)]   # object with a gene subset
  obj <- RunDecant(subset(cs$obj, features = keep), raw = cs$raw)
  expect_equal(nrow(obj[["Decant"]]), length(keep))
})

test_that("prefixed barcodes from merge(add.cell.ids=) are matched", {
  cs <- make_case(prefix = "S1")
  obj <- RunDecant(cs$obj, raw = cs$raw)
  expect_true("Decant" %in% Assays(obj))
})

test_that("multi-sample merged object with split v5 layers", {
  a <- make_case(seed = 1); b <- make_case(seed = 2)
  m <- merge(a$obj, b$obj, add.cell.ids = c("A", "B"))
  m$sample <- ifelse(startsWith(colnames(m), "A_"), "A", "B")
  expect_gt(length(Layers(m[["RNA"]])), 1)                 # counts.1, counts.2
  out <- RunDecant(m, raw = list(A = a$raw, B = b$raw), sample_col = "sample",
                   clusters = "type")
  expect_equal(ncol(out[["Decant"]]), ncol(m))
  expect_true("hierarchical_ambient" %in% Misc(out, "decant")$modules)
  expect_equal(ncol(Misc(out, "decant")$ambient), 2)
})

test_that("unsafe inputs fail loudly", {
  cs <- make_case()
  norm <- suppressWarnings(CreateSeuratObject(cs$sim$observed / 7.3))
  expect_error(RunDecant(norm, raw = cs$raw), "raw UMI counts")
  expect_error(RunDecant(cs$obj), "exactly one")
  wrong <- cs$raw; colnames(wrong) <- paste0("zz", colnames(wrong))
  expect_error(RunDecant(cs$obj, raw = wrong), "could not match")
  fewer <- cs$raw[-(1:5), ]
  expect_error(RunDecant(cs$obj, raw = fewer), "not in the raw matrix")
  obj <- RunDecant(cs$obj, raw = cs$raw)
  expect_error(RunDecant(obj, raw = cs$raw), "already exists")
})

test_that("raw given as a Cell Ranger directory (symbols, as Read10X uses)", {
  cs <- make_case()
  d <- tempfile(); dir.create(d)
  Matrix::writeMM(cs$raw, file.path(d, "matrix.mtx"))
  writeLines(colnames(cs$raw), file.path(d, "barcodes.tsv"))
  utils::write.table(data.frame(paste0("ENSG", seq_len(nrow(cs$raw))), rownames(cs$raw),
                                "Gene Expression"),
                     file.path(d, "features.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  from_dir <- RunDecant(cs$obj, raw = d)
  from_mat <- RunDecant(cs$obj, raw = cs$raw)
  expect_equal(LayerData(from_dir[["Decant"]], layer = "counts"),
               LayerData(from_mat[["Decant"]], layer = "counts"))
})

test_that("DecantDE runs end to end from a Seurat object and flags leakage", {
  skip_if_not_installed("edgeR")
  st <- simulate_ambient_study(seed = 4, G = 500, n_cells = 400, n_per = 4)
  rownames(st$X) <- paste0("GENE_", seq_len(nrow(st$X)))           # Seurat renames "_" -> "-"
  emp <- lapply(st$empties, function(e) { rownames(e) <- rownames(st$X); e })
  bc <- paste0("S", st$sample, "x", seq_len(ncol(st$X)))
  colnames(st$X) <- bc
  raw <- lapply(seq_along(emp), function(s) {
    e <- emp[[s]]; colnames(e) <- paste0("S", s, "e", seq_len(ncol(e)))
    methods::as(cbind(st$X[, st$sample == s], e), "CsparseMatrix") })
  names(raw) <- paste0("d", seq_along(raw))
  obj <- suppressWarnings(CreateSeuratObject(methods::as(st$X, "CsparseMatrix")))
  obj$donor <- paste0("d", st$sample)
  obj$status <- as.character(st$cond)[st$sample]
  obj$type <- st$label
  res <- DecantDE(obj, cell_type = "R", sample_col = "donor", formula = ~ status,
                  raw = raw, celltype_col = "type")
  expect_true(all(c("PValue", "FDR", "ambient_driven") %in% names(res)))
  gi <- as.integer(sub("GENE-", "", res$gene))
  expect_equal(sum(res$FDR[gi %in% st$de_abund] < 0.1), 0)
  expect_equal(nrow(attr(res, "samples")), 8)
  ## design variable not constant within sample
  obj$bad <- sample(c("a", "b"), ncol(obj), TRUE)
  expect_error(DecantDE(obj, "R", "donor", ~ bad, raw = raw, celltype_col = "type"),
               "not constant")
})
