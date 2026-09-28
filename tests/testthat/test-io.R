test_that("read_10x_counts + split_droplets round-trip", {
  set.seed(1)
  G <- 50; cells <- 30; empt <- 200
  cellm <- matrix(rpois(G * cells, 20), G)
  emptm <- matrix(rpois(G * empt, 0.5), G)
  raw <- methods::as(methods::as(cbind(cellm, emptm), "CsparseMatrix"), "generalMatrix")
  d <- tempfile(); dir.create(d)
  Matrix::writeMM(raw, file.path(d, "matrix.mtx"))
  bcs <- paste0("BC", seq_len(ncol(raw)))
  writeLines(bcs, file.path(d, "barcodes.tsv"))
  utils::write.table(data.frame(paste0("ENSG", 1:G), paste0("G", 1:G), "Gene Expression"),
                     file.path(d, "features.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  x <- read_10x_counts(d)
  expect_equal(dim(x), dim(raw))
  expect_equal(rownames(x)[1], "ENSG1")

  sp <- split_droplets(x, cells = bcs[1:cells])
  expect_equal(ncol(sp$cells), cells)
  expect_true(all(colSums(sp$empties) <= 100 & colSums(sp$empties) >= 1))
  expect_error(split_droplets(x, cells = "nope"), "not found")
})
