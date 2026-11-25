## Streaming input (I5): the soup statistics and the cells must match the
## in-memory path exactly, whatever the chunk size.
write_raw_10x <- function(seed = 1, G = 60, cells = 40, empt = 400, gz = TRUE) {
  set.seed(seed)
  cellm <- matrix(stats::rpois(G * cells, 15), G)
  emptm <- matrix(stats::rpois(G * empt, 0.3), G)
  emptm[, 1:20] <- 0                                     # zero-count barcodes
  raw <- methods::as(methods::as(cbind(emptm[, 1:200], cellm, emptm[, -(1:200)]), "CsparseMatrix"), "generalMatrix")
  d <- tempfile(); dir.create(d)
  Matrix::writeMM(raw, file.path(d, "matrix.mtx"))
  bcs <- paste0("BC", seq_len(ncol(raw)))
  writeLines(bcs, file.path(d, "barcodes.tsv"))
  ftype <- c(rep("Gene Expression", G - 3), rep("Antibody Capture", 3))
  utils::write.table(data.frame(paste0("ENSG", 1:G), paste0("G", 1:G), ftype),
                     file.path(d, "features.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  if (gz) for (f in c("matrix.mtx", "barcodes.tsv", "features.tsv")) {
    x <- readBin(file.path(d, f), "raw", file.size(file.path(d, f)))
    con <- gzfile(file.path(d, paste0(f, ".gz")), "wb"); writeBin(x, con); close(con)
    unlink(file.path(d, f))
  }
  list(dir = d, cells = bcs[200 + seq_len(cells)])
}

test_that("read_10x_streaming matches read_10x_counts + split_droplets to 1e-8", {
  w <- write_raw_10x()
  ref <- split_droplets(read_10x_counts(w$dir), cells = w$cells)
  for (ch in c(37, 1e6)) {
    st <- read_10x_streaming(w$dir, cells = w$cells, chunk_lines = ch)
    expect_equal(st$summary, ref$summary)
    expect_equal(as.matrix(st$cells), as.matrix(ref$cells), tolerance = 0)
    expect_equal(dimnames(st$cells), dimnames(ref$cells))
    sref <- soup_summary(ref$empties)
    for (f in c("gene_sums", "gene_sumsq", "gene_detected", "total_umis", "n_droplets"))
      expect_equal(st$soup[[f]], sref[[f]], tolerance = 1e-8, label = f)
    expect_equal(ambient_global(st$soup), ambient_global(ref$empties), tolerance = 1e-12)
  }
  ## without cell calls, the UMI threshold (here every cell has > 500 UMIs)
  st2 <- read_10x_streaming(w$dir, chunk_lines = 100)
  expect_equal(colnames(st2$cells), colnames(split_droplets(read_10x_counts(w$dir))$cells))
})

test_that("decant() gives the same result from a soup summary as from the empties", {
  sim <- small_sim(rho_mean = 0.1)
  cl <- quick_labels(sim$observed, k = 4)
  a <- decant(as_sparse(sim$observed), as_sparse(sim$empty), clusters = cl)
  b <- decant(as_sparse(sim$observed), soup_summary(sim$empty), clusters = cl)
  expect_equal(b$rho, a$rho)
  expect_equal(as.matrix(b$corrected), as.matrix(a$corrected))
  ## multi-sample: matrices and summaries mixed
  s <- rep(1:2, length.out = ncol(sim$observed))
  e <- list(as_sparse(sim$empty[, 1:700]), as_sparse(sim$empty[, 701:1500]))
  m1 <- decant(as_sparse(sim$observed), e, sample_of = s, clusters = cl)
  m2 <- decant(as_sparse(sim$observed), list(soup_summary(e[[1]]), e[[2]]), sample_of = s, clusters = cl)
  expect_equal(as.matrix(m2$corrected), as.matrix(m1$corrected))
  expect_error(decant(sim$observed, soup_summary(sim$empty), clusters = cl, basis = diag(3)), "soup summaries")
})

test_that("soup_dispersion recovers NB dispersion from moments", {
  set.seed(2)
  G <- 300; n <- 5000; mu <- stats::rgamma(G, 2, 4); alpha <- 0.3
  e <- matrix(stats::rnbinom(G * n, size = 1 / alpha, mu = mu), G)
  d <- soup_dispersion(soup_summary(e))
  ok <- mu > 0.2
  expect_lt(abs(stats::median(attr(d, "raw")[ok]) - alpha), 0.05)
  expect_lt(abs(stats::median(d[ok]) - alpha), 0.05)
})

test_that("read_10x_streaming does not depend on the entry order of the file", {
  w <- write_raw_10x(seed = 3, gz = FALSE)
  f <- file.path(w$dir, "matrix.mtx")
  l <- readLines(f); hdr <- grep("^%", l); body <- l[-c(hdr, max(hdr) + 1)]
  set.seed(9); writeLines(c(l[c(hdr, max(hdr) + 1)], sample(body)), f)
  ref <- split_droplets(read_10x_counts(w$dir), cells = w$cells)
  st <- read_10x_streaming(w$dir, cells = w$cells, chunk_lines = 50)
  expect_true(methods::validObject(st$cells))
  expect_equal(as.matrix(st$cells), as.matrix(ref$cells), tolerance = 0)
  expect_equal(st$soup$gene_sums, soup_summary(ref$empties)$gene_sums)
})
