#!/usr/bin/env Rscript
## validation/benchmark_pbmc_six_tools.R
## Decant against SoupX, DecontX, CellBender, scAR and AmbiDose on the real
## 10x PBMC 4k data bundled with SoupX. No ground truth; see BENCHMARK.md.
##
##   Rscript validation/benchmark_pbmc_six_tools.R export     write the dataset for the Python tools
##   <venv>/bin/python validation/py/run_python_tools.py validation/data/bench_pbmc
##   Rscript validation/benchmark_pbmc_six_tools.R score      run the R tools and score everything
##
## The SoupX data ships the summed empty-droplet counts per gene, not the
## droplets. AmbiDose and CellBender need droplets, so 40,000 empty droplets
## are RECONSTRUCTED by assigning each of the 2.7M soup UMIs to a droplet at
## random. Their pooled profile is exactly the real soup; their droplet-to-
## droplet variation is not real. Every tool that takes empties gets the same
## reconstructed set, so that no tool is fed a different soup.

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))
source(file.path("validation", "comparators.R")); source(file.path("validation", "python_tool_io.R"))
mode <- commandArgs(trailingOnly = TRUE)[1]
DIR <- file.path("validation", "data", "bench_pbmc"); OUT <- file.path("validation", "results")

CACHE <- file.path("validation", "data", "SoupX")
if (!dir.exists(CACHE))
  try(system2("git", c("clone", "-q", "--depth", "1", "https://github.com/constantAmateur/SoupX", CACHE)))
e <- new.env()
## The same objects ship with the CRAN SoupX package; use them when
## github.com cannot be reached.
if (file.exists(file.path(CACHE, "data", "PBMC_sc.RData"))) {
  load(file.path(CACHE, "data", "PBMC_sc.RData"), e); load(file.path(CACHE, "data", "PBMC_metaData.RData"), e)
} else utils::data(list = c("PBMC_sc", "PBMC_metaData"), package = "SoupX", envir = e)
md <- e$PBMC_metaData[colnames(e$PBMC_sc$toc), ]; ok <- md$Annotation != "?"
cells <- methods::as(methods::as(e$PBMC_sc$toc, "CsparseMatrix"), "dgCMatrix")[, ok]; md <- md[ok, ]
colnames(cells) <- paste0(colnames(cells), "-1")
soup_counts <- e$PBMC_sc$soupProfile$counts
set.seed(1); N_EMPTY <- 40000
empties <- Matrix::sparseMatrix(i = rep(seq_along(soup_counts), soup_counts),
                                j = sample(N_EMPTY, sum(soup_counts), replace = TRUE), x = 1,
                                dims = c(nrow(cells), N_EMPTY),
                                dimnames = list(rownames(cells), paste0("EMPTY", seq_len(N_EMPTY), "-1")))
stopifnot(all(Matrix::rowSums(empties) == soup_counts))

if (identical(mode, "export")) {
  export_dataset(DIR, cells, empties, md$Cluster)
  cat("exported to", DIR, ": cells", ncol(cells), "empties", ncol(empties),
      "median empty UMIs", stats::median(Matrix::colSums(empties)), "\n")
  quit(status = 0)
}

soup <- soup_counts / sum(soup_counts); names(soup) <- rownames(cells)
T_c <- Matrix::colSums(cells); AMBIENT_BOUND <- 0.15
timed <- function(f) { t <- system.time(r <- f())[["elapsed"]]; r$seconds <- t; r }
res <- list(
  Decant = timed(function() run_decant(cells, empties, md$Cluster)),
  `Decant (bayes, previous default)` = timed(function() run_decant(cells, empties, md$Cluster, correction = "bayes")),
  SoupX = timed(function() run_soupx(cells, empties, md$Cluster)),
  DecontX = timed(function() run_decontx(cells, empties, md$Cluster)),
  CellBender = read_python_tool(DIR, "cellbender", cells),
  scAR = read_python_tool(DIR, "scar", cells),
  AmbiDose = read_python_tool(DIR, "ambidose", cells))
versions <- c(Decant = as.character(utils::packageVersion("Decant")),
              `Decant (bayes, previous default)` = as.character(utils::packageVersion("Decant")),
              SoupX = as.character(utils::packageVersion("SoupX")),
              DecontX = as.character(utils::packageVersion("decontX")),
              CellBender = res$CellBender$version %||% NA, scAR = res$scAR$version %||% NA,
              AmbiDose = res$AmbiDose$version %||% NA)

sets <- list(immunoglobulin = list(genes = c("IGKC", "IGHM", "IGLC2"), home = "B"),
             monocyte = list(genes = c("LYZ", "S100A8", "S100A9"), home = "MNP"),
             nk_cytotoxic = list(genes = c("GNLY", "NKG7"), home = c("NK", "T_CD8")))
rows <- lapply(names(res), function(nm) {
  r <- res[[nm]]
  if (is.null(r)) return(data.frame(method = nm, version = versions[[nm]], status = "not run"))
  if (!is.null(r$error)) return(data.frame(method = nm, version = versions[[nm]], status = paste("error:", substr(r$error, 1, 120))))
  corr <- methods::as(methods::as(r$corrected, "CsparseMatrix"), "dgCMatrix")
  row <- data.frame(method = nm, version = versions[[nm]], status = "ok",
                    rho_median = if (is.null(r$rho)) NA else stats::median(r$rho),
                    total_removed = 1 - sum(corr) / sum(cells),
                    ## count integrity: UMIs in entries where the corrected value exceeds the observed
                    counts_added = sum(pmax(corr - cells, 0)),
                    entries_increased = sum((corr - cells)@x > 1e-8),
                    seconds = r$seconds)
  for (s in names(sets)) {
    g <- sets[[s]]$genes; home <- md$Annotation %in% sets[[s]]$home
    row[[paste0(s, "_foreign_removed")]] <- 1 - sum(corr[g, !home]) / sum(cells[g, !home])
    row[[paste0(s, "_native_removed")]] <- 1 - sum(corr[g, home]) / sum(cells[g, home])
    ## Split the foreign counts by whether ambient RNA could have produced
    ## them: an entry is "incompatible" when its count is beyond the 99.9%
    ## Poisson bound of a deliberately generous ambient expectation
    ## (contamination 0.15, about twice any tool's estimate here). The bound
    ## uses no tool's own rho.
    X <- as.matrix(cells[g, !home]); C <- as.matrix(corr[g, !home])
    inc <- X > 0 & stats::ppois(X - 1, outer(soup[g], AMBIENT_BOUND * T_c[!home]), lower.tail = FALSE) < 0.001
    row[[paste0(s, "_compatible_removed")]] <- 1 - sum(C[!inc]) / sum(X[!inc])
    row[[paste0(s, "_incompatible_removed")]] <- 1 - sum(C[inc]) / sum(X[inc])
    row[[paste0(s, "_incompatible_share")]] <- sum(X[inc]) / sum(X)
  }
  row$note <- r$note %||% ""
  row
})
allc <- unique(unlist(lapply(rows, names)))
tab <- do.call(rbind, lapply(rows, function(r) { r[setdiff(allc, names(r))] <- NA; r[allc] }))
print(t(tab), quote = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, file.path(OUT, "pbmc_six_tools.csv"), row.names = FALSE)
