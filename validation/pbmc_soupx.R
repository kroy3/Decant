#!/usr/bin/env Rscript
## validation/pbmc_soupx.R
## Real-data check on the 10x PBMC 4k data bundled with SoupX (Young & Behjati
## 2020): 2,170 cells, 33,694 genes, a soup profile from 2.7M empty-droplet
## UMIs, and the authors' cluster annotation. Fetched from the SoupX GitHub
## repository, so it needs no host beyond github.com.
##
## There is NO ground truth here. What can be checked:
##   1. rho against SoupX's own estimates on the same data, as printed in its
##      tutorial (inst/doc/pbmcTutorial.html): 0.06 from autoEstCont, 7.69%
##      from the immunoglobulin genes in non-B cells.
##   2. Whether rho * soup predicts the marker counts seen in cell types that
##      should not express them (immunoglobulins outside B cells, LYZ /
##      S100A8 / S100A9 outside monocytes, GNLY / NKG7 outside NK and CD8 T).
##   3. How much of those counts each correction rule removes, and how much
##      it removes from the cell type that does express them.
##
## Run from the repo root:  Rscript validation/pbmc_soupx.R

if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
suppressPackageStartupMessages(library(Matrix))

CACHE <- file.path("validation", "data", "SoupX")
OUT   <- file.path("validation", "results")
if (!dir.exists(CACHE))
  system2("git", c("clone", "-q", "--depth", "1", "https://github.com/constantAmateur/SoupX", CACHE))
e <- new.env()
load(file.path(CACHE, "data", "PBMC_sc.RData"), e)
load(file.path(CACHE, "data", "PBMC_metaData.RData"), e)
sc <- e$PBMC_sc; md <- e$PBMC_metaData

cells <- methods::as(methods::as(sc$toc, "CsparseMatrix"), "dgCMatrix")
md <- md[colnames(cells), ]
ok <- md$Annotation != "?"
cells <- cells[, ok]; md <- md[ok, ]
## the SoupChannel keeps the summed empty-droplet counts per gene, not the
## droplets; one pooled column is all decant() needs for a single sample
empties <- Matrix::Matrix(matrix(sc$soupProfile$counts, ncol = 1,
                                 dimnames = list(rownames(cells), "pooled")), sparse = TRUE)

res <- withCallingHandlers(decant(cells, empties, clusters = md$Cluster),
                           warning = function(w) { message("decant warning: ", conditionMessage(w))
                                                   invokeRestart("muffleWarning") })
cat(sprintf("cells %d, genes %d, soup UMIs %d\n", ncol(cells), nrow(cells), sum(empties)))
cat(sprintf("Decant rho: median %.3f, mean %.3f, UMI-weighted %.3f   (SoupX tutorial: 0.06 auto, 0.077 IG genes)\n",
            stats::median(res$rho), mean(res$rho), stats::weighted.mean(res$rho, colSums(cells))))
cat(sprintf("fraction of UMIs removed: %.3f   soup_mismatch: %s\n",
            1 - sum(res$corrected) / sum(cells), isTRUE(res$soup_mismatch)))
print(res$rho_diagnostics[, c("cluster", "n_cells", "rho", "n_absent_genes", "low_power", "borrowed")],
      row.names = FALSE, digits = 3)

## ---- lineage-restricted genes: removed where foreign, kept where native ----
## For each marker set: counts in the cell types that should not express it
## ("foreign") and in the type that does ("native"). `expected_ambient` is
## rho * UMIs * soup summed over the foreign cells, as a share of what is
## observed there: near 1 means the ambient model accounts for those counts.
sets <- list(immunoglobulin = list(genes = c("IGKC", "IGHM", "IGLC2"), home = "B"),
             monocyte = list(genes = c("LYZ", "S100A8", "S100A9"), home = "MNP"),
             nk_cytotoxic = list(genes = c("GNLY", "NKG7"), home = c("NK", "T_CD8")))
soup <- empties[, 1] / sum(empties)
T_c <- Matrix::colSums(cells)
rows <- list()
for (m in c("bayes", "redistribute", "subtract")) {
  corr <- if (m == "bayes") res$corrected
          else suppressWarnings(decant(cells, empties, clusters = md$Cluster, correction = m))$corrected
  for (nm in names(sets)) {
    g <- intersect(sets[[nm]]$genes, rownames(cells)); home <- md$Annotation %in% sets[[nm]]$home
    rows[[length(rows) + 1]] <- data.frame(
      correction = m, set = nm,
      foreign_umis = sum(cells[g, !home]),
      expected_ambient = sum(res$rho[!home] * T_c[!home]) * sum(soup[g]) / sum(cells[g, !home]),
      foreign_removed = 1 - sum(corr[g, !home]) / sum(cells[g, !home]),
      native_removed = 1 - sum(corr[g, home]) / sum(cells[g, home]),
      total_removed = 1 - sum(corr) / sum(cells))
  }
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE, digits = 3)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, file.path(OUT, "pbmc_soupx_markers.csv"), row.names = FALSE)
utils::write.csv(cbind(res$rho_diagnostics, rho_cell_median = stats::median(res$rho)),
                 file.path(OUT, "pbmc_soupx_rho.csv"), row.names = FALSE)
