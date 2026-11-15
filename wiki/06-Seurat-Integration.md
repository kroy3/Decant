# Chapter 6. Seurat Integration

`RunDecant()` corrects a Seurat object in one step. It finds the empty
droplets, matches genes and barcodes to the object, runs `decant()`, and
stores the result as a new assay.

## 6.1 A complete example

```r
library(Decant)
library(Seurat)
library(Matrix)

sim <- simulate_experiment(seed = 1)
colnames(sim$empty) <- paste0("empty", seq_len(ncol(sim$empty)))
raw <- Matrix(cbind(sim$observed, sim$empty), sparse = TRUE)

obj <- CreateSeuratObject(counts = Matrix(sim$observed, sparse = TRUE))
Idents(obj) <- factor(sim$labels)

obj <- RunDecant(obj, raw = raw)
obj
#> An object of class Seurat
#> 1200 features across 1200 samples within 2 assays
#> Active assay: RNA (600 features, 0 variable features)
#>  1 layer present: counts
#>  1 other assay present: Decant
```

- The first block builds a raw matrix and a Seurat object from simulated
  data. In practice you already have both.
- `Idents(obj)` holds the clusters. `RunDecant()` uses them by default.
- `RunDecant(obj, raw = raw)` reads the empty droplets from `raw`, corrects
  the counts of the cells in `obj`, and adds an assay called `Decant`. The
  original `RNA` assay is not changed.

## 6.2 With a Cell Ranger folder

```r
# Requires your own data
obj <- RunDecant(obj, raw = "sample/outs/raw_feature_bc_matrix")
```

- `raw` can be a path. `RunDecant()` reads gene symbols and converts them
  the same way Seurat does (underscores become dashes), so genes match.
- Barcodes are matched even if Seurat added a prefix (from
  `merge(add.cell.ids = ...)`) or a suffix (from integration).
- Any cell that cannot be matched stops the run with examples, rather than
  being dropped silently.

## 6.3 What is added to the object

```r
summary(obj$decant_rho)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02389 0.14950 0.19933 0.20574 0.25223 0.55303 
names(Misc(obj, "decant"))
#> [1] "rho_diagnostics" "soup_mismatch"   "ambient"         "modules"        
#> [5] "source_assay"   
```

| Location | Content |
|---|---|
| assay `Decant` | Corrected counts |
| `obj$decant_rho` | Contamination fraction of each cell, in `meta.data` |
| `Misc(obj, "decant")` | Per-cluster diagnostics, soup profile, steps that ran |

`decant_rho` can be plotted like any metadata column, for example
`FeaturePlot(obj, "decant_rho")` or `VlnPlot(obj, "decant_rho")`.

## 6.4 Continue the analysis on corrected counts

```r
obj <- NormalizeData(obj, assay = "Decant")
DefaultAssay(obj) <- "Decant"
```

- The `Decant` assay holds counts, so normalise it before use.
- Setting it as the default assay makes `FindVariableFeatures()`,
  `ScaleData()`, `RunPCA()` and `FindMarkers()` use the corrected counts.

## 6.5 Several samples in one object

```r
# Requires your own data
obj <- RunDecant(obj,
                 raw = list(donor1 = "donor1/outs/raw_feature_bc_matrix",
                            donor2 = "donor2/outs/raw_feature_bc_matrix"),
                 sample_col = "donor")
```

- `raw` is a named list with one raw matrix per sample.
- `sample_col` names the `meta.data` column that says which sample each
  cell came from. Its values must match the names of `raw`.
- Each sample gets its own soup, shared information between samples is
  used as described in [Chapter 8](08-Multi-Sample-Experiments), and Seurat
  v5 split layers are joined automatically.

## 6.6 Options

| Argument | Default | Purpose |
|---|---|---|
| `raw` | `NULL` | Raw matrix, path, or named list per sample |
| `empties` | `NULL` | Empty droplets directly, instead of `raw` |
| `assay` | `"RNA"` | Assay holding raw counts |
| `clusters` | `NULL` | A `meta.data` column name or a vector; `NULL` uses `Idents()` |
| `sample_col` | `NULL` | Sample column, required for several samples |
| `new_assay` | `"Decant"` | Name of the new assay |
| `set_default` | `FALSE` | Make the new assay the default |
| `empty_umi_range` | `c(1, 100)` | UMI range that defines an empty droplet |
| `...` | | Passed to `decant()`, for example `correction` or `rho` |

`RunDecant()` refuses normalised or SCTransformed data: Decant works on raw
UMI counts only.

---

**Previous:** [Chapter 5, Running Decant](05-Running-Decant) | **Next:** [Chapter 7, Estimation and Correction](07-Estimation-and-Correction)
