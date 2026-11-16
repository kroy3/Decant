# Chapter 5. Running Decant

[`decant()`](https://kroy3.github.io/Decant/reference/decant.md) is the
main function. It estimates the soup, estimates each cell’s
contamination, and removes it.

## 5.1 The basic call

``` r

library(Decant)
sim <- simulate_experiment(seed = 1)
res <- decant(sim$observed, sim$empty, clusters = sim$labels)
res
#> <decant> 600 genes x 1200 cells
#>   modules : cluster_rho, correction_bayes_cell
#>   rho     : median 0.199 (10%-90%: 0.117-0.305)
```

- First argument: the cell counts, genes × cells.
- Second argument: the empty droplets, genes × droplets.
- `clusters`: one label per cell.
- The printout lists the steps that ran and summarises rho.

## 5.2 Arguments

| Argument | Default | Purpose |
|----|----|----|
| `obs_total` | (required) | Cell counts, genes × cells, dense or sparse |
| `empties` | (required) | Empty droplets; a matrix, or a list of matrices for several samples ([Chapter 8](https://kroy3.github.io/Decant/articles/08-Multi-Sample-Experiments.md)) |
| `clusters` | `NULL` | Cluster labels; if missing, k-means with `k` clusters is run |
| `k` | `8` | Number of k-means clusters when `clusters` is missing |
| `rho` | `NULL` | Your own contamination fraction, one value or one per cell; skips estimation |
| `correction` | `"bayes_cell"` | Correction rule ([Chapter 7](https://kroy3.github.io/Decant/articles/07-Estimation-and-Correction.md)) |
| `sample_of` | `NULL` | Sample of each cell, when `empties` is a list |
| `obs_unspliced`, `obs_spliced`, `empties_unspliced`, `empties_spliced` | `NULL` | Velocity layers ([Chapter 10](https://kroy3.github.io/Decant/articles/10-Splice-Aware-and-Genotype-Modules.md)) |
| `allelic` | `NULL` | Genotype data for contamination estimation ([Chapter 10](https://kroy3.github.io/Decant/articles/10-Splice-Aware-and-Genotype-Modules.md)) |
| `verbose` | `FALSE` | Print a short summary while running |

## 5.3 The result object

``` r

names(res)
#> [1] "rho"                "rho_diagnostics"    "soup_mismatch"     
#> [4] "rho_unidentifiable" "ambient"            "clusters"          
#> [7] "corrected"          "modules"           
```

| Element | Content |
|----|----|
| `corrected` | Corrected counts, same genes, cells and class as the input |
| `rho` | Estimated contamination fraction of each cell |
| `rho_diagnostics` | One row per cluster: its rho and how well it was determined |
| `soup_mismatch` | `TRUE` if much of the soup is in genes no cell contains ([Chapter 11](https://kroy3.github.io/Decant/articles/11-Diagnostics-and-Troubleshooting.md)) |
| `rho_unidentifiable` | `TRUE` if contamination could not be estimated from the data; supply `rho` ([Chapter 11.4](https://kroy3.github.io/Decant/articles/11-Diagnostics-and-Troubleshooting.md)) |
| `ambient` | The soup profile, genes × samples |
| `clusters` | The cluster labels used |
| `modules` | The steps that ran |

## 5.4 Per-cluster diagnostics

``` r

res$rho_diagnostics[, c("cluster", "n_cells", "rho", "n_absent_genes", "absent_soup_frac", "low_power")]
#>   cluster n_cells       rho n_absent_genes absent_soup_frac low_power
#> 1       1     222 0.2008547            188       0.82722628     FALSE
#> 2       2     202 0.2032755            214       0.84945162     FALSE
#> 3       3     201 0.2130331            226       0.86184901     FALSE
#> 4       4     179 0.2195561            239       0.85055138     FALSE
#> 5       5     196 0.2099970            185       0.73371592     FALSE
#> 6       6     200 0.1895710             97       0.05638816     FALSE
```

- `rho`: the cluster’s contamination estimate.
- `n_absent_genes`: how many genes the cluster was found not to express.
  These genes carry the evidence for rho.
- `absent_soup_frac`: the share of the soup those genes cover. The
  larger it is, the better rho is determined. Cluster 6 relies on only
  5.6% of the soup and is the least certain.
- `low_power`: `TRUE` when the evidence is too thin, or when the
  cluster’s estimate was implausible (above 0.5). Such a cluster borrows
  its rho from the other clusters, and Decant prints a warning.
- Two further columns, not shown here: `dispersion` (how widely absent
  genes scatter around rho, from gene to gene) and `implausible`.

## 5.5 Sparse input

``` r

library(Matrix)
sp <- Matrix(sim$observed, sparse = TRUE)
res_sp <- decant(sp, Matrix(sim$empty, sparse = TRUE), clusters = sim$labels)
class(res_sp$corrected)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```

- Real data are sparse. Decant works directly on sparse `dgCMatrix`
  objects and returns the same class, without ever building a dense
  genes × cells matrix.

## 5.6 Supplying your own contamination fraction

``` r

res_fixed <- decant(sim$observed, sim$empty, rho = 0.2, clusters = sim$labels)
res_fixed$modules
#> [1] "user_rho"              "correction_bayes_cell"
```

- `rho = 0.2` applies 20% contamination to every cell and skips the
  estimation step. You can also pass one value per cell.
- Use this when you have an independent estimate, for example from
  genotypes or from a species-mixing control. The `bayes_cell` rule
  still needs `clusters`.

## 5.7 Without clusters

``` r

res_k <- decant(sim$observed, sim$empty)
res_k$modules
#> [1] "kmeans_clusters"       "cluster_rho"           "correction_bayes_cell"
```

- Without `clusters`, Decant clusters the cells itself with k-means
  (`kmeans_clusters` in `modules`). Supply your own clusters for real
  analyses.

------------------------------------------------------------------------

**Previous:** [Chapter 4, Preparing Your
Data](https://kroy3.github.io/Decant/articles/04-Preparing-Your-Data.md)
\| **Next:** [Chapter 6, Seurat
Integration](https://kroy3.github.io/Decant/articles/06-Seurat-Integration.md)
