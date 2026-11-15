# Chapter 4. Preparing Your Data

Decant needs two matrices from each sample: the **cells** and the
**empty droplets**. Both come from the raw (unfiltered) count matrix.

## 4.1 Read a 10x Genomics matrix

``` r

# Requires your own data
raw <- read_10x_counts("sample/outs/raw_feature_bc_matrix")
dim(raw)
```

- [`read_10x_counts()`](https://kroy3.github.io/Decant/reference/read_10x_counts.md)
  reads the three files Cell Ranger writes (`matrix.mtx.gz`,
  `features.tsv.gz`, `barcodes.tsv.gz`) and returns a sparse matrix with
  genes in rows and barcodes in columns.
- By default rows are named by Ensembl gene ID (`gene_column = 1`),
  which is unique. Use `gene_column = 2` for gene symbols, which some
  genes share.
- Only “Gene Expression” features are kept. Antibody or CRISPR features
  in the same run are dropped (`feature_type = NULL` keeps everything).
- A raw matrix usually has several hundred thousand columns; this is
  expected.

## 4.2 Separate cells from empty droplets

``` r

library(Decant)
library(Matrix)
sim <- simulate_experiment(seed = 1)
colnames(sim$empty) <- paste0("empty", seq_len(ncol(sim$empty)))
raw <- Matrix(cbind(sim$observed, sim$empty), sparse = TRUE)
dim(raw)
#> [1]  600 5200

d <- split_droplets(raw, cells = colnames(sim$observed))
d$summary
#>   n_cells n_empties n_dropped soup_umis
#>      1200      4000         0    100021
```

- The first four lines build a raw matrix from simulated data, so this
  example runs anywhere. With real data you would use the `raw` from
  4.1.
- `split_droplets(raw, cells = ...)` returns a list with `cells` (genes
  × cells), `empties` (genes × empty droplets) and `summary`.
- `cells` is the vector of cell barcodes. Use the barcodes of Cell
  Ranger’s filtered matrix, or the cells you kept after your own quality
  control.
- `summary` reports how many droplets went into each group and how many
  UMIs the soup estimate rests on (`soup_umis`). More than about 50,000
  soup UMIs gives a stable profile.

### How droplets are assigned

| Droplet | Rule | Used as |
|----|----|----|
| Listed in `cells` | always | cell |
| Total UMIs within `empty_umi_range` (default 1 to 100) | not a cell | empty droplet |
| Anything else | neither | dropped |

Droplets between the two groups are dropped on purpose: they may contain
a damaged cell, and treating them as soup would bias the profile.

### When you do not have a list of cells

``` r

d2 <- split_droplets(raw)
d2$summary
#>   n_cells n_empties n_dropped soup_umis
#>      1200      4000         0    100021
```

- Without `cells`, droplets with at least `cell_min_umi` (default 500)
  UMIs are treated as cells. This is a crude rule; a proper cell caller
  such as Cell Ranger, EmptyDrops or CellBender is recommended for real
  data.

## 4.3 Clusters

Decant estimates contamination cluster by cluster ([Chapter
2.5](https://kroy3.github.io/Decant/articles/02-Key-Concepts.md)), so it
needs one label per cell.

``` r

# Requires your own data
clusters <- as.character(Idents(seurat_object))
```

- Any clustering works: Seurat, Scanpy (imported into R), or cell-type
  annotations.
- Prefer a moderate number of reasonably large clusters (for example
  Seurat at resolution 0.3 to 0.8). Very small clusters have too few
  counts to tell a weakly expressed gene from an absent one, which
  pushes the estimate up.
- If you supply no clusters, Decant runs a quick k-means clustering
  (`k = 8`) itself. This is fine for exploration but less accurate.

## 4.4 Checklist

The matrix holds **raw UMI counts**, not normalised or scaled values.

Cells and empty droplets have the **same genes in the same order**.
([`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md)
guarantees this.)

The cluster vector has **one entry per cell**, in column order.

There are enough empty droplets: tens of thousands is typical.

------------------------------------------------------------------------

**Previous:** [Chapter 3, Quick
Start](https://kroy3.github.io/Decant/articles/03-Quick-Start.md) \|
**Next:** [Chapter 5, Running
Decant](https://kroy3.github.io/Decant/articles/05-Running-Decant.md)
