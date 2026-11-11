# Default clustering used when the caller supplies none.

k-means on log-normalised top-variance genes. Sparse-safe (only an n_hvg
x cells block is densified). For real analyses pass your own Seurat /
SingleCellExperiment clusters instead.

## Usage

``` r
quick_labels(observed, k = 6, n_hvg = 200, seed = 1)
```

## Arguments

- observed:

  genes x cells counts.

- k:

  number of clusters.

- n_hvg:

  number of top-variance genes used.

- seed:

  RNG seed for k-means.

## Value

integer cluster labels.
