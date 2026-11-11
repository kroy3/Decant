# Derive a program basis (cluster mean profiles) from the observed cells.

Derive a program basis (cluster mean profiles) from the observed cells.

## Usage

``` r
cluster_basis(observed, k = 6, n_hvg = 200, seed = 1, clusters = NULL)
```

## Arguments

- observed:

  genes x cells counts.

- k:

  number of k-means clusters (ignored if `clusters` is given).

- n_hvg:

  top-variance genes used for clustering.

- seed:

  RNG seed.

- clusters:

  optional cluster labels (e.g. from Seurat / SCE); when given, no
  clustering is run.

## Value

genes x K matrix, columns are probability vectors.
