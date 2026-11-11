# Contamination-aware differential expression on pseudobulk counts.

Contamination-aware differential expression on pseudobulk counts.

## Usage

``` r
ambient_de(
  counts,
  empties = NULL,
  design,
  coef = ncol(design),
  soup = NULL,
  min_count = 10,
  naive = TRUE,
  maxit = 50,
  tol = 1e-08
)
```

## Arguments

- counts:

  genes x samples pseudobulk counts of ONE cell type.

- empties:

  list (one per sample, same order as `counts` columns) of genes x
  droplets empty-droplet matrices for the same genes; OR supply `soup`
  directly.

- design:

  sample-level design matrix (rows = samples).

- coef:

  column of `design` to test (index or name).

- soup:

  optional precomputed
  [`soup_eb()`](https://kroy3.github.io/Decant/reference/soup_eb.md)
  result.

- min_count:

  genes are tested only if their soup-removed pseudobulk passes edgeR's
  `filterByExpr` at this `min.count`.

- naive:

  also run standard edgeR QL on the raw pseudobulk and report which of
  its hits the contamination-aware test does not support.

- maxit, tol:

  Fisher-scoring controls.

## Value

data.frame per tested gene: logFC (native-expression log2 fold change),
F, PValue, FDR, ambient_frac (mean fraction of the gene's counts
attributed to ambient), and when `naive = TRUE` naive_PValue, naive_FDR
and `ambient_driven` (naive FDR \< 0.05 but contamination-aware FDR \>=
0.05). Attribute `"samples"`: per-sample ambient fraction and library
factors.
