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
  robust = TRUE,
  joint_rho = TRUE,
  soup_design = TRUE,
  lfc = 0,
  max_ambient_frac = 0.75,
  huber_k = 2,
  robust_iter = 2,
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

- robust:

  down-weight outlier samples per gene (Huber weights on Pearson
  residuals, with the dispersion re-estimated under those weights).
  Guards against the false positives NB tests produce in large cohorts
  when a few samples carry extreme counts.

- joint_rho:

  refine each sample's ambient fraction on genes that look absent from
  the cell type across ALL samples pooled. A single sample cannot tell a
  weakly expressed gene from an absent one, which biases its ambient
  fraction upward; the pooled test can.

- soup_design:

  shrink each sample's soup toward the design-fitted soup (see
  [`soup_eb()`](https://kroy3.github.io/Decant/reference/soup_eb.md))
  instead of the overall mean.

- lfc:

  test against `|log2 fold change| <= lfc` instead of zero (as limma's
  TREAT). In cohorts of a hundred or more samples per group, residual
  biases well under 1% become statistically detectable; a small
  threshold such as `log2(1.1)` asks for effects that matter.

- max_ambient_frac:

  genes whose counts in this cell type are, on average, more than this
  fraction ambient are not tested (they are listed in attribute
  `"ambient_dominated"`). Their native signal is too small for the
  ambient estimate's error.

- huber_k:

  Huber tuning constant, in Pearson-residual units. The default 2 kept
  null genes calibrated with and without outliers in
  [`gate_ambient_de_cohort()`](https://kroy3.github.io/Decant/reference/gate_ambient_de_cohort.md);
  smaller values over-trim clean data and inflate false positives,
  larger ones let outliers through.

- robust_iter:

  number of re-weighting rounds.

- maxit, tol:

  Fisher-scoring controls.

## Value

data.frame per tested gene: logFC (native-expression log2 fold change),
F, PValue, FDR, ambient_frac (mean fraction of the gene's counts
attributed to ambient), and when `naive = TRUE` naive_PValue, naive_FDR
and `ambient_driven` (naive FDR \< 0.05 but contamination-aware FDR \>=
0.05). Attribute `"samples"`: per-sample ambient fraction and library
factors. Attribute `"ambient_dominated"`: genes skipped as mostly soup.
