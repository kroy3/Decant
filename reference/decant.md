# Run Decant.

Run Decant.

## Usage

``` r
decant(
  obs_total,
  empties,
  sample_of = NULL,
  clusters = NULL,
  k = 8,
  rho = NULL,
  correction = c("bayes", "subtract", "redistribute"),
  obs_unspliced = NULL,
  obs_spliced = NULL,
  empties_unspliced = NULL,
  empties_spliced = NULL,
  allelic = NULL,
  basis = NULL,
  verbose = FALSE
)
```

## Arguments

- obs_total:

  genes x cells total counts (base matrix or dgCMatrix).

- empties:

  EITHER a genes x droplets matrix (single sample) OR a list of such
  matrices, one per sample, which triggers hierarchical pooling (Gap 3).
  See
  [`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md)
  to obtain these from a raw 10x matrix.

- sample_of:

  length-cells sample index (required if empties is a list).

- clusters:

  optional length-cells cluster labels (e.g. from Seurat/SCE). Used by
  the rho estimator. If NULL, k-means with `k` clusters is run.

- k:

  number of clusters when `clusters` is NULL.

- rho:

  optional user-supplied per-cell (or scalar) contamination fraction;
  skips estimation entirely.

- correction:

  correction rule, see
  [`correct_counts()`](https://kroy3.github.io/Decant/reference/correct_counts.md).
  The default is "bayes". On simulated data the three rules are within a
  few percent of each other
  ([`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md));
  on real PBMC data (`validation/pbmc_soupx.R`) "bayes" removed 83 to
  94% of marker counts in cell types that do not express them, against
  19 to 57% for "subtract", at the same loss in the cell types that do.

- obs_unspliced, obs_spliced:

  genes x cells layers; if both given with empty layers, correction is
  splice-aware (Gap 1).

- empties_unspliced, empties_spliced:

  layer-resolved empties (matrix or list per sample), required for
  splice-aware mode.

- allelic:

  optional list(own, other, donor, abund) for genotype-aware rho (Gap
  5); overrides count-based rho when present.

- basis:

  optional genes x K program basis; if given, returns a lysis diagnostic
  (Gap 4).

- verbose:

  print a short run summary.

## Value

an object of class `decant`: a list with \$corrected (same class as the
input; + \$corrected_unspliced/\$corrected_spliced in splice mode),
\$rho, \$rho_diagnostics (per-cluster, from
[`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md)),
\$ambient (genes x samples), \$clusters, \$lysis (if basis), \$modules.
