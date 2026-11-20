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
  correction = c("bayes_cell", "bayes", "subtract", "redistribute", "bayes_absent"),
  cap_test = c("entry", "cluster", "hybrid"),
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
  The default, "bayes_cell", has the lowest error on simulated ground
  truth
  ([`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md)).
  On real PBMC data (`validation/benchmark_pbmc_six_tools.R`) it removes
  88 to 95% of marker counts that ambient RNA could explain in cell
  types that do not express them, and leaves most of the counts that are
  too high in a cell to be ambient (likely doublets), which "bayes"
  removed.

- cap_test:

  for "bayes_cell": "entry", "cluster" or "hybrid", how the test that
  sets aside counts too high to be ambient is applied; see
  [`correct_counts()`](https://kroy3.github.io/Decant/reference/correct_counts.md).

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
