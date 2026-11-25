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
  correction = c("bayes_absent", "bayes_cell", "bayes", "subtract", "redistribute"),
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
  Each may instead be a
  [`soup_summary()`](https://kroy3.github.io/Decant/reference/soup_summary.md)
  (O(genes) statistics of the empty droplets), e.g. from
  [`read_10x_streaming()`](https://kroy3.github.io/Decant/reference/read_10x_streaming.md),
  which never loads the empties. See
  [`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md)
  to obtain matrices from a raw 10x matrix.

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
  The default, "bayes_absent", decides per gene and cluster whether a
  gene is expressed at all and removes all counts of genes that are not;
  where that call is unreliable (overlapping components) it falls back
  to "bayes_cell". It has the lowest error on simulated ground truth
  ([`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md)).
  On two unused 10x human/mouse mixtures it removed 0.891 of ambient
  signal (bayes_cell 0.531, DecontX 0.655) and kept 0.9918 of native
  signal (bayes_cell 0.9945); on the hardest mixture seen (hgmm 10k) it
  kept 2.5 points less than bayes_cell. On PBMC 4k its call never fired
  and the output matched bayes_cell. It missed its pre-registered retain
  criterion by 0.02 points against DecontX and was made the default by
  maintainer decision (`validation/BENCHMARK.md` 1h). Use "bayes_cell"
  to keep the most native signal.

- cap_test:

  for "bayes_cell" ("bayes_absent" always uses "hybrid"): "entry",
  "cluster" or "hybrid", how the test that sets aside counts too high to
  be ambient is applied; see
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
\$identifiability (per cluster: cosine and Gram condition number kappa
between the cluster and the soup, Jensen-Shannon divergence, Fisher
information per UMI and its effective number of genes, rho's standard
error and 95% CI width, and whether it is weak, why, and what was done;
NULL when rho was supplied), \$ambient (genes x samples), \$clusters,
\$lysis (if basis), \$modules.
