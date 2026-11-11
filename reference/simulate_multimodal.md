# Multi-sample, splice-resolved contaminated experiment with ground truth.

Multi-sample, splice-resolved contaminated experiment with ground truth.

## Usage

``` r
simulate_multimodal(
  n_genes = 600,
  n_types = 6,
  n_cells = 1500,
  n_samples = 6,
  empty_per_sample = NULL,
  soup_bias = 3,
  rho_mean = 0.2,
  rho_conc = 25,
  between_sample_var = 0.15,
  splice_distinct = 0.6,
  lib_mean = 3000,
  seed = 1
)
```

## Arguments

- n_genes, n_types, n_cells:

  dimensions.

- n_samples:

  number of samples; each gets its own soup that shares a global
  structure but deviates (between_sample_var), and its own empty-droplet
  pool.

- empty_per_sample:

  vector (length n_samples) of empties per sample; set some low to
  stress hierarchical pooling.

- soup_bias, rho_mean, rho_conc, between_sample_var, lib_mean, seed:

  simulation knobs.

- splice_distinct:

  in 0-1: how differently ambient vs native split into
  unspliced/spliced. 0 = identical (negative control for the splice
  module).

## Value

list of layers (genes x cells) for spliced/unspliced observed and truth,
per-sample empties (with splice layers), sample/label vectors, and the
true soup profiles.
