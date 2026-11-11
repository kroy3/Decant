# Cluster-pooled, absent-gene rho estimator with per-cell shrinkage.

Cluster-pooled, absent-gene rho estimator with per-cell shrinkage.

## Usage

``` r
estimate_rho_cluster(
  observed,
  soup,
  clusters,
  alpha = 0.01,
  min_absent_soup = 0.02,
  max_iter = 100
)
```

## Arguments

- observed:

  genes x cells counts (base matrix or dgCMatrix).

- soup:

  length-G ambient profile (e.g. from
  [`ambient_global()`](https://kroy3.github.io/Decant/reference/ambient_global.md)).

- clusters:

  length-cells cluster labels. Prefer a moderate number of reasonably
  large clusters: very small clusters lose the power to reject
  weakly-expressed genes and bias rho upward.

- alpha:

  one-sided Poisson level for calling a gene "above soup" (i.e. natively
  expressed) in a cluster.

- min_absent_soup:

  a cluster whose absent set covers less than this fraction of soup mass
  is flagged `low_power`.

- max_iter:

  iteration cap for the absent-set fixed point.

## Value

numeric vector of per-cell rho, with attribute `"clusters"`: a
data.frame of per-cluster diagnostics (n_cells, rho, n_absent_genes,
absent_soup_frac, absent_counts, low_power, borrowed, rho_two_sided),
and attribute `"soup_mismatch"` (logical). Low-power clusters take the
evidence-weighted rho of the other clusters (`borrowed = TRUE`).
`soup_mismatch = TRUE` means the soup profile looks inconsistent with
the cells and rho is probably UNDERestimated – inspect the empties or
supply rho yourself.
