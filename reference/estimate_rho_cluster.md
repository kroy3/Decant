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
  max_iter = 100,
  min_expected = 5,
  floor_rho = 0.005
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

- min_expected:

  a gene can join the absent set only if its expected ambient count in
  the cluster, at the current rho, is at least this. Below it the
  Poisson test cannot tell a weakly expressed gene from an absent one,
  and on real data (tens of thousands of low-count genes) letting those
  genes in drove rho to several times its true value.

- floor_rho:

  the smallest contamination considered possible. Genes whose total
  count in the cells is significantly below this fraction of their soup
  expectation are treated as soup-only and dropped from the soup.

## Value

numeric vector of per-cell rho, with attribute `"clusters"`: a
data.frame of per-cluster diagnostics (n_cells, rho, n_absent_genes,
absent_soup_frac, absent_counts, low_power, borrowed, rho_two_sided),
and attributes `"soup_mismatch"` (logical), `"soup_only"` (logical per
gene) and `"soup_only_frac"` (their share of the soup). Low-power
clusters take the evidence-weighted rho of the other clusters
(`borrowed = TRUE`). `soup_mismatch = TRUE` means more than 5% of the
soup was in genes no cell contains; those genes are dropped, but the
empties deserve a look.
