# Legacy v0 per-cell rho estimator (superseded).

Median of observed/soup ratios over soup-rich genes. Because those genes
include ones the cell expresses natively, it overestimates rho 1.5-4x on
ground truth (see
[`gate_rho()`](https://kroy3.github.io/Decant/reference/gate_rho.md)).
Kept only so earlier benchmark numbers remain reproducible; use
[`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md)
instead.

## Usage

``` r
estimate_rho(observed, soup, low_expr_q = 0.5)
```

## Arguments

- observed:

  genes x cells.

- soup:

  length-G ambient profile.

- low_expr_q:

  genes below this soup-quantile are ignored as estimators.

## Value

per-cell rho.
