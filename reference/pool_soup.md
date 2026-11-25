# Empirical-Bayes-style shrinkage of per-sample soup profiles.

Empirical-Bayes-style shrinkage of per-sample soup profiles.

## Usage

``` r
pool_soup(empties_total)
```

## Arguments

- empties_total:

  list of genes x droplets matrices (or
  [`soup_summary()`](https://kroy3.github.io/Decant/reference/soup_summary.md)
  objects), one per sample.

## Value

list with \$pooled (genes x samples) and \$independent (genes x
samples).
