# Empirical-Bayes soup profiles per sample.

Each sample's soup is shrunk, gene by gene, toward the cross-sample mean
by `w = tau2 / (tau2 + v)`, where `v` is the multinomial sampling
variance of that sample's estimate and `tau2` the estimated
between-sample variance of the true soups. Well-measured samples keep
their own soup; poorly measured ones borrow. (A fixed-rule pooling was
tried and rejected: it erased real sample-specific soup differences and
inflated false positives.)

## Usage

``` r
soup_eb(empties, design = NULL)
```

## Arguments

- empties:

  list of genes x droplets empty-droplet matrices, one per sample (same
  genes, same order).

- design:

  optional sample-level design matrix. Each sample's soup is then shrunk
  toward the design-fitted soup instead of the overall mean, so real
  soup differences between conditions are not shrunk away.

## Value

list(soup = genes x samples, var = genes x samples posterior variance,
umis = total soup UMIs per sample).
