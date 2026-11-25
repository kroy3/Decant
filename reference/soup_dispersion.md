# Per-gene negative-binomial dispersion of the empty droplets, from moments.

For each gene, the method-of-moments estimate
`alpha = max(0, (var - mean) / mean^2)` of per-droplet counts, then
shrunk toward a trend in log mean (a running median over genes ordered
by mean, weighted by detection), so that rarely detected genes borrow
strength. This ignores droplet-to-droplet differences in library size,
which inflate `alpha`; in empty droplets (a narrow UMI range) the effect
is small.

## Usage

``` r
soup_dispersion(soup, min_detected = 20)
```

## Arguments

- soup:

  a
  [`soup_summary()`](https://kroy3.github.io/Decant/reference/soup_summary.md)
  (or an empty-droplet matrix).

- min_detected:

  genes detected in fewer droplets take the trend value.

## Value

numeric vector of per-gene dispersions (NA where the gene has no
counts), with attribute `"raw"` (unshrunk estimates).
