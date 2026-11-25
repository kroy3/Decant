# Summary statistics of empty droplets (a streamed soup).

An O(genes) summary that
[`decant()`](https://kroy3.github.io/Decant/reference/decant.md),
[`ambient_global()`](https://kroy3.github.io/Decant/reference/ambient_global.md)
and
[`pool_soup()`](https://kroy3.github.io/Decant/reference/pool_soup.md)
accept in place of an empty-droplet matrix. Build it from a matrix with
`soup_summary()` or, without ever loading the empties, with
[`read_10x_streaming()`](https://kroy3.github.io/Decant/reference/read_10x_streaming.md).

## Usage

``` r
soup_summary(empty)
```

## Arguments

- empty:

  genes x droplets matrix (dense or sparse).

## Value

an object of class `decant_soup`: a list with `gene_sums`, `gene_sumsq`
(sum of squared counts), `gene_detected` (droplets with a count),
`total_umis`, `n_droplets` and `genes` (row names).
