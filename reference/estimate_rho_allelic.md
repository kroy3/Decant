# Allelic rho estimate. Fraction of mismatched-allele reads = rho \* P(soup read is from another donor). With pool abundances a_d, for a cell of donor d that probability is (1 - a_d). Invert to recover rho.

Allelic rho estimate. Fraction of mismatched-allele reads = rho \*
P(soup read is from another donor). With pool abundances a_d, for a cell
of donor d that probability is (1 - a_d). Invert to recover rho.

## Usage

``` r
estimate_rho_allelic(allelic)
```

## Arguments

- allelic:

  list with own, other, donor, abund (see
  [`simulate_allelic()`](https://kroy3.github.io/Decant/reference/simulate_allelic.md)).

## Value

per-cell rho.
