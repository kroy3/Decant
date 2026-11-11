# Structured soup profile (lysis diagnostic).

Structured soup profile (lysis diagnostic).

## Usage

``` r
ambient_structured(empty, basis)
```

## Arguments

- empty:

  genes x droplets background matrix.

- basis:

  genes x K matrix of cell-state program profiles (columns need not be
  normalised; they are normalised internally). Supply cluster means from
  your filtered cells – see
  [`cluster_basis()`](https://kroy3.github.io/Decant/reference/cluster_basis.md).

## Value

list with \$soup (length-G probability vector) and \$lysis (length-K
non-negative weights summing to 1; the interpretable output telling you
which populations are contaminating the rest).
