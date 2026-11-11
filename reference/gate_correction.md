# GATE: which mass-conserving correction rule is closest to truth?

Scores "subtract" (v0), "redistribute" (SoupX-style) and "bayes" with
the ESTIMATED rho (not oracle), on easy and stress simulations. The
default in
[`decant()`](https://kroy3.github.io/Decant/reference/decant.md) is
whichever has the lowest mean `l1_error` here.

## Usage

``` r
gate_correction(rho_grid = c(0.05, 0.2), seeds = 1:3, verbose = TRUE)
```

## Arguments

- rho_grid:

  true mean contamination levels.

- seeds:

  replicate seeds.

- verbose:

  print the table.

## Value

invisibly, a data.frame with attribute `"winner"`.
