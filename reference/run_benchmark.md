# Global vs structured soup benchmark.

Global vs structured soup benchmark.

## Usage

``` r
run_benchmark(
  soup_bias_grid = c(0, 1, 2, 4),
  reps = 3,
  rho_mean = 0.2,
  n_genes = 600,
  n_types = 6,
  n_cells = 1200,
  n_empty = 4000
)
```

## Arguments

- soup_bias_grid:

  vector of soup-bias values to test (0 = uniform soup).

- reps:

  replicate simulations per setting (different seeds).

- rho_mean, n_genes, n_types, n_cells, n_empty:

  passed to
  [`simulate_experiment()`](https://kroy3.github.io/Decant/reference/simulate_experiment.md).

## Value

data.frame, one row per (method, soup_bias, rep).
