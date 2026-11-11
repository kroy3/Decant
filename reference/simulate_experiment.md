# Simulate a contaminated droplet experiment with known ground truth.

Simulate a contaminated droplet experiment with known ground truth.

## Usage

``` r
simulate_experiment(
  n_genes = 600,
  n_types = 6,
  n_cells = 1200,
  soup_bias = 3,
  rho_mean = 0.2,
  rho_conc = 25,
  n_empty = 4000,
  lib_mean = 3000,
  seed = 1,
  shared_frac = 0,
  overdispersion = 0,
  type_prob = NULL,
  lib_sdlog = 0.4
)
```

## Arguments

- n_genes, n_types, n_cells:

  dimensions of the experiment.

- soup_bias:

  controls how non-uniform the soup is. 0 = every type lyses equally
  (soup == global average of cells, SoupX's assumption holds exactly);
  higher values = a few fragile types dominate the soup (the regime
  where a structured model should help, and where global models should
  fail).

- rho_mean:

  mean contamination fraction per cell (snRNA-seq tends to run higher
  than scRNA-seq; try 0.1 vs 0.3).

- rho_conc:

  Beta concentration of per-cell rho (higher = less spread).

- n_empty:

  number of empty droplets profiled (the soup observation).

- lib_mean:

  mean library size.

- seed:

  RNG seed.

- shared_frac:

  in \[0,1): fraction of every type's profile drawn from one shared
  "housekeeping" program. Higher = fewer genes that are cleanly absent
  from any type, which is exactly what makes rho hard to identify. 0
  (default) reproduces the original, easy simulator.

- overdispersion:

  cell-level gamma noise on the expression profile (variance of the
  multiplicative factor). 0 (default) = pure multinomial.

- type_prob:

  optional length-n_types relative abundances (e.g. to create rare
  types). NULL (default) = equal abundances.

- lib_sdlog:

  log-scale sd of library sizes.

## Value

list with \$observed, \$truth (clean own counts), \$empty, \$soup_true,
\$lysis_true, \$rho_true, \$labels, \$profiles. Matrices are genes x
cells.

## Details

The three stress knobs (`shared_frac`, `overdispersion`, `type_prob`)
consume no random draws when left at their defaults, so results for a
given seed are identical to earlier versions.
