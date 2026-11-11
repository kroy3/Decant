# GATE: is the cluster/absent-gene rho estimator calibrated, and better than the legacy marker-ratio estimator?

Pass requires, in every (setting, rho) cell: lower RMSE than legacy, and
\|mean bias\| \< max(25% of true rho, 0.015). Negative control: at rho =
0 the estimator must not hallucinate contamination (mean estimate \<
0.01). The soup-mismatch tripwire must not fire on any of these clean
simulations.

## Usage

``` r
gate_rho(
  rho_grid = c(0, 0.02, 0.05, 0.1, 0.2, 0.3),
  seeds = 1:3,
  verbose = TRUE
)
```

## Arguments

- rho_grid:

  true mean contamination levels to test.

- seeds:

  replicate seeds.

- verbose:

  print the table.

## Value

invisibly, a data.frame of results with attribute `"pass"`.
