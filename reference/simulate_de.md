# Self-contained pseudobulk DE simulator with three gene classes: null : no difference anywhere true_de : differs in NATIVE expression between conditions (real signal) ambient_de: native identical, but SOUP differs between conditions (a trap)

Self-contained pseudobulk DE simulator with three gene classes: null :
no difference anywhere true_de : differs in NATIVE expression between
conditions (real signal) ambient_de: native identical, but SOUP differs
between conditions (a trap)

## Usage

``` r
simulate_de(
  n_genes = 400,
  n_true = 40,
  n_amb = 40,
  n_rep = 5,
  rho_bar = 0.2,
  fc = 2.5,
  amb_fc = 3,
  seed = 1
)
```

## Arguments

- n_genes, n_true, n_amb, n_rep:

  gene-class sizes and replicates per condition.

- rho_bar:

  mean contamination.

- fc, amb_fc:

  native and ambient fold changes.

- seed:

  RNG seed.
