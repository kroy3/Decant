# Simulate a pooled multi-genotype experiment with allele-resolved SNP genes.

Simulate a pooled multi-genotype experiment with allele-resolved SNP
genes.

## Usage

``` r
simulate_allelic(
  n_cells = 1500,
  n_donors = 8,
  n_snp = 300,
  rho_mean = 0.2,
  rho_conc = 25,
  snp_lib = 800,
  seed = 1
)
```

## Arguments

- n_cells, n_donors, n_snp:

  dimensions.

- rho_mean, rho_conc:

  Beta parameters of per-cell contamination.

- snp_lib:

  mean SNP-covering reads per cell.

- seed:

  RNG seed.

## Value

list with per-cell own/other allele counts and true rho.
