# Simulate a case/control single-cell study with sample-specific soup.

Simulate a case/control single-cell study with sample-specific soup.

## Usage

``` r
simulate_ambient_study(
  n_per = 6,
  G = 1000,
  n_cells = 700,
  rho_ctrl = 0.05,
  rho_dis = 0.12,
  bcv = 0.25,
  n_de_rare = 40,
  fc_rare = 2,
  n_de_abund = 80,
  fc_abund = 3,
  frac = c(A = 0.55, B = 0.2, C = 0.15, R = 0.1),
  n_empty = 2000,
  n_shared = 0,
  seed = 1
)
```

## Arguments

- n_per:

  samples per condition.

- G:

  genes.

- n_cells:

  cells per sample.

- rho_ctrl, rho_dis:

  mean contamination in control / disease samples (disease tissue is
  often more fragile).

- bcv:

  biological coefficient of variation between samples.

- n_de_rare, fc_rare:

  true disease genes in the rare type and their fold.

- n_de_abund, fc_abund:

  disease genes in the abundant type A (all chosen among genes the rare
  type does not express).

- frac:

  cell-type proportions; must include types "A" and "R".

- n_empty:

  empty droplets per sample.

- n_shared:

  genes that are genuinely disease genes in BOTH A and R.

- seed:

  RNG seed.

## Value

list: X (genes x cells), label, sample, cond (per sample), empties (list
per sample), de_rare, de_abund (leak-only genes), shared, rho_true (per
sample), prof.
