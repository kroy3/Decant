# Simulate a large case/control cohort at the pseudobulk level.

The same ambient mechanism as
[`simulate_ambient_study()`](https://kroy3.github.io/Decant/reference/simulate_ambient_study.md),
but drawing the rare type's pseudobulk directly so that hundreds of
samples are cheap, and with optional departures from the negative
binomial: per-gene outlier samples (`out_p`, `out_fold`) or heavy-tailed
biological noise (`heavy_df`). These are the conditions under which Li
et al. (2022, Genome Biology 23:79) found NB-based DE to lose
false-positive control in population-scale data.

## Usage

``` r
simulate_ambient_cohort(
  n_per = 50,
  G = 2000,
  rho_ctrl = 0.05,
  rho_dis = 0.12,
  bcv = 0.25,
  n_de_rare = 40,
  fc_rare = 1.5,
  n_de_abund = 80,
  fc_abund = 3,
  frac = c(A = 0.55, B = 0.2, C = 0.15, R = 0.1),
  out_p = 0,
  out_fold = 10,
  heavy_df = Inf,
  n_batch = 0,
  batch_fc = 1.5,
  batch_imbalance = 0.7,
  depth = 175000,
  soup_umis = 50000,
  seed = 1
)
```

## Arguments

- n_per:

  samples per condition.

- G:

  genes.

- rho_ctrl, rho_dis:

  mean contamination in control / disease samples.

- bcv:

  biological coefficient of variation between samples.

- n_de_rare, fc_rare:

  true disease genes in the rare type and their fold.

- n_de_abund, fc_abund:

  disease genes in the abundant type A.

- frac:

  cell-type proportions; must include types "A" and "R".

- out_p:

  probability that a gene's expression in a sample (in any cell type) is
  an outlier.

- out_fold:

  median fold of an outlier.

- heavy_df:

  if finite, biological noise is log-t with this many degrees of freedom
  instead of gamma.

- n_batch, batch_fc:

  number of rare-type genes shifted by a batch, and their fold. The
  genes are not disease genes.

- batch_imbalance:

  probability that a disease sample is in batch 1 (a control sample is
  with one minus this), so batch confounds condition.

- depth:

  mean pseudobulk depth of the rare type per sample.

- soup_umis:

  empty-droplet UMIs per sample.

- seed:

  RNG seed.

## Value

list: Y (genes x samples pseudobulk of the rare type), empties (list per
sample), cond, batch, de_rare, de_abund, batch_genes, rho_true.
