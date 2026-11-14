# GATE: does contamination-aware DE stay calibrated in large cohorts whose counts are not negative binomial?

Simulated case/control cohorts
([`simulate_ambient_cohort()`](https://kroy3.github.io/Decant/reference/simulate_ambient_cohort.md))
with `n_per` samples per group, in four cases: clean NB-like counts,
per-gene outlier samples, heavy-tailed biological noise, and a batch
that confounds condition and shifts 200 genes (the model-based methods
adjust for it; a rank test cannot). Compared against ambient_de without
outlier down-weighting, naive edgeR (robust QL), a naive Wilcoxon
rank-sum test, and a Wilcoxon test on ambient-subtracted counts. Pass
requires of ambient_de, averaged over seeds, in every case: at most 1
leaked gene called, null false-positive rate \<= 0.075, false discovery
proportion \<= 0.15, and power no more than 0.1 below the best
comparator that itself passes the first three criteria.

## Usage

``` r
gate_ambient_de_cohort(seeds = 1:4, n_per = 50, verbose = TRUE)
```

## Arguments

- seeds:

  replicate seeds.

- n_per:

  samples per condition.

- verbose:

  print the table.

## Value

invisibly, a data.frame with attribute `"pass"`.
