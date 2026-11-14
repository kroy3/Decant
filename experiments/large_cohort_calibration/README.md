# Large-cohort calibration of contamination-aware DE

Li et al. (2022, Genome Biology 23:79) reported that negative-binomial DE
tests lose false-positive control in population-scale cohorts, because a few
samples with extreme counts shift a group mean while the NB dispersion
understates their variance. `ambient_de()` is an NB quasi-likelihood test, so
it was checked for the same failure, and it had it.

`run.R` sweeps 6, 15, 50 and 100 samples per group over three cases from
`simulate_ambient_cohort()`: clean counts, 1% outlier entries (median 10-fold),
and heavy-tailed (log-t, 3 df) biological noise. Contamination is about 5% in
controls and 12% in disease. Means of 3 seeds are in `results.csv`.

Methods:

- `ambient_de`: with Huber down-weighting of outlier samples (the default).
- `ambient_de_nonrobust`: `robust = FALSE`, the behaviour before this change.
- `naive_edgeR`: edgeR quasi-likelihood on raw pseudobulk, `robust = TRUE`.
- `naive_wilcoxon`: Wilcoxon rank-sum on TMM-normalised CPM of raw pseudobulk.
- `subtract_wilcoxon`: subtract the ambient estimate `ambient_de` uses, then
  Wilcoxon.

## Null false-positive rate (nominal 0.05)

| per group | case | ambient_de | non-robust | naive edgeR | naive Wilcoxon | subtract + Wilcoxon |
|---|---|---|---|---|---|---|
| 6 | clean | 0.051 | 0.049 | 0.131 | 0.117 | 0.044 |
| 6 | outliers | 0.068 | 0.124 | 0.129 | 0.101 | 0.041 |
| 6 | heavy tail | 0.058 | 0.086 | 0.118 | 0.087 | 0.039 |
| 15 | clean | 0.055 | 0.052 | 0.198 | 0.184 | 0.046 |
| 15 | outliers | 0.055 | 0.200 | 0.260 | 0.185 | 0.045 |
| 15 | heavy tail | 0.061 | 0.109 | 0.189 | 0.147 | 0.044 |
| 50 | clean | 0.056 | 0.050 | 0.402 | 0.393 | 0.048 |
| 50 | outliers | 0.050 | 0.338 | 0.430 | 0.333 | 0.054 |
| 50 | heavy tail | 0.049 | 0.155 | 0.303 | 0.261 | 0.050 |
| 100 | clean | 0.063 | 0.058 | 0.574 | 0.552 | 0.058 |
| 100 | outliers | 0.062 | 0.407 | 0.497 | 0.483 | 0.067 |
| 100 | heavy tail | 0.036 | 0.190 | 0.439 | 0.412 | 0.055 |

## Power (FDR 0.1, 40 true DE genes at 1.5-fold)

| per group | case | ambient_de | subtract + Wilcoxon |
|---|---|---|---|
| 6 | clean | 0.18 | 0.00 |
| 6 | outliers | 0.14 | 0.00 |
| 15 | clean | 0.75 | 0.50 |
| 15 | outliers | 0.73 | 0.59 |
| 15 | heavy tail | 0.38 | 0.18 |
| 50 | all | 0.98-0.99 | 0.99-1.00 |
| 100 | all | 0.98-1.00 | 1.00 |

No method that models or subtracts the ambient term called a leaked gene in
any run. The naive methods called 12 to 52 of 80.

## What this shows

1. Without down-weighting, `ambient_de` fails as Li et al. describe: up to
   41% of null genes at p < 0.05 with 100 per group.
2. A rank test is not a fix for contamination. Naive Wilcoxon is as inflated
   as naive edgeR whenever contamination differs between groups, and the
   inflation grows with sample size.
3. Subtracting the ambient estimate and running Wilcoxon is as good as
   `ambient_de` from 50 per group up. Below that `ambient_de` has more power.
   A rank test also cannot adjust for covariates; that case is not tested
   here.

## Limits

- With 100 per group the clean and outlier cases drift to 0.06 for both
  ambient-aware methods, and the false discovery proportion reaches 0.14 to
  0.17. The per-sample contamination estimate is biased upward by about one
  to two points, and that bias does not shrink as samples are added, so it
  becomes detectable. Expect this to worsen in larger cohorts.
- With 6 per group and outliers the rate is 0.068. Outliers cannot be
  identified with so few samples.
- At small sample sizes the false discovery proportion is high (0.23 to 0.64
  with 6 per group) on few calls.
- `huber_k = 2` was chosen on this simulator: 1.345 gave 0.07 on clean data,
  4 let outliers through.
- The simulator generates data from the additive ambient model that
  `ambient_de` fits. Nothing here is real data.
