# Large-cohort calibration of contamination-aware DE

Li et al. (2022, Genome Biology 23:79) reported that negative-binomial DE
tests lose false-positive control in population-scale cohorts, because a few
samples with extreme counts shift a group mean while the NB dispersion
understates their variance. `ambient_de()` is an NB quasi-likelihood test, so
it was checked for the same failure, and it had it.

`run.R` sweeps 6, 15, 50, 100 and 200 samples per group over four cases from
`simulate_ambient_cohort()`: clean counts, 1% outlier entries (median 10-fold),
heavy-tailed (log-t, 3 df) biological noise, and a batch that is confounded
with condition and shifts 200 genes. Contamination is about 5% in controls and
12% in disease. Means of 3 seeds are in `results.csv`.

Methods:

- `ambient_de`: current defaults.
- `ambient_de_lfc1.1`: the same, tested against a 10% minimum fold change.
- `ambient_de_nonrobust`: `robust = FALSE`.
- `naive_edgeR`: edgeR quasi-likelihood on raw pseudobulk, `robust = TRUE`.
- `naive_wilcoxon`: Wilcoxon rank-sum on TMM-normalised CPM of raw pseudobulk.
- `subtract_wilcoxon`: subtract the ambient estimate `ambient_de` uses, then
  Wilcoxon.

In the batch case the model-based methods have batch in the design; the rank
tests cannot.

## Null false-positive rate (nominal 0.05)

| per group | case | ambient_de | non-robust | naive edgeR | naive Wilcoxon | subtract + Wilcoxon |
|---|---|---|---|---|---|---|
| 6 | clean | 0.050 | 0.046 | 0.137 | 0.117 | 0.036 |
| 6 | outliers | 0.073 | 0.133 | 0.123 | 0.095 | 0.047 |
| 6 | heavy tail | 0.055 | 0.074 | 0.106 | 0.078 | 0.039 |
| 6 | batch | 0.050 | 0.048 | 0.144 | 0.120 | 0.040 |
| 15 | clean | 0.049 | 0.044 | 0.201 | 0.185 | 0.041 |
| 15 | outliers | 0.045 | 0.201 | 0.248 | 0.160 | 0.042 |
| 15 | heavy tail | 0.055 | 0.101 | 0.185 | 0.142 | 0.052 |
| 15 | batch | 0.060 | 0.058 | 0.225 | 0.226 | 0.066 |
| 50 | clean | 0.056 | 0.052 | 0.382 | 0.369 | 0.054 |
| 50 | outliers | 0.053 | 0.335 | 0.444 | 0.366 | 0.053 |
| 50 | heavy tail | 0.044 | 0.146 | 0.303 | 0.267 | 0.047 |
| 50 | batch | 0.050 | 0.044 | 0.451 | 0.462 | 0.141 |
| 100 | clean | 0.054 | 0.049 | 0.574 | 0.550 | 0.050 |
| 100 | outliers | 0.047 | 0.390 | 0.515 | 0.533 | 0.047 |
| 100 | heavy tail | 0.034 | 0.187 | 0.414 | 0.386 | 0.050 |
| 100 | batch | 0.058 | 0.053 | 0.612 | 0.614 | 0.212 |
| 200 | clean | 0.064 | 0.057 | 0.759 | 0.738 | 0.053 |
| 200 | outliers | 0.047 | 0.418 | 0.552 | 0.697 | 0.055 |
| 200 | heavy tail | 0.046 | 0.272 | 0.566 | 0.587 | 0.053 |
| 200 | batch | 0.056 | 0.050 | 0.764 | 0.741 | 0.213 |

## Power and false discovery proportion (FDR 0.1, 40 true DE genes at 1.5-fold)

| per group | case | ambient_de power | FDP | subtract + Wilcoxon power | FDP | lfc 1.1 power | FDP |
|---|---|---|---|---|---|---|---|
| 6 | clean | 0.18 | 0.11 | 0.00 | 0.00 | 0.11 | 0.05 |
| 6 | outliers | 0.15 | 0.27 | 0.00 | 0.00 | 0.07 | 0.33 |
| 6 | heavy tail | 0.05 | 0.41 | 0.00 | 0.00 | 0.04 | 0.40 |
| 15 | clean | 0.86 | 0.07 | 0.61 | 0.05 | 0.59 | 0.04 |
| 15 | outliers | 0.65 | 0.07 | 0.47 | 0.02 | 0.38 | 0.00 |
| 15 | heavy tail | 0.44 | 0.27 | 0.18 | 0.05 | 0.16 | 0.42 |
| 15 | batch | 0.73 | 0.11 | 0.60 | 0.09 | 0.59 | 0.02 |
| 50 | clean | 1.00 | 0.16 | 1.00 | 0.11 | 1.00 | 0.00 |
| 50 | batch | 0.99 | 0.12 | 1.00 | 0.64 | 0.99 | 0.00 |
| 100 | clean | 1.00 | 0.16 | 1.00 | 0.12 | 1.00 | 0.00 |
| 100 | batch | 1.00 | 0.17 | 1.00 | 0.84 | 1.00 | 0.00 |
| 200 | clean | 1.00 | 0.14 | 1.00 | 0.10 | 1.00 | 0.00 |
| 200 | heavy tail | 1.00 | 0.16 | 1.00 | 0.09 | 1.00 | 0.00 |
| 200 | batch | 0.99 | 0.10 | 1.00 | 0.85 | 0.99 | 0.00 |

`ambient_de` called no leaked gene in any run. The naive methods called 12 to
68 of 80. Subtract + Wilcoxon called up to 4 in the batch case.

## What this shows

1. Without down-weighting, `ambient_de` fails as Li et al. describe: up to
   42% of null genes at p < 0.05 with 200 per group.
2. A rank test is not a fix for contamination. Naive Wilcoxon is as inflated
   as naive edgeR whenever contamination differs between groups, and the
   inflation grows with sample size (0.74 at 200 per group).
3. In a two-group design with 50 or more per group, subtracting the ambient
   estimate and running Wilcoxon is as good as `ambient_de`, and its false
   discovery proportion is lower (about 0.10 against 0.14 to 0.16).
4. `ambient_de` is the better choice below about 50 per group (power 0.86
   against 0.61 at 15) and whenever a covariate has to be adjusted for: with
   a confounding batch the rank test's false discovery proportion is 0.64 to
   0.85.
5. A 10% minimum fold change removes essentially all false discoveries from
   50 per group up, at no cost in power for 1.5-fold effects.

## Limits

- `ambient_de`'s false discovery proportion at a nominal 0.10 is 0.10 to 0.17
  from 50 per group up. The null false-positive rate is near nominal, so these
  are a handful of genes with very small residual shifts; `lfc` removes them.
- With 6 per group and outliers the null rate is 0.073 and the false
  discovery proportion 0.27 to 0.41 on few calls. Outliers cannot be
  identified with so few samples.
- The outlier weighting costs about half a point on clean data (0.050 to
  0.064 against 0.044 to 0.057 without it).
- `huber_k = 2` and `max_ambient_frac = 0.75` were chosen on this simulator.
- The simulator generates data from the additive ambient model that
  `ambient_de` fits. Nothing here is real data.
