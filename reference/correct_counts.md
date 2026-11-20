# Remove ambient counts from each cell (mass-conserving).

Remove ambient counts from each cell (mass-conserving).

## Usage

``` r
correct_counts(
  observed,
  soup,
  rho,
  method = c("subtract", "redistribute", "bayes", "bayes_cell", "bayes_absent"),
  clusters = NULL,
  max_iter = 50,
  drop_zeros = TRUE,
  cap_q = 0.001,
  cap_test = c("entry", "cluster", "hybrid"),
  soup_umis = NULL
)
```

## Arguments

- observed:

  genes x cells counts (base matrix or dgCMatrix).

- soup:

  length-G ambient probability vector (one soup for all cells; for
  per-sample soups use
  [`decant()`](https://kroy3.github.io/Decant/reference/decant.md)).

- rho:

  length-cells contamination fractions.

- method:

  one of "subtract", "redistribute", "bayes", "bayes_cell",
  "bayes_absent" (see Details).

- clusters:

  cluster labels, required for the two bayes methods.

- max_iter:

  iterations for "redistribute".

- drop_zeros:

  for sparse input, drop entries corrected to exactly 0. FALSE keeps the
  input's sparsity pattern (explicit zeros).

- cap_q:

  for "bayes_cell": an entry is treated as not explainable by ambient
  RNA when its count exceeds the `1 - cap_q` Poisson quantile of its
  ambient expectation.

- cap_test:

  for "bayes_cell": how the "not explainable by ambient" test is
  applied. "entry" tests each entry on its own; "cluster" corrects it
  for the number of cells of its cluster it is applied to (Bonferroni):
  an entry is set aside only if a count that high would be unexpected
  among all those cells. With "entry", any single count of a gene whose
  expected ambient count is below `cap_q` is set aside and kept,
  although thousands of cells make some such counts expected; on the 20k
  species mixture this kept most ambient counts of rarely contaminating
  genes. "cluster" removes more but also keeps less of the cells' own
  rare counts: on held-out mixtures removal 0.519 vs 0.444, retain 0.974
  vs 0.986 on the 10k mixture (`validation/BENCHMARK.md`, section 1e),
  so it is not the default. "hybrid" (used by "bayes_absent") sets an
  entry aside if it passes the "cluster" test, or the "entry" test with
  at least 5 UMIs: a single UMI against a tiny expectation is weak
  evidence on its own, but five are not, so doublet-level counts stay
  protected.

- soup_umis:

  for "bayes_absent": total UMIs the soup profile was estimated from. If
  given, the absence test compares each gene's count in a cluster with
  its count in the empty droplets (a binomial given their total), which
  carries the soup's own sampling noise; otherwise the soup is taken as
  known.

## Value

corrected counts, same class as `observed` (dense or dgCMatrix).
Guaranteed 0 \<= corrected \<= observed element-wise.

## Details

See the header of `R/correct.R` for the rules and
[`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md)
for how they compare against ground truth.
