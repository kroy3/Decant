# Remove ambient counts from each cell (mass-conserving).

Remove ambient counts from each cell (mass-conserving).

## Usage

``` r
correct_counts(
  observed,
  soup,
  rho,
  method = c("subtract", "redistribute", "bayes"),
  clusters = NULL,
  max_iter = 50,
  drop_zeros = TRUE
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

  one of "subtract", "redistribute", "bayes" (see Details).

- clusters:

  cluster labels, required for method = "bayes".

- max_iter:

  iterations for "redistribute".

- drop_zeros:

  for sparse input, drop entries corrected to exactly 0. FALSE keeps the
  input's sparsity pattern (explicit zeros).

## Value

corrected counts, same class as `observed` (dense or dgCMatrix).
Guaranteed 0 \<= corrected \<= observed element-wise.

## Details

See the header of `R/correct.R` for the three rules and
[`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md)
for how they compare against ground truth.
