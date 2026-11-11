# Score decontamination on a species-mixing experiment (full matrix).

Doublets are excluded from all metrics.

## Usage

``` r
score_species_mix(
  cells,
  empties,
  species,
  corrected,
  rho_hat = NULL,
  truth = NULL
)
```

## Arguments

- cells, empties, species:

  as for
  [`species_mix_truth()`](https://kroy3.github.io/Decant/reference/species_mix_truth.md).

- corrected:

  corrected counts (same shape as `cells`).

- rho_hat:

  per-cell rho estimate (NULL to skip rho metrics).

- truth:

  optional precomputed
  [`species_mix_truth()`](https://kroy3.github.io/Decant/reference/species_mix_truth.md)
  result.

## Value

one-row data.frame: rho bias / RMSE / correlation, the fraction of
cross-species (known-contaminant) counts removed, and own-species
removal relative to the expected own-species contamination (1 = right
amount, above 1 = over-correction destroying real signal).
