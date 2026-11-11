# Hidden-species view of a mixture: one species' cells, only its genes.

The estimator sees no cross-species genes, so it must work as on
ordinary single-species data. Truth is rescaled to what is estimable
from these genes: contamination among the kept genes' counts, rho_hidden
= rho_true \* s_own \* T / own.

## Usage

``` r
species_mix_hidden(cells, empties, species, keep_species = NULL, truth = NULL)
```

## Arguments

- cells, empties, species:

  as for
  [`species_mix_truth()`](https://kroy3.github.io/Decant/reference/species_mix_truth.md).

- keep_species:

  which genome to keep (default: the majority genome).

- truth:

  optional precomputed
  [`species_mix_truth()`](https://kroy3.github.io/Decant/reference/species_mix_truth.md)
  result.

## Value

list(cells, empties, rho_true, barcodes).
