# Ground-truth contamination from a two-species mixture.

Ground-truth contamination from a two-species mixture.

## Usage

``` r
species_mix_truth(cells, empties, species, max_minor = 0.25)
```

## Arguments

- cells:

  genes x cells counts (both genomes' genes).

- empties:

  genes x droplets counts from empty droplets.

- species:

  length-nrow genome label per gene (see
  [`gene_species()`](https://kroy3.github.io/Decant/reference/gene_species.md)).

- max_minor:

  cells whose minority-species fraction exceeds this are called doublets
  and excluded. Contamination alone rarely exceeds ~0.2 of the VISIBLE
  (cross-species) fraction.

## Value

data.frame per cell: barcode, species, total, own, other, other_frac,
rho_true, doublet. Attribute `"soup_share"`: soup fraction per genome.
