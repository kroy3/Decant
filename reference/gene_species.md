# Species of each gene from a species-prefixed feature name.

Cell Ranger names features in multi-genome references as `<genome>_<id>`
(e.g. `GRCh38_ENSG...`, `mm10___ENSMUSG...`).

## Usage

``` r
gene_species(genes)
```

## Arguments

- genes:

  character vector of feature names.

## Value

character vector of genome labels (the prefix before the first "\_").
