# Read a 10x Genomics matrix directory into a sparse matrix.

Read a 10x Genomics matrix directory into a sparse matrix.

## Usage

``` r
read_10x_counts(dir, gene_column = 1, feature_type = "Gene Expression")
```

## Arguments

- dir:

  path to a Cell Ranger `raw_feature_bc_matrix` (or filtered) directory
  containing matrix.mtx, features.tsv (or genes.tsv) and barcodes.tsv,
  optionally gzipped.

- gene_column:

  which column of features.tsv to use as row names (1 = Ensembl ID, 2 =
  symbol). IDs are the safe default: symbols repeat.

- feature_type:

  keep only features of this type (column 3 of features.tsv), e.g. "Gene
  Expression". NULL keeps all.

## Value

genes x barcodes dgCMatrix.
