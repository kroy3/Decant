# Read a raw 10x matrix in a streaming pass, keeping only the cells.

Reads `matrix.mtx(.gz)` in chunks, twice: first to total each barcode,
then to keep the cells' entries and reduce the empty droplets to a
[`soup_summary()`](https://kroy3.github.io/Decant/reference/soup_summary.md).
Neither the raw matrix nor the empty-droplet matrix is ever held in
memory; peak memory is the cells plus one chunk. Cells and empties are
defined as in
[`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md).

## Usage

``` r
read_10x_streaming(
  dir,
  cells = NULL,
  cell_min_umi = 500,
  empty_umi_range = c(1, 100),
  gene_column = 1,
  feature_type = "Gene Expression",
  chunk_lines = 2e+06
)
```

## Arguments

- dir:

  a Cell Ranger `raw_feature_bc_matrix` directory (matrix.mtx,
  features.tsv or genes.tsv, barcodes.tsv; optionally gzipped).

- cells:

  optional character vector of cell barcodes (recommended).

- cell_min_umi:

  used only when `cells` is NULL.

- empty_umi_range:

  inclusive UMI range for empty droplets.

- gene_column, feature_type:

  as in
  [`read_10x_counts()`](https://kroy3.github.io/Decant/reference/read_10x_counts.md).

- chunk_lines:

  matrix entries read per chunk.

## Value

list(cells = genes x cells dgCMatrix, soup = a `decant_soup`, summary =
the same counts as
[`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md)).
