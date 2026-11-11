# Split a raw droplet matrix into cells and empty droplets.

Cells are either the barcodes you pass (e.g. from Cell Ranger's filtered
matrix, EmptyDrops or CellBender calls – recommended) or, failing that,
droplets with at least `cell_min_umi` counts. Empties are droplets whose
total lies in `empty_umi_range` (SoupX's default soup range is 0-100
UMIs) and that are not cells. Droplets in between are ambiguous and are
dropped from both sets rather than risk calling a cell "soup".

## Usage

``` r
split_droplets(
  raw,
  cells = NULL,
  cell_min_umi = 500,
  empty_umi_range = c(1, 100)
)
```

## Arguments

- raw:

  genes x barcodes raw matrix (see
  [`read_10x_counts()`](https://kroy3.github.io/Decant/reference/read_10x_counts.md)).

- cells:

  optional character vector of cell barcodes.

- cell_min_umi:

  UMI threshold used only when `cells` is NULL.

- empty_umi_range:

  inclusive UMI range for empty droplets. The lower bound defaults to 1
  because zero-count barcodes carry no information.

## Value

list(cells = genes x cells, empties = genes x droplets, summary).
