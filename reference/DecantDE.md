# Contamination-aware DE for one cell type of a multi-sample Seurat object.

Pseudobulks the chosen cell type per sample, measures each sample's soup
from its own empty droplets, and tests a sample-level design with
[`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md).
Unlike correcting counts and then testing, this stops an abundant cell
type's condition-dependent genes from being reported as DE in the target
type when ambient RNA differs between samples.

## Usage

``` r
DecantDE(
  object,
  cell_type,
  sample_col,
  formula,
  raw = NULL,
  empties = NULL,
  celltype_col = NULL,
  coef = NULL,
  assay = "RNA",
  min_cells = 10,
  empty_umi_range = c(1, 100),
  ...
)
```

## Arguments

- object:

  a Seurat object with raw UMI counts and cells from several samples.

- cell_type:

  the cell type to test (a value of `celltype_col`).

- sample_col:

  meta.data column naming each cell's sample.

- formula:

  sample-level design, e.g. `~ condition` or `~ batch + condition`.
  Every variable must be a meta.data column that is constant within each
  sample.

- raw:

  named list, one entry per sample (names = values of `sample_col`), of
  raw droplet matrices or Cell Ranger `raw_feature_bc_matrix` paths. Or
  supply `empties` instead.

- empties:

  named list of empty-droplet matrices per sample.

- celltype_col:

  meta.data column of cell-type labels; NULL uses `Idents(object)`.

- coef:

  design column to test; default the last one (the last term of
  `formula`).

- assay:

  assay holding raw counts.

- min_cells:

  samples with fewer target cells are dropped (with a warning).

- empty_umi_range:

  UMI range defining empty droplets.

- ...:

  passed to
  [`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md).

## Value

a data.frame as from
[`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md),
with `ambient_driven` marking genes a standard edgeR pseudobulk test
calls DE but the contamination- aware test does not. Attributes
`"samples"` (per-sample ambient fraction, soup depth, cells) and
`"design"`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- DecantDE(obj, cell_type = "Microglia", sample_col = "donor",
                formula = ~ diagnosis, raw = raw_paths_by_donor)
head(res)
subset(res, ambient_driven)   # standard-pipeline hits explained by soup
} # }
```
