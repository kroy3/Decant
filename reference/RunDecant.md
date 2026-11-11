# Run Decant on a Seurat object.

Run Decant on a Seurat object.

## Usage

``` r
RunDecant(
  object,
  raw = NULL,
  empties = NULL,
  assay = "RNA",
  clusters = NULL,
  sample_col = NULL,
  new_assay = "Decant",
  set_default = FALSE,
  empty_umi_range = c(1, 100),
  ...
)
```

## Arguments

- object:

  a Seurat object holding raw UMI counts.

- raw:

  the raw (unfiltered) droplet matrix: a path to a Cell Ranger
  `raw_feature_bc_matrix` directory or a genes x barcodes matrix. For a
  multi-sample object, a named list of these, one per sample, with names
  matching the values of `meta.data[[sample_col]]`. Either `raw` or
  `empties` is required.

- empties:

  alternatively, the empty-droplet matrix (or named list per sample)
  directly.

- assay:

  assay holding the raw counts.

- clusters:

  cluster labels: the name of a meta.data column, a vector with one
  label per cell, or NULL to use `Idents(object)`. k-means is used only
  if the identities are a single level.

- sample_col:

  meta.data column naming each cell's sample; required when
  `raw`/`empties` is a list.

- new_assay:

  name of the assay that receives the corrected counts.

- set_default:

  make `new_assay` the default assay.

- empty_umi_range:

  UMI range defining empty droplets (see
  [`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md)).

- ...:

  passed to
  [`decant()`](https://kroy3.github.io/Decant/reference/decant.md) (e.g.
  `correction`, `rho`, `k`).

## Value

the Seurat object with: a new assay (`new_assay`) of corrected counts,
meta.data column `decant_rho`, and `Misc(object, "decant")` holding the
per-cluster rho diagnostics, the soup profile(s), the modules that ran
and the soup-mismatch flag. Run `NormalizeData()` on the new assay
before downstream use.

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- RunDecant(obj, raw = "sample/outs/raw_feature_bc_matrix")
obj <- NormalizeData(obj, assay = "Decant")
} # }
```
