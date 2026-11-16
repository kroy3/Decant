# Chapter 1. Installation

Decant is an R package. It needs R 4.1 or later and four packages that
come with every R installation: Matrix, methods, stats and utils.

## 1.1 Install Decant

``` r

install.packages("remotes")
remotes::install_github("kroy3/Decant")
```

- `install.packages("remotes")` installs the helper package that can
  install packages from GitHub. Skip it if you already have `remotes`.
- `remotes::install_github("kroy3/Decant")` downloads Decant from its
  GitHub repository, builds it and installs it into your R library.

## 1.2 Optional packages

Some features need additional packages. Decant works without them and
tells you when one is missing.

| Package | Needed for | Install with |
|----|----|----|
| SeuratObject (5.0 or later) | [`RunDecant()`](https://kroy3.github.io/Decant/reference/RunDecant.md), [`DecantDE()`](https://kroy3.github.io/Decant/reference/DecantDE.md) | `install.packages("Seurat")` |
| edgeR, limma | [`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md), [`DecantDE()`](https://kroy3.github.io/Decant/reference/DecantDE.md) | `BiocManager::install(c("edgeR", "limma"))` |

``` r

install.packages("Seurat")
install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma"))
```

- The first line installs Seurat, which also installs SeuratObject.
- `BiocManager` is the installer for Bioconductor packages.
- The last line installs edgeR and limma, which the differential
  expression functions use for normalisation and testing.

## 1.3 Check the installation

``` r

library(Decant)
packageVersion("Decant")
#> [1] '0.4.0'
```

- [`library(Decant)`](https://github.com/kroy3/Decant) loads the
  package. No message is printed when it loads successfully.
- `packageVersion("Decant")` prints the installed version. Include it
  when you report a problem.

## 1.4 What you need from your experiment

Decant measures the ambient RNA from **empty droplets**, so it needs the
**raw (unfiltered)** count matrix, not only the filtered one. For 10x
Genomics data this is the `raw_feature_bc_matrix` folder written by Cell
Ranger, next to `filtered_feature_bc_matrix`. If you only kept the
filtered matrix, you can still use Decant with an approximate soup (see
[Chapter
11.6](https://kroy3.github.io/Decant/articles/11-Diagnostics-and-Troubleshooting.md)),
but the result is less reliable.

------------------------------------------------------------------------

**Previous:**
[Preface](https://kroy3.github.io/Decant/articles/00-Preface.md) \|
**Next:** [Chapter 2, Key
Concepts](https://kroy3.github.io/Decant/articles/02-Key-Concepts.md)
