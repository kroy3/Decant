# Validation scripts

Run from the repo root. Results go to `validation/results/`.

| script | ground truth | status |
|---|---|---|
| `de_comparators.R` | simulated case/control study | run; results in `results/de_comparators.csv` |
| `species_mix.R` | 10x human/mouse mixtures | **not yet run** (needs `cf.10xgenomics.com`) |
| `kidney_mixture.R` | Janssen et al. 2023 genotype mixture | **not yet run** (needs `zenodo.org`) |

## Comparator wrappers (`comparators.R`)

Official packages only. A tool that is not installed is reported as
"not installed" and never replaced by a reimplementation.

| tool | install | wrapper status |
|---|---|---|
| SoupX | `install.packages("SoupX")` | run (SoupX 1.6.2, R 4.3.3) |
| DecontX | `BiocManager::install("decontX")` | written, **not yet run** |
| FastCAR | `remotes::install_github("LungCellAtlas/FastCAR")` | written against the FastCAR source, **not yet run** |
| CellBender | `pip install cellbender`, plus R package `hdf5r` | written against the v0.3 docs, **not yet run** |

Do not quote a number from a wrapper marked "not yet run" until its output
has been checked on a case with a known answer.

`species_mix.R` still carries its own SoupX and DecontX wrappers. Pointing it
at `comparators.R` is left for when it is first run on real data, because
FastCAR and CellBender return no per-cell rho and its scoring expects one.

## Reproduce the DE comparator table

```bash
Rscript validation/de_comparators.R 1 2 3 4
```
