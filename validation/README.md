# Validation scripts

Run from the repo root. Results go to `validation/results/`.

| script | ground truth | status |
|---|---|---|
| `benchmark_differential_expression.R` | simulated case/control study | run; results in `results/differential_expression.csv` |
| `validate_pbmc_markers.R` | none (10x PBMC 4k bundled with SoupX; compared with SoupX's published estimates and marker specificity) | run; results in `results/pbmc_*.csv` |
| `benchmark_pbmc_six_tools.R` + `py/run_python_tools.py` | none (PBMC 4k); six tools | run; results in `results/pbmc_six_tools.csv` |
| `benchmark_species_mixture.R` | species mixing; the two scores of Cargnelli et al. (2026) | self-test run on a synthetic mixture; `low_complexity` run; results in `results/species_mixture_scores.csv` (BENCHMARK.md section 1c) |
| `benchmark_species_mixture_rho.R` | cross-species reads give each cell's contamination | run; results in `results/species_mixture_rho.csv` |
| `benchmark_negative_control.R` | Smart-seq2 negative control (E-MTAB-5061), no ambient RNA | run; results in `results/negative_control_scores*.csv` |
| `validate_species_mixture.R` | 10x human/mouse mixtures | **not yet run** |
| `validate_genotype_mixture.R` | Janssen et al. 2023 genotype mixture | **not yet run** (needs `zenodo.org`) |

## Comparator wrappers (`comparators.R`)

Official packages only. A tool that is not installed is reported as
"not installed" and never replaced by a reimplementation.

| tool | install | wrapper status |
|---|---|---|
| SoupX | `install.packages("SoupX")` | run (SoupX 1.6.2, R 4.3.3) |
| DecontX | `BiocManager::install("decontX")` | run (decontX 1.10.0) |
| FastCAR | `remotes::install_github("LungCellAtlas/FastCAR")` | written against the FastCAR source, **not yet run** |
| CellBender | `pip install cellbender`, plus R package `hdf5r` | written against the v0.3 docs, **not yet run** |

Do not quote a number from a wrapper marked "not yet run" until its output
has been checked on a case with a known answer.

`validate_species_mixture.R` still carries its own SoupX and DecontX wrappers. Pointing it
at `comparators.R` is left for when it is first run on real data, because
FastCAR and CellBender return no per-cell rho and its scoring expects one.

## Reproduce the DE comparator table

```bash
Rscript validation/benchmark_differential_expression.R 1 2 3 4
```
