# Decant

**Ambient RNA decontamination for single-cell and single-nucleus RNA-seq,
with a differential expression test that models ambient RNA instead of
subtracting it.**

Every droplet experiment captures RNA that did not come from the cell in the
droplet: transcripts released by lysed cells float in the suspension and end
up in every droplet. This *ambient RNA*, or *soup*, puts the markers of
abundant cell types into every other cell. It blurs rare populations, and it
can make differential expression report the soup instead of biology. Decant
measures the soup from your own empty droplets, estimates how much of each
cell is ambient, and removes it.

[![Documentation Status](https://readthedocs.org/projects/decant/badge/?version=latest)](https://decant.readthedocs.io/en/latest/)

Source code: [github.com/kroy3/Decant](https://github.com/kroy3/Decant) ·
Changes: [NEWS](https://github.com/kroy3/Decant/blob/main/NEWS.md) ·
Issues: [GitHub issues](https://github.com/kroy3/Decant/issues)

---

## Install

```r
# install.packages("remotes")
remotes::install_github("kroy3/Decant")
```

Decant is an R package. Seurat is optional; edgeR is used for differential
expression. See [Installation](guide/01-Installation.md).

## Quick start

You need the **raw (unfiltered)** matrix: its empty droplets measure the
ambient RNA.

```r
library(Decant)
raw   <- read_10x_counts("sample/outs/raw_feature_bc_matrix")
cells <- colnames(read_10x_counts("sample/outs/filtered_feature_bc_matrix"))
d     <- split_droplets(raw, cells = cells)

res <- decant(d$cells, d$empties, clusters = my_clusters)
res$corrected   # corrected counts (sparse), never above the observed counts
res$rho         # contamination fraction per cell
```

From Seurat:

```r
obj <- RunDecant(obj, raw = "sample/outs/raw_feature_bc_matrix")   # adds a "Decant" assay
```

The full walk-through is in [Quick Start](guide/03-Quick-Start.md).

## What Decant does

| Function | Purpose | Status |
|---|---|---|
| `decant()`, `RunDecant()` | Estimate contamination per cell and correct counts | stable |
| `correct_counts(method = "bayes_absent")` | Default correction: decides per gene and cluster whether a gene is expressed at all, and removes all counts of genes that are not; keeps counts too high to be ambient (doublets, rare expressers). Removes the most ambient RNA on species mixtures, at a small cost in retained native signal; no gain shown yet in real tissue | default |
| `correct_counts(method = "bayes_cell")` | Removes each count's expected ambient share only; keeps the most native signal | stable, previous default |
| `ambient_de()`, `DecantDE()` | Differential expression with each sample's ambient RNA inside the negative binomial model: `mu = native + ambient` | experimental |
| `chirp()` | Ligand-receptor calls that do not mistake soup for signalling | experimental, simulation only |

Three rules hold throughout:

1. **The soup is measured, not guessed.** It comes from your empty droplets.
2. **Counts are never added.** Every corrected count lies between zero and
   the observed count, which is checked on every run.
3. **Every default is tested.** A default is enabled only after it passes a
   test against data with a known answer. Tests that went against Decant
   are published with the rest.

## How it compares

All numbers are per-gene scores from the benchmark of Cargnelli et al.
(2026): **retain** is the share of a cell's own signal kept and **removal**
the share of ambient signal removed (1 is best). They are measured on 10x
human/mouse mixtures, where cross-species counts give the truth. All tools
were run on the same inputs.

| | Decant default (`bayes_absent`) | Decant `bayes_cell` | DecontX | SoupX |
|---|---|---|---|---|
| Removal, two mixtures never used in design | **0.891** | 0.531 | 0.655 | 0.529 |
| Retain, same mixtures | 0.9918 | **0.9945** | 0.9920 | 0.9931 |
| Retain, hardest mixture seen (hgmm 10k) | 0.961 | **0.986** | not run here | not run here |
| Retain, Smart-seq2 control with no ambient RNA | **0.998** | **0.998** | 0.952 | 0.989 |

**What this means in practice:**

- **On species mixtures, the default removes the most of any tool run here.** It keeps about as
  much as DecontX: 0.02 points less on the pre-registered test. That missed
  the pre-stated bar; it was made the default by the maintainer's decision,
  and the record says so.
- **It costs native signal on hard data.** Against `bayes_cell` it keeps
  0.3 to 0.6 points less on most mixtures and 2.5 points less on the
  hardest. If keeping every native count matters most, use the previous
  default:

  ```r
  res <- decant(d$cells, d$empties, clusters = my_clusters, correction = "bayes_cell")
  ```

- **Its gain is shown only on species mixtures.** In real tissue its
  gene-level call almost never applies, and it then behaves as `bayes_cell`.
  On PBMC 4k it never applied and the output matched. In brain (six GSE272082
  libraries, pre-registered) it applied in 1 of 42 clusters, and microglia
  lost the same share of genes foreign to them as with `bayes_cell` (0.092
  against 0.091). Oligodendrocytes and astrocytes lost 2.7 points more
  foreign signal, from its cap rule, not its gene-level call; part of that is
  counts too high to be ambient.
- **Tiny clusters:** the one brain cluster where it applied had 9 cells, and
  it called 17,543 genes absent there, which looks like a degenerate fit.
  Avoid it on clusters of a few cells, or merge them first.

- **Neither rule removes counts that cannot be ambient.** On real PBMC data
  both keep counts far above what ambient RNA could produce, which DecontX,
  AmbiDose and scAR remove.

## Limits

Read these before relying on the output.

- **Contamination estimates can run low.** On one of four mixtures (hgmm
  10k) Decant underestimates contamination and ranks cells poorly; DecontX
  does better there.
- **A single cell type gives no estimate.** When the soup looks like the
  cells, contamination cannot be estimated from expression. Decant warns
  and asks for `rho`; the number it reports then is not reliable.
- **Skewed empty droplets are not detected.** Empty droplets that
  over-represent genes the cells also express make the estimate too low.
- **`ambient_de()` is experimental.** On real data it is checked for
  calibration only: with shuffled labels it reports hits a little more often
  than standard edgeR. That it removes ambient-driven false hits is shown
  only in simulation. Two pre-registered real-data tests were null, because
  the standard test found no such hits to remove. Report its hits together
  with a standard pseudobulk test.
- **`chirp()` is untested on real data.**
- **Ground truth comes from cell lines.** It is human/mouse mixtures and one
  negative control. The paper's genotype mixture and medium- and
  high-complexity datasets have not been run.

## User guide

1. [Installation](guide/01-Installation.md)
2. [Key Concepts](guide/02-Key-Concepts.md)
3. [Quick Start](guide/03-Quick-Start.md)
4. [Preparing Your Data](guide/04-Preparing-Your-Data.md)
5. [Running Decant](guide/05-Running-Decant.md)
6. [Seurat Integration](guide/06-Seurat-Integration.md)
7. [Estimation and Correction](guide/07-Estimation-and-Correction.md)
8. [Multi-Sample Experiments](guide/08-Multi-Sample-Experiments.md)
9. [Differential Expression](guide/09-Differential-Expression.md)
10. [Splice-Aware and Genotype Modules](guide/10-Splice-Aware-and-Genotype-Modules.md)
11. [Diagnostics and Troubleshooting](guide/11-Diagnostics-and-Troubleshooting.md)
12. [Validation and Benchmarks](guide/12-Validation-and-Benchmarks.md)
13. [Function Reference](guide/13-Function-Reference.md)
14. [Cell-Cell Communication](guide/14-Cell-Cell-Communication.md)

## Validation

The validation is published in full in the source repository:

- [`validation/BENCHMARK.md`](https://github.com/kroy3/Decant/blob/main/validation/BENCHMARK.md):
  every benchmark and the changes that failed.
- [`validation/preregistration/`](https://github.com/kroy3/Decant/tree/main/validation/preregistration):
  tests whose criteria were committed before their data were seen.
- [`validation/SCORECARD.md`](https://github.com/kroy3/Decant/blob/main/validation/SCORECARD.md):
  where Decant is ahead, level or behind.

## Citation and license

Benchmark metrics follow Cargnelli et al. (2026). License: see
[LICENSE](https://github.com/kroy3/Decant/blob/main/LICENSE).
