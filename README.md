# Decant

[![Documentation Status](https://readthedocs.org/projects/decant/badge/?version=latest)](https://decant.readthedocs.io/en/latest/)

Ambient RNA decontamination for single-cell and single-nucleus RNA-seq, with a
differential expression test that models ambient RNA instead of subtracting
it. No other ambient RNA tool we found offers one.

**Documentation:** [decant.readthedocs.io](https://decant.readthedocs.io/en/latest/).
The user guide covers installation, the Seurat workflow, multi-sample data,
differential expression and diagnostics, chapter by chapter (also on the
[GitHub wiki](https://github.com/kroy3/Decant/wiki)).

## What is new

1. **Contamination-aware differential expression** (`ambient_de()`,
   `DecantDE()`). When disease tissue releases more ambient RNA than controls,
   an abundant cell type's disease genes leak into every other cell type and
   are reported as differentially expressed. Existing tools correct counts and
   then test. Decant puts each sample's measured ambient RNA, and its
   uncertainty, inside the negative binomial model:
   `mu = native expression + ambient`. It supports covariates and is robust to
   outlier samples in large cohorts. **Status: experimental.** With shuffled
   labels on three real datasets it now rarely reports false hits, a little
   more often than standard edgeR; it has not yet been shown to beat a
   standard test on real data (see [Limits](#limits)).
2. **Correction that decides whether a gene is expressed at all**
   (`bayes_absent`, the default). For each gene in each cluster it computes
   the probability that the cluster does not express the gene, from the
   gene's count in the cluster against its count in the empty droplets, and
   removes all counts of genes it does not express. Counts far above what
   ambient RNA could produce in a cell (doublets, rare expressers) keep their
   molecules. On two 10x mixtures never used in its design it removed the
   most ambient RNA of every tool run here.
3. **Never adds counts.** Corrected counts are bounded by the observed counts,
   checked on every run.
4. **Ambient-aware ligand-receptor calls** (`chirp()`). A ligand made by the
   cell type that dominates the soup is detected in every cluster, so the
   usual detection rule makes every cluster look like a sender. `chirp()`
   calls a gene in a cluster only when its count is well above what ambient
   RNA would put there, and scores pairs on the excess. **Status:
   experimental, simulation only** (see [Limits](#limits)).
5. **Every default is gated.** Each component must pass a falsification test
   against ground truth in CI. Components that failed are kept off and
   documented.

## How it compares

### Correction against experimental ground truth

Two datasets from the benchmark of Cargnelli et al. (2026), scored
with the paper's metrics. Retain is the share of the cell's own signal that
is kept, and removal the share of ambient signal removed, both per gene, 1
is best. Decant, SoupX and DecontX were run on the same inputs. CellBender's
numbers are the paper's own, not a run here.

| | Decant `bayes_absent` (default) | Decant `bayes_cell` | SoupX | DecontX | CellBender (published) |
|---|---|---|---|---|---|
| 2 mixtures never used in design (hgmm 1k): retain | 0.9918 | **0.9945** | 0.9931 | 0.9920 | not run |
| 2 mixtures never used in design (hgmm 1k): **removal** | **0.891** | 0.531 | 0.529 | 0.655 | not run |
| mixtures 1k, 5k, 10k v3 (seen in design): retain | 0.982 | **0.992** | 0.982 | 0.984 | not run |
| mixtures 1k, 5k, 10k v3 (seen in design): **removal** | **0.862** | 0.444 | 0.679 | 0.637 | not run |
| mixture 20k: retain / removal | not rerun | 0.989 / 0.616 | 0.992 / 0.345 | 0.983 / 0.722 | 0.996 / 0.729 |
| Smart-seq2, no ambient RNA: retain | 0.998 | 0.998 | 0.989 | 0.952 | 0.797 |

**The default removes the most ambient RNA of every tool run here, at a
small cost in the cells' own signal.** On the two mixtures that judged it
(downloaded only after the method was frozen and the criteria were
committed) it removed 0.891 of ambient signal against DecontX's 0.655, and
kept 0.9918 against DecontX's 0.9920. That missed the pre-stated bar of
keeping at least as much as DecontX by 0.02 points; it was made the default
by the maintainer's decision anyway (`validation/BENCHMARK.md`, 1h).

What it costs, against the previous default `bayes_cell`:

- 0.3 to 0.6 points of native signal on most mixtures, and 2.5 points on the
  hardest (hgmm 10k: 0.961 against 0.986).
- Its gain is shown only where clusters differ sharply. On real blood data
  (PBMC 4k) the clusters share most genes, its gene-level call never
  applied, and its output matched `bayes_cell`.

Use `correction = "bayes_cell"` when keeping every native count matters
more than removing the soup. Results on further mixtures and the history
of the contamination estimate are in
[`validation/BENCHMARK.md`](validation/BENCHMARK.md), sections 1c to 1h.

### Correction: six tools on real 10x PBMC 4k data

Same cells, soup and clusters for every tool. "Ambient-compatible" foreign
counts are marker counts in cell types that do not express the marker, at a
level ambient RNA could explain; "not ambient" are those too high to be
ambient, mostly B-cell doublets. A good tool removes the first and keeps the
second.

| | Decant | SoupX | DecontX | CellBender¹ | AmbiDose | scAR |
|---|---|---|---|---|---|---|
| ambient-compatible foreign counts removed | 88–95% | 76–94% | 97–100% | 54–98% | 67–79% | 99–100% |
| counts that are not ambient removed (lower is better) | 14–37% | 12–31% | 96–100% | 4–68% | 82–100% | 58–82% |
| native marker signal lost | 1–4% | 1–2% | 2–6% | 0–1% | 0–6% | 18–26% |
| counts added | 0 | 0 | 0 | 0 | 0 | 727,390 |
| run time (4 CPU cores) | 11 s | 12 s | 34 s | 55 min | 73 s | 15 min |

¹ Reduced CPU configuration; CellBender is designed for a GPU.

Decant removes more ambient-compatible counts than SoupX and AmbiDose on all
three marker sets and than CellBender (as run here) on two of three. It keeps
the counts that are not ambient, which DecontX, AmbiDose and scAR remove, and
it is as fast as the fastest tool. DecontX removes 5 to 9 points
more of the ambient-compatible counts. Details and caveats:
[`validation/BENCHMARK.md`](validation/BENCHMARK.md).

### Differential expression (simulated case/control studies)

| | leaked genes called (of 80), 6 vs 6 | null false-positive rate, 6 vs 6 | null false-positive rate, 50 vs 50, confounding batch |
|---|---|---|---|
| **`ambient_de`** | **0** | **0.05** | **0.05** |
| SoupX, then edgeR | 1.5 | 0.04 | not run |
| naive edgeR | 53 | 0.22 | 0.45 |
| naive Wilcoxon | not run | not run | 0.46 |
| subtract ambient, then Wilcoxon | not run | not run | 0.14 |

Nominal false-positive rate 0.05. In the 6 vs 6 study, correcting with
Decant's `bayes` rule and then running edgeR was as good as `ambient_de`;
`ambient_de` is the option that stays valid with covariates such as batch.

**Calibration on real data.** Condition labels shuffled between donors so
that there is nothing to find, 24 shuffles per cell type. A calibrated test
reports a hit in a few shuffles, and few genes when it does.

| data, cell types | `ambient_de`: shuffles with a hit | most hits in one shuffle | naive edgeR: shuffles with a hit |
|---|---|---|---|
| GSE272082 brain, AD vs control; microglia, oligodendrocytes | 4 of 48 | 9 | 1 of 48 |
| GSE190015 brain, nine controls; microglia, oligodendrocytes | 5 of 48 | 1 | 3 of 48 |
| GSE273597 pancreas, T1D vs control; acinar | 1 of 24 | 1 | 0 of 24 |
| GSE273597 pancreas; ductal | 4 of 24 (before this release: 5, up to **2,300** genes) | 1 | 0 of 24 |
| GSE273597 pancreas; alpha | 6 of 24 | 5 | 4 of 24 |

On the pancreas data, for 6 of 23 donors the ductal ambient fraction came out
at 0.86 to 0.96 (the others 0.02 to 0.16): a single cell type's pseudobulk
cannot always identify its own ambient fraction. The model treated most of
those donors' counts as soup, and splits that separated them gave over a
thousand false hits. `ambient_de` now replaces a value above `max_rho` (0.5)
with the median of the other samples, adds their spread to its uncertainty,
and warns. That removed the false hits; it changed nothing on the two brain
datasets. The pancreas data were used to find the problem, so they do not
show the guard works elsewhere. `ambient_de` is slightly more liberal than
naive edgeR overall (hits in 20 of 168 shuffles against 8, never more than 9
genes). Details: [`validation/BENCHMARK.md`](validation/BENCHMARK.md),
sections 3b and 3c.

An earlier problem, the cause of hits in most shuffles on the brain data
(82 of 96 tests, up to 574 genes), was fixed by using edgeR's trended
dispersion (section 3b).

## Install

```r
remotes::install_github("kroy3/Decant")
```

## Quick start

The raw (unfiltered) matrix is required: its empty droplets measure the
ambient RNA.

```r
library(Decant)
raw <- read_10x_counts("sample/outs/raw_feature_bc_matrix")
cells <- colnames(read_10x_counts("sample/outs/filtered_feature_bc_matrix"))
d <- split_droplets(raw, cells = cells)
res <- decant(d$cells, d$empties, clusters = my_clusters)
res$corrected   # corrected counts, sparse
res$rho         # per-cell contamination
```

From Seurat:

```r
obj <- RunDecant(obj, raw = "sample/outs/raw_feature_bc_matrix")  # new "Decant" assay

res <- DecantDE(obj, cell_type = "Microglia", sample_col = "donor",
                formula = ~ diagnosis, raw = raw_paths_by_donor)
subset(res, ambient_driven)   # hits a standard pipeline reports that the soup explains
```

## Limits

- Ground truth so far is four human/mouse mixtures and one negative control,
  all from cell lines. On the 10k mixture Decant underestimates contamination
  (0.021 against 0.062) and ranks cells poorly; DecontX does better there.
  The genotype mixture and the paper's medium- and high-complexity datasets
  are not yet run.
- When the soup resembles the cells (a single cell type), contamination
  cannot be estimated from expression. Decant then warns and asks for `rho`;
  the number it reports in that case is not reliable.
- Memory: about 9.5 GB for 11,000 cells and 1.9 million empty droplets.
- **`ambient_de()` is experimental.** On real data it has been checked for
  calibration only (three datasets, shuffled labels); before this release it
  failed in pancreas ductal cells. That it removes
  ambient-driven hits a standard test reports is shown only in simulation,
  from a simulator that follows the model's own mechanism. Two pre-registered
  real-data tests (pancreas, T1D vs control; brain regions within donors)
  could not show it: the standard test found no ambient-driven hits to
  remove (`validation/BENCHMARK.md`, 3c and 3d). Until there is real-data evidence, report `ambient_de` hits together
  with a standard pseudobulk test. Analyses made with earlier versions
  should be rerun.
- **`chirp()` is experimental and untested on real data.** On simulated
  data it made no false sender calls where the 10% detection rule made
  608-944 and correcting with `bayes_cell` first still left 70-388. The
  simulator matches its model, CellChat and CellPhoneDB were not run, and
  it misses most pairs expressed at 5-20 per 10k. It relies on `rho`; with
  `rho` at half its true value it made up to 100 false calls.
- Empty droplets that over-represent genes the cells also express are not
  detected and make the contamination estimate too low.

The full scorecard, design notes and history are in
[`validation/SCORECARD.md`](validation/SCORECARD.md); changes are in
[`NEWS.md`](NEWS.md).

## Development

```r
devtools::test()                              # unit tests
Rscript scripts/run_falsification_gates.R     # falsification gates (also run in CI)
Rscript scripts/check_user_guide.R            # run every example in the user guide
Rscript scripts/build_user_guide.R            # copy wiki/ into the website articles
# the GitHub wiki is updated from wiki/ by .github/workflows/wiki.yaml;
# Read the Docs rebuilds from main (.readthedocs.yaml, mkdocs.yml, rtd/)
```

## License

MIT.
