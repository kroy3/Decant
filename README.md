# Decant

Ambient RNA decontamination for single-cell and single-nucleus RNA-seq, with a
differential expression test that models ambient RNA instead of subtracting
it. No other ambient RNA tool we found offers one.

**Documentation:** the [User Guide](https://github.com/kroy3/Decant/wiki)
covers installation, the Seurat workflow, multi-sample data, differential
expression and diagnostics, chapter by chapter.

## What is new

1. **Contamination-aware differential expression** (`ambient_de()`,
   `DecantDE()`). When disease tissue releases more ambient RNA than controls,
   an abundant cell type's disease genes leak into every other cell type and
   are reported as differentially expressed. Existing tools correct counts and
   then test. Decant puts each sample's measured ambient RNA, and its
   uncertainty, inside the negative binomial model:
   `mu = native expression + ambient`. It supports covariates and is robust to
   outlier samples in large cohorts. **Status: experimental.** On real data
   with shuffled labels it was calibrated on two brain datasets but not on a
   pancreas dataset, and it has not yet been shown to beat a standard test
   on real data (see [Limits](#limits)).
2. **Correction that keeps what cannot be ambient** (`bayes_cell`). Counts far
   above what ambient RNA could produce in a cell (doublets, rare expressers)
   keep their molecules; only the expected ambient share is removed.
3. **Never adds counts.** Corrected counts are bounded by the observed counts,
   checked on every run.
4. **Every default is gated.** Each component must pass a falsification test
   against ground truth in CI. Components that failed are kept off and
   documented.

## How it compares

### Correction against experimental ground truth

Two datasets from the benchmark of Cargnelli et al. (2026), scored
with the paper's metrics. Retain is the share of the cell's own signal that
is kept, and removal the share of ambient signal removed, both per gene, 1
is best. Decant, SoupX and DecontX were run on the same inputs. CellBender's
numbers are the paper's own, not a run here.

| | Decant | SoupX | DecontX | CellBender (published) |
|---|---|---|---|---|
| human/mouse mixture: retain | 0.989 | 0.992 | 0.983 | 0.996 |
| human/mouse mixture: **removal** | **0.616** | 0.345 | **0.722** | 0.729 |
| Smart-seq2, no ambient RNA: retain | 0.998 | 0.989 | 0.952 | 0.797 |

Decant still removes about 10 points less ambient RNA than DecontX or
CellBender on the mixture, and keeps slightly more of the cells' own counts
than DecontX. It removes almost nothing where there is nothing to remove.
Its contamination estimate now matches the cross-species truth (median 0.027
against 0.029; before this release 0.013). That mixture was used to find the
estimator's flaw; results on held-out mixtures are in
[`validation/BENCHMARK.md`](validation/BENCHMARK.md), sections 1c and 1d.

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
| GSE273597 pancreas; ductal, alpha | **11 of 48** | **2,300** | 4 of 48 |

**`ambient_de` failed on pancreas ductal and alpha cells.** For 6 of 23
donors, the ductal ambient fraction came out at 0.86 to 0.96 (the others 0.02
to 0.16): a single cell type's pseudobulk cannot always identify its own
ambient fraction. The model then treated most of those donors' counts as
soup, and any split that separated them gave over a thousand false hits.
This release replaces such values (above 0.5) with the median of the other
samples and warns; its effect on these data is in
[`validation/BENCHMARK.md`](validation/BENCHMARK.md), section 3c. The
pancreas data were used to find that problem, so they cannot show the fix
works elsewhere.

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
  calibration only (three datasets, shuffled labels) and failed in two
  pancreas cell types before the guard in this release. That it removes
  ambient-driven hits a standard test reports is shown only in simulation,
  from a simulator that follows the model's own mechanism. A pre-registered
  real-data test (pancreas, T1D vs control) could not show it: the standard
  test found no ambient-driven hits to remove (`validation/BENCHMARK.md`,
  3c). Until there is real-data evidence, report `ambient_de` hits together
  with a standard pseudobulk test. Analyses made with earlier versions
  should be rerun.
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
# the GitHub wiki is updated from wiki/ by .github/workflows/wiki.yaml
```

## License

MIT.
