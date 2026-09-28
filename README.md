# Decant

A benchmark-first toolkit for ambient RNA decontamination research in
single-cell and single-nucleus RNA-seq. Decant is assembled from components,
and a component becomes a default ONLY after it passes a falsification gate
against ground truth. This is not a finished CellBender replacement and does not
claim to be one. Its value is that every capability in it has been shown to earn
its place, and the things that did not are documented as failures.

> "Decant" is a working placeholder name. Check CRAN, Bioconductor, and GitHub
> for collisions before any release.

## Install and run

```r
# install.packages("remotes")
remotes::install_github("kroy3/Decant")
```

On real 10x data (the RAW matrix is required: the empty droplets are the soup
measurement):

```r
library(Decant)
raw  <- read_10x_counts("sample/outs/raw_feature_bc_matrix")
cell_bcs <- colnames(read_10x_counts("sample/outs/filtered_feature_bc_matrix"))
d    <- split_droplets(raw, cells = cell_bcs)       # cells vs empty droplets
res  <- decant(d$cells, d$empties, clusters = my_seurat_clusters)
res                       # summary, including any warnings about rho reliability
res$corrected             # corrected counts, still a sparse dgCMatrix
res$rho_diagnostics       # per-cluster rho, absent-gene power, borrowing
```

Pass your own clusters (Seurat / SCE) when you have them; otherwise k-means is
run. Inputs may be dense matrices or `Matrix::dgCMatrix`. Sparse inputs stay
sparse end to end: a 20k-gene x 20k-cell sparse matrix (229 MB) runs in
~10-20 s at ~1.1 GB peak, where v0.1 needed a 3 GB dense intermediate.

## The design rule

Earlier work in this repo established, against ground truth, that putting a
smarter model on the same count matrix (a fancier soup profile) does not move
the sensitivity/specificity frontier. So every module here brings in ORTHOGONAL
information to break the ambient-vs-native identifiability problem, and each is
gated. Run the scorecard yourself (~1 min):

```bash
Rscript run_all_gates.R
```

## Scorecard (from `run_all_gates.R`, v0.2)

### Core components (every run depends on these)

| Component | Gate result | Default |
|-----------|-------------|---------|
| **Contamination fraction (rho)**: cluster-pooled absent-gene estimator with empirical-Bayes per-cell shrinkage (`estimate_rho_cluster`) | **PASS** on 12 settings (easy + stress simulator, rho 0 to 0.3, k-means clusters, not oracle). Mean bias within 25% everywhere; negative control holds (rho=0 gives <0.01). Easy: RMSE 0.006-0.030, per-cell r = 0.94-0.98. Stress: RMSE 0.03-0.08, r = 0.50-0.81. | ON |
| Legacy marker-ratio rho (v0.1 `estimate_rho`) | **FAIL.** Overestimated rho 1.5-4x on easy data (worst at the low rho typical of scRNA-seq) and reported 0.57-0.71 when the truth was 0.0-0.3 under stress. | kept only for reproducibility |
| **Correction rule** | `subtract` (v0) has the lowest mean L1 error (0.071), vs `redistribute` (SoupX-style, 0.073) and `bayes` (0.076). The alternatives remove more contamination but destroy more signal; none dominates. | `subtract`; others via `correction=` |

What the rho fix buys end to end (`compare_methods.R`, mean of 4 seeds,
L1 error to the clean truth, lower is better):

| true rho | no correction | v0.1 pipeline | v0.2 `decant()` |
|---------:|--------------:|--------------:|----------------:|
| 0.1 | 0.109 | **0.147 (worse than doing nothing)** | 0.042 |
| 0.3 | 0.426 | 0.176 | 0.101 |

At rho=0.1 the v0.1 pipeline destroyed 12.5% of real signal through
over-correction.

### Modules

| Gap | Module | Gate result | Default |
|-----|--------|-------------|---------|
| 1 | Splice-aware layer decontamination | PASS. Cuts unspliced-layer error by 25% (0.044 vs 0.059) when ambient is splice-distinct. Ties exactly when it is not (negative control holds). Outputs corrected spliced/unspliced layers for RNA velocity. | ON |
| 3 | Hierarchical multi-sample soup pooling | PASS. Reduces soup error for empty-poor samples; neutral for empty-rich. | ON |
| 5 | Allelic / genotype-aware rho | PASS. Near-ground-truth rho in pooled designs. (Gate contrast is soft, 0.015 vs 0.017 RMSE; the absolute accuracy is the real point.) | ON when allelic data present |
| 4 | Structured-soup lysis diagnostic | PASS as QC only, with an oracle basis. It is a diagnostic readout, never a corrector (that use was falsified). | ON as diagnostic |
| 2 | Decontamination-aware differential expression | FAIL. With the benchmark fixed (see below), the covariate model has 90% false positives on ambient-trap genes vs 8% for naive subtract-then-test, and 16% on null genes. | OFF (experimental) |

## Corrections to v0.1 claims

This release re-ran every gate on the shipped estimator. Three earlier claims
did not survive, and they are corrected here rather than quietly edited:

1. **"Splice-aware halves unspliced-layer error."** The 2x came partly from the
   legacy rho over-correcting the baseline more heavily. With a calibrated rho
   the gain is 25%, and the negative control initially FAILED because each
   layer was clamped separately (losing removal mass twice). Fixed by clamping
   the total once and routing overflow between layers. The control now ties.
2. **"DE: the problem is real and severe."** The v0.1 DE benchmark was broken.
   It renormalised condition B after up-regulating the true-DE genes, and
   compared plain CPM, so naive DE called 58% of NULL genes significant. With
   median-of-ratios normalisation the benchmark is calibrated (7% null FP), and
   naive subtract-then-test has only 8% FP on ambient-trap genes *when given the
   true ambient load*. The remaining risk is ambient **estimation** error, which
   this benchmark does not yet model. The gate now prints its own calibration
   line and declares itself invalid if that line fails.
3. **"Default rho estimator overestimates ~2x."** It was 1.5-4x, and
   catastrophically worse under realistic heterogeneity. Replaced (see above).

## Known limitations (stated, not hidden)

- **Junk in the empty droplets makes rho collapse (detected, NOT fixed).** If
  the empties contain material no cell contains (debris, cell types filtered
  out of the cell set), those genes sit far below the expected soup level and
  the estimator collapses toward 0. In testing, 10% junk took a true rho of 0.2
  to 0.005. A robust fix is an open problem: a two-sided variant resists junk
  but is biased up on clean data. Decant runs both and warns
  (`soup_mismatch`, "rho is probably UNDERestimated") when they disagree
  sharply. On the gate grid it has zero false alarms across the 12 clean
  settings. It caught 20 of 22 per-seed junk runs where rho was underestimated
  by more than 0.05. The 2 misses (rho 0.1, 30% junk) are the worst kind: heavy
  junk collapses the two-sided variant too, so there is no disagreement to see.
  Treat a silent result as "not flagged", not "clean".
- **Weak identifiability when few genes are cleanly absent.** When cell types
  share most of their expression program, rho is not identifiable from counts
  alone and the estimator biases UP (~+10% under stress). Clusters with
  almost no absent genes borrow the evidence-weighted rho of the others and are
  flagged `low_power`. Per-cell resolution is limited (r = 0.5-0.8 under
  stress).
- **Over-clustering biases rho up.** Very small clusters lack the power to
  reject weakly-expressed genes. Prefer a moderate number of reasonably large
  clusters.
- **The DecontX-like comparator is not a fair stand-in for DecontX.** It
  removes only ~8-10% of contamination here, which says more about the
  reimplementation than about DecontX. There is no SoupX reimplementation at
  all. Nothing here is a claim against either tool.
- The lysis diagnostic was validated with an oracle basis; real k-means
  clustering degrades it.
- **All results are on simulated data.** Simulations can be rigged. These are
  sufficiency tests and negative controls, not evidence against CellBender. The
  stress simulator exists because the original one made rho estimation too
  easy.

## The non-negotiable next step

Nothing here is a real-world claim until it is re-run on experimental ground
truth: species-mixing (human/mouse, e.g. 10x hgmm) and genotype-mixing
(demuxlet/souporcell-style) datasets, where cross-species or wrong-genotype
reads measure contamination directly. They should be scored with this same
metric suite, with the official SoupX, DecontX (celda), and CellBender as
comparators. `read_10x_counts()` and `split_droplets()` now get such data in;
the comparator wrappers are next.

## Development

```r
devtools::test()            # unit tests: guarantees, sparse==dense, calibration, I/O
devtools::check()           # R CMD check (clean as of v0.2)
Rscript run_all_gates.R     # the scientific scorecard
```

CI runs `R CMD check` on Linux/macOS/Windows and the gate scorecard on every PR,
and fails if the rho gate stops passing.

## License

MIT. See LICENSE.
