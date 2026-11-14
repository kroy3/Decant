# Decant (development version)

* **Bug fix, found on real data.** `estimate_rho_cluster()` overestimated rho
  several-fold on real 10x data (median 0.22 on PBMC 4k, where SoupX gives
  0.06 to 0.077). Low-count genes passed the absent-gene test for lack of
  power. The absent set now requires an expected ambient count of at least
  `min_expected = 5`; the PBMC estimate is 0.063 (median).
* `decant()` now defaults to `correction = "bayes"`: on the PBMC data it
  removed 83 to 94% of marker counts from cell types that do not express
  them, against 19 to 57% for `subtract`.
* Genes present in the empties but absent from every cell are dropped from
  the soup (`floor_rho`), and reported in the mismatch warning.
* The two-sided soup/cell mismatch tripwire is retired: it raised a false
  alarm on the PBMC data. `soup_mismatch` now means more than 5% of the soup
  was in soup-only genes. Empties that over-represent genes the cells do
  contain are no longer flagged; this is a known limitation.
* A cluster that cannot identify its own rho borrows the weighted median of
  the others, not a count-weighted mean.
* New `validation/pbmc_soupx.R`, the first validation run on real data.
* New `validation/benchmark_soupx.R` and `validation/BENCHMARK.md`: Decant
  against the official SoupX on real PBMC and simulated data, with the
  literature and what is and is not new. With the `bayes` default,
  correcting with Decant and then running edgeR matches `ambient_de` on the
  simulated 6 vs 6 study.
* `ambient_de()` gains `robust = TRUE` (default): outlier samples are
  down-weighted per gene by Huber weights, the dispersion is re-estimated
  under them, and the residual degrees of freedom are reduced to match.
  Without it, cohorts of 50 or more samples per group with a few extreme
  counts gave 16 to 41% false positives on null genes (the failure Li et al.
  2022 report for NB tests). `DecantDE()` inherits the default.
* `ambient_de()` estimates the ambient term more accurately:
  - `joint_rho = TRUE`: each sample's ambient fraction is refined on genes
    that look absent across all samples pooled, which removes an upward bias
    of about one point; its sampling variance now enters the model.
  - `soup_design = TRUE`, and `soup_eb(design = )`: soups are shrunk toward
    the design-fitted soup, not a condition-blind mean.
* `ambient_de(max_ambient_frac = 0.75)`: genes whose counts are mostly soup
  are not tested and are listed in attribute `"ambient_dominated"`.
  Previously most such genes were dropped only as a side effect of the
  over-estimated ambient fraction.
* `ambient_de(lfc = )`: test against a minimum fold change (as limma's
  TREAT), for cohorts large enough that sub-1% biases become significant.
* New `simulate_ambient_cohort()` and `gate_ambient_de_cohort()`: large
  cohorts with outlier samples, heavy-tailed noise, or a confounding batch.
* `run_all_gates.R` exits non-zero when a required gate fails; CI runs each
  gate once.
* `experiments/large_cohort_calibration/`: sweep over 6 to 200 samples per
  group, with its limits.
* New `validation/de_comparators.R`: contamination-aware DE against
  correct-then-test with the official SoupX, DecontX, FastCAR and CellBender
  (each run only if installed). Wrappers live in `validation/comparators.R`.
* README: the claim that correcting counts first does not stop DE leakage is
  narrowed. It holds for Decant's default `subtract` rule; official SoupX
  followed by edgeR removes almost all leakage on the simulator.
* New vignette, "Contamination-aware differential expression".
* pkgdown site configuration and deployment workflow.
* `DESCRIPTION`: `URL`, `BugReports`, `biocViews`, `VignetteBuilder`.
