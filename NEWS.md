# Decant (development version)

* `ambient_de()` gains `robust = TRUE` (default): outlier samples are
  down-weighted per gene by Huber weights and the dispersion is re-estimated
  under them. Without it, cohorts of 50 or more samples per group with a few
  extreme counts gave 16 to 41% false positives on null genes (the failure
  Li et al. 2022 report for NB tests). `DecantDE()` inherits the default.
* New `simulate_ambient_cohort()` and `gate_ambient_de_cohort()`: large
  cohorts with outlier samples or heavy-tailed noise. The gate runs in CI.
* `experiments/large_cohort_calibration/`: sweep over 6 to 100 samples per
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
