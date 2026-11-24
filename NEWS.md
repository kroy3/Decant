# Decant (development version)

* **Identifiability diagnostics and a classed warning.** `decant()` returns
  `$identifiability`: per cluster, the weighted cosine and Gram condition
  number between its profile and the soup, their Jensen-Shannon divergence,
  the Fisher information about rho per UMI and how many genes carry it, and
  rho's standard error and 95% interval width (from the absent-gene
  estimator's negative binomial variance). A cluster is flagged `weak` when
  its estimate rests on thin evidence (interval wider than 0.05, almost no
  absent genes, or an implausible own estimate), with the reason and what
  was done; the warning has class `decant_weak_identifiability` (and
  `decant_unidentifiable` when nothing can be estimated) and carries the
  table. Closeness to the soup alone (`collinear`: cosine > 0.9) is
  reported but does not flag a cluster by default: in simulation such a
  cluster with enough absent genes was estimated within 0.005 of the truth.
  `weak_rule = "any"` flags it too. Known limitation, found with these
  diagnostics: a borrowed cluster's per-cell values still use its own
  absent genes and can run high (simulation: cluster mean 0.26, borrowed
  0.101, truth 0.103); not changed here.

* **Large raw matrices without loading the empty droplets.**
  `read_10x_streaming()` reads a raw 10x `matrix.mtx(.gz)` in two chunked
  passes and returns the cells as a sparse matrix plus a `soup_summary()`:
  per-gene sums, sums of squares and detection counts of the empty droplets
  (O(genes)). `decant()`, `ambient_global()` and `pool_soup()` accept such a
  summary wherever they took an empty-droplet matrix, with identical
  results. `soup_dispersion()` gives per-gene negative binomial dispersions
  of the empties from those moments. On hgmm 10k (6.8M barcodes, 1.9M empty
  droplets) the read peaks at 1.8 GB instead of 4.3 GB, and the soup
  summary takes 12 MB. Measured on that channel, the empties were never the
  main memory cost: the cells (112k genes x 11k cells, 55M non-zeros) are.

* **Correction is about 3x faster and uses half the memory.** The
  cell-aware rules run over column blocks, group the absence likelihood by
  distinct count pairs, and compute each tail probability once. On hgmm 10k
  `decant()` takes 55 s instead of 193 s and the correction step peaks at
  2.2 GB instead of 4.7 GB above its input; output is unchanged (largest
  difference 2.6e-10).

* **`bayes_absent` no longer fits degenerate absence calls in tiny
  clusters.** The call now needs at least 50 genes with 5 or more expected
  ambient counts in the cluster (and in the empties), and a cluster whose
  share of genes called absent exceeds every other cluster's by more than
  0.10 keeps the `bayes_cell` rule. In simulation a 1-cell cluster kept
  0.37 of its native counts before and 0.83 now (`bayes_cell`: 0.83). The
  9-cell brain cluster of 1i is the case this targets; it has not been
  re-run.

* **The bound 0 <= corrected <= raw is now asserted, not clipped.** Every
  rule satisfies it by construction (checked over 64.5M entries of the gate
  simulations before the change); a violation is now an error of class
  `decant_invariant_violation`. New property tests cover random negative
  binomial data across contamination levels and degenerate clusters, and
  doublet-level counts.

* **Real-tissue test of the default (pre-registered): no gain.** In six brain
  libraries (GSE272082), scored against genes foreign to each cell type in an
  independent C1 reference, `bayes_absent` removed the same share of foreign
  signal from microglia as `bayes_cell` (0.092 against 0.091). Its gene-level
  call applied in 1 of 42 clusters, so its large gain is shown only on
  species mixtures. Oligodendrocytes and astrocytes lost 2.7 points more
  foreign signal, from the cap rule. The one cluster where the call applied
  had 9 cells and 17,543 genes called absent, a degenerate fit; avoid the
  default on tiny clusters (`validation/BENCHMARK.md` 1i).

* **The default correction is now `bayes_absent`** (was `bayes_cell`).
  It removes the most ambient RNA of every tool run here: 0.891 on two
  unused 10x mixtures, against `bayes_cell` 0.531 and DecontX 0.655.
  It missed its pre-registered retain criterion by 0.02 points against
  DecontX and was made the default by maintainer decision, recorded in
  `validation/BENCHMARK.md` 1h. Costs: 0.3 to 0.6 points of native signal
  against `bayes_cell` on most mixtures, 2.5 on the hardest (hgmm 10k). On
  PBMC 4k its output matches `bayes_cell`. Results from earlier versions
  change; pass `correction = "bayes_cell"` to reproduce them.
  `gate_correction()` now includes `bayes_absent` and fails if the default
  is not the rule with the lowest error (now 0.0628 against 0.0631, a narrow
  margin).

* **`bayes_absent` now uses `cap_test = "hybrid"` and no longer removes
  doublet-level counts (PBMC identical to `bayes_cell`).** Judged on two
  unused 10x mixtures under a pre-registration: removal 0.891 against
  DecontX 0.655, retain 0.99183 against 0.99201. The retain criterion fails
  by 0.02 points, so it stays opt-in; the default is unchanged
  (`validation/BENCHMARK.md` 1h).

* **New opt-in correction `correction = "bayes_absent"`, judged on fresh
  data: removes the most, fails two pre-stated criteria.** On three 10x
  mixtures never used for design it removes 0.885 of ambient signal (DecontX
  0.761, SoupX 0.614, default `bayes_cell` 0.566). But its mean retain is
  0.9837 vs DecontX 0.9847, and on PBMC 4k it removes more counts too high
  to be ambient (likely doublets), so it is not the default. The negative
  control (0.998), native markers and the simulation gate pass
  (`validation/BENCHMARK.md` 1g, `validation/check_bayes_absent.R`).

* **New, experimental: `chirp()`, ambient-aware ligand-receptor calls.** A
  ligand made by the cell type that dominates the soup is detected in every
  cluster, so a detection rule (>= 10% of cells) reports every cluster as a
  sender. `chirp()` calls a gene in a cluster only when its count is
  significantly above `margin` times the ambient count expected there, and
  scores pairs on the excess. Simulated ground truth only (2 seeds, rho
  0.05-0.3): the detection rule called 600-950 absent pairs, `chirp()` none,
  and it found every pair expressed at >= 20 per 10k. It misses most pairs at
  5-20 per 10k. With rho at 2/3 of its value it still made no false calls; at
  1/2 it made up to 100. Since the contamination estimate is known to run
  low on some real data (1f below), check `ambient_share_ligand` before
  trusting a call. Not yet tested on real data. Correcting counts with
  `bayes_cell` and then applying the detection rule left 70-388 of those
  false calls on the same simulations; `subtract` left all of them, since
  its fractional counts still count as detected.

* **Second pre-registered real-data DE test: null.** Paired brain regions
  (GSE272082, hippocampus or entorhinal cortex vs prefrontal cortex, within
  donors): standard edgeR called no gene foreign to microglia, so there was
  nothing for `ambient_de` to remove (`validation/BENCHMARK.md` 3d). In
  microglia EC vs PFC, `ambient_de` had a one-gene hit in 8 of 15 label
  swaps, above the pre-stated limit.
* `ambient_de()`: robust weights are kept finite and positive (a NaN or
  infinite Pearson residual made edgeR stop on real data).

* **Contamination estimate on the 10k mixture: diagnosed, not fixed.** The
  absent-gene fixed point collapses there (1 to 6 absent genes kept of
  5,000+) because its overdispersion term starts at zero and its low start
  settles on the soup's most abundant genes. A fallback start that fixed the
  10k mixture failed on held-out data (hgmm 5k median 43% low, Spearman
  0.49) and turned the negative control's correct near-zero into NA, so the
  estimator is unchanged (`validation/BENCHMARK.md` 1f).
  `validate_species_mixture.R` now records a method that declines (rho NA).

* **Ambient removal: held-out results, and a change that did not pass.** On
  held-out human/mouse mixtures (hgmm 1k, 5k, 10k) Decant removes 0.444 of
  ambient signal against DecontX 0.637 and SoupX 0.679, while keeping the
  most (retain 0.992). On the 20k mixture the loss was traced to
  `bayes_cell`'s cap, which keeps every single count of a gene with a tiny
  ambient expectation. Correcting that test for the number of cells
  (`cap_test = "cluster"`, opt-in) raised held-out removal to 0.519 but
  dropped retain to 0.974 on the 10k mixture, so it is not the default
  (`validation/BENCHMARK.md` 1e). The larger part of the held-out gap comes
  from contamination estimates that are too low on those mixtures.

* **`ambient_de()` guards against implausible ambient fractions.** A single
  cell type's pseudobulk cannot always identify its own ambient fraction. On
  GSE273597 (pancreas) six donors' ductal values came out at 0.86 to 0.96 and
  shuffled labellings gave up to 2,300 false hits. A sample above the new
  `max_rho` (0.5) now takes the median of the plausible samples, with their
  spread added to its uncertainty, and `ambient_de` warns
  (`attr(res, "samples")$ambient_frac_imputed`). Ductal shuffles: 4 of 24 with
  a hit, at most 1 gene. GSE190015 and GSE272082 unchanged. The pancreas data
  found the problem, so they do not validate the guard elsewhere.
* **Pre-registered real-data test (`validation/BENCHMARK.md` 3c): null.** On
  T1D vs non-diabetic pancreas, naive edgeR found no ambient-driven hits to
  remove. `ambient_de` is now described as experimental in the README.

* **`ambient_de()` is now calibrated on real data.** It fitted each gene with
  edgeR's tagwise NB dispersion and then ran a quasi-likelihood F-test, which
  estimates the gene's extra variance again: the variance was counted twice
  and the test was liberal. It now uses the trended dispersion, as edgeR's
  quasi-likelihood pipeline does (`dispersion = "trended"`; `"tagwise"`
  restores the old behaviour). With condition labels shuffled between donors
  on two real datasets, the old default reported genes in 82 of 96 tests (up
  to 574 genes), the new one in 9 (at most 9 genes), naive edgeR in 4.
  GSE190015 (nine controls) was held out from the diagnosis. The simulation
  gates still pass with unchanged power. Re-run analyses made with earlier
  versions. `joint_rho`, first suspected, was not the cause.
  New: `validation/validate_de_calibration.R`, `validation/BENCHMARK.md` 3b.

* **Contamination estimate no longer halved on real data.** Genes a cluster
  does not express scatter around rho by about 25-30% from gene to gene, far
  more than Poisson noise. `estimate_rho_cluster()` tested them against a
  Poisson, rejected the absent genes above the current rho and kept those
  below, and settled near the bottom of their spread. It now uses a negative
  binomial whose extra variance is estimated from the lower tail of the
  absent genes (`overdispersed = TRUE`). On a 20k human/mouse mixture the
  median rho went from 0.013 to 0.027 (truth 0.029) and ambient removal from
  0.415 to 0.616. On further mixtures (`validation/BENCHMARK.md`, 1d) it
  has the lowest error of Decant, SoupX and DecontX on 1k and 5k cells, but
  still underestimates on 10k cells, where DecontX ranks cells better.
* **rho = 0 trap fixed.** With sparse counts the starting value could be 0,
  which the fixed point cannot leave; on a 1k-cell mixture 7 of 8 clusters
  returned exactly 0. The start is now at least `floor_rho`.
* **Unidentifiable contamination is reported.** When the soup resembles the
  cells (a single cell type), no gene is absent and rho cannot be estimated;
  the estimator returned 0.5-0.94 for a true 0.004 without a warning.
  Clusters above `max_rho = 0.5` are flagged `implausible` and borrow rho from
  the others; when most cells are in such clusters `decant()` warns that rho
  is not identifiable and should be supplied (`$rho_unidentifiable`).
* `$rho_diagnostics` gains `dispersion` and `implausible` columns.
* The placeholder-name note is removed from DESCRIPTION: no CRAN or
  Bioconductor package is called Decant.
* New user guide in 13 chapters (`wiki/`, also published as website
  articles). Every example runs on the bundled simulators, and
  `scripts/check_user_guide.R` checks that each shown output is current.
* Scripts renamed for clarity: repository-level scripts moved to `scripts/`
  (`run_falsification_gates.R`, `compare_correction_methods.R`,
  `demo_structured_soup.R`, `stress_test_empty_mismatch.R`); validation
  scripts follow `benchmark_*` (Decant against other tools) and `validate_*`
  (Decant against ground truth); result files renamed to match.
* New correction rule `"bayes_cell"`, now the default in `decant()`. `"bayes"`
  scaled every cell of a cluster by one ratio, so a cell whose count of a
  gene was far too high to be ambient (a doublet) lost most of it and
  inflated the cluster's native profile. `"bayes_cell"` removes only the
  expected ambient from such entries and builds the profile from the rest.
  It has the lowest error in `gate_correction()` (0.065, against 0.071 to
  0.075). On real PBMC data it removes 88 to 95% of ambient-compatible
  foreign marker counts (previously 84 to 94%) and 14 to 37% of the
  incompatible ones (previously 82 to 98%), at the same native loss.
* `validation/benchmark_pbmc_six_tools.R` splits foreign marker counts into those
  ambient RNA could and could not have produced. DecontX removes both
  (96 to 100% of the incompatible ones).
* New `validation/benchmark_pbmc_six_tools.R` and `validation/py/run_python_tools.py`:
  Decant against SoupX, DecontX, CellBender, scAR and AmbiDose on real PBMC
  data. DecontX removes more foreign marker signal than Decant; scAR adds
  counts and removes 44% of UMIs; Decant is the fastest. See
  `validation/BENCHMARK.md`.
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
* New `validation/validate_pbmc_markers.R`, the first validation run on real data.
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
* `scripts/run_falsification_gates.R` exits non-zero when a required gate fails; CI runs each
  gate once.
* `experiments/large_cohort_calibration/`: sweep over 6 to 200 samples per
  group, with its limits.
* New `validation/benchmark_differential_expression.R`: contamination-aware DE against
  correct-then-test with the official SoupX, DecontX, FastCAR and CellBender
  (each run only if installed). Wrappers live in `validation/comparators.R`.
* README: the claim that correcting counts first does not stop DE leakage is
  narrowed. It holds for Decant's default `subtract` rule; official SoupX
  followed by edgeR removes almost all leakage on the simulator.
* New vignette, "Contamination-aware differential expression".
* pkgdown site configuration and deployment workflow.
* `DESCRIPTION`: `URL`, `BugReports`, `biocViews`, `VignetteBuilder`.
