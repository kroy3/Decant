# Decant against existing tools

What has been measured, against what, and what has not. Reproduce with
`validation/benchmark_all.R` (six tools, real PBMC), `validation/benchmark_soupx.R`
and `Rscript validation/de_comparators.R 1 2 3 4`.

Six tools have been run on one real dataset (section 1a): Decant, SoupX
1.6.2, DecontX 1.10.0, CellBender 0.4.0, scAR (scvi-tools 1.4.2) and AmbiDose
0.5.15, all official releases installed from PyPI or from their GitHub
source. FastCAR and CellClear have not been run. Sections 2 and 3 compare
against SoupX only.

## 1a. Six tools on real data: 10x PBMC 4k (no ground truth)

2,150 annotated cells, 33,694 genes, soup from 2.7M empty-droplet UMIs (the
data bundled with SoupX). Every tool gets the same cells and clusters.
"Foreign" is a marker set counted in cell types that should not express it;
"native" is the same set in the type that does.

| | Decant | SoupX | DecontX | CellBender (CPU budget) | AmbiDose | scAR |
|---|---|---|---|---|---|---|
| rho, median | 0.063 | 0.061 | 0.078 | not reported | 0.185 | not reported |
| UMIs removed | 7.9% | 6.1% | 11.6% | 3.0% | 21.3% | 43.8% |
| immunoglobulin: foreign removed | 86% | 65% | 100% | 45% | 82% | 89% |
| monocyte: foreign removed | 94% | 93% | 100% | 97% | 68% | 100% |
| NK: foreign removed | 83% | 65% | 97% | 47% | 80% | 92% |
| immunoglobulin: native removed | 1.3% | 1.1% | 1.6% | 0.0% | 0.4% | 26.3% |
| monocyte: native removed | 4.2% | 2.2% | 5.9% | 0.6% | 5.8% | 19.7% |
| NK: native removed | 2.0% | 2.0% | 3.2% | 0.4% | 1.1% | 17.7% |
| counts added (UMIs) | 0 | 0 | 0 | 0 | 0 | 727,390 |
| run time, 4 CPU cores | 5 s | 14 s | 32 s | 55 min | 73 s | 15 min |

Reading it:

- **DecontX removes the most foreign signal** (97 to 100%), ahead of Decant
  (83 to 94%) on all three sets. It also removes more in total and loses
  somewhat more native signal.
- **CellBender, as run here, is the most conservative**: almost no native
  loss, but only 45 to 47% of two of the three foreign sets removed.
- **Decant, SoupX and DecontX form a middle group**: 6 to 12% of UMIs
  removed, native loss under 6%, no counts added. Within it Decant is the
  fastest, ahead of SoupX on foreign removal and behind DecontX.
- **AmbiDose estimates three times the contamination of the others** (0.185)
  and removes 21% of UMIs, yet removes less monocyte signal than Decant. It
  flagged 65% of its own per-cell estimates as "interpret with caution".
- **scAR removes 44% of all UMIs, loses 18 to 26% of native marker signal,
  and adds 727,390 counts** across 666,033 entries where the cell had fewer.
  Its output is one random sample from its model, which accounts for part of
  the added counts.

Limits that matter for this table:

- There is no ground truth. Which total removal is right (3%, 8%, 12%) cannot
  be decided from marker sets; "removes more foreign signal" is only better if
  those counts are ambient. Under Decant's own estimate the ambient model
  predicts 92 to 106% of them, which is circular as evidence for Decant.
- **The empty droplets are reconstructed.** The SoupX data ships the summed
  empty-droplet counts, not the droplets. 40,000 droplets were made by
  assigning each soup UMI to a droplet at random, so the pooled soup is real
  but the droplets are not. This is most likely to handicap CellBender, which
  models droplets directly; AmbiDose also uses them.
- **CellBender was not run at its defaults.** With 150 epochs and cells plus
  10,000 droplets it completed 28 epochs in two hours on 4 CPU cores. The row
  above is 100 epochs with cells plus 3,000 droplets. It needs a GPU to be
  run as intended.
- One dataset, one tissue, clusters taken from the SoupX annotation.

## 1b. Decant and SoupX on the real soup profile (no reconstructed droplets)

2,150 annotated cells, 33,694 genes, soup from 2.7M empty-droplet UMIs (the
data bundled with SoupX). Both tools were given the same cells, soup and
clusters. "Foreign" is a marker set counted in cell types that should not
express it; "native" is the same set in the type that does.

| | Decant | SoupX 1.6.2 |
|---|---|---|
| rho, median / mean | 0.063 / 0.076 (per cell) | 0.061 (one value for all cells) |
| UMIs removed | 7.9% | 6.1% |
| immunoglobulin (IGKC, IGHM, IGLC2): foreign removed | 86% | 65% |
| immunoglobulin: native removed | 1.3% | 1.1% |
| monocyte (LYZ, S100A8, S100A9): foreign removed | 94% | 93% |
| monocyte: native removed | 4.2% | 2.2% |
| NK (GNLY, NKG7): foreign removed | 83% | 65% |
| NK: native removed | 2.0% | 2.0% |
| counts added to any entry | 0 | 0 |
| run time | 7 s | 13 s |

Under Decant's rho and soup, the ambient model predicts 94%, 106% and 92% of
the foreign counts in the three sets, so removing most of them is the
expected behaviour, not over-correction. Decant removes more of the foreign
signal for two of three sets and about the same for the third. It also
removes twice as much native monocyte signal (4.2% against 2.2%), which
follows from its higher rho in the monocyte clusters; whether that is
over-correction cannot be decided without ground truth.

This is one dataset from one tissue. Before the fix in this release Decant
estimated rho at 0.22 on it, and would have lost this comparison badly.

## 2. Simulated data with ground truth

`simulate_experiment()`, easy and stress settings, k-means clusters given to
both tools, mean of 3 seeds. L1 error is the distance of the corrected counts
from the true counts, as a share of true counts (lower is better).

| setting | true rho | Decant rho | SoupX rho | Decant L1 | SoupX L1 | SoupX runs that failed |
|---|---|---|---|---|---|---|
| easy | 0.00 | 0.000 | 0.016 | 0.000 | 0.016 | 0 of 3 |
| easy | 0.05 | 0.052 | 0.045 | 0.023 | 0.032 | 0 of 3 |
| easy | 0.10 | 0.102 | 0.084 | 0.049 | 0.055 | 0 of 3 |
| easy | 0.20 | 0.203 | not reported | 0.085 | not reported | 1 of 3 |
| stress | 0.00 | 0.004 | 0.011 | 0.003 | 0.011 | 0 of 3 |
| stress | 0.05 | 0.062 | 0.075 | 0.054 | 0.068 | 0 of 3 |
| stress | 0.10 | 0.110 | not reported | 0.082 | not reported | 1 of 3 |
| stress | 0.20 | 0.220 | not reported | 0.139 | not reported | 2 of 3 |

**This table is weak evidence and should not be quoted as Decant beating
SoupX.** In every run with contamination, SoupX's default marker search found
no markers; its own fallback (`soupQuantile = 0.5`) was applied once, and four
runs still failed. The simulated profiles do not have the marker structure
SoupX is built to find, and the simulator was written alongside Decant. What
the table does show is that Decant recovers rho on data generated from its
own model, including a true rho of zero.

## 3. Differential expression (simulated case/control)

See `results/de_comparators.csv`. Six control and six disease samples;
disease samples carry more ambient RNA; 80 genes that are disease genes only
in an abundant cell type can leak into the rare type being tested.

Mean of 4 seeds. Every corrector gets the true cell-type labels.

| case | method | leaked genes called | null FP | power | FDP |
|---|---|---|---|---|---|
| more soup in disease | naive edgeR | 52.8 / 80 | 0.223 | 0.91 | 0.82 |
| | `ambient_de` | 0 / 80 | 0.048 | 0.92 | 0.13 |
| | Decant `bayes` (default) + edgeR | 0 / 80 | 0.045 | 0.89 | 0.11 |
| | Decant `redistribute` + edgeR | 6.0 / 80 | 0.081 | 0.90 | 0.31 |
| | SoupX + edgeR | 1.5 / 80 | 0.043 | 0.85 | 0.13 |
| equal soup, shifted composition | naive edgeR | 12.5 / 80 | 0.060 | 0.92 | 0.35 |
| | `ambient_de` | 0 / 80 | 0.051 | 0.94 | 0.11 |
| | Decant `bayes` (default) + edgeR | 0 / 80 | 0.051 | 0.94 | 0.10 |
| | Decant `redistribute` + edgeR | 0.5 / 80 | 0.060 | 0.94 | 0.12 |
| | SoupX + edgeR | 0.25 / 80 | 0.042 | 0.89 | 0.09 |
| genes truly DE in both types | naive edgeR | 30.8 / 60 | 0.227 | 0.92 | 0.71 |
| | `ambient_de` | 0 / 60 | 0.057 | 0.97 | 0.11 |
| | Decant `bayes` (default) + edgeR | 0 / 60 | 0.045 | 0.95 | 0.09 |
| | Decant `redistribute` + edgeR | 9.0 / 60 | 0.114 | 0.94 | 0.43 |
| | SoupX + edgeR | 2.3 / 60 | 0.051 | 0.90 | 0.14 |

**Correcting with Decant's default rule and then running edgeR matches
`ambient_de`**: no leaked gene, calibrated, slightly lower false discovery
proportion, 0 to 3 points less power. SoupX followed by edgeR is close
behind, with 4 to 7 points less power and up to 2.3 leaked genes. On this
simulation the gain from modelling the ambient term inside the test, over a
good correction, is within noise. SoupX needed its `soupQuantile = 0.5`
fallback on every sample.

## 4. What each tool does

| | Decant | SoupX | DecontX | CellBender | scAR | AmbiDose |
|---|---|---|---|---|---|---|
| uses empty droplets | yes | yes | optional | yes (all droplets) | yes | yes |
| contamination estimate | per cell, automatic | one per channel, automatic or manual | per cell | per cell | per cell | per cell |
| can add counts | no (checked on every run) | no (0 here) | no (0 here) | no (0 here) | yes (727k UMIs here) | no (0 here) |
| needs GPU | no | no | no | in practice | in practice | no |
| contamination-aware DE test | yes (`ambient_de`) | no | no | no | no | no |
| run here | yes | yes | yes | yes, reduced configuration | yes | yes |

## 5. What the recent literature supports, and what it does not

These points come from search summaries of the papers. The full texts were
not reachable from where this was written, so details should be checked
against the papers before being cited.

- **Count integrity matters.** A 2026 preprint benchmarking six tools on
  mixtures, tissue and well-plate data reports that scAR and CellClear
  restructure count matrices, including creating counts, and argues that
  count integrity should be a primary criterion
  (bioRxiv 10.64898/2026.04.08.717130). This supports Decant's hard
  constraint that a corrected count never exceeds the observed one. It does
  not distinguish Decant from SoupX, DecontX or CellBender, which the same
  preprint does not accuse of adding counts.
- **No tool wins everywhere.** A 2026 preprint concludes that no single
  method is best across datasets and metrics, with CellBender, DecontX and
  SoupX generally performing well (bioRxiv 10.64898/2026.01.13.699237). A 2026
  Nature Communications benchmark reports under-correction of highly
  contaminating genes by some tools and over-correction of weakly
  contaminating genes by others (10.1038/s41467-026-77458-4). Neither
  includes Decant. Both imply that a claim to outperform needs several
  datasets with ground truth, which Decant does not yet have.
- **Ambient RNA distorts differential expression.** FastCAR was written for
  this problem and notes that contamination can differ between patients and
  controls (Berg et al. 2023, BMC Genomics, 10.1186/s12864-023-09822-3). A
  2025 PLOS ONE study reports many false-positive DE genes from uncorrected
  ambient contamination (10.1371/journal.pone.0332440). Caglayan et al. (2022,
  Neuron) showed ambient RNA producing misread cell types in brain
  single-nucleus data. These support the problem `ambient_de()` addresses.
  They all correct counts and then test; none models the ambient term inside
  the test.
- **Large-cohort DE.** Li et al. (2022, Genome Biology,
  10.1186/s13059-022-02648-4) showed NB-based tests losing false-positive
  control with outlier samples. `ambient_de()` had that failure and now
  down-weights outliers.
- **A close competitor exists.** AmbiDose (bioRxiv 10.64898/2026.09.27.754849,
  a Python package) describes per-cell ambient doses from empty droplets,
  protection of lineage genes, and never inflating counts. That overlaps
  with Decant's correction step almost entirely, so "per-cell, conservative,
  never adds counts" is not a novelty claim Decant can make.

## 6. Where Decant is ahead, level, behind, or unknown

| claim | status | evidence |
|---|---|---|
| Removes more foreign marker signal than SoupX at similar rho | ahead on 2 of 3 marker sets, one real dataset | sections 1a, 1b |
| Preserves native signal as well as SoupX | level on 2 of 3 sets, behind on monocytes | sections 1a, 1b |
| Never adds counts | level with SoupX, DecontX, CellBender, AmbiDose; ahead of scAR (727k UMIs added) | section 1a |
| Per-cell contamination estimate | ahead of SoupX (one value per channel); level with DecontX, CellBender, AmbiDose by design; accuracy of the per-cell values untested on real data | section 4 |
| DE without leakage, power kept | level with Decant-correct-then-edgeR; slightly ahead of SoupX-then-edgeR (4 to 7 points of power); simulated only | section 3 |
| DE with covariates in large cohorts | ahead of subtract-then-rank-test in simulation | `experiments/large_cohort_calibration/` |
| Foreign marker removal against DecontX | behind on all three sets | section 1a |
| Foreign marker removal against CellBender, AmbiDose | ahead on 2 of 3 and 3 of 3 sets; CellBender in a reduced configuration | section 1a |
| Native signal kept, against scAR | far ahead (1 to 4% lost against 18 to 26%) | section 1a |
| Speed | fastest of the six on this dataset | section 1a |
| Accuracy on data with experimental ground truth | unknown, not run | species and genotype mixture scripts exist, data unreachable |
| Robust to empties that over-represent some genes | behind: not detected | README, Known limitations |

## 7. Novelty

What appears to be new, as far as the searches above show:

1. **The ambient term inside the DE model.** `ambient_de()` fits
   `mu = N exp(X beta) + E` with `E` the per-sample ambient expectation and
   its uncertainty in the variance, instead of correcting counts and testing
   afterwards. No tool found in the searches does this. It is new, but
   section 3 shows it has not yet been shown to be *better* than correcting
   well and then testing.
2. **Falsification gates in CI.** Every default has a gate that can fail the
   build, and results that went against the package are kept in the README.

What is not new: absent-gene estimation of rho (SoupX), per-cell
contamination (DecontX, CellBender, AmbiDose), mass-conserving correction
(SoupX and others), posterior-based removal (DecontX, CellBender).
