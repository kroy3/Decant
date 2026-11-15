# Decant against existing tools

What has been measured, against what, and what has not. Reproduce with
`validation/benchmark_pbmc_six_tools.R` (six tools, real PBMC), `validation/benchmark_soupx.R`
and `Rscript validation/benchmark_differential_expression.R 1 2 3 4`.

Six tools have been run on one real dataset (section 1a): Decant, SoupX
1.6.2, DecontX 1.10.0, CellBender 0.4.0, scAR (scvi-tools 1.4.2) and AmbiDose
0.5.15, all official releases installed from PyPI or from their GitHub
source. FastCAR and CellClear have not been run. Section 1c scores Decant,
SoupX and DecontX on two datasets with experimental ground truth from
Cargnelli et al. (2026). Sections 2 and 3 compare against SoupX
only.

## 1a. Six tools on real data: 10x PBMC 4k (no ground truth)

2,150 annotated cells, 33,694 genes, soup from 2.7M empty-droplet UMIs (the
data bundled with SoupX). Every tool gets the same cells and clusters.
"Foreign" is a marker set counted in cell types that should not express it;
"native" is the same set in the type that does.

| | Decant | SoupX | DecontX | CellBender (CPU budget) | AmbiDose | scAR |
|---|---|---|---|---|---|---|
| rho, median | 0.063 | 0.061 | 0.078 | not reported | 0.185 | not reported |
| UMIs removed | 7.9% | 6.1% | 11.6% | 3.0% | 21.3% | 43.8% |
| immunoglobulin: foreign removed | 70% | 65% | 100% | 45% | 82% | 89% |
| monocyte: foreign removed | 93% | 93% | 100% | 97% | 68% | 100% |
| NK: foreign removed | 79% | 65% | 97% | 47% | 80% | 92% |
| immunoglobulin: native removed | 1.3% | 1.1% | 1.6% | 0.0% | 0.4% | 26.3% |
| monocyte: native removed | 4.1% | 2.2% | 5.9% | 0.6% | 5.8% | 19.7% |
| NK: native removed | 2.0% | 2.0% | 3.2% | 0.4% | 1.1% | 17.7% |
| counts added (UMIs) | 0 | 0 | 0 | 0 | 0 | 727,390 |
| run time, 4 CPU cores | 11 s | 12 s | 34 s | 55 min | 73 s | 15 min |

### Foreign counts are not all ambient

"Foreign removed" treats every marker count in the wrong cell type as
contamination. It is not. Fifteen non-B cells (14 annotated as monocytes)
hold 397 of the 1,188 foreign immunoglobulin UMIs. They carry 16.7 UMIs per
10,000 of CD79A, CD79B and MS4A1, against 0.7 in other non-B cells and 41.7
in B cells, and have 1.5 times the UMIs of their neighbours: most likely
B-cell doublets. Those molecules were in the droplet with the cell. An
ambient tool that removes them is not removing ambient RNA.

So the foreign counts are split. An entry is **incompatible** with ambient
when its count is beyond the 99.9% Poisson bound of a deliberately generous
ambient expectation (contamination 0.15, about twice any tool's estimate
here, so the bound uses no tool's own rho). The share of foreign UMIs that
is incompatible is 31% for immunoglobulins, 2% for monocyte markers and 18%
for NK markers. A correct tool removes the compatible counts and leaves the
incompatible ones.

| | Decant | Decant, previous `bayes` rule | SoupX | DecontX | CellBender (CPU budget) | AmbiDose | scAR |
|---|---|---|---|---|---|---|---|
| immunoglobulin: compatible removed | 95% | 88% | 88% | 100% | 63% | 75% | 99% |
| monocyte: compatible removed | 95% | 94% | 94% | 100% | 98% | 67% | 100% |
| NK: compatible removed | 88% | 84% | 76% | 97% | 54% | 79% | 100% |
| immunoglobulin: incompatible removed | 14% | 82% | 12% | 100% | 4% | 96% | 66% |
| monocyte: incompatible removed | 16% | 98% | 31% | 100% | 68% | 100% | 82% |
| NK: incompatible removed | 37% | 82% | 14% | 96% | 18% | 82% | 58% |

Reading both tables:

- **DecontX removes everything foreign, including the counts that cannot be
  ambient** (96 to 100% of them). Its lead on "foreign removed" is partly
  that. On the compatible counts it is still ahead of Decant: by 5 points on
  two sets and 9 on the NK set.
- **Decant's previous `bayes` rule had the same fault to a lesser degree**: it
  scaled every cell of a cluster by one ratio and removed 82 to 98% of the
  incompatible counts. The new default, `bayes_cell`, removes 14 to 37% of
  them and more of the compatible counts than before (88 to 95%), at the
  same native loss and the same total removal.
- **SoupX and CellBender also leave the incompatible counts**, and remove
  less of the compatible ones than Decant on two of three sets (SoupX) or on
  two of three (CellBender, which is ahead on the monocyte set).
- **AmbiDose and scAR** remove most incompatible counts; scAR also adds
  727,390 counts and loses 18 to 26% of native marker signal.
- **Decant and SoupX are the fastest** (11 and 12 s), and no tool except scAR
  adds counts.

Whether the remaining 5 to 12 points between Decant and DecontX on the
compatible counts are real ambient RNA that Decant misses, or low native
expression that DecontX removes, cannot be decided here. Decant removes what
its rho and the soup predict; in the clusters where counts remain, the
observed foreign counts exceed that prediction by 1.3 to 2.3 times.

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

## 1c. Ground truth: two datasets from Cargnelli et al. (2026)

Two datasets from "Benchmarking computational decontamination of ambient RNA"
(bioRxiv, doi:10.64898/2026.01.13.699237), scored with the paper's two
metrics. For each gene with input count x > 0 (summed over a sample's cells
of one species) and corrected count y:

- **endogenous retain score** = mean over the cell's own species' genes of max(0, 1 - |y - x| / x). 1 means nothing endogenous was removed.
- **ambient removal score** = mean over the other species' genes of min(1, |x - y| / x). 1 means all ambient counts were removed.

Both are unweighted means over genes. The tools that ran here were Decant,
SoupX 1.6.2 and DecontX 1.10.0 (R 4.3.3, Seurat 5.0.1, DropletUtils 1.22.0).
The other rows are the paper's published numbers, from the paper's own
pipeline. They are context only, not a run here.

### Low-complexity human/mouse mixture (10x 20k HGMM, 3' HT)

Reproduce with `Rscript validation/benchmark_species_mixture.R low_complexity`.
Preprocessing follows the paper: emptyDrops (lower = 500, FDR <= 0.05), cells
with at least 75% of UMIs from one species kept, and Seurat clusters at
resolution 0.1. This gives 17,619 cells (9,030 human, 8,589 mouse) and 6
clusters. The soup is the 2,952,920 non-cell barcodes with 1 to 100 UMIs,
which is SoupX's default range. Every tool gets the same cells, soup and
clusters. The paper's other QC steps are not reproduced.

| | retain | removal | run time |
|---|---|---|---|
| **Decant** | **0.989** | **0.616** | 194 s |
| Decant before the estimator fix (0.4.0) | 0.991 | 0.415 | 136 s |
| SoupX, run here | 0.992 | 0.345 | 216 s |
| DecontX, run here | 0.983 | 0.722 | 134 s |
| *published:* SoupX | 0.988 | 0.422 | |
| *published:* DecontX | 0.980 | 0.624 | |
| *published:* CellBender | 0.996 | 0.729 | |
| *published:* scAR | 0.470 | 0.959 | |
| *published:* CellClear | 0.876 | 0.023 | |
| *published:* FastCAR / scCDC | 1.000 | 0.000 / 0.002 | |

**Decant still removes less ambient RNA than DecontX and CellBender, but the
gap fell from about 30 points to about 10.** It removes 0.616 against
DecontX's 0.722 and the published CellBender's 0.729, while keeping more of
the cells' own counts than DecontX (0.989 against 0.983). The change comes
from the contamination estimate (`validation/benchmark_species_mixture_rho.R`).
Cross-species reads give each cell's contamination directly: in a human
cell, mouse reads can only be ambient.

| | rho, median | rho, mean | Spearman with cross-species rho |
|---|---|---|---|
| cross-species (truth) | 0.029 | 0.061 | 1 |
| **Decant** | **0.027** | **0.057** | **0.985** |
| Decant before the fix | 0.013 | 0.036 | 0.81 |
| SoupX | 0.014 | 0.014 | one value per channel |
| DecontX | 0.024 | 0.093 | 0.95 |

**What was wrong.** Every gene the estimator called absent was truly foreign,
but it kept only the low-ratio part of them (1,584 of 10,549 usable foreign
genes in the largest cluster). Absent genes scatter around rho by about
25-30% from gene to gene, far more than Poisson noise: their variance was 4
to 200 times the Poisson variance. The Poisson test rejected absent genes
above the current rho and kept those below, so the estimate settled near
the bottom of their spread. The test is now a negative binomial whose extra
variance is estimated from the lower tail of the absent genes (10th to 25th
percentile), which expressed genes cannot reach. With no extra variance it
is the Poisson test. A first version that used the median instead failed
`gate_rho` (weakly expressed genes inflated the variance and rho ran 26-36%
high in its stress setting) and was replaced.

**This dataset was used to find the problem**, so it is a development set.
The held-out mixtures in section 1d test the fix on data it was not built
on.

Caveats:
- SoupX here scores 0.345 against the paper's 0.422. The setup differs: no
  extra QC, and SoupX ran without `forceAccept`. The comparison between
  tools within this run is fairer than comparing one of these rows with a
  published row.
- SoupX and DecontX receive the empties as per-gene totals, not as millions
  of droplet columns. Building SoupX's channel and DecontX's background over
  2.95M droplets ran out of memory in 16 GB. Both tools read their empties
  only through row sums (`SoupX::estimateSoup`; decontX's
  `eta_tilda <- rowSums(countsBackground)`), so the results are the same. This
  was checked on the negative control below: per-donor scores were identical
  to the full-matrix run.
- Species mixtures have unusually clean marker genes, which helps the
  absent-gene estimators (Decant, SoupX) more than it would in a
  single-species tissue.

### Negative control: Smart-seq2 human islets (E-MTAB-5061)

Reproduce with `Rscript validation/benchmark_negative_control.R`. Each cell was
sorted into its own well (Segerstolpe et al. 2016), so there is no ambient
RNA. The ideal correction is none, and every gene counts as endogenous.
Following the paper, the donors are HP1506401, HP1507101 and HP1508501T2D,
each run separately. Cells are wells with 100,000 to 500,000 counts and
2,500 to 10,500 detected genes. The "empty" boundaries are set at 1,000 to
100,000 counts. SoupX gets `soupRange = c(1000, 1e5)` with
`forceAccept = TRUE`. DecontX gets every non-cell well as background, as in
the paper, which passes the raw matrix. Decant gets the same soup wells as
SoupX. Clusters are computed on all three donors with Harmony by donor, 10
PCs and resolution 0.3. Scores are means of the three donors.

| | retain | counts removed | rho, mean |
|---|---|---|---|
| **Decant** | **0.998** | **0.17%** | **0.0020** |
| Decant, soup = pooled expression of the donor's cells (no non-cell wells needed) | 0.998 | 0.28% | 0.0036 |
| SoupX, run here | 0.989 | 1.0% | 0.010 |
| DecontX, run here | 0.952 | 5.0% | 0.050 |
| *published:* SoupX | 0.988 | | |
| *published:* DecontX | 0.833 | | |
| *published:* CellBender | 0.797 | | |
| *published:* scAR | 0.401 | | |

Decant estimates near-zero contamination and removes almost nothing, which
is the correct answer here. SoupX matches its published score, which
suggests the setup is close to the paper's. DecontX does not match: 0.952
here against 0.833 published.

Caveats, and why this table should not be read alone:
- **Doing nothing scores 1.0 here.** This control penalises over-removal
  only. Read together with the species mixture, it shows Decant is
  conservative: it is the best tool here for the same reason it is behind
  DecontX and CellBender on removal above.
- SoupX's 0.010 is the lower bound of its search grid, not an estimate, so
  SoupX cannot score higher than about 0.99 on any dataset without ambient
  RNA. For one donor (HP1506401), SoupX's defaults found no marker genes, and
  the wrapper retried once with `soupQuantile = 0.5`, as the SoupX error
  message suggests. The paper does not report doing this.
- Counts are the authors' processed RefSeq count matrix from EBI
  (`pancreas_refseq_rpkms_counts_3514sc.txt`), not the paper's own STARsolo
  re-alignment. The paper's thresholds applied to these counts keep 129, 155
  and 131 cells. The DecontX gap with the paper is most likely due to this
  difference, but that has not been checked.
- The 1k to 100k "soup" wells (124 to 168 per donor) include 150 wells that
  the original authors rated as good cells, only shallower. The soup is
  therefore close to cell expression, as it is in the paper's setup. This
  makes the control easier for estimators that look for genes absent from
  the soup.

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

See `results/differential_expression.csv`. Six control and six disease samples;
disease samples carry more ambient RNA; 80 genes that are disease genes only
in an abundant cell type can leak into the rare type being tested.

Mean of 4 seeds. Every corrector gets the true cell-type labels. The Decant
correction arm uses the `bayes` rule; it has not been rerun with the current
default, `bayes_cell`.

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
| Removal of ambient-compatible foreign counts, against DecontX | behind by 5 to 9 points on all three sets | section 1a |
| Leaving counts that cannot be ambient (likely doublets) | ahead of DecontX, AmbiDose, scAR; level with SoupX and CellBender | section 1a |
| Removal of ambient-compatible foreign counts, against SoupX, CellBender, AmbiDose | ahead on 3 of 3, 2 of 3 and 3 of 3 sets; CellBender in a reduced configuration | section 1a |
| Native signal kept, against scAR | far ahead (1 to 4% lost against 18 to 26%) | section 1a |
| Speed | level with SoupX, faster than the other four on this dataset | section 1a |
| Ambient removal with experimental ground truth (species mixture) | behind DecontX and the published CellBender (0.42 against 0.72 to 0.73); ahead of SoupX run here (0.35), level with published SoupX (0.42); contamination estimated at about half its true level | section 1c |
| Endogenous counts kept with experimental ground truth | level with SoupX on the species mixture (0.991 against 0.992), ahead of DecontX (0.983); best of those run on the Smart-seq2 negative control (0.998 against SoupX 0.989, DecontX 0.952) | section 1c |
| Genotype mixture (kidney) and medium/high-complexity sets of the same paper | unknown, not run | `validate_genotype_mixture.R` exists |
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
