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

## 1d. Ground truth on further human/mouse mixtures (hgmm 1k, 5k, 10k)

Reproduce with `Rscript validation/validate_species_mixture.R`. Three public
10x v3 mixtures of human (HEK293T) and mouse (NIH3T3) cell lines. Doublets
are excluded; the soup is the barcodes with 1 to 100 UMIs; clusters are
k-means (k = 8) on the cells. Every tool gets the same inputs. Truth is each
cell's contamination from its cross-species reads.

**Full test (both genomes).** Mean rho, root-mean-square error against the
per-cell truth, and Pearson correlation with it:

| Dataset | True mean | Decant | Decant, 0.4.0 | DecontX | SoupX |
|---|---|---|---|---|---|
| 1k (1,045 cells) | 0.014 | **0.011** (error 0.023, r 0.88) | 0.001 (0.034, 0.68) | 0.011 (0.032, 0.65) | 0.010 (0.040) |
| 5k (5,247 cells) | 0.040 | **0.030** (0.023, 0.95) | 0.029 (0.024, 0.94) | 0.027 (0.031, 0.91) | 0.016 (0.065) |
| 10k (10,699 cells) | 0.062 | 0.021 (0.101, 0.18) | 0.021 (0.105, 0.11) | 0.026 (0.080, **0.67**) | 0.010 (0.108) |

- On 1k and 5k Decant has the lowest error of the three tools.
- **On 10k Decant fails.** Seven of eight clusters have too few absent genes
  to estimate rho and borrow it from the one that does, so per-cell values
  barely track the truth (r 0.18). DecontX ranks the cells much better. All
  tools underestimate on this dataset.

**Hidden test (the majority species' cells and genes only: mouse for 1k
and 5k, human for 10k).** This removes the other genome, so the estimator
must work as on single-species data. These are
single cell lines: the soup has the same composition as the cells, no gene
is absent, and contamination is not identifiable from expression.

| Dataset | True mean | Decant | Decant, 0.4.0 | DecontX | SoupX |
|---|---|---|---|---|---|
| 1k | 0.004 | 0.13, warns | 0.64, silent | 0.35 | 0.026 |
| 5k | 0.018 | 0.22, warns | 0.25, silent | 0.39 | failed |
| 10k | 0.054 | 0.24, warns | 0.37, silent | 0.09 | 0.010 |

Decant now reports that rho is not identifiable here (`$rho_unidentifiable`)
instead of returning 25-64% without a warning. Its estimate is still wrong
and should not be used; the warning asks for `rho` to be supplied.

**Which data shaped the method.** The 20k mixture (1c) exposed the
underestimate and was used to design the overdispersion-aware test. The 1k
set exposed two further failures, the rho = 0 trap and the silent
implausible estimates, and their fixes were designed on it. The 10k set was
used to try a third change (estimating the spread from the lowest-ratio
genes when the first pass leaves too few): it brought 10k to 0.046 (r 0.44)
but drove PBMC to rho 0.99 (79% of UMIs removed) and the Smart-seq2 control
from 0.002 to 0.015, so it was rejected. **Only the 5k set and PBMC were not
used to make any design decision**; they are the cleanest evidence above.

**Memory.** `decant()` peaks near 9.5 GB on the 10k channel (10,699 cells,
1.9 million empty droplets), about the same as before these changes.

## 1e. Ambient removal on held-out mixtures, and a change that failed

Script: `benchmark_removal_heldout.R` (the per-gene retain/removal scorer of
section 1c, on the hgmm 1k, 5k and 10k mixtures; clusters by the same Seurat
recipe as the 20k mixture). Results: `results/removal_heldout_current.csv`,
`results/removal_heldout_cap_cluster.csv`. The 20k mixture is the
development set: used below to diagnose, not to judge.

**Held-out, current default (`bayes_cell`).** Decant removes the least and
keeps the most:

| mixture | Decant retain | Decant removal | DecontX retain | DecontX removal | SoupX retain | SoupX removal |
|---|---|---|---|---|---|---|
| 1k | 0.998 | 0.489 | 0.995 | 0.653 | 0.994 | 0.767 |
| 5k | 0.993 | 0.464 | 0.988 | 0.656 | 0.988 | 0.609 |
| 10k | 0.986 | 0.378 | 0.969 | 0.602 | 0.963 | 0.662 |
| mean | 0.992 | **0.444** | 0.984 | **0.637** | 0.982 | **0.679** |

**On held-out data the removal gap to DecontX is 19 points, not the 10 seen
on the 20k mixture.**

**Diagnosis (20k, development set).**

| rule | rho | retain | removal |
|---|---|---|---|
| subtract | estimated / true | 0.997 / 0.996 | 0.048 / 0.049 |
| redistribute | estimated / true | 0.994 / 0.994 | 0.067 / 0.069 |
| bayes | estimated / true | 0.983 / 0.982 | 0.718 / 0.741 |
| bayes_cell (default) | estimated / true | 0.989 / 0.988 | 0.616 / 0.639 |
| bayes_cell, no cap | estimated / true | 0.983 / 0.982 | 0.718 / 0.741 |

- Estimation costs about 2 points on the 20k mixture (per-cell rho is close
  to the truth there: median 0.027 vs 0.029, Spearman 0.985).
- Allocation costs about 10: `bayes_cell`'s cap. An entry whose count
  exceeds the 99.9% Poisson quantile of its ambient expectation keeps all but
  that expectation. For a gene whose expected ambient count per cell is below
  0.001 the quantile is 0, so every single count is kept. Foreign genes with
  1 or 2 counts had 0.14 to 0.18 removal, those with over 100 counts 0.83.
- The cap is also what keeps real signal: without it retain falls from
  0.989 to 0.983. The single counts it keeps of foreign genes and of the
  cell's own rarely expressed genes look the same to a per-entry test.

**The change tried** (`cap_test = "cluster"`): correct the cap's test for the
number of cells of the cluster it is applied to (Bonferroni; no new
constant). On the 20k mixture it gave retain 0.983, removal 0.717: it removes
the cap's effect on both sides.

| held-out | retain 1k / 5k / 10k | removal 1k / 5k / 10k | mean removal |
|---|---|---|---|
| current | 0.998 / 0.993 / 0.986 | 0.489 / 0.464 / 0.378 | 0.444 |
| `cap_test = "cluster"` | 0.997 / 0.990 / **0.974** | 0.534 / 0.540 / 0.484 | 0.519 |
| DecontX | 0.995 / 0.988 / 0.969 | 0.653 / 0.656 / 0.602 | 0.637 |

**It fails the pre-stated criteria** (mean removal within 0.03 of DecontX:
it is 0.118 behind; retain >= 0.985 on every mixture: 0.974 on 10k). The
default is unchanged; `cap_test = "cluster"` stays available as a documented
trade-off (more removal, less retained signal). The negative-control and
PBMC checks were not run, as they could not change the decision.

**What was learned.** On the held-out mixtures Decant's contamination
estimate is low (mean 0.0105 vs truth 0.014 on 1k, 0.0295 vs 0.040 on 5k,
0.021 vs 0.062 on 10k; section 1d), unlike on the 20k mixture. The removal
gap there is mostly estimation, not allocation: the allocation fix recovered
7.5 of 19 points. The next step is the estimator on these designs, not the
correction rule.

## 1f. Contamination estimate on the 10k mixture: diagnosis, and a change that failed

The 10k mixture is development data from here on; held-out data judge.

**Diagnosis (hgmm 10k, Seurat clusters).** The absent-gene set collapses.
Each cluster has 5,000 to 13,000 genes that are truly absent and well
covered (the other species' genes), but the fixed point kept 1 to 6. Four of
five clusters were flagged low-power and borrowed one value (0.027) from a
104-cell cluster; truths 0.026 to 0.066, per-cell Spearman 0.66.

- The overdispersion term is estimated from the current absent set. On the
  first pass that set is empty, so it is 0 and the test is Poisson. Absent
  genes scatter 2 to 4 times more than Poisson (log-ratio SD 0.43 to 0.76
  against 0.21 to 0.30), so 62 to 77% of them were rejected and the spread
  never became estimable.
- The empty droplets do not match the cells' ambient profile: in the human
  cluster the typical mouse gene sits at twice the true rho (median ratio
  0.052 against 0.026); only the count-weighted average is right, because
  the most abundant soup genes sit low.
- The fixed point only rejects genes above the current rho, so it settles on
  whichever self-consistent value is nearest its start. Started at the median
  ratio of well-covered genes instead of the lowest decile, it found
  6,000 to 13,000 absent genes and rho within 10% of each cluster's
  count-weighted truth.

**The change tried.** As the default, the high start failed `gate_rho`
(about 0.68 whatever the truth in its stress setting: it pulls expressed
genes into the absent set). Used only when the low start's absent set has
collapsed, and accepted only if well covered and at most 0.5, it passed
`gate_rho`; on the 10k mixture (development) median 0.036 vs truth 0.038,
Spearman 0.98. With it, a result was also declared not identifiable (rho NA,
counts not corrected) when most cells sit in clusters without usable
evidence, instead of only when all do.

| check | current | candidate | bar |
|---|---|---|---|
| hgmm 5k full, held out: median vs truth 0.0256 | 0.0129 (-50%) | 0.0145 (**-43%**) | within 25% |
| hgmm 5k full, held out: per-cell Spearman | 0.35 | **0.49** | >= 0.6 |
| hidden-species tests 1k / 5k / 10k (truth 0.004 / 0.018 / 0.054) | 0.135 / 0.221 / 0.239 | declined (NA) | an honest answer |
| PBMC 4k median (SoupX 0.06 to 0.077) | 0.067 | 0.067 | |
| 20k median, Spearman | 0.027, 0.985 | unchanged | no regression |
| 1k full | 0.0105 | unchanged | no regression |
| Smart-seq2 negative control (truth 0) | 0.0014 | **NA** | no regression |
| `gate_rho` | pass | pass | pass |

**It fails** (5k median and Spearman; the negative control's correct answer
became NA). The estimator is unchanged. Two lessons: the mean and Pearson
correlation of rho flatter an estimator (on 5k the candidate's mean was
within 18% and Pearson 0.975, driven by a minority of heavily contaminated
cells), so medians and rank correlations are the measures to report; and the
"not identifiable" rule needs to tell a sample with no contamination from one
whose contamination cannot be seen, which the share of low-power cells does
not.

## 1g. `correction = "bayes_absent"`: judged on fresh mixtures

Scripts: `benchmark_removal_heldout.R` (`DECANT_HELDOUT=fresh`) and
`check_bayes_absent.R`. Results: `results/removal_heldout_fresh_*.csv`,
`results/bayes_absent_checks_*.csv`. The algorithm was frozen (commit
38a9dc4) and the criteria fixed before hgmm_6k, hgmm_12k and
10k_hgmm_3p_nextgem_Chromium_X were downloaded. None of the three was used
for any design decision.

| fresh mixture | bayes_absent retain / removal | bayes_cell (default) | DecontX | SoupX |
|---|---|---|---|---|
| hgmm_6k | 0.982 / 0.892 | 0.989 / 0.597 | 0.984 / 0.800 | 0.973 / 0.842 |
| hgmm_12k | 0.985 / 0.872 | 0.989 / 0.453 | 0.983 / 0.680 | 0.986 / 0.522 |
| 10k nextgem X | 0.984 / 0.891 | 0.989 / 0.648 | 0.986 / 0.804 | 0.990 / 0.477 |
| mean | 0.9837 / **0.885** | 0.9894 / 0.566 | **0.9847** / 0.761 | 0.9829 / 0.614 |

| pre-stated criterion | result |
|---|---|
| mean removal >= DecontX | **pass** (0.885 vs 0.761; best on every mixture) |
| mean retain >= DecontX | **fail** (0.9837 vs 0.9847) |
| Smart-seq2 negative control retain >= 0.99 | pass (0.998, same as bayes_cell) |
| PBMC native markers no worse than bayes_cell | pass (1.3 / 4.7 / 2.1% lost, unchanged) |
| PBMC counts that cannot be ambient no worse than bayes_cell | **fail**: removed 17 / 33 / 60% (Ig / monocyte / NK) vs 13 / 16 / 32% |
| simulation gate, counts never added | pass (mean L1 0.0627 vs 0.0631; 0 fabricated) |

**It fails two criteria, so the default stays `bayes_cell`.** `bayes_absent`
ships as opt-in for analyses that need the most ambient RNA removed. Its
costs are about 0.5 points less retain than the default and more likely
doublets removed. Suspected cause of the doublet loss (not tested): it always
uses the cluster-level cap test (`cap_test = "cluster"`), which protects fewer
single high counts than the per-entry test. (Correction: three 10x human/mouse
mixtures were still unused at this point; see 1h. The numbers in this section
are for the version frozen at 38a9dc4, before the hybrid cap.)

## 1h. `bayes_absent` with the hybrid cap: judged on the unused mixtures

Pre-registration: `preregistration/bayes_absent_hybrid.md`, committed with the
frozen code (db6f672) before hgmm_1k, hgmm_1k_v2 and hgmm_100 were downloaded.

The change was developed on data already seen (the hgmm v3 mixtures and PBMC
4k). The PBMC doublet-level loss came entirely from the per-cluster cap test.
`cap_test = "hybrid"` also sets aside entries that pass the per-entry test
with at least 5 UMIs. On seen data, removal was the same and PBMC identical to
`bayes_cell`; retain was unchanged, and was predicted to fail again.

hgmm_1k and hgmm_1k_v2 are distinct libraries (5 shared barcodes of about
1,000), so both are primary. hgmm_100 is reported but not scored.

| unused mixture | bayes_absent retain / removal | bayes_cell (default) | DecontX | SoupX |
|---|---|---|---|---|
| hgmm_1k | 0.9934 / 0.899 | 0.9961 / 0.487 | 0.9933 / 0.674 | 0.9918 / 0.634 |
| hgmm_1k_v2 | 0.9903 / 0.882 | 0.9930 / 0.575 | 0.9907 / 0.635 | 0.9945 / 0.423 |
| **mean (primary)** | 0.99183 / **0.891** | **0.99452** / 0.531 | 0.99201 / 0.655 | 0.99314 / 0.529 |
| hgmm_100 (not scored) | 0.9982 / 0.813 | 0.9990 / 0.347 | 0.9981 / 0.583 | failed (no markers) |

| pre-stated criterion | result |
|---|---|
| 1. mean removal >= DecontX | **pass** (0.891 vs 0.655; best on every mixture) |
| 2. mean retain >= DecontX | **fail** (0.99183 vs 0.99201, 0.02 points) |
| 3. Smart-seq2 negative control retain >= 0.99 | pass (0.998) |
| 4. gate: no counts added, L1 within 0.002 of bayes_cell | pass (0.0628 vs 0.0631; 0 added) |
| 5. PBMC non-ambient counts no worse than bayes_cell (regression check, not held out) | pass (12.7 / 16.1 / 32.1%, identical) |

**Criterion 2 fails, as predicted, so `bayes_absent` stays opt-in and the
default stays `bayes_cell`.** The doublet problem of 1g is fixed. The
remaining cost is retain: on six fresh or unused mixtures it keeps
0.2 to 0.7 points less of the cells' own counts than the default, and about
as much as DecontX (behind by 0.1 and 0.02 points in the two tests), while
removing 0.89 of ambient signal against DecontX's 0.66 to 0.76. Part of that
retain loss is within-species ambient RNA, which this scorer cannot tell from
native counts (on hgmm 10k, half of the native-species genes called absent
in a cluster were 3-fold enriched in another cluster of the same species).
That is a property of the scorer, stated here, not a pass. No 10x human/mouse
mixture is left unused, so nothing further is tuned on these.

**Decision after the result: `bayes_absent` made the default anyway.** The
pre-registered rule said it stays opt-in. The maintainer overrode that rule,
after seeing the result above, because the default removed the least ambient
RNA of the tools compared, while `bayes_absent` removed the most for 0.02
points of retain against DecontX. This is recorded as a maintainer decision,
not as a pass. What was known when it was made:

- Retain cost against `bayes_cell`: 0.27 points on the unused mixtures, 0.57
  on the fresh ones (1g, earlier cap), 2.5 on hgmm 10k v3 (0.961 against
  0.986, seen in design).
- On PBMC 4k the absence posterior was not used in any of the 11 clusters
  (the separation guard: components overlap), and the output matched
  `bayes_cell`. Its gain is shown only on species mixtures, where absent
  genes are sharply separated. Whether it fires, and helps, on real tissue
  with very distinct cell types is untested.
- `gate_correction()` now scores `bayes_absent` and fails CI unless the
  default has the lowest mean L1 error: 0.0628 against 0.0631 for
  `bayes_cell`. The margin is narrow.

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

Mean of 4 seeds, run before `ambient_de` switched to the trended dispersion
(section 3b); its gate (`gate_ambient_de`) gives similar numbers after the
switch. Every corrector gets the true cell-type labels. The Decant
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

## 3b. Calibration of `ambient_de` on real data

Script: `validate_de_calibration.R`; results:
`results/de_calibration_gse272082_pfc.csv`,
`results/de_calibration_gse190015_ctrl.csv`.

A calibrated test, given condition labels that carry no effect, reports a
hit in few shuffles and about 1% of genes at p < 0.01. Labels were shuffled
between donors 24 times (seed 7; the same shuffles for every method), and
microglia and oligodendrocytes were tested by pseudobulk per donor. Cell
types come from a fixed recipe (Seurat clustering, clusters labelled by
marker genes).

Data:

- **GSE272082** (sporadic early-onset Alzheimer's disease, 10x Multiome,
  gene expression used), prefrontal cortex, 5 control and 4 patient donors,
  about 31,600 nuclei, 340,000 to 620,000 empty droplets per sample. At most
  2 of the 4 patients are in the shuffled patient arm, so these shuffles are
  close to, but not exactly, null. This dataset was used to find and
  diagnose the problem.
- **GSE190015**, held out: nine control donors (occipital cortex, one
  batch), about 69,000 nuclei (barcodes with >= 500 UMIs), shuffled 4 vs 5.
  No disease in either arm: every hit is false. Not used until the fix was
  made.

| data | cell type | method | shuffles with a hit (of 24) | max hits | median share p < 0.01 |
|---|---|---|---|---|---|
| GSE272082 | microglia | `ambient_de` | 2 | 9 | 0.009 |
| | | `ambient_de`, tagwise (before the fix) | 21 | 90 | 0.016 |
| | | naive edgeR | 1 | 1 | 0.007 |
| | oligodendrocytes | `ambient_de` | 2 | 1 | 0.005 |
| | | `ambient_de`, tagwise (before the fix) | 19 | 574 | 0.011 |
| | | naive edgeR | 0 | 0 | 0.004 |
| GSE190015 | microglia | `ambient_de` | 1 | 1 | 0.009 |
| | | `ambient_de`, tagwise (before the fix) | 21 | 125 | 0.016 |
| | | naive edgeR | 1 | 1 | 0.007 |
| | oligodendrocytes | `ambient_de` | 4 | 1 | 0.007 |
| | | `ambient_de`, tagwise (before the fix) | 21 | 176 | 0.014 |
| | | naive edgeR | 2 | 1 | 0.006 |

**Cause.** `ambient_de` fitted each gene with edgeR's tagwise dispersion (the
gene's own estimate, shrunk) and then applied a quasi-likelihood F-test,
which estimates the gene's extra variance again from the residual deviance.
The gene's own variance was absorbed twice, the residual variance came out
too small, and the F-test was liberal. edgeR's quasi-likelihood pipeline
uses the trended dispersion for this reason; `ambient_de` now does too
(`dispersion = "trended"`). The diagnosis on GSE272082, microglia: with the
ambient term set to zero, the tagwise version still reported hits in 22 of 24
shuffles (so the ambient model was not the cause); with the trended
dispersion, 2 of 24, with or without the ambient term. Turning off the Huber
down-weighting made no difference.

`joint_rho` was suspected first. Turning it off also removes every hit, but
by breaking the test: without it, each sample's ambient fraction comes from
that sample's pseudobulk alone, which for a single cell type is not
identifiable (0.78 to 0.83 for some donors against 0.02 to 0.06 for the
others), most of their counts are treated as ambient, and the share of genes
at p < 0.01 falls to 0.1 to 0.3%. It stays on. (Since the ambient-fraction
guard of section 3c, `joint_rho = FALSE` replaces those values with the median
of the other samples instead; one GSE272082 oligodendrocyte shuffle then gives
91 hits.)

The simulation gates pass with the trended dispersion, with power unchanged
(`gate_ambient_de`: null false-positive rate 0.045 to 0.054, power 0.90 to
0.975, no leaked gene; `gate_ambient_de_cohort` passes).

On the true GSE272082 labels, `ambient_de` now reports 1 gene in microglia
and 2 in oligodendrocytes (before the fix 20 and 292); naive edgeR 0 and 3.

What this does and does not show: `ambient_de` is now about as well
calibrated as naive edgeR on real data, and slightly more liberal (9 against
4 tests with a hit, of 96). It does not show that `ambient_de` removes
ambient-driven hits on real data: on these donors naive edgeR found almost
nothing to remove. That remains a simulation result.

This dataset was first used for a pre-registered test of whether the genes
`ambient_de` flags as ambient-driven have closed chromatin in the tested cell
type (from the ATAC side of the same nuclei). That test could not be run:
naive edgeR found 0 genes in microglia and 3 in oligodendrocytes, so there
was nothing to flag.

## 3c. Pre-registered test on type 1 diabetes pancreas (GSE273597)

Pre-registration: `preregistration/t1d_leakage.md` (committed before any DE,
with two amendments recorded before any DE, each with its reason). Script:
`validate_de_leakage_t1d.R`. Results: `results/t1d_leakage_*.csv` (as
pre-registered) and `results/t1d_leakage_*_guard.csv` (rerun with the
ambient-fraction guard, below).

Design: whole-pancreas snRNA-seq, 11 non-diabetic vs 12 T1D donors
(Melton et al. 2025); acinar, ductal and alpha pseudobulks. T1D destroys
beta cells, and soup insulin varies 61-fold between libraries; insulin in
acinar, alpha and ductal pseudobulks tracks it (Spearman 0.74 to 0.81). A
gene called DE in a cell type that does not express it (>= 100-fold below
another type in Segerstolpe Smart-seq2 islets) counts as an ambient-driven
false positive.

**Primary result: null.** Naive edgeR called no gene at all in acinar cells,
foreign or not (pre-stated precondition: at least 5 foreign genes). Between
donors, soup composition varies as much within each group as between them,
so 11 vs 12 donors did not separate. Every method called 0 or 1 foreign
genes in every test, including the contrast with long-standing T1D only
(5 donors). These data do not show ambient-driven DE from a standard test,
and so cannot show `ambient_de` removing it.

**Calibration failed as pre-registered**, in ductal cells: in 5 of 24
shuffled labellings `ambient_de` reported 1,106 to 2,300 genes (naive edgeR
0), and 1,064 on the true labels. Cause: for 6 of 23 donors the ductal
ambient fraction came out at 0.86 to 0.96 (others 0.02 to 0.16). In the five
bad shuffles, setting the ambient term to zero, or replacing those six values
with the median of the others, removed every hit; widening their variance
100-fold, turning off the Huber weights or `joint_rho`, did not.

**Guard** (`max_rho = 0.5`, the plausibility limit of
`estimate_rho_cluster`): after the joint refinement, a sample above it takes
the median of the plausible samples, with their between-sample variance
added to its uncertainty; `ambient_de` warns and names it. A first version
also replaced implausible single-sample start values before the joint
refinement; on held-out GSE190015 oligodendrocytes that produced a shuffle
with 76 hits (it fired on four samples the joint step corrects by itself),
so it was removed.

| cell type (24 shuffles each) | `ambient_de` as pre-registered | with the guard | naive edgeR |
|---|---|---|---|
| acinar | 1 shuffle with a hit, max 1 gene | 1, max 1 | 0 |
| ductal | 5, max 2,300 | 4, max 1 | 0 |
| alpha | 6, max 5 | 6, max 5 | 4, max 3 |

With the guard, `ambient_de` reports 0 genes in ductal cells on the true
labels (1,064 before). GSE190015 and GSE272082 are unchanged by it with
default settings. Alpha cells still exceed the pre-stated limit of 5 shuffles
with a hit (6), though never more than 5 genes; naive edgeR also has hits in 4
of 24 there (17 donors, often fewer than 100 alpha cells each).

What this does not show: the pancreas data were used to find the problem, so
they do not show that the guard works on other data; and there is still no
real dataset on which `ambient_de` beats a standard test.

## 3d. Pre-registered test: paired brain regions (GSE272082)

Pre-registration: `preregistration/region_leakage.md` (committed before any
DE; Amendment 1, a crash fix in `ambient_de`'s robust weights, made before
any result was seen). Script: `validate_de_region_leakage.R`. Results:
`results/region_leakage_*.csv`.

Design: six donors with prefrontal cortex (PFC), entorhinal cortex (EC) and
hippocampus (HIP); microglia (primary) and oligodendrocytes; `~ donor +
region`. Truth: genes foreign to each type in Darmanis et al. 2015 (C1 single
cells; microglia 1,435 genes including SNAP25, SYT1, GAD1; oligodendrocytes
690, without canonical neuronal genes because the reference's
oligodendrocytes carry neuronal RNA at about 10%). One donor's PFC had 9
microglia, below the minimum of 20, so the microglia tests have 5 donors.

**Primary result: null.** Microglia, HIP vs PFC: standard edgeR called no
gene at all (precondition: at least 5 foreign genes), and neither did any
other method.

| test | standard edgeR: foreign / other | `ambient_de` | Decant, SoupX, DecontX then edgeR |
|---|---|---|---|
| microglia HIP vs PFC (primary) | 0 / 0 | 0 / 0 | 0 / 0 each |
| microglia EC vs PFC | 0 / 0 | 0 / 0 | 0 / 0 each |
| oligodendrocytes HIP vs PFC | 3 / 37 | 2 / 45 | 0 / 0, 0 / 0, 0 / 1 |
| oligodendrocytes EC vs PFC | 0 / 0 | 0 / 0 | 0 / 0 each |

Calibration, region labels swapped within donors:

| test | `ambient_de`: relabelings with a hit (max genes) | standard edgeR (max genes) |
|---|---|---|
| microglia HIP | 0 of 15 | 0 of 15 |
| microglia EC | **8 of 15** (1) | 0 of 15 |
| oligodendrocytes HIP | 2 of 31 (1) | 3 of 31 (**107**) |
| oligodendrocytes EC | 3 of 31 (2) | 2 of 31 (**393**) |

- The precondition failed, so the test says nothing about whether
  `ambient_de` removes ambient-driven hits. Region differences in the soup
  did not make a standard test call genes foreign to microglia here.
- Microglia EC vs PFC exceeds the pre-stated calibration limit (6 of 31,
  i.e. about 3 of 15): `ambient_de` has a hit in 8 of 15 relabelings, one gene
  each time.
- Not pre-registered as an endpoint, but noted: in two oligodendrocyte
  relabelings standard edgeR reported 107 and 393 genes; `ambient_de` never
  more than 2.

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
| DE calibration on real data (shuffled labels) | slightly more liberal than naive edgeR: hits in 20 of 168 shuffles against 8, at most 9 genes (three datasets, one held out). Before the dispersion fix and the ambient-fraction guard, far behind (up to 2,300 genes) | sections 3b, 3c |
| DE: removes ambient-driven hits on real data | unknown: two pre-registered tests (pancreas, paired brain regions) were null, because the standard test found no ambient-driven hits to remove | sections 3c, 3d |
| DE with covariates in large cohorts | ahead of subtract-then-rank-test in simulation | `experiments/large_cohort_calibration/` |
| Removal of ambient-compatible foreign counts, against DecontX | behind by 5 to 9 points on all three sets | section 1a |
| Leaving counts that cannot be ambient (likely doublets) | ahead of DecontX, AmbiDose, scAR; level with SoupX and CellBender | section 1a |
| Removal of ambient-compatible foreign counts, against SoupX, CellBender, AmbiDose | ahead on 3 of 3, 2 of 3 and 3 of 3 sets; CellBender in a reduced configuration | section 1a |
| Native signal kept, against scAR | far ahead (1 to 4% lost against 18 to 26%) | section 1a |
| Speed | level with SoupX, faster than the other four on this dataset | section 1a |
| Ambient removal with experimental ground truth (20k species mixture) | behind DecontX and the published CellBender (0.62 against 0.72 to 0.73); well ahead of SoupX run here (0.35) and published (0.42) | section 1c |
| Ambient removal with the default (`bayes_absent`) on unused mixtures | **ahead of every tool run**: 0.891 against DecontX 0.655, SoupX 0.529; retain 0.9918 against DecontX 0.9920 (missed the pre-registered bar by 0.02 points; made default by maintainer decision) | section 1h |
| Ambient removal with `bayes_cell` (default until 1h) on held-out mixtures (hgmm 1k, 5k, 10k) | behind DecontX by 19 points (0.44 against 0.64) and SoupX (0.68); keeps the most (retain 0.992 against 0.984, 0.982) | section 1e |
| Native signal kept by the default on hard data | behind `bayes_cell` by 2.5 points on hgmm 10k (0.961 against 0.986) | section 1h |
| Default's gain on real tissue | unknown: on PBMC 4k its gene-level call never applied (output identical to `bayes_cell`) | section 1h |
| Contamination estimate on held-out species mixtures | lowest error of Decant, SoupX and DecontX on 1k and 5k; behind DecontX on 10k (per-cell correlation 0.18 against 0.67), where all tools underestimate | section 1d |
| Endogenous counts kept with experimental ground truth | level with SoupX on the species mixture (0.989 against 0.992), ahead of DecontX (0.983); best of those run on the Smart-seq2 negative control (0.998 against SoupX 0.989, DecontX 0.952) | section 1c |
| Contamination not identifiable (single cell type) | Decant warns (DecontX returns 9-39% silently; SoupX failed on one set); the number Decant reports is still wrong | section 1d |
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
