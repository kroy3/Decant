# Pre-registration: ambient-driven DE in type 1 diabetes pancreas

Written before any differential expression was computed on these data. The
commit that adds this file is the timestamp. Changes after that are listed
as amendments at the end, with the reason, and are made before the result
they could affect is seen where possible.

## Question

On real case/control data, does a standard pseudobulk test report genes as
differentially expressed in a cell type that does not express them, because
the ambient RNA differs between conditions? And does `ambient_de()` (or
correcting counts first) stop it?

## Data

GSE273597 (Melton et al. 2025, *Science Advances*, doi
10.1126/sciadv.ady0080): whole-pancreas snRNA-seq (10x), nPOD donors, 33
libraries from 32 donors. Raw matrices (`raw_feature_bc_matrix.h5`) for all.
Donor status from the paper's Supplementary Table 1
(`validation/sciadv.ady0080_tables_s1_to_s17.zip`), extracted to
`validation/t1d_donor_labels.csv`: 11 non-diabetic (ND), 9 autoantibody
positive, 7 early T1D, 5 long-standing T1D. Donor 6282 has two libraries;
they are pooled.

Why this dataset: T1D destroys beta cells, so the ambient RNA of T1D donors
carries far less beta-cell RNA. A label-free check (before this file) found
soup insulin varies 61-fold between libraries (10th to 90th percentile) and
insulin in acinar, alpha and ductal pseudobulks tracks it (Spearman 0.79,
0.81, 0.74; log-log slopes 1.19, 0.90, 0.81).

## Cells and cell types (fixed)

- Cells: barcodes with >= 500 UMIs and >= 300 detected genes. Empty
  droplets: 1 to 100 UMIs (`split_droplets` defaults).
- Cell types from marker genes, per cell: score = mean log1p(CP10k) over each
  set (beta INS; alpha GCG; delta SST; acinar PRSS1, CPA1, CTRB2, CELA3A;
  ductal KRT19, CFTR). A cell takes the top set if its score exceeds the
  second by more than 1; otherwise it is unassigned.
- Pseudobulk per donor and cell type. A donor enters a cell type's test only
  with >= 20 cells of that type.

## Truth (independent of these data)

Segerstolpe et al. 2016 (E-MTAB-5061), Smart-seq2 islet cells, non-diabetic
donors only, authors' cell-type labels. Plate-based: essentially no ambient
RNA, different platform, different donors.

For a tested type T, a gene is **foreign to T** if, in Segerstolpe:
- some other type S has >= 100 CPM (pooled counts), and
- T has <= 1/100 of that type's CPM, and
- the gene is detected (count > 0) in < 5% of T cells.

A foreign gene called DE in T (FDR < 0.05, either direction) counts as an
**ambient-driven false positive**. Genes not foreign to T are not scored as
right or wrong (T1D also changes non-beta cells for real).

## Tests

Primary: acinar cells, ND (11 donors) vs T1D (12 donors), design `~ condition`.
Secondary: ductal and alpha cells; ND vs long-standing T1D only.

Methods (FDR < 0.05 throughout):
1. naive edgeR quasi-likelihood on raw pseudobulk (`ambient_de(naive = TRUE)`);
2. `ambient_de()` with default settings;
3. correct counts then the same edgeR test: Decant (`decant()`, default
   `bayes_cell`), SoupX, DecontX, each run per library with the empty
   droplets and the marker cell types as clusters.

Calibration: 24 shuffles of the condition labels among the 23 donors (seed 7),
`ambient_de` and naive edgeR, as in `validate_de_calibration.R`.

## Pre-stated outcomes (primary test)

- **Precondition:** naive edgeR calls >= 5 foreign genes. Otherwise the test
  is reported as null (no ambient-driven DE to remove).
- **Support:** `ambient_de` calls at most half as many foreign genes as naive
  edgeR, and its shuffled-label calibration stays within what was seen on
  GSE272082 and GSE190015 (hits in at most 5 of 24 shuffles).
- **Against:** `ambient_de` calls more than half as many foreign genes as
  naive edgeR.
- Correct-then-test methods are reported alongside, without a threshold.
- For transparency, the number of non-foreign hits per method is reported
  (unscored; a method that removes everything also removes foreign hits).

## Known weaknesses

- Alpha cells can contain real insulin (bihormonal cells, doublets): acinar
  is the primary test for that reason.
- Segerstolpe cells are whole cells (cytoplasmic RNA), these are nuclei;
  specificity of most marker genes carries over, but some genes may differ.
- 23 donors; effects of T1D on exocrine tissue are real and unscored.

## Amendments

**1 (before any DE was computed).** Method 1 is changed from
`ambient_de(naive = TRUE)` to a standalone standard pipeline: `filterByExpr`
on the raw pseudobulk, TMM, `estimateDisp`, `glmQLFit`, `glmQLFTest`. The
naive test inside `ambient_de()` only covers genes that pass `ambient_de`'s
own filter, which skips genes whose counts are mostly ambient: exactly the
genes most likely to leak. Using it would understate the leakage of the
standard test. The correct-then-test methods use the same standalone
pipeline. Genes that `ambient_de` declines to test as mostly ambient count as
not called, and their number is reported.
