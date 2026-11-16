# Pre-registration: ambient-driven DE between brain regions, paired by donor

Written before any differential expression was computed on these data. The
commit that adds this file is the timestamp. Amendments are listed at the end,
each made before the result it could affect, with its reason.

## Question

Does a standard pseudobulk test report genes as differentially expressed in
microglia (or oligodendrocytes) between two brain regions of the same donors
because the regions' ambient RNA differs in neuronal composition, and does
`ambient_de()` (or correcting counts first) stop it?

## Data

GSE272082 (10x Multiome, gene expression used). Six donors with prefrontal
cortex (PFC), entorhinal cortex (EC) and hippocampus (HIP) libraries:
NIH01/02/03, NIH04/05/06, NIH10/11/12, NIH13/14/15, NIH16/17/18, NIH28/29/30
(PFC/EC/HIP in that order; donors identified by sex, age and post-mortem
interval in GEO). Disease status is a donor property and is absorbed by the
donor term.

## Cells, cell types, pseudobulk (fixed)

- Cells: Cell Ranger ARC filtered barcodes with >= 300 detected genes.
  Empty droplets: raw barcodes with 1 to 100 UMIs.
- Cell types, per cell, from marker scores (mean log1p CP10k over each set):
  microglia P2RY12, CSF1R, CX3CR1; oligodendrocyte MBP, PLP1, MOBP;
  astrocyte GFAP, AQP4, SLC1A2; OPC PDGFRA, VCAN; excitatory SLC17A7, SATB2;
  inhibitory GAD1, GAD2. A cell takes the top set if its score exceeds the
  second by more than 1; otherwise it is unassigned.
- Pseudobulk per library and cell type. A library enters if it has >= 20
  cells of the type; a donor enters a contrast only if both of its regions do.

## Truth (independent of these data)

Darmanis et al. 2015 (PNAS; GSE67835): single cells from adult human
temporal cortex captured in Fluidigm C1 chambers (not droplets), authors'
labels: 131 neurons, 62 astrocytes, 38 oligodendrocytes, 18 OPC, 16 microglia,
20 endothelial cells.

A gene is **foreign to type T** if another type has >= 100 CPM (pooled
counts) and T has <= 1/100 of that.

Reference check, done before this rule was fixed:
- Microglia are clean for neuronal genes: SNAP25 0, SYT1 4 CPM against
  2,579 and 2,398 in neurons (>= 600-fold). They carry oligodendrocyte RNA
  (PLP1 3,567 CPM, about 10% of oligodendrocytes), so PLP1 and MBP are not
  foreign to microglia.
- Oligodendrocytes carry neuronal RNA at about 10% of the neuronal level
  (SNAP25 252, SYT1 269, GAD1 134 CPM; detected in 47 to 68% of them): C1
  capture is not free of ambient RNA. Under the 100-fold rule no canonical
  neuronal gene is foreign to oligodendrocytes. The rule is kept, not
  loosened: the oligodendrocyte test is therefore conservative and probably
  underpowered, and is secondary.
- Truth sets: microglia 1,435 genes (including SNAP25, SYT1, GAD1, CAMK2A,
  STMN2, GAP43, NEFL); oligodendrocytes 690 genes.

A foreign gene called DE in T (FDR < 0.05, either direction) counts as an
ambient-driven false positive. Other genes are not scored: microglia and
oligodendrocytes do differ between regions.

## Tests

Primary: **microglia, HIP vs PFC**, design `~ donor + region`.
Secondary: microglia EC vs PFC; oligodendrocytes HIP vs PFC and EC vs PFC.

Methods (FDR < 0.05 throughout):
1. standard edgeR quasi-likelihood on the raw pseudobulk: `filterByExpr`,
   TMM, `estimateDisp`, `glmQLFit`, `glmQLFTest` on the region coefficient;
2. `ambient_de()` with default settings, empty droplets as per-library
   per-gene totals;
3. correct counts then the same edgeR test: Decant (`decant()`, defaults),
   SoupX, DecontX, each per library with the marker cell types as clusters.

Calibration: region labels swapped within donors. With 6 donors there are 31
distinct relabelings besides the true one (a relabeling and its mirror image
give the same test); all 31 are run for `ambient_de` and standard edgeR.

## Pre-stated outcomes (primary test)

- **Precondition:** standard edgeR calls >= 5 foreign genes. Otherwise the
  test is reported as null.
- **Support:** `ambient_de` calls at most half as many foreign genes as
  standard edgeR, and reports a hit in at most 6 of the 31 relabelings.
- **Against:** `ambient_de` calls more than half as many foreign genes as
  standard edgeR.
- Calibration failure (hits in more than 6 relabelings) is reported as such,
  whatever the foreign-gene counts.
- Correct-then-test methods are reported alongside without a threshold, and
  non-foreign hits per method are reported (unscored).

## README

If the test supports `ambient_de`, the result goes in the README. Otherwise
the README must not imply a real-data advantage.

## Known weaknesses

- Six donors.
- The reference is temporal cortex: hippocampus-specific neuronal genes may
  not appear among its neurons and so not in the truth set (conservative).
- The reference itself carries trace ambient RNA (see above).
- Microglia differ between regions for real; only foreign genes are scored.
