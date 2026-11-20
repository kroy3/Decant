# Pre-registration: does the default correction remove more ambient RNA in real tissue?

Written before any correction was run on these data for this question. The
commit that adds this file, together with the script, is the timestamp.
Amendments are listed at the end, each made before the result it could
affect, with its reason.

## Question

`bayes_absent` became the default because it removed far more ambient RNA
than `bayes_cell` on human/mouse mixtures (0.891 against 0.531; BENCHMARK.md
1h). On PBMC 4k its gene-level call never applied, because the clusters share
most genes, and the output matched `bayes_cell`. Does the gain appear in real
tissue whose cell types differ sharply, here brain?

## Data

GSE272082 (10x Multiome; gene expression used): the six prefrontal cortex
libraries NIH01, NIH04, NIH10, NIH13, NIH16 and NIH28. These were used before
only for differential expression tests (BENCHMARK.md 3b and 3d), never in
designing or tuning a correction.

- Cells: Cell Ranger ARC filtered barcodes with >= 300 detected genes.
- Empty droplets: raw barcodes with 1 to 100 UMIs.
- Cell types and clusters: the marker-score recipe of
  `preregistration/region_leakage.md`, fixed there:
  - microglia: P2RY12, CSF1R, CX3CR1;
  - oligodendrocyte: MBP, PLP1, MOBP;
  - astrocyte: GFAP, AQP4, SLC1A2;
  - OPC: PDGFRA, VCAN;
  - excitatory neuron: SLC17A7, SATB2;
  - inhibitory neuron: GAD1, GAD2;
  - unassigned: cells whose top score does not exceed the second by more
    than 1.

  These labels are the clusters given to every method.

## Truth (independent of these data)

Darmanis et al. 2015 (GSE67835): adult temporal cortex captured in Fluidigm
C1 chambers, not droplets, with the authors' labels. CPM is pooled per type.

- A gene is **foreign** to type T if another type has >= 100 CPM and T has
  <= 1/100 of that. These are the same rules as in `region_leakage.md`.
- A gene is **native-specific** to T if T has >= 100 CPM and at least 10
  times every other type.

Scored types:
- **microglia**, primary: the cleanest truth; neuronal genes are at least
  600-fold below neurons;
- **oligodendrocytes** and **astrocytes**, secondary.

The reference carries trace ambient RNA (oligodendrocytes about 10% of
neuronal levels), so its foreign sets are conservative.

## Methods

All are run per library, with the same cells, empty droplets and labels:
- `decant()` with `correction = "bayes_absent"` (the default);
- `decant()` with `correction = "bayes_cell"`;
- SoupX and DecontX, through `validation/comparators.R`.

## Metrics, per library and scored type T (cells labelled T)

- **foreign removed:** 1 - (corrected foreign counts) / (observed foreign
  counts). These counts are ambient RNA or doublets; higher is better.
- **native removed:** the same over native-specific genes. Some ambient share
  is correct to remove (about rho), so lower is better.
- **doublet-level kept:** the share kept of foreign entries too high to be
  ambient, which are likely doublets. An entry is too high if its count is
  improbable under ambient alone: P(X >= x) < 0.001 at a contamination bound
  of 0.15, as in the PBMC benchmark.

The reported figure is the mean over the six libraries.

The script also records, for every library and cluster, whether
`bayes_absent`'s gene-level call was used (the separation guard) and how
many genes it called absent.

## Pre-stated outcomes (microglia)

- **The gain holds in real tissue:** `bayes_absent` foreign removed >= 
  `bayes_cell` + 0.05, and native removed <= `bayes_cell` + 0.01.
- **No gain (falls back):** foreign removed within 0.01 of `bayes_cell`.
  Then the default's advantage stays shown only on species mixtures, and the
  README says so.
- **Harmful:** native removed > `bayes_cell` + 0.01, whatever the removal.
  Then the README states the cost.
- Anything else is reported as in between, with its numbers.
- Doublet-level kept is reported for every method, without a threshold.
- SoupX and DecontX are reported alongside, without a threshold.
- Oligodendrocytes and astrocytes are reported with the same three outcomes
  but are secondary.

Whatever the outcome, `bayes_absent` is not tuned on these data.

## Known weaknesses

- **There is no exact truth.** Foreign counts in microglia are taken as
  ambient or doublet, but microglia may carry phagocytosed neuronal RNA.
- **The reference is temporal cortex** and these data are prefrontal cortex.
- **Labels from three marker genes leave some mixed cells.** They are a
  fixed recipe, not tuned here.

## Amendments
