# Chapter 2. Key Concepts

This chapter introduces the ideas the rest of the guide relies on. There is
no code here.

## 2.1 Ambient RNA

During tissue dissociation some cells break open and release their RNA into
the cell suspension. When the suspension is loaded onto a droplet device,
this free RNA is captured in every droplet together with the intact cell.
After sequencing, each cell's counts are a mixture:

> observed counts = the cell's own RNA + ambient RNA

Ambient RNA is dominated by the most abundant and most fragile cell types.
In a brain sample, for example, neuronal transcripts can appear in microglia;
in a blood sample, haemoglobin and immunoglobulin genes appear in every cell.

## 2.2 The soup

The **soup** (or **ambient profile**) is the composition of the ambient RNA:
for every gene, the fraction of ambient molecules that belong to that gene.
It is the same for every droplet in a sample, because all droplets draw from
the same suspension.

## 2.3 Empty droplets

Most droplets in a run contain no cell. They still capture ambient RNA, so
their counts are a direct measurement of the soup. A typical 10x run has
hundreds of thousands of empty droplets with between 1 and 100 UMIs each.
Decant adds them up to estimate the soup.

## 2.4 The contamination fraction (rho)

The **contamination fraction**, written **rho** (ρ), is the share of a
cell's counts that are ambient. A rho of 0.05 means 5% of the cell's UMIs
came from the soup. Typical values are 0.02 to 0.10 for single-cell data
and higher, up to 0.30, for single-nucleus data. Rho differs from cell to
cell: smaller cells and nuclei carry relatively more ambient RNA.

## 2.5 How Decant estimates rho

Some genes are not expressed by a given cell type at all. Haemoglobin, for
example, is not made by T cells. Every haemoglobin count in a T cell must
therefore be ambient. Decant finds, for each cluster, the genes that the
cluster does not express (the **absent genes**) and asks how much soup would
be needed to explain the counts it sees for them. That amount is rho. The
cluster estimate is then refined for each cell.

This is why Decant needs **clusters**: it looks for absent genes cluster by
cluster. Clusters from your usual Seurat analysis work well.

## 2.6 How Decant removes contamination

Once rho and the soup are known, the expected ambient count of every gene in
every cell is `rho × cell total × soup fraction`. Decant removes that
expectation, gene by gene, with two guarantees:

- a corrected count is never negative, and
- a corrected count is never larger than the observed count.

Decant offers four rules for distributing the removal; the default,
`bayes_cell`, is explained in [Chapter 7](07-Estimation-and-Correction).

## 2.7 Ambient RNA and differential expression

In a case/control study, samples can carry different amounts of soup. If
diseased tissue is more fragile, its samples contain more ambient RNA, and
that RNA is rich in the disease genes of the abundant cell types. A
standard test on a rare cell type then reports those genes as changed in the
rare type. Decant's `ambient_de()` puts each sample's measured soup into
the statistical model, so the test separates real expression changes from
changes in the soup ([Chapter 9](09-Differential-Expression)).

## 2.8 Glossary

| Term | Meaning |
|---|---|
| Ambient RNA, soup | Cell-free RNA captured in every droplet |
| Soup profile | Fraction of ambient molecules per gene |
| Empty droplet | A droplet without a cell; it measures the soup |
| rho (ρ) | Fraction of a cell's counts that are ambient |
| Absent gene | A gene a cluster does not express; its counts in that cluster are ambient |
| UMI | Unique molecular identifier; one UMI is one captured molecule |
| Pseudobulk | Counts of one cell type summed per sample |
| Gate | An automated test against a known answer that a component must pass |

---

**Previous:** [Chapter 1, Installation](01-Installation) | **Next:** [Chapter 3, Quick Start](03-Quick-Start)
