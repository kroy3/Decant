# Chapter 11. Diagnostics and Troubleshooting

Decant reports when its estimate is uncertain. This chapter explains each
message and what to do about it.

## 11.1 Three checks after every run

```r
library(Decant)
sim <- simulate_experiment(seed = 1)
res <- decant(sim$observed, sim$empty, clusters = sim$labels)

quantile(res$rho, c(0.1, 0.5, 0.9))
#>       10%       50%       90%
#> 0.1083049 0.1922303 0.2939032
res$soup_mismatch
#> [1] FALSE
any(res$rho_diagnostics$low_power)
#> [1] FALSE
```

1. **Is rho plausible?** Single-cell data typically show 0.02 to 0.10 and
   single-nucleus data up to 0.30. A median far outside this range deserves
   a closer look.
2. **`soup_mismatch`** should be `FALSE`.
3. **`low_power`** should be `FALSE` for all clusters, or at least for the
   large ones.

## 11.2 Warning: soup/cell mismatch

> *X% of the soup is in N gene(s) the cells do not contain ... rho is
> UNDERestimated.*

**Meaning.** A noticeable share of the soup (more than 5%) is in genes that
no cell expresses. The empty droplets contain material that is not ambient
RNA from these cells, for example debris, red blood cells, or a cell type
that was removed by filtering.

**What Decant does.** It removes those genes from the soup before
estimating rho.

**What to do.**
- Look at the top soup genes: `head(sort(res$ambient[, 1], decreasing = TRUE), 20)`.
- If they belong to a cell type you filtered out, consider keeping those
  cells.
- If rho still looks too low, supply an independent estimate with
  `rho = ...` ([Chapter 5.6](05-Running-Decant)).

## 11.3 Warning: low-power clusters

> *N cluster(s) have almost no genes that are cleanly absent, so rho is
> weakly identified there and likely OVERestimated.*

**Meaning.** In these clusters almost every gene in the soup is also
expressed, so there is little evidence about the soup alone. This happens
with very small clusters and with cell types whose expression resembles the
soup.

**What Decant does.** It gives these clusters the rho of the other
clusters, weighted by their evidence (`borrowed = TRUE` in
`res$rho_diagnostics`).

**What to do.** Merge very small clusters, or cluster at a lower
resolution, and run again.

## 11.4 Common errors

| Message | Cause | Fix |
|---|---|---|
| `empties needs nrow(obs_total) genes` | Cells and empty droplets have different genes | Create both with `split_droplets()` from the same raw matrix |
| `length(clusters) must equal ncol(obs_total)` | The cluster vector does not match the cells | Use one label per cell, in column order |
| `` `sample_of` is required when `empties` is a list `` | Several soups but no sample assignment | Pass `sample_of` ([Chapter 8](08-Multi-Sample-Experiments)) |
| `need at least 2 cells` | Too few cells | Check the cell filtering |
| `RunDecant()` refuses the object | The assay holds normalised or SCTransformed values | Use the assay with raw UMI counts (`assay = "RNA"`) |

## 11.5 Frequently asked questions

**Do I need the raw matrix?**
It is strongly recommended, because the empty droplets measure the soup.
With filtered data only, you can approximate the soup by the pooled
expression of all cells, as SoupX does in the same situation:

```r
pooled <- Matrix::Matrix(Matrix::rowSums(sim$observed), ncol = 1, sparse = TRUE,
                         dimnames = list(rownames(sim$observed), "pooled"))
res_pooled <- decant(sim$observed, pooled, clusters = sim$labels)
```

- `pooled` is a one-column matrix holding the summed counts of all cells.
  Decant treats it as a single, very deep empty droplet.
- This soup is biased towards the abundant cell types, so treat the result
  as approximate.

**How many empty droplets do I need?**
The soup should rest on at least about 50,000 UMIs (`soup_umis` in the
`split_droplets()` summary). A standard 10x run has far more.

**Should I correct before or after quality control?**
Remove empty and damaged barcodes first (cell calling), then run Decant,
then do the rest of the analysis (normalisation, clustering, markers) on the
corrected counts.

**Can I use Decant with single-nucleus data?**
Yes. Contamination is usually higher, which makes correction more
important. For velocity, see [Chapter 10](10-Splice-Aware-and-Genotype-Modules).

**Why are some counts still present for marker genes in the wrong cell
type?**
Decant removes what ambient RNA can explain. Counts far above that level,
often from doublets, are left in place on purpose: they are not ambient RNA
and should be handled by doublet detection.

---

**Previous:** [Chapter 10, Splice-Aware and Genotype Modules](10-Splice-Aware-and-Genotype-Modules) | **Next:** [Chapter 12, Validation and Benchmarks](12-Validation-and-Benchmarks)
