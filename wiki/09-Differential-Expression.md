# Chapter 9. Differential Expression

Correcting counts and then running a standard test is one way to deal with
ambient RNA in differential expression. Decant offers a second way:
`ambient_de()` keeps the raw counts and models each sample's soup inside the
test.

## 9.1 The problem

Suppose disease samples carry more ambient RNA than control samples, and
that ambient RNA is rich in the disease genes of an abundant cell type, A. A
pseudobulk test on a rare cell type, R, sees those genes rise in R's disease
samples and reports them as differentially expressed in R. They are not.
Their rise is entirely due to the soup.

## 9.2 The model

For each gene `g` and sample `s`, `ambient_de()` fits

```
counts[g, s]  ~  NegativeBinomial( native expression[g, s]  +  rho[s] × total[s] × soup[g, s] )
```

- The **native expression** term depends on the design (for example
  control versus disease) and is what is tested.
- The **ambient** term is measured: `soup` from the sample's own empty
  droplets and `rho` estimated from the data. It is added to the mean, not
  subtracted from the counts, and its uncertainty is carried into the test.
- The test is a quasi-likelihood F-test, as in edgeR.

## 9.3 Example

```r
library(Decant)
st <- simulate_ambient_study(seed = 3, G = 600, n_cells = 500)
table(st$cond)
#>
#> ctrl  dis
#>    6    6
```

- `simulate_ambient_study()` creates 6 control and 6 disease samples. Cell
  type A is abundant and has 80 disease genes; cell type R is rare and has
  40 of its own. Disease samples carry more ambient RNA than controls.

### Pseudobulk the cell type of interest

```r
is_R <- st$label == "R"
samples <- sort(unique(st$sample))
Y <- vapply(samples, function(s)
  rowSums(st$X[, is_R & st$sample == s, drop = FALSE]), numeric(nrow(st$X)))
rownames(Y) <- rownames(st$X)
dim(Y)
#> [1] 600  12
```

- `is_R` selects the cells of the rare type.
- `vapply()` sums those cells' counts within each sample, giving one column
  per sample: the pseudobulk matrix `Y`, genes × samples.

### Run the test

```r
design <- model.matrix(~ cond, data = data.frame(cond = st$cond))
res <- ambient_de(Y, st$empties, design)
head(res[, c("gene", "logFC", "FDR", "ambient_frac", "naive_FDR", "ambient_driven")], 3)
#>     gene    logFC          FDR ambient_frac    naive_FDR ambient_driven
#> 117 g218 1.461912 1.129595e-07   0.02752255 3.990660e-07          FALSE
#> 173 g309 1.226158 2.978136e-07   0.01134339 1.797587e-06          FALSE
#> 189 g338 1.183623 4.248957e-06   0.01707194 9.899770e-06          FALSE
```

- `design` is a standard model matrix with one row per sample. The last
  column (`conddis`, disease versus control) is tested by default.
- `st$empties` is a list of empty-droplet matrices, one per sample, in the
  same order as the columns of `Y`.
- Results are ordered by p-value. The key columns are:

| Column | Meaning |
|---|---|
| `logFC` | log2 fold change of native expression |
| `PValue`, `FDR` | Contamination-aware p-value and false discovery rate |
| `ambient_frac` | Average share of the gene's counts that are ambient |
| `naive_FDR` | FDR of a standard edgeR test on the same pseudobulk |
| `ambient_driven` | `TRUE` if the standard test calls the gene and the contamination-aware test does not |

### Genes explained by the soup

```r
sum(res$ambient_driven)
#> [1] 12
```

- Twelve genes would be reported by a standard pipeline but are explained
  by differences in the soup between conditions. These are the calls
  `ambient_de()` protects you from.

### Per-sample estimates

```r
head(attr(res, "samples"), 3)
#>   sample ambient_frac soup_umis norm_factor
#> 1      1   0.05692409     49848   0.9663427
#> 2      2   0.05808700     49966   1.0214589
#> 3      3   0.04450122     50234   0.9746253
```

- `ambient_frac` is the estimated ambient fraction of the cell type in each
  sample, and `soup_umis` the number of UMIs the sample's soup rests on.

## 9.4 From a Seurat object

```r
# Requires your own data
res <- DecantDE(obj, cell_type = "Microglia", sample_col = "donor",
                formula = ~ diagnosis, raw = raw_paths_by_donor)
subset(res, ambient_driven)
```

- `DecantDE()` builds the pseudobulk, reads each sample's empty droplets and
  calls `ambient_de()`.
- `cell_type` is the identity to test (`Idents(obj)` by default, or a
  column given in `celltype_col`).
- `sample_col` names the sample of each cell; `formula` is the
  sample-level design. Every variable in `formula` must be constant within
  a sample.
- `raw` is a named list with one raw matrix or Cell Ranger path per sample.

## 9.5 Useful options

| Argument | Default | Purpose |
|---|---|---|
| `coef` | last column | Design column to test |
| `robust` | `TRUE` | Down-weight outlier samples gene by gene; important in large cohorts |
| `lfc` | `0` | Test against a minimum fold change, for example `log2(1.1)` |
| `max_ambient_frac` | `0.75` | Skip genes whose counts are mostly ambient; they are listed in `attr(res, "ambient_dominated")` |
| `naive` | `TRUE` | Also run the standard test, to fill `naive_FDR` and `ambient_driven` |
| `joint_rho` | `TRUE` | Refine each sample's ambient fraction on genes absent from all samples pooled. Leave it on: without it, the ambient fraction of a single cell type can be wildly wrong |
| `dispersion` | `"trended"` | NB dispersion in the fit. `"tagwise"` was the earlier default and is not calibrated on real data (9.6) |

## 9.6 Limits

- **Calibration on real data.** On two real single-nucleus datasets
  (GSE272082, Alzheimer's vs control; GSE190015, nine controls), condition
  labels were shuffled between donors 24 times, so there was nothing to
  find. In microglia and oligodendrocytes, `ambient_de()` reported a gene in
  9 of the 96 tests, never more than 9 genes; naive edgeR in 4. Earlier
  versions used `dispersion = "tagwise"` and reported genes in 82 of 96
  tests, up to 574; rerun analyses made with them. To check any method on
  your own data, shuffle your condition labels between samples and count
  the hits.
- On real data only calibration has been checked. That `ambient_de()`
  removes ambient-driven hits is shown in simulation.
- With very few empty droplets per sample (a few hundred), the soup is
  uncertain and some leakage returns.
- With three samples per group, power is lower than with six.
- Correcting counts with a good tool and then testing performs
  similarly on the same simulations; `ambient_de()` has the advantage when
  the design includes covariates such as batch.

---

**Previous:** [Chapter 8, Multi-Sample Experiments](08-Multi-Sample-Experiments) | **Next:** [Chapter 10, Splice-Aware and Genotype Modules](10-Splice-Aware-and-Genotype-Modules)
