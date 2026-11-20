# Chapter 7. Estimation and Correction

`decant()` runs three steps in order. Each is also available as its own
function, which is useful when you want to inspect or replace one of them.

| Step | Function | Output |
|---|---|---|
| 1. Soup profile | `ambient_global()` | Fraction of ambient RNA per gene |
| 2. Contamination | `estimate_rho_cluster()` | rho per cell, with per-cluster diagnostics |
| 3. Correction | `correct_counts()` | Corrected counts |

## 7.1 Step by step

```r
library(Decant)
sim <- simulate_experiment(seed = 1)

soup <- ambient_global(sim$empty)
round(sum(soup), 6)
#> [1] 1

rho <- estimate_rho_cluster(sim$observed, soup, clusters = sim$labels)
round(median(rho), 3)
#> [1] 0.199

corrected <- correct_counts(sim$observed, soup, rho,
                            method = "bayes_cell", clusters = sim$labels)
```

- `ambient_global()` sums the empty droplets gene by gene and divides by
  the total, giving a soup profile that adds up to 1.
- `estimate_rho_cluster()` finds the genes each cluster does not express
  and measures how much soup explains their counts
  ([Chapter 2.5](02-Key-Concepts)). It returns one value per cell.
- `correct_counts()` removes `rho × cell total × soup` from each cell,
  using the rule given in `method`.

## 7.2 How rho is estimated

For each cluster, Decant repeats two steps until they agree:

1. Using the current rho, test every gene: is its count in the cluster
   clearly above what the soup alone would give? Genes that are not are
   called **absent**.
2. Re-estimate rho from the absent genes only: the ratio of their observed
   counts to their expected soup counts.

Five safeguards make the estimate reliable on real data:

- **Gene-to-gene spread is allowed for.** In real data, absent genes do not
  sit exactly on the line `rho × soup`: they scatter around it by about
  25-30%, because soup and capture differ slightly from gene to gene. Step 1
  therefore uses a negative binomial test whose extra spread is measured from
  the lowest absent genes, which expressed genes can never reach. A plain
  Poisson test would reject the absent genes that fall above the line and
  keep those below, and the estimate would settle at about half the true
  value. With no extra spread the test is the Poisson test.
- **Enough expected counts.** A gene enters the absent set only if enough
  ambient counts are expected for it (`min_expected = 5`). Without this,
  thousands of rarely detected genes look absent simply for lack of counts,
  and rho is overestimated.
- **Per-cell shrinkage.** Each cell's rho is pulled towards its cluster's
  rho, more strongly for cells with few counts. Small cells therefore do not
  get extreme values.
- **Weak clusters borrow.** A cluster without enough absent genes
  (`low_power`) borrows rho from the other clusters instead of guessing.
- **Implausible values are not trusted.** A cluster whose estimate exceeds
  0.5 (`implausible`) would be mostly soup, which happens when the soup looks
  like the cluster's own expression. It borrows rho too. If most cells are in
  such clusters, rho cannot be estimated from the data at all, and `decant()`
  warns ([Chapter 11.4](11-Diagnostics-and-Troubleshooting)).

## 7.3 The five correction rules

| Rule | How it removes contamination | When to use |
|---|---|---|
| `bayes_absent` (default) | First decides, gene by gene and cluster by cluster, whether the cluster expresses the gene at all; a gene it does not express loses all its counts. Otherwise as `bayes_cell`. Where the two groups of genes overlap, it does not decide and behaves as `bayes_cell` | General use; removes the most |
| `bayes_cell` | Assigns each count to cell or soup by how well each explains it, cell by cell. Counts far too high to be ambient keep all but their ambient expectation | Keeping the most native signal |
| `bayes` | The same, but one ratio for the whole cluster | Reproducing earlier results |
| `redistribute` | Removes exactly `rho × total` counts per cell, in soup proportion (SoupX's approach) | Comparing with SoupX |
| `subtract` | Subtracts the expectation gene by gene and stops at zero | Simplest reference; removes too little |

```r
for (m in c("subtract", "redistribute", "bayes", "bayes_cell", "bayes_absent")) {
  r <- decant(sim$observed, sim$empty, clusters = sim$labels, correction = m)
  e <- score_correction(sim$observed, r$corrected, sim$truth)$l1_error
  cat(sprintf("%-13s %.3f\n", m, e))
}
#> subtract      0.067
#> redistribute  0.063
#> bayes         0.047
#> bayes_cell    0.049
#> bayes_absent  0.048
```

- The loop runs Decant once per rule and prints the error against the
  truth (lower is better; no correction scores 0.245).
- On this one simulation the three Bayesian rules are close. Across all
  simulated settings in `gate_correction()`, `bayes_absent` has the lowest
  average error, by a small margin over `bayes_cell`.
- The difference shows on real mixtures of human and mouse cells, where the
  truth is known. On two 10x mixtures never used to design it,
  `bayes_absent` removed 0.891 of the ambient signal (`bayes_cell` 0.531,
  DecontX 0.655), and kept 0.9918 of the cells' own signal (`bayes_cell`
  0.9945). On the hardest mixture it kept 2.5 points less than
  `bayes_cell`. That is why it is the default, and why `bayes_cell` stays
  available: use `correction = "bayes_cell"` when keeping every native
  count matters more than removing the soup.
- On real blood data (PBMC 4k) the clusters share most of their genes, the
  two groups overlap, and `bayes_absent` gives the same result as
  `bayes_cell`. Both keep the counts of likely doublets, which `bayes`
  removes.

## 7.4 The guarantee

Every rule satisfies `0 ≤ corrected ≤ observed` for every gene in every
cell. Decant checks this after each run and stops with an error if it were
ever violated. Some other tools return corrected counts that are larger
than the observed ones; Decant cannot.

## 7.5 Tuning parameters

The defaults were chosen on test data and rarely need changing.

| Function | Argument | Default | Effect |
|---|---|---|---|
| `estimate_rho_cluster()` | `alpha` | `0.01` | Stringency for calling a gene expressed; smaller calls fewer genes expressed |
| `estimate_rho_cluster()` | `min_expected` | `5` | Minimum expected ambient count for an absent gene |
| `estimate_rho_cluster()` | `floor_rho` | `0.005` | Smallest contamination considered possible; also the lowest starting value |
| `estimate_rho_cluster()` | `max_rho` | `0.5` | Largest contamination treated as identifiable |
| `estimate_rho_cluster()` | `overdispersed` | `TRUE` | Allow for gene-to-gene spread of absent genes; `FALSE` gives the earlier Poisson test |
| `correct_counts()` | `cap_q` | `0.001` | For `bayes_cell` and `bayes_absent`: how unusual a count must be before it is treated as not ambient |

---

**Previous:** [Chapter 6, Seurat Integration](06-Seurat-Integration) | **Next:** [Chapter 8, Multi-Sample Experiments](08-Multi-Sample-Experiments)
