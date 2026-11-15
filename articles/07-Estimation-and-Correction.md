# Chapter 7. Estimation and Correction

[`decant()`](https://kroy3.github.io/Decant/reference/decant.md) runs
three steps in order. Each is also available as its own function, which
is useful when you want to inspect or replace one of them.

| Step | Function | Output |
|----|----|----|
| 1\. Soup profile | [`ambient_global()`](https://kroy3.github.io/Decant/reference/ambient_global.md) | Fraction of ambient RNA per gene |
| 2\. Contamination | [`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md) | rho per cell, with per-cluster diagnostics |
| 3\. Correction | [`correct_counts()`](https://kroy3.github.io/Decant/reference/correct_counts.md) | Corrected counts |

## 7.1 Step by step

``` r

library(Decant)
sim <- simulate_experiment(seed = 1)

soup <- ambient_global(sim$empty)
round(sum(soup), 6)
#> [1] 1

rho <- estimate_rho_cluster(sim$observed, soup, clusters = sim$labels)
round(median(rho), 3)
#> [1] 0.192

corrected <- correct_counts(sim$observed, soup, rho,
                            method = "bayes_cell", clusters = sim$labels)
```

- [`ambient_global()`](https://kroy3.github.io/Decant/reference/ambient_global.md)
  sums the empty droplets gene by gene and divides by the total, giving
  a soup profile that adds up to 1.
- [`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md)
  finds the genes each cluster does not express and measures how much
  soup explains their counts ([Chapter
  2.5](https://kroy3.github.io/Decant/articles/02-Key-Concepts.md)). It
  returns one value per cell.
- [`correct_counts()`](https://kroy3.github.io/Decant/reference/correct_counts.md)
  removes `rho × cell total × soup` from each cell, using the rule given
  in `method`.

## 7.2 How rho is estimated

For each cluster, Decant repeats two steps until they agree:

1.  Using the current rho, test every gene: is its count in the cluster
    clearly above what the soup alone would give? Genes that are not are
    called **absent**.
2.  Re-estimate rho from the absent genes only: the ratio of their
    observed counts to their expected soup counts.

Three safeguards make the estimate reliable on real data:

- A gene enters the absent set only if enough ambient counts are
  expected for it (`min_expected = 5`). Without this, thousands of
  rarely detected genes look absent simply for lack of counts, and rho
  is overestimated.
- Each cell’s rho is pulled towards its cluster’s rho, more strongly for
  cells with few counts. Small cells therefore do not get extreme
  values.
- A cluster without enough absent genes (`low_power`) borrows rho from
  the other clusters instead of guessing.

## 7.3 The four correction rules

| Rule | How it removes contamination | When to use |
|----|----|----|
| `bayes_cell` (default) | Assigns each count to cell or soup by how well each explains it, cell by cell. Counts far too high to be ambient keep all but their ambient expectation | General use |
| `bayes` | The same, but one ratio for the whole cluster | Reproducing earlier results |
| `redistribute` | Removes exactly `rho × total` counts per cell, in soup proportion (SoupX’s approach) | Comparing with SoupX |
| `subtract` | Subtracts the expectation gene by gene and stops at zero | Simplest reference; removes too little |

``` r

for (m in c("subtract", "redistribute", "bayes", "bayes_cell")) {
  r <- decant(sim$observed, sim$empty, clusters = sim$labels, correction = m)
  e <- score_correction(sim$observed, r$corrected, sim$truth)$l1_error
  cat(sprintf("%-13s %.3f\n", m, e))
}
#> subtract      0.070
#> redistribute  0.062
#> bayes         0.047
#> bayes_cell    0.049
```

- The loop runs Decant once per rule and prints the error against the
  truth (lower is better; no correction scores 0.245).
- On this one simulation `bayes` and `bayes_cell` are close. Across all
  simulated settings in
  [`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md),
  `bayes_cell` has the lowest average error. On real blood data it also
  keeps the counts of likely doublets, which `bayes` removes. That is
  why it is the default.

## 7.4 The guarantee

Every rule satisfies `0 ≤ corrected ≤ observed` for every gene in every
cell. Decant checks this after each run and stops with an error if it
were ever violated. Some other tools return corrected counts that are
larger than the observed ones; Decant cannot.

## 7.5 Tuning parameters

The defaults were chosen on test data and rarely need changing.

| Function | Argument | Default | Effect |
|----|----|----|----|
| [`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md) | `alpha` | `0.01` | Stringency for calling a gene expressed; smaller calls fewer genes expressed |
| [`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md) | `min_expected` | `5` | Minimum expected ambient count for an absent gene |
| [`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md) | `floor_rho` | `0.005` | Smallest contamination considered possible |
| [`correct_counts()`](https://kroy3.github.io/Decant/reference/correct_counts.md) | `cap_q` | `0.001` | For `bayes_cell`: how unusual a count must be before it is treated as not ambient |

------------------------------------------------------------------------

**Previous:** [Chapter 6, Seurat
Integration](https://kroy3.github.io/Decant/articles/06-Seurat-Integration.md)
\| **Next:** [Chapter 8, Multi-Sample
Experiments](https://kroy3.github.io/Decant/articles/08-Multi-Sample-Experiments.md)
