# Chapter 3. Quick Start

This chapter runs a complete analysis on simulated data in which the true
contamination is known, so you can see Decant recover it.

## 3.1 Simulate an experiment

```r
library(Decant)
sim <- simulate_experiment(seed = 1)
dim(sim$observed)
#> [1]  600 1200
dim(sim$empty)
#> [1]  600 4000
```

- `simulate_experiment(seed = 1)` creates a contaminated experiment with
  600 genes, 1,200 cells of six cell types, and 4,000 empty droplets. The
  `seed` makes the result identical every time.
- `sim$observed` holds the contaminated cell counts (genes × cells).
- `sim$empty` holds the empty droplets (genes × droplets).
- The simulation also stores the answer: `sim$truth` (the clean counts),
  `sim$rho_true` (each cell's true contamination) and `sim$labels` (the
  cell types).

## 3.2 Run Decant

```r
res <- decant(sim$observed, sim$empty, clusters = sim$labels)
res
#> <decant> 600 genes x 1200 cells
#>   modules : cluster_rho, correction_bayes_cell
#>   rho     : median 0.199 (10%-90%: 0.117-0.305)
```

- `decant()` takes the cell counts, the empty droplets and the cluster
  labels, and returns a result object.
- Printing the result shows its size, the steps that ran (`modules`) and a
  summary of the estimated contamination: the median cell is 19.9%
  ambient, and 80% of cells lie between 11.7% and 30.5%.

## 3.3 Compare with the truth

```r
summary(res$rho)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02389 0.14950 0.19933 0.20574 0.25223 0.55303 
summary(sim$rho_true)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0148  0.1371  0.1907  0.1970  0.2473  0.5317 
cor(res$rho, sim$rho_true)
#> [1] 0.9760583
```

- `res$rho` is the estimated contamination of every cell.
- The estimates sit slightly above the truth here (median 0.199 against
  0.191), and the correlation of 0.98 shows Decant ranks cells correctly.

## 3.4 Score the correction

```r
score_correction(sim$observed, sim$observed, sim$truth)$l1_error
#> [1] 0.2451949
score_correction(sim$observed, res$corrected, sim$truth)$l1_error
#> [1] 0.04877535
```

- `score_correction()` compares a matrix with the clean truth. Its
  `l1_error` adds up every count that is wrong, whether contamination left
  behind or real signal removed, as a fraction of the true total.
- Without correction (first line) 24.5% of counts are wrong. After Decant
  the error falls to 4.9%.

## 3.5 Use the corrected counts

```r
corrected <- res$corrected
all(corrected <= sim$observed)
#> [1] TRUE
```

- `res$corrected` has the same genes, cells and class as the input and can
  be used wherever you would use raw counts.
- The check confirms that no count was increased.

## 3.6 With your own data

```r
# Requires your own data
library(Decant)
raw   <- read_10x_counts("sample/outs/raw_feature_bc_matrix")
cells <- colnames(read_10x_counts("sample/outs/filtered_feature_bc_matrix"))
d     <- split_droplets(raw, cells = cells)
res   <- decant(d$cells, d$empties, clusters = my_clusters)
```

- `read_10x_counts()` reads a Cell Ranger folder into a sparse matrix.
- The filtered matrix is read only for its column names, the barcodes Cell
  Ranger called as cells.
- `split_droplets()` separates the raw matrix into cells and empty droplets.
- `my_clusters` is a vector with one cluster label per cell, in the same
  order as the columns of `d$cells`, for example from Seurat.

[Chapter 4](04-Preparing-Your-Data) explains each of these steps.

---

**Previous:** [Chapter 2, Key Concepts](02-Key-Concepts) | **Next:** [Chapter 4, Preparing Your Data](04-Preparing-Your-Data)
