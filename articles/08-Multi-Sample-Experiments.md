# Chapter 8. Multi-Sample Experiments

Each sample (each 10x channel) has its own suspension and therefore its
own soup. When cells from several samples are analysed together, Decant
gives every sample its own soup profile.

## 8.1 Why samples need separate soups

The soup reflects the cells that lysed in that particular sample. Two
donors, or a control and a treated sample, can have quite different
soups. Using one pooled soup for all of them over-corrects some genes in
one sample and under-corrects them in another.

At the same time, a sample with few empty droplets gives a noisy soup.
Decant handles both problems with **partial pooling**: a sample with
many empty droplets keeps its own soup, and a sample with few is pulled
towards the average of all samples.

## 8.2 Running several samples

``` r

library(Decant)
mm <- simulate_multimodal(seed = 1)
empties <- lapply(mm$empties, function(e) e$unspliced + e$spliced)
sapply(empties, ncol)
#> [1] 3000  300 3000  300 3000  300

res <- decant(mm$obs_total, empties,
              sample_of = mm$sample_of, clusters = mm$labels)
res
#> <decant> 600 genes x 1500 cells
#>   modules : hierarchical_ambient, cluster_rho, correction_bayes_absent
#>   rho     : median 0.203 (10%-90%: 0.115-0.315)
```

- [`simulate_multimodal()`](https://kroy3.github.io/Decant/reference/simulate_multimodal.md)
  simulates six samples. Half have 3,000 empty droplets, half only 300.
- The simulation stores spliced and unspliced counts separately; the
  second line adds them to get total counts per empty droplet.
- `empties` is a **list with one matrix per sample**. Passing a list is
  what switches on per-sample soups (`hierarchical_ambient` in
  `modules`).
- `sample_of` gives the sample of each cell, as a number matching the
  position in the list, or as a name matching `names(empties)`.

## 8.3 Per-sample results

``` r

dim(res$ambient)
#> [1] 600   6
round(rbind(estimated = tapply(res$rho, mm$sample_of, median),
            true      = tapply(mm$rho_true, mm$sample_of, median)), 3)
#>               1     2     3     4     5     6
#> estimated 0.208 0.214 0.200 0.194 0.204 0.204
#> true      0.197 0.195 0.186 0.183 0.191 0.195
```

- `res$ambient` now has one soup column per sample.
- The table compares the median estimated and true contamination in each
  sample. The estimates are one to two percentage points above the
  truth, and no worse for the samples with only 300 empty droplets.

## 8.4 The pooled soups directly

``` r

pooled <- pool_soup(empties)
names(pooled)
#> [1] "pooled"      "independent" "lambda"      "mass"
pooled$mass
#> [1] 75064  7432 74846  7554 75017  7575
round(pooled$lambda, 2)
#> [1] 0.65 0.15 0.64 0.15 0.65 0.16
```

- [`pool_soup()`](https://kroy3.github.io/Decant/reference/pool_soup.md)
  is the function
  [`decant()`](https://kroy3.github.io/Decant/reference/decant.md) uses
  internally.
- `independent` holds each sample’s own soup and `pooled` the soups
  after partial pooling.
- `mass` is the number of UMIs in each sample’s empty droplets.
- `lambda` is the weight each sample’s own soup receives; the rest comes
  from the average of all samples. Samples with ten times fewer soup
  UMIs (2, 4 and 6) rely on their own soup for only 15% and borrow the
  rest.

## 8.5 With Seurat

[`RunDecant()`](https://kroy3.github.io/Decant/reference/RunDecant.md)
does the same when given one raw matrix per sample and a sample column
([Chapter
6.5](https://kroy3.github.io/Decant/articles/06-Seurat-Integration.md)).

------------------------------------------------------------------------

**Previous:** [Chapter 7, Estimation and
Correction](https://kroy3.github.io/Decant/articles/07-Estimation-and-Correction.md)
\| **Next:** [Chapter 9, Differential
Expression](https://kroy3.github.io/Decant/articles/09-Differential-Expression.md)
