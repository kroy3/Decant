# Contamination-aware differential expression

## The problem

In a case/control single-cell study, samples differ in how much ambient
RNA they carry. Fragile disease tissue lyses more, so its soup is larger
and is enriched for the disease genes of whichever cell type is
abundant. Those reads land in every droplet. A pseudobulk test on a rare
cell type then reports the abundant type’s disease genes as
differentially expressed in the rare type.

[`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md)
does not subtract the soup. It puts each sample’s measured soup inside
the model as a known additive term:

``` math
y_{gs} \sim \mathrm{NB}\big(N_s \exp(x_s^\top \beta_g) + \rho_s T_s\,\mathrm{soup}_{gs}\big)
```

and tests $`\beta_g`$ with a quasi-likelihood F-test.

## A simulated study with known truth

[`simulate_ambient_study()`](https://kroy3.github.io/Decant/reference/simulate_ambient_study.md)
makes 6 control and 6 disease samples. Cell type `A` is abundant and has
80 disease genes. Cell type `R` is rare and has 40 disease genes of its
own. Disease samples carry about 12% ambient RNA, controls about 5%.

``` r

library(Decant)
st <- simulate_ambient_study(seed = 3, G = 600, n_cells = 500)
table(st$cond)
#> 
#> ctrl  dis 
#>    6    6
round(tapply(st$rho_true, st$cond, mean), 3)
#>  ctrl   dis 
#> 0.052 0.111
```

Pseudobulk the rare type, one column per sample:

``` r

is_R <- st$label == "R"
samples <- sort(unique(st$sample))
Y <- vapply(samples, function(s)
  rowSums(st$X[, is_R & st$sample == s, drop = FALSE]), numeric(nrow(st$X)))
rownames(Y) <- rownames(st$X)
dim(Y)
#> [1] 600  12
```

## Run the test

`empties` is a list with one empty-droplet matrix per sample, in the
same order as the columns of `Y`.

``` r

design <- model.matrix(~ st$cond)
res <- ambient_de(Y, st$empties, design)
head(res[, c("gene", "logFC", "PValue", "FDR", "ambient_frac", "naive_FDR", "ambient_driven")])
#>     gene     logFC       PValue          FDR ambient_frac    naive_FDR
#> 117 g218  1.462029 2.746077e-10 8.952210e-08   0.02752370 2.085851e-07
#> 173 g309  1.222120 6.475971e-09 1.055583e-06   0.01134446 7.910796e-06
#> 189 g338  1.183664 3.708674e-08 4.030093e-06   0.01707288 7.910796e-06
#> 84  g158 -1.092676 1.740945e-07 1.418870e-05   0.01294913 1.077668e-05
#> 99  g184 -1.182222 4.557539e-07 2.971515e-05   0.01613194 7.910796e-06
#> 47  g121 -1.033911 6.203022e-07 3.069116e-05   0.01149664 1.295151e-05
#>     ambient_driven
#> 117          FALSE
#> 173          FALSE
#> 189          FALSE
#> 84           FALSE
#> 99           FALSE
#> 47           FALSE
```

`naive_FDR` is standard edgeR on the same pseudobulk. `ambient_driven`
marks genes that edgeR calls and the contamination-aware test does not.

## Score against the truth

``` r

gi <- as.integer(sub("g", "", res$gene))
leaked <- gi %in% st$de_abund   # A's disease genes: any call in R is false
true_de <- gi %in% st$de_rare   # R's own disease genes
c(leaked_called_naive  = sum(res$naive_FDR[leaked] < 0.1),
  leaked_called_decant = sum(res$FDR[leaked] < 0.1),
  power_naive  = mean(res$naive_FDR[true_de] < 0.1),
  power_decant = mean(res$FDR[true_de] < 0.1))
#>  leaked_called_naive leaked_called_decant          power_naive 
#>            0.0000000            0.0000000            0.8333333 
#>         power_decant 
#>            0.8611111
```

The per-sample ambient fractions the test estimated:

``` r

attr(res, "samples")
#>    sample ambient_frac soup_umis norm_factor
#> 1       1   0.05691182     49848   0.9663453
#> 2       2   0.05807062     49966   1.0214628
#> 3       3   0.04450713     50234   0.9746262
#> 4       4   0.03639673     49819   1.0010121
#> 5       5   0.07353880     50063   1.0171349
#> 6       6   0.03910505     50089   0.9802863
#> 7       7   0.08618538     49997   1.0165516
#> 8       8   0.09673501     50090   0.9952828
#> 9       9   0.14177055     49851   0.9904143
#> 10     10   0.09325414     49814   0.9953381
#> 11     11   0.13165005     49671   1.0026047
#> 12     12   0.11165293     50086   1.0414611
```

## On a Seurat object

[`DecantDE()`](https://kroy3.github.io/Decant/reference/DecantDE.md)
does the pseudobulking and reads each donor’s raw matrix:

``` r

res <- DecantDE(obj, cell_type = "Microglia", sample_col = "donor",
                formula = ~ diagnosis, raw = raw_paths_by_donor)
subset(res, ambient_driven)
```

## What this does not show

- The simulator is built around the additive-soup mechanism the test
  assumes. This vignette shows the method does what it is designed to
  do. It is not real-data evidence.
- Correcting counts first and then testing is not always worse. In
  `validation/de_comparators.R`, official SoupX followed by edgeR
  removes almost all leakage on this simulator, at some cost in power.
  See the README for the table.
- With very few empty droplets, or three samples per group, the test
  leaks or loses power. The README lists the measured limits.

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Decant_0.4.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.52        rlang_1.3.0       xfun_0.62        
#>  [5] otel_0.3.0        textshaping_1.0.5 jsonlite_2.0.0    statmod_1.5.2    
#>  [9] htmltools_0.5.9   ragg_1.5.2        sass_0.4.10       locfit_1.5-9.12  
#> [13] rmarkdown_2.32    grid_4.6.2        evaluate_1.0.5    jquerylib_0.1.4  
#> [17] fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5   compiler_4.6.2   
#> [21] fs_2.1.0          limma_3.68.5      edgeR_4.10.7      lattice_0.23-1   
#> [25] systemfonts_1.3.2 digest_0.6.39     R6_2.6.1          splines_4.6.2    
#> [29] Matrix_1.7-6      bslib_0.12.0      tools_4.6.2       pkgdown_2.2.1    
#> [33] cachem_1.1.0      desc_1.4.3
```
