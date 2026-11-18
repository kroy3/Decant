# Ambient-aware ligand-receptor inference between clusters ("chirp").

For each gene and cluster, compares the observed cluster total with the
count ambient RNA is expected to contribute,
`sum_c rho_c * T_c * soup_g` (per sample when `soup` has one column per
sample). A gene counts as natively expressed in a cluster when the
observed total is significantly above `margin` times that expectation
(Poisson upper tail, Benjamini-Hochberg across genes x clusters). Pair
scores use native expression, the observed total minus the expected
ambient count, per 10,000 native UMIs.

## Usage

``` r
chirp(
  observed,
  pairs,
  fit = NULL,
  soup = NULL,
  rho = NULL,
  clusters = NULL,
  sample_of = NULL,
  alpha = 0.05,
  margin = 2,
  min_native = 1
)
```

## Arguments

- observed:

  genes x cells raw counts (matrix or dgCMatrix), uncorrected.

- pairs:

  data.frame with columns `ligand` and `receptor` (gene names matching
  `rownames(observed)`). A multi-subunit receptor or ligand can be given
  as genes joined by `"+"`; the complex is scored by its least expressed
  subunit and called only if every subunit is called.

- fit:

  optional `decant` object; supplies `soup`, `rho` and `clusters` when
  those are not given.

- soup:

  ambient profile: a length-genes vector, or genes x samples matrix with
  `sample_of`.

- rho:

  per-cell (or scalar) contamination fraction.

- clusters:

  length-cells cluster labels.

- sample_of:

  length-cells sample index into the columns of `soup` (required when
  `soup` has more than one column).

- alpha:

  FDR level for the native-expression call.

- margin:

  multiplier on the expected ambient count before testing. The default,
  2, kept false calls at zero in simulations with rho underestimated by
  a third; it does not protect if rho is half its true value. Lower it
  for more sensitivity only if rho is trusted.

- min_native:

  minimum native expression, in counts per 10,000 native UMIs, for a
  call.

## Value

data.frame, one row per pair x sender x receiver, sorted by score:
`ligand`, `receptor`, `sender`, `receiver`, `score` (geometric mean of
native ligand in sender and native receptor in receiver), `call` (both
sides natively expressed), `naive_score` (same on raw counts),
`naive_call` (both sides detected in at least 10% of cells, the usual
rule), `ambient_share_ligand` (expected ambient share of the ligand
counts in the sender) and `p_ligand`, `p_receptor` (BH-adjusted). Rows
with `naive_call & !call` are the interactions ambient RNA would have
produced. Attribute `"native"`: genes x clusters native expression.

## Examples

``` r
sim <- simulate_experiment(n_genes = 200, n_types = 4, n_cells = 300, seed = 1)
g <- rownames(sim$observed)
pairs <- data.frame(ligand = g[1:5], receptor = g[6:10])
res <- chirp(sim$observed, pairs, soup = ambient_global(sim$empty),
             rho = sim$rho_true, clusters = sim$labels)
head(res[res$call, ])
#>   ligand receptor sender receiver    score call naive_score naive_call
#> 1  gene5   gene10      1        1 629.7999 TRUE    512.4324       TRUE
#> 2  gene2    gene7      1        1 612.6521 TRUE    498.0305       TRUE
#> 3  gene3    gene8      1        1 607.6516 TRUE    493.6902       TRUE
#> 4  gene1    gene6      1        1 606.9207 TRUE    493.0299       TRUE
#> 5  gene4    gene9      1        1 605.7259 TRUE    492.0313       TRUE
#> 6  gene5   gene10      2        1 100.9971 TRUE     83.2499       TRUE
#>   ambient_share_ligand      p_ligand p_receptor
#> 1         0.0007390784  0.000000e+00          0
#> 2         0.0001889548  0.000000e+00          0
#> 3         0.0003020906  0.000000e+00          0
#> 4         0.0004197138  0.000000e+00          0
#> 5         0.0002685611  0.000000e+00          0
#> 6         0.0282350389 1.544711e-255          0
```
