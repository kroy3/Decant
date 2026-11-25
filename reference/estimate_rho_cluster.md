# Cluster-pooled, absent-gene rho estimator with per-cell shrinkage.

Cluster-pooled, absent-gene rho estimator with per-cell shrinkage.

## Usage

``` r
estimate_rho_cluster(
  observed,
  soup,
  clusters,
  alpha = 0.01,
  min_absent_soup = 0.02,
  max_iter = 100,
  min_expected = 5,
  floor_rho = 0.005,
  overdispersed = TRUE,
  max_rho = 0.5,
  weak_cosine = 0.9,
  weak_kappa = 30,
  weak_ci_width = 0.05,
  weak_rule = c("evidence", "any")
)
```

## Arguments

- observed:

  genes x cells counts (base matrix or dgCMatrix).

- soup:

  length-G ambient profile (e.g. from
  [`ambient_global()`](https://kroy3.github.io/Decant/reference/ambient_global.md)).

- clusters:

  length-cells cluster labels. Prefer a moderate number of reasonably
  large clusters: very small clusters lose the power to reject
  weakly-expressed genes and bias rho upward.

- alpha:

  one-sided Poisson level for calling a gene "above soup" (i.e. natively
  expressed) in a cluster.

- min_absent_soup:

  a cluster whose absent set covers less than this fraction of soup mass
  is flagged `low_power`.

- max_iter:

  iteration cap for the absent-set fixed point.

- min_expected:

  a gene can join the absent set only if its expected ambient count in
  the cluster, at the current rho, is at least this. Below it the
  Poisson test cannot tell a weakly expressed gene from an absent one,
  and on real data (tens of thousands of low-count genes) letting those
  genes in drove rho to several times its true value.

- floor_rho:

  the smallest contamination considered possible. Genes whose total
  count in the cells is significantly below this fraction of their soup
  expectation are treated as soup-only and dropped from the soup.

- overdispersed:

  test genes against a negative binomial whose extra gene-to-gene
  variance is estimated from the data, instead of a Poisson. On real
  data, genes a cluster does not express scatter around rho by about
  25-30% (soup and capture differ slightly gene by gene), far more than
  Poisson noise. A Poisson test then rejects the absent genes that fall
  above the current rho and keeps those below, and the fixed point
  settles near the bottom of their spread: on a human/mouse mixture this
  halved rho. The variance is estimated from the lower tail of the
  absent genes, which cannot contain expressed genes. With no extra
  variance the test reduces to the Poisson one.

- max_rho:

  the largest contamination treated as identifiable. A cluster whose
  estimate exceeds it is flagged (`implausible`) and borrows rho from
  the other clusters, like a low-power cluster. Estimates above about
  0.5 mean the cells would be mostly soup, which points to a soup that
  resembles the cells (a single cell type, or a cluster that dominates
  the soup) rather than to real contamination: on a mouse-only cell line
  with 0.4% true contamination the estimator returned 0.5 to 0.94.

- weak_cosine, weak_kappa:

  thresholds above which a cluster's profile counts as close to the soup
  (`collinear`): the weighted cosine between the cluster's raw profile
  and the soup (0.9: a two-profile fit's variance inflation 1 / (1 -
  c^2) then exceeds 5) and the Gram condition number (30, i.e. a cosine
  of 0.935, so it only matters if `weak_cosine` is raised).

- weak_ci_width:

  a cluster whose rho has a 95% interval wider than this is weakly
  identified (0.05: wider than typical contamination itself, 0.01-0.1).

- weak_rule:

  "evidence" (default): a cluster is `weak` when the estimator's own
  evidence is thin (interval wider than `weak_ci_width`, almost no
  absent genes, or an implausible own estimate). "any": also when it is
  `collinear`. A collinear cluster whose absent genes still cover enough
  soup is usually estimated well, so "any" raises false alarms on
  soup-dominant cell types.

## Value

numeric vector of per-cell rho, with attribute `"clusters"`: a
data.frame of per-cluster diagnostics (n_cells, rho, n_absent_genes,
absent_soup_frac, absent_counts, low_power, borrowed, rho_two_sided,
dispersion: the extra squared coefficient of variation used,
implausible: rho above `max_rho`; identifiability: cosine, kappa, js,
fisher_umi, n_eff, rho_se, ci_width, collinear, and weak / weak_reason /
action_taken), attribute `"unidentifiable"` (no cluster had usable
evidence, or more than half the cells are in implausible clusters; rho
should then be supplied), and attributes `"soup_mismatch"` (logical),
`"soup_only"` (logical per gene) and `"soup_only_frac"` (their share of
the soup). Low-power clusters take the evidence-weighted rho of the
other clusters (`borrowed = TRUE`). `soup_mismatch = TRUE` means more
than 5% of the soup was in genes no cell contains; those genes are
dropped, but the empties deserve a look.
