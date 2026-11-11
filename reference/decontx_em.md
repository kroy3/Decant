# DecontX-style decontamination (Bayesian-mixture EM, reimplementation).

Model (Yang et al. 2020): each cell belongs to a cluster k. Its observed
counts are a mixture of native expression phi_k and a contamination
distribution eta_k (expression originating from OTHER clusters in the
experiment). Per-cell contamination theta_c and the profiles are fit by
EM. Crucially, this uses NO empty droplets – contamination is inferred
from population structure alone. Returns a decontaminated (native) count
matrix.

## Usage

``` r
decontx_em(observed, labels, iters = 30, eps = 1e-10)
```

## Arguments

- observed:

  genes x cells counts.

- labels:

  integer cluster labels per cell (length = ncol).

- iters:

  EM iterations.

- eps:

  numerical floor.

## Value

list with \$corrected (genes x cells) and \$theta (per-cell
contamination).
