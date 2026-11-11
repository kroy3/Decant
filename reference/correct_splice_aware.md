# Splice-aware layer decontamination.

Splice-aware layer decontamination.

## Usage

``` r
correct_splice_aware(obs_u, obs_s, emp_u, emp_s, rho, drop_zeros = TRUE)
```

## Arguments

- obs_u, obs_s:

  genes x cells observed unspliced / spliced counts.

- emp_u, emp_s:

  genes x droplets empty-droplet unspliced / spliced counts.

- rho:

  per-cell contamination fraction (estimate upstream).

- drop_zeros:

  for sparse input, drop entries corrected to exactly 0.

## Value

list with \$unspliced, \$spliced (corrected, mass-conserving) and the
per-gene ambient unspliced fraction used.
