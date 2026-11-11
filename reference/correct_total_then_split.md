# Baseline to beat: correct the TOTAL, then split the removal by each cell's own observed unspliced/spliced ratio (what you get if you decontaminate total and naively apportion to layers).

Baseline to beat: correct the TOTAL, then split the removal by each
cell's own observed unspliced/spliced ratio (what you get if you
decontaminate total and naively apportion to layers).

## Usage

``` r
correct_total_then_split(obs_u, obs_s, emp_u, emp_s, rho)
```

## Arguments

- obs_u, obs_s:

  genes x cells observed unspliced / spliced counts.

- emp_u, emp_s:

  genes x droplets empty-droplet unspliced / spliced counts.

- rho:

  per-cell contamination fraction (estimate upstream).

## Value

list with \$unspliced, \$spliced.
