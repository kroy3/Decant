# Score a corrected matrix against ground truth.

Score a corrected matrix against ground truth.

## Usage

``` r
score_correction(observed, corrected, truth)
```

## Arguments

- observed, corrected, truth:

  genes x cells matrices (truth = clean own counts).

## Value

named list of scalar metrics. `l1_error` = sum\|corrected - truth\| /
sum(truth): residual contamination and destroyed signal on one scale, so
a method cannot look good by trading one for the other.
