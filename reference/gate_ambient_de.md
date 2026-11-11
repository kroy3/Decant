# GATE: does contamination-aware DE stop ambient leakage while staying calibrated and powerful?

Simulated 6 vs 6 case/control studies in which disease tissue carries
more ambient RNA. Tests DE in a rare cell type; the abundant type's
disease genes leak in through the soup. Pass requires, averaged over
seeds, in every case: at most 1 leaked gene called, null false-positive
rate \<= 0.075, false discovery proportion \<= 0.15, and power within
0.1 of naive edgeR (including on genes truly DE in both types).

## Usage

``` r
gate_ambient_de(seeds = 1:2, verbose = TRUE)
```

## Arguments

- seeds:

  replicate seeds.

- verbose:

  print the table.

## Value

invisibly, a data.frame with attribute `"pass"`.
