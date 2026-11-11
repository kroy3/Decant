# GATE: false-positive rate on ambient-trap genes (should be low for Decant, inflated for naive) while preserving power on true_de genes. Also reports the NULL-gene false-positive rate: if that is not near alpha the simulator or normalisation is broken and the gate result is meaningless.

GATE: false-positive rate on ambient-trap genes (should be low for
Decant, inflated for naive) while preserving power on true_de genes.
Also reports the NULL-gene false-positive rate: if that is not near
alpha the simulator or normalisation is broken and the gate result is
meaningless.

## Usage

``` r
gate_de(seeds = 1:4, alpha = 0.05)
```

## Arguments

- seeds:

  replicate seeds.

- alpha:

  significance level.
