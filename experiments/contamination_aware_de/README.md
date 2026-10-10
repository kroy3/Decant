# Contamination-aware DE: falsification experiment

Question: when ambient RNA differs between samples, do standard pseudobulk DE
pipelines report an abundant cell type's disease genes as DE in a rare cell
type, and does modelling the measured soup inside the test fix it?

Design (`sim.R`): 6 control vs 6 disease samples, biological between-sample
variation (BCV 0.25), an abundant type A with 80 disease genes (3x), a rare
type R (10% of cells) with 40 true disease genes (2x), contamination ~5% in
controls vs ~12% in fragile disease tissue, a per-sample soup measured from
that sample's empty droplets.

Methods (`methods.R`), all testing type R:
- `naive`: edgeR QL on raw pseudobulk.
- `corrected`: Decant correction per sample, then edgeR QL.
- `caware`: negative binomial with mean `N_s * exp(b0 + b1*x_s) + E_s`, where
  `E_s = sum over R cells of rho_c * T_c` times the sample's soup profile is a
  KNOWN additive term (identity scale, not a log covariate). LRT on b1,
  dispersion borrowed from edgeR on the corrected counts.

Run from this directory: `Rscript run1.R`, `Rscript run2.R`
(set `DECANT_ROOT` if not run from inside the repo).

## Results (mean of 3 seeds)

`fp_leak` = fraction of A's disease genes called DE in R (all false).

| case | naive | corrected | caware |
|---|---|---|---|
| base | 1.00 | 0.99 | 0.00 |
| negative control (equal mean rho) | 0.59 | 0.18 | 0.00 |
| rho over-estimated 30% | 1.00 | 0.99 | 0.00 |
| rho under-estimated 30% | 1.00 | 0.99 | 0.52 |
| only 150 empty droplets | 1.00 | 0.99 | 0.46 |
| genes truly DE in both A and R | 1.00 | 0.98 | 0.02 (power on them 0.92) |

Other-null false-positive rate at p < 0.05: caware 0.045-0.09 vs naive
0.20-0.24. Power on true R genes: ~0.93 for all methods.

## What this does and does not show

- Shows: correcting counts first does NOT stop leakage (0.99); the residual
  after subtraction still tracks condition. Modelling the soup in the test does.
- Weakness: caware trusts the ambient load. Under-estimating it, or a noisy
  soup from few empties, brings leakage back. Next: (1) re-calibrate the
  ambient load per sample from the target type's own pseudobulk, (2) carry the
  soup's sampling variance into the likelihood, (3) pool soups across samples.
- Simulated data only; the contamination model matches the test's assumption.
  The k-means arm was uninformative here (types were perfectly separable).
