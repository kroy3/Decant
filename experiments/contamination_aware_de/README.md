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

Run from this directory: `Rscript run1.R` (v1), `Rscript run2.R` (v1-v4 stress grid)
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

## Iteration: making the model robust (`methods2.R`, `run2.R`)

v1 trusted the ambient load and broke when it was under-estimated or the soup
was noisy. Changes, each kept only if it helped:

- **Self-calibration (kept):** re-estimate each sample's ambient fraction from
  the target type's own pseudobulk with the absent-gene estimator. Removes
  the dependence on upstream cell-level rho: leakage under a 30% rho
  under-estimate went 0.52 -> 0.01.
- **Fixed-rule soup pooling (rejected):** `pool_soup()` shrinks every sample
  toward the global soup by a fixed rule, erasing real sample-specific soup
  differences; null false positives rose to 0.107.
- **Empirical-Bayes soup (kept):** per gene, shrink by sampling variance vs
  between-sample variance, and carry the remaining soup variance into the
  likelihood. Few-empties leakage 0.46 -> 0.23 without harming rich samples.
- **TMM composition normalisation on the soup-removed counts (kept):** fixed
  null-gene inflation (0.13 -> 0.07) when a few genes change strongly.

| case (fp_leak) | naive | corrected | v1 | **v4** |
|---|---|---|---|---|
| base | 1.00 | 0.99 | 0.00 | 0.012 |
| negative control | 0.59 | 0.18 | 0.00 | 0.008 |
| rho under-estimated 30% | 1.00 | 0.99 | 0.52 | 0.012 |
| only 150 empty droplets | 1.00 | 0.99 | 0.46 | 0.23 |
| genes truly DE in both A and R | 1.00 | 0.98 | 0.02 | 0.043 |

v4: null FP rate 0.065-0.071 at nominal 0.05 (slightly liberal; LRT with
plug-in dispersion), FDP 0.13-0.15 at a 0.10 target, power 0.98-0.99, and
power 1.00 on genes truly DE in both types.

## What this does and does not show

- Shows: correcting counts first does NOT stop leakage (0.99); the residual
  after subtraction still tracks condition. Modelling the soup in the test does.
- Remaining weaknesses: (1) extreme soup scarcity (150 empties, ~3.7k soup
  UMIs) still leaks 23%; (2) slightly liberal tests (0.07 vs 0.05); a
  quasi-likelihood F-test, as edgeR uses, is the likely fix; (3) one test per
  gene via optim is slow for whole transcriptomes.
- Simulated data only; the contamination model matches the test's assumption.
  The k-means arm was uninformative here (types were perfectly separable).
