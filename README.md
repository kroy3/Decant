# Decant

A benchmark-first toolkit for ambient RNA decontamination research in
single-cell and single-nucleus RNA-seq. Decant is assembled from modules, and a
module is wired into the default pipeline ONLY after it passes a falsification
gate against ground truth. This is not a finished CellBender replacement and
does not claim to be one. Its value is that every capability in it has been shown
to earn its place, and the things that did not are documented as failures.

> "Decant" is a working placeholder name. Check CRAN, Bioconductor, and GitHub
> for collisions before any release.

## The design rule

Earlier work in this repo established, against ground truth, that putting a
smarter model on the same count matrix (a fancier soup profile) does not move the
sensitivity/specificity frontier. So every module here brings in ORTHOGONAL
information to break the ambient-vs-native identifiability problem, and each is
gated. Run the scorecard yourself:

```bash
Rscript run_all_gates.R
```

## Module scorecard (from run_all_gates.R)

| Gap | Module | Gate result | Default |
|-----|--------|-------------|---------|
| 1 | Splice-aware layer decontamination | PASS. Halves unspliced-layer error when ambient is splice-distinct; correctly neutral when it is not (negative control holds). Uniquely outputs corrected spliced/unspliced layers for RNA velocity. | ON |
| 3 | Hierarchical multi-sample soup pooling | PASS. Reduces soup error for empty-poor samples; neutral for empty-rich. | ON |
| 5 | Allelic / genotype-aware rho | PASS. Near-ground-truth rho in pooled designs. (Gate contrast was soft; the absolute accuracy is the real point.) | ON when allelic data present |
| 4 | Structured-soup lysis diagnostic | PASS as QC only. Recovers which populations seed the soup; it is a diagnostic readout, never a corrector (that use was falsified earlier). | ON as diagnostic |
| 2 | Decontamination-aware differential expression | FAIL. Made ambient-driven DE false positives worse than naive (collinearity between the ambient covariate and the condition). The problem is real and severe; this fix is not the answer. | OFF (experimental, kept for reproducibility) |

## The assembled pipeline

`decant()` switches each passing capability on when the data that justifies it is
present (layer matrices trigger splice-aware mode; a list of per-sample empties
triggers pooling; an allelic table triggers genotype-aware rho; a program basis
adds the lysis diagnostic). The mass-conservation guarantee (never fabricate
counts) holds end to end and is asserted at the end of every run.

```r
res <- decant(obs_total, empties = list_of_per_sample_empties,
              sample_of = cell_sample_index,
              obs_unspliced = U, obs_spliced = S,
              empties_unspliced = eU, empties_spliced = eS,
              basis = program_basis)
res$corrected            # corrected total
res$corrected_unspliced  # corrected layers (for velocity)
res$lysis                # which populations seed the soup
res$modules              # what actually ran
```

## Known limitations (stated, not hidden)

- The default marker-based rho estimator is poorly calibrated and can overestimate
  contamination roughly 2x on multi-sample data, causing over-correction. Use the
  allelic estimator where genotype data exists, or replace this estimator. This is
  the single weakest component and the clearest target for the next iteration.
- The lysis diagnostic was validated with an oracle basis; real k-means clustering
  degrades it.
- All results here are on simulated data. Simulations can be rigged; these are
  sufficiency tests and negative controls, not evidence against CellBender.

## The non-negotiable next step

Nothing here is a real-world claim until it is re-run on experimental ground
truth: species-mixing and genotype-mixing datasets, scored with this same metric
suite, with the official SoupX, DecontX (celda), and CellBender as comparators.
The harness is built so those wrappers drop in beside the simulator.

## License

MIT. See LICENSE.
