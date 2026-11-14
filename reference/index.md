# Package index

## All functions

- [`DecantDE()`](https://kroy3.github.io/Decant/reference/DecantDE.md) :
  Contamination-aware DE for one cell type of a multi-sample Seurat
  object.

- [`RunDecant()`](https://kroy3.github.io/Decant/reference/RunDecant.md)
  : Run Decant on a Seurat object.

- [`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md)
  : Contamination-aware differential expression on pseudobulk counts.

- [`ambient_global()`](https://kroy3.github.io/Decant/reference/ambient_global.md)
  : Global soup profile: normalised aggregate of the empty-droplet pool.

- [`ambient_structured()`](https://kroy3.github.io/Decant/reference/ambient_structured.md)
  : Structured soup profile (lysis diagnostic).

- [`cluster_basis()`](https://kroy3.github.io/Decant/reference/cluster_basis.md)
  : Derive a program basis (cluster mean profiles) from the observed
  cells.

- [`correct_counts()`](https://kroy3.github.io/Decant/reference/correct_counts.md)
  : Remove ambient counts from each cell (mass-conserving).

- [`correct_splice_aware()`](https://kroy3.github.io/Decant/reference/correct_splice_aware.md)
  : Splice-aware layer decontamination.

- [`correct_total_then_split()`](https://kroy3.github.io/Decant/reference/correct_total_then_split.md)
  : Baseline to beat: correct the TOTAL, then split the removal by each
  cell's own observed unspliced/spliced ratio (what you get if you
  decontaminate total and naively apportion to layers).

- [`de_ambient_aware()`](https://kroy3.github.io/Decant/reference/de_ambient_aware.md)
  : Decant: test condition on OBSERVED data with the ambient load as a
  covariate, so condition-correlated ambient is absorbed instead of
  mistaken for signal.

- [`de_naive()`](https://kroy3.github.io/Decant/reference/de_naive.md) :
  Naive: subtract ambient point estimate, then test condition per gene.

- [`decant()`](https://kroy3.github.io/Decant/reference/decant.md) : Run
  Decant.

- [`decant_de_experimental()`](https://kroy3.github.io/Decant/reference/decant_de_experimental.md)
  :

  EXPERIMENTAL and OFF by default: failed its gate (made ambient-driven
  DE false positives worse via covariate collinearity). Kept only so the
  negative result is reproducible. Superseded by
  [`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md)
  /
  [`DecantDE()`](https://kroy3.github.io/Decant/reference/DecantDE.md).

- [`decontx_em()`](https://kroy3.github.io/Decant/reference/decontx_em.md)
  : DecontX-style decontamination (Bayesian-mixture EM,
  reimplementation).

- [`estimate_rho()`](https://kroy3.github.io/Decant/reference/estimate_rho.md)
  : Legacy v0 per-cell rho estimator (superseded).

- [`estimate_rho_allelic()`](https://kroy3.github.io/Decant/reference/estimate_rho_allelic.md)
  : Allelic rho estimate. Fraction of mismatched-allele reads = rho \*
  P(soup read is from another donor). With pool abundances a_d, for a
  cell of donor d that probability is (1 - a_d). Invert to recover rho.

- [`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md)
  : Cluster-pooled, absent-gene rho estimator with per-cell shrinkage.

- [`gate_allelic()`](https://kroy3.github.io/Decant/reference/gate_allelic.md)
  : GATE: allelic rho error vs a marker-style proxy (here: assuming
  uniform pool, the naive 1/n_donors correction a non-allelic method
  would implicitly make).

- [`gate_ambient_de()`](https://kroy3.github.io/Decant/reference/gate_ambient_de.md)
  : GATE: does contamination-aware DE stop ambient leakage while staying
  calibrated and powerful?

- [`gate_ambient_de_cohort()`](https://kroy3.github.io/Decant/reference/gate_ambient_de_cohort.md)
  : GATE: does contamination-aware DE stay calibrated in large cohorts
  whose counts are not negative binomial?

- [`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md)
  : GATE: which mass-conserving correction rule is closest to truth?

- [`gate_de()`](https://kroy3.github.io/Decant/reference/gate_de.md) :
  GATE: false-positive rate on ambient-trap genes (should be low for
  Decant, inflated for naive) while preserving power on true_de genes.
  Also reports the NULL-gene false-positive rate: if that is not near
  alpha the simulator or normalisation is broken and the gate result is
  meaningless.

- [`gate_diagnose()`](https://kroy3.github.io/Decant/reference/gate_diagnose.md)
  : GATE: with a clean (oracle) program basis, do recovered lysis
  weights track the true lysis weights? Reported as rank correlation;
  clustering noise lowers this in real use.

- [`gate_hier()`](https://kroy3.github.io/Decant/reference/gate_hier.md)
  : GATE: for empty-poor samples, is the pooled soup closer to truth
  than the independent estimate? For empty-rich samples, are they about
  equal (pooling should not hurt)?

- [`gate_rho()`](https://kroy3.github.io/Decant/reference/gate_rho.md) :
  GATE: is the cluster/absent-gene rho estimator calibrated, and better
  than the legacy marker-ratio estimator?

- [`gate_splice()`](https://kroy3.github.io/Decant/reference/gate_splice.md)
  : GATE: does splice-aware beat total-then-split on per-layer recovery,
  and does it correctly NOT help when ambient is not splice-distinct
  (negative control)?

- [`gene_species()`](https://kroy3.github.io/Decant/reference/gene_species.md)
  : Species of each gene from a species-prefixed feature name.

- [`pool_soup()`](https://kroy3.github.io/Decant/reference/pool_soup.md)
  : Empirical-Bayes-style shrinkage of per-sample soup profiles.

- [`print(`*`<decant>`*`)`](https://kroy3.github.io/Decant/reference/print.decant.md)
  : Print a decant result.

- [`quick_labels()`](https://kroy3.github.io/Decant/reference/quick_labels.md)
  : Default clustering used when the caller supplies none.

- [`read_10x_counts()`](https://kroy3.github.io/Decant/reference/read_10x_counts.md)
  : Read a 10x Genomics matrix directory into a sparse matrix.

- [`run_benchmark()`](https://kroy3.github.io/Decant/reference/run_benchmark.md)
  : Global vs structured soup benchmark.

- [`score_correction()`](https://kroy3.github.io/Decant/reference/score_correction.md)
  : Score a corrected matrix against ground truth.

- [`score_species_mix()`](https://kroy3.github.io/Decant/reference/score_species_mix.md)
  : Score decontamination on a species-mixing experiment (full matrix).

- [`simulate_allelic()`](https://kroy3.github.io/Decant/reference/simulate_allelic.md)
  : Simulate a pooled multi-genotype experiment with allele-resolved SNP
  genes.

- [`simulate_ambient_cohort()`](https://kroy3.github.io/Decant/reference/simulate_ambient_cohort.md)
  : Simulate a large case/control cohort at the pseudobulk level.

- [`simulate_ambient_study()`](https://kroy3.github.io/Decant/reference/simulate_ambient_study.md)
  : Simulate a case/control single-cell study with sample-specific soup.

- [`simulate_de()`](https://kroy3.github.io/Decant/reference/simulate_de.md)
  : Self-contained pseudobulk DE simulator with three gene classes: null
  : no difference anywhere true_de : differs in NATIVE expression
  between conditions (real signal) ambient_de: native identical, but
  SOUP differs between conditions (a trap)

- [`simulate_experiment()`](https://kroy3.github.io/Decant/reference/simulate_experiment.md)
  : Simulate a contaminated droplet experiment with known ground truth.

- [`simulate_multimodal()`](https://kroy3.github.io/Decant/reference/simulate_multimodal.md)
  : Multi-sample, splice-resolved contaminated experiment with ground
  truth.

- [`soup_eb()`](https://kroy3.github.io/Decant/reference/soup_eb.md) :
  Empirical-Bayes soup profiles per sample.

- [`soup_profile_error()`](https://kroy3.github.io/Decant/reference/soup_profile_error.md)
  : Soup-profile recovery error.

- [`species_mix_hidden()`](https://kroy3.github.io/Decant/reference/species_mix_hidden.md)
  : Hidden-species view of a mixture: one species' cells, only its
  genes.

- [`species_mix_truth()`](https://kroy3.github.io/Decant/reference/species_mix_truth.md)
  : Ground-truth contamination from a two-species mixture.

- [`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md)
  : Split a raw droplet matrix into cells and empty droplets.

- [`summarise_benchmark()`](https://kroy3.github.io/Decant/reference/summarise_benchmark.md)
  : Aggregate to mean +/- sd per (method, soup_bias).
