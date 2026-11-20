# Chapter 13. Function Reference

Every exported function, grouped by task. Type `?function_name` in R for
the full documentation.

## 13.1 Main workflow

| Function | Purpose | Chapter |
|----|----|----|
| [`decant()`](https://kroy3.github.io/Decant/reference/decant.md) | Estimate and remove ambient RNA | [5](https://kroy3.github.io/Decant/articles/05-Running-Decant.md) |
| [`RunDecant()`](https://kroy3.github.io/Decant/reference/RunDecant.md) | [`decant()`](https://kroy3.github.io/Decant/reference/decant.md) for a Seurat object | [6](https://kroy3.github.io/Decant/articles/06-Seurat-Integration.md) |
| [`read_10x_counts()`](https://kroy3.github.io/Decant/reference/read_10x_counts.md) | Read a Cell Ranger matrix folder | [4](https://kroy3.github.io/Decant/articles/04-Preparing-Your-Data.md) |
| [`split_droplets()`](https://kroy3.github.io/Decant/reference/split_droplets.md) | Separate cells from empty droplets | [4](https://kroy3.github.io/Decant/articles/04-Preparing-Your-Data.md) |

## 13.2 Individual steps

| Function | Purpose | Chapter |
|----|----|----|
| [`ambient_global()`](https://kroy3.github.io/Decant/reference/ambient_global.md) | Soup profile from empty droplets | [7](https://kroy3.github.io/Decant/articles/07-Estimation-and-Correction.md) |
| [`pool_soup()`](https://kroy3.github.io/Decant/reference/pool_soup.md) | Per-sample soups with partial pooling | [8](https://kroy3.github.io/Decant/articles/08-Multi-Sample-Experiments.md) |
| [`estimate_rho_cluster()`](https://kroy3.github.io/Decant/reference/estimate_rho_cluster.md) | Contamination per cell (default estimator) | [7](https://kroy3.github.io/Decant/articles/07-Estimation-and-Correction.md) |
| [`estimate_rho()`](https://kroy3.github.io/Decant/reference/estimate_rho.md) | Older marker-ratio estimator, kept for comparison |  |
| [`estimate_rho_allelic()`](https://kroy3.github.io/Decant/reference/estimate_rho_allelic.md) | Contamination from genotypes | [10](https://kroy3.github.io/Decant/articles/10-Splice-Aware-and-Genotype-Modules.md) |
| [`correct_counts()`](https://kroy3.github.io/Decant/reference/correct_counts.md) | Remove contamination with a chosen rule | [7](https://kroy3.github.io/Decant/articles/07-Estimation-and-Correction.md) |
| [`correct_splice_aware()`](https://kroy3.github.io/Decant/reference/correct_splice_aware.md) | Correct spliced and unspliced layers | [10](https://kroy3.github.io/Decant/articles/10-Splice-Aware-and-Genotype-Modules.md) |
| [`correct_total_then_split()`](https://kroy3.github.io/Decant/reference/correct_total_then_split.md) | Reference method for splice-aware correction |  |
| [`quick_labels()`](https://kroy3.github.io/Decant/reference/quick_labels.md) | Fallback k-means clustering | [5](https://kroy3.github.io/Decant/articles/05-Running-Decant.md) |
| [`ambient_structured()`](https://kroy3.github.io/Decant/reference/ambient_structured.md), [`cluster_basis()`](https://kroy3.github.io/Decant/reference/cluster_basis.md) | Experimental soup model built from cell programs; off by default |  |

## 13.3 Differential expression

| Function | Purpose | Chapter |
|----|----|----|
| [`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md) | Contamination-aware test on pseudobulk counts | [9](https://kroy3.github.io/Decant/articles/09-Differential-Expression.md) |
| [`DecantDE()`](https://kroy3.github.io/Decant/reference/DecantDE.md) | [`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md) for a Seurat object | [9](https://kroy3.github.io/Decant/articles/09-Differential-Expression.md) |
| [`soup_eb()`](https://kroy3.github.io/Decant/reference/soup_eb.md) | Per-sample soups shrunk by empirical Bayes, used by [`ambient_de()`](https://kroy3.github.io/Decant/reference/ambient_de.md) |  |
| [`de_naive()`](https://kroy3.github.io/Decant/reference/de_naive.md), [`de_ambient_aware()`](https://kroy3.github.io/Decant/reference/de_ambient_aware.md), [`decant_de_experimental()`](https://kroy3.github.io/Decant/reference/decant_de_experimental.md) | Earlier approaches; [`decant_de_experimental()`](https://kroy3.github.io/Decant/reference/decant_de_experimental.md) failed its gate and is kept only for reproducibility |  |

[`chirp()`](https://kroy3.github.io/Decant/reference/chirp.md) \|
Ambient-awareligand-receptor calls between clusters (experimental) \|
[14](https://kroy3.github.io/Decant/articles/14-Cell-Cell-Communication.md)
\|

## 13.4 Simulation

| Function | Simulates |
|----|----|
| [`simulate_experiment()`](https://kroy3.github.io/Decant/reference/simulate_experiment.md) | One contaminated sample with known truth |
| [`simulate_multimodal()`](https://kroy3.github.io/Decant/reference/simulate_multimodal.md) | Several samples with spliced and unspliced layers |
| [`simulate_allelic()`](https://kroy3.github.io/Decant/reference/simulate_allelic.md) | Pooled donors with allele counts |
| [`simulate_ambient_study()`](https://kroy3.github.io/Decant/reference/simulate_ambient_study.md) | A case/control study for differential expression |
| [`simulate_ambient_cohort()`](https://kroy3.github.io/Decant/reference/simulate_ambient_cohort.md) | Large cohorts, with outliers and batch effects |
| [`simulate_de()`](https://kroy3.github.io/Decant/reference/simulate_de.md) | Earlier differential expression simulator |

## 13.5 Evaluation

| Function | Purpose |
|----|----|
| [`score_correction()`](https://kroy3.github.io/Decant/reference/score_correction.md) | Compare corrected counts with the truth |
| [`soup_profile_error()`](https://kroy3.github.io/Decant/reference/soup_profile_error.md) | Compare a soup estimate with the true soup |
| [`run_benchmark()`](https://kroy3.github.io/Decant/reference/run_benchmark.md), [`summarise_benchmark()`](https://kroy3.github.io/Decant/reference/summarise_benchmark.md) | Run and summarise simulated benchmarks |
| [`gate_rho()`](https://kroy3.github.io/Decant/reference/gate_rho.md), [`gate_correction()`](https://kroy3.github.io/Decant/reference/gate_correction.md), [`gate_hier()`](https://kroy3.github.io/Decant/reference/gate_hier.md), [`gate_splice()`](https://kroy3.github.io/Decant/reference/gate_splice.md), [`gate_allelic()`](https://kroy3.github.io/Decant/reference/gate_allelic.md), [`gate_diagnose()`](https://kroy3.github.io/Decant/reference/gate_diagnose.md), [`gate_de()`](https://kroy3.github.io/Decant/reference/gate_de.md), [`gate_ambient_de()`](https://kroy3.github.io/Decant/reference/gate_ambient_de.md), [`gate_ambient_de_cohort()`](https://kroy3.github.io/Decant/reference/gate_ambient_de_cohort.md) | Falsification gates ([Chapter 12](https://kroy3.github.io/Decant/articles/12-Validation-and-Benchmarks.md)) |
| [`species_mix_truth()`](https://kroy3.github.io/Decant/reference/species_mix_truth.md), [`species_mix_hidden()`](https://kroy3.github.io/Decant/reference/species_mix_hidden.md), [`score_species_mix()`](https://kroy3.github.io/Decant/reference/score_species_mix.md), [`gene_species()`](https://kroy3.github.io/Decant/reference/gene_species.md) | Ground truth and scoring for species mixtures |
| [`decontx_em()`](https://kroy3.github.io/Decant/reference/decontx_em.md) | Simplified DecontX-style model used in simulation comparisons; for real comparisons use the official DecontX |

------------------------------------------------------------------------

**Previous:** [Chapter 12, Validation and
Benchmarks](https://kroy3.github.io/Decant/articles/12-Validation-and-Benchmarks.md)
\| **Next:** [Chapter 14, Cell-Cell
Communication](https://kroy3.github.io/Decant/articles/14-Cell-Cell-Communication.md)
