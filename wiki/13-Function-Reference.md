# Chapter 13. Function Reference

Every exported function, grouped by task. Type `?function_name` in R for
the full documentation.

## 13.1 Main workflow

| Function | Purpose | Chapter |
|---|---|---|
| `decant()` | Estimate and remove ambient RNA | [5](05-Running-Decant) |
| `RunDecant()` | `decant()` for a Seurat object | [6](06-Seurat-Integration) |
| `read_10x_counts()` | Read a Cell Ranger matrix folder | [4](04-Preparing-Your-Data) |
| `split_droplets()` | Separate cells from empty droplets | [4](04-Preparing-Your-Data) |

## 13.2 Individual steps

| Function | Purpose | Chapter |
|---|---|---|
| `ambient_global()` | Soup profile from empty droplets | [7](07-Estimation-and-Correction) |
| `pool_soup()` | Per-sample soups with partial pooling | [8](08-Multi-Sample-Experiments) |
| `estimate_rho_cluster()` | Contamination per cell (default estimator) | [7](07-Estimation-and-Correction) |
| `estimate_rho()` | Older marker-ratio estimator, kept for comparison | |
| `estimate_rho_allelic()` | Contamination from genotypes | [10](10-Splice-Aware-and-Genotype-Modules) |
| `correct_counts()` | Remove contamination with a chosen rule | [7](07-Estimation-and-Correction) |
| `correct_splice_aware()` | Correct spliced and unspliced layers | [10](10-Splice-Aware-and-Genotype-Modules) |
| `correct_total_then_split()` | Reference method for splice-aware correction | |
| `quick_labels()` | Fallback k-means clustering | [5](05-Running-Decant) |
| `ambient_structured()`, `cluster_basis()` | Experimental soup model built from cell programs; off by default | |

## 13.3 Differential expression

| Function | Purpose | Chapter |
|---|---|---|
| `ambient_de()` | Contamination-aware test on pseudobulk counts | [9](09-Differential-Expression) |
| `DecantDE()` | `ambient_de()` for a Seurat object | [9](09-Differential-Expression) |
| `soup_eb()` | Per-sample soups shrunk by empirical Bayes, used by `ambient_de()` | |
| `de_naive()`, `de_ambient_aware()`, `decant_de_experimental()` | Earlier approaches; `decant_de_experimental()` failed its gate and is kept only for reproducibility | |

| `chirp()` | Ambient-aware ligand-receptor calls between clusters (experimental) | [14](14-Cell-Cell-Communication) |

## 13.4 Simulation

| Function | Simulates |
|---|---|
| `simulate_experiment()` | One contaminated sample with known truth |
| `simulate_multimodal()` | Several samples with spliced and unspliced layers |
| `simulate_allelic()` | Pooled donors with allele counts |
| `simulate_ambient_study()` | A case/control study for differential expression |
| `simulate_ambient_cohort()` | Large cohorts, with outliers and batch effects |
| `simulate_de()` | Earlier differential expression simulator |

## 13.5 Evaluation

| Function | Purpose |
|---|---|
| `score_correction()` | Compare corrected counts with the truth |
| `soup_profile_error()` | Compare a soup estimate with the true soup |
| `run_benchmark()`, `summarise_benchmark()` | Run and summarise simulated benchmarks |
| `gate_rho()`, `gate_correction()`, `gate_hier()`, `gate_splice()`, `gate_allelic()`, `gate_diagnose()`, `gate_de()`, `gate_ambient_de()`, `gate_ambient_de_cohort()` | Falsification gates ([Chapter 12](12-Validation-and-Benchmarks)) |
| `species_mix_truth()`, `species_mix_hidden()`, `score_species_mix()`, `gene_species()` | Ground truth and scoring for species mixtures |
| `decontx_em()` | Simplified DecontX-style model used in simulation comparisons; for real comparisons use the official DecontX |

---

**Previous:** [Chapter 12, Validation and Benchmarks](12-Validation-and-Benchmarks) | **Next:** [Chapter 14, Cell-Cell Communication](14-Cell-Cell-Communication)
