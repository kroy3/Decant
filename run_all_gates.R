#!/usr/bin/env Rscript
## run_all_gates.R -- the honest scorecard. Each module must pass its own
## falsification gate to be enabled by default in Decant. Failures are reported,
## not hidden.
for (f in list.files("Decant/R", pattern = "\\.R$", full.names = TRUE)) source(f)

cat("=========================================================\n")
cat(" DECANT MODULE GATES (each must earn its place)\n")
cat("=========================================================\n\n")

cat("[GAP 1] Splice-aware layer decontamination\n");      gate_splice();   cat("\n")
cat("[GAP 3] Hierarchical multi-sample soup pooling\n");  gate_hier();     cat("\n")
cat("[GAP 5] Allelic / genotype-aware rho estimation\n"); gate_allelic();  cat("\n")
cat("[GAP 4] Structured-soup lysis diagnostic\n");        gate_diagnose(); cat("\n")
cat("[GAP 2] Decontamination-aware differential expression\n"); gate_de(); cat("\n")

cat("=========================================================\n")
cat(" Modules that fail their gate ship OFF by default.\n")
cat("=========================================================\n")
