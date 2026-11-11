# Changelog

## Decant (development version)

- New `validation/de_comparators.R`: contamination-aware DE against
  correct-then-test with the official SoupX, DecontX, FastCAR and
  CellBender (each run only if installed). Wrappers live in
  `validation/comparators.R`.
- README: the claim that correcting counts first does not stop DE
  leakage is narrowed. It holds for Decant’s default `subtract` rule;
  official SoupX followed by edgeR removes almost all leakage on the
  simulator.
- New vignette, “Contamination-aware differential expression”.
- pkgdown site configuration and deployment workflow.
- `DESCRIPTION`: `URL`, `BugReports`, `biocViews`, `VignetteBuilder`.
