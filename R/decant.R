## decant.R -- the assembled pipeline.
## Composes ONLY the modules that passed their gate. Each capability switches on
## when the data that justifies it is present, and the mass-conservation
## guarantee holds end to end. The DE module is intentionally excluded by default
## (it failed its gate); it is reachable as decant_de_experimental() with a loud
## warning, not wired into the main path.

#' Run Decant.
#'
#' @param obs_total genes x cells total counts (required).
#' @param empties EITHER a genes x droplets matrix (single sample) OR a list of
#'   such matrices, one per sample, which triggers hierarchical pooling [GAP 3].
#' @param sample_of length-cells sample index (required if empties is a list).
#' @param obs_unspliced,obs_spliced genes x cells layers; if both given with
#'   per-sample empty layers, correction is splice-aware [GAP 1].
#' @param empties_unspliced,empties_spliced layer-resolved empties (matrix or
#'   list per sample), required for splice-aware mode.
#' @param allelic optional list(own, other, donor, abund) for genotype-aware rho
#'   [GAP 5]; overrides marker-based rho when present.
#' @param basis optional genes x K program basis; if given, returns a lysis
#'   diagnostic [GAP 4].
#' @return list: $corrected (+ $corrected_unspliced/$corrected_spliced if splice
#'   mode), $rho, $ambient (per sample), $lysis (if basis), $modules (what ran).
decant <- function(obs_total, empties, sample_of = NULL,
                   obs_unspliced = NULL, obs_spliced = NULL,
                   empties_unspliced = NULL, empties_spliced = NULL,
                   allelic = NULL, basis = NULL) {
  modules <- character(0)
  multi <- is.list(empties)
  G <- nrow(obs_total); N <- ncol(obs_total)

  ## ---- ambient profile (hierarchical pooling if multi-sample) [GAP 3] ----
  if (multi) {
    stopifnot(!is.null(sample_of))
    pp <- pool_soup(empties)
    soup_by_sample <- pp$pooled
    modules <- c(modules, "hierarchical_ambient")
  } else {
    soup_by_sample <- matrix(ambient_global(empties), ncol = 1)
    sample_of <- rep(1L, N)
  }
  soup_for_cell <- soup_by_sample[, sample_of, drop = FALSE]   # G x N

  ## ---- per-cell rho (allelic if available) [GAP 5] ----
  if (!is.null(allelic)) {
    rho <- estimate_rho_allelic(allelic)
    modules <- c(modules, "allelic_rho")
  } else {
    rho <- estimate_rho(obs_total, rowMeans(soup_by_sample))
    modules <- c(modules, "marker_rho")
  }

  ## ---- correction ----
  splice_mode <- !is.null(obs_unspliced) && !is.null(obs_spliced) &&
                 !is.null(empties_unspliced) && !is.null(empties_spliced)
  out <- list(rho = rho, ambient = soup_by_sample)

  if (splice_mode) {
    eu <- if (is.list(empties_unspliced)) do.call(cbind, empties_unspliced) else empties_unspliced
    es <- if (is.list(empties_spliced))   do.call(cbind, empties_spliced)   else empties_spliced
    sa <- correct_splice_aware(obs_unspliced, obs_spliced, eu, es, rho)
    out$corrected_unspliced <- sa$unspliced
    out$corrected_spliced   <- sa$spliced
    out$corrected <- sa$unspliced + sa$spliced
    modules <- c(modules, "splice_aware_correction")
  } else {
    ## per-cell correction with that cell's (possibly pooled) soup
    T_c <- colSums(obs_total)
    expected <- sweep(soup_for_cell, 2, rho * T_c, "*")
    corr <- obs_total - expected; corr[corr < 0] <- 0
    corr <- pmin(corr, obs_total)
    out$corrected <- corr
    modules <- c(modules, "global_correction")
  }

  ## ---- lysis diagnostic [GAP 4] ----
  if (!is.null(basis)) {
    st <- ambient_structured(if (multi) do.call(cbind, empties) else empties, basis)
    out$lysis <- st$lysis
    modules <- c(modules, "lysis_diagnostic")
  }

  ## hard guarantee check
  stopifnot(all(out$corrected >= 0), all(out$corrected <= obs_total + 1e-6))
  out$modules <- modules
  out
}

#' EXPERIMENTAL and OFF by default: failed its gate (made ambient-driven DE
#' false positives worse via covariate collinearity). Kept only so the negative
#' result is reproducible. Do not use for real inference.
decant_de_experimental <- function(...) {
  warning("decant_de_experimental FAILED its benchmark gate (higher false-positive ",
          "rate than naive). Provided for reproducibility only; do not use.")
  de_ambient_aware(...)
}
