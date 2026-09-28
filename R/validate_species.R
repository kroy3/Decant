## validate_species.R
## Experimental ground truth from species-mixing (e.g. 10x human/mouse "hgmm")
## experiments. In a human cell every mouse read is contamination, so the
## cross-species reads measure rho directly -- no simulator involved.
##
## For a human cell c with T_c total counts and m_c mouse counts, and a soup
## whose mouse share is s_mouse (measured from the empty droplets):
##     rho_c = (m_c / T_c) / s_mouse
## i.e. the visible cross-species contamination scaled up by the part of the
## soup that is invisible (human-derived soup landing in a human cell).
##
## CAUTION -- species mixing is an EASY case for absent-gene estimators: the
## other genome is a perfect set of absent genes. A method that scores well only
## on the full matrix has shown little. species_mix_hidden() removes that crutch:
## it keeps one species' cells and ONLY that species' genes, so the estimator
## must find rho from within-species structure (as on real single-species
## data), while the truth still comes from the hidden cross-species reads.

#' Species of each gene from a species-prefixed feature name.
#'
#' Cell Ranger names features in multi-genome references as
#' `<genome>_<id>` (e.g. `GRCh38_ENSG...`, `mm10___ENSMUSG...`).
#'
#' @param genes character vector of feature names.
#' @return character vector of genome labels (the prefix before the first "_").
#' @export
gene_species <- function(genes) {
  sp <- sub("_.*$", "", genes)
  if (length(unique(sp)) != 2)
    stop("expected exactly two genome prefixes, found: ",
         paste(utils::head(unique(sp), 5), collapse = ", "), call. = FALSE)
  sp
}

#' Ground-truth contamination from a two-species mixture.
#'
#' @param cells genes x cells counts (both genomes' genes).
#' @param empties genes x droplets counts from empty droplets.
#' @param species length-nrow genome label per gene (see [gene_species()]).
#' @param max_minor cells whose minority-species fraction exceeds this are
#'   called doublets and excluded. Contamination alone rarely exceeds ~0.2 of
#'   the VISIBLE (cross-species) fraction.
#' @return data.frame per cell: barcode, species, total, own, other,
#'   other_frac, rho_true, doublet. Attribute `"soup_share"`: soup fraction
#'   per genome.
#' @export
species_mix_truth <- function(cells, empties, species, max_minor = 0.25) {
  cells <- .as_counts(cells, "cells"); empties <- .as_counts(empties, "empties")
  if (length(species) != nrow(cells)) stop("length(species) must equal nrow(cells)", call. = FALSE)
  g <- sort(unique(species))
  by_sp <- vapply(g, function(s) .col_sums(cells[species == s, , drop = FALSE]), numeric(ncol(cells)))
  if (is.null(dim(by_sp))) by_sp <- matrix(by_sp, nrow = 1, dimnames = list(NULL, g))
  soup <- vapply(g, function(s) sum(empties[species == s, , drop = FALSE]), numeric(1))
  soup_share <- soup / sum(soup)
  total <- rowSums(by_sp)
  call <- g[max.col(by_sp, ties.method = "first")]
  own <- by_sp[cbind(seq_along(call), match(call, g))]
  other <- total - own
  other_frac <- other / pmax(total, 1)
  other_share <- soup_share[ifelse(call == g[1], g[2], g[1])]
  out <- data.frame(barcode = colnames(cells), species = call, total = total,
                    own = own, other = other, other_frac = other_frac,
                    rho_true = pmin(other_frac / other_share, 1),
                    doublet = other_frac > max_minor,
                    stringsAsFactors = FALSE, row.names = NULL)
  attr(out, "soup_share") <- soup_share
  out
}

#' Score decontamination on a species-mixing experiment (full matrix).
#'
#' Doublets are excluded from all metrics.
#'
#' @param cells,empties,species as for [species_mix_truth()].
#' @param corrected corrected counts (same shape as `cells`).
#' @param rho_hat per-cell rho estimate (NULL to skip rho metrics).
#' @param truth optional precomputed [species_mix_truth()] result.
#' @return one-row data.frame: rho bias / RMSE / correlation, the fraction of
#'   cross-species (known-contaminant) counts removed, and own-species removal
#'   relative to the expected own-species contamination (1 = right amount,
#'   above 1 = over-correction destroying real signal).
#' @export
score_species_mix <- function(cells, empties, species, corrected, rho_hat = NULL,
                              truth = NULL) {
  if (is.null(truth)) truth <- species_mix_truth(cells, empties, species)
  keep <- !truth$doublet
  share <- attr(truth, "soup_share")
  cells <- .as_counts(cells); corrected <- .as_counts(corrected)
  g <- names(share)
  removed_by_sp <- vapply(g, function(s)
    .col_sums(cells[species == s, , drop = FALSE]) -
      .col_sums(corrected[species == s, , drop = FALSE]), numeric(ncol(cells)))
  own_idx <- match(truth$species, g)
  own_removed <- removed_by_sp[cbind(seq_along(own_idx), own_idx)]
  other_removed <- rowSums(removed_by_sp) - own_removed
  own_expected <- truth$rho_true * share[truth$species] * truth$total
  out <- data.frame(
    n_cells = sum(keep),
    rho_true_mean = mean(truth$rho_true[keep]),
    cross_removed = sum(other_removed[keep]) / max(sum(truth$other[keep]), 1),
    own_removed_ratio = sum(own_removed[keep]) / max(sum(own_expected[keep]), 1e-12))
  if (!is.null(rho_hat)) {
    r <- as.numeric(rho_hat)[keep]; t <- truth$rho_true[keep]
    out$rho_hat_mean <- mean(r)
    out$rho_bias <- mean(r) - mean(t)
    out$rho_rmse <- sqrt(mean((r - t)^2))
    out$rho_cor <- if (stats::sd(r) > 0) stats::cor(r, t) else NA_real_
  }
  out
}

#' Hidden-species view of a mixture: one species' cells, only its genes.
#'
#' The estimator sees no cross-species genes, so it must work as on ordinary
#' single-species data. Truth is rescaled to what is estimable from these
#' genes: contamination among the kept genes' counts,
#'     rho_hidden = rho_true * s_own * T / own.
#'
#' @param cells,empties,species as for [species_mix_truth()].
#' @param keep_species which genome to keep (default: the majority genome).
#' @param truth optional precomputed [species_mix_truth()] result.
#' @return list(cells, empties, rho_true, barcodes).
#' @export
species_mix_hidden <- function(cells, empties, species, keep_species = NULL, truth = NULL) {
  if (is.null(truth)) truth <- species_mix_truth(cells, empties, species)
  if (is.null(keep_species)) keep_species <- names(which.max(table(truth$species[!truth$doublet])))
  share <- attr(truth, "soup_share")
  sel <- which(truth$species == keep_species & !truth$doublet)
  gsel <- species == keep_species
  rho_hidden <- truth$rho_true[sel] * share[keep_species] * truth$total[sel] / truth$own[sel]
  list(cells = cells[gsel, sel, drop = FALSE], empties = empties[gsel, , drop = FALSE],
       rho_true = pmin(rho_hidden, 1), barcodes = truth$barcode[sel],
       species = keep_species)
}
