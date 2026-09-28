#' Decant: benchmark-first ambient RNA decontamination
#'
#' Start with [decant()] for correction, [split_droplets()] /
#' [read_10x_counts()] to get real data in, and the `gate_*()` functions to
#' re-derive every default from ground truth.
#'
#' @keywords internal
#' @import Matrix
#' @importFrom methods as
#' @importFrom stats aggregate cor kmeans lm median quantile rbeta rbinom
#'   rgamma rlnorm rmultinom rpois runif sd var ppois
"_PACKAGE"
