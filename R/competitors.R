## competitors.R
## Faithful, clearly-labelled REIMPLEMENTATIONS of competing algorithms, for
## controlled benchmarking only. These are NOT the official packages (which
## would not compile in this sandbox). They reproduce the published algorithm
## structure so the comparison is about the METHOD, not a wrapper call.
##
## Before publishing any comparison, re-run against the official SoupX (CRAN) and
## DecontX (Bioconductor 'celda') packages. These reimplementations are for
## internal method development.

#' DecontX-style decontamination (Bayesian-mixture EM, reimplementation).
#'
#' Model (Yang et al. 2020): each cell belongs to a cluster k. Its observed
#' counts are a mixture of native expression phi_k and a contamination
#' distribution eta_k (expression originating from OTHER clusters in the
#' experiment). Per-cell contamination theta_c and the profiles are fit by EM.
#' Crucially, this uses NO empty droplets -- contamination is inferred from
#' population structure alone. Returns a decontaminated (native) count matrix.
#'
#' @param observed genes x cells counts.
#' @param labels integer cluster labels per cell (length = ncol).
#' @param iters EM iterations.
#' @return list with $corrected (genes x cells) and $theta (per-cell contamination).
decontx_em <- function(observed, labels, iters = 30, eps = 1e-10) {
  G <- nrow(observed); N <- ncol(observed)
  K <- length(unique(labels)); ulab <- sort(unique(labels))
  cs <- colSums(observed); cs[cs == 0] <- 1

  ## init native profiles phi_k = normalised cluster sums
  phi <- vapply(ulab, function(k) {
    m <- rowSums(observed[, labels == k, drop = FALSE]); m / sum(m)
  }, numeric(G))
  sizes <- vapply(ulab, function(k) sum(labels == k), numeric(1))

  ## init contamination profiles eta_k = normalised aggregate of OTHER clusters
  eta <- vapply(seq_len(K), function(ki) {
    w <- sizes; w[ki] <- 0
    m <- as.numeric(phi %*% w); m / sum(m)
  }, numeric(G))

  theta <- rep(0.1, N)
  lab_idx <- match(labels, ulab)

  for (it in seq_len(iters)) {
    native <- matrix(0, G, N)
    new_theta <- numeric(N)
    for (c in seq_len(N)) {
      k <- lab_idx[c]
      x <- observed[, c]
      num <- theta[c] * eta[, k]
      den <- num + (1 - theta[c]) * phi[, k] + eps
      r <- num / den                      # responsibility: contamination
      contam_counts <- sum(x * r)
      new_theta[c] <- contam_counts / cs[c]
      native[, c] <- x * (1 - r)          # native counts (mass-conserving)
    }
    theta <- pmin(pmax(new_theta, 0), 0.99)
    ## M-step: refresh phi from native counts, eta from updated phi
    phi <- vapply(seq_len(K), function(ki) {
      m <- rowSums(native[, lab_idx == ki, drop = FALSE]); s <- sum(m)
      if (s == 0) phi[, ki] else m / s
    }, numeric(G))
    eta <- vapply(seq_len(K), function(ki) {
      w <- sizes; w[ki] <- 0
      m <- as.numeric(phi %*% w); m / sum(m)
    }, numeric(G))
  }
  list(corrected = native, theta = theta)
}

#' Simple k-means labels for feeding DecontX-style EM (kept separate so the
#' clustering choice is explicit and shared across methods that need it).
quick_labels <- function(observed, k = 6, n_hvg = 200, seed = 1) {
  set.seed(seed)
  cs <- colSums(observed); cs[cs == 0] <- 1
  logn <- log1p(sweep(observed, 2, cs, "/") * 1e4)
  v <- apply(logn, 1, var)
  hvg <- order(v, decreasing = TRUE)[seq_len(min(n_hvg, nrow(logn)))]
  kmeans(t(logn[hvg, ]), centers = k, nstart = 5, iter.max = 50)$cluster
}
