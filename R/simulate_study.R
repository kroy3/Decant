## simulate_study.R
## Multi-sample case/control study with sample-specific ambient RNA, for
## testing whether DE in a rare cell type is contaminated by an abundant cell
## type's disease genes. Same honesty caveat as the other simulators: it is
## built around the ambient mechanism the contamination-aware test models, so
## it tests sufficiency and failure modes, not real-world superiority.

#' Simulate a case/control single-cell study with sample-specific soup.
#'
#' @param n_per samples per condition.
#' @param G genes.
#' @param n_cells cells per sample.
#' @param rho_ctrl,rho_dis mean contamination in control / disease samples
#'   (disease tissue is often more fragile).
#' @param bcv biological coefficient of variation between samples.
#' @param n_de_rare,fc_rare true disease genes in the rare type and their fold.
#' @param n_de_abund,fc_abund disease genes in the abundant type A (all
#'   chosen among genes the rare type does not express).
#' @param frac cell-type proportions; must include types "A" and "R".
#' @param n_empty empty droplets per sample.
#' @param n_shared genes that are genuinely disease genes in BOTH A and R.
#' @param seed RNG seed.
#' @return list: X (genes x cells), label, sample, cond (per sample),
#'   empties (list per sample), de_rare, de_abund (leak-only genes), shared,
#'   rho_true (per sample), prof.
#' @export
simulate_ambient_study <- function(n_per = 6, G = 1000, n_cells = 700,
                                   rho_ctrl = 0.05, rho_dis = 0.12, bcv = 0.25,
                                   n_de_rare = 40, fc_rare = 2, n_de_abund = 80, fc_abund = 3,
                                   frac = c(A = 0.55, B = 0.2, C = 0.15, R = 0.10),
                                   n_empty = 2000, n_shared = 0, seed = 1) {
  if (!all(c("A", "R") %in% names(frac))) stop("`frac` must name types 'A' and 'R'", call. = FALSE)
  set.seed(seed)
  K <- length(frac); types <- names(frac)
  prof <- sapply(seq_len(K), function(k) {
    b <- stats::rgamma(G, 0.4, 1); idx <- ((k - 1) * 40 + 1):(k * 40)
    b[idx] <- b[idx] + 25; b / sum(b) })
  colnames(prof) <- types
  a_rank <- order(prof[, "A"] / (prof[, "R"] + 1e-6), decreasing = TRUE)
  de_abund <- a_rank[seq_len(n_de_abund)]
  cand <- setdiff(order(prof[, "R"], decreasing = TRUE)[1:400], de_abund)
  de_rare <- sample(cand, n_de_rare)
  dir_rare <- sample(c(-1, 1), n_de_rare, TRUE)
  shared <- if (n_shared > 0) de_abund[seq_len(n_shared)] else integer(0)
  if (n_shared > 0) {
    prof[shared, "R"] <- prof[shared, "A"] * 0.5
    prof <- sweep(prof, 2, colSums(prof), "/")
    de_rare <- c(de_rare, shared); dir_rare <- c(dir_rare, rep(1, n_shared))
  }
  cond <- rep(c("ctrl", "dis"), each = n_per); S <- length(cond)
  cells <- emp <- lab <- smp <- vector("list", S); rho_true <- numeric(S)
  for (s in seq_len(S)) {
    dis <- cond[s] == "dis"
    P <- prof * matrix(stats::rgamma(G * K, 1 / bcv^2, 1 / bcv^2), G, K)
    if (dis) {
      P[de_abund, "A"] <- P[de_abund, "A"] * fc_abund
      P[de_rare, "R"] <- P[de_rare, "R"] * fc_rare^dir_rare
    }
    P <- sweep(P, 2, colSums(P), "/")
    rho_s <- min(max(stats::rnorm(1, if (dis) rho_dis else rho_ctrl, 0.015), 0.005), 0.5)
    rho_true[s] <- rho_s
    soup <- as.numeric(P %*% frac); soup <- soup / sum(soup)
    lt <- sample(types, n_cells, TRUE, prob = frac)
    libs <- round(stats::rlnorm(n_cells, log(2500), 0.4))
    cells[[s]] <- vapply(seq_len(n_cells), function(c) {
      r <- stats::rbeta(1, rho_s * 40, (1 - rho_s) * 40); no <- stats::rbinom(1, libs[c], 1 - r)
      stats::rmultinom(1, no, P[, lt[c]])[, 1] + stats::rmultinom(1, libs[c] - no, soup)[, 1]
    }, numeric(G))
    lab[[s]] <- lt; smp[[s]] <- rep(s, n_cells)
    emp[[s]] <- vapply(stats::rpois(n_empty, 25), function(L) stats::rmultinom(1, L, soup)[, 1],
                       numeric(G))
  }
  X <- do.call(cbind, cells)
  rownames(X) <- paste0("g", seq_len(G)); colnames(X) <- paste0("c", seq_len(ncol(X)))
  emp <- lapply(emp, function(e) { rownames(e) <- rownames(X); e })
  list(X = X, label = unlist(lab), sample = unlist(smp), cond = factor(cond), empties = emp,
       de_rare = de_rare, de_abund = setdiff(de_abund, shared), shared = shared,
       rho_true = rho_true, prof = prof)
}

#' Simulate a large case/control cohort at the pseudobulk level.
#'
#' The same ambient mechanism as [simulate_ambient_study()], but drawing the
#' rare type's pseudobulk directly so that hundreds of samples are cheap, and
#' with optional departures from the negative binomial: per-gene outlier
#' samples (`out_p`, `out_fold`) or heavy-tailed biological noise
#' (`heavy_df`). These are the conditions under which Li et al. (2022,
#' Genome Biology 23:79) found NB-based DE to lose false-positive control in
#' population-scale data.
#'
#' @param n_per samples per condition.
#' @param G genes.
#' @param rho_ctrl,rho_dis mean contamination in control / disease samples.
#' @param bcv biological coefficient of variation between samples.
#' @param n_de_rare,fc_rare true disease genes in the rare type and their fold.
#' @param n_de_abund,fc_abund disease genes in the abundant type A.
#' @param frac cell-type proportions; must include types "A" and "R".
#' @param out_p probability that a gene's expression in a sample (in any cell
#'   type) is an outlier.
#' @param out_fold median fold of an outlier.
#' @param heavy_df if finite, biological noise is log-t with this many degrees
#'   of freedom instead of gamma.
#' @param n_batch,batch_fc number of rare-type genes shifted by a batch, and
#'   their fold. The genes are not disease genes.
#' @param batch_imbalance probability that a disease sample is in batch 1
#'   (a control sample is with one minus this), so batch confounds condition.
#' @param depth mean pseudobulk depth of the rare type per sample.
#' @param soup_umis empty-droplet UMIs per sample.
#' @param seed RNG seed.
#' @return list: Y (genes x samples pseudobulk of the rare type), empties
#'   (list per sample), cond, batch, de_rare, de_abund, batch_genes, rho_true.
#' @export
simulate_ambient_cohort <- function(n_per = 50, G = 2000, rho_ctrl = 0.05, rho_dis = 0.12,
                                    bcv = 0.25, n_de_rare = 40, fc_rare = 1.5,
                                    n_de_abund = 80, fc_abund = 3,
                                    frac = c(A = 0.55, B = 0.2, C = 0.15, R = 0.10),
                                    out_p = 0, out_fold = 10, heavy_df = Inf,
                                    n_batch = 0, batch_fc = 1.5, batch_imbalance = 0.7,
                                    depth = 175000, soup_umis = 50000, seed = 1) {
  if (!all(c("A", "R") %in% names(frac))) stop("`frac` must name types 'A' and 'R'", call. = FALSE)
  set.seed(seed)
  K <- length(frac); types <- names(frac)
  prof <- sapply(seq_len(K), function(k) {
    b <- stats::rgamma(G, 0.4, 1); idx <- ((k - 1) * 40 + 1):(k * 40)
    b[idx] <- b[idx] + 25; b / sum(b) })
  colnames(prof) <- types
  a_rank <- order(prof[, "A"] / (prof[, "R"] + 1e-6), decreasing = TRUE)
  de_abund <- a_rank[seq_len(n_de_abund)]
  cand <- setdiff(order(prof[, "R"], decreasing = TRUE)[seq_len(min(800, G))], de_abund)
  de_rare <- sample(cand, n_de_rare); dir_rare <- sample(c(-1, 1), n_de_rare, TRUE)
  cond <- rep(c("ctrl", "dis"), each = n_per); S <- length(cond)
  ## a batch that is over-represented in disease and shifts `n_batch` genes
  batch_genes <- if (n_batch > 0) sample(setdiff(cand, de_rare), n_batch) else integer(0)
  dir_batch <- sample(c(-1, 1), n_batch, TRUE)
  batch <- stats::rbinom(S, 1, ifelse(cond == "dis", batch_imbalance, 1 - batch_imbalance))
  Y <- matrix(0, G, S); emp <- vector("list", S); rho_true <- numeric(S)
  for (s in seq_len(S)) {
    dis <- cond[s] == "dis"
    noise <- if (is.finite(heavy_df)) exp(bcv * matrix(stats::rt(G * K, heavy_df), G, K))
             else matrix(stats::rgamma(G * K, 1 / bcv^2, 1 / bcv^2), G, K)
    P <- prof * noise
    if (out_p > 0) {
      hit <- matrix(stats::runif(G * K) < out_p, G, K)
      P[hit] <- P[hit] * exp(stats::rnorm(sum(hit), log(out_fold), 0.5))
    }
    if (dis) {
      P[de_abund, "A"] <- P[de_abund, "A"] * fc_abund
      P[de_rare, "R"] <- P[de_rare, "R"] * fc_rare^dir_rare
    }
    if (n_batch > 0 && batch[s] == 1) P[batch_genes, "R"] <- P[batch_genes, "R"] * batch_fc^dir_batch
    P <- sweep(P, 2, colSums(P), "/")
    rho_s <- min(max(stats::rnorm(1, if (dis) rho_dis else rho_ctrl, 0.015), 0.005), 0.5)
    rho_true[s] <- rho_s
    soup <- as.numeric(P %*% frac); soup <- soup / sum(soup)
    T_s <- depth * stats::rlnorm(1, 0, 0.3)
    Y[, s] <- stats::rpois(G, T_s * (1 - rho_s) * P[, "R"]) + stats::rpois(G, T_s * rho_s * soup)
    emp[[s]] <- matrix(stats::rmultinom(1, soup_umis, soup), ncol = 1)
  }
  rownames(Y) <- paste0("g", seq_len(G)); colnames(Y) <- paste0("s", seq_len(S))
  emp <- lapply(emp, function(e) { rownames(e) <- rownames(Y); e })
  list(Y = Y, empties = emp, cond = factor(cond), batch = factor(batch), de_rare = de_rare,
       de_abund = de_abund, batch_genes = batch_genes, rho_true = rho_true)
}
