## Falsification simulation for contamination-aware DE.
## Cell-level counts -> per-sample empties -> rare-type pseudobulk.
sim_study <- function(n_per = 6, G = 1000, n_cells = 700, seed = 1,
                      rho_ctrl = 0.05, rho_dis = 0.12, bcv = 0.25,
                      n_de_rare = 40, fc_rare = 2, n_de_abund = 80, fc_abund = 3,
                      frac = c(A = 0.55, B = 0.2, C = 0.15, R = 0.10), n_empty = 2000, n_shared = 0) {
  set.seed(seed)
  K <- length(frac); types <- names(frac)
  prof <- sapply(seq_len(K), function(k) {
    b <- rgamma(G, 0.4, 1); idx <- ((k - 1) * 40 + 1):(k * 40); b[idx] <- b[idx] + 25; b / sum(b) })
  colnames(prof) <- types
  ## abundant-type disease genes: A-specific markers + A-high genes (not expressed in R)
  a_rank <- order(prof[, "A"] / (prof[, "R"] + 1e-6), decreasing = TRUE)
  de_abund <- a_rank[1:n_de_abund]
  ## rare-type true DE: genes reasonably expressed in R, not in the A set
  cand <- setdiff(order(prof[, "R"], decreasing = TRUE)[1:400], de_abund)
  de_rare <- sample(cand, n_de_rare)
  dir_rare <- sample(c(-1, 1), n_de_rare, TRUE)
  ## genes that GENUINELY change in both A and R (true R DE that also leaks from A)
  shared <- if (n_shared > 0) de_abund[seq_len(n_shared)] else integer(0)
  if (n_shared > 0) {
    prof[shared, "R"] <- prof[shared, "A"] * 0.5            # R expresses them natively too
    prof <- sweep(prof, 2, colSums(prof), "/")
    de_rare <- c(de_rare, shared); dir_rare <- c(dir_rare, rep(1, n_shared))
  }
  cond <- rep(c("ctrl", "dis"), each = n_per); S <- length(cond)
  cells <- list(); emp <- list(); lab <- list(); smp <- list(); rho_true <- numeric(S)
  for (s in seq_len(S)) {
    dis <- cond[s] == "dis"
    P <- prof * matrix(rgamma(G * K, 1 / bcv^2, 1 / bcv^2), G, K)       # sample biological var
    if (dis) {
      P[de_abund, "A"] <- P[de_abund, "A"] * fc_abund
      P[de_rare, "R"] <- P[de_rare, "R"] * fc_rare^dir_rare
    }
    P <- sweep(P, 2, colSums(P), "/")
    rho_s <- min(max(rnorm(1, if (dis) rho_dis else rho_ctrl, 0.015), 0.005), 0.5)
    rho_true[s] <- rho_s
    soup <- as.numeric(P %*% frac); soup <- soup / sum(soup)
    lt <- sample(types, n_cells, TRUE, prob = frac)
    libs <- round(rlnorm(n_cells, log(2500), 0.4))
    X <- vapply(seq_len(n_cells), function(c) {
      r <- rbeta(1, rho_s * 40, (1 - rho_s) * 40); no <- rbinom(1, libs[c], 1 - r)
      rmultinom(1, no, P[, lt[c]])[, 1] + rmultinom(1, libs[c] - no, soup)[, 1] }, numeric(G))
    cells[[s]] <- X; lab[[s]] <- lt; smp[[s]] <- rep(s, n_cells)
    emp[[s]] <- vapply(rpois(n_empty, 25), function(L) rmultinom(1, L, soup)[, 1], numeric(G))
  }
  X <- do.call(cbind, cells); rownames(X) <- paste0("g", seq_len(G))
  emp <- lapply(emp, function(e) { rownames(e) <- rownames(X); e })
  list(X = X, label = unlist(lab), sample = unlist(smp), cond = factor(cond), empties = emp,
       de_rare = de_rare, de_abund = setdiff(de_abund, shared), shared = shared, rho_true = rho_true, prof = prof)
}
