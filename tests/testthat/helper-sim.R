small_sim <- function(rho_mean = 0.1, seed = 42, ...) {
  simulate_experiment(n_genes = 300, n_types = 4, n_cells = 400, n_empty = 1500,
                      rho_mean = rho_mean, seed = seed, ...)
}
as_sparse <- function(x) methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix")
