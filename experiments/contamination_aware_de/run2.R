source("sim.R"); source("methods.R")
cases <- list(
  base            = list(sim = list(), run = list()),
  neg_control     = list(sim = list(rho_dis = 0.05), run = list()),
  kmeans_clusters = list(sim = list(), run = list(clusters = "kmeans")),
  rho_under30     = list(sim = list(), run = list(rho_scale = 0.7)),
  rho_over30      = list(sim = list(), run = list(rho_scale = 1.3)),
  few_empties     = list(sim = list(n_empty = 150), run = list()),
  shared_true_de  = list(sim = list(n_shared = 20), run = list()))
all <- list()
for (nm in names(cases)) for (s in 1:3) {
  st <- do.call(sim_study, c(cases[[nm]]$sim, list(seed = s)))
  r <- do.call(run_all, c(list(st), cases[[nm]]$run))
  all[[length(all) + 1]] <- cbind(case = nm, seed = s, score(st, r))
}
out <- do.call(rbind, all)
agg <- aggregate(cbind(fp_leak, fp_other_null, power, n_disc, fdp) ~ case + method, out, mean)
sh <- aggregate(power_shared ~ case + method, out[out$case == "shared_true_de", ], mean)
agg <- agg[order(match(agg$case, names(cases)), agg$method), ]
print(agg, digits = 2, row.names = FALSE); cat("\nshared genes (true DE in R that also change in A):\n"); print(sh, digits = 2, row.names = FALSE)
