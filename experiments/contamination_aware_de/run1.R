source("sim.R"); source("methods.R")
res <- list()
for (s in 1:3) {
  st <- sim_study(seed = s); r <- run_all(st)
  cat(sprintf("seed %d  true rho ctrl=%.3f dis=%.3f | est R-cells ctrl=%.3f dis=%.3f\n", s,
      mean(st$rho_true[1:6]), mean(st$rho_true[7:12]), mean(r$rho_est[1:6]), mean(r$rho_est[7:12])))
  res[[s]] <- cbind(seed = s, score(st, r))
}
out <- do.call(rbind, res)
print(aggregate(cbind(fp_leak, fp_other_null, power, n_disc, fdp, leak_tested) ~ method, out, mean), digits = 3)
