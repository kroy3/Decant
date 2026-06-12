## The one property that survived every experiment: the correction can never
## fabricate counts. This test fails loudly if a future edit breaks that.
for (f in list.files("../R", pattern="\\.R$", full.names=TRUE)) source(f)
sim  <- simulate_experiment(n_cells = 300, n_empty = 1000, seed = 42)
soup <- ambient_global(sim$empty)
rho  <- estimate_rho(sim$observed, soup)
corr <- correct_counts(sim$observed, soup, rho)
stopifnot(all(corr >= 0))                 # never negative
stopifnot(all(corr <= sim$observed))      # never adds counts
stopifnot(sum(pmax(corr - sim$observed, 0)) == 0)
cat("PASS: mass-conservation guarantee holds (no fabricated counts)\n")
