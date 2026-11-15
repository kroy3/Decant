#!/usr/bin/env Rscript
## scripts/check_user_guide.R
## Runs every R code block in the user guide (wiki/*.md) and checks that the
## output shown in the guide (lines starting with "#>") is what the code
## prints. Installation commands, and blocks whose first line starts with
## "# Requires", are skipped. Blocks of one chapter run in order in one
## session, so later blocks can use earlier objects.
## Usage, from the repository root:  Rscript scripts/check_user_guide.R

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
})
norm <- function(x) gsub("[[:space:]]+", " ", trimws(paste(x, collapse = "\n")))

blocks_of <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^```r\\s*$", lines); ends <- grep("^```\\s*$", lines)
  lapply(starts, function(s) lines[(s + 1):(min(ends[ends > s]) - 1)])
}

run_block <- function(block, env) {
  if (grepl("^# Requires", block[1]) || any(grepl("install(\\.packages|_github|::install)", block)))
    return(list(status = "skipped"))
  expected <- sub("^#> ?", "", grep("^#>", block, value = TRUE))
  code <- grep("^#>", block, value = TRUE, invert = TRUE)
  got <- utils::capture.output({
    for (e in parse(text = code)) {
      v <- withVisible(eval(e, env))
      if (v$visible) print(v$value)
    }
  }, type = "output")
  if (!length(expected)) return(list(status = "ran"))
  ok <- identical(norm(expected), norm(got))
  list(status = if (ok) "match" else "MISMATCH", expected = expected, got = got)
}

files <- sort(list.files("wiki", pattern = "^[0-9].*\\.md$", full.names = TRUE))
bad <- 0
for (f in files) {
  env <- new.env(parent = globalenv())
  set.seed(1)
  res <- tryCatch(suppressMessages(suppressWarnings(lapply(blocks_of(f), run_block, env = env))),
                  error = function(e) { bad <<- bad + 1; cat("ERROR in", f, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  st <- vapply(res, `[[`, "", "status")
  cat(sprintf("%-50s %s\n", basename(f), paste(names(table(st)), table(st), sep = "=", collapse = " ")))
  for (r in res[st == "MISMATCH"]) {
    bad <- bad + 1
    cat("  expected:\n", paste0("    ", r$expected, "\n"), "  got:\n", paste0("    ", r$got, "\n"), sep = "")
  }
}
if (bad) { cat(bad, "problem(s)\n"); quit(status = 1) } else cat("All code blocks run and match their shown output.\n")
