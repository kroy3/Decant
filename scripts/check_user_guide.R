#!/usr/bin/env Rscript
## scripts/check_user_guide.R
## Runs every R code block in the user guide (wiki/*.md) and checks that the
## output shown in the guide (lines starting with "#>") is what the code
## prints. Installation commands, and blocks whose first line starts with
## "# Requires", are skipped. Blocks of one chapter run in order in one
## session, so later blocks can use earlier objects.
## Usage, from the repository root:
##   Rscript scripts/check_user_guide.R            check only
##   Rscript scripts/check_user_guide.R --update   also rewrite the "#>" lines of
##                                                 mismatching blocks with the real
##                                                 output; then re-read the prose
##                                                 around them, which may quote
##                                                 the old numbers

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Decant)
})
update <- "--update" %in% commandArgs(trailingOnly = TRUE)
norm <- function(x) gsub("[[:space:]]+", " ", trimws(paste(x, collapse = "\n")))

block_ranges <- function(lines) {
  starts <- grep("^```r\\s*$", lines); ends <- grep("^```\\s*$", lines)
  lapply(starts, function(s) c(s + 1, min(ends[ends > s]) - 1))
}

## Returns the status, the shown and the real output, and the block rewritten
## with the real output placed after the line that produced it.
run_block <- function(block, env) {
  if (grepl("^# Requires", block[1]) || any(grepl("install(\\.packages|_github|::install)", block)))
    return(list(status = "skipped"))
  expected <- sub("^#> ?", "", grep("^#>", block, value = TRUE))
  code <- grep("^#>", block, value = TRUE, invert = TRUE)
  exprs <- parse(text = code, keep.source = TRUE)
  last_line <- vapply(attr(exprs, "srcref"), function(r) r[3], 0)
  out_after <- vector("list", length(code))
  for (i in seq_along(exprs)) {
    o <- utils::capture.output({
      v <- withVisible(eval(exprs[[i]], env))
      if (v$visible) print(v$value)
    }, type = "output")
    out_after[[last_line[i]]] <- c(out_after[[last_line[i]]], o)
  }
  got <- unlist(out_after)
  rewritten <- unlist(lapply(seq_along(code), function(j)
    c(code[j], if (length(out_after[[j]])) paste0("#> ", out_after[[j]]))))
  rewritten <- sub("^#> $", "#>", rewritten)
  if (!length(expected)) return(list(status = "ran"))
  ok <- identical(norm(expected), norm(got))
  list(status = if (ok) "match" else "MISMATCH", expected = expected, got = got, rewritten = rewritten)
}

files <- sort(list.files("wiki", pattern = "^[0-9].*\\.md$", full.names = TRUE))
bad <- 0
for (f in files) {
  env <- new.env(parent = globalenv())
  set.seed(1)
  lines <- readLines(f, warn = FALSE); ranges <- block_ranges(lines)
  res <- tryCatch(suppressMessages(suppressWarnings(
    lapply(ranges, function(r) run_block(lines[r[1]:r[2]], env)))),
    error = function(e) { bad <<- bad + 1; cat("ERROR in", f, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  st <- vapply(res, `[[`, "", "status")
  cat(sprintf("%-50s %s\n", basename(f), paste(names(table(st)), table(st), sep = "=", collapse = " ")))
  for (r in res[st == "MISMATCH"]) {
    bad <- bad + 1
    cat("  expected:\n", paste0("    ", r$expected, "\n"), "  got:\n", paste0("    ", r$got, "\n"), sep = "")
  }
  if (update && any(st == "MISMATCH")) {
    for (b in rev(which(st == "MISMATCH"))) {
      r <- ranges[[b]]
      lines <- c(lines[seq_len(r[1] - 1)], res[[b]]$rewritten, lines[(r[2] + 1):length(lines)])
    }
    writeLines(lines, f)
    cat("  updated", sum(st == "MISMATCH"), "block(s) in", basename(f), "\n")
  }
}
if (bad) {
  cat(bad, "problem(s)", if (update) "(blocks updated; re-run to confirm)" else "", "\n")
  quit(status = 1)
} else cat("All code blocks run and match their shown output.\n")
