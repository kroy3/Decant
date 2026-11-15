#!/usr/bin/env Rscript
## scripts/build_user_guide.R
## The user guide is written once, as GitHub wiki pages in wiki/. This script
## copies those pages into vignettes/articles/ as pkgdown articles, so the
## same text also appears on the package website. Run it after editing
## wiki/, from the repository root:
##   Rscript scripts/build_user_guide.R
## Wiki links such as (03-Quick-Start) become (03-Quick-Start.html); the wiki
## home page becomes the article "00-Preface".

out_dir <- file.path("vignettes", "articles")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pages <- sort(list.files("wiki", pattern = "^(Home|[0-9]{2}-.*)\\.md$", full.names = TRUE))

slug_of <- function(path) {
  b <- sub("\\.md$", "", basename(path))
  if (b == "Home") "00-Preface" else b
}

for (p in pages) {
  lines <- readLines(p, warn = FALSE)
  title <- sub("^# +", "", lines[grep("^# ", lines)[1]])
  body <- lines[-grep("^# ", lines)[1]]
  body <- gsub("\\]\\(Home\\)", "](00-Preface.html)", body)
  body <- gsub("\\]\\(([0-9]{2}-[A-Za-z-]+)\\)", "](\\1.html)", body)
  yaml <- c("---", sprintf("title: \"%s\"", gsub("\"", "\\\\\"", title)), "---", "",
            "<!-- Generated from wiki/ by scripts/build_user_guide.R. Edit the wiki page, not this file. -->", "")
  writeLines(c(yaml, body), file.path(out_dir, paste0(slug_of(p), ".Rmd")))
}
cat(sprintf("Wrote %d articles to %s\n", length(pages), out_dir))
