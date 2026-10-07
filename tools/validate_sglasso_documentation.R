#!/usr/bin/env Rscript
# Read-only checks of the local documentation build and its recorded source protection.
site <- normalizePath(yaml::read_yaml("_pkgdown.yml")$destination, mustWork = TRUE)
validation <- file.path("output/validation",
  sub("^sglasso_", "sglasso_documentation_", basename(site)), "build")
hash <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
record <- readRDS(file.path(validation, "BUILD_ACCEPTED.rds"))
stopifnot(isTRUE(record$accepted), record$article_count == 5L,
          identical(record$production_runs, 0L), !record$reference_examples_executed)
protected <- record$protected_files
stopifnot(identical(unname(vapply(protected$file, hash, character(1))), protected$sha256))
m <- read.csv(file.path(validation, "SITE_MANIFEST.csv"), stringsAsFactors = FALSE)
paths <- file.path(site, m$file)
stopifnot(all(file.exists(paths)), all(file.info(paths)$size == m$bytes),
          identical(unname(vapply(paths, hash, character(1))), m$sha256))
cff <- yaml::read_yaml("CITATION.cff")$`preferred-citation`
stopifnot(cff$year == 2026, cff$start == 1L, cff$end == 23L,
          cff$`date-published` == "2026-09-17", cff$status == "advance-online",
          is.null(cff$volume), is.null(cff$issue))
cit <- utils::readCitationFile("inst/CITATION", meta = as.list(read.dcf("DESCRIPTION")[1L, ]))[[1L]]
stopifnot(cit$pages == "1--23", cit$doi == cff$doi,
          grepl("17 September 2026", cit$note, fixed = TRUE),
          identical(hash("CITATION.cff"), hash(file.path(site, "CITATION.cff"))))
for (p in c("index.html", "authors.html", "articles/reproducibility.html", "news/index.html")) {
  txt <- gsub("[[:space:]]+", " ", xml2::xml_text(xml2::read_html(file.path(site, p))))
  stopifnot(grepl("17 September 2026", txt, fixed = TRUE),
            grepl("10.1080/00031305.2026.2709494", txt, fixed = TRUE))
}
files <- list.files(site, pattern = "[.]html$", recursive = TRUE, full.names = TRUE)
checked <- 0L
for (f in files) {
  html <- xml2::read_html(f)
  for (attr in c("href", "src")) {
    nodes <- xml2::xml_find_all(html, if (attr == "href") ".//a[@href]|.//link[@href]" else ".//img[@src]|.//script[@src]")
    urls <- xml2::xml_attr(nodes, attr)
    for (u in urls) {
      if (!nzchar(u) || grepl("^([a-zA-Z]+:|//|#)", u)) next
      rel <- utils::URLdecode(sub("[?#].*$", "", u))
      path <- normalizePath(file.path(dirname(f), rel), mustWork = FALSE)
      if (dir.exists(path)) path <- file.path(path, "index.html")
      if (!file.exists(path)) stop("Broken local link in ", basename(f), ": ", u)
      checked <- checked + 1L
    }
  }
  imgs <- xml2::xml_find_all(html, ".//img")
  stopifnot(!anyNA(xml2::xml_attr(imgs, "alt")))
}
public <- c("README.Rmd", "README.md", "CITATION.cff", "NEWS.md", "_pkgdown.yml",
            list.files("vignettes/articles", pattern = "[.]Rmd$", full.names = TRUE))
stopifnot(!any(vapply(public, function(f) any(grepl("/Users/|/arf/|Documents/Codex|gwrs|GWR-ENET",
  readLines(f, warn = FALSE), ignore.case = TRUE)), logical(1))))
cat("VERIFIED:", length(files), "HTML pages;", checked,
    "local links/assets; manifest;", nrow(protected), "protected files unchanged.\n")
cat("Five guides; no identifying local paths; no gwrs content; no production runs.\n")
