#!/usr/bin/env Rscript
# Local documentation build only. Reference examples are intentionally disabled.
# Use an already installed, source-matching package; never install or compile here.
# Rscript --vanilla tools/build_sglasso_documentation.R --library=/path/to/library
args <- commandArgs(trailingOnly = TRUE)
finish_site <- "--finish-site" %in% args
citation_refresh <- "--citation-refresh" %in% args
if (citation_refresh && finish_site) stop("Choose one build mode.")
arg <- grep("^--library=", args, value = TRUE)
if (length(arg) != 1L) stop("Provide exactly one --library=/path/to/installed/library.")
lib <- normalizePath(sub("^--library=", "", arg), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
for (p in c("sglasso", "pkgdown", "rmarkdown", "knitr", "yaml", "digest", "xml2")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing installed package: ", p)
}
stopifnot(utils::packageVersion("sglasso") == "1.2.0")
stopifnot(normalizePath(dirname(find.package("sglasso"))) == lib)
root <- normalizePath(".")
config <- yaml::read_yaml("_pkgdown.yml")
site <- file.path(root, config$destination)
validation <- file.path(root, "output/validation",
                        sub("^sglasso_", "sglasso_documentation_", basename(site)), "build")
if (!finish_site && (file.exists(site) || file.exists(validation))) {
  stop("Site/validation path already exists. Preserve it; use an approved new version.")
}
if (finish_site) stopifnot(dir.exists(site),
                          !file.exists(file.path(validation, "BUILD_ACCEPTED.rds")))
if (citation_refresh) {
  previous_arg <- grep("^--previous-site=", args, value = TRUE)
  if (length(previous_arg) != 1L) stop("Provide one --previous-site=path for a citation-only refresh.")
  previous <- normalizePath(sub("^--previous-site=", "", previous_arg), mustWork = TRUE)
  previous_validation <- file.path(root, "output/validation",
    sub("^sglasso_", "sglasso_documentation_", basename(previous)), "build")
  previous_record <- readRDS(file.path(previous_validation, "BUILD_ACCEPTED.rds"))
  stopifnot(isTRUE(previous_record$accepted), previous != site)
  previous_manifest <- read.csv(file.path(previous_validation, "SITE_MANIFEST.csv"))
  previous_paths <- file.path(previous, previous_manifest$file)
  stopifnot(all(file.info(previous_paths)$size == previous_manifest$bytes),
    identical(unname(vapply(previous_paths, function(f)
      digest::digest(file = f, algo = "sha256", serialize = FALSE), character(1))),
      previous_manifest$sha256))
  dir.create(site, recursive = TRUE)
  copied <- list.files(previous, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  copied <- copied[basename(copied) != "CITATION.cff"]
  stopifnot(all(file.copy(copied, site, recursive = TRUE, overwrite = FALSE)))
}
hash <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
core <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
          list.files("src", pattern = "[.](cpp|h|hpp)$", full.names = TRUE))
legacy <- c("man/figures/README-example-1.png", "man/figures/README-example-2.png",
            "docs/help/gfortran_installation_guide.html")
protected <- c(core, legacy)
before <- vapply(protected, hash, character(1))
for (f in c("tools/render_sglasso_gallery.R", "tools/build_sglasso_documentation.R")) parse(f)
cff <- yaml::read_yaml("CITATION.cff")
stopifnot(cff$`cff-version` == "1.2.0", cff$type == "software",
          cff$version == "1.2.0", length(cff$authors) == 2L,
          cff$authors[[1]]$`family-names` == "Yüzbaşı",
          cff$`preferred-citation`$doi == "10.1080/00031305.2026.2709494")
form <- yaml::read_yaml(".github/ISSUE_TEMPLATE/bug-report.yml")
stopifnot(length(form$body) >= 4L)
articles <- list.files("vignettes/articles", pattern = "[.]Rmd$", full.names = TRUE)
stopifnot(length(articles) == 5L)
for (f in articles) {
  metadata <- rmarkdown::yaml_front_matter(f)
  stopifnot(nzchar(metadata$title))
}
gallery <- "output/validation/sglasso_documentation_v1/gallery"
m <- utils::read.csv(file.path(gallery, "MANIFEST.csv"), stringsAsFactors = FALSE)
gp <- file.path(gallery, m$file)
stopifnot(all(file.info(gp)$size == m$bytes),
          identical(unname(vapply(gp, hash, character(1))), m$sha256))
assets <- c("sglasso-quickstart-path.svg", "sglasso-external-roc.svg")
dest <- file.path("man/figures", assets)
for (i in seq_along(dest)) {
  if (file.exists(dest[i])) {
    stopifnot(identical(hash(dest[i]), hash(file.path(gallery, assets[i]))))
  } else stopifnot(file.copy(file.path(gallery, assets[i]), dest[i], overwrite = FALSE))
}
dir.create(validation, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = file.path(root, "output/validation/sglasso_documentation_v1/dependency_cache"))
if (!finish_site) {
  rmarkdown::render("README.Rmd", output_format = rmarkdown::github_document(html_preview = FALSE),
                    quiet = TRUE, envir = new.env(parent = globalenv()))
  if (citation_refresh) {
    # Reuse validated pages; this guide has no executable model-fitting chunks.
    pkgdown::build_home(root, preview = FALSE, quiet = FALSE)
    pkgdown::build_article("articles/reproducibility", pkg = root,
                          new_process = FALSE, quiet = FALSE)
  } else pkgdown::build_site(pkg = root, examples = FALSE, run_dont_run = FALSE,
                      install = FALSE, devel = FALSE, new_process = FALSE,
                      preview = FALSE, quiet = FALSE)
}
if (finish_site || citation_refresh) {
  # Resume only the nonnumerical stages after all articles/reference pages exist.
  stopifnot(all(file.exists(file.path(site, "articles", paste0(
    tools::file_path_sans_ext(basename(articles)), ".html")))),
    file.exists(file.path(site, "reference/sglasso.html")))
  pkgdown::build_news(root, preview = FALSE)
  pkg <- pkgdown::as_pkgdown(root)
  pkgdown:::build_sitemap(pkg)
  pkgdown:::build_redirects(pkg)
  pkgdown:::build_search(pkg)
  stopifnot(file.copy("pkgdown/extra.css", file.path(site, "extra.css"), overwrite = TRUE))
}
# README links must work both on GitHub (source files) and on the built website.
index_file <- file.path(site, "index.html")
index <- xml2::read_html(index_file)
links <- xml2::xml_find_all(index, ".//a[@href]")
for (link in links) {
  href <- xml2::xml_attr(link, "href")
  if (grepl("^vignettes/articles/.*[.]Rmd$", href)) {
    href <- sub("[.]Rmd$", ".html", sub("^vignettes/articles/", "articles/", href))
  } else if (href %in% c("tools/render_sglasso_gallery.R", "logistic_paper_codes/README.txt")) {
    href <- paste0("https://github.com/byuzbasi/sglasso/blob/main/", href)
  } else if (href == "paper_codes/") {
    href <- "https://github.com/byuzbasi/sglasso/tree/main/paper_codes"
  } else if (href %in% c("NEWS.md", "NEWS.html")) href <- "news/index.html"
  xml2::xml_set_attr(link, "href", href)
}
xml2::write_html(index, index_file)
if (file.exists(file.path(site, "CITATION.cff"))) {
  stopifnot(identical(hash("CITATION.cff"), hash(file.path(site, "CITATION.cff"))))
} else stopifnot(file.copy("CITATION.cff", file.path(site, "CITATION.cff"), overwrite = FALSE))
stopifnot(file.exists(file.path(site, "index.html")),
          all(file.exists(file.path(site, "articles", paste0(
            tools::file_path_sans_ext(basename(articles)), ".html")))))
after <- vapply(protected, hash, character(1))
stopifnot(identical(before, after))
saveRDS(list(schema = "sglasso_documentation_build_v1", accepted = TRUE,
             reference_examples_executed = FALSE, production_runs = 0L,
             citation_only_refresh = citation_refresh,
             protected_files = data.frame(file = protected, sha256 = before),
             article_count = length(articles), R = R.version.string,
             package_versions = vapply(c("sglasso", "pkgdown", "rmarkdown", "knitr"),
               function(p) as.character(utils::packageVersion(p)), character(1))),
        file.path(validation, "BUILD_ACCEPTED.rds"))
files <- list.files(site, recursive = TRUE, full.names = TRUE)
utils::write.csv(data.frame(file = substring(files, nchar(site) + 2L),
  bytes = file.info(files)$size, sha256 = vapply(files, hash, character(1))),
  file.path(validation, "SITE_MANIFEST.csv"), row.names = FALSE)
cat("VERIFIED: five articles; core and original graphics unchanged; no production run.\n")
cat("Local site:", site, "\n")
