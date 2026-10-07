#!/usr/bin/env Rscript
# One entry point. The default verifies the supplied evidence without fitting.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
stopifnot(length(script) == 1L)
# R's non-trailing --file argument encodes ASCII spaces as ~+~ on some builds.
# Keep a literal existing filename intact; decode only when lookup failed.
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed=TRUE)
bundle <- normalizePath(dirname(script), mustWork=TRUE)
source(file.path(bundle, "adapter/portable.R"))
portable_main(bundle, commandArgs(TRUE))
