portable_assert <- function (ok, message) 
if (!isTRUE(ok)) stop(message, call. = FALSE)

portable_hash <- function (p) 
digest::digest(file = p, algo = "sha256")

portable_inventory <- function (bundle) 
{
    f <- sort(list.files(bundle, recursive = TRUE, all.files = TRUE, no.. = TRUE))
    f <- f[!dir.exists(file.path(bundle, f)) & !f %in% c("MANIFEST.csv", "MANIFEST.sha256")]
    data.frame(file = f, bytes = unname(file.info(file.path(bundle, f))$size), sha256 = unname(vapply(file.path(bundle, 
        f), portable_hash, "")))
}

portable_require <- function (packages) 
{
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    portable_assert(!length(missing), paste("Missing installed dependencies (no automatic installation):", paste(missing, 
        collapse = ", ")))
}

portable_sources <- function (bundle) 
{
    f <- portable_inventory(bundle)
    f <- f[startsWith(f$file, "research/"), ]
    f$file <- substring(f$file, 10L)
    f
}

portable_simulation <- function (bundle, out, stage) 
{
    root <- file.path(bundle, "research")
    scope <- new.env(parent = globalenv())
    sys.source(file.path(root, "R/logistic_sglasso_workflow_v19.R"), scope)
    e <- scope$lsg_load_v19(root)
    base <- e$lsg_configuration_v7
    e$lsg_configuration_v7 <- function(stage) {
        cfg <- base(stage)
        if (stage == "smoke") 
            cfg$alpha_grid <- cfg$d_grid <- cfg$benchmark_alpha_grid <- c(0, 0.5, 1)
        else portable_assert(identical(cfg, paper_configuration(bundle, "simulation")$configuration), "Production simulation settings changed")
        cfg
    }
    e$lsg_source_files_v7 <- function() portable_sources(bundle)$file
    e$lsg_output_v7 <- function(root, version) file.path(out, version)
    environment(e$lsg_run_v19)$lsg_operational_preflight_v15 <- function(root, stage, cores, max_seconds, task_limit, 
        interval) {
        portable_assert(cores >= 1L && cores <= 56L && max_seconds > 0 && interval >= 1 && interval <= 60, "Invalid execution budget")
    }
    e
}

portable_external <- function (bundle, out, stage, version) 
{
    root <- file.path(bundle, "research")
    source(file.path(root, "R/corsiv_sz_analysis_v1.R"), local = environment())
    b <- csa_context(root)
    b$gse_load_data <- function(root, stage) portable_data(bundle, stage)
    b$gse_source_inventory <- function(root) portable_sources(bundle)
    b$gse_output <- function(root, version) file.path(out, version)
    e <- b$allb_load_environment_v7(root)
    parent.env(e) <- b
    spec <- b$gse_specification(root, e, stage, version)
    if (stage == "production") {
        original <- paper_configuration(bundle, "external")
        portable_assert(identical(spec$configuration, original$configuration) && identical(spec$fold, original$fold), 
            "Production external settings/folds changed")
    }
    list(b = b, e = e, spec = spec, data = b$gse_load_data(root, stage))
}

portable_output <- function (bundle, path) 
{
    portable_assert(nzchar(path), "An explicit output directory is required")
    portable_assert(startsWith(path, "/"), "Use an absolute output directory")
    portable_assert(!grepl("(^|/)\\.\\.(/|$)", path), "Parent traversal refused")
    existing <- path
    while (!dir.exists(existing)) existing <- dirname(existing)
    canonical <- file.path(normalizePath(existing), substring(path, nchar(existing) + 2L))
    if (dir.exists(path)) 
        canonical <- normalizePath(path)
    portable_assert(!identical(canonical, bundle) && !startsWith(canonical, paste0(bundle, "/")), "Output must be outside the immutable supplement")
    canonical
}

