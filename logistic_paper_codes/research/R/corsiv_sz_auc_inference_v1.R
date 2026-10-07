# Post-estimation only. No fitting, threshold selection, or solver loading.
cai_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}
cai_hash <- function(path) digest::digest(file = path, algo = "sha256")
cai_save <- function(x, path) {
  cai_assert(!file.exists(path), paste("Refusing overwrite:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("cai_atomic_", dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp))
  saveRDS(x, tmp)
  cai_assert(file.rename(tmp, path), "Atomic RDS rename failed")
}
cai_same_or_save <- function(x, path) {
  if (file.exists(path)) cai_assert(identical(x, readRDS(path)),
    paste("Existing output differs:", path)) else cai_save(x, path)
}
cai_manifest <- function(out, files) {
  paths <- file.path(out, files)
  data.frame(file = files, bytes = unname(file.info(paths)$size),
    sha256 = unname(vapply(paths, cai_hash, "")))
}
cai_check_manifest <- function(out, m) {
  cai_assert(nrow(m) > 0L && !anyDuplicated(m$file) &&
    !any(grepl("(^/|(^|/)\\.\\.(/|$))", m$file)), "Invalid manifest paths")
  for (i in seq_len(nrow(m))) {
    p <- file.path(out, m$file[i])
    cai_assert(file.exists(p) && !dir.exists(p) && file.info(p)$size == m$bytes[i] &&
      identical(cai_hash(p), m$sha256[i]), paste("Checksum failure:", p))
  }
  invisible(TRUE)
}
cai_inputs <- function(root) {
  out <- file.path(root, "outputs/study/corsiv_sz_external_v1")
  m <- readRDS(file.path(out, "CHECKPOINT_MANIFEST.rds"))
  cai_check_manifest(out, m)
  spec <- readRDS(file.path(out, "specification.rds"))
  complete <- readRDS(file.path(out, "COMPLETE.rds"))
  release <- readRDS(file.path(root, "release/corsiv_sz_v1/LOCAL_VALIDATED_corsiv_sz_external_v1.rds"))
  cai_assert(complete$validated_units == 6L && complete$method_rows == 6L &&
    identical(complete$scientific_signature, spec$scientific_signature) &&
    identical(release$production_signature, spec$scientific_signature), "Parent study not verified")
  methods <- c("Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)",
    "Logistic Group Elastic Net (adelie)", "Logistic Group Lasso (grpreg)",
    "Logistic Group MCP (grpreg)", "Logistic Group SCAD (grpreg)")
  z <- readRDS(file.path(out, "final/predictions.rds"))
  reported <- readRDS(file.path(out, "final/response_summary.rds"))
  ids <- spec$data_identity$test_sample_id
  cai_assert(nrow(z) == 4050L && length(ids) == 675L && !anyDuplicated(ids) &&
    setequal(unique(z$method), methods), "Unexpected predictions or methods")
  parts <- lapply(methods, function(method) {
    p <- z[z$method == method, , drop = FALSE]
    cai_assert(nrow(p) == 675L && !anyDuplicated(p$sample_id) &&
      setequal(p$sample_id, ids), "Unpaired test identities")
    p[match(ids, p$sample_id), , drop = FALSE]
  })
  y <- parts[[1L]]$y
  cai_assert(all(y %in% 0:1) && sum(y) == 353L && sum(y == 0L) == 322L &&
    all(vapply(parts, function(p) identical(p$y, y) &&
      all(is.finite(p$probability)) && all(p$probability >= 0 & p$probability <= 1), logical(1))),
    "Response alignment, class counts, or probability validation failed")
  rocs <- lapply(parts, function(p) pROC::roc(y, p$probability,
    levels = c(0, 1), direction = "<", percent = FALSE, quiet = TRUE, na.rm = FALSE))
  names(rocs) <- methods
  auc <- vapply(rocs, function(r) as.numeric(pROC::auc(r)), 0.0)
  cai_assert(identical(as.character(reported$method), methods) &&
    max(abs(auc-reported$auc)) < 1e-12, "pROC AUC differs from frozen reported AUC")
  list(rocs = rocs, methods = methods, auc = auc, ids = ids, y = y,
    identity = list(parent_signature = spec$scientific_signature,
      manifest_sha256 = cai_hash(file.path(out, "CHECKPOINT_MANIFEST.rds")),
      predictions_sha256 = cai_hash(file.path(out, "final/predictions.rds"))))
}
cai_spec <- function(root, input, stage) {
  files <- c("R/corsiv_sz_auc_inference_v1.R", "scripts/113_bootstrap_corsiv_sz_v1.R",
    "tests/test_corsiv_sz_auc_inference_v1.R", "CORSIV_SZ_AUC_INFERENCE_PROTOCOL_V1.md",
    "manuscript/corsiv_sz_auc_inference_methods_v1.tex", "manuscript/corsiv_sz_auc_inference_v1.bib")
  x <- list(schema = "corsiv_auc_inference_v1", stage = stage,
    run_id = if (stage == "smoke") "corsiv_sz_auc_inference_smoke_v1" else "corsiv_sz_auc_inference_b10000_v1",
    boot_n = if (stage == "smoke") 25L else 10000L, conf_level = .95,
    seed_ci = 1010601L, seed_comparison = 1010602L,
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    stratified = TRUE, paired = TRUE, alternative = "two.sided", direction = "<",
    levels = c(0, 1), percent = FALSE, adjustment = "holm",
    multiplicity_family = "five_sglasso_contrasts_separately_per_test_method",
    interpretation = "exploratory_conditional_on_trained_models_and_class_counts",
    model_fits = 0L, methods = input$methods, input_identity = input$identity,
    sources = cai_manifest(root, files),
    runtime = list(R = R.version.string, platform = R.version$platform,
      packages = vapply(c("pROC", "digest", "jsonlite", "knitr"),
        function(p) as.character(utils::packageVersion(p)), "")))
  x$signature <- digest::digest(x, algo = "sha256")
  x
}
cai_units <- function() c(paste0("auc_", 1:6), paste0("compare_", 2:6))
cai_compute <- function(unit, input, spec) {
  do.call(RNGkind, as.list(spec$rng_kind))
  index <- as.integer(sub("^.*_", "", unit))
  if (startsWith(unit, "auc_")) {
    # Identical seed preserves the same stratified resamples across all six CIs.
    set.seed(spec$seed_ci)
    pROC::ci.auc(input$rocs[[index]], method = "bootstrap",
      conf.level = spec$conf_level, boot.n = spec$boot_n,
      boot.stratified = TRUE, reuse.auc = TRUE)
  } else {
    set.seed(spec$seed_comparison)
    list(bootstrap = pROC::roc.test(input$rocs[[1]], input$rocs[[index]],
      method = "bootstrap", paired = TRUE, alternative = "two.sided",
      boot.n = spec$boot_n, boot.stratified = TRUE, reuse.auc = TRUE),
      delong = pROC::roc.test(input$rocs[[1]], input$rocs[[index]],
        method = "delong", paired = TRUE, alternative = "two.sided",
        conf.level = spec$conf_level, reuse.auc = TRUE))
  }
}
cai_validate_value <- function(value, unit, input, spec) {
  index <- as.integer(sub("^.*_", "", unit))
  if (startsWith(unit, "auc_")) {
    cai_assert(inherits(value, "ci.auc") && length(value) == 3L &&
      all(is.finite(value)) && all(value >= 0 & value <= 1) && all(diff(value) >= 0) &&
      identical(attr(value, "method"), "bootstrap") &&
      attr(value, "boot.n") == spec$boot_n && isTRUE(attr(value, "boot.stratified")) &&
      attr(value, "conf.level") == spec$conf_level &&
      abs(as.numeric(attr(value, "auc"))-input$auc[index]) < 1e-12,
      "Invalid pROC bootstrap AUC interval")
  } else {
    for (name in c("bootstrap", "delong")) {
      v <- value[[name]]
      cai_assert(inherits(v, "htest") && is.finite(v$p.value) &&
        v$p.value >= 0 && v$p.value <= 1 && grepl("correlated", v$method) &&
        identical(v$alternative, "two.sided") &&
        max(abs(as.numeric(v$estimate)-input$auc[c(1L,index)])) < 1e-12,
        paste("Invalid paired test:", name))
    }
    cai_assert(value$bootstrap$parameter[["boot.n"]] == spec$boot_n &&
      value$bootstrap$parameter[["boot.stratified"]] == 1 &&
      is.finite(as.numeric(value$bootstrap$statistic)) &&
      length(value$delong$conf.int) == 2L && all(is.finite(value$delong$conf.int)),
      "Invalid test controls or DeLong interval")
  }
  invisible(TRUE)
}
cai_path <- function(out, unit) file.path(out, "shards", paste0(unit, ".rds"))
cai_valid <- function(path, unit, input, spec) {
  if (!file.exists(path) || !file.exists(paste0(path, ".receipt.rds"))) return(FALSE)
  tryCatch({
    receipt <- readRDS(paste0(path, ".receipt.rds")); z <- readRDS(path)
    identical(cai_hash(path), receipt$sha256) && file.info(path)$size == receipt$bytes &&
      identical(z$signature, spec$signature) && identical(z$unit, unit) &&
      isTRUE(z$validated) && isTRUE(cai_validate_value(z$value, unit, input, spec))
  }, error = function(e) FALSE)
}
cai_progress <- function(out, spec, status, started, phase, current = NULL, durations = numeric()) {
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  completed <- sum(status %in% c("passed", "resumed")); failed <- sum(status == "failed")
  remaining <- names(status)[status == "pending"]
  type <- function(x) sub("_.*$", "", x)
  estimates <- vapply(remaining, function(u) {
    z <- durations[type(names(durations)) == type(u)]
    if (length(z)) mean(z) else NA_real_
  }, 0.0)
  eta <- if (length(estimates) && all(is.finite(estimates))) sum(estimates) else NULL
  now <- Sys.time(); deadline <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME", "")))
  p <- list(run_id = spec$run_id, slurm_job_id = Sys.getenv("SLURM_JOB_ID", "local"),
    phase = phase, work_unit = "six_auc_intervals_plus_five_paired_comparisons",
    total = length(status), completed = completed, running = as.integer(!is.null(current)),
    failed = failed, pending = length(status)-completed-failed-as.integer(!is.null(current)),
    percent_completed = 100*completed/length(status), elapsed_seconds = elapsed,
    throughput_units_per_second = if(elapsed>0) sum(status == "passed")/elapsed else 0,
    eta_seconds = eta, estimated_completion_utc = if(is.null(eta)) NULL else format(now+eta,tz="UTC",usetz=TRUE),
    eta_status = if(is.null(eta)) "estimating_or_unknown" else "current_invocation_type_specific",
    eta_scope = "remaining_resampling_units_including_current_excludes_final_aggregation_verification",
    slurm_remaining_seconds = if(is.finite(deadline)) max(0,deadline-as.numeric(now)) else NULL,
    heartbeat_utc = format(now,tz="UTC",usetz=TRUE), current_work_unit = current,
    most_recently_completed_work_unit = if(completed) tail(names(status)[status %in% c("passed","resumed")],1L) else NULL)
  tmp <- tempfile("cai_progress_",out)
  jsonlite::write_json(p,tmp,auto_unbox=TRUE,null="null",pretty=TRUE)
  cai_assert(file.rename(tmp,file.path(out,"progress.json")),"Progress rename failed")
  fields <- c("heartbeat_utc","phase","total","completed","running","failed","pending",
    "percent_completed","elapsed_seconds","throughput_units_per_second")
  row <- as.data.frame(p[fields]); row$eta_seconds <- if(is.null(eta)) NA_real_ else eta
  hist <- file.path(out,"progress.tsv")
  write.table(row,hist,sep="\t",row.names=FALSE,quote=FALSE,col.names=!file.exists(hist),append=file.exists(hist))
  cat(p$heartbeat_utc,phase,completed,"/",length(status),"validated; unit:",current,
    "ETA:",if(is.null(eta)) "unknown" else round(eta),"seconds\n"); flush.console()
}
cai_one <- function(unit, out, input, spec, work = function() cai_compute(unit,input,spec)) {
  path <- cai_path(out, unit)
  if(cai_valid(path,unit,input,spec)) return(invisible("resumed"))
  cai_assert(!file.exists(path) && !file.exists(paste0(path,".receipt.rds")), "Invalid existing shard; preserve and inspect")
  started <- Sys.time(); warnings <- character()
  value <- withCallingHandlers(work(), warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
  })
  cai_assert(!length(warnings),paste("pROC warning; inspect before inference:",paste(warnings,collapse="; ")))
  cai_validate_value(value,unit,input,spec)
  cai_save(list(signature=spec$signature,unit=unit,validated=TRUE,value=value,
    runtime_seconds=as.numeric(difftime(Sys.time(),started,units="secs"))),path)
  cai_save(list(bytes=file.info(path)$size,sha256=cai_hash(path)),paste0(path,".receipt.rds"))
  invisible("passed")
}
cai_tables <- function(out,input,spec) {
  get <- function(u) readRDS(cai_path(out,u))$value
  ci <- do.call(rbind,lapply(1:6,function(i) {
    z <- get(paste0("auc_",i))
    data.frame(method=input$methods[i],auc=input$auc[i],lower=z[1],upper=z[3],
      bootstrap_median=z[2],confidence_level=spec$conf_level,boot_n=spec$boot_n)
  }))
  comparisons <- do.call(rbind,lapply(2:6,function(i) {
    z <- get(paste0("compare_",i))
    data.frame(reference=input$methods[1],comparator=input$methods[i],
      auc_difference=input$auc[1]-input$auc[i],bootstrap_z=as.numeric(z$bootstrap$statistic),
      bootstrap_p=z$bootstrap$p.value,delong_p=z$delong$p.value,
      delong_difference_lower=unname(z$delong$conf.int[1]),
      delong_difference_upper=unname(z$delong$conf.int[2]))
  }))
  comparisons$bootstrap_p_holm <- p.adjust(comparisons$bootstrap_p,method="holm")
  comparisons$delong_p_holm <- p.adjust(comparisons$delong_p,method="holm")
  rownames(ci) <- rownames(comparisons) <- NULL
  list(auc_intervals=ci,paired_comparisons=comparisons)
}
cai_finalize <- function(out,input,spec) {
  tabs <- cai_tables(out,input,spec)
  for(n in names(tabs)) cai_same_or_save(tabs[[n]],file.path(out,"final",paste0(n,".rds")))
  # Manuscript-ready table material, generated from the same stored results.
  for(n in names(tabs)) {
    path <- file.path(out,"final",paste0(n,".tex"))
    lines <- as.character(knitr::kable(tabs[[n]],format="latex",booktabs=TRUE,digits=6,row.names=FALSE))
    if(file.exists(path)) cai_assert(identical(readLines(path,warn=FALSE),strsplit(lines,"\n",fixed=TRUE)[[1]]),"Table differs")
    else writeLines(lines,path)
  }
  files <- c("specification.rds",list.files(file.path(out,"shards"),full.names=FALSE),
    list.files(file.path(out,"final"),full.names=FALSE))
  files <- c(files[1],paste0("shards/",list.files(file.path(out,"shards"))),
    paste0("final/",list.files(file.path(out,"final"))))
  cai_same_or_save(cai_manifest(out,files),file.path(out,"MANIFEST.rds"))
  cai_same_or_save(list(signature=spec$signature,units=11L,auc_rows=6L,comparison_rows=5L,
    model_fits=0L,manifest_sha256=cai_hash(file.path(out,"MANIFEST.rds"))),file.path(out,"COMPLETE.rds"))
}
cai_verify <- function(out,input,spec) {
  cai_assert(identical(readRDS(file.path(out,"specification.rds")),spec),"Specification changed")
  cai_assert(all(vapply(cai_units(),function(u) cai_valid(cai_path(out,u),u,input,spec),logical(1))),"Incomplete inference units")
  m <- readRDS(file.path(out,"MANIFEST.rds")); cai_check_manifest(out,m)
  done <- readRDS(file.path(out,"COMPLETE.rds"))
  cai_assert(identical(done$signature,spec$signature) && done$units==11L && done$model_fits==0L &&
    identical(done$manifest_sha256,cai_hash(file.path(out,"MANIFEST.rds"))),"Completion validation failed")
  expected <- cai_tables(out,input,spec)
  for(n in names(expected)) cai_assert(identical(expected[[n]],readRDS(file.path(out,"final",paste0(n,".rds")))),"Summary reconstruction failed")
  invisible(TRUE)
}
